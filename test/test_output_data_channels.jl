using Test
import HDF5
import OWENS

const OUTPUT_DATASET_NAMES = [
    "t",
    "aziHist",
    "OmegaHist",
    "OmegaDotHist",
    "gbHist",
    "gbDotHist",
    "gbDotDotHist",
    "FReactionHist",
    "FTwrBsHist",
    "genTorque",
    "genPower",
    "torqueDriveShaft",
    "uHist",
    "uHist_prp",
    "epsilon_x_hist",
    "epsilon_y_hist",
    "epsilon_z_hist",
    "kappa_x_hist",
    "kappa_y_hist",
    "kappa_z_hist",
    "massOwens",
    "stress_U",
    "SF_ult_U",
    "SF_buck_U",
    "stress_L",
    "SF_ult_L",
    "SF_buck_L",
    "stress_TU",
    "SF_ult_TU",
    "SF_buck_TU",
    "stress_TL",
    "SF_ult_TL",
    "SF_buck_TL",
    "topstrainout_blade_U",
    "topstrainout_blade_L",
    "topstrainout_tower_U",
    "topstrainout_tower_L",
    "topDamage_blade_U",
    "topDamage_blade_L",
    "topDamage_tower_U",
    "topDamage_tower_L",
]

@testset "outputData channel registry" begin
    channels = OWENS.output_data_channels()
    @test channels isa Vector{OWENS.ResultChannel}
    @test length(channels) == 41
    @test OWENS.output_data_channel_names() == OUTPUT_DATASET_NAMES
    @test [channel.name for channel in channels] == OUTPUT_DATASET_NAMES
    @test length(unique(OWENS.output_data_channel_names())) == 41

    time_channel = OWENS.output_data_channel("t")
    @test time_channel.units == "s"
    @test time_channel.dimensions == ["time"]
    @test time_channel.frame == "scalar"
    @test time_channel.association == "time"
    @test time_channel.source == "OWENS runtime"
    @test time_channel.sign_convention == "positive elapsed simulation time"

    reaction_channel = OWENS.output_data_channel("FReactionHist")
    @test reaction_channel.units == "N or N*m by DOF"
    @test reaction_channel.dimensions == ["time", "structural_dof"]
    @test reaction_channel.frame == "structural mesh"
    @test reaction_channel.sign_convention == "unknown"

    stress_channel = OWENS.output_data_channel("stress_U")
    @test stress_channel.units == "Pa"
    @test stress_channel.dimensions ==
          ["time", "span_station", "chord_station", "stress_component"]
    @test stress_channel.association == "blade upper laminate"

    @test_throws KeyError OWENS.output_data_channel("not_a_channel")
end

@testset "outputData writes channel metadata" begin
    mktempdir() do dir
        data_output = joinpath(dir, "unit.out")
        inputs = (; dataOutputFilename = data_output)

        vector = [1.0, 2.0, 3.0]
        dof_history = reshape(collect(1.0:18.0), 3, 6)
        strain_history = reshape(collect(1.0:24.0), 4, 2, 3)
        stress_history = reshape(collect(1.0:72.0), 3, 2, 4, 3)
        safety_factor = reshape(collect(1.0:24.0), 3, 2, 4)
        topstrain = reshape(collect(1.0:216.0), 3, 2, 4, 9)
        damage = reshape(collect(1.0:8.0), 2, 4)

        OWENS.outputData(;
            inputs,
            t = vector,
            aziHist = vector .+ 0.1,
            OmegaHist = vector .+ 0.2,
            OmegaDotHist = vector .+ 0.3,
            gbHist = vector .+ 0.4,
            gbDotHist = vector .+ 0.5,
            gbDotDotHist = vector .+ 0.6,
            FReactionHist = dof_history,
            FTwrBsHist = dof_history .+ 10.0,
            genTorque = vector .+ 0.7,
            genPower = vector .+ 0.8,
            torqueDriveShaft = vector .+ 0.9,
            uHist = dof_history .+ 20.0,
            uHist_prp = dof_history .+ 30.0,
            epsilon_x_hist = strain_history .+ 1.0,
            epsilon_y_hist = strain_history .+ 2.0,
            epsilon_z_hist = strain_history .+ 3.0,
            kappa_x_hist = strain_history .+ 4.0,
            kappa_y_hist = strain_history .+ 5.0,
            kappa_z_hist = strain_history .+ 6.0,
            massOwens = 12.5,
            stress_U = stress_history .+ 1.0,
            SF_ult_U = safety_factor .+ 1.0,
            SF_buck_U = safety_factor .+ 2.0,
            stress_L = stress_history .+ 2.0,
            SF_ult_L = safety_factor .+ 3.0,
            SF_buck_L = safety_factor .+ 4.0,
            stress_TU = stress_history .+ 3.0,
            SF_ult_TU = safety_factor .+ 5.0,
            SF_buck_TU = safety_factor .+ 6.0,
            stress_TL = stress_history .+ 4.0,
            SF_ult_TL = safety_factor .+ 7.0,
            SF_buck_TL = safety_factor .+ 8.0,
            topstrainout_blade_U = topstrain .+ 1.0,
            topstrainout_blade_L = topstrain .+ 2.0,
            topstrainout_tower_U = topstrain .+ 3.0,
            topstrainout_tower_L = topstrain .+ 4.0,
            topDamage_blade_U = damage .+ 1.0,
            topDamage_blade_L = damage .+ 2.0,
            topDamage_tower_U = damage .+ 3.0,
            topDamage_tower_L = damage .+ 4.0,
        )

        h5_file = joinpath(dir, "unit.h5")
        @test isfile(h5_file)

        HDF5.h5open(h5_file, "r") do file
            @test sort(keys(file)) == sort(OUTPUT_DATASET_NAMES)
            @test HDF5.read(file["t"]) == vector
            @test HDF5.read(file["massOwens"]) == 12.5
            @test size(HDF5.read(file["FReactionHist"])) == (3, 6)
            @test size(HDF5.read(file["epsilon_x_hist"])) == (4, 2, 3)
            @test size(HDF5.read(file["stress_U"])) == (3, 2, 4, 3)
            @test size(HDF5.read(file["topstrainout_blade_U"])) == (3, 2, 4, 9)
            @test size(HDF5.read(file["topDamage_blade_U"])) == (2, 4)

            for name in OUTPUT_DATASET_NAMES
                channel = OWENS.output_data_channel(name)
                attrs = HDF5.attrs(file[name])
                @test attrs["owens_channel_name"] == name
                @test attrs["units"] == channel.units
                @test attrs["dimensions"] == join(channel.dimensions, ",")
                @test attrs["frame"] == channel.frame
                @test attrs["association"] == channel.association
                @test attrs["source"] == channel.source
                @test attrs["sign_convention"] == channel.sign_convention
                @test attrs["description"] == channel.description
            end
        end
    end
end

@testset "outputData summary" begin
    mktempdir() do dir
        data_output = joinpath(dir, "unit.out")
        inputs = (; dataOutputFilename = data_output)

        vector = [1.0, 2.0, 3.0]
        dof_history = reshape(collect(1.0:18.0), 3, 6)
        strain_history = reshape(collect(1.0:24.0), 4, 2, 3)
        stress_history = reshape(collect(1.0:72.0), 3, 2, 4, 3)
        safety_factor = reshape(collect(1.0:24.0), 3, 2, 4)
        topstrain = reshape(collect(1.0:216.0), 3, 2, 4, 9)
        damage = reshape(collect(1.0:8.0), 2, 4)

        OWENS.outputData(;
            inputs,
            t = vector,
            aziHist = vector .+ 0.1,
            OmegaHist = vector .+ 0.2,
            OmegaDotHist = vector .+ 0.3,
            gbHist = vector .+ 0.4,
            gbDotHist = vector .+ 0.5,
            gbDotDotHist = vector .+ 0.6,
            FReactionHist = dof_history,
            FTwrBsHist = dof_history .+ 10.0,
            genTorque = vector .+ 0.7,
            genPower = vector .+ 0.8,
            torqueDriveShaft = vector .+ 0.9,
            uHist = dof_history .+ 20.0,
            uHist_prp = dof_history .+ 30.0,
            epsilon_x_hist = strain_history .+ 1.0,
            epsilon_y_hist = strain_history .+ 2.0,
            epsilon_z_hist = strain_history .+ 3.0,
            kappa_x_hist = strain_history .+ 4.0,
            kappa_y_hist = strain_history .+ 5.0,
            kappa_z_hist = strain_history .+ 6.0,
            massOwens = 12.5,
            stress_U = stress_history .+ 1.0,
            SF_ult_U = safety_factor .+ 1.0,
            SF_buck_U = safety_factor .+ 2.0,
            stress_L = stress_history .+ 2.0,
            SF_ult_L = safety_factor .+ 3.0,
            SF_buck_L = safety_factor .+ 4.0,
            stress_TU = stress_history .+ 3.0,
            SF_ult_TU = safety_factor .+ 5.0,
            SF_buck_TU = safety_factor .+ 6.0,
            stress_TL = stress_history .+ 4.0,
            SF_ult_TL = safety_factor .+ 7.0,
            SF_buck_TL = safety_factor .+ 8.0,
            topstrainout_blade_U = topstrain .+ 1.0,
            topstrainout_blade_L = topstrain .+ 2.0,
            topstrainout_tower_U = topstrain .+ 3.0,
            topstrainout_tower_L = topstrain .+ 4.0,
            topDamage_blade_U = damage .+ 1.0,
            topDamage_blade_L = damage .+ 2.0,
            topDamage_tower_U = damage .+ 3.0,
            topDamage_tower_L = damage .+ 4.0,
        )

        h5_file = joinpath(dir, "unit.h5")
        summary = OWENS.output_data_summary(h5_file)
        @test length(summary) == 41
        @test [row.name for row in summary] == OUTPUT_DATASET_NAMES
        @test all(row.present for row in summary)
        @test all(row.registered for row in summary)
        @test all(row.has_channel_attrs for row in summary)
        @test all(isempty(row.attr_mismatches) for row in summary)
        @test summary[1].name == "t"
        @test summary[1].shape == [3]
        @test summary[1].ndims == 1
        @test summary[1].eltype == "Float64"
        @test summary[1].units == "s"
        @test summary[1].dimensions == ["time"]
        @test summary[8].name == "FReactionHist"
        @test summary[8].shape == [3, 6]
        @test summary[8].units == "N or N*m by DOF"
        @test summary[21].name == "massOwens"
        @test summary[21].shape == Int[]
        @test summary[21].ndims == 0
        @test summary[22].name == "stress_U"
        @test summary[22].shape == [3, 2, 4, 3]

        partial_file = joinpath(dir, "partial.h5")
        HDF5.h5open(partial_file, "w") do file
            HDF5.write(file, "t", vector)
            HDF5.write(file, "custom_channel", [4.0, 5.0])
        end

        partial_summary = OWENS.output_data_summary(partial_file)
        @test length(partial_summary) == 41
        @test partial_summary[1].name == "t"
        @test partial_summary[1].present === true
        @test partial_summary[1].shape == [3]
        @test partial_summary[1].has_channel_attrs === false
        @test partial_summary[1].attr_mismatches == [
            "missing:owens_channel_name",
            "missing:units",
            "missing:dimensions",
            "missing:frame",
            "missing:association",
            "missing:source",
            "missing:sign_convention",
            "missing:description",
        ]
        @test partial_summary[2].name == "aziHist"
        @test partial_summary[2].present === false
        @test ismissing(partial_summary[2].shape)
        @test ismissing(partial_summary[2].ndims)
        @test ismissing(partial_summary[2].eltype)
        @test partial_summary[2].units == "rad"
        @test partial_summary[2].dimensions == ["time"]

        with_extra = OWENS.output_data_summary(partial_file; include_unregistered = true)
        @test length(with_extra) == 42
        @test with_extra[end].name == "custom_channel"
        @test with_extra[end].present === true
        @test with_extra[end].registered === false
        @test with_extra[end].shape == [2]
        @test with_extra[end].eltype == "Float64"
        @test ismissing(with_extra[end].units)
        @test ismissing(with_extra[end].dimensions)

        t_only = OWENS.output_data_summary(partial_file; channels = ["t"])
        @test length(t_only) == 1
        @test t_only[1].name == "t"
        @test t_only[1].present === true
        @test OWENS.output_data_summary(partial_file; channels = "t") == t_only

        @test_throws KeyError OWENS.output_data_summary(
            partial_file;
            channels = ["not_a_channel"],
        )
        @test_throws ArgumentError OWENS.output_data_summary(joinpath(dir, "missing.h5"))
    end
end

@testset "outputData channel metrics" begin
    mktempdir() do dir
        reference_file = joinpath(dir, "reference.h5")
        candidate_file = joinpath(dir, "candidate.h5")

        reference_reactions = reshape(collect(1.0:6.0), 2, 3)
        HDF5.h5open(reference_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0, 3.0])
            HDF5.write(file, "FReactionHist", reference_reactions)
            HDF5.write(file, "OmegaHist", [0.0, 0.0])
            HDF5.write(file, "OmegaDotHist", [0.0, 0.0])
            HDF5.write(file, "massOwens", 12.5)
        end
        HDF5.h5open(candidate_file, "w") do file
            HDF5.write(file, "t", [2.0, 1.0, 5.0])
            HDF5.write(file, "FReactionHist", reference_reactions .+ 0.5)
            HDF5.write(file, "OmegaHist", [0.0, 0.0])
            HDF5.write(file, "OmegaDotHist", [1.0, 0.0])
            HDF5.write(file, "genPower", [10.0, 20.0, 30.0])
            HDF5.write(file, "massOwens", 13.0)
        end

        metrics = OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["t", "FReactionHist", "massOwens"],
        )
        @test length(metrics) == 3
        @test [row.name for row in metrics] == ["t", "FReactionHist", "massOwens"]
        @test keys(metrics[1]) == (
            :name,
            :comparable,
            :passed,
            :status,
            :units,
            :dimensions,
            :reference_shape,
            :candidate_shape,
            :n,
            :bias,
            :rmse,
            :max_abs_error,
            :mean_abs_error,
            :reference_rms,
            :max_abs_reference,
            :relative_rmse,
            :tolerance,
        )

        time_metric = metrics[1]
        @test time_metric.comparable === true
        @test time_metric.passed === false
        @test time_metric.status == "fail"
        @test time_metric.units == "s"
        @test time_metric.dimensions == ["time"]
        @test time_metric.reference_shape == [3]
        @test time_metric.candidate_shape == [3]
        @test time_metric.n == 3
        @test time_metric.bias == 2.0 / 3.0
        @test time_metric.rmse == sqrt(2.0)
        @test time_metric.max_abs_error == 2.0
        @test time_metric.mean_abs_error == 4.0 / 3.0
        @test time_metric.reference_rms == sqrt(14.0 / 3.0)
        @test time_metric.max_abs_reference == 3.0
        @test time_metric.relative_rmse == sqrt(2.0) / sqrt(14.0 / 3.0)
        @test time_metric.tolerance == 1.0e-6 * 3.0

        reaction_metric = metrics[2]
        @test reaction_metric.comparable === true
        @test reaction_metric.passed === false
        @test reaction_metric.status == "fail"
        @test reaction_metric.units == "N or N*m by DOF"
        @test reaction_metric.dimensions == ["time", "structural_dof"]
        @test reaction_metric.reference_shape == [2, 3]
        @test reaction_metric.candidate_shape == [2, 3]
        @test reaction_metric.n == 6
        @test reaction_metric.bias == 0.5
        @test reaction_metric.rmse == 0.5
        @test reaction_metric.max_abs_error == 0.5
        @test reaction_metric.mean_abs_error == 0.5
        @test reaction_metric.reference_rms == sqrt(91.0 / 6.0)
        @test reaction_metric.max_abs_reference == 6.0
        @test reaction_metric.relative_rmse == 0.5 / sqrt(91.0 / 6.0)
        @test reaction_metric.tolerance == 1.0e-6 * 6.0

        scalar_metric = metrics[3]
        @test scalar_metric.comparable === true
        @test scalar_metric.passed === false
        @test scalar_metric.status == "fail"
        @test scalar_metric.reference_shape == Int[]
        @test scalar_metric.candidate_shape == Int[]
        @test scalar_metric.n == 1
        @test scalar_metric.bias == 0.5
        @test scalar_metric.rmse == 0.5
        @test scalar_metric.max_abs_error == 0.5
        @test scalar_metric.mean_abs_error == 0.5
        @test scalar_metric.reference_rms == 12.5
        @test scalar_metric.max_abs_reference == 12.5
        @test scalar_metric.relative_rmse == 0.5 / 12.5
        @test scalar_metric.tolerance == 1.0e-6 * 12.5

        relaxed_metrics = OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["t", "FReactionHist", "massOwens"],
            atol = 2.0,
        )
        @test [row.passed for row in relaxed_metrics] == [true, true, true]
        @test [row.status for row in relaxed_metrics] == ["pass", "pass", "pass"]
        @test [row.tolerance for row in relaxed_metrics] == [2.0 + 1.0e-6 * 3.0, 2.0 + 1.0e-6 * 6.0, 2.0 + 1.0e-6 * 12.5]

        string_channel_metric = OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = "t",
            atol = 2.0,
        )
        @test length(string_channel_metric) == 1
        @test string_channel_metric[1] == relaxed_metrics[1]

        zero_scale_metrics = OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["OmegaHist", "OmegaDotHist"],
        )
        @test zero_scale_metrics[1].status == "pass"
        @test zero_scale_metrics[1].relative_rmse == 0.0
        @test zero_scale_metrics[1].tolerance == 0.0
        @test zero_scale_metrics[2].status == "fail"
        @test zero_scale_metrics[2].relative_rmse == Inf
        @test zero_scale_metrics[2].max_abs_error == 1.0
        @test zero_scale_metrics[2].tolerance == 0.0

        missing_reference_file = joinpath(dir, "missing_reference.h5")
        HDF5.h5open(missing_reference_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0, 3.0])
        end
        missing_reference_metric = OWENS.output_data_channel_metrics(
            missing_reference_file,
            candidate_file;
            channels = ["genPower"],
        )
        @test missing_reference_metric[1].comparable === false
        @test missing_reference_metric[1].passed === false
        @test missing_reference_metric[1].status == "missing_reference"
        @test ismissing(missing_reference_metric[1].reference_shape)
        @test missing_reference_metric[1].candidate_shape == [3]
        @test missing_reference_metric[1].n == 0
        @test ismissing(missing_reference_metric[1].rmse)
        @test ismissing(missing_reference_metric[1].tolerance)

        missing_candidate_metric = OWENS.output_data_channel_metrics(
            reference_file,
            missing_reference_file;
            channels = ["massOwens"],
        )
        @test missing_candidate_metric[1].comparable === false
        @test missing_candidate_metric[1].passed === false
        @test missing_candidate_metric[1].status == "missing_candidate"
        @test missing_candidate_metric[1].reference_shape == Int[]
        @test ismissing(missing_candidate_metric[1].candidate_shape)

        missing_both_metric = OWENS.output_data_channel_metrics(
            missing_reference_file,
            missing_reference_file;
            channels = ["genTorque"],
        )
        @test missing_both_metric[1].comparable === false
        @test missing_both_metric[1].passed === false
        @test missing_both_metric[1].status == "missing_reference_and_candidate"
        @test ismissing(missing_both_metric[1].reference_shape)
        @test ismissing(missing_both_metric[1].candidate_shape)

        shape_mismatch_file = joinpath(dir, "shape_mismatch.h5")
        HDF5.h5open(shape_mismatch_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0])
        end
        shape_mismatch_metric = OWENS.output_data_channel_metrics(
            reference_file,
            shape_mismatch_file;
            channels = ["t"],
        )
        @test shape_mismatch_metric[1].comparable === false
        @test shape_mismatch_metric[1].passed === false
        @test shape_mismatch_metric[1].status == "shape_mismatch"
        @test shape_mismatch_metric[1].reference_shape == [3]
        @test shape_mismatch_metric[1].candidate_shape == [2]
        @test ismissing(shape_mismatch_metric[1].max_abs_error)

        nonnumeric_file = joinpath(dir, "nonnumeric.h5")
        HDF5.h5open(nonnumeric_file, "w") do file
            HDF5.write(file, "genTorque", "not numeric")
        end
        nonnumeric_metric = OWENS.output_data_channel_metrics(
            nonnumeric_file,
            nonnumeric_file;
            channels = ["genTorque"],
        )
        @test nonnumeric_metric[1].comparable === false
        @test nonnumeric_metric[1].passed === false
        @test nonnumeric_metric[1].status == "non_numeric"
        @test nonnumeric_metric[1].reference_shape == Int[]
        @test nonnumeric_metric[1].candidate_shape == Int[]

        nonfinite_file = joinpath(dir, "nonfinite.h5")
        HDF5.h5open(nonfinite_file, "w") do file
            HDF5.write(file, "t", [1.0, NaN])
        end
        nonfinite_metric = OWENS.output_data_channel_metrics(
            nonfinite_file,
            nonfinite_file;
            channels = ["t"],
        )
        @test nonfinite_metric[1].comparable === false
        @test nonfinite_metric[1].passed === false
        @test nonfinite_metric[1].status == "non_finite"
        @test nonfinite_metric[1].reference_shape == [2]
        @test nonfinite_metric[1].candidate_shape == [2]

        empty_file = joinpath(dir, "empty.h5")
        HDF5.h5open(empty_file, "w") do file
            HDF5.write(file, "t", Float64[])
        end
        empty_metric =
            OWENS.output_data_channel_metrics(empty_file, empty_file; channels = ["t"])
        @test empty_metric[1].comparable === false
        @test empty_metric[1].passed === false
        @test empty_metric[1].status == "empty"
        @test empty_metric[1].reference_shape == [0]
        @test empty_metric[1].candidate_shape == [0]

        @test_throws KeyError OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["not_a_channel"],
        )
        @test_throws ArgumentError OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["t"],
            atol = -1.0,
        )
        @test_throws ArgumentError OWENS.output_data_channel_metrics(
            reference_file,
            candidate_file;
            channels = ["t"],
            rtol = -1.0,
        )
        @test_throws ArgumentError OWENS.output_data_channel_metrics(
            joinpath(dir, "does_not_exist_reference.h5"),
            candidate_file,
        )
        @test_throws ArgumentError OWENS.output_data_channel_metrics(
            reference_file,
            joinpath(dir, "does_not_exist_candidate.h5"),
        )
    end
end
