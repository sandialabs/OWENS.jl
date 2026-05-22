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
    @test stress_channel.dimensions == ["time", "span_station", "chord_station", "stress_component"]
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
            t=vector,
            aziHist=vector .+ 0.1,
            OmegaHist=vector .+ 0.2,
            OmegaDotHist=vector .+ 0.3,
            gbHist=vector .+ 0.4,
            gbDotHist=vector .+ 0.5,
            gbDotDotHist=vector .+ 0.6,
            FReactionHist=dof_history,
            FTwrBsHist=dof_history .+ 10.0,
            genTorque=vector .+ 0.7,
            genPower=vector .+ 0.8,
            torqueDriveShaft=vector .+ 0.9,
            uHist=dof_history .+ 20.0,
            uHist_prp=dof_history .+ 30.0,
            epsilon_x_hist=strain_history .+ 1.0,
            epsilon_y_hist=strain_history .+ 2.0,
            epsilon_z_hist=strain_history .+ 3.0,
            kappa_x_hist=strain_history .+ 4.0,
            kappa_y_hist=strain_history .+ 5.0,
            kappa_z_hist=strain_history .+ 6.0,
            massOwens=12.5,
            stress_U=stress_history .+ 1.0,
            SF_ult_U=safety_factor .+ 1.0,
            SF_buck_U=safety_factor .+ 2.0,
            stress_L=stress_history .+ 2.0,
            SF_ult_L=safety_factor .+ 3.0,
            SF_buck_L=safety_factor .+ 4.0,
            stress_TU=stress_history .+ 3.0,
            SF_ult_TU=safety_factor .+ 5.0,
            SF_buck_TU=safety_factor .+ 6.0,
            stress_TL=stress_history .+ 4.0,
            SF_ult_TL=safety_factor .+ 7.0,
            SF_buck_TL=safety_factor .+ 8.0,
            topstrainout_blade_U=topstrain .+ 1.0,
            topstrainout_blade_L=topstrain .+ 2.0,
            topstrainout_tower_U=topstrain .+ 3.0,
            topstrainout_tower_L=topstrain .+ 4.0,
            topDamage_blade_U=damage .+ 1.0,
            topDamage_blade_L=damage .+ 2.0,
            topDamage_tower_U=damage .+ 3.0,
            topDamage_tower_L=damage .+ 4.0,
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
