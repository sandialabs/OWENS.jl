using Test
import HDF5
import OWENS
import OrderedCollections

@testset "Output data validation report" begin
    mktempdir() do dir
        reference_file = joinpath(dir, "reference.h5")
        candidate_file = joinpath(dir, "candidate.h5")
        report_file = joinpath(dir, "reports", "validation.yml")

        HDF5.h5open(reference_file, "w") do file
            HDF5.write(file, "t", [0.0, 1.0, 2.0])
            HDF5.write(file, "genPower", [10.0, 20.0, 30.0])
            HDF5.write(file, "OmegaHist", [0.0, 0.0])
            HDF5.write(file, "OmegaDotHist", [1.0, NaN])
        end
        HDF5.h5open(candidate_file, "w") do file
            HDF5.write(file, "t", [0.0, 1.05, 1.95])
            HDF5.write(file, "genPower", [10.0, 25.0, 30.0])
            HDF5.write(file, "OmegaHist", [0.0, 0.0])
            HDF5.write(file, "OmegaDotHist", [1.0, 2.0])
        end

        tolerances = OrderedCollections.OrderedDict{String,Any}(
            "t" => (; atol = 0.1),
            "genPower" => Dict(:rtol => 0.1),
        )
        report = OWENS.build_output_data_validation_report(
            reference_file,
            candidate_file;
            channels = ["t", "genPower", "OmegaHist", "OmegaDotHist", "genTorque"],
            tolerances,
            case_id = "validation-report-unit",
            metadata = Dict(:source => :unit, :active => true),
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        @test collect(keys(report)) == [
            "schema_version",
            "case_id",
            "created_at_utc",
            "status",
            "comparison",
            "reference",
            "candidate",
            "summary",
            "metadata",
            "channels",
        ]
        @test report["schema_version"] == "owens-output-data-validation/v1"
        @test report["case_id"] == "validation-report-unit"
        @test report["created_at_utc"] == "2026-05-20T00:00:00.000Z"
        @test report["status"] == "failed"
        @test report["comparison"]["kind"] == "outputData_channel_metrics"
        @test report["comparison"]["default_atol"] == 0.0
        @test report["comparison"]["default_rtol"] == 1.0e-6
        @test report["comparison"]["note"] ==
              "Elementwise comparison of selected registered outputData channels; no time alignment, interpolation, phase correction, or whole-revolution windowing is applied."
        @test report["reference"]["path"] == abspath(reference_file)
        @test report["reference"]["bytes"] == stat(reference_file).size
        @test report["reference"]["sha256"] == OWENS.file_sha256(reference_file)
        @test report["reference"]["role"] == "reference_output_data"
        @test report["candidate"]["path"] == abspath(candidate_file)
        @test report["candidate"]["bytes"] == stat(candidate_file).size
        @test report["candidate"]["sha256"] == OWENS.file_sha256(candidate_file)
        @test report["candidate"]["role"] == "candidate_output_data"
        @test report["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "channels_requested" => 5,
            "channels_comparable" => 3,
            "channels_noncomparable" => 2,
            "channels_passed" => 2,
            "channels_failed" => 3,
        )
        @test collect(keys(report["metadata"])) == ["active", "source"]
        @test report["metadata"]["active"] === true
        @test report["metadata"]["source"] == "unit"

        @test length(report["channels"]) == 5
        @test [row["name"] for row in report["channels"]] == ["t", "genPower", "OmegaHist", "OmegaDotHist", "genTorque"]
        @test collect(keys(report["channels"][1])) == [
            "name",
            "passed",
            "comparable",
            "status",
            "units",
            "dimensions",
            "reference_shape",
            "candidate_shape",
            "n",
            "bias",
            "rmse",
            "max_abs_error",
            "mean_abs_error",
            "reference_rms",
            "max_abs_reference",
            "relative_rmse",
            "tolerance",
            "atol",
            "rtol",
        ]

        time_row = report["channels"][1]
        @test time_row["passed"] === true
        @test time_row["comparable"] === true
        @test time_row["status"] == "pass"
        @test time_row["units"] == "s"
        @test time_row["dimensions"] == ["time"]
        @test time_row["reference_shape"] == [3]
        @test time_row["candidate_shape"] == [3]
        @test time_row["n"] == 3
        @test time_row["bias"] == 0.0
        time_difference = [0.0, 1.05 - 1.0, 1.95 - 2.0]
        @test time_row["rmse"] == sqrt(sum(abs2, time_difference) / 3.0)
        @test time_row["max_abs_error"] == maximum(abs.(time_difference))
        @test time_row["mean_abs_error"] == sum(abs.(time_difference)) / 3.0
        @test time_row["reference_rms"] == sqrt(5.0 / 3.0)
        @test time_row["max_abs_reference"] == 2.0
        @test time_row["relative_rmse"] ==
              sqrt(sum(abs2, time_difference) / 3.0) / sqrt(5.0 / 3.0)
        @test time_row["tolerance"] == 0.1 + 1.0e-6 * 2.0
        @test time_row["atol"] == 0.1
        @test time_row["rtol"] == 1.0e-6

        power_row = report["channels"][2]
        @test power_row["passed"] === false
        @test power_row["comparable"] === true
        @test power_row["status"] == "fail"
        @test power_row["n"] == 3
        @test power_row["bias"] == 5.0 / 3.0
        @test power_row["rmse"] == sqrt(25.0 / 3.0)
        @test power_row["max_abs_error"] == 5.0
        @test power_row["mean_abs_error"] == 5.0 / 3.0
        @test power_row["reference_rms"] == sqrt(1400.0 / 3.0)
        @test power_row["max_abs_reference"] == 30.0
        @test power_row["relative_rmse"] == sqrt(25.0 / 3.0) / sqrt(1400.0 / 3.0)
        @test power_row["tolerance"] == 3.0
        @test power_row["atol"] == 0.0
        @test power_row["rtol"] == 0.1

        omega_row = report["channels"][3]
        @test omega_row["passed"] === true
        @test omega_row["status"] == "pass"
        @test omega_row["relative_rmse"] == 0.0
        @test omega_row["tolerance"] == 0.0

        nonfinite_row = report["channels"][4]
        @test nonfinite_row["passed"] === false
        @test nonfinite_row["comparable"] === false
        @test nonfinite_row["status"] == "non_finite"
        @test nonfinite_row["reference_shape"] == [2]
        @test nonfinite_row["candidate_shape"] == [2]
        @test isnothing(nonfinite_row["rmse"])
        @test isnothing(nonfinite_row["tolerance"])
        @test nonfinite_row["atol"] == 0.0
        @test nonfinite_row["rtol"] == 1.0e-6

        missing_row = report["channels"][5]
        @test missing_row["passed"] === false
        @test missing_row["comparable"] === false
        @test missing_row["status"] == "missing_reference_and_candidate"
        @test isnothing(missing_row["reference_shape"])
        @test isnothing(missing_row["candidate_shape"])
        @test missing_row["n"] == 0
        @test isnothing(missing_row["max_abs_error"])

        written_report = OWENS.write_output_data_validation_report(
            report_file,
            reference_file,
            candidate_file;
            channels = ["t", "genPower", "OmegaHist", "OmegaDotHist", "genTorque"],
            tolerances,
            case_id = "validation-report-unit",
            metadata = Dict(:source => :unit, :active => true),
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        loaded_report = OWENS.read_output_data_validation_report(report_file)
        @test isfile(report_file)
        @test written_report["schema_version"] == "owens-output-data-validation/v1"
        @test loaded_report["schema_version"] == "owens-output-data-validation/v1"
        @test loaded_report["case_id"] == "validation-report-unit"
        @test loaded_report["status"] == "failed"
        @test loaded_report["summary"]["channels_requested"] == 5
        @test loaded_report["channels"][1]["status"] == "pass"
        @test loaded_report["channels"][2]["status"] == "fail"
        @test loaded_report["channels"][4]["status"] == "non_finite"

        @test_throws KeyError OWENS.build_output_data_validation_report(
            reference_file,
            candidate_file;
            channels = ["not_a_channel"],
        )
        @test_throws ArgumentError OWENS.build_output_data_validation_report(
            reference_file,
            candidate_file;
            channels = ["t"],
            tolerances = ["not" => "a dict"],
        )
        @test_throws ArgumentError OWENS.build_output_data_validation_report(
            reference_file,
            candidate_file;
            channels = ["t"],
            tolerances = Dict("t" => (; atol = -1.0)),
        )
        @test_throws ArgumentError OWENS.read_output_data_validation_report(
            joinpath(dir, "missing.yml"),
        )
    end
end
