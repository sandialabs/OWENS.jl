using Test
import HDF5
import OWENS
import OrderedCollections

@testset "File provenance health rows" begin
    mktempdir() do dir
        input_file = joinpath(dir, "input.txt")
        write(input_file, "alpha\n")
        record = OWENS.file_provenance(input_file; root = dir, role = "input")

        ok = OWENS.verify_file_provenance(record; root = dir)
        @test collect(keys(ok)) == [
            "status",
            "issues",
            "path",
            "resolved_path",
            "role",
            "expected_bytes",
            "actual_bytes",
            "expected_sha256",
            "actual_sha256",
        ]
        @test ok["status"] == "ok"
        @test ok["issues"] == String[]
        @test ok["path"] == "input.txt"
        @test ok["resolved_path"] == abspath(input_file)
        @test ok["role"] == "input"
        @test ok["expected_bytes"] == 6
        @test ok["actual_bytes"] == 6
        @test ok["expected_sha256"] ==
              "b6a98d9ce9a2d9149288fa3df42d377c3e42737afdcdaf714e33c0a100b51060"
        @test ok["actual_sha256"] ==
              "b6a98d9ce9a2d9149288fa3df42d377c3e42737afdcdaf714e33c0a100b51060"

        write(input_file, "changed\n")
        modified = OWENS.verify_file_provenance(record; root = dir)
        @test modified["status"] == "modified"
        @test modified["issues"] == ["bytes mismatch", "sha256 mismatch"]
        @test modified["expected_bytes"] == 6
        @test modified["actual_bytes"] == 8
        @test modified["expected_sha256"] ==
              "b6a98d9ce9a2d9149288fa3df42d377c3e42737afdcdaf714e33c0a100b51060"
        @test modified["actual_sha256"] ==
              "7f8b1dfc466b6249f06cbe55c9174df2578e7754da793fded244ef5cba2a38f1"

        rm(input_file)
        missing = OWENS.verify_file_provenance(record; root = dir)
        @test missing["status"] == "missing"
        @test missing["issues"] == ["missing file"]
        @test missing["expected_bytes"] == 6
        @test isnothing(missing["actual_bytes"])
        @test isnothing(missing["actual_sha256"])

        invalid = OWENS.verify_file_provenance(
            Dict("path" => 42, "bytes" => true, "sha256" => "bad", "role" => 7);
            root = dir,
        )
        @test invalid["status"] == "invalid_record"
        @test invalid["issues"] == [
            "path must be a string",
            "bytes must be a non-negative integer",
            "sha256 must be a lowercase 64-character SHA-256 digest",
            "role must be a string when present",
        ]
        @test isnothing(invalid["path"])
        @test isnothing(invalid["resolved_path"])
        @test invalid["expected_bytes"] === true
        @test invalid["expected_sha256"] == "bad"

        missing_fields = OWENS.verify_file_provenance(Dict{String,Any}(); root = dir)
        @test missing_fields["status"] == "invalid_record"
        @test missing_fields["issues"] ==
              ["path is required", "bytes is required", "sha256 is required"]

        non_dictionary = OWENS.verify_file_provenance("not a manifest record"; root = dir)
        @test non_dictionary["status"] == "invalid_record"
        @test non_dictionary["issues"] == ["record must be a dictionary"]
    end
end

@testset "Run manifest health summary" begin
    mktempdir() do dir
        input_file = joinpath(dir, "input.yml")
        output_file = joinpath(dir, "output.h5")
        generated_file = joinpath(dir, "generated.vtp")
        manifest_file = joinpath(dir, "run_manifest.yml")

        write(input_file, "input: unit\n")
        HDF5.h5open(output_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0, 3.0])
        end
        write(generated_file, "generated: unit\n")

        manifest = OWENS.build_run_manifest(;
            run_id = "health-unit",
            run_name = "health unit",
            project_root = dir,
            solver = "unit-solver",
            input_files = [input_file],
            output_files = [output_file],
            generated_files = [generated_file],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        OWENS.write_run_manifest(manifest_file, manifest)

        health = OWENS.run_manifest_health(manifest)
        @test collect(keys(health)) == [
            "schema_version",
            "status",
            "manifest_path",
            "root",
            "manifest_issues",
            "summary",
            "inputs",
            "outputs",
            "generated",
        ]
        @test health["schema_version"] == "owens-run-health/v1"
        @test health["status"] == "ok"
        @test isnothing(health["manifest_path"])
        @test health["root"] == abspath(dir)
        @test health["manifest_issues"] == String[]
        @test health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 3,
            "modified" => 0,
            "missing" => 0,
            "invalid_record" => 0,
        )
        @test length(health["inputs"]) == 1
        @test length(health["outputs"]) == 1
        @test length(health["generated"]) == 1
        @test health["inputs"][1]["status"] == "ok"
        @test health["outputs"][1]["status"] == "ok"
        @test health["generated"][1]["status"] == "ok"
        @test haskey(health["outputs"][1], "output_data_summary")
        @test length(health["outputs"][1]["output_data_summary"]) == 41
        @test health["outputs"][1]["output_data_summary"][1]["name"] == "t"
        @test health["outputs"][1]["output_data_summary"][1]["present"] === true
        @test health["outputs"][1]["output_data_summary"][1]["shape"] == [3]
        @test health["outputs"][1]["output_data_summary"][1]["attr_mismatches"] == [
            "missing:owens_channel_name",
            "missing:units",
            "missing:dimensions",
            "missing:frame",
            "missing:association",
            "missing:source",
            "missing:sign_convention",
            "missing:description",
        ]

        loaded_health = OWENS.run_manifest_health(manifest_file)
        @test loaded_health["status"] == "ok"
        @test loaded_health["manifest_path"] == abspath(manifest_file)
        @test loaded_health["summary"] == health["summary"]

        no_project_root_manifest = deepcopy(manifest)
        delete!(no_project_root_manifest, "project_root")
        path_root_health = OWENS.run_manifest_health(
            no_project_root_manifest;
            manifest_path = manifest_file,
        )
        @test path_root_health["root"] == abspath(dir)

        cwd_root_health = cd(dir) do
            OWENS.run_manifest_health(no_project_root_manifest)
        end
        cwd_root_expected = cd(dir) do
            OWENS._canonical_abs_path(pwd())
        end
        @test cwd_root_health["root"] == cwd_root_expected

        no_output_summary = OWENS.run_manifest_health(manifest; summarize_outputs = false)
        @test no_output_summary["outputs"][1]["status"] == "ok"
        @test !haskey(no_output_summary["outputs"][1], "output_data_summary")

        stale_root_manifest = deepcopy(manifest)
        stale_root_manifest["project_root"] = joinpath(dir, "stale")
        stale_health = OWENS.run_manifest_health(stale_root_manifest)
        @test stale_health["status"] == "attention"
        @test stale_health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 0,
            "modified" => 0,
            "missing" => 3,
            "invalid_record" => 0,
        )
        override_health = OWENS.run_manifest_health(stale_root_manifest; root = dir)
        @test override_health["status"] == "ok"
        @test override_health["root"] == abspath(dir)
        @test override_health["summary"] == health["summary"]

        write(input_file, "input: changed\n")
        rm(generated_file)
        drift_health = OWENS.run_manifest_health(manifest)
        @test drift_health["status"] == "attention"
        @test drift_health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 1,
            "modified" => 1,
            "missing" => 1,
            "invalid_record" => 0,
        )
        @test drift_health["inputs"][1]["status"] == "modified"
        @test drift_health["inputs"][1]["issues"] == ["bytes mismatch", "sha256 mismatch"]
        @test drift_health["outputs"][1]["status"] == "ok"
        @test drift_health["generated"][1]["status"] == "missing"

        invalid_manifest = deepcopy(manifest)
        push!(
            invalid_manifest["inputs"],
            Dict("path" => 42, "bytes" => true, "sha256" => "bad"),
        )
        invalid_health = OWENS.run_manifest_health(invalid_manifest)
        @test invalid_health["status"] == "attention"
        @test invalid_health["summary"]["invalid_record"] == 1
        @test invalid_health["inputs"][2]["status"] == "invalid_record"
        @test invalid_health["manifest_issues"][1] == "inputs[2].path must be a string"
    end
end
