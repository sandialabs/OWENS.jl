using Test
import OWENS
import Dates
import OrderedCollections

const ALPHA_SHA256 = "b6a98d9ce9a2d9149288fa3df42d377c3e42737afdcdaf714e33c0a100b51060"
const GAMMA_SHA256 = "ae9a6306a205417afddd14316cc1d0d5e04a98f1be10865dce643925ee070ce2"
const PROJECT_SHA256 = "a9a9ed3918e077096e2af7e826fef92cd55faf11207764efb6a77a4d09ecb890"
const INPUT_SHA256 = "ddb019638b1605fc8a6bcd8ea9a8e578ea56fd91a44ff0c5ac7f450f1bde51d4"
const OUTPUT_SHA256 = "d1bf324bfe2e9bfa7067ff1f27010f6d98e3bb3dd0104864f8632e11678bfa36"
const GENERATED_SHA256 = "90e6df81ce3177ae762f91b5d27efa667cc2efcdc9c7e7310a45fb5e581255a7"

@testset "File provenance records" begin
    mktempdir() do dir
        input_file = joinpath(dir, "input.txt")
        write(input_file, "alpha\n")

        record = OWENS.file_provenance(input_file; root = dir, role = :airfoil)
        @test record isa OrderedCollections.OrderedDict{String,Any}
        @test collect(keys(record)) == ["path", "bytes", "sha256", "role"]
        @test record["path"] == "input.txt"
        @test record["bytes"] == 6
        @test record["sha256"] == ALPHA_SHA256
        @test record["role"] == "airfoil"
        @test OWENS.file_sha256(input_file) == ALPHA_SHA256

        write(input_file, "gamma\n")
        @test OWENS.file_sha256(input_file) == GAMMA_SHA256

        missing_file = joinpath(dir, "missing.txt")
        @test_throws ArgumentError OWENS.file_sha256(missing_file)
        @test_throws ArgumentError OWENS.file_provenance(missing_file; root = dir)
    end
end

@testset "Run manifest YAML round trip" begin
    mktempdir() do dir
        project_file = joinpath(dir, "project.yml")
        input_file = joinpath(dir, "input.yml")
        output_file = joinpath(dir, "output.h5")
        generated_file = joinpath(dir, "generated.vtp")
        manifest_file = joinpath(dir, "runs", "run_manifest.yml")

        write(project_file, "project: unit\n")
        write(input_file, "input: unit\n")
        write(output_file, "output: unit\n")
        write(generated_file, "generated: unit\n")

        parameters = OrderedCollections.OrderedDict{String,Any}(
            "rpm" => 12.5,
            "flags" => [:tip_loss, :dynamic_stall],
            "nested" => Dict(:b => 2, :a => 1),
        )
        metadata = Dict(:case => "unit", :active => true)

        manifest = OWENS.build_run_manifest(;
            run_id = "run-unit-001",
            run_name = "manifest unit test",
            project_root = dir,
            solver = "unit-solver",
            project_file = project_file,
            input_files = [input_file],
            output_files = [output_file],
            generated_files = [generated_file],
            parameters,
            metadata,
            warnings = ["low residual warning"],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        @test manifest isa OrderedCollections.OrderedDict{String,Any}
        @test collect(keys(manifest)) == [
            "schema_version",
            "run_id",
            "run_name",
            "created_at_utc",
            "project_root",
            "solver",
            "status",
            "julia",
            "packages",
            "git",
            "parameters",
            "metadata",
            "inputs",
            "outputs",
            "generated",
            "warnings",
        ]
        @test manifest["schema_version"] == "owens-run-manifest/v1"
        @test manifest["run_id"] == "run-unit-001"
        @test manifest["run_name"] == "manifest unit test"
        @test manifest["created_at_utc"] == "2026-05-20T00:00:00.000Z"
        @test manifest["project_root"] == abspath(dir)
        @test manifest["solver"] == "unit-solver"
        @test manifest["status"] == "complete"
        @test manifest["julia"]["version"] == string(VERSION)
        @test manifest["julia"]["threads"] isa Int
        @test manifest["packages"]["OWENS"] isa Union{Nothing,String}
        @test manifest["git"]["available"] isa Bool
        @test manifest["git"]["root"] isa String
        @test manifest["warnings"] == ["low residual warning"]

        @test length(manifest["inputs"]) == 2
        @test manifest["inputs"][1]["path"] == "project.yml"
        @test manifest["inputs"][1]["bytes"] == 14
        @test manifest["inputs"][1]["sha256"] == PROJECT_SHA256
        @test manifest["inputs"][1]["role"] == "project"
        @test manifest["inputs"][2]["path"] == "input.yml"
        @test manifest["inputs"][2]["bytes"] == 12
        @test manifest["inputs"][2]["sha256"] == INPUT_SHA256
        @test manifest["inputs"][2]["role"] == "input"

        @test length(manifest["outputs"]) == 1
        @test manifest["outputs"][1]["path"] == "output.h5"
        @test manifest["outputs"][1]["bytes"] == 13
        @test manifest["outputs"][1]["sha256"] == OUTPUT_SHA256
        @test manifest["outputs"][1]["role"] == "output"

        @test length(manifest["generated"]) == 1
        @test manifest["generated"][1]["path"] == "generated.vtp"
        @test manifest["generated"][1]["bytes"] == 16
        @test manifest["generated"][1]["sha256"] == GENERATED_SHA256
        @test manifest["generated"][1]["role"] == "generated"

        @test manifest["parameters"]["rpm"] == 12.5
        @test manifest["parameters"]["flags"] == ["tip_loss", "dynamic_stall"]
        @test collect(keys(manifest["parameters"]["nested"])) == ["a", "b"]
        @test manifest["parameters"]["nested"]["a"] == 1
        @test manifest["parameters"]["nested"]["b"] == 2
        @test collect(keys(manifest["metadata"])) == ["active", "case"]
        @test manifest["metadata"]["active"] === true
        @test manifest["metadata"]["case"] == "unit"
        @test OWENS.run_manifest_issues(manifest) == String[]
        @test OWENS.validate_run_manifest(manifest) === manifest

        written = OWENS.write_run_manifest(manifest_file, manifest)
        loaded = OWENS.read_run_manifest(manifest_file)
        @test isfile(manifest_file)
        @test written["schema_version"] == "owens-run-manifest/v1"
        @test loaded["schema_version"] == "owens-run-manifest/v1"
        @test loaded["run_id"] == "run-unit-001"
        @test loaded["inputs"][1]["sha256"] == PROJECT_SHA256
        @test loaded["outputs"][1]["sha256"] == OUTPUT_SHA256
        @test loaded["generated"][1]["sha256"] == GENERATED_SHA256
        @test loaded["parameters"]["flags"] == ["tip_loss", "dynamic_stall"]
        @test OWENS.run_manifest_issues(manifest_file) == String[]
        @test OWENS.validate_run_manifest(manifest_file)["run_id"] == "run-unit-001"

        @test_throws ArgumentError OWENS.read_run_manifest(joinpath(dir, "missing.yml"))
    end
end

@testset "Run manifest defaults and value normalization" begin
    mktempdir() do dir
        input_file = joinpath(dir, "input.txt")
        manifest_file = joinpath(dir, "manifest.yml")
        write(input_file, "alpha\n")

        manifest = OWENS.write_run_manifest(
            manifest_file;
            run_name="defaults unit test",
            project_root=dir,
            input_files=[input_file],
            parameters=Dict(:release => v"1.2.3", :date => Dates.Date(2026, 5, 22)),
        )
        loaded = OWENS.read_run_manifest(manifest_file)

        @test startswith(manifest["run_id"], "run-")
        @test endswith(manifest["created_at_utc"], "Z")
        @test manifest["inputs"][1]["path"] == "input.txt"
        @test manifest["parameters"]["date"] == "2026-05-22"
        @test manifest["parameters"]["release"] == "1.2.3"
        @test loaded["run_id"] == manifest["run_id"]
        @test loaded["parameters"]["date"] == "2026-05-22"
    end

    @test OWENS._manifest_value(v"2.0.1") == "2.0.1"
    @test OWENS._manifest_value(Dates.Date(2026, 5, 22)) == "2026-05-22"
    @test OWENS._manifest_value((:alpha, v"1.0.0")) == ["alpha", "1.0.0"]
end

@testset "Run manifest validation diagnostics" begin
    mktempdir() do dir
        input_file = joinpath(dir, "input.yml")
        output_file = joinpath(dir, "output.h5")
        write(input_file, "input: unit\n")
        write(output_file, "output: unit\n")

        manifest = OWENS.build_run_manifest(;
            run_id = "validation-unit",
            run_name = "validation unit",
            project_root = dir,
            solver = "unit-solver",
            input_files = [input_file],
            output_files = [output_file],
            warnings = ["unit warning"],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        @test OWENS.run_manifest_issues(manifest) == String[]

        missing_schema = deepcopy(manifest)
        delete!(missing_schema, "schema_version")
        @test OWENS.run_manifest_issues(missing_schema) ==
              ["missing required key: schema_version"]

        old_schema = deepcopy(manifest)
        old_schema["schema_version"] = "owens-run-manifest/v0"
        @test OWENS.run_manifest_issues(old_schema) ==
              ["schema_version must equal owens-run-manifest/v1"]

        malformed = deepcopy(manifest)
        malformed["run_id"] = 42
        malformed["created_at_utc"] = false
        malformed["julia"] = "Julia 1.12"
        malformed["parameters"] = ["rpm"]
        malformed["inputs"] = [
            Dict(
                "path" => 42,
                "bytes" => true,
                "sha256" => uppercase(INPUT_SHA256),
                "role" => 7,
            ),
            "not a record",
        ]
        malformed["outputs"] =
            [Dict("path" => "output.h5", "bytes" => -1, "sha256" => "abc")]
        malformed["generated"] = [Dict("path" => "generated.vtp", "bytes" => 1)]
        malformed["warnings"] = [3]
        @test OWENS.run_manifest_issues(malformed) == [
            "run_id must be a string",
            "created_at_utc must be a string",
            "julia must be a dictionary",
            "parameters must be a dictionary",
            "inputs[1].path must be a string",
            "inputs[1].bytes must be a non-negative integer",
            "inputs[1].sha256 must be a lowercase 64-character SHA-256 digest",
            "inputs[1].role must be a string when present",
            "inputs[2] must be a dictionary",
            "outputs[1].bytes must be a non-negative integer",
            "outputs[1].sha256 must be a lowercase 64-character SHA-256 digest",
            "generated[1].sha256 is required",
            "warnings[1] must be a string",
        ]
        @test_throws ArgumentError OWENS.validate_run_manifest(malformed)

        missing_keys = OrderedCollections.OrderedDict{String,Any}()
        @test OWENS.run_manifest_issues(missing_keys) == [
            "missing required key: schema_version",
            "missing required key: run_id",
            "missing required key: run_name",
            "missing required key: created_at_utc",
            "missing required key: project_root",
            "missing required key: solver",
            "missing required key: status",
            "missing required key: julia",
            "missing required key: packages",
            "missing required key: git",
            "missing required key: parameters",
            "missing required key: metadata",
            "missing required key: inputs",
            "missing required key: outputs",
            "missing required key: generated",
            "missing required key: warnings",
        ]
    end
end

@testset "Git state outside repository" begin
    mktempdir() do dir
        git_state = OWENS.collect_git_state(dir)
        @test git_state isa OrderedCollections.OrderedDict{String,Any}
        @test git_state["available"] === false
        @test git_state["root"] == abspath(dir)
        @test isnothing(git_state["branch"])
        @test isnothing(git_state["commit"])
        @test isnothing(git_state["dirty"])
    end
end

@testset "Git metadata parsing" begin
    mktempdir() do dir
        repo_root = joinpath(dir, "repo")
        git_dir = joinpath(repo_root, ".git")
        ref_commit = "0123456789abcdef0123456789abcdef01234567"
        mkpath(joinpath(git_dir, "refs", "heads"))
        mkpath(joinpath(repo_root, "subdir"))
        write(joinpath(git_dir, "HEAD"), "ref: refs/heads/main\n")
        write(joinpath(git_dir, "refs", "heads", "main"), "$ref_commit\n")

        @test OWENS._git_repo_root(joinpath(repo_root, "subdir")) == abspath(repo_root)
        @test OWENS._git_dir(repo_root) == git_dir
        @test OWENS._git_head(git_dir) == ("main", ref_commit)
        @test OWENS._git_ref_commit(git_dir, "refs/heads/main") == ref_commit

        git_state = OWENS.collect_git_state(repo_root)
        @test git_state["available"] === true
        @test git_state["root"] == abspath(repo_root)
        @test git_state["branch"] == "main"
        @test git_state["commit"] == ref_commit
        @test git_state["dirty"] isa Union{Bool,Nothing}

        packed_git_dir = joinpath(dir, "packed.git")
        packed_commit = "fedcba9876543210fedcba9876543210fedcba98"
        mkpath(packed_git_dir)
        write(joinpath(packed_git_dir, "HEAD"), "ref: refs/heads/packed\n")
        write(
            joinpath(packed_git_dir, "packed-refs"),
            "# pack-refs with: peeled fully-peeled sorted\n$packed_commit refs/heads/packed\n",
        )
        @test OWENS._git_head(packed_git_dir) == ("packed", packed_commit)
        @test OWENS._packed_ref_commit(packed_git_dir, "refs/heads/missing") === nothing

        worktree_root = joinpath(dir, "worktree")
        worktree_git_dir = joinpath(dir, "actual.git")
        detached_commit = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        mkpath(worktree_root)
        mkpath(worktree_git_dir)
        write(joinpath(worktree_root, ".git"), "gitdir: ../actual.git\n")
        write(joinpath(worktree_git_dir, "HEAD"), "$detached_commit\n")
        @test OWENS._git_dir(worktree_root) == worktree_git_dir
        @test OWENS._git_head(worktree_git_dir) == (nothing, detached_commit)
    end
end
