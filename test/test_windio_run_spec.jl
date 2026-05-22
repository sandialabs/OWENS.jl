using Test
import OWENS

const WINDIO_SPEC_MODEL_SHA256 = "5fdc1fb3c0b14924ab13cbe75816f20ada7685ccb6ad4ca8a5103251233ddaf0"
const WINDIO_SPEC_DESIGN_SHA256 = "8c6ed05c7c0f22c45fc5acea73c206ab9ca0b1b7d62b7fbb6b246cb3f080e496"
const WINDIO_SPEC_OUTPUT_SHA256 = "a02b48c12e33e850387f1890f7aaa11c845f5d62f783919978585d2d52c9c7aa"

@testset "WindIO run spec paths" begin
    mktempdir() do dir
        model_file = joinpath(dir, "modeling_options.yml")
        windio_file = joinpath(dir, "design.yml")
        run_path = joinpath(dir, "run")
        write(model_file, "OWENS_Options:\n  numTS: 2\n")
        write(windio_file, "name: unit\n")
        mkdir(run_path)

        spec = OWENS.windio_run_spec(model_file, windio_file, run_path)
        @test spec isa OWENS.WindIORunSpec
        @test spec.modeling_options_file == abspath(model_file)
        @test spec.windio_file == abspath(windio_file)
        @test spec.run_path == abspath(run_path)

        created_run_path = joinpath(dir, "created_run")
        created_spec = OWENS.windio_run_spec(
            model_file,
            windio_file,
            created_run_path;
            create_run_path = true,
        )
        @test isdir(created_run_path)
        @test created_spec.run_path == abspath(created_run_path)

        @test_throws ArgumentError OWENS.windio_run_spec(
            joinpath(dir, "missing_model.yml"),
            windio_file,
            run_path,
        )
        @test_throws ArgumentError OWENS.windio_run_spec(
            model_file,
            joinpath(dir, "missing_design.yml"),
            run_path,
        )
        @test_throws ArgumentError OWENS.windio_run_spec(
            model_file,
            windio_file,
            joinpath(dir, "missing_run"),
        )
    end
end

@testset "WindIO run manifest export" begin
    mktempdir() do dir
        model_file = joinpath(dir, "modeling_options.yml")
        windio_file = joinpath(dir, "design.yml")
        run_path = joinpath(dir, "run")
        output_file = joinpath(run_path, "output.h5")
        manifest_file = joinpath(run_path, "run_manifest.yml")
        write(model_file, "OWENS_Options:\n  numTS: 2\n")
        write(windio_file, "name: unit\n")
        mkdir(run_path)
        write(output_file, "output: done\n")

        spec = OWENS.windio_run_spec(model_file, windio_file, run_path)
        manifest = OWENS.build_windio_run_manifest(
            spec;
            run_id = "windio-unit-001",
            run_name = "windio unit",
            output_files = [output_file],
            parameters = Dict(:case => "manifest"),
            metadata = Dict(:source => "test"),
            warnings = ["unit warning"],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        @test manifest["schema_version"] == "owens-run-manifest/v1"
        @test manifest["run_id"] == "windio-unit-001"
        @test manifest["run_name"] == "windio unit"
        @test manifest["solver"] == "runOWENSWINDIO"
        @test manifest["project_root"] == abspath(run_path)
        @test manifest["status"] == "complete"
        @test manifest["created_at_utc"] == "2026-05-20T00:00:00.000Z"
        @test manifest["parameters"]["case"] == "manifest"
        @test manifest["metadata"]["source"] == "test"
        @test manifest["warnings"] == ["unit warning"]

        @test length(manifest["inputs"]) == 2
        @test manifest["inputs"][1]["path"] ==
              relpath(abspath(model_file), abspath(run_path))
        @test manifest["inputs"][1]["role"] == "modeling_options"
        @test manifest["inputs"][1]["sha256"] == WINDIO_SPEC_MODEL_SHA256
        @test manifest["inputs"][2]["path"] ==
              relpath(abspath(windio_file), abspath(run_path))
        @test manifest["inputs"][2]["role"] == "windio"
        @test manifest["inputs"][2]["sha256"] == WINDIO_SPEC_DESIGN_SHA256

        @test length(manifest["outputs"]) == 1
        @test manifest["outputs"][1]["path"] == "output.h5"
        @test manifest["outputs"][1]["role"] == "output"
        @test manifest["outputs"][1]["sha256"] == WINDIO_SPEC_OUTPUT_SHA256

        written = OWENS.write_windio_run_manifest(
            manifest_file,
            spec;
            run_id = "windio-unit-002",
            output_files = [output_file],
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        loaded = OWENS.read_run_manifest(manifest_file)
        @test written["run_id"] == "windio-unit-002"
        @test loaded["run_id"] == "windio-unit-002"
        @test loaded["inputs"][1]["sha256"] == WINDIO_SPEC_MODEL_SHA256
        @test loaded["inputs"][2]["sha256"] == WINDIO_SPEC_DESIGN_SHA256
        @test loaded["outputs"][1]["sha256"] == WINDIO_SPEC_OUTPUT_SHA256
    end
end

@testset "WindIO run script export" begin
    mktempdir() do dir
        model_file = joinpath(dir, "modeling options.yml")
        windio_file = joinpath(dir, "design data.yml")
        run_path = joinpath(dir, "run path")
        script_file = joinpath(dir, "scripts", "run_windio.jl")
        write(model_file, "OWENS_Options:\n  numTS: 2\n")
        write(windio_file, "name: unit\n")
        mkdir(run_path)

        spec = OWENS.windio_run_spec(model_file, windio_file, run_path)
        expected_script = join(
            [
                "# Generated by OWENS.render_windio_run_script.",
                "using OWENS",
                "",
                "modeling_options_file = $(repr(abspath(model_file)))",
                "windio_file = $(repr(abspath(windio_file)))",
                "run_path = $(repr(abspath(run_path)))",
                "",
                "OWENS.runOWENSWINDIO(modeling_options_file, windio_file, run_path)",
                "",
            ],
            "\n",
        )

        @test OWENS.render_windio_run_script(spec) == expected_script
        written_script = OWENS.write_windio_run_script(script_file, spec)
        @test written_script == expected_script
        @test read(script_file, String) == expected_script
    end
end
