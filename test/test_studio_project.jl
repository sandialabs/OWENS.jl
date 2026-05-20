using Test
import HDF5
import OWENS
import OrderedCollections
import YAML

include(joinpath(@__DIR__, "..", "app", "OWENS_APP", "src", "OWENS_APP.jl"))

const STUDIO_MODEL_SHA256 = "5fdc1fb3c0b14924ab13cbe75816f20ada7685ccb6ad4ca8a5103251233ddaf0"
const STUDIO_WINDIO_SHA256 = "8c6ed05c7c0f22c45fc5acea73c206ab9ca0b1b7d62b7fbb6b246cb3f080e496"
const STUDIO_RM2_MODEL_SHA256 = "df24a053994c15fa83dcab09846d1401b14f892478333875971a320de9d4e94a"
const STUDIO_RM2_WINDIO_SHA256 = "18fbfb761fe866e18d6fb24ed6f5800c26f7dcca225cf7f0e859729a23e74c3c"

@testset "OWENS Studio project manifest and health" begin
    mktempdir() do dir
        model_file = joinpath(dir, "modeling_options.yml")
        windio_file = joinpath(dir, "design.yml")
        run_dir = joinpath(dir, "runs", "run-001")
        output_file = joinpath(run_dir, "output.h5")
        manifest_file = joinpath(run_dir, "run_manifest.yml")
        project_file = joinpath(dir, "owens_project.yml")
        html_file = joinpath(dir, "workbench.html")
        mkpath(run_dir)
        write(model_file, "OWENS_Options:\n  numTS: 2\n")
        write(windio_file, "name: unit\n")
        HDF5.h5open(output_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0, 3.0])
        end

        OWENS.write_run_manifest(
            manifest_file;
            run_id = "run-001",
            run_name = "studio run",
            project_root = run_dir,
            solver = "unit-solver",
            input_files = [model_file, windio_file],
            output_files = [output_file],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        project = OWENS.build_studio_project(
            dir;
            project_id = "studio-unit",
            name = "Studio Unit",
            description = "Headless workbench unit fixture",
            modeling_options_file = model_file,
            windio_file,
            run_manifests = [manifest_file],
            metadata = Dict(:source => :unit, :active => true),
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        @test collect(keys(project)) == [
            "schema_version",
            "project_id",
            "name",
            "description",
            "created_at_utc",
            "updated_at_utc",
            "root",
            "files",
            "runs",
            "metadata",
        ]
        @test project["schema_version"] == "owens-studio-project/v1"
        @test project["project_id"] == "studio-unit"
        @test project["name"] == "Studio Unit"
        @test project["description"] == "Headless workbench unit fixture"
        @test project["created_at_utc"] == "2026-05-20T00:00:00.000Z"
        @test project["updated_at_utc"] == "2026-05-20T00:00:00.000Z"
        @test project["root"] == abspath(dir)
        @test length(project["files"]) == 2
        @test project["files"][1]["path"] == "modeling_options.yml"
        @test project["files"][1]["role"] == "modeling_options"
        @test project["files"][1]["sha256"] == STUDIO_MODEL_SHA256
        @test project["files"][2]["path"] == "design.yml"
        @test project["files"][2]["role"] == "windio"
        @test project["files"][2]["sha256"] == STUDIO_WINDIO_SHA256
        @test length(project["runs"]) == 1
        @test project["runs"][1]["path"] == joinpath("runs", "run-001", "run_manifest.yml")
        @test project["runs"][1]["role"] == "run_manifest"
        @test collect(keys(project["metadata"])) == ["active", "source"]
        @test project["metadata"]["active"] === true
        @test project["metadata"]["source"] == "unit"
        @test OWENS.studio_project_issues(project) == String[]
        @test OWENS.validate_studio_project(project) === project

        written_project = OWENS.write_studio_project(project_file, project)
        loaded_project = OWENS.read_studio_project(project_file)
        @test isfile(project_file)
        @test written_project["project_id"] == "studio-unit"
        @test loaded_project["project_id"] == "studio-unit"
        @test OWENS.studio_project_issues(project_file) == String[]

        health = OWENS.studio_project_health(project_file)
        @test collect(keys(health)) == [
            "schema_version",
            "status",
            "project_path",
            "root",
            "project_id",
            "name",
            "project_issues",
            "summary",
            "files",
            "runs",
        ]
        @test health["schema_version"] == "owens-studio-workbench/v1"
        @test health["status"] == "ok"
        @test health["project_path"] == abspath(project_file)
        @test health["root"] == abspath(dir)
        @test health["project_id"] == "studio-unit"
        @test health["name"] == "Studio Unit"
        @test health["project_issues"] == String[]
        @test health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 3,
            "modified" => 0,
            "missing" => 0,
            "invalid_record" => 0,
        )
        @test [row["status"] for row in health["files"]] == ["ok", "ok"]
        @test health["runs"][1]["status"] == "ok"
        @test health["runs"][1]["run_manifest_health"]["status"] == "ok"
        @test health["runs"][1]["run_manifest_health"]["outputs"][1]["output_data_summary"][1]["name"] ==
              "t"
        @test health["runs"][1]["run_manifest_health"]["outputs"][1]["output_data_summary"][1]["present"] ===
              true

        no_run_summary = OWENS.studio_project_health(project; summarize_runs = false)
        @test no_run_summary["runs"][1]["status"] == "ok"
        @test !haskey(no_run_summary["runs"][1], "run_manifest_health")

        html = OWENS.render_studio_workbench_html(health)
        @test occursin("<title>Studio Unit - OWENS Studio</title>", html)
        @test occursin("<span class=\"status ok\">ok</span>", html)
        @test occursin("<h3>Project Files</h3>", html)
        @test occursin("modeling_options.yml", html)
        @test occursin("run_manifest.yml", html)
        @test occursin("No schema issues.", html)
        written_html = OWENS.write_studio_workbench_html(html_file, project_file)
        @test written_html == read(html_file, String)
        @test occursin("OWENS Studio", written_html)

        stale_project = deepcopy(project)
        stale_project["root"] = joinpath(dir, "stale")
        stale_health = OWENS.studio_project_health(stale_project)
        @test stale_health["status"] == "attention"
        @test stale_health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 0,
            "modified" => 0,
            "missing" => 3,
            "invalid_record" => 0,
        )
        override_health = OWENS.studio_project_health(stale_project; root = dir)
        @test override_health["status"] == "ok"
        @test override_health["summary"] == health["summary"]

        malformed = deepcopy(project)
        malformed["schema_version"] = "owens-studio-project/v0"
        malformed["files"] = [Dict("path" => 7, "bytes" => true, "sha256" => "bad")]
        @test OWENS.studio_project_issues(malformed) == [
            "schema_version must equal owens-studio-project/v1",
            "files[1].path must be a string",
            "files[1].bytes must be a non-negative integer",
            "files[1].sha256 must be a lowercase 64-character SHA-256 digest",
        ]
        @test_throws ArgumentError OWENS.validate_studio_project(malformed)
    end
end

@testset "OWENS Studio template projects" begin
    @test OWENS.studio_project_template_names() == ["blank", "rm2"]

    mktempdir() do dir
        target = joinpath(dir, "rm2-studio")
        created = OWENS.create_studio_project_template(
            target;
            template = "rm2",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        @test collect(keys(created)) == [
            "template",
            "project_file",
            "project",
            "run_manifest_file",
            "run_manifest",
            "script_file",
            "script",
        ]
        @test created["template"] == "rm2"
        @test created["project_file"] == abspath(joinpath(target, "owens_project.yml"))
        @test created["run_manifest_file"] ==
              abspath(joinpath(target, "runs", "rm2", "run_manifest.yml"))
        @test created["script_file"] ==
              abspath(joinpath(target, "runs", "rm2", "run_rm2_windio.jl"))
        @test isfile(created["project_file"])
        @test isfile(created["run_manifest_file"])
        @test isfile(created["script_file"])

        project = created["project"]
        @test project["schema_version"] == "owens-studio-project/v1"
        @test project["project_id"] == "rm2"
        @test project["name"] == "RM2 VAWT Template"
        @test project["root"] == abspath(target)
        @test project["metadata"] == OrderedCollections.OrderedDict{String,Any}(
            "generated_script" => joinpath("runs", "rm2", "run_rm2_windio.jl"),
            "template" => "rm2",
            "template_description" => "RM2 VAWT WindIO project",
        )
        @test project["files"][1]["path"] ==
              joinpath("inputs", "modeling_options_OWENS_RM2.yml")
        @test project["files"][1]["role"] == "modeling_options"
        @test project["files"][1]["sha256"] == STUDIO_RM2_MODEL_SHA256
        @test project["files"][2]["path"] == joinpath("inputs", "WINDIO_RM2.yaml")
        @test project["files"][2]["role"] == "windio"
        @test project["files"][2]["sha256"] == STUDIO_RM2_WINDIO_SHA256
        @test project["runs"][1]["path"] == joinpath("runs", "rm2", "run_manifest.yml")
        @test project["runs"][1]["role"] == "run_manifest"

        manifest = created["run_manifest"]
        @test manifest["schema_version"] == "owens-run-manifest/v1"
        @test manifest["run_id"] == "rm2-template"
        @test manifest["run_name"] == "RM2 Studio Template"
        @test manifest["solver"] == "runOWENSWINDIO"
        @test manifest["project_root"] == abspath(joinpath(target, "runs", "rm2"))
        @test manifest["status"] == "created"
        @test manifest["metadata"] == OrderedCollections.OrderedDict{String,Any}(
            "template" => "rm2",
            "template_description" => "RM2 VAWT WindIO project",
        )
        @test manifest["inputs"][1]["path"] ==
              joinpath("..", "..", "inputs", "modeling_options_OWENS_RM2.yml")
        @test manifest["inputs"][1]["sha256"] == STUDIO_RM2_MODEL_SHA256
        @test manifest["inputs"][2]["path"] ==
              joinpath("..", "..", "inputs", "WINDIO_RM2.yaml")
        @test manifest["inputs"][2]["sha256"] == STUDIO_RM2_WINDIO_SHA256
        @test length(manifest["generated"]) == 1
        @test manifest["generated"][1]["path"] == "run_rm2_windio.jl"
        @test manifest["generated"][1]["role"] == "generated"

        @test occursin("OWENS.runOWENSWINDIO", created["script"])
        @test occursin(
            repr(abspath(joinpath(target, "inputs", "modeling_options_OWENS_RM2.yml"))),
            created["script"],
        )
        @test read(created["script_file"], String) == created["script"]

        health = OWENS.studio_project_health(created["project_file"])
        @test health["status"] == "ok"
        @test health["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 3,
            "modified" => 0,
            "missing" => 0,
            "invalid_record" => 0,
        )
        @test health["runs"][1]["run_manifest_health"]["summary"] ==
              OrderedCollections.OrderedDict{String,Any}(
            "records" => 3,
            "ok" => 3,
            "modified" => 0,
            "missing" => 0,
            "invalid_record" => 0,
        )

        @test_throws ArgumentError OWENS.create_studio_project_template(
            target;
            template = "rm2",
        )
        @test_throws ArgumentError OWENS.create_studio_project_template(
            joinpath(dir, "bad-template");
            template = "missing",
        )

        blank_target = joinpath(dir, "blank-studio")
        blank = OWENS.create_studio_project_template(
            blank_target;
            template = "blank",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        @test blank["template"] == "blank"
        @test blank["run_manifest_file"] === nothing
        @test blank["script_file"] === nothing
        @test blank["project"]["files"] == OrderedCollections.OrderedDict{String,Any}[]
        @test blank["project"]["runs"] == OrderedCollections.OrderedDict{String,Any}[]
        @test OWENS.studio_project_health(blank["project_file"])["status"] == "ok"
    end
end

@testset "OWENS Studio app services" begin
    mktempdir() do dir
        model_file = joinpath(dir, "modeling_options.yml")
        windio_file = joinpath(dir, "design.yml")
        run_dir = joinpath(dir, "run")
        output_file = joinpath(run_dir, "output.h5")
        manifest_file = joinpath(run_dir, "run_manifest.yml")
        project_file = joinpath(dir, "owens_project.yml")
        html_file = joinpath(dir, "studio.html")
        mkpath(run_dir)
        write(model_file, "OWENS_Options:\n  numTS: 2\n")
        write(windio_file, "name: unit\n")
        HDF5.h5open(output_file, "w") do file
            HDF5.write(file, "t", [1.0, 2.0, 3.0])
        end
        OWENS.write_run_manifest(
            manifest_file;
            run_id = "run-app",
            run_name = "app run",
            project_root = run_dir,
            solver = "unit-solver",
            input_files = [model_file, windio_file],
            output_files = [output_file],
            status = "complete",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        OWENS.write_studio_project(
            project_file,
            dir;
            project_id = "studio-app",
            name = "Studio App",
            modeling_options_file = model_file,
            windio_file,
            run_manifests = [manifest_file],
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        manifest_health = OWENS_APP.inspect_run_manifest(manifest_file)
        @test manifest_health["status"] == "ok"
        @test manifest_health["summary"]["records"] == 3

        output_summary = OWENS_APP.inspect_output_data(output_file; channels = ["t"])
        @test output_summary["path"] == abspath(output_file)
        @test length(output_summary["channels"]) == 1
        @test output_summary["channels"][1]["name"] == "t"
        @test output_summary["channels"][1]["shape"] == [3]
        @test output_summary["channels"][1]["units"] == "s"

        prepared = OWENS_APP.prepare_windio_run(
            model_file,
            windio_file,
            run_dir;
            run_id = "windio-app",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        @test prepared["spec"]["modeling_options_file"] == abspath(model_file)
        @test prepared["spec"]["windio_file"] == abspath(windio_file)
        @test prepared["spec"]["run_path"] == abspath(run_dir)
        @test occursin("OWENS.runOWENSWINDIO", prepared["script"])
        @test prepared["manifest"]["run_id"] == "windio-app"
        @test prepared["manifest"]["solver"] == "runOWENSWINDIO"

        project_health = OWENS_APP.inspect_studio_project(project_file)
        @test project_health["status"] == "ok"
        @test project_health["name"] == "Studio App"
        html = OWENS_APP.write_studio_project_workbench(html_file, project_health)
        @test html == read(html_file, String)
        @test occursin("Studio App", html)

        manifest_cli =
            OWENS_APP.real_main(["manifest-health", manifest_file]; io = IOBuffer())
        @test manifest_cli["status"] == "ok"
        summary_cli = OWENS_APP.real_main(["output-summary", output_file]; io = IOBuffer())
        @test summary_cli["channels"][1]["name"] == "t"
        windio_cli = OWENS_APP.real_main(
            ["windio-script", model_file, windio_file, run_dir];
            io = IOBuffer(),
        )
        @test occursin("OWENS.runOWENSWINDIO", windio_cli["script"])
        template_cli = OWENS_APP.real_main(
            ["project-template", "rm2", joinpath(dir, "template-cli")];
            io = IOBuffer(),
        )
        @test template_cli["template"] == "rm2"
        @test template_cli["project_status"] == "ok"
        @test template_cli["project_health"]["summary"]["records"] == 3
        @test isfile(template_cli["project_file"])
        @test isfile(template_cli["run_manifest_file"])
        @test isfile(template_cli["script_file"])
        project_cli = OWENS_APP.real_main(["project-health", project_file]; io = IOBuffer())
        @test project_cli["status"] == "ok"
        project_html_cli =
            OWENS_APP.real_main(["project-html", project_file, html_file]; io = IOBuffer())
        @test project_html_cli["output_html"] == abspath(html_file)
        @test project_html_cli["project_status"] == "ok"

        @test_throws ArgumentError OWENS_APP.real_main(String[]; io = IOBuffer())
        @test_throws ArgumentError OWENS_APP.real_main(["bad-command"]; io = IOBuffer())
    end
end
