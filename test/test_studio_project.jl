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
            "metadata",
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
        @test health["metadata"] == OrderedCollections.OrderedDict{String,Any}(
            "active" => true,
            "source" => "unit",
        )
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
    catalog = OWENS.studio_project_template_catalog()
    @test collect(keys(catalog)) == ["schema_version", "templates"]
    @test catalog["schema_version"] == "owens-studio-template-catalog/v1"
    @test catalog["templates"] == OrderedCollections.OrderedDict{String,Any}[
        OrderedCollections.OrderedDict{String,Any}(
            "template" => "blank",
            "title" => "Blank OWENS Studio Project",
            "description" => "Blank OWENS Studio project",
            "turbine_type" => "custom",
            "solver_path" => nothing,
            "creates_generated_script" => false,
            "creates_run_manifest" => false,
        ),
        OrderedCollections.OrderedDict{String,Any}(
            "template" => "rm2",
            "title" => "RM2 VAWT Template",
            "description" => "RM2 VAWT WindIO project",
            "turbine_type" => "VAWT",
            "solver_path" => "runOWENSWINDIO",
            "creates_generated_script" => true,
            "creates_run_manifest" => true,
        ),
    ]
    @test OWENS.studio_example_project_names() == ["rm2"]
    example_catalog = OWENS.studio_example_project_catalog()
    @test collect(keys(example_catalog)) == ["schema_version", "examples"]
    @test example_catalog["schema_version"] == "owens-studio-example-catalog/v1"
    @test length(example_catalog["examples"]) == 1
    @test example_catalog["examples"][1]["example"] == "rm2"
    @test example_catalog["examples"][1]["project_relative_path"] ==
          joinpath("examples", "gui", "rm2", "owens_project.yml")
    @test example_catalog["examples"][1]["available"] === true
    @test isfile(example_catalog["examples"][1]["project_file"])

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
        @test OWENS.studio_project_generated_script_path(created["project_file"]) ==
              created["script_file"]
        @test OWENS.read_studio_project_generated_script(created["project_file"]) ==
              created["script"]
        bundle = OWENS.write_studio_workbench_bundle(
            joinpath(dir, "rm2-bundle"),
            created["project_file"],
        )
        @test bundle["schema_version"] == "owens-studio-bundle/v1"
        @test bundle["project_file"] == created["project_file"]
        @test bundle["project_status"] == "ok"
        @test isdir(bundle["bundle_dir"])
        @test isfile(bundle["index_html"])
        @test isfile(bundle["health_file"])
        @test isfile(bundle["script_file"])
        @test bundle["bytes"]["index_html"] == stat(bundle["index_html"]).size
        @test bundle["bytes"]["health_file"] == stat(bundle["health_file"]).size
        @test bundle["bytes"]["script_file"] == stat(bundle["script_file"]).size
        @test read(bundle["script_file"], String) == created["script"]
        bundle_html = read(bundle["index_html"], String)
        @test occursin("Health YAML", bundle_html)
        @test occursin("Generated Julia", bundle_html)
        @test occursin("health.yml", bundle_html)
        @test occursin("generated_script.jl", bundle_html)
        bundle_health = YAML.load_file(
            bundle["health_file"];
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_health["status"] == "ok"
        @test bundle_health["metadata"]["generated_script"] ==
              joinpath("runs", "rm2", "run_rm2_windio.jl")

        health = OWENS.studio_project_health(created["project_file"])
        @test health["status"] == "ok"
        @test health["metadata"]["generated_script"] ==
              joinpath("runs", "rm2", "run_rm2_windio.jl")
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
        @test OWENS.studio_project_generated_script_path(blank["project_file"]) === nothing
        @test OWENS.read_studio_project_generated_script(
            blank["project_file"];
            required = false,
        ) === nothing
        @test_throws ArgumentError OWENS.read_studio_project_generated_script(
            blank["project_file"],
        )
    end
end

@testset "OWENS Studio GUI fixtures" begin
    fixture_project =
        normpath(joinpath(@__DIR__, "..", "examples", "gui", "rm2", "owens_project.yml"))
    fixture_root = dirname(fixture_project)
    fixture_run_manifest = joinpath(fixture_root, "runs", "rm2", "run_manifest.yml")
    fixture_script = joinpath(fixture_root, "runs", "rm2", "run_rm2_windio.jl")

    @test isfile(fixture_project)
    @test isfile(fixture_run_manifest)
    @test isfile(fixture_script)
    @test OWENS.studio_project_issues(fixture_project) == String[]

    health = OWENS.studio_project_health(fixture_project)
    @test health["schema_version"] == "owens-studio-workbench/v1"
    @test health["status"] == "ok"
    @test health["project_id"] == "rm2-gui-fixture"
    @test health["root"] == fixture_root
    @test health["summary"] == OrderedCollections.OrderedDict{String,Any}(
        "records" => 3,
        "ok" => 3,
        "modified" => 0,
        "missing" => 0,
        "invalid_record" => 0,
    )
    @test [row["path"] for row in health["files"]] == [
        joinpath("..", "..", "RM2", "modeling_options_OWENS_RM2.yml"),
        joinpath("..", "..", "RM2", "WINDIO_RM2.yaml"),
    ]
    @test health["runs"][1]["resolved_path"] == fixture_run_manifest
    @test health["runs"][1]["run_manifest_health"]["status"] == "ok"
    @test health["runs"][1]["run_manifest_health"]["root"] == fixture_root
    @test health["runs"][1]["run_manifest_health"]["summary"] ==
          OrderedCollections.OrderedDict{String,Any}(
        "records" => 3,
        "ok" => 3,
        "modified" => 0,
        "missing" => 0,
        "invalid_record" => 0,
    )
    @test OWENS.studio_project_generated_script_path(fixture_project) == fixture_script
    @test occursin(
        "rm2_input_root = normpath(joinpath(@__DIR__",
        OWENS.read_studio_project_generated_script(fixture_project),
    )

    open_payload = OWENS_APP.open_studio_project(fixture_project)
    @test open_payload["schema_version"] == "owens-studio-open/v1"
    @test open_payload["project_status"] == "ok"
    @test open_payload["generated_script"]["path"] == fixture_script
    @test open_payload["generated_script"]["relative_path"] ==
          joinpath("runs", "rm2", "run_rm2_windio.jl")
    @test open_payload["generated_script"]["available"] === true
    @test open_payload["generated_script"]["bytes"] == stat(fixture_script).size
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

        template_catalog = OWENS_APP.list_studio_project_templates()
        @test template_catalog["schema_version"] == "owens-studio-template-catalog/v1"
        @test [row["template"] for row in template_catalog["templates"]] == ["blank", "rm2"]
        example_catalog = OWENS_APP.list_studio_example_projects()
        @test example_catalog["schema_version"] == "owens-studio-example-catalog/v1"
        @test [row["example"] for row in example_catalog["examples"]] == ["rm2"]
        route_catalog = OWENS_APP.studio_route_catalog()
        @test route_catalog["schema_version"] == "owens-studio-route-catalog/v1"
        @test [row["name"] for row in route_catalog["routes"]] == [
            "route_catalog",
            "template_catalog",
            "example_catalog",
            "project_open",
            "project_health",
            "project_workbench",
            "project_script",
            "project_bundle",
            "create_template_project",
        ]
        @test [row["method"] for row in route_catalog["routes"]] == ["GET", "GET", "GET", "GET", "GET", "GET", "GET", "POST", "POST"]
        @test route_catalog["routes"][4]["required_params"] == ["project_path"]
        @test route_catalog["routes"][4]["optional_params"] == ["summarize_runs"]
        @test route_catalog["routes"][6]["content_type"] == "text/html; charset=utf-8"
        @test route_catalog["routes"][8]["required_params"] ==
              ["project_path", "output_dir"]

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
        templates_cli = OWENS_APP.real_main(["project-templates"]; io = IOBuffer())
        @test templates_cli["schema_version"] == "owens-studio-template-catalog/v1"
        @test templates_cli["templates"][2]["solver_path"] == "runOWENSWINDIO"
        examples_cli = OWENS_APP.real_main(["project-examples"]; io = IOBuffer())
        @test examples_cli["schema_version"] == "owens-studio-example-catalog/v1"
        @test examples_cli["examples"][1]["example"] == "rm2"
        routes_cli = OWENS_APP.real_main(["project-routes"]; io = IOBuffer())
        @test routes_cli["schema_version"] == "owens-studio-route-catalog/v1"
        @test routes_cli["routes"][1]["path"] == "/api/routes"
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
        open_payload = OWENS_APP.open_studio_project(template_cli["project_file"])
        @test open_payload["schema_version"] == "owens-studio-open/v1"
        @test open_payload["project_file"] == template_cli["project_file"]
        @test open_payload["project_status"] == "ok"
        @test open_payload["generated_script"]["path"] == template_cli["script_file"]
        @test open_payload["generated_script"]["relative_path"] ==
              joinpath("runs", "rm2", "run_rm2_windio.jl")
        @test open_payload["generated_script"]["available"] === true
        @test open_payload["generated_script"]["sha256"] ==
              OWENS.file_sha256(template_cli["script_file"])
        @test [row["route"] for row in open_payload["actions"]] == ["project_health", "project_workbench", "project_script", "project_bundle"]
        @test open_payload["routes"]["schema_version"] == "owens-studio-route-catalog/v1"
        @test open_payload["templates"]["schema_version"] ==
              "owens-studio-template-catalog/v1"
        @test open_payload["examples"]["schema_version"] ==
              "owens-studio-example-catalog/v1"
        open_cli = OWENS_APP.real_main(
            ["project-open", template_cli["project_file"]];
            io = IOBuffer(),
        )
        @test open_cli["schema_version"] == "owens-studio-open/v1"
        @test open_cli["generated_script"]["available"] === true
        script_cli = OWENS_APP.real_main(
            ["project-script", template_cli["project_file"]];
            io = IOBuffer(),
        )
        @test script_cli["script_file"] == template_cli["script_file"]
        @test occursin("OWENS.runOWENSWINDIO", script_cli["script"])
        bundle_cli = OWENS_APP.real_main(
            ["project-bundle", template_cli["project_file"], joinpath(dir, "bundle-cli")];
            io = IOBuffer(),
        )
        @test bundle_cli["schema_version"] == "owens-studio-bundle/v1"
        @test bundle_cli["project_status"] == "ok"
        @test isfile(bundle_cli["index_html"])
        @test isfile(bundle_cli["health_file"])
        @test isfile(bundle_cli["script_file"])
        @test isfile(bundle_cli["open_file"])
        @test bundle_cli["bytes"]["index_html"] == stat(bundle_cli["index_html"]).size
        @test bundle_cli["bytes"]["open_file"] == stat(bundle_cli["open_file"]).size
        bundle_index = read(bundle_cli["index_html"], String)
        @test occursin("Open Payload", bundle_index)
        @test occursin("open.yml", bundle_index)
        bundle_open = YAML.load_file(
            bundle_cli["open_file"];
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_open["schema_version"] == "owens-studio-open/v1"
        @test bundle_open["generated_script"]["available"] === true
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

@testset "OWENS Studio app route handlers" begin
    mktempdir() do dir
        created = OWENS.create_studio_project_template(
            joinpath(dir, "rm2-route");
            template = "rm2",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        project_file = created["project_file"]

        routes_response = OWENS_APP.studio_routes_route()
        @test routes_response.status == 200
        @test routes_response.content_type == "application/x-yaml; charset=utf-8"
        routes_payload = YAML.load(
            routes_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test routes_payload["schema_version"] == "owens-studio-route-catalog/v1"
        @test routes_payload["routes"][1]["handler"] == "studio_routes_route"
        @test routes_payload["routes"][end]["handler"] == "studio_project_template_route"

        templates_response = OWENS_APP.studio_project_templates_route()
        @test templates_response.status == 200
        @test templates_response.content_type == "application/x-yaml; charset=utf-8"
        templates_payload = YAML.load(
            templates_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test templates_payload["schema_version"] == "owens-studio-template-catalog/v1"
        @test templates_payload["templates"][1]["template"] == "blank"
        @test templates_payload["templates"][2]["template"] == "rm2"

        examples_response = OWENS_APP.studio_project_examples_route()
        @test examples_response.status == 200
        @test examples_response.content_type == "application/x-yaml; charset=utf-8"
        examples_payload = YAML.load(
            examples_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test examples_payload["schema_version"] == "owens-studio-example-catalog/v1"
        @test examples_payload["examples"][1]["example"] == "rm2"
        @test examples_payload["examples"][1]["available"] === true
        dispatch_examples = OWENS_APP.dispatch_studio_route("/api/examples"; method = "GET")
        @test dispatch_examples.status == 200
        @test occursin("owens-studio-example-catalog/v1", dispatch_examples.body)

        open_response = OWENS_APP.studio_project_open_route(project_file)
        @test open_response.status == 200
        @test open_response.content_type == "application/x-yaml; charset=utf-8"
        open_payload = YAML.load(
            open_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test open_payload["schema_version"] == "owens-studio-open/v1"
        @test open_payload["project_status"] == "ok"
        @test open_payload["generated_script"]["available"] === true
        @test open_payload["actions"][3]["route"] == "project_script"
        @test open_payload["actions"][3]["enabled"] === true

        dispatch_open = OWENS_APP.dispatch_studio_route(
            "/api/project/open";
            method = "GET",
            params = (; project_path = project_file),
        )
        @test dispatch_open.status == 200
        @test dispatch_open.content_type == "application/x-yaml; charset=utf-8"
        dispatch_open_payload = YAML.load(
            dispatch_open.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test dispatch_open_payload["project_file"] == project_file
        @test dispatch_open_payload["routes"]["routes"][4]["name"] == "project_open"

        health_response = OWENS_APP.studio_project_health_route(project_file)
        @test health_response isa OWENS_APP.StudioRouteResponse
        @test health_response.status == 200
        @test health_response.content_type == "application/x-yaml; charset=utf-8"
        health_payload = YAML.load(
            health_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test health_payload["schema_version"] == "owens-studio-workbench/v1"
        @test health_payload["status"] == "ok"
        @test health_payload["summary"]["records"] == 3
        @test health_payload["runs"][1]["run_manifest_health"]["summary"]["ok"] == 3

        html_response = OWENS_APP.studio_project_workbench_route(project_file)
        @test html_response.status == 200
        @test html_response.content_type == "text/html; charset=utf-8"
        @test occursin(
            "<title>RM2 VAWT Template - OWENS Studio</title>",
            html_response.body,
        )
        @test occursin("run_manifest.yml", html_response.body)
        @test occursin("Generated Script", html_response.body)
        @test occursin("run_rm2_windio.jl", html_response.body)

        script_response = OWENS_APP.studio_project_script_route(project_file)
        @test script_response.status == 200
        @test script_response.content_type == "text/plain; charset=utf-8"
        @test script_response.body == created["script"]
        @test occursin("OWENS.runOWENSWINDIO", script_response.body)

        bundle_response = OWENS_APP.studio_project_bundle_route(
            project_file,
            joinpath(dir, "route-bundle"),
        )
        @test bundle_response.status == 200
        @test bundle_response.content_type == "application/x-yaml; charset=utf-8"
        bundle_payload = YAML.load(
            bundle_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_payload["schema_version"] == "owens-studio-bundle/v1"
        @test bundle_payload["project_status"] == "ok"
        @test isfile(bundle_payload["index_html"])
        @test isfile(bundle_payload["health_file"])
        @test isfile(bundle_payload["script_file"])
        @test isfile(bundle_payload["open_file"])
        @test bundle_payload["bytes"]["index_html"] ==
              stat(bundle_payload["index_html"]).size
        @test bundle_payload["bytes"]["open_file"] == stat(bundle_payload["open_file"]).size
        bundle_payload_index = read(bundle_payload["index_html"], String)
        @test occursin("Open Payload", bundle_payload_index)
        @test occursin("open.yml", bundle_payload_index)
        bundle_open_payload = YAML.load_file(
            bundle_payload["open_file"];
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_open_payload["project_file"] == project_file
        @test bundle_open_payload["actions"][4]["route"] == "project_bundle"

        dispatch_bundle = OWENS_APP.dispatch_studio_route(
            "project_bundle";
            method = "POST",
            params = Dict(
                :project_path => project_file,
                :output_dir => joinpath(dir, "dispatch-bundle"),
                :include_script => false,
            ),
        )
        @test dispatch_bundle.status == 200
        dispatch_bundle_payload = YAML.load(
            dispatch_bundle.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test isfile(dispatch_bundle_payload["open_file"])
        @test dispatch_bundle_payload["script_file"] === nothing

        bad_method = OWENS_APP.dispatch_studio_route(
            "project_open";
            method = "POST",
            params = (; project_path = project_file),
        )
        @test bad_method.status == 405
        @test occursin("not allowed", bad_method.body)
        missing_param = OWENS_APP.dispatch_studio_route("project_open")
        @test missing_param.status == 400
        @test occursin("Missing Studio route parameter: project_path", missing_param.body)
        missing_route = OWENS_APP.dispatch_studio_route("missing_route")
        @test missing_route.status == 404
        @test occursin("Unknown Studio route", missing_route.body)

        template_response = OWENS_APP.studio_project_template_route(
            joinpath(dir, "created-from-route");
            template = "blank",
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        @test template_response.status == 200
        @test template_response.content_type == "application/x-yaml; charset=utf-8"
        template_payload = YAML.load(
            template_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test template_payload["template"] == "blank"
        @test template_payload["project_status"] == "ok"
        @test template_payload["run_manifest_file"] === nothing
        @test isfile(template_payload["project_file"])

        error_response = OWENS_APP.studio_project_health_route(joinpath(dir, "missing.yml"))
        @test error_response.status == 400
        @test error_response.content_type == "application/x-yaml; charset=utf-8"
        error_payload = YAML.load(
            error_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test error_payload["status"] == "error"
        @test occursin("Cannot read missing Studio project", error_payload["message"])
    end
end
