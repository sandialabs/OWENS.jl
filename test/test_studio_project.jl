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

_portable_path(path::AbstractString) = replace(path, '\\' => '/')
_portable_paths(paths) = [_portable_path(path) for path in paths]

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
        @test occursin(
            "meta name=\"owens-studio-style\" content=\"owens-studio-shared-style/v1\"",
            html,
        )
        @test occursin("--ink: #17202a;", html)
        @test occursin("<span class=\"status ok\">ok</span>", html)
        @test occursin("<h3>Project Files</h3>", html)
        @test occursin("modeling_options.yml", html)
        @test occursin("run_manifest.yml", html)
        @test occursin("No schema issues.", html)
        @test !occursin("href=\"#\"", html)
        @test occursin("href=\"/workbench?project_path=", html)
        @test occursin("data-capability=\"geometry_editor\"", html)
        @test occursin("aria-disabled=\"true\"", html)
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

        rm(output_file)
        missing_output_health = OWENS.studio_project_health(project)
        @test missing_output_health["status"] == "attention"
        @test missing_output_health["runs"][1]["status"] == "ok"
        @test missing_output_health["runs"][1]["run_manifest_health"]["status"] ==
              "attention"
        missing_output_remediation = missing_output_health["runs"][1]["run_manifest_health"]["outputs"][1]["remediation"]
        @test missing_output_remediation["schema_version"] ==
              "owens-run-artifact-remediation/v1"
        @test missing_output_remediation["code"] == "missing_run_artifact"
        @test missing_output_remediation["field"] == "outputs[1].path"
        missing_output_html = OWENS.render_studio_workbench_html(missing_output_health)
        @test occursin("missing_run_artifact", missing_output_html)
        @test occursin("outputs[1].path", missing_output_html)
        @test occursin("Re-run the case", missing_output_html)

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
    home_html = OWENS.render_studio_home_html()
    @test occursin("<title>OWENS Studio</title>", home_html)
    @test occursin(
        "meta name=\"owens-studio-style\" content=\"owens-studio-shared-style/v1\"",
        home_html,
    )
    @test occursin("grid-template-columns: 220px minmax(0, 1fr);", home_html)
    @test occursin("Project Gallery", home_html)
    @test occursin("Example Projects", home_html)
    @test occursin("RM2 GUI Fixture", home_html)
    @test occursin(joinpath("examples", "gui", "rm2", "owens_project.yml"), home_html)
    @test occursin("New Project Templates", home_html)
    @test occursin("Blank OWENS Studio Project", home_html)
    @test occursin("RM2 VAWT Template", home_html)
    @test !occursin("href=\"#\"", home_html)
    @test occursin("href=\"/\"", home_html)
    @test occursin("data-capability=\"project_workbench\"", home_html)
    @test occursin("Project workflow is unavailable in this context.", home_html)

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

@testset "OWENS Studio capability and input summaries" begin
    @test OWENS.studio_gui_capability_ids()[1:5] == [
        "project_manifest_health",
        "template_catalog",
        "generated_script_export",
        "input_file_summary",
        "route_contracts",
    ]
    capabilities = OWENS.studio_gui_capability_catalog()
    @test capabilities["schema_version"] == "owens-studio-capability-catalog/v1"
    @test capabilities["capabilities"][4]["id"] == "input_file_summary"
    @test capabilities["capabilities"][4]["status"] == "implemented"
    @test "OWENS_APP" in capabilities["capabilities"][4]["surfaces"]
    hawt_capability_index =
        findfirst(
            row -> row["id"] == "hawt_aeroelastic_workflow",
            capabilities["capabilities"],
        )
    @test hawt_capability_index !== nothing
    hawt_capability = capabilities["capabilities"][hawt_capability_index]
    @test hawt_capability["status"] == "experimental"
    @test occursin("validation-gated", hawt_capability["description"])
    @test "OWENSOpenFASTWrappers" in hawt_capability["surfaces"]
    @test OWENS.studio_gui_quality_gate_ids()[1:4] == [
        "service_contract_tests",
        "file_provenance",
        "generated_script_equivalence",
        "physical_validation",
    ]
    quality_gates = OWENS.studio_gui_quality_gate_catalog()
    @test quality_gates["schema_version"] == "owens-studio-quality-gates/v1"
    @test quality_gates["summary"]["gates"] == length(quality_gates["gates"])
    @test quality_gates["summary"]["required_for_done"] == length(quality_gates["gates"])
    @test quality_gates["gates"][1]["status"] == "implemented"
    @test "structured error tests" in quality_gates["gates"][1]["evidence_required"]
    browser_gate_index =
        findfirst(row -> row["id"] == "browser_interaction_tests", quality_gates["gates"])
    @test browser_gate_index !== nothing
    @test quality_gates["gates"][browser_gate_index]["status"] == "planned"

    mktempdir() do dir
        model_file = joinpath(dir, "modeling_options.yml")
        windio_file = joinpath(dir, "design.yml")
        project_file = joinpath(dir, "owens_project.yml")
        write(model_file, "OWENS_Options:\n  numTS: 2\n  delta_t: 0.1\n")
        write(windio_file, "name: editor-unit\nassembly:\n  turbine_class: VAWT\n")
        OWENS.write_studio_project(
            project_file,
            dir;
            project_id = "editor-unit",
            name = "Editor Unit",
            modeling_options_file = model_file,
            windio_file,
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )

        summary = OWENS.studio_project_input_summary(project_file; include_text = true)
        @test summary["schema_version"] == "owens-studio-input-summary/v1"
        @test summary["project_id"] == "editor-unit"
        @test summary["capability_gates"] == OrderedCollections.OrderedDict{String,Any}[]
        @test summary["summary"] == OrderedCollections.OrderedDict{String,Any}(
            "files" => 2,
            "editable" => 2,
            "parse_errors" => 0,
            "text_included" => 2,
        )
        @test [row["role"] for row in summary["files"]] == ["modeling_options", "windio"]
        @test [row["format"] for row in summary["files"]] == ["yaml", "yaml"]
        @test all(row["editable"] === true for row in summary["files"])
        @test summary["files"][1]["top_level_keys"] == ["OWENS_Options"]
        @test summary["files"][2]["top_level_keys"] == ["name", "assembly"]
        @test summary["files"][1]["validation_status"] == "ok"
        @test summary["files"][1]["validation_blocking"] === false
        @test isempty(summary["files"][1]["validation_issues"])
        @test occursin("numTS: 2", summary["files"][1]["text"])
        session = OWENS.studio_project_session_summary(project_file)
        @test session["schema_version"] == "owens-studio-session/v1"
        @test session["session_state"] == "clean"
        @test session["dirty"] === false
        @test session["reload_required"] === false
        @test session["summary"]["save_conflicts"] == 0
        @test length(session["file_states"]) == 2
        @test all(row["needs_reload"] === false for row in session["file_states"])

        conflict_model_file = joinpath(dir, "conflict_modeling_options.yml")
        conflict_windio_file = joinpath(dir, "conflict_windio.yml")
        conflict_project_file = joinpath(dir, "conflict_project.yml")
        write(conflict_model_file, "OWENS_Options:\n  numTS: 2\n")
        write(conflict_windio_file, "name: conflict\nassembly:\n  turbine_class: VAWT\n")
        OWENS.write_studio_project(
            conflict_project_file,
            dir;
            project_id = "conflict-unit",
            name = "Conflict Unit",
            modeling_options_file = conflict_model_file,
            windio_file = conflict_windio_file,
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        conflict_clean = OWENS.studio_project_session_summary(conflict_project_file)
        @test conflict_clean["session_state"] == "clean"
        write(conflict_model_file, "OWENS_Options:\n  numTS: 99\n")
        conflict_dirty = OWENS.studio_project_session_summary(conflict_project_file)
        @test conflict_dirty["session_state"] == "needs_reload"
        @test conflict_dirty["dirty"] === true
        @test conflict_dirty["reload_required"] === true
        @test conflict_dirty["summary"]["external_changes"] == 1
        @test conflict_dirty["summary"]["save_conflicts"] == 1
        @test conflict_dirty["save_conflicts"][1]["role"] == "modeling_options"
        @test conflict_dirty["save_conflicts"][1]["needs_reload"] === true
        editor_html = OWENS.render_studio_project_editor_html(summary)
        @test occursin("<title>Editor Unit - OWENS Studio Editor</title>", editor_html)
        @test occursin(
            "meta name=\"owens-studio-style\" content=\"owens-studio-shared-style/v1\"",
            editor_html,
        )
        @test occursin(
            "grid-template-columns: 220px minmax(0, 1fr) 320px;",
            editor_html,
        )
        @test occursin("name=\"expected_sha256\"", editor_html)
        @test occursin(
            "<label for=\"studio-input-modeling_options\">modeling_options input text</label>",
            editor_html,
        )
        @test occursin("id=\"studio-input-modeling_options-help\"", editor_html)
        @test occursin(
            "aria-describedby=\"studio-input-modeling_options-help\"",
            editor_html,
        )
        @test occursin("aria-invalid=\"false\"", editor_html)
        @test occursin("textarea:focus", editor_html)
        @test occursin("name=\"text\"", editor_html)
        @test occursin("OWENS_Options", editor_html)
        @test !occursin("href=\"#\"", editor_html)
        @test occursin("href=\"/workbench?project_path=", editor_html)
        @test occursin("data-capability=\"run_workflow\"", editor_html)
        @test occursin("Simulation<span class=\"nav-status\">planned</span>", editor_html)
        editor_file = joinpath(dir, "editor.html")
        written_editor = OWENS.write_studio_project_editor_html(editor_file, summary)
        @test written_editor == read(editor_file, String)
        @test occursin("action=\"/api/project/input\"", written_editor)

        truncated = OWENS.studio_project_input_summary(
            project_file;
            include_text = true,
            max_text_bytes = 1,
        )
        @test truncated["summary"]["text_included"] == 0
        @test all(row["text_truncated"] === true for row in truncated["files"])
        @test all(row["parse_status"] == "skipped_size_limit" for row in truncated["files"])
        @test all(row["validation_status"] == "not_validated" for row in truncated["files"])
        @test all(row["validation_blocking"] === false for row in truncated["files"])
        @test all(length(row["validation_issues"]) == 1 for row in truncated["files"])
        @test truncated["files"][1]["validation_issues"][1]["remediation_action"] ==
              "increase_max_text_bytes_or_validate_externally"
        @test occursin(
            "max_text_bytes=1",
            truncated["files"][1]["validation_issues"][1]["message"],
        )
        saved_project = OWENS.save_studio_project_input_text(
            project_file,
            "modeling_options",
            "OWENS_Options:\n  numTS: 4\n  delta_t: 0.05\n";
            expected_sha256 = summary["files"][1]["actual_sha256"],
            updated_at_utc = "2026-05-21T00:00:00.000Z",
        )
        @test saved_project["updated_at_utc"] == "2026-05-21T00:00:00.000Z"
        @test saved_project["files"][1]["sha256"] == OWENS.file_sha256(model_file)
        @test saved_project["files"][1]["sha256"] != summary["files"][1]["actual_sha256"]
        @test saved_project["last_input_save"]["schema_version"] ==
              "owens-studio-input-save-provenance/v1"
        save_info = saved_project["last_input_save"]
        @test save_info["role"] == "modeling_options"
        @test save_info["atomic_write"] === true
        @test save_info["validation"]["schema_version"] ==
              "owens-studio-input-validation/v1"
        @test save_info["validation"]["status"] == "ok"
        @test save_info["validation"]["blocking"] === false
        @test isempty(save_info["validation"]["issues"])
        @test save_info["before"]["sha256"] == summary["files"][1]["actual_sha256"]
        @test save_info["after"]["sha256"] == OWENS.file_sha256(model_file)
        @test save_info["backup"]["sha256"] == summary["files"][1]["actual_sha256"]
        @test occursin(".owens-studio-history", save_info["backup"]["path"])
        backup_path = joinpath(dir, save_info["backup"]["path"])
        @test isfile(backup_path)
        @test occursin("numTS: 2", read(backup_path, String))
        saved_summary =
            OWENS.studio_project_input_summary(project_file; include_text = true)
        @test saved_summary["summary"]["parse_errors"] == 0
        @test occursin("numTS: 4", saved_summary["files"][1]["text"])
        history_dir = joinpath(dir, ".owens-studio-history")
        history_count = length(readdir(history_dir))
        current_sha = OWENS.file_sha256(model_file)
        @test_throws ArgumentError OWENS.save_studio_project_input_text(
            project_file,
            "modeling_options",
            "OWENS_Options:\n  bad: [\n";
            expected_sha256 = current_sha,
        )
        @test OWENS.file_sha256(model_file) == current_sha
        @test length(readdir(history_dir)) == history_count
        parse_validation = OWENS.validate_studio_project_input_text(
            "modeling_options",
            model_file,
            "OWENS_Options:\n  bad: [\n",
        )
        @test parse_validation["status"] == "error"
        @test parse_validation["blocking"] === true
        parse_issue = parse_validation["issues"][1]
        @test parse_issue["schema_version"] == "owens-studio-validation-issue/v1"
        @test parse_issue["severity"] == "parse_error"
        @test parse_issue["role"] == "modeling_options"
        @test parse_issue["path"] == model_file
        @test parse_issue["field"] == "document"
        @test parse_issue["yaml_path"] == "document"
        @test parse_issue["remediation_action"] == "fix_yaml_syntax"
        @test occursin("Fix the YAML syntax", parse_issue["suggested_fix"])
        schema_validation = OWENS.validate_studio_project_input_text(
            "modeling_options",
            model_file,
            "Not_OWENS_Options:\n  numTS: 5\n",
        )
        @test schema_validation["status"] == "error"
        @test schema_validation["blocking"] === true
        schema_issue = schema_validation["issues"][1]
        @test schema_issue["schema_version"] == "owens-studio-validation-issue/v1"
        @test schema_issue["severity"] == "schema_error"
        @test schema_issue["role"] == "modeling_options"
        @test schema_issue["field"] == "OWENS_Options"
        @test schema_issue["yaml_path"] == "OWENS_Options"
        @test schema_issue["documentation"] == "docs/src/getting_started.md"
        @test schema_issue["remediation_action"] == "add_owens_options_mapping"
        @test occursin("Add a top-level OWENS_Options", schema_issue["suggested_fix"])
        @test occursin("time stepping", schema_issue["physical_implication"])
        @test_throws ArgumentError OWENS.save_studio_project_input_text(
            project_file,
            "modeling_options",
            "OWENS_Options:\n  numTS: 5\n";
            expected_sha256 = summary["files"][1]["actual_sha256"],
        )
        @test_throws ArgumentError OWENS.studio_project_input_summary(
            project_file;
            max_text_bytes = -1,
        )

        hawt_windio_file = joinpath(dir, "hawt_design.yml")
        hawt_project_file = joinpath(dir, "hawt_project.yml")
        write(hawt_windio_file, "name: hawt-editor-unit\nassembly:\n  turbine_class: HAWT\n")
        OWENS.write_studio_project(
            hawt_project_file,
            dir;
            project_id = "hawt-editor-unit",
            name = "HAWT Editor Unit",
            modeling_options_file = model_file,
            windio_file = hawt_windio_file,
            created_at_utc = "2026-05-20T00:00:00.000Z",
        )
        hawt_summary =
            OWENS.studio_project_input_summary(hawt_project_file; include_text = true)
        @test length(hawt_summary["capability_gates"]) == 1
        @test hawt_summary["files"][2]["validation_status"] == "warning"
        @test hawt_summary["files"][2]["validation_blocking"] === false
        @test length(hawt_summary["files"][2]["validation_issues"]) == 1
        @test hawt_summary["files"][2]["validation_issues"][1]["yaml_path"] ==
              "assembly.turbine_class"
        hawt_gate = hawt_summary["capability_gates"][1]
        @test hawt_gate["capability"] == "hawt_aeroelastic_workflow"
        @test hawt_gate["status"] == "experimental"
        @test hawt_gate["severity"] == "warning"
        @test hawt_gate["source_role"] == "windio"
        @test hawt_gate["detected_value"] == "HAWT"
        @test occursin("validation-gated", hawt_gate["message"])
        hawt_validation = OWENS.validate_studio_project_input_text(
            "windio",
            hawt_windio_file,
            read(hawt_windio_file, String),
        )
        @test hawt_validation["status"] == "warning"
        @test hawt_validation["blocking"] === false
        hawt_issue = hawt_validation["issues"][1]
        @test hawt_issue["schema_version"] == "owens-studio-validation-issue/v1"
        @test hawt_issue["severity"] == "unsupported_feature"
        @test hawt_issue["role"] == "windio"
        @test hawt_issue["field"] == "assembly.turbine_class"
        @test hawt_issue["yaml_path"] == "assembly.turbine_class"
        @test hawt_issue["remediation_action"] == "acknowledge_experimental_hawt"
        @test occursin("validation workflow", hawt_issue["suggested_fix"])
        hawt_editor_html = OWENS.render_studio_project_editor_html(hawt_summary)
        @test occursin("warning validation", hawt_editor_html)
        @test occursin("data-validation-severity=\"unsupported_feature\"", hawt_editor_html)
        @test occursin("data-validation-field=\"assembly.turbine_class\"", hawt_editor_html)
        @test occursin("id=\"studio-input-windio-issues\"", hawt_editor_html)
        @test occursin(
            "aria-describedby=\"studio-input-windio-help studio-input-windio-issues\"",
            hawt_editor_html,
        )
        @test occursin("Suggested fix", hawt_editor_html)
        hawt_open_payload = OWENS_APP.open_studio_project(hawt_project_file)
        @test hawt_open_payload["inputs"]["capability_gates"] ==
              hawt_summary["capability_gates"]
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
    @test _portable_paths([row["path"] for row in health["files"]]) == ["../../RM2/modeling_options_OWENS_RM2.yml", "../../RM2/WINDIO_RM2.yaml"]
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
        capability_catalog = OWENS_APP.list_studio_gui_capabilities()
        @test capability_catalog["schema_version"] == "owens-studio-capability-catalog/v1"
        @test capability_catalog["capabilities"][1]["id"] == "project_manifest_health"
        @test capability_catalog["capabilities"][4]["id"] == "input_file_summary"
        quality_gate_catalog = OWENS_APP.list_studio_quality_gates()
        @test quality_gate_catalog["schema_version"] == "owens-studio-quality-gates/v1"
        @test quality_gate_catalog["summary"]["partial"] >= 1
        @test quality_gate_catalog["gates"][2]["id"] == "file_provenance"
        doctor = OWENS_APP.studio_doctor(; output_dir = dir)
        @test doctor["schema_version"] == "owens-studio-doctor/v1"
        @test doctor["status"] == "ok"
        @test doctor["output_dir"] == abspath(dir)
        @test doctor["summary"]["checks"] == length(doctor["checks"])
        @test doctor["summary"]["failed"] == 0
        @test doctor["summary"]["routes"] > 0
        @test doctor["quality_gates"]["schema_version"] == "owens-studio-quality-gates/v1"
        @test any(
            row -> row["name"] == "output_dir_writable" && row["passed"] === true,
            doctor["checks"],
        )
        @test occursin("studio-home", doctor["commands"]["home_html"])
        example_catalog = OWENS_APP.list_studio_example_projects()
        @test example_catalog["schema_version"] == "owens-studio-example-catalog/v1"
        @test [row["example"] for row in example_catalog["examples"]] == ["rm2"]
        route_catalog = OWENS_APP.studio_route_catalog()
        @test route_catalog["schema_version"] == "owens-studio-route-catalog/v1"
        @test [row["name"] for row in route_catalog["routes"]] == [
            "route_catalog",
            "studio_home",
            "capability_catalog",
            "quality_gate_catalog",
            "template_catalog",
            "example_catalog",
            "project_open",
            "project_inputs",
            "project_editor",
            "project_input_save",
            "project_session",
            "project_health",
            "project_workbench",
            "project_script",
            "project_bundle",
            "create_template_project",
            "studio_doctor",
        ]
        @test [row["method"] for row in route_catalog["routes"]] == [
            "GET",
            "GET",
            "GET",
            "GET",
            "GET",
            "GET",
            "GET",
            "GET",
            "GET",
            "POST",
            "GET",
            "GET",
            "GET",
            "GET",
            "POST",
            "POST",
            "GET",
        ]
        @test route_catalog["routes"][2]["path"] == "/"
        @test route_catalog["contract"]["schema_version"] ==
              "owens-studio-route-contract/v1"
        @test route_catalog["contract"]["request_schema_version"] ==
              "owens-studio-route-request/v1"
        @test route_catalog["contract"]["response_schema_version"] ==
              "owens-studio-route-response/v1"
        @test route_catalog["contract"]["error_schema_version"] == "owens-studio-error/v1"
        @test route_catalog["contract"]["json_supported"] === false
        @test occursin("YAML-first", route_catalog["contract"]["json_status_reason"])
        @test route_catalog["contract"]["compatibility"]["route_order_is_not_contractual"] ===
              true
        @test route_catalog["contract"]["security"]["workspace_root_param"] ==
              "workspace_root"
        @test route_catalog["contract"]["security"]["workspace_root_enforced_path_params"] ==
              ["project_path", "output_dir", "target"]
        @test route_catalog["contract"]["security"]["workspace_root_rejects_allow_external"] ===
              true
        @test route_catalog["contract"]["security"]["local_only_server_binding_required_when_http_is_added"] ===
              true
        @test route_catalog["contract"]["security"]["csrf_or_non_browser_token_required_when_http_is_added"] ===
              true
        @test all(
            row["request_schema"]["schema_version"] == "owens-studio-route-request/v1" for
            row in route_catalog["routes"]
        )
        @test all(
            row["response_schema"]["schema_version"] == "owens-studio-route-response/v1" for
            row in route_catalog["routes"]
        )
        @test all(
            row["response_schema"]["error_schema_version"] == "owens-studio-error/v1" for
            row in route_catalog["routes"]
        )
        route_by_name = OrderedCollections.OrderedDict(
            row["name"] => row for row in route_catalog["routes"]
        )
        @test route_catalog["routes"][7]["required_params"] == ["project_path"]
        @test route_catalog["routes"][7]["optional_params"] == ["summarize_runs"]
        @test route_by_name["project_open"]["request_schema"]["params"]["required"][1] ==
              OrderedCollections.OrderedDict{String,Any}(
                  "name" => "project_path",
                  "type" => "string",
                  "required" => true,
                  "description" => "Path to an OWENS Studio project manifest.",
              )
        @test route_catalog["routes"][8]["required_params"] == ["project_path"]
        @test route_catalog["routes"][8]["optional_params"] ==
              ["include_text", "max_text_bytes"]
        @test route_by_name["project_inputs"]["request_schema"]["params"]["optional"][1]["name"] ==
              "include_text"
        @test route_by_name["project_inputs"]["request_schema"]["params"]["optional"][1]["type"] ==
              "boolean"
        @test route_catalog["routes"][9]["required_params"] == ["project_path"]
        @test route_catalog["routes"][9]["content_type"] == "text/html; charset=utf-8"
        @test route_catalog["routes"][10]["required_params"] ==
              ["project_path", "role", "text"]
        @test route_catalog["routes"][10]["optional_params"] ==
              ["expected_sha256", "allow_external", "updated_at_utc"]
        @test route_catalog["routes"][11]["required_params"] == ["project_path"]
        @test route_catalog["routes"][11]["optional_params"] ==
              ["include_text", "max_text_bytes"]
        @test route_catalog["routes"][13]["content_type"] == "text/html; charset=utf-8"
        @test route_catalog["routes"][15]["required_params"] ==
              ["project_path", "output_dir"]
        @test route_by_name["quality_gate_catalog"]["response_schema"]["payload_schema_version"] ==
              "owens-studio-quality-gates/v1"
        @test route_by_name["project_inputs"]["response_schema"]["payload_schema_version"] ==
              "owens-studio-input-summary/v1"
        @test route_by_name["project_session"]["response_schema"]["payload_schema_version"] ==
              "owens-studio-session/v1"
        @test route_by_name["project_workbench"]["response_schema"]["body_kind"] ==
              "html_document"
        @test route_by_name["project_workbench"]["response_schema"]["payload_schema_version"] ===
              nothing
        @test route_by_name["project_script"]["response_schema"]["body_kind"] == "plain_text"
        @test route_by_name["create_template_project"]["response_schema"]["payload_schema_version"] ==
              "owens-studio-template-create/v1"
        @test route_by_name["studio_doctor"]["optional_params"] == ["output_dir"]
        @test route_by_name["studio_doctor"]["response_schema"]["payload_schema_version"] ==
              "owens-studio-doctor/v1"

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
        project_inputs =
            OWENS_APP.inspect_studio_project_inputs(project_file; include_text = true)
        @test project_inputs["schema_version"] == "owens-studio-input-summary/v1"
        @test project_inputs["summary"]["editable"] == 2
        @test project_inputs["files"][1]["top_level_keys"] == ["OWENS_Options"]
        @test occursin("numTS: 2", project_inputs["files"][1]["text"])
        project_session = OWENS_APP.inspect_studio_project_session(project_file)
        @test project_session["schema_version"] == "owens-studio-session/v1"
        @test project_session["session_state"] == "clean"
        @test project_session["summary"]["save_conflicts"] == 0
        html = OWENS_APP.write_studio_project_workbench(html_file, project_health)
        @test html == read(html_file, String)
        @test occursin("Studio App", html)

        manifest_cli =
            OWENS_APP.real_main(["manifest-health", manifest_file]; io = IOBuffer())
        @test manifest_cli["status"] == "ok"
        saved_input = OWENS_APP.save_studio_project_input(
            project_file,
            "modeling_options",
            "OWENS_Options:\n  numTS: 3\n";
            expected_sha256 = project_inputs["files"][1]["actual_sha256"],
            updated_at_utc = "2026-05-21T00:00:00.000Z",
        )
        @test saved_input["schema_version"] == "owens-studio-input-save/v1"
        @test saved_input["project_status"] == "attention"
        @test saved_input["project"]["updated_at_utc"] == "2026-05-21T00:00:00.000Z"
        @test saved_input["project"]["files"][1]["sha256"] == OWENS.file_sha256(model_file)
        @test saved_input["save"]["atomic_write"] === true
        @test saved_input["save"]["before"]["sha256"] ==
              project_inputs["files"][1]["actual_sha256"]
        @test saved_input["save"]["after"]["sha256"] == OWENS.file_sha256(model_file)
        @test saved_input["save"]["backup"]["sha256"] ==
              project_inputs["files"][1]["actual_sha256"]
        @test saved_input["save"]["validation"]["status"] == "ok"
        @test saved_input["health"]["runs"][1]["run_manifest_health"]["status"] ==
              "attention"
        @test saved_input["health"]["runs"][1]["run_manifest_health"]["inputs"][1]["remediation"]["code"] ==
              "modified_run_artifact"
        @test occursin("numTS: 3", saved_input["inputs"]["files"][1]["text"])
        home_cli = OWENS_APP.real_main(["studio-home", html_file]; io = IOBuffer())
        @test home_cli["output_html"] == abspath(html_file)
        @test home_cli["bytes"] == stat(html_file).size
        @test occursin("Project Gallery", read(html_file, String))
        doctor_cli = OWENS_APP.real_main(["studio-doctor", dir]; io = IOBuffer())
        @test doctor_cli["schema_version"] == "owens-studio-doctor/v1"
        @test doctor_cli["status"] == "ok"
        @test doctor_cli["output_dir"] == abspath(dir)
        capabilities_cli = OWENS_APP.real_main(["project-capabilities"]; io = IOBuffer())
        @test capabilities_cli["schema_version"] == "owens-studio-capability-catalog/v1"
        quality_gates_cli =
            OWENS_APP.real_main(["project-quality-gates"]; io = IOBuffer())
        @test quality_gates_cli["schema_version"] == "owens-studio-quality-gates/v1"
        @test quality_gates_cli["summary"]["gates"] == length(quality_gates_cli["gates"])
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
        @test template_cli["schema_version"] == "owens-studio-template-create/v1"
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
        @test [row["route"] for row in open_payload["actions"]] == [
            "project_health",
            "project_inputs",
            "project_editor",
            "project_input_save",
            "project_workbench",
            "project_script",
            "project_bundle",
            "capability_catalog",
        ]
        @test open_payload["capabilities"]["schema_version"] ==
              "owens-studio-capability-catalog/v1"
        @test open_payload["inputs"]["schema_version"] == "owens-studio-input-summary/v1"
        @test open_payload["inputs"]["summary"]["editable"] == 2
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
        inputs_cli = OWENS_APP.real_main(
            ["project-inputs", template_cli["project_file"]];
            io = IOBuffer(),
        )
        @test inputs_cli["schema_version"] == "owens-studio-input-summary/v1"
        @test inputs_cli["summary"]["editable"] == 2
        limited_inputs_cli = OWENS_APP.real_main(
            ["project-inputs", template_cli["project_file"], "true", "1"];
            io = IOBuffer(),
        )
        @test limited_inputs_cli["summary"]["text_included"] == 0
        @test all(
            row["parse_status"] == "skipped_size_limit" for
            row in limited_inputs_cli["files"]
        )
        @test limited_inputs_cli["files"][1]["validation_issues"][1]["remediation_action"] ==
              "increase_max_text_bytes_or_validate_externally"
        session_cli = OWENS_APP.real_main(
            ["project-session", template_cli["project_file"]];
            io = IOBuffer(),
        )
        @test session_cli["schema_version"] == "owens-studio-session/v1"
        @test session_cli["session_state"] == "clean"
        @test session_cli["summary"]["save_conflicts"] == 0
        limited_session_cli = OWENS_APP.real_main(
            ["project-session", template_cli["project_file"], "true", "1"];
            io = IOBuffer(),
        )
        @test limited_session_cli["inputs"]["summary"]["text_included"] == 0
        @test all(
            row["parse_status"] == "skipped_size_limit" for
            row in limited_session_cli["inputs"]["files"]
        )
        editor_cli = OWENS_APP.real_main(
            [
                "project-editor-html",
                template_cli["project_file"],
                joinpath(dir, "template-editor.html"),
            ];
            io = IOBuffer(),
        )
        @test editor_cli["output_html"] == abspath(joinpath(dir, "template-editor.html"))
        @test editor_cli["bytes"] == stat(editor_cli["output_html"]).size
        @test occursin("OWENS Studio Editor", read(editor_cli["output_html"], String))
        save_text_file = joinpath(dir, "updated_modeling_options.yml")
        write(save_text_file, "OWENS_Options:\n  numTS: 6\n")
        save_cli = OWENS_APP.real_main(
            [
                "project-save-input",
                template_cli["project_file"],
                "modeling_options",
                save_text_file,
            ];
            io = IOBuffer(),
        )
        @test save_cli["schema_version"] == "owens-studio-input-save/v1"
        @test save_cli["project_status"] == "attention"
        @test save_cli["health"]["runs"][1]["run_manifest_health"]["inputs"][1]["remediation"]["code"] ==
              "modified_run_artifact"
        @test occursin("numTS: 6", save_cli["inputs"]["files"][1]["text"])
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
        @test bundle_cli["project_status"] == "attention"
        @test isfile(bundle_cli["index_html"])
        @test isfile(bundle_cli["editor_html"])
        @test isfile(bundle_cli["health_file"])
        @test isfile(bundle_cli["inputs_file"])
        @test isfile(bundle_cli["script_file"])
        @test isfile(bundle_cli["open_file"])
        @test bundle_cli["bytes"]["index_html"] == stat(bundle_cli["index_html"]).size
        @test bundle_cli["bytes"]["editor_html"] == stat(bundle_cli["editor_html"]).size
        @test bundle_cli["bytes"]["inputs_file"] == stat(bundle_cli["inputs_file"]).size
        @test bundle_cli["bytes"]["open_file"] == stat(bundle_cli["open_file"]).size
        bundle_index = read(bundle_cli["index_html"], String)
        @test occursin("Open Payload", bundle_index)
        @test occursin("Editor", bundle_index)
        @test occursin("open.yml", bundle_index)
        bundle_editor = read(bundle_cli["editor_html"], String)
        @test occursin("OWENS Studio Editor", bundle_editor)
        @test occursin("project_path", bundle_editor)
        bundle_open = YAML.load_file(
            bundle_cli["open_file"];
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_open["schema_version"] == "owens-studio-open/v1"
        @test bundle_open["generated_script"]["available"] === true
        project_cli = OWENS_APP.real_main(["project-health", project_file]; io = IOBuffer())
        @test project_cli["status"] == "attention"
        project_html_cli =
            OWENS_APP.real_main(["project-html", project_file, html_file]; io = IOBuffer())
        @test project_html_cli["output_html"] == abspath(html_file)
        @test project_html_cli["project_status"] == "attention"

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
        @test routes_payload["routes"][end - 1]["handler"] ==
              "studio_project_template_route"
        @test routes_payload["routes"][end]["handler"] == "studio_doctor_route"

        home_response = OWENS_APP.studio_home_route()
        @test home_response.status == 200
        @test home_response.content_type == "text/html; charset=utf-8"
        @test occursin("Project Gallery", home_response.body)
        @test occursin("RM2 GUI Fixture", home_response.body)
        dispatch_home = OWENS_APP.dispatch_studio_route("/"; method = "GET")
        @test dispatch_home.status == 200
        @test occursin("New Project Templates", dispatch_home.body)

        doctor_response = OWENS_APP.studio_doctor_route(output_dir = dir)
        @test doctor_response.status == 200
        @test doctor_response.content_type == "application/x-yaml; charset=utf-8"
        doctor_payload = YAML.load(
            doctor_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test doctor_payload["schema_version"] == "owens-studio-doctor/v1"
        @test doctor_payload["status"] == "ok"
        @test doctor_payload["output_dir"] == abspath(dir)
        dispatch_doctor = OWENS_APP.dispatch_studio_route(
            "/api/doctor";
            method = "GET",
            params = Dict("output_dir" => dir, "workspace_root" => dir),
        )
        @test dispatch_doctor.status == 200
        @test occursin("owens-studio-doctor/v1", dispatch_doctor.body)

        capabilities_response = OWENS_APP.studio_gui_capabilities_route()
        @test capabilities_response.status == 200
        @test capabilities_response.content_type == "application/x-yaml; charset=utf-8"
        capabilities_payload = YAML.load(
            capabilities_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test capabilities_payload["schema_version"] == "owens-studio-capability-catalog/v1"
        @test capabilities_payload["capabilities"][4]["id"] == "input_file_summary"
        dispatch_capabilities =
            OWENS_APP.dispatch_studio_route("/api/capabilities"; method = "GET")
        @test dispatch_capabilities.status == 200
        @test occursin("owens-studio-capability-catalog/v1", dispatch_capabilities.body)

        quality_gates_response = OWENS_APP.studio_quality_gates_route()
        @test quality_gates_response.status == 200
        @test quality_gates_response.content_type == "application/x-yaml; charset=utf-8"
        quality_gates_payload = YAML.load(
            quality_gates_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test quality_gates_payload["schema_version"] == "owens-studio-quality-gates/v1"
        @test quality_gates_payload["summary"]["gates"] ==
              length(quality_gates_payload["gates"])
        dispatch_quality_gates =
            OWENS_APP.dispatch_studio_route("/api/quality-gates"; method = "GET")
        @test dispatch_quality_gates.status == 200
        @test occursin("owens-studio-quality-gates/v1", dispatch_quality_gates.body)

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
        @test open_payload["actions"][6]["route"] == "project_script"
        @test open_payload["actions"][6]["enabled"] === true
        @test open_payload["inputs"]["summary"]["editable"] == 2
        @test open_payload["session"]["schema_version"] == "owens-studio-session/v1"
        @test open_payload["session"]["session_state"] == "clean"

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
        @test dispatch_open_payload["routes"]["routes"][7]["name"] == "project_open"
        sandboxed_open = OWENS_APP.dispatch_studio_route(
            "project_open";
            method = "GET",
            params = Dict("project_path" => project_file, "workspace_root" => dir),
        )
        @test sandboxed_open.status == 200
        sandboxed_open_payload = YAML.load(
            sandboxed_open.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test sandboxed_open_payload["project_file"] == project_file
        @test sandboxed_open_payload["project_status"] == "ok"

        inputs_response =
            OWENS_APP.studio_project_inputs_route(project_file; include_text = true)
        @test inputs_response.status == 200
        @test inputs_response.content_type == "application/x-yaml; charset=utf-8"
        inputs_payload = YAML.load(
            inputs_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test inputs_payload["schema_version"] == "owens-studio-input-summary/v1"
        @test inputs_payload["summary"]["editable"] == 2
        @test inputs_payload["summary"]["text_included"] == 2
        @test !occursin("runOWENSWINDIO", inputs_payload["files"][1]["text"])
        route_model_file = inputs_payload["files"][1]["resolved_path"]
        session_response =
            OWENS_APP.studio_project_session_route(project_file; include_text = true)
        @test session_response.status == 200
        @test session_response.content_type == "application/x-yaml; charset=utf-8"
        session_payload = YAML.load(
            session_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test session_payload["schema_version"] == "owens-studio-session/v1"
        @test session_payload["session_state"] == "clean"
        @test session_payload["summary"]["save_conflicts"] == 0
        @test session_payload["inputs"]["summary"]["text_included"] == 2

        dispatch_session = OWENS_APP.dispatch_studio_route(
            "/api/project/session";
            method = "GET",
            params = Dict("project_path" => project_file, "workspace_root" => dir),
        )
        @test dispatch_session.status == 200
        @test occursin("owens-studio-session/v1", dispatch_session.body)

        dispatch_inputs = OWENS_APP.dispatch_studio_route(
            "project_inputs";
            method = "GET",
            params = Dict(
                "project_path" => project_file,
                "include_text" => true,
                "max_text_bytes" => 1,
            ),
        )
        @test dispatch_inputs.status == 200
        @test occursin("text_truncated: true", dispatch_inputs.body)

        editor_response = OWENS_APP.studio_project_editor_route(project_file)
        @test editor_response.status == 200
        @test editor_response.content_type == "text/html; charset=utf-8"
        @test occursin("OWENS Studio Editor", editor_response.body)
        @test occursin("name=\"expected_sha256\"", editor_response.body)
        @test occursin("modeling_options", editor_response.body)
        dispatch_editor = OWENS_APP.dispatch_studio_route(
            "/project/edit";
            method = "GET",
            params = Dict("project_path" => project_file),
        )
        @test dispatch_editor.status == 200
        @test occursin("name=\"text\"", dispatch_editor.body)

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

        save_response = OWENS_APP.studio_project_input_save_route(
            project_file,
            "modeling_options",
            "OWENS_Options:\n  numTS: 7\n";
            expected_sha256 = inputs_payload["files"][1]["actual_sha256"],
            updated_at_utc = "2026-05-21T00:00:00.000Z",
        )
        @test save_response.status == 200
        @test save_response.content_type == "application/x-yaml; charset=utf-8"
        save_payload = YAML.load(
            save_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test save_payload["schema_version"] == "owens-studio-input-save/v1"
        @test save_payload["project_status"] == "attention"
        @test save_payload["project"]["updated_at_utc"] == "2026-05-21T00:00:00.000Z"
        @test save_payload["save"]["atomic_write"] === true
        @test save_payload["save"]["before"]["sha256"] ==
              inputs_payload["files"][1]["actual_sha256"]
        @test save_payload["save"]["after"]["sha256"] == OWENS.file_sha256(
            joinpath(dirname(project_file), save_payload["save"]["after"]["path"]),
        )
        @test save_payload["save"]["validation"]["status"] == "ok"
        @test save_payload["health"]["runs"][1]["run_manifest_health"]["inputs"][1]["remediation"]["code"] ==
              "modified_run_artifact"
        @test occursin("numTS: 7", save_payload["inputs"]["files"][1]["text"])

        stale_save = OWENS_APP.dispatch_studio_route(
            "project_input_save";
            method = "POST",
            params = Dict(
                "project_path" => project_file,
                "role" => "modeling_options",
                "text" => "OWENS_Options:\n  numTS: 8\n",
                "expected_sha256" => inputs_payload["files"][1]["actual_sha256"],
            ),
        )
        @test stale_save.status == 400
        @test occursin("Refusing to overwrite modified Studio input", stale_save.body)
        current_route_sha = OWENS.file_sha256(route_model_file)
        invalid_save = OWENS_APP.dispatch_studio_route(
            "project_input_save";
            method = "POST",
            params = Dict(
                "project_path" => project_file,
                "role" => "modeling_options",
                "text" => "OWENS_Options:\n  bad: [\n",
                "expected_sha256" => current_route_sha,
            ),
        )
        @test invalid_save.status == 400
        @test occursin("Refusing to save invalid Studio input", invalid_save.body)
        @test OWENS.file_sha256(route_model_file) == current_route_sha
        sandbox_external_save = OWENS_APP.dispatch_studio_route(
            "project_input_save";
            method = "POST",
            params = Dict(
                "project_path" => project_file,
                "role" => "modeling_options",
                "text" => "OWENS_Options:\n  numTS: 9\n",
                "expected_sha256" => current_route_sha,
                "allow_external" => true,
                "workspace_root" => dir,
            ),
        )
        @test sandbox_external_save.status == 403
        sandbox_external_payload = YAML.load(
            sandbox_external_save.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test sandbox_external_payload["schema_version"] == "owens-studio-error/v1"
        @test sandbox_external_payload["code"] == "workspace_boundary_violation"
        @test sandbox_external_payload["route"] == "project_input_save"
        @test occursin("allow_external=true", sandbox_external_payload["message"])

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
        @test bundle_payload["project_status"] == "attention"
        @test isfile(bundle_payload["index_html"])
        @test isfile(bundle_payload["editor_html"])
        @test isfile(bundle_payload["health_file"])
        @test isfile(bundle_payload["inputs_file"])
        @test isfile(bundle_payload["script_file"])
        @test isfile(bundle_payload["open_file"])
        @test bundle_payload["bytes"]["index_html"] ==
              stat(bundle_payload["index_html"]).size
        @test bundle_payload["bytes"]["editor_html"] ==
              stat(bundle_payload["editor_html"]).size
        @test bundle_payload["bytes"]["open_file"] == stat(bundle_payload["open_file"]).size
        bundle_payload_index = read(bundle_payload["index_html"], String)
        @test occursin("Open Payload", bundle_payload_index)
        @test occursin("editor.html", bundle_payload_index)
        @test occursin("open.yml", bundle_payload_index)
        @test occursin("name=\"role\"", read(bundle_payload["editor_html"], String))
        bundle_open_payload = YAML.load_file(
            bundle_payload["open_file"];
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test bundle_open_payload["project_file"] == project_file
        @test bundle_open_payload["actions"][7]["route"] == "project_bundle"

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
        @test template_payload["schema_version"] == "owens-studio-template-create/v1"
        @test template_payload["template"] == "blank"
        @test template_payload["project_status"] == "ok"
        @test template_payload["run_manifest_file"] === nothing
        @test isfile(template_payload["project_file"])
        mktempdir(dirname(dir)) do outside_dir
            outside_project = OWENS.create_studio_project_template(
                joinpath(outside_dir, "outside-rm2");
                template = "rm2",
                created_at_utc = "2026-05-20T00:00:00.000Z",
            )
            outside_open = OWENS_APP.dispatch_studio_route(
                "project_open";
                method = "GET",
                params = Dict(
                    "project_path" => outside_project["project_file"],
                    "workspace_root" => dir,
                ),
            )
            @test outside_open.status == 403
            outside_open_payload = YAML.load(
                outside_open.body;
                dicttype = OrderedCollections.OrderedDict{String,Any},
            )
            @test outside_open_payload["schema_version"] == "owens-studio-error/v1"
            @test outside_open_payload["code"] == "workspace_boundary_violation"
            @test outside_open_payload["route"] == "project_open"
            @test occursin("project_path", outside_open_payload["message"])

            outside_template = OWENS_APP.dispatch_studio_route(
                "create_template_project";
                method = "POST",
                params = Dict(
                    "target" => joinpath(outside_dir, "blocked-template"),
                    "template" => "blank",
                    "workspace_root" => dir,
                ),
            )
            @test outside_template.status == 403
            outside_template_payload = YAML.load(
                outside_template.body;
                dicttype = OrderedCollections.OrderedDict{String,Any},
            )
            @test outside_template_payload["code"] == "workspace_boundary_violation"
            @test outside_template_payload["route"] == "create_template_project"
            @test occursin("target", outside_template_payload["message"])
        end

        error_response = OWENS_APP.studio_project_health_route(joinpath(dir, "missing.yml"))
        @test error_response.status == 400
        @test error_response.content_type == "application/x-yaml; charset=utf-8"
        error_payload = YAML.load(
            error_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test error_payload["status"] == "error"
        @test error_payload["schema_version"] == "owens-studio-error/v1"
        @test error_payload["code"] == "invalid_request"
        @test error_payload["severity"] == "error"
        @test error_payload["exception_type"] == "ArgumentError"
        @test error_payload["route"] === nothing
        @test occursin("project path", error_payload["suggested_fix"])
        @test error_payload["developer_detail"] == error_payload["message"]
        @test occursin("Cannot read missing Studio project", error_payload["message"])
        editor_error =
            OWENS_APP.studio_project_editor_route(joinpath(dir, "missing-editor.yml"))
        @test editor_error.status == 400
        @test editor_error.content_type == "text/html; charset=utf-8"
        @test occursin("OWENS Studio Error", editor_error.body)
        @test occursin("Suggested fix", editor_error.body)
        @test occursin("Developer detail", editor_error.body)
        @test occursin("data-error-schema=\"owens-studio-error/v1\"", editor_error.body)
        @test occursin("data-error-code=\"invalid_request\"", editor_error.body)
        @test occursin("data-error-route=\"project_editor\"", editor_error.body)
        @test occursin("Cannot read missing Studio project", editor_error.body)
        workbench_error = OWENS_APP.dispatch_studio_route(
            "/workbench";
            method = "GET",
            params = Dict("project_path" => joinpath(dir, "missing-workbench.yml")),
        )
        @test workbench_error.status == 400
        @test workbench_error.content_type == "text/html; charset=utf-8"
        @test occursin("data-error-route=\"project_workbench\"", workbench_error.body)
        @test occursin("Suggested fix", workbench_error.body)
        missing_route_response =
            OWENS_APP.dispatch_studio_route("/api/not-a-route"; method = "GET")
        @test missing_route_response.status == 404
        missing_route_payload = YAML.load(
            missing_route_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test missing_route_payload["schema_version"] == "owens-studio-error/v1"
        @test missing_route_payload["code"] == "route_not_found"
        @test missing_route_payload["route"] == "/api/not-a-route"
        @test occursin("/api/routes", missing_route_payload["suggested_fix"])
        wrong_method_response =
            OWENS_APP.dispatch_studio_route("project_open"; method = "POST")
        @test wrong_method_response.status == 405
        wrong_method_payload = YAML.load(
            wrong_method_response.body;
            dicttype = OrderedCollections.OrderedDict{String,Any},
        )
        @test wrong_method_payload["schema_version"] == "owens-studio-error/v1"
        @test wrong_method_payload["code"] == "method_not_allowed"
        @test wrong_method_payload["route"] == "project_open"
        @test occursin("Use GET", wrong_method_payload["suggested_fix"])
    end
end
