export inspect_run_manifest,
    inspect_output_data,
    prepare_windio_run,
    list_studio_gui_capabilities,
    list_studio_quality_gates,
    list_studio_project_templates,
    list_studio_example_projects,
    create_studio_template_project,
    open_studio_project,
    inspect_studio_project,
    inspect_studio_project_inputs,
    inspect_studio_project_session,
    studio_doctor,
    save_studio_project_input,
    inspect_studio_project_script,
    write_studio_home,
    write_studio_project_bundle,
    write_studio_project_editor,
    write_studio_project_workbench

const STUDIO_OPEN_SCHEMA_VERSION = "owens-studio-open/v1"
const STUDIO_INPUT_SAVE_SCHEMA_VERSION = "owens-studio-input-save/v1"
const STUDIO_TEMPLATE_CREATE_SCHEMA_VERSION = "owens-studio-template-create/v1"
const STUDIO_DOCTOR_SCHEMA_VERSION = "owens-studio-doctor/v1"

function inspect_run_manifest(
    path::AbstractString;
    root = nothing,
    summarize_outputs::Bool = true,
)
    return OWENS.run_manifest_health(path; root, summarize_outputs)
end

function inspect_output_data(
    path::AbstractString;
    channels = OWENS.output_data_channel_names(),
    include_unregistered::Bool = false,
)
    return OrderedCollections.OrderedDict{String,Any}(
        "path" => abspath(path),
        "channels" => [
            _service_value(row) for
            row in OWENS.output_data_summary(path; channels, include_unregistered)
        ],
    )
end

function prepare_windio_run(
    modeling_options_file::AbstractString,
    windio_file::AbstractString,
    run_path::AbstractString;
    create_run_path::Bool = false,
    build_manifest::Bool = true,
    manifest_kwargs...,
)
    spec =
        OWENS.windio_run_spec(modeling_options_file, windio_file, run_path; create_run_path)
    script = OWENS.render_windio_run_script(spec)
    manifest =
        build_manifest ? OWENS.build_windio_run_manifest(spec; manifest_kwargs...) : nothing
    return OrderedCollections.OrderedDict{String,Any}(
        "spec" => OrderedCollections.OrderedDict{String,Any}(
            "modeling_options_file" => spec.modeling_options_file,
            "windio_file" => spec.windio_file,
            "run_path" => spec.run_path,
        ),
        "script" => script,
        "manifest" => manifest,
    )
end

function list_studio_gui_capabilities()
    return OWENS.studio_gui_capability_catalog()
end

function list_studio_quality_gates()
    return OWENS.studio_gui_quality_gate_catalog()
end

function list_studio_project_templates()
    return OWENS.studio_project_template_catalog()
end

function list_studio_example_projects()
    return OWENS.studio_example_project_catalog()
end

function create_studio_template_project(
    target::AbstractString;
    template::AbstractString = "rm2",
    overwrite::Bool = false,
    created_at_utc = nothing,
)
    created =
        OWENS.create_studio_project_template(target; template, overwrite, created_at_utc)
    health = OWENS.studio_project_health(created["project_file"])
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_TEMPLATE_CREATE_SCHEMA_VERSION,
        "template" => created["template"],
        "project_file" => created["project_file"],
        "run_manifest_file" => created["run_manifest_file"],
        "script_file" => created["script_file"],
        "project_status" => health["status"],
        "project_health" => health,
    )
end

function open_studio_project(path::AbstractString; summarize_runs::Bool = true)
    health = inspect_studio_project(path; summarize_runs)
    script = _studio_script_artifact(path, health["root"])
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_OPEN_SCHEMA_VERSION,
        "project_file" => abspath(path),
        "project_status" => health["status"],
        "project_id" => health["project_id"],
        "name" => health["name"],
        "root" => health["root"],
        "metadata" => health["metadata"],
        "generated_script" => script,
        "actions" => _studio_open_actions(script["available"]),
        "capabilities" => list_studio_gui_capabilities(),
        "inputs" => inspect_studio_project_inputs(path),
        "session" => inspect_studio_project_session(path),
        "templates" => list_studio_project_templates(),
        "examples" => list_studio_example_projects(),
        "routes" => studio_route_catalog(),
        "health" => health,
    )
end

function inspect_studio_project(
    path::AbstractString;
    root = nothing,
    summarize_runs::Bool = true,
)
    return OWENS.studio_project_health(path; root, summarize_runs)
end

function inspect_studio_project_inputs(
    path::AbstractString;
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    return OWENS.studio_project_input_summary(path; include_text, max_text_bytes)
end

function inspect_studio_project_session(
    path::AbstractString;
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    return OWENS.studio_project_session_summary(path; include_text, max_text_bytes)
end

function studio_doctor(; output_dir::AbstractString = pwd())
    checks = OrderedCollections.OrderedDict{String,Any}[]
    templates = list_studio_project_templates()
    examples = list_studio_example_projects()
    capabilities = list_studio_gui_capabilities()
    routes = studio_route_catalog()
    output_dir_abs = abspath(output_dir)

    _push_studio_doctor_check!(
        checks,
        "julia_version",
        true,
        "Julia $(VERSION) is running with $(Threads.nthreads()) thread(s).",
        severity = "info",
    )
    _push_studio_doctor_check!(
        checks,
        "owens_loaded",
        true,
        "OWENS is loaded and Studio service functions are callable.",
    )
    _push_studio_doctor_check!(
        checks,
        "yaml_loaded",
        isdefined(YAML, :write),
        "YAML serialization is available for Studio command payloads.",
    )
    _push_studio_doctor_check!(
        checks,
        "hdf5_available",
        !isnothing(Base.find_package("HDF5")),
        "HDF5 is visible to Julia for result summary workflows.",
        severity = "warning",
        suggested_fix = "Instantiate the OWENS project environment if HDF5 result summaries are needed.",
    )
    _push_studio_doctor_check!(
        checks,
        "template_catalog",
        "rm2" in [row["template"] for row in templates["templates"]],
        "Built-in Studio templates are available.",
    )
    _push_studio_doctor_check!(
        checks,
        "example_catalog",
        all(row -> row["available"] === true, examples["examples"]),
        "Committed Studio example projects are available.",
        severity = "warning",
        suggested_fix = "Check the examples/gui directory if example projects are missing.",
    )
    _push_studio_doctor_check!(
        checks,
        "route_catalog",
        length(routes["routes"]) > 0,
        "Dependency-light Studio route contracts are available.",
    )
    _push_studio_doctor_check!(
        checks,
        "output_dir_writable",
        _studio_output_dir_writable(output_dir_abs),
        "Diagnostic output directory is writable: $output_dir_abs",
        suggested_fix = "Choose a writable output directory for generated Studio HTML and bundles.",
    )

    has_error = any(checks) do row
        row["passed"] === false && row["severity"] == "error"
    end
    has_failure = any(checks) do row
        row["passed"] === false
    end
    status = has_error ? "error" : has_failure ? "attention" : "ok"

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_DOCTOR_SCHEMA_VERSION,
        "status" => status,
        "output_dir" => output_dir_abs,
        "julia" => OrderedCollections.OrderedDict{String,Any}(
            "version" => string(VERSION),
            "threads" => Threads.nthreads(),
            "project" => Base.active_project(),
        ),
        "commands" => OrderedCollections.OrderedDict{String,Any}(
            "diagnose" => "studio-doctor $(output_dir_abs)",
            "home_html" => "studio-home $(joinpath(output_dir_abs, "owens_studio_home.html"))",
            "create_rm2_template" =>
                "project-template rm2 $(joinpath(output_dir_abs, "owens-rm2-studio"))",
        ),
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "checks" => length(checks),
            "passed" => count(row -> row["passed"] === true, checks),
            "failed" => count(row -> row["passed"] === false, checks),
            "templates" => length(templates["templates"]),
            "examples" => length(examples["examples"]),
            "capabilities" => length(capabilities["capabilities"]),
            "routes" => length(routes["routes"]),
        ),
        "checks" => checks,
        "templates" => templates,
        "examples" => examples,
        "capabilities" => capabilities,
        "quality_gates" => list_studio_quality_gates(),
        "routes" => routes,
    )
end

function save_studio_project_input(
    project_path::AbstractString,
    role::AbstractString,
    text::AbstractString;
    expected_sha256 = nothing,
    allow_external::Bool = false,
    updated_at_utc = nothing,
)
    project = OWENS.save_studio_project_input_text(
        project_path,
        role,
        text;
        expected_sha256,
        allow_external,
        updated_at_utc,
    )
    health = inspect_studio_project(project_path)
    inputs = inspect_studio_project_inputs(project_path; include_text = true)
    save = get(project, "last_input_save", nothing)
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_INPUT_SAVE_SCHEMA_VERSION,
        "project_file" => abspath(project_path),
        "role" => string(role),
        "save" => save,
        "project_status" => health["status"],
        "project" => project,
        "health" => health,
        "inputs" => inputs,
    )
end

function _studio_script_artifact(project_path::AbstractString, root::AbstractString)
    script_path = OWENS.studio_project_generated_script_path(project_path)
    if isnothing(script_path)
        return OrderedCollections.OrderedDict{String,Any}(
            "path" => nothing,
            "relative_path" => nothing,
            "available" => false,
            "bytes" => nothing,
            "sha256" => nothing,
        )
    end

    available = isfile(script_path)
    return OrderedCollections.OrderedDict{String,Any}(
        "path" => script_path,
        "relative_path" => relpath(script_path, root),
        "available" => available,
        "bytes" => available ? stat(script_path).size : nothing,
        "sha256" => available ? OWENS.file_sha256(script_path) : nothing,
    )
end

function _studio_open_actions(script_available::Bool)
    return OrderedCollections.OrderedDict{String,Any}[
        _studio_open_action("project_health", "Inspect project health", true),
        _studio_open_action("project_inputs", "Inspect editable inputs", true),
        _studio_open_action("project_editor", "Edit project inputs", true),
        _studio_open_action("project_input_save", "Save editable input", true),
        _studio_open_action("project_workbench", "Render workbench HTML", true),
        _studio_open_action("project_script", "View generated Julia", script_available),
        _studio_open_action("project_bundle", "Write static workbench bundle", true),
        _studio_open_action("capability_catalog", "Inspect GUI capabilities", true),
    ]
end

function _studio_open_action(route::AbstractString, label::AbstractString, enabled::Bool)
    return OrderedCollections.OrderedDict{String,Any}(
        "route" => string(route),
        "label" => string(label),
        "enabled" => enabled,
    )
end

function inspect_studio_project_script(path::AbstractString)
    script_path = OWENS.studio_project_generated_script_path(path)
    script = OWENS.read_studio_project_generated_script(path)
    return OrderedCollections.OrderedDict{String,Any}(
        "project_file" => abspath(path),
        "script_file" => script_path,
        "script" => script,
    )
end

function write_studio_home(output_html::AbstractString)
    return OWENS.write_studio_home_html(output_html)
end

function write_studio_project_bundle(
    output_dir::AbstractString,
    project_path::AbstractString;
    include_script::Bool = true,
)
    bundle = OWENS.write_studio_workbench_bundle(output_dir, project_path; include_script)
    open_file = joinpath(bundle["bundle_dir"], "open.yml")
    YAML.write_file(open_file, open_studio_project(project_path))
    html = OWENS.render_studio_workbench_html(
        inspect_studio_project(project_path);
        health_href = basename(bundle["health_file"]),
        script_href = isnothing(bundle["script_file"]) ? nothing :
                      basename(bundle["script_file"]),
        open_href = basename(open_file),
        editor_href = basename(bundle["editor_html"]),
    )
    open(bundle["index_html"], "w") do io
        write(io, html)
    end

    bundle["open_file"] = open_file
    bundle["bytes"]["index_html"] = stat(bundle["index_html"]).size
    bundle["bytes"]["open_file"] = stat(open_file).size
    return bundle
end

function write_studio_project_editor(
    output_html::AbstractString,
    project_path::AbstractString,
)
    inputs = inspect_studio_project_inputs(project_path; include_text = true)
    return OWENS.write_studio_project_editor_html(output_html, inputs)
end

function write_studio_project_workbench(output_html::AbstractString, project_or_health;)
    return OWENS.write_studio_workbench_html(output_html, project_or_health)
end

function _push_studio_doctor_check!(
    checks::Vector{OrderedCollections.OrderedDict{String,Any}},
    name::AbstractString,
    passed::Bool,
    message::AbstractString;
    severity::AbstractString = passed ? "info" : "error",
    suggested_fix::AbstractString = "",
)
    push!(
        checks,
        OrderedCollections.OrderedDict{String,Any}(
            "name" => string(name),
            "passed" => passed,
            "severity" => string(severity),
            "message" => string(message),
            "suggested_fix" => passed ? "" : string(suggested_fix),
        ),
    )
    return checks
end

function _studio_output_dir_writable(output_dir::AbstractString)
    try
        mkpath(output_dir)
        probe = tempname(output_dir; cleanup = false)
        open(probe, "w") do io
            write(io, "owens-studio-doctor\n")
        end
        isfile(probe) || return false
        rm(probe; force = true)
        return true
    catch
        return false
    end
end

function _service_value(value)
    if value === missing
        return nothing
    elseif value isa NamedTuple
        return OrderedCollections.OrderedDict{String,Any}(
            string(key) => _service_value(getfield(value, key)) for key in keys(value)
        )
    elseif value isa AbstractDict
        normalized = OrderedCollections.OrderedDict{String,Any}()
        keys_iter =
            value isa OrderedCollections.OrderedDict ? collect(keys(value)) :
            sort(collect(keys(value)); by = string)
        for key in keys_iter
            normalized[string(key)] = _service_value(value[key])
        end
        return normalized
    elseif value isa AbstractVector || value isa Tuple
        return [_service_value(item) for item in value]
    elseif value isa Symbol
        return string(value)
    elseif value isa VersionNumber
        return string(value)
    elseif value isa AbstractString ||
           value isa Number ||
           value isa Bool ||
           isnothing(value)
        return value
    else
        return string(value)
    end
end
