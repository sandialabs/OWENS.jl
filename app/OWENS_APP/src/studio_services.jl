export inspect_run_manifest,
    inspect_output_data,
    prepare_windio_run,
    list_studio_project_templates,
    list_studio_example_projects,
    create_studio_template_project,
    open_studio_project,
    inspect_studio_project,
    inspect_studio_project_script,
    write_studio_project_bundle,
    write_studio_project_workbench

const STUDIO_OPEN_SCHEMA_VERSION = "owens-studio-open/v1"

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
        _studio_open_action("project_workbench", "Render workbench HTML", true),
        _studio_open_action("project_script", "View generated Julia", script_available),
        _studio_open_action("project_bundle", "Write static workbench bundle", true),
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
    )
    open(bundle["index_html"], "w") do io
        write(io, html)
    end

    bundle["open_file"] = open_file
    bundle["bytes"]["index_html"] = stat(bundle["index_html"]).size
    bundle["bytes"]["open_file"] = stat(open_file).size
    return bundle
end

function write_studio_project_workbench(output_html::AbstractString, project_or_health;)
    return OWENS.write_studio_workbench_html(output_html, project_or_health)
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
