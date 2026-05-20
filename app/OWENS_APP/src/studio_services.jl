export inspect_run_manifest,
    inspect_output_data,
    prepare_windio_run,
    create_studio_template_project,
    inspect_studio_project,
    write_studio_project_workbench

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

function inspect_studio_project(
    path::AbstractString;
    root = nothing,
    summarize_runs::Bool = true,
)
    return OWENS.studio_project_health(path; root, summarize_runs)
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
