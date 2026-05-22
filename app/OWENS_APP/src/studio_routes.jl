export studio_route_catalog,
    StudioRouteResponse,
    dispatch_studio_route,
    studio_routes_route,
    studio_home_route,
    studio_project_templates_route,
    studio_project_examples_route,
    studio_project_open_route,
    studio_project_health_route,
    studio_project_workbench_route,
    studio_project_script_route,
    studio_project_bundle_route,
    studio_project_template_route

struct StudioRouteResponse
    status::Int
    content_type::String
    body::String
end

const STUDIO_ROUTE_CATALOG_SCHEMA_VERSION = "owens-studio-route-catalog/v1"

"""
    studio_route_catalog()

Return the dependency-light route map that a Genie shell should wrap.
"""
function studio_route_catalog()
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_CATALOG_SCHEMA_VERSION,
        "routes" => OrderedCollections.OrderedDict{String,Any}[
            _studio_route_record(
                "route_catalog",
                "GET",
                "/api/routes",
                "studio_routes_route",
                "application/x-yaml; charset=utf-8",
                "List OWENS Studio route contracts.",
            ),
            _studio_route_record(
                "studio_home",
                "GET",
                "/",
                "studio_home_route",
                "text/html; charset=utf-8",
                "Render the Studio project chooser.",
            ),
            _studio_route_record(
                "template_catalog",
                "GET",
                "/api/templates",
                "studio_project_templates_route",
                "application/x-yaml; charset=utf-8",
                "List built-in project templates.",
            ),
            _studio_route_record(
                "example_catalog",
                "GET",
                "/api/examples",
                "studio_project_examples_route",
                "application/x-yaml; charset=utf-8",
                "List committed Studio example projects.",
            ),
            _studio_route_record(
                "project_open",
                "GET",
                "/api/project/open",
                "studio_project_open_route",
                "application/x-yaml; charset=utf-8",
                "Open a project and return workbench bootstrap data.",
                required_params = ["project_path"],
                optional_params = ["summarize_runs"],
            ),
            _studio_route_record(
                "project_health",
                "GET",
                "/api/project/health",
                "studio_project_health_route",
                "application/x-yaml; charset=utf-8",
                "Inspect Studio project health.",
                required_params = ["project_path"],
                optional_params = ["summarize_runs"],
            ),
            _studio_route_record(
                "project_workbench",
                "GET",
                "/workbench",
                "studio_project_workbench_route",
                "text/html; charset=utf-8",
                "Render the Studio workbench HTML.",
                required_params = ["project_path"],
            ),
            _studio_route_record(
                "project_script",
                "GET",
                "/api/project/script",
                "studio_project_script_route",
                "text/plain; charset=utf-8",
                "Return the generated Julia driver.",
                required_params = ["project_path"],
            ),
            _studio_route_record(
                "project_bundle",
                "POST",
                "/api/project/bundle",
                "studio_project_bundle_route",
                "application/x-yaml; charset=utf-8",
                "Write a static workbench bundle.",
                required_params = ["project_path", "output_dir"],
                optional_params = ["include_script"],
            ),
            _studio_route_record(
                "create_template_project",
                "POST",
                "/api/project/template",
                "studio_project_template_route",
                "application/x-yaml; charset=utf-8",
                "Create a project from a built-in template.",
                required_params = ["target"],
                optional_params = ["template", "overwrite", "created_at_utc"],
            ),
        ],
    )
end

"""
    dispatch_studio_route(route; method=nothing, params=Dict())

Resolve a Studio route by catalog name or path and call the matching
dependency-light route handler. This keeps future Genie/HTTP glue thin: the web
layer should translate query/body data into `params` and leave OWENS project
semantics here.
"""
function dispatch_studio_route(
    route::AbstractString;
    method = nothing,
    params = OrderedCollections.OrderedDict{String,Any}(),
)
    record = _studio_route_record_for(route)
    if isnothing(record)
        return _studio_route_error_response(
            ArgumentError("Unknown Studio route: $route");
            status = 404,
        )
    end

    requested_method = isnothing(method) ? record["method"] : uppercase(string(method))
    if requested_method != record["method"]
        return _studio_route_error_response(
            ArgumentError(
                "Method $requested_method is not allowed for Studio route $(record["name"])",
            );
            status = 405,
        )
    end

    try
        return _dispatch_studio_route_name(record["name"], _studio_route_params(params))
    catch err
        return _studio_route_error_response(err)
    end
end

"""
    studio_routes_route()

Return the route catalog as a route response.
"""
function studio_routes_route()
    return _studio_yaml_route_response() do
        studio_route_catalog()
    end
end

"""
    studio_home_route()

Return the static OWENS Studio project chooser HTML.
"""
function studio_home_route()
    try
        return StudioRouteResponse(
            200,
            "text/html; charset=utf-8",
            OWENS.render_studio_home_html(),
        )
    catch err
        return _studio_route_error_response(err)
    end
end

"""
    studio_project_open_route(project_path; summarize_runs=true)

Return the workbench bootstrap payload for a Studio project.
"""
function studio_project_open_route(
    project_path::AbstractString;
    summarize_runs::Bool = true,
)
    return _studio_yaml_route_response() do
        open_studio_project(project_path; summarize_runs)
    end
end

"""
    studio_project_templates_route()

Return the built-in OWENS Studio project template catalog as a route response.
"""
function studio_project_templates_route()
    return _studio_yaml_route_response() do
        list_studio_project_templates()
    end
end

"""
    studio_project_examples_route()

Return committed OWENS Studio example projects as a route response.
"""
function studio_project_examples_route()
    return _studio_yaml_route_response() do
        list_studio_example_projects()
    end
end

"""
    studio_project_health_route(project_path; summarize_runs=true)

Return a dependency-light route response for the project health API. Genie or
another web shell can wrap this response without owning OWENS project semantics.
"""
function studio_project_health_route(
    project_path::AbstractString;
    summarize_runs::Bool = true,
)
    return _studio_yaml_route_response() do
        inspect_studio_project(project_path; summarize_runs)
    end
end

"""
    studio_project_workbench_route(project_path)

Return a route response containing the static OWENS Studio workbench HTML.
"""
function studio_project_workbench_route(project_path::AbstractString)
    try
        html = OWENS.render_studio_workbench_html(project_path)
        return StudioRouteResponse(200, "text/html; charset=utf-8", html)
    catch err
        return _studio_route_error_response(err)
    end
end

"""
    studio_project_script_route(project_path)

Return a route response containing the generated Julia script referenced by a
Studio project.
"""
function studio_project_script_route(project_path::AbstractString)
    try
        script = OWENS.read_studio_project_generated_script(project_path)
        return StudioRouteResponse(200, "text/plain; charset=utf-8", script)
    catch err
        return _studio_route_error_response(err)
    end
end

"""
    studio_project_bundle_route(project_path, output_dir)

Write the static workbench bundle and return its file manifest as a route
response.
"""
function studio_project_bundle_route(
    project_path::AbstractString,
    output_dir::AbstractString;
    include_script::Bool = true,
)
    return _studio_yaml_route_response() do
        write_studio_project_bundle(output_dir, project_path; include_script)
    end
end

"""
    studio_project_template_route(target; template="rm2", overwrite=false)

Create a Studio project template and return the resulting health payload as a
route response.
"""
function studio_project_template_route(
    target::AbstractString;
    template::AbstractString = "rm2",
    overwrite::Bool = false,
    created_at_utc = nothing,
)
    return _studio_yaml_route_response() do
        create_studio_template_project(target; template, overwrite, created_at_utc)
    end
end

function _studio_yaml_route_response(build_payload::Function)
    try
        return StudioRouteResponse(
            200,
            "application/x-yaml; charset=utf-8",
            _studio_yaml_body(build_payload()),
        )
    catch err
        return _studio_route_error_response(err)
    end
end

function _studio_route_error_response(err; status::Integer = 400)
    payload = OrderedCollections.OrderedDict{String,Any}(
        "status" => "error",
        "message" => sprint(showerror, err),
    )
    return StudioRouteResponse(
        Int(status),
        "application/x-yaml; charset=utf-8",
        _studio_yaml_body(payload),
    )
end

function _studio_yaml_body(payload)
    io = IOBuffer()
    YAML.write(io, payload)
    write(io, "\n")
    return String(take!(io))
end

function _studio_route_record(
    name::AbstractString,
    method::AbstractString,
    path::AbstractString,
    handler::AbstractString,
    content_type::AbstractString,
    description::AbstractString,
    ;
    required_params = String[],
    optional_params = String[],
)
    return OrderedCollections.OrderedDict{String,Any}(
        "name" => string(name),
        "method" => string(method),
        "path" => string(path),
        "handler" => string(handler),
        "content_type" => string(content_type),
        "description" => string(description),
        "required_params" => string.(required_params),
        "optional_params" => string.(optional_params),
    )
end

function _studio_route_record_for(route::AbstractString)
    route_id = string(route)
    for record in studio_route_catalog()["routes"]
        if route_id == record["name"] || route_id == record["path"]
            return record
        end
    end
    return nothing
end

function _dispatch_studio_route_name(name::AbstractString, params::AbstractDict)
    if name == "route_catalog"
        return studio_routes_route()
    elseif name == "studio_home"
        return studio_home_route()
    elseif name == "template_catalog"
        return studio_project_templates_route()
    elseif name == "example_catalog"
        return studio_project_examples_route()
    elseif name == "project_open"
        return studio_project_open_route(
            _studio_route_required(params, "project_path");
            summarize_runs = _studio_route_optional(params, "summarize_runs", true),
        )
    elseif name == "project_health"
        return studio_project_health_route(
            _studio_route_required(params, "project_path");
            summarize_runs = _studio_route_optional(params, "summarize_runs", true),
        )
    elseif name == "project_workbench"
        return studio_project_workbench_route(
            _studio_route_required(params, "project_path"),
        )
    elseif name == "project_script"
        return studio_project_script_route(_studio_route_required(params, "project_path"))
    elseif name == "project_bundle"
        return studio_project_bundle_route(
            _studio_route_required(params, "project_path"),
            _studio_route_required(params, "output_dir");
            include_script = _studio_route_optional(params, "include_script", true),
        )
    elseif name == "create_template_project"
        return studio_project_template_route(
            _studio_route_required(params, "target");
            template = _studio_route_optional(params, "template", "rm2"),
            overwrite = _studio_route_optional(params, "overwrite", false),
            created_at_utc = _studio_route_optional(params, "created_at_utc", nothing),
        )
    end

    throw(ArgumentError("Unhandled Studio route: $name"))
end

function _studio_route_params(params)
    if isnothing(params)
        return OrderedCollections.OrderedDict{String,Any}()
    elseif params isa AbstractDict || params isa NamedTuple
        normalized = OrderedCollections.OrderedDict{String,Any}()
        for (key, value) in pairs(params)
            normalized[string(key)] = value
        end
        return normalized
    end

    throw(ArgumentError("Studio route params must be a dictionary or named tuple"))
end

function _studio_route_required(params::AbstractDict, key::AbstractString)
    haskey(params, key) || throw(ArgumentError("Missing Studio route parameter: $key"))
    return params[key]
end

function _studio_route_optional(params::AbstractDict, key::AbstractString, default)
    return haskey(params, key) ? params[key] : default
end
