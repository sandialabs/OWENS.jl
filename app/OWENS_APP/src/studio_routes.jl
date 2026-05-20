export studio_route_catalog,
    StudioRouteResponse,
    studio_routes_route,
    studio_project_templates_route,
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
                "template_catalog",
                "GET",
                "/api/templates",
                "studio_project_templates_route",
                "application/x-yaml; charset=utf-8",
                "List built-in project templates.",
            ),
            _studio_route_record(
                "project_health",
                "GET",
                "/api/project/health",
                "studio_project_health_route",
                "application/x-yaml; charset=utf-8",
                "Inspect Studio project health.",
            ),
            _studio_route_record(
                "project_workbench",
                "GET",
                "/workbench",
                "studio_project_workbench_route",
                "text/html; charset=utf-8",
                "Render the Studio workbench HTML.",
            ),
            _studio_route_record(
                "project_script",
                "GET",
                "/api/project/script",
                "studio_project_script_route",
                "text/plain; charset=utf-8",
                "Return the generated Julia driver.",
            ),
            _studio_route_record(
                "project_bundle",
                "POST",
                "/api/project/bundle",
                "studio_project_bundle_route",
                "application/x-yaml; charset=utf-8",
                "Write a static workbench bundle.",
            ),
            _studio_route_record(
                "create_template_project",
                "POST",
                "/api/project/template",
                "studio_project_template_route",
                "application/x-yaml; charset=utf-8",
                "Create a project from a built-in template.",
            ),
        ],
    )
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
    studio_project_templates_route()

Return the built-in OWENS Studio project template catalog as a route response.
"""
function studio_project_templates_route()
    return _studio_yaml_route_response() do
        list_studio_project_templates()
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

function _studio_route_error_response(err)
    payload = OrderedCollections.OrderedDict{String,Any}(
        "status" => "error",
        "message" => sprint(showerror, err),
    )
    return StudioRouteResponse(
        400,
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
)
    return OrderedCollections.OrderedDict{String,Any}(
        "name" => string(name),
        "method" => string(method),
        "path" => string(path),
        "handler" => string(handler),
        "content_type" => string(content_type),
        "description" => string(description),
    )
end
