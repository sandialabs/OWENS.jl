export StudioRouteResponse,
    studio_project_health_route,
    studio_project_workbench_route,
    studio_project_script_route,
    studio_project_template_route

struct StudioRouteResponse
    status::Int
    content_type::String
    body::String
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
