export studio_route_catalog,
    StudioRouteResponse,
    dispatch_studio_route,
    studio_routes_route,
    studio_home_route,
    studio_doctor_route,
    studio_gui_capabilities_route,
    studio_quality_gates_route,
    studio_project_templates_route,
    studio_project_examples_route,
    studio_project_open_route,
    studio_project_inputs_route,
    studio_project_editor_route,
    studio_project_input_save_route,
    studio_project_session_route,
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

struct StudioRouteWorkspaceBoundaryError <: Exception
    message::String
end

Base.showerror(io::IO, err::StudioRouteWorkspaceBoundaryError) = print(io, err.message)

const STUDIO_ROUTE_CATALOG_SCHEMA_VERSION = "owens-studio-route-catalog/v1"
const STUDIO_ROUTE_CONTRACT_SCHEMA_VERSION = "owens-studio-route-contract/v1"
const STUDIO_ROUTE_REQUEST_SCHEMA_VERSION = "owens-studio-route-request/v1"
const STUDIO_ROUTE_RESPONSE_SCHEMA_VERSION = "owens-studio-route-response/v1"
const STUDIO_ROUTE_ERROR_SCHEMA_VERSION = "owens-studio-error/v1"

const STUDIO_ROUTE_PARAM_SPECS = OrderedCollections.OrderedDict{String,Any}(
    "allow_external" => ("boolean", "Whether project input saves may write files outside the project root."),
    "created_at_utc" => ("string", "Optional ISO-8601 UTC timestamp for deterministic project creation."),
    "expected_sha256" => ("string", "Optional optimistic-lock SHA-256 for the file being saved."),
    "include_script" => ("boolean", "Whether a static bundle should include the generated Julia script."),
    "include_text" => ("boolean", "Whether input summaries should include bounded source text."),
    "max_text_bytes" => ("integer", "Maximum bytes of source text to include per input file."),
    "output_dir" => ("string", "Directory where generated bundle artifacts should be written."),
    "overwrite" => ("boolean", "Whether a project template may replace an existing target directory."),
    "project_path" => ("string", "Path to an OWENS Studio project manifest."),
    "role" => ("string", "Editable Studio input role such as modeling_options or windio."),
    "summarize_runs" => ("boolean", "Whether project health should inspect referenced run manifests."),
    "target" => ("string", "Target directory for a generated Studio project template."),
    "template" => ("string", "Built-in Studio project template identifier."),
    "text" => ("string", "Replacement source text for the editable input file."),
    "updated_at_utc" => ("string", "Optional ISO-8601 UTC timestamp for deterministic input saves."),
    "workspace_root" => ("string", "Optional route-dispatch sandbox root. When supplied, filesystem path parameters must stay inside this root."),
)

const STUDIO_ROUTE_PAYLOAD_SCHEMAS = OrderedCollections.OrderedDict{String,Any}(
    "route_catalog" => STUDIO_ROUTE_CATALOG_SCHEMA_VERSION,
    "studio_home" => nothing,
    "studio_doctor" => STUDIO_DOCTOR_SCHEMA_VERSION,
    "capability_catalog" => "owens-studio-capability-catalog/v1",
    "quality_gate_catalog" => "owens-studio-quality-gates/v1",
    "template_catalog" => "owens-studio-template-catalog/v1",
    "example_catalog" => "owens-studio-example-catalog/v1",
    "project_open" => STUDIO_OPEN_SCHEMA_VERSION,
    "project_inputs" => OWENS.STUDIO_INPUT_SUMMARY_SCHEMA_VERSION,
    "project_editor" => nothing,
    "project_input_save" => STUDIO_INPUT_SAVE_SCHEMA_VERSION,
    "project_session" => OWENS.STUDIO_SESSION_SCHEMA_VERSION,
    "project_health" => OWENS.STUDIO_WORKBENCH_SCHEMA_VERSION,
    "project_workbench" => nothing,
    "project_script" => nothing,
    "project_bundle" => OWENS.STUDIO_WORKBENCH_BUNDLE_SCHEMA_VERSION,
    "create_template_project" => STUDIO_TEMPLATE_CREATE_SCHEMA_VERSION,
)

"""
    studio_route_catalog()

Return the dependency-light route map that a Genie shell should wrap.
"""
function studio_route_catalog()
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_CATALOG_SCHEMA_VERSION,
        "contract" => _studio_route_contract_record(),
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
                "capability_catalog",
                "GET",
                "/api/capabilities",
                "studio_gui_capabilities_route",
                "application/x-yaml; charset=utf-8",
                "List OWENS Studio capability status.",
            ),
            _studio_route_record(
                "quality_gate_catalog",
                "GET",
                "/api/quality-gates",
                "studio_quality_gates_route",
                "application/x-yaml; charset=utf-8",
                "List OWENS Studio quality gates and current evidence.",
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
                "project_inputs",
                "GET",
                "/api/project/inputs",
                "studio_project_inputs_route",
                "application/x-yaml; charset=utf-8",
                "Inspect project input files for editor panes.",
                required_params = ["project_path"],
                optional_params = ["include_text", "max_text_bytes"],
            ),
            _studio_route_record(
                "project_editor",
                "GET",
                "/project/edit",
                "studio_project_editor_route",
                "text/html; charset=utf-8",
                "Render the Studio project input editor HTML.",
                required_params = ["project_path"],
            ),
            _studio_route_record(
                "project_input_save",
                "POST",
                "/api/project/input",
                "studio_project_input_save_route",
                "application/x-yaml; charset=utf-8",
                "Save editable project input text and refresh project provenance.",
                required_params = ["project_path", "role", "text"],
                optional_params = ["expected_sha256", "allow_external", "updated_at_utc"],
            ),
            _studio_route_record(
                "project_session",
                "GET",
                "/api/project/session",
                "studio_project_session_route",
                "application/x-yaml; charset=utf-8",
                "Inspect active project session state and save-conflict status.",
                required_params = ["project_path"],
                optional_params = ["include_text", "max_text_bytes"],
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
            _studio_route_record(
                "studio_doctor",
                "GET",
                "/api/doctor",
                "studio_doctor_route",
                "application/x-yaml; charset=utf-8",
                "Run local OWENS Studio diagnostics.",
                optional_params = ["output_dir"],
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
            code = "route_not_found",
            route,
            suggested_fix = "Use /api/routes or project-routes to inspect available Studio routes.",
        )
    end

    requested_method = isnothing(method) ? record["method"] : uppercase(string(method))
    if requested_method != record["method"]
        return _studio_route_error_response(
            ArgumentError(
                "Method $requested_method is not allowed for Studio route $(record["name"])",
            );
            status = 405,
            code = "method_not_allowed",
            route = record["name"],
            suggested_fix = "Use $(record["method"]) for Studio route $(record["name"]).",
        )
    end

    try
        route_params = _studio_route_params(params)
        _studio_enforce_workspace_root!(record["name"], route_params)
        return _dispatch_studio_route_name(record["name"], route_params)
    catch err
        return _studio_route_error_response(
            err;
            status = _studio_route_error_status(err),
            route = record["name"],
        )
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
        return _studio_html_route_error_response(err; route = "studio_home")
    end
end

"""
    studio_doctor_route(; output_dir=pwd())

Return local OWENS Studio diagnostics as a route response.
"""
function studio_doctor_route(; output_dir::AbstractString = pwd())
    return _studio_yaml_route_response() do
        studio_doctor(; output_dir)
    end
end

"""
    studio_gui_capabilities_route()

Return the GUI capability catalog as a route response.
"""
function studio_gui_capabilities_route()
    return _studio_yaml_route_response() do
        list_studio_gui_capabilities()
    end
end

"""
    studio_quality_gates_route()

Return the Studio quality-gate catalog as a route response.
"""
function studio_quality_gates_route()
    return _studio_yaml_route_response() do
        list_studio_quality_gates()
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
    studio_project_inputs_route(project_path; include_text=false, max_text_bytes=200_000)

Return project input-file summaries for Studio editor panes.
"""
function studio_project_inputs_route(
    project_path::AbstractString;
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    return _studio_yaml_route_response() do
        inspect_studio_project_inputs(project_path; include_text, max_text_bytes)
    end
end

"""
    studio_project_editor_route(project_path)

Return a route response containing the OWENS Studio project input editor HTML.
"""
function studio_project_editor_route(project_path::AbstractString)
    try
        html = OWENS.render_studio_project_editor_html(project_path)
        return StudioRouteResponse(200, "text/html; charset=utf-8", html)
    catch err
        return _studio_html_route_error_response(err; route = "project_editor")
    end
end

"""
    studio_project_input_save_route(project_path, role, text; kwargs...)

Save editable project input text and return refreshed project health.
"""
function studio_project_input_save_route(
    project_path::AbstractString,
    role::AbstractString,
    text::AbstractString;
    expected_sha256 = nothing,
    allow_external::Bool = false,
    updated_at_utc = nothing,
)
    return _studio_yaml_route_response() do
        save_studio_project_input(
            project_path,
            role,
            text;
            expected_sha256,
            allow_external,
            updated_at_utc,
        )
    end
end

"""
    studio_project_session_route(project_path; include_text=false, max_text_bytes=200_000)

Return active project session state for dirty/reload/save-conflict decisions.
"""
function studio_project_session_route(
    project_path::AbstractString;
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    return _studio_yaml_route_response() do
        inspect_studio_project_session(project_path; include_text, max_text_bytes)
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
        return _studio_html_route_error_response(err; route = "project_workbench")
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

function _studio_route_error_response(
    err;
    status::Integer = 400,
    code = nothing,
    route = nothing,
    suggested_fix = nothing,
)
    payload = _studio_route_error_payload(
        err;
        status,
        code,
        route,
        suggested_fix,
    )
    return StudioRouteResponse(
        Int(status),
        "application/x-yaml; charset=utf-8",
        _studio_yaml_body(payload),
    )
end

function _studio_html_route_error_response(
    err;
    status::Integer = _studio_route_error_status(err),
    code = nothing,
    route = nothing,
    suggested_fix = nothing,
)
    payload = _studio_route_error_payload(
        err;
        status,
        code,
        route,
        suggested_fix,
    )
    return StudioRouteResponse(
        Int(status),
        "text/html; charset=utf-8",
        _studio_error_html_body(payload),
    )
end

function _studio_route_error_payload(
    err;
    status::Integer = 400,
    code = nothing,
    route = nothing,
    suggested_fix = nothing,
)
    message = sprint(showerror, err)
    error_code = isnothing(code) ? _studio_route_error_code(err, status) : string(code)
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_ERROR_SCHEMA_VERSION,
        "status" => "error",
        "code" => error_code,
        "severity" => "error",
        "message" => message,
        "route" => isnothing(route) ? nothing : string(route),
        "exception_type" => string(nameof(typeof(err))),
        "suggested_fix" => isnothing(suggested_fix) ?
                           _studio_route_error_suggested_fix(error_code, message) :
                           string(suggested_fix),
        "developer_detail" => message,
    )
end

function _studio_route_error_code(err, status::Integer)
    status == 404 && return "route_not_found"
    status == 405 && return "method_not_allowed"
    err isa StudioRouteWorkspaceBoundaryError && return "workspace_boundary_violation"
    err isa ArgumentError && return "invalid_request"
    err isa KeyError && return "missing_parameter"
    return "route_error"
end

function _studio_route_error_status(err)
    err isa StudioRouteWorkspaceBoundaryError && return 403
    return 400
end

function _studio_route_error_suggested_fix(code::AbstractString, message::AbstractString)
    if code == "missing_parameter"
        return "Provide all required parameters listed in /api/routes."
    elseif code == "invalid_request"
        return "Check the route parameters, project path, and editable input contents."
    elseif code == "workspace_boundary_violation"
        return "Choose a project, target, or output path inside the configured workspace_root, or omit workspace_root only for trusted local scripts."
    elseif code == "route_not_found"
        return "Use /api/routes or project-routes to inspect available Studio routes."
    elseif code == "method_not_allowed"
        return "Use the HTTP method listed for this route in /api/routes."
    end
    return "Inspect developer_detail and retry after correcting the request."
end

function _studio_yaml_body(payload)
    io = IOBuffer()
    YAML.write(io, payload)
    write(io, "\n")
    return String(take!(io))
end

function _studio_error_html_body(payload::AbstractDict)
    schema = _studio_route_html_escape(get(payload, "schema_version", STUDIO_ROUTE_ERROR_SCHEMA_VERSION))
    code = _studio_route_html_escape(get(payload, "code", "route_error"))
    route = _studio_route_html_escape(get(payload, "route", ""))
    message = _studio_route_html_escape(get(payload, "message", "Unknown Studio route error."))
    suggested = _studio_route_html_escape(get(payload, "suggested_fix", "Inspect the request and retry."))
    detail = _studio_route_html_escape(get(payload, "developer_detail", get(payload, "message", "")))
    exception = _studio_route_html_escape(get(payload, "exception_type", "Exception"))
    severity = _studio_route_html_escape(get(payload, "severity", "error"))
    return """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>OWENS Studio Error</title>
  <style>
    :root {
      color-scheme: light;
      --ink: #17202a;
      --muted: #5c6670;
      --line: #d9dee4;
      --bad: #9b1c1c;
      --panel: #f8fafc;
    }
    body {
      margin: 0;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: #ffffff;
    }
    main {
      max-width: 860px;
      margin: 0 auto;
      padding: 36px 24px;
    }
    h1 {
      margin: 0 0 12px;
      font-size: 28px;
    }
    .status {
      display: inline-block;
      margin-bottom: 18px;
      padding: 4px 8px;
      border: 1px solid #f3b6b6;
      border-radius: 6px;
      color: var(--bad);
      background: #fff1f2;
      font-size: 13px;
      font-weight: 650;
    }
    section {
      border-top: 1px solid var(--line);
      padding: 16px 0;
    }
    h2 {
      margin: 0 0 8px;
      font-size: 16px;
    }
    p, code {
      overflow-wrap: anywhere;
    }
    p {
      margin: 0;
      line-height: 1.45;
    }
    .detail {
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--panel);
      white-space: pre-wrap;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
    }
    .muted {
      color: var(--muted);
    }
  </style>
</head>
<body>
  <main data-error-schema="$schema" data-error-code="$code" data-error-route="$route">
    <h1>OWENS Studio Error</h1>
    <span class="status">$severity</span>
    <section>
      <h2>What happened</h2>
      <p>$message</p>
    </section>
    <section>
      <h2>Suggested fix</h2>
      <p>$suggested</p>
    </section>
    <section>
      <h2>Developer detail</h2>
      <p class="detail">$detail</p>
      <p class="muted">Exception: <code>$exception</code></p>
    </section>
  </main>
</body>
</html>
"""
end

function _studio_route_html_escape(value)
    return replace(
        string(value),
        "&" => "&amp;",
        "<" => "&lt;",
        ">" => "&gt;",
        "\"" => "&quot;",
        "'" => "&#39;",
    )
end

function _studio_route_contract_record()
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_CONTRACT_SCHEMA_VERSION,
        "request_schema_version" => STUDIO_ROUTE_REQUEST_SCHEMA_VERSION,
        "response_schema_version" => STUDIO_ROUTE_RESPONSE_SCHEMA_VERSION,
        "error_schema_version" => STUDIO_ROUTE_ERROR_SCHEMA_VERSION,
        "transport" => "dependency_light_function_dispatch",
        "primary_serialization" => "application/x-yaml; charset=utf-8",
        "json_supported" => false,
        "json_status_reason" =>
            "The current OWENS Studio shell is dependency-light and YAML-first; JSON should be added when the HTTP server boundary is implemented so both formats share one schema validator.",
        "compatibility" => OrderedCollections.OrderedDict{String,Any}(
            "route_names_are_stable_ids" => true,
            "route_order_is_not_contractual" => true,
            "unknown_response_fields_are_allowed" => true,
            "error_payload_schema_version" => STUDIO_ROUTE_ERROR_SCHEMA_VERSION,
        ),
        "security" => OrderedCollections.OrderedDict{String,Any}(
            "workspace_root_param" => "workspace_root",
            "workspace_root_required" => false,
            "workspace_root_enforced_path_params" => ["project_path", "output_dir", "target"],
            "workspace_root_rejects_allow_external" => true,
            "local_only_server_binding_required_when_http_is_added" => true,
            "csrf_or_non_browser_token_required_when_http_is_added" => true,
        ),
    )
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
        "request_schema" =>
            _studio_route_request_schema(string.(required_params), string.(optional_params)),
        "response_schema" => _studio_route_response_schema(name, content_type),
    )
end

function _studio_route_request_schema(
    required_params::AbstractVector{String},
    optional_params::AbstractVector{String},
)
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_REQUEST_SCHEMA_VERSION,
        "params" => OrderedCollections.OrderedDict{String,Any}(
            "required" => [_studio_route_param_schema(name, true) for name in required_params],
            "optional" => [_studio_route_param_schema(name, false) for name in optional_params],
        ),
        "body" => OrderedCollections.OrderedDict{String,Any}(
            "encoding" => "query_or_form_params",
            "required_for_methods" => ["POST"],
        ),
    )
end

function _studio_route_param_schema(name::AbstractString, required::Bool)
    spec = get(STUDIO_ROUTE_PARAM_SPECS, string(name), ("string", "Undocumented route parameter."))
    return OrderedCollections.OrderedDict{String,Any}(
        "name" => string(name),
        "type" => string(spec[1]),
        "required" => required,
        "description" => string(spec[2]),
    )
end

function _studio_route_response_schema(name::AbstractString, content_type::AbstractString)
    route_name = string(name)
    content = string(content_type)
    payload_schema = get(STUDIO_ROUTE_PAYLOAD_SCHEMAS, route_name, nothing)
    body_kind =
        startswith(content, "text/html") ? "html_document" :
        startswith(content, "text/plain") ? "plain_text" : "yaml_mapping"
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_ROUTE_RESPONSE_SCHEMA_VERSION,
        "success_status" => 200,
        "content_type" => content,
        "body_kind" => body_kind,
        "payload_schema_version" => payload_schema,
        "error_schema_version" => STUDIO_ROUTE_ERROR_SCHEMA_VERSION,
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
    elseif name == "studio_doctor"
        return studio_doctor_route(
            output_dir = _studio_route_optional(params, "output_dir", pwd()),
        )
    elseif name == "capability_catalog"
        return studio_gui_capabilities_route()
    elseif name == "quality_gate_catalog"
        return studio_quality_gates_route()
    elseif name == "template_catalog"
        return studio_project_templates_route()
    elseif name == "example_catalog"
        return studio_project_examples_route()
    elseif name == "project_open"
        return studio_project_open_route(
            _studio_route_required(params, "project_path");
            summarize_runs = _studio_route_optional(params, "summarize_runs", true),
        )
    elseif name == "project_inputs"
        return studio_project_inputs_route(
            _studio_route_required(params, "project_path");
            include_text = _studio_route_optional(params, "include_text", false),
            max_text_bytes = _studio_route_optional(params, "max_text_bytes", 200_000),
        )
    elseif name == "project_editor"
        return studio_project_editor_route(_studio_route_required(params, "project_path"))
    elseif name == "project_input_save"
        return studio_project_input_save_route(
            _studio_route_required(params, "project_path"),
            _studio_route_required(params, "role"),
            _studio_route_required(params, "text");
            expected_sha256 = _studio_route_optional(params, "expected_sha256", nothing),
            allow_external = _studio_route_optional(params, "allow_external", false),
            updated_at_utc = _studio_route_optional(params, "updated_at_utc", nothing),
        )
    elseif name == "project_session"
        return studio_project_session_route(
            _studio_route_required(params, "project_path");
            include_text = _studio_route_optional(params, "include_text", false),
            max_text_bytes = _studio_route_optional(params, "max_text_bytes", 200_000),
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

function _studio_enforce_workspace_root!(route_name::AbstractString, params::AbstractDict)
    haskey(params, "workspace_root") || return params
    root = _studio_route_abs_path(params["workspace_root"])
    isdir(root) ||
        throw(ArgumentError("Studio workspace_root does not exist or is not a directory: $root"))

    if get(params, "allow_external", false) === true
        throw(
            StudioRouteWorkspaceBoundaryError(
                "Studio route $(route_name) cannot set allow_external=true when workspace_root is active.",
            ),
        )
    end

    for key in ("project_path", "output_dir", "target")
        haskey(params, key) || continue
        value = params[key]
        isnothing(value) && continue
        path = _studio_route_abs_path(value)
        _studio_route_path_within_root(path, root) ||
            throw(
                StudioRouteWorkspaceBoundaryError(
                    "Studio route $(route_name) path parameter $key is outside workspace root: $path",
                ),
            )
    end

    delete!(params, "workspace_root")
    return params
end

function _studio_route_abs_path(path)
    path_string = string(path)
    return ispath(path_string) ? realpath(path_string) : normpath(abspath(path_string))
end

function _studio_route_path_within_root(path::AbstractString, root::AbstractString)
    normalized_path = normpath(path)
    normalized_root = normpath(root)
    normalized_path == normalized_root && return true
    separator = Sys.iswindows() ? "\\" : "/"
    return startswith(normalized_path * separator, normalized_root * separator)
end

function _studio_route_required(params::AbstractDict, key::AbstractString)
    haskey(params, key) || throw(ArgumentError("Missing Studio route parameter: $key"))
    return params[key]
end

function _studio_route_optional(params::AbstractDict, key::AbstractString, default)
    return haskey(params, key) ? params[key] : default
end
