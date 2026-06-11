export build_studio_project,
    write_studio_project,
    read_studio_project,
    studio_project_issues,
    validate_studio_project,
    studio_project_health,
    studio_project_generated_script_path,
    read_studio_project_generated_script,
    studio_project_input_summary,
    studio_project_session_summary,
    validate_studio_project_input_text,
    refresh_studio_project,
    save_studio_project_input_text,
    render_studio_project_editor_html,
    render_studio_home_html,
    write_studio_home_html,
    write_studio_project_editor_html,
    render_studio_workbench_html,
    write_studio_workbench_html,
    write_studio_workbench_bundle

const STUDIO_PROJECT_SCHEMA_VERSION = "owens-studio-project/v1"
const STUDIO_WORKBENCH_SCHEMA_VERSION = "owens-studio-workbench/v1"
const STUDIO_WORKBENCH_BUNDLE_SCHEMA_VERSION = "owens-studio-bundle/v1"
const STUDIO_INPUT_SUMMARY_SCHEMA_VERSION = "owens-studio-input-summary/v1"
const STUDIO_SESSION_SCHEMA_VERSION = "owens-studio-session/v1"
const STUDIO_INPUT_VALIDATION_SCHEMA_VERSION = "owens-studio-input-validation/v1"
const STUDIO_INPUT_VALIDATION_ISSUE_SCHEMA_VERSION =
    "owens-studio-validation-issue/v1"
const STUDIO_INPUT_VALIDATION_DOC = "docs/src/getting_started.md"
const STUDIO_INPUT_VALIDATION_SEVERITIES = (
    "parse_error",
    "schema_error",
    "physical_error",
    "unsupported_feature",
    "warning",
    "info",
)
const STUDIO_INPUT_BLOCKING_SEVERITIES =
    Set(["parse_error", "schema_error", "physical_error"])
const STUDIO_NAV_ITEMS = OrderedCollections.OrderedDict{String,Any}[
    OrderedCollections.OrderedDict(
        "label" => "Home",
        "route" => "/",
        "capability" => "project_gallery",
        "status" => "available",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Project",
        "route" => "/workbench",
        "capability" => "project_workbench",
        "status" => "available",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Geometry",
        "capability" => "geometry_editor",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Airfoils",
        "capability" => "airfoil_polar_manager",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Structure",
        "capability" => "structural_composites_editor",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Environment",
        "capability" => "environment_editor",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Controls",
        "capability" => "controls_editor",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Simulation",
        "capability" => "run_workflow",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Validation",
        "capability" => "validation_workflow",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Results",
        "capability" => "results_visualization",
        "status" => "planned",
    ),
    OrderedCollections.OrderedDict(
        "label" => "Reports",
        "capability" => "report_builder",
        "status" => "planned",
    ),
]
const STUDIO_PROJECT_REQUIRED_KEYS = [
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
const STUDIO_EDITABLE_INPUT_ROLES = Set(["modeling_options", "windio"])

"""
    build_studio_project(root; kwargs...)

Build a deterministic OWENS Studio project manifest. The manifest stores
project-level file references and run-manifest references with SHA-256
provenance records so the GUI can detect stale or missing artifacts before a
run or validation report is trusted.
"""
function build_studio_project(
    root::AbstractString;
    project_id::AbstractString = "owens-studio-project",
    name::AbstractString = "OWENS Studio Project",
    description::AbstractString = "",
    modeling_options_file = nothing,
    windio_file = nothing,
    run_manifests = String[],
    metadata = OrderedCollections.OrderedDict{String,Any}(),
    created_at_utc = nothing,
    updated_at_utc = nothing,
)
    project_root = abspath(root)
    isdir(project_root) ||
        throw(ArgumentError("Project root does not exist: $project_root"))
    created = isnothing(created_at_utc) ? _utc_timestamp() : string(created_at_utc)
    updated = isnothing(updated_at_utc) ? created : string(updated_at_utc)

    files = OrderedCollections.OrderedDict{String,Any}[]
    _push_studio_file_record!(
        files,
        modeling_options_file,
        project_root,
        "modeling_options",
    )
    _push_studio_file_record!(files, windio_file, project_root, "windio")

    runs = OrderedCollections.OrderedDict{String,Any}[]
    for run_manifest in run_manifests
        _push_studio_file_record!(runs, run_manifest, project_root, "run_manifest")
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_PROJECT_SCHEMA_VERSION,
        "project_id" => string(project_id),
        "name" => string(name),
        "description" => string(description),
        "created_at_utc" => created,
        "updated_at_utc" => updated,
        "root" => project_root,
        "files" => files,
        "runs" => runs,
        "metadata" => _studio_project_value(metadata),
    )
end

"""
    write_studio_project(path, project)
    write_studio_project(path, root; kwargs...)

Write an OWENS Studio project manifest and return the manifest written.
"""
function write_studio_project(path::AbstractString, project::AbstractDict)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    normalized = _studio_project_value(project)
    YAML.write_file(path, normalized)
    return normalized
end

function write_studio_project(path::AbstractString, root::AbstractString; kwargs...)
    project = build_studio_project(root; kwargs...)
    write_studio_project(path, project)
    return project
end

"""
    read_studio_project(path)

Read an OWENS Studio project manifest using string-keyed ordered dictionaries.
"""
function read_studio_project(path::AbstractString)
    isfile(path) || throw(ArgumentError("Cannot read missing Studio project: $path"))
    return YAML.load_file(path; dicttype = OrderedCollections.OrderedDict{String,Any})
end

"""
    studio_project_issues(project_or_path)

Return deterministic schema and type diagnostics for an OWENS Studio project
manifest. This validates project shape only; use `studio_project_health` to
recompute file hashes and run-manifest health.
"""
studio_project_issues(path::AbstractString) =
    studio_project_issues(read_studio_project(path))

function studio_project_issues(project::AbstractDict)
    issues = String[]
    for key in STUDIO_PROJECT_REQUIRED_KEYS
        haskey(project, key) || push!(issues, "missing required key: $key")
    end

    if haskey(project, "schema_version") &&
       project["schema_version"] != STUDIO_PROJECT_SCHEMA_VERSION
        push!(issues, "schema_version must equal $STUDIO_PROJECT_SCHEMA_VERSION")
    end

    for key in
        ("project_id", "name", "description", "created_at_utc", "updated_at_utc", "root")
        _require_studio_string!(issues, project, key)
    end

    for key in ("files", "runs")
        _require_studio_vector!(issues, project, key)
    end
    _require_studio_dict!(issues, project, "metadata")

    for section in ("files", "runs")
        if haskey(project, section) && project[section] isa AbstractVector
            for (index, record) in enumerate(project[section])
                for issue in _studio_file_record_issues(record)
                    push!(issues, "$section[$index].$issue")
                end
            end
        end
    end

    return issues
end

"""
    validate_studio_project(project_or_path)

Return a valid Studio project, or throw `ArgumentError` with all schema
diagnostics.
"""
validate_studio_project(path::AbstractString) =
    validate_studio_project(read_studio_project(path))

function validate_studio_project(project::AbstractDict)
    issues = studio_project_issues(project)
    isempty(issues) ||
        throw(ArgumentError("Invalid OWENS Studio project:\n- " * join(issues, "\n- ")))
    return project
end

"""
    studio_project_health(project_or_path; root=nothing, summarize_runs=true)

Return project-health rows for Studio file references and run manifests. Healthy
run-manifest references can include nested `run_manifest_health` summaries, so
the first GUI health panel can display stale project files, stale runs, and
missing result-channel metadata through one API.
Relative project roots are resolved against the project manifest file when a
path is supplied, which keeps committed GUI fixtures relocatable across clones.
"""
studio_project_health(path::AbstractString; kwargs...) =
    studio_project_health(read_studio_project(path); project_path = path, kwargs...)

function studio_project_health(
    project::AbstractDict;
    project_path = nothing,
    root = nothing,
    summarize_runs::Bool = true,
)
    project_issues = studio_project_issues(project)
    project_root = _studio_project_root(project, project_path, root)
    file_rows = _studio_file_rows(project, "files", project_root)
    run_rows = _studio_file_rows(project, "runs", project_root)
    if summarize_runs
        for row in run_rows
            if row["status"] == "ok"
                row["run_manifest_health"] = run_manifest_health(row["resolved_path"])
            end
        end
    end

    rows = vcat(file_rows, run_rows)
    nested_runs_ok =
        !summarize_runs ||
        all(
            row ->
                !haskey(row, "run_manifest_health") ||
                    row["run_manifest_health"]["status"] == "ok",
            run_rows,
        )
    status =
        isempty(project_issues) && all(row["status"] == "ok" for row in rows) &&
        nested_runs_ok ? "ok" :
        "attention"

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_WORKBENCH_SCHEMA_VERSION,
        "status" => status,
        "project_path" => isnothing(project_path) ? nothing : abspath(string(project_path)),
        "root" => project_root,
        "project_id" => get(project, "project_id", nothing),
        "name" => get(project, "name", nothing),
        "metadata" => _studio_project_value(
            get(project, "metadata", OrderedCollections.OrderedDict{String,Any}()),
        ),
        "project_issues" => project_issues,
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "records" => length(rows),
            "ok" => count(row -> row["status"] == "ok", rows),
            "modified" => count(row -> row["status"] == "modified", rows),
            "missing" => count(row -> row["status"] == "missing", rows),
            "invalid_record" => count(row -> row["status"] == "invalid_record", rows),
        ),
        "files" => file_rows,
        "runs" => run_rows,
    )
end

"""
    studio_project_generated_script_path(project_or_path; root=nothing)

Return the generated Julia script path recorded in a Studio project manifest, or
`nothing` when the project has no generated-script metadata. Relative script
paths are resolved against the project root.
"""
studio_project_generated_script_path(path::AbstractString; kwargs...) =
    studio_project_generated_script_path(
        read_studio_project(path);
        project_path = path,
        kwargs...,
    )

function studio_project_generated_script_path(
    project::AbstractDict;
    project_path = nothing,
    root = nothing,
)
    project_root = _studio_project_root(project, project_path, root)
    metadata = get(project, "metadata", OrderedCollections.OrderedDict{String,Any}())
    metadata isa AbstractDict || return nothing
    script_ref = get(metadata, "generated_script", nothing)
    script_ref isa AbstractString || return nothing
    isempty(script_ref) && return nothing

    return isabspath(script_ref) ? normpath(script_ref) :
           normpath(joinpath(project_root, script_ref))
end

"""
    read_studio_project_generated_script(project_or_path; required=true)

Read the generated Julia script referenced by a Studio project manifest. When
`required=false`, projects without generated-script metadata return `nothing`
instead of throwing.
"""
read_studio_project_generated_script(path::AbstractString; kwargs...) =
    read_studio_project_generated_script(
        read_studio_project(path);
        project_path = path,
        kwargs...,
    )

function read_studio_project_generated_script(
    project::AbstractDict;
    project_path = nothing,
    root = nothing,
    required::Bool = true,
)
    script_path = studio_project_generated_script_path(project; project_path, root)
    if isnothing(script_path)
        required || return nothing
        throw(ArgumentError("Studio project does not define metadata.generated_script"))
    end
    isfile(script_path) ||
        throw(ArgumentError("Generated Studio script does not exist: $script_path"))

    return read(script_path, String)
end

"""
    studio_project_input_summary(project_or_path; include_text=false, max_text_bytes=200_000)

Return a versioned summary of Studio project input files for GUI editor panes.
The summary preserves provenance health, resolves paths, identifies editable
YAML/text inputs, lists YAML top-level keys, and can include bounded file text
when `include_text=true`.
"""
studio_project_input_summary(path::AbstractString; kwargs...) =
    studio_project_input_summary(read_studio_project(path); project_path = path, kwargs...)

function studio_project_input_summary(
    project::AbstractDict;
    project_path = nothing,
    root = nothing,
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    max_text_bytes < 0 && throw(ArgumentError("max_text_bytes must be non-negative"))
    health = studio_project_health(project; project_path, root, summarize_runs = false)
    input_rows = [
        _studio_input_file_summary(
            row;
            include_text,
            max_text_bytes,
            project_root = health["root"],
        ) for row in health["files"]
    ]
    capability_gates = _studio_project_capability_gates(input_rows)

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_INPUT_SUMMARY_SCHEMA_VERSION,
        "project_path" => get(health, "project_path", nothing),
        "root" => health["root"],
        "project_id" => health["project_id"],
        "name" => health["name"],
        "capability_gates" => capability_gates,
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "files" => length(input_rows),
            "editable" => count(row -> row["editable"] === true, input_rows),
            "parse_errors" => count(row -> row["parse_status"] == "error", input_rows),
            "text_included" =>
                count(row -> get(row, "text_included", false) === true, input_rows),
        ),
        "files" => input_rows,
    )
end

"""
    studio_project_session_summary(project_or_path; include_text=false, max_text_bytes=200_000)

Return a versioned, read-only Studio session snapshot for an active project.
The snapshot classifies file provenance, editable input parse state, and
external-change conflicts so the GUI can decide whether it is safe to save,
reload, or launch without relying on hidden frontend state.
"""
studio_project_session_summary(path::AbstractString; kwargs...) =
    studio_project_session_summary(read_studio_project(path); project_path = path, kwargs...)

function studio_project_session_summary(
    project::AbstractDict;
    project_path = nothing,
    root = nothing,
    include_text::Bool = false,
    max_text_bytes::Integer = 200_000,
)
    health = studio_project_health(project; project_path, root, summarize_runs = false)
    inputs = studio_project_input_summary(
        project;
        project_path,
        root = health["root"],
        include_text,
        max_text_bytes,
    )
    file_states = [_studio_session_file_state(row) for row in inputs["files"]]
    external_changes = [row for row in file_states if row["external_change"] === true]
    parse_errors = [row for row in file_states if row["parse_error"] === true]
    unavailable = [row for row in file_states if row["available"] === false]
    save_conflicts = [
        row for row in file_states if
        row["editable"] === true && row["external_change"] === true
    ]
    reload_required =
        !isempty(external_changes) || !isempty(parse_errors) || !isempty(unavailable)
    session_state =
        reload_required ? "needs_reload" :
        health["status"] == "ok" ? "clean" : "attention"

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_SESSION_SCHEMA_VERSION,
        "project_path" => get(health, "project_path", nothing),
        "root" => health["root"],
        "project_id" => health["project_id"],
        "name" => health["name"],
        "session_state" => session_state,
        "dirty" => reload_required,
        "reload_required" => reload_required,
        "active_project" => OrderedCollections.OrderedDict{String,Any}(
            "project_path" => get(health, "project_path", nothing),
            "root" => health["root"],
            "project_id" => health["project_id"],
            "name" => health["name"],
        ),
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "files" => length(file_states),
            "external_changes" => length(external_changes),
            "parse_errors" => length(parse_errors),
            "unavailable" => length(unavailable),
            "save_conflicts" => length(save_conflicts),
        ),
        "file_states" => file_states,
        "save_conflicts" => save_conflicts,
        "inputs" => inputs,
    )
end

"""
    refresh_studio_project(path; updated_at_utc=nothing)

Recompute provenance records for a Studio project manifest after inputs or run
artifacts have changed, write the refreshed manifest, and return it. Missing
tracked files raise `ArgumentError` so the GUI does not silently trust a stale
project.
"""
function refresh_studio_project(path::AbstractString; updated_at_utc = nothing)
    project = validate_studio_project(path)
    project_root = _studio_project_root(project, path, nothing)
    refreshed = deepcopy(project)
    refreshed["updated_at_utc"] =
        isnothing(updated_at_utc) ? _utc_timestamp() : string(updated_at_utc)
    refreshed["files"] = _refresh_studio_records(project, "files", project_root)
    refreshed["runs"] = _refresh_studio_records(project, "runs", project_root)
    write_studio_project(path, refreshed)
    return refreshed
end

"""
    save_studio_project_input_text(project_path, role, text; kwargs...)

Write new text for an editable Studio project input role, refresh the project
manifest provenance, and return the refreshed project. By default this refuses
to write outside the project root, rejects invalid YAML candidates, and can
enforce an `expected_sha256` optimistic-lock check. The returned project
includes a non-persistent `last_input_save` record with atomic-write, backup,
validation, and before/after provenance details.
"""
function save_studio_project_input_text(
    project_path::AbstractString,
    role::AbstractString,
    text::AbstractString;
    expected_sha256 = nothing,
    allow_external::Bool = false,
    updated_at_utc = nothing,
)
    role_key = string(role)
    role_key in STUDIO_EDITABLE_INPUT_ROLES ||
        throw(ArgumentError("Studio input role is not editable: $role"))
    project = validate_studio_project(project_path)
    project_root = _studio_project_root(project, project_path, nothing)
    record = _studio_project_record_for_role(project, "files", role_key)
    isnothing(record) &&
        throw(ArgumentError("Studio project does not define editable input role: $role"))
    input_path = _resolve_manifest_file_path(record["path"], project_root)
    isfile(input_path) ||
        throw(ArgumentError("Cannot edit missing Studio input file: $input_path"))
    if !allow_external && !_studio_path_within_root(input_path, project_root)
        throw(
            ArgumentError(
                "Refusing to edit Studio input outside project root: $input_path",
            ),
        )
    end

    before_record = file_provenance(input_path; root = project_root, role = role_key)
    if !isnothing(expected_sha256) && string(expected_sha256) != before_record["sha256"]
        throw(
            ArgumentError(
                "Refusing to overwrite modified Studio input $role: expected SHA-256 $(expected_sha256), found $(before_record["sha256"])",
            ),
        )
    end
    validation = validate_studio_project_input_text(role_key, input_path, text)
    if validation["blocking"] === true
        throw(
            ArgumentError(
                "Refusing to save invalid Studio input $role: $(_studio_validation_message(validation))",
            ),
        )
    end

    updated = isnothing(updated_at_utc) ? _utc_timestamp() : string(updated_at_utc)
    backup_path = _backup_studio_input_file(input_path, project_root, role_key, updated)
    _write_studio_text_atomically(input_path, text)

    refreshed = refresh_studio_project(project_path; updated_at_utc = updated)
    after_record = file_provenance(input_path; root = project_root, role = role_key)
    backup_record =
        file_provenance(backup_path; root = project_root, role = "$(role_key)_backup")
    refreshed["last_input_save"] = OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => "owens-studio-input-save-provenance/v1",
        "role" => role_key,
        "path" => before_record["path"],
        "atomic_write" => true,
        "validation" => validation,
        "backup" => backup_record,
        "before" => before_record,
        "after" => after_record,
    )
    return refreshed
end

"""
    validate_studio_project_input_text(role, path, text)

Validate candidate text for an editable Studio project input before saving it.
The returned payload uses severity levels `parse_error`, `schema_error`,
`physical_error`, `unsupported_feature`, `warning`, and `info`; parse/schema/
physical issues are blocking. This is intentionally a save gate, not a full
solver validation report.
"""
function validate_studio_project_input_text(
    role::AbstractString,
    path::AbstractString,
    text::AbstractString,
)
    role_key = string(role)
    issues = OrderedCollections.OrderedDict{String,Any}[]
    format = _studio_input_file_format(path)
    top_level_keys = String[]

    if format == "yaml"
        parsed = nothing
        try
            parsed = YAML.load(text; dicttype = OrderedCollections.OrderedDict{String,Any})
        catch err
            _push_studio_validation_issue!(
                issues,
                "parse_error",
                "Candidate text is not valid YAML: $(sprint(showerror, err))",
                role = role_key,
                path = path,
                field = "document",
                yaml_path = "document",
                physical_implication = "OWENS cannot build a trustworthy input model from malformed YAML.",
                suggested_fix = "Fix the YAML syntax near the parser error and retry the save.",
                remediation_action = "fix_yaml_syntax",
            )
        end

        if isempty(issues)
            if parsed isa AbstractDict
                top_level_keys = String.(collect(Base.keys(parsed)))
                _validate_studio_yaml_schema!(issues, role_key, parsed, top_level_keys)
            else
                _push_studio_validation_issue!(
                    issues,
                    "schema_error",
                    "Candidate YAML root must be a mapping for Studio input role $role_key.",
                    role = role_key,
                    path = path,
                    field = "document",
                    yaml_path = "document",
                    physical_implication = "OWENS expects structured key-value input and cannot derive solver settings from a scalar or sequence root.",
                    suggested_fix = "Replace the YAML root with a mapping that contains the required Studio input keys.",
                    remediation_action = "replace_yaml_root_with_mapping",
                )
            end
        end
    else
        _push_studio_validation_issue!(
            issues,
            "info",
            "Input format $format is treated as plain text; YAML schema checks were skipped.",
            role = role_key,
            path = path,
            field = role_key,
            yaml_path = nothing,
            physical_implication = "No YAML-specific solver schema checks were applied to this file.",
            suggested_fix = "No action is required unless this file should be a YAML Studio input.",
            remediation_action = "no_action_plain_text_input",
        )
    end

    blocking = any(issue -> issue["blocking"] === true, issues)
    status = blocking ? "error" : (isempty(issues) ? "ok" : "warning")
    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_INPUT_VALIDATION_SCHEMA_VERSION,
        "role" => role_key,
        "format" => format,
        "status" => status,
        "blocking" => blocking,
        "top_level_keys" => top_level_keys,
        "issues" => issues,
    )
end

"""
    render_studio_project_editor_html(project_or_inputs; kwargs...)

Render a dependency-light OWENS Studio project editor page. The page exposes
one save form per editable input role and includes the SHA-256 optimistic-lock
token that `save_studio_project_input_text` expects.
"""
function render_studio_project_editor_html(
    project_or_inputs;
    save_action::AbstractString = "/api/project/input",
    workbench_href = nothing,
    inputs_href = nothing,
)
    inputs = _studio_input_summary_input(project_or_inputs)
    title_text = string(get(inputs, "name", "OWENS Studio Project"))
    title = _html_escape(title_text)
    project_path = get(inputs, "project_path", nothing)
    page_css = """
    main {
      padding: 22px 28px;
      min-width: 0;
    }
    aside {
      border-left: 1px solid var(--line);
      background: var(--panel);
      padding: 22px 18px;
    }
    .toolbar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 18px;
    }
    .toolbar h2 {
      font-size: 24px;
      margin: 0;
    }
    .editor {
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      margin-bottom: 18px;
      overflow: hidden;
    }
    .editor header {
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 12px;
      padding: 14px 16px;
      border-bottom: 1px solid var(--line);
      background: #f8fafc;
    }
    .editor h3 {
      font-size: 17px;
      margin: 0 0 6px;
    }
    .editor p, aside p {
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
      margin: 0 0 8px;
      overflow-wrap: anywhere;
    }
    .pill {
      display: inline-flex;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 4px 8px;
      font-size: 12px;
      color: var(--muted);
      background: #fff;
      white-space: nowrap;
    }
    .pill.ok { color: var(--ok); background: #eef8f2; }
    .pill.bad { color: var(--bad); background: #fff1f2; }
    .pill.attention { color: var(--attention); background: #fff7ed; }
    .editor form {
      padding: 14px 16px 16px;
    }
    textarea {
      width: 100%;
      min-height: 260px;
      resize: vertical;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 10px 12px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      line-height: 1.45;
      color: var(--ink);
      background: #fff;
    }
    label {
      display: block;
      margin: 0 0 6px;
      font-weight: 650;
      font-size: 13px;
    }
    .field-help {
      margin: 0 0 8px;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
      overflow-wrap: anywhere;
    }
    textarea[readonly] {
      background: #f3f5f7;
      color: var(--muted);
    }
    textarea:focus,
    button:focus-visible,
    nav .nav-item:focus-visible,
    aside a:focus-visible {
      outline: 3px solid #74a8d8;
      outline-offset: 2px;
    }
    button {
      margin-top: 10px;
      border: 1px solid var(--blue);
      border-radius: 6px;
      background: var(--blue);
      color: #fff;
      padding: 8px 12px;
      font-weight: 650;
    }
    button:disabled {
      border-color: var(--line);
      background: #d9dee4;
      color: var(--muted);
    }
    .meta {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 8px;
    }
    aside a {
      color: var(--blue);
    }
    .empty {
      color: var(--muted);
      margin: 0;
    }
    .validation-issues {
      display: grid;
      gap: 8px;
      margin: 10px 0 0;
    }
    .validation-issue {
      border: 1px solid var(--line);
      border-left: 4px solid var(--attention);
      border-radius: 6px;
      padding: 10px 12px;
      background: #fffaf0;
    }
    .validation-issue.blocking {
      border-left-color: var(--bad);
      background: #fff7f7;
    }
    .validation-issue p {
      margin: 4px 0 0;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
      overflow-wrap: anywhere;
    }
    .validation-issue .issue-title {
      margin-top: 0;
      color: var(--ink);
      font-weight: 650;
    }
$(_studio_shell_mobile_css())"""

    return """
<!doctype html>
<html lang="en">
$(_studio_document_head("$title_text - OWENS Studio Editor", page_css))
<body>
  <div class="studio">
    <nav>
      <h1>OWENS Studio</h1>
      $(_studio_nav_html("Project"; project_path))
    </nav>
    <main>
      <div class="toolbar">
        <h2>$title</h2>
      </div>
      $(_studio_input_editor_forms_html(inputs, string(save_action)))
    </main>
    <aside>
      <h3>Project</h3>
      <p>$(_html_escape(string(project_path)))</p>
      <h3>Root</h3>
      <p>$(_html_escape(string(get(inputs, "root", ""))))</p>
      <h3>Input Summary</h3>
      $(_studio_input_summary_cards_html(inputs["summary"]))
      $(_studio_editor_artifact_links_html(workbench_href, inputs_href))
    </aside>
  </div>
</body>
</html>
"""
end

"""
    render_studio_workbench_html(project_or_health; health_href=nothing, script_href=nothing, open_href=nothing, editor_href=nothing)

Render a dependency-light OWENS Studio workbench shell as static HTML. This is
the first GUI slice: a project health view that later Genie routes can serve
without changing the underlying health API.
"""
function render_studio_workbench_html(
    project_or_health;
    health_href = nothing,
    script_href = nothing,
    open_href = nothing,
    editor_href = nothing,
)
    health = _studio_health_input(project_or_health)
    title_text = string(get(health, "name", "OWENS Studio"))
    title = _html_escape(title_text)
    status = _html_escape(string(health["status"]))
    page_css = """
    main {
      padding: 22px 28px;
      min-width: 0;
    }
    aside {
      border-left: 1px solid var(--line);
      background: var(--panel);
      padding: 22px 18px;
    }
    .toolbar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 18px;
    }
    .toolbar h2 {
      font-size: 24px;
      margin: 0;
    }
    .status {
      display: inline-flex;
      align-items: center;
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 5px 10px;
      font-size: 13px;
      font-weight: 650;
      text-transform: uppercase;
      letter-spacing: 0;
      background: #fff;
    }
    .status.ok { color: var(--ok); }
    .status.attention { color: var(--attention); }
    .cards {
      display: grid;
      grid-template-columns: repeat(5, minmax(90px, 1fr));
      gap: 10px;
      margin: 0 0 22px;
    }
    .metric {
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 12px;
      background: #fff;
    }
    .metric strong {
      display: block;
      font-size: 24px;
      line-height: 1.1;
    }
    .metric span {
      color: var(--muted);
      font-size: 12px;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      margin-bottom: 22px;
    }
    th, td {
      text-align: left;
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }
    th {
      background: #f3f5f7;
      color: var(--muted);
      font-weight: 650;
    }
    tr:last-child td { border-bottom: 0; }
    .issue {
      color: var(--bad);
      margin: 0 0 7px;
      overflow-wrap: anywhere;
    }
    .empty {
      color: var(--muted);
      margin: 0;
    }
$(_studio_shell_mobile_css())
    @media (max-width: 980px) {
      .cards { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
    }"""

    return """
<!doctype html>
<html lang="en">
$(_studio_document_head("$title_text - OWENS Studio", page_css))
<body>
  <div class="studio">
    <nav>
      <h1>OWENS Studio</h1>
      $(_studio_nav_html("Project"; project_path = get(health, "project_path", nothing)))
    </nav>
    <main>
      <div class="toolbar">
        <h2>$title</h2>
        <span class="status $status">$status</span>
      </div>
      $(_studio_metric_cards_html(health["summary"]))
      <h3>Project Files</h3>
      $(_studio_records_table_html(health["files"]))
      <h3>Run Manifests</h3>
      $(_studio_records_table_html(health["runs"]))
    </main>
    <aside>
      <h3>Project Health</h3>
      $(_studio_issues_html(health["project_issues"]))
      <h3>Root</h3>
      <p>$(_html_escape(string(health["root"])))</p>
      <h3>Project Manifest</h3>
      <p>$(_html_escape(string(get(health, "project_path", nothing))))</p>
      $(_studio_artifact_links_html(health_href, script_href, open_href, editor_href))
      $(_studio_generated_script_html(health))
    </aside>
  </div>
</body>
</html>
"""
end

"""
    render_studio_home_html(; templates=studio_project_template_catalog(), examples=studio_example_project_catalog())

Render the dependency-light OWENS Studio project chooser. The generated HTML is
static so the same surface can be written to disk, served by the app route
layer, or wrapped later by Genie.
"""
function render_studio_home_html(;
    templates = studio_project_template_catalog(),
    examples = studio_example_project_catalog(),
)
    page_css = """
    main {
      padding: 24px 30px;
      min-width: 0;
    }
    .toolbar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 22px;
    }
    .toolbar h2 {
      font-size: 24px;
      margin: 0;
    }
    .grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 12px;
      margin-bottom: 26px;
    }
    .card {
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      background: #fff;
      min-height: 132px;
    }
    .card h4 {
      font-size: 16px;
      margin: 0 0 8px;
    }
    .card p {
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
      margin: 0 0 10px;
      overflow-wrap: anywhere;
    }
    .meta {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
    }
    .pill {
      border: 1px solid var(--line);
      border-radius: 999px;
      padding: 4px 8px;
      font-size: 12px;
      color: var(--muted);
      background: var(--panel);
    }
    .pill.ok { color: var(--green); background: #eef8f2; }
    .empty {
      color: var(--muted);
      margin: 0;
    }
$(_studio_shell_mobile_css(include_aside = false, max_width = 840))"""

    return """
<!doctype html>
<html lang="en">
$(_studio_document_head("OWENS Studio", page_css; columns = "220px minmax(0, 1fr)"))
<body>
  <div class="studio">
    <nav>
      <h1>OWENS Studio</h1>
      $(_studio_nav_html("Home"))
    </nav>
    <main>
      <div class="toolbar">
        <h2>Project Gallery</h2>
      </div>
      <h3>Example Projects</h3>
      $(_studio_example_cards_html(examples["examples"]))
      <h3>New Project Templates</h3>
      $(_studio_template_cards_html(templates["templates"]))
    </main>
  </div>
</body>
</html>
"""
end

"""
    write_studio_home_html(path; kwargs...)

Write the static OWENS Studio project chooser and return the rendered HTML.
"""
function write_studio_home_html(path::AbstractString; kwargs...)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    html = render_studio_home_html(; kwargs...)
    open(path, "w") do io
        write(io, html)
    end

    return html
end

"""
    write_studio_workbench_html(path, project_or_health)

Write the static OWENS Studio workbench shell and return the rendered HTML.
"""
function write_studio_workbench_html(path::AbstractString, project_or_health)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    html = render_studio_workbench_html(project_or_health)
    open(path, "w") do io
        write(io, html)
    end

    return html
end

"""
    write_studio_project_editor_html(path, project_or_inputs; kwargs...)

Write the static OWENS Studio project editor shell and return the rendered
HTML.
"""
function write_studio_project_editor_html(
    path::AbstractString,
    project_or_inputs;
    kwargs...,
)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    html = render_studio_project_editor_html(project_or_inputs; kwargs...)
    open(path, "w") do io
        write(io, html)
    end

    return html
end

"""
    write_studio_workbench_bundle(output_dir, project_path; include_script=true)

Write a static OWENS Studio workbench bundle containing `index.html`,
`health.yml`, and, when available, `generated_script.jl`. This provides a
server-free GUI artifact for review and gives the future web shell the exact
files it should serve.
"""
function write_studio_workbench_bundle(
    output_dir::AbstractString,
    project_path::AbstractString;
    include_script::Bool = true,
)
    bundle_dir = abspath(output_dir)
    mkpath(bundle_dir)

    health = studio_project_health(project_path)
    inputs = studio_project_input_summary(project_path; include_text = true)
    health_file = joinpath(bundle_dir, "health.yml")
    inputs_file = joinpath(bundle_dir, "inputs.yml")
    YAML.write_file(health_file, health)
    YAML.write_file(inputs_file, inputs)

    script_file = nothing
    if include_script
        script = read_studio_project_generated_script(project_path; required = false)
        if !isnothing(script)
            script_file = joinpath(bundle_dir, "generated_script.jl")
            open(script_file, "w") do io
                write(io, script)
            end
        end
    end

    index_file = joinpath(bundle_dir, "index.html")
    html = render_studio_workbench_html(
        health;
        health_href = basename(health_file),
        script_href = isnothing(script_file) ? nothing : basename(script_file),
        editor_href = "editor.html",
    )
    open(index_file, "w") do io
        write(io, html)
    end

    editor_file = joinpath(bundle_dir, "editor.html")
    editor_html = render_studio_project_editor_html(
        inputs;
        workbench_href = basename(index_file),
        inputs_href = basename(inputs_file),
    )
    open(editor_file, "w") do io
        write(io, editor_html)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_WORKBENCH_BUNDLE_SCHEMA_VERSION,
        "bundle_dir" => bundle_dir,
        "project_file" => abspath(project_path),
        "index_html" => index_file,
        "editor_html" => editor_file,
        "health_file" => health_file,
        "inputs_file" => inputs_file,
        "script_file" => script_file,
        "project_status" => health["status"],
        "bytes" => OrderedCollections.OrderedDict{String,Any}(
            "index_html" => stat(index_file).size,
            "editor_html" => stat(editor_file).size,
            "health_file" => stat(health_file).size,
            "inputs_file" => stat(inputs_file).size,
            "script_file" => isnothing(script_file) ? nothing : stat(script_file).size,
        ),
    )
end

function _push_studio_file_record!(
    records,
    path,
    root::AbstractString,
    role::AbstractString,
)
    if !isnothing(path)
        push!(records, file_provenance(string(path); root, role))
    end
    return records
end

function _studio_project_root(project::AbstractDict, project_path, root)
    if !isnothing(root)
        return _canonical_abs_path(string(root))
    elseif haskey(project, "root") && project["root"] isa AbstractString
        project_root = project["root"]
        if isabspath(project_root) || isnothing(project_path)
            return _canonical_abs_path(project_root)
        end
        return _canonical_abs_path(
            joinpath(dirname(abspath(string(project_path))), project_root),
        )
    elseif !isnothing(project_path)
        return _canonical_abs_path(dirname(abspath(string(project_path))))
    end

    return _canonical_abs_path(pwd())
end

function _backup_studio_input_file(
    input_path::AbstractString,
    project_root::AbstractString,
    role::AbstractString,
    timestamp::AbstractString,
)
    history_dir = joinpath(project_root, ".owens-studio-history")
    mkpath(history_dir)
    safe_timestamp = _studio_safe_filename_label(timestamp)
    safe_role = _studio_safe_filename_label(role)
    backup_base = joinpath(history_dir, "$safe_timestamp-$safe_role-$(basename(input_path))")
    backup_path = _unique_studio_backup_path(backup_base)
    cp(input_path, backup_path; force = false)
    return backup_path
end

function _write_studio_text_atomically(path::AbstractString, text::AbstractString)
    parent = dirname(path)
    mkpath(parent)
    temp_path = tempname(parent; cleanup = false)
    try
        open(temp_path, "w") do io
            write(io, text)
            _flush_studio_file(io)
        end
        mv(temp_path, path; force = true)
    finally
        ispath(temp_path) && rm(temp_path; force = true)
    end
    return path
end

function _flush_studio_file(io)
    flush(io)
    try
        raw_fd = fd(io)
        ccall(:fsync, Cint, (Cint,), reinterpret(Cint, raw_fd))
    catch
        # Some Julia streams or filesystems do not expose fsync cleanly. The
        # flush above is still required; fsync is best-effort for local disks.
    end
    return nothing
end

function _studio_safe_filename_label(value)
    label = replace(string(value), r"[^A-Za-z0-9_.-]+" => "-")
    label = strip(label, ['-', '.', '_'])
    return isempty(label) ? "value" : label
end

function _unique_studio_backup_path(path::AbstractString)
    candidate = path
    counter = 1
    while ispath(candidate)
        candidate = "$path.$counter"
        counter += 1
    end
    return candidate
end

function _studio_file_rows(
    project::AbstractDict,
    section::AbstractString,
    root::AbstractString,
)
    records = get(project, section, OrderedCollections.OrderedDict{String,Any}[])
    records isa AbstractVector || return OrderedCollections.OrderedDict{String,Any}[]
    return [verify_file_provenance(record; root) for record in records]
end

function _studio_health_input(project_or_health)
    if project_or_health isa AbstractString
        return studio_project_health(project_or_health)
    elseif project_or_health isa AbstractDict &&
           get(project_or_health, "schema_version", nothing) ==
           STUDIO_WORKBENCH_SCHEMA_VERSION
        return project_or_health
    elseif project_or_health isa AbstractDict
        return studio_project_health(project_or_health)
    end

    throw(
        ArgumentError("Expected a Studio project path, project manifest, or health report"),
    )
end

function _studio_file_record_issues(record)
    if !(record isa AbstractDict)
        return String["record must be a dictionary"]
    end

    health = verify_file_provenance(record)
    return health["status"] == "invalid_record" ? String.(health["issues"]) : String[]
end

function _require_studio_string!(
    issues::Vector{String},
    project::AbstractDict,
    key::AbstractString,
)
    haskey(project, key) || return issues
    project[key] isa AbstractString || push!(issues, "$key must be a string")
    return issues
end

function _require_studio_vector!(
    issues::Vector{String},
    project::AbstractDict,
    key::AbstractString,
)
    haskey(project, key) || return issues
    project[key] isa AbstractVector || push!(issues, "$key must be a vector")
    return issues
end

function _require_studio_dict!(
    issues::Vector{String},
    project::AbstractDict,
    key::AbstractString,
)
    haskey(project, key) || return issues
    project[key] isa AbstractDict || push!(issues, "$key must be a dictionary")
    return issues
end

const STUDIO_SHARED_STYLE_VERSION = "owens-studio-shared-style/v1"

function _studio_document_head(
    title::AbstractString,
    page_css::AbstractString;
    columns::AbstractString = "220px minmax(0, 1fr) 320px",
)
    return """
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="owens-studio-style" content="$STUDIO_SHARED_STYLE_VERSION">
  <title>$(_html_escape(title))</title>
  <style>
$(_studio_shared_shell_css(; columns))
$page_css
  </style>
</head>"""
end

function _studio_shared_shell_css(; columns::AbstractString = "220px minmax(0, 1fr) 320px")
    return """
    :root {
      color-scheme: light;
      --ink: #17202a;
      --muted: #5c6670;
      --line: #d9dee4;
      --panel: #f6f8fa;
      --ok: #116b3a;
      --attention: #9a4f00;
      --bad: #9b1c1c;
      --blue: #1f5e99;
      --green: #116b3a;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: #ffffff;
    }
    .studio {
      display: grid;
      grid-template-columns: $columns;
      min-height: 100vh;
    }
    nav {
      border-right: 1px solid var(--line);
      padding: 18px 14px;
      background: #f9fafb;
    }
    nav h1 {
      font-size: 18px;
      margin: 0 0 18px;
    }
    nav .nav-item {
      display: block;
      padding: 9px 10px;
      color: var(--ink);
      text-decoration: none;
      border-radius: 6px;
      font-size: 14px;
    }
    nav .nav-item.active {
      background: #e7f0f8;
      color: var(--blue);
      font-weight: 650;
    }
    nav .nav-item.disabled {
      color: var(--muted);
      cursor: not-allowed;
    }
    nav .nav-status {
      display: inline-block;
      margin-left: 6px;
      font-size: 11px;
      color: var(--muted);
    }"""
end

function _studio_shell_mobile_css(; include_aside::Bool = true, max_width::Integer = 980)
    border_selector = include_aside ? "nav, aside" : "nav"
    border_rule =
        include_aside ? "border: 0; border-bottom: 1px solid var(--line);" :
        "border-right: 0; border-bottom: 1px solid var(--line);"
    return """
    @media (max-width: $(max_width)px) {
      .studio { grid-template-columns: 1fr; }
      $border_selector { $border_rule }
    }"""
end

function _studio_nav_html(
    active_item::AbstractString = "Project";
    project_path = nothing,
)
    return join(
        [
            _studio_nav_item_html(row, active_item; project_path) for
            row in STUDIO_NAV_ITEMS
        ],
        "\n      ",
    )
end

function _studio_nav_item_html(
    row::AbstractDict,
    active_item::AbstractString;
    project_path = nothing,
)
    label = string(row["label"])
    status = string(get(row, "status", "planned"))
    capability = string(get(row, "capability", lowercase(label)))
    active = label == active_item
    href = _studio_nav_href(row; project_path)
    class = "nav-item$(active ? " active" : "")"
    label_html = _html_escape(label)
    capability_html = _html_escape(capability)

    if !isnothing(href)
        aria_current = active ? " aria-current=\"page\"" : ""
        return "<a class=\"$class\" href=\"$(_html_escape(href))\" data-capability=\"$capability_html\" data-status=\"$(_html_escape(status))\"$aria_current>$label_html</a>"
    end

    title = status == "planned" ?
            "$label workflow is planned and not available in this Studio shell yet." :
            "$label workflow is unavailable in this context."
    return "<span class=\"$class disabled\" aria-disabled=\"true\" title=\"$(_html_escape(title))\" data-capability=\"$capability_html\" data-status=\"$(_html_escape(status))\">$label_html<span class=\"nav-status\">$(_html_escape(status))</span></span>"
end

function _studio_nav_href(row::AbstractDict; project_path = nothing)
    status = string(get(row, "status", "planned"))
    status == "available" || return nothing
    route = get(row, "route", nothing)
    route isa AbstractString || return nothing
    if route == "/workbench"
        isnothing(project_path) && return nothing
        return "$route?project_path=$(_studio_query_escape(project_path))"
    end
    return string(route)
end

function _studio_query_escape(value)
    io = IOBuffer()
    for byte in codeunits(string(value))
        if (UInt8('A') <= byte <= UInt8('Z')) ||
           (UInt8('a') <= byte <= UInt8('z')) ||
           (UInt8('0') <= byte <= UInt8('9')) ||
           byte in (UInt8('-'), UInt8('_'), UInt8('.'), UInt8('~'), UInt8('/'))
            write(io, byte)
        else
            write(io, '%')
            write(io, uppercase(string(byte; base = 16, pad = 2)))
        end
    end
    return String(take!(io))
end

function _studio_metric_cards_html(summary::AbstractDict)
    keys = ["records", "ok", "modified", "missing", "invalid_record"]
    labels = ["Records", "OK", "Modified", "Missing", "Invalid"]
    cards = [
        "<div class=\"metric\"><strong>$(_html_escape(string(get(summary, key, 0))))</strong><span>$(label)</span></div>"
        for (key, label) in zip(keys, labels)
    ]
    return "<section class=\"cards\">\n        " *
           join(cards, "\n        ") *
           "\n      </section>"
end

function _studio_example_cards_html(rows::AbstractVector)
    if isempty(rows)
        return "<p class=\"empty\">No example projects.</p>"
    end

    cards = [_studio_example_card_html(row) for row in rows]
    return "<section class=\"grid\">\n        " *
           join(cards, "\n        ") *
           "\n      </section>"
end

function _studio_template_cards_html(rows::AbstractVector)
    if isempty(rows)
        return "<p class=\"empty\">No project templates.</p>"
    end

    cards = [_studio_template_card_html(row) for row in rows]
    return "<section class=\"grid\">\n        " *
           join(cards, "\n        ") *
           "\n      </section>"
end

function _studio_example_card_html(row::AbstractDict)
    available = get(row, "available", false) === true
    project_path = string(get(row, "project_relative_path", get(row, "project_file", "")))
    status_class = available ? " ok" : ""
    status_text = available ? "available" : "missing"
    return """
<article class="card">
  <h4>$(_html_escape(string(get(row, "title", get(row, "example", "Example")))))</h4>
  <p>$(_html_escape(string(get(row, "description", ""))))</p>
  <p>$(_html_escape(project_path))</p>
  <div class="meta">
    <span class="pill$status_class">$(_html_escape(status_text))</span>
    <span class="pill">$(_html_escape(string(get(row, "turbine_type", ""))))</span>
    <span class="pill">$(_html_escape(string(get(row, "template", ""))))</span>
  </div>
</article>"""
end

function _studio_template_card_html(row::AbstractDict)
    generated = get(row, "creates_generated_script", false) === true
    manifest = get(row, "creates_run_manifest", false) === true
    return """
<article class="card">
  <h4>$(_html_escape(string(get(row, "title", get(row, "template", "Template")))))</h4>
  <p>$(_html_escape(string(get(row, "description", ""))))</p>
  <div class="meta">
    <span class="pill">$(_html_escape(string(get(row, "template", ""))))</span>
    <span class="pill">$(_html_escape(string(get(row, "turbine_type", ""))))</span>
    <span class="pill$(generated ? " ok" : "")">$(_html_escape(generated ? "script" : "no script"))</span>
    <span class="pill$(manifest ? " ok" : "")">$(_html_escape(manifest ? "manifest" : "no manifest"))</span>
  </div>
</article>"""
end

function _studio_records_table_html(rows::AbstractVector)
    if isempty(rows)
        return "<p class=\"empty\">No records.</p>"
    end

    body = join([_studio_record_row_html(row) for row in rows], "\n        ")
    return """
<table>
  <thead>
    <tr><th>Status</th><th>Role</th><th>Path</th><th>Issues</th></tr>
  </thead>
  <tbody>
        $body
  </tbody>
</table>"""
end

function _studio_record_row_html(row::AbstractDict)
    issue_parts = String.(get(row, "issues", String[]))
    if haskey(row, "remediation")
        push!(issue_parts, _studio_remediation_summary(row["remediation"]))
    end
    append!(issue_parts, _studio_nested_run_remediation_summaries(row))
    issue_text = isempty(issue_parts) ? "" : join(issue_parts, "; ")
    return "<tr><td>$(_html_escape(string(row["status"])))</td><td>$(_html_escape(string(get(row, "role", ""))))</td><td>$(_html_escape(string(get(row, "path", ""))))</td><td>$(_html_escape(issue_text))</td></tr>"
end

function _studio_nested_run_remediation_summaries(row::AbstractDict)
    health = get(row, "run_manifest_health", nothing)
    health isa AbstractDict || return String[]
    summaries = String[]
    for section in ("inputs", "outputs", "generated")
        records = get(health, section, OrderedCollections.OrderedDict{String,Any}[])
        records isa AbstractVector || continue
        for record in records
            record isa AbstractDict || continue
            haskey(record, "remediation") || continue
            push!(summaries, _studio_remediation_summary(record["remediation"]))
        end
    end
    return summaries
end

function _studio_remediation_summary(remediation)
    remediation isa AbstractDict || return string(remediation)
    parts = String[]
    for key in ("code", "field", "suggested_fix")
        value = get(remediation, key, nothing)
        value isa AbstractString && !isempty(value) || continue
        push!(parts, value)
    end
    return isempty(parts) ? string(remediation) : join(parts, " | ")
end

function _studio_issues_html(issues::AbstractVector)
    if isempty(issues)
        return "<p class=\"empty\">No schema issues.</p>"
    end

    return join(
        ["<p class=\"issue\">$(_html_escape(string(issue)))</p>" for issue in issues],
        "\n",
    )
end

function _studio_artifact_links_html(health_href, script_href, open_href, editor_href)
    links = String[]
    if open_href isa AbstractString && !isempty(open_href)
        push!(links, "<a href=\"$(_html_escape(open_href))\">Open Payload</a>")
    end
    if editor_href isa AbstractString && !isempty(editor_href)
        push!(links, "<a href=\"$(_html_escape(editor_href))\">Editor</a>")
    end
    if health_href isa AbstractString && !isempty(health_href)
        push!(links, "<a href=\"$(_html_escape(health_href))\">Health YAML</a>")
    end
    if script_href isa AbstractString && !isempty(script_href)
        push!(links, "<a href=\"$(_html_escape(script_href))\">Generated Julia</a>")
    end
    isempty(links) && return ""

    return "<h3>Artifacts</h3>\n      <p>" * join(links, " &middot; ") * "</p>"
end

function _studio_generated_script_html(health::AbstractDict)
    metadata = get(health, "metadata", OrderedCollections.OrderedDict{String,Any}())
    metadata isa AbstractDict || return ""
    script = get(metadata, "generated_script", nothing)
    script isa AbstractString || return ""
    isempty(script) && return ""

    return "<h3>Generated Script</h3>\n      <p>$(_html_escape(script))</p>"
end

function _studio_input_file_summary(
    row::AbstractDict;
    include_text::Bool,
    max_text_bytes::Integer,
    project_root::AbstractString,
)
    resolved_path = get(row, "resolved_path", nothing)
    role = get(row, "role", nothing)
    format =
        resolved_path isa AbstractString ? _studio_input_file_format(resolved_path) :
        "unknown"
    available = resolved_path isa AbstractString && isfile(resolved_path)
    editable =
        available &&
        role isa AbstractString &&
        role in STUDIO_EDITABLE_INPUT_ROLES &&
        format in ("yaml", "text") &&
        _studio_path_within_root(resolved_path, project_root)
    parse = _studio_input_parse_summary(resolved_path, format, available, max_text_bytes)
    validation =
        _studio_input_validation_summary(role, resolved_path, format, available, max_text_bytes)

    summary = OrderedCollections.OrderedDict{String,Any}(
        "role" => role,
        "path" => get(row, "path", nothing),
        "resolved_path" => resolved_path,
        "status" => get(row, "status", nothing),
        "issues" => get(row, "issues", String[]),
        "expected_bytes" => get(row, "expected_bytes", nothing),
        "actual_bytes" => get(row, "actual_bytes", nothing),
        "expected_sha256" => get(row, "expected_sha256", nothing),
        "actual_sha256" => get(row, "actual_sha256", nothing),
        "format" => format,
        "editable" => editable,
        "parse_status" => parse["status"],
        "parse_message" => parse["message"],
        "top_level_keys" => parse["top_level_keys"],
        "validation_status" => validation["status"],
        "validation_blocking" => validation["blocking"],
        "validation_issues" => validation["issues"],
        "line_count" => available ? _studio_input_line_count(resolved_path) : nothing,
        "text_included" => false,
        "text_truncated" => false,
    )

    if include_text && available && format in ("yaml", "text")
        actual_bytes = stat(resolved_path).size
        if actual_bytes <= max_text_bytes
            summary["text"] = read(resolved_path, String)
            summary["text_included"] = true
        else
            summary["text_truncated"] = true
            summary["text_limit_bytes"] = Int(max_text_bytes)
        end
    end

    return summary
end

function _studio_session_file_state(row::AbstractDict)
    expected_sha = get(row, "expected_sha256", nothing)
    actual_sha = get(row, "actual_sha256", nothing)
    status = string(get(row, "status", "unknown"))
    parse_status = string(get(row, "parse_status", "unknown"))
    available =
        get(row, "resolved_path", nothing) isa AbstractString &&
        get(row, "actual_bytes", nothing) !== nothing
    external_change =
        status != "ok" ||
        (!isnothing(expected_sha) && !isnothing(actual_sha) && expected_sha != actual_sha)
    parse_error = parse_status == "error"
    needs_reload = external_change || parse_error || !available
    return OrderedCollections.OrderedDict{String,Any}(
        "role" => get(row, "role", nothing),
        "path" => get(row, "path", nothing),
        "resolved_path" => get(row, "resolved_path", nothing),
        "editable" => get(row, "editable", false) === true,
        "available" => available,
        "status" => status,
        "parse_status" => parse_status,
        "external_change" => external_change,
        "parse_error" => parse_error,
        "needs_reload" => needs_reload,
        "expected_sha256" => expected_sha,
        "actual_sha256" => actual_sha,
        "issues" => get(row, "issues", String[]),
    )
end

function _studio_project_capability_gates(input_rows::AbstractVector)
    gates = OrderedCollections.OrderedDict{String,Any}[]
    for row in input_rows
        get(row, "role", nothing) == "windio" || continue
        get(row, "parse_status", nothing) == "ok" || continue
        path = get(row, "resolved_path", nothing)
        path isa AbstractString && isfile(path) || continue
        detected = _studio_windio_turbine_kind(path)
        isnothing(detected) && continue
        _studio_hawt_kind(detected) || continue
        push!(
            gates,
            OrderedCollections.OrderedDict{String,Any}(
                "capability" => "hawt_aeroelastic_workflow",
                "status" => "experimental",
                "severity" => "warning",
                "source_role" => "windio",
                "source_path" => get(row, "path", nothing),
                "detected_value" => detected,
                "message" =>
                    "WindIO declares a HAWT/axial-flow turbine. OWENS HAWT setup is experimental and validation-gated; unsupported production workflows should not run without explicit checks.",
            ),
        )
    end
    return gates
end

function _studio_windio_turbine_kind(path::AbstractString)
    parsed = YAML.load_file(path; dicttype = OrderedCollections.OrderedDict{String,Any})
    parsed isa AbstractDict || return nothing
    return _studio_windio_turbine_kind(parsed)
end

function _studio_windio_turbine_kind(parsed::AbstractDict)
    record = _studio_windio_turbine_kind_record(parsed)
    return isnothing(record) ? nothing : record.value
end

function _studio_windio_turbine_kind_record(parsed::AbstractDict)
    for key_path in (
        ("assembly", "turbine_class"),
        ("assembly", "turbine_type"),
        ("turbine", "turbine_class"),
        ("turbine", "turbine_type"),
        ("turbine_class",),
        ("turbine_type",),
    )
        value = _studio_nested_value(parsed, key_path)
        value isa AbstractString || continue
        isempty(strip(value)) || return (
            value = String(value),
            yaml_path = join(string.(key_path), "."),
        )
    end
    return nothing
end

function _validate_studio_yaml_schema!(
    issues::Vector{OrderedCollections.OrderedDict{String,Any}},
    role::AbstractString,
    parsed::AbstractDict,
    top_level_keys::AbstractVector{String},
)
    if role == "modeling_options"
        if !_studio_has_key(parsed, "OWENS_Options")
            _push_studio_validation_issue!(
                issues,
                "schema_error",
                "Modeling-options YAML must contain a top-level OWENS_Options mapping.",
                role = role,
                field = "OWENS_Options",
                yaml_path = "OWENS_Options",
                physical_implication = "Without OWENS_Options, Studio cannot identify solver time stepping, coupling, or runtime options.",
                suggested_fix = "Add a top-level OWENS_Options mapping before saving.",
                remediation_action = "add_owens_options_mapping",
            )
        elseif !(_studio_dict_value(parsed, "OWENS_Options") isa AbstractDict)
            _push_studio_validation_issue!(
                issues,
                "schema_error",
                "Modeling-options top-level OWENS_Options value must be a mapping.",
                role = role,
                field = "OWENS_Options",
                yaml_path = "OWENS_Options",
                physical_implication = "Studio cannot read solver options from a scalar or sequence OWENS_Options value.",
                suggested_fix = "Change OWENS_Options to a YAML mapping of option names to values.",
                remediation_action = "convert_owens_options_to_mapping",
            )
        end
    elseif role == "windio"
        if isempty(top_level_keys)
            _push_studio_validation_issue!(
                issues,
                "schema_error",
                "WindIO YAML must contain at least one top-level mapping key.",
                role = role,
                field = "document",
                yaml_path = "document",
                physical_implication = "Studio cannot infer turbine geometry, assembly metadata, or backend capability from an empty WindIO document.",
                suggested_fix = "Add the required WindIO project and turbine mappings before saving.",
                remediation_action = "add_windio_top_level_mappings",
            )
        end
        detected_record = _studio_windio_turbine_kind_record(parsed)
        detected = isnothing(detected_record) ? nothing : detected_record.value
        if detected isa AbstractString && _studio_hawt_kind(detected)
            _push_studio_validation_issue!(
                issues,
                "unsupported_feature",
                "WindIO declares a HAWT/axial-flow turbine. OWENS HAWT setup is experimental and validation-gated.",
                role = role,
                field = detected_record.yaml_path,
                yaml_path = detected_record.yaml_path,
                physical_implication = "HAWT setup paths are not yet acceptance-validated across frames, controls, and aerodynamic backends.",
                suggested_fix = "Use a validated VAWT project for production work, or run this HAWT case only through an explicit experimental/validation workflow.",
                remediation_action = "acknowledge_experimental_hawt",
            )
        end
    else
        _push_studio_validation_issue!(
            issues,
            "warning",
            "No role-specific Studio schema checks are available for input role $role.",
            role = role,
            field = role,
            yaml_path = nothing,
            physical_implication = "Studio cannot provide field-level solver checks for this input role yet.",
            suggested_fix = "Review this file manually and add a role-specific validator before relying on GUI validation.",
            remediation_action = "add_role_specific_validator",
        )
    end

    return issues
end

function _push_studio_validation_issue!(
    issues::Vector{OrderedCollections.OrderedDict{String,Any}},
    severity::AbstractString,
    message::AbstractString,
    ;
    role = nothing,
    path = nothing,
    field = nothing,
    yaml_path = nothing,
    physical_implication = nothing,
    suggested_fix = nothing,
    documentation = STUDIO_INPUT_VALIDATION_DOC,
    remediation_action = nothing,
)
    severity in STUDIO_INPUT_VALIDATION_SEVERITIES ||
        throw(ArgumentError("Unsupported Studio validation severity: $severity"))
    push!(
        issues,
        OrderedCollections.OrderedDict{String,Any}(
            "schema_version" => STUDIO_INPUT_VALIDATION_ISSUE_SCHEMA_VERSION,
            "severity" => string(severity),
            "blocking" => severity in STUDIO_INPUT_BLOCKING_SEVERITIES,
            "message" => string(message),
            "role" => _studio_optional_string(role),
            "path" => _studio_optional_string(path),
            "field" => _studio_optional_string(field),
            "yaml_path" => _studio_optional_string(yaml_path),
            "physical_implication" => _studio_optional_string(physical_implication),
            "suggested_fix" => _studio_optional_string(suggested_fix),
            "documentation" => _studio_optional_string(documentation),
            "remediation_action" => _studio_optional_string(remediation_action),
        ),
    )
    return issues
end

_studio_optional_string(value) = isnothing(value) ? nothing : string(value)

function _studio_validation_message(validation::AbstractDict)
    issues = get(validation, "issues", OrderedCollections.OrderedDict{String,Any}[])
    blocking = [
        string(get(issue, "message", "invalid input"))
        for issue in issues
        if get(issue, "blocking", false) === true
    ]
    return isempty(blocking) ? "validation failed" : join(blocking, "; ")
end

function _studio_has_key(data::AbstractDict, key::AbstractString)
    return haskey(data, key) || haskey(data, Symbol(key))
end

function _studio_dict_value(data::AbstractDict, key::AbstractString)
    if haskey(data, key)
        return data[key]
    elseif haskey(data, Symbol(key))
        return data[Symbol(key)]
    end
    return nothing
end

function _studio_nested_value(data, key_path)
    value = data
    for key in key_path
        value isa AbstractDict || return nothing
        if haskey(value, key)
            value = value[key]
        elseif haskey(value, Symbol(key))
            value = value[Symbol(key)]
        else
            return nothing
        end
    end
    return value
end

function _studio_hawt_kind(value::AbstractString)
    normalized = uppercase(replace(value, r"[^A-Za-z0-9]+" => ""))
    return occursin("HAWT", normalized) ||
           occursin("AXIALFLOW", normalized) ||
           normalized == "AXIAL"
end

function _studio_input_file_format(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    if extension in (".yml", ".yaml")
        return "yaml"
    elseif extension in (".jl", ".txt", ".csv", ".dat", ".inp", ".md")
        return "text"
    elseif extension == ".h5"
        return "hdf5"
    end

    return isempty(extension) ? "unknown" : extension[2:end]
end

function _studio_input_parse_summary(
    path,
    format::AbstractString,
    available::Bool,
    max_text_bytes::Integer,
)
    if !available
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "missing",
            "message" => "file is not available",
            "top_level_keys" => String[],
        )
    elseif format != "yaml"
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "not_parsed",
            "message" => "not a YAML file",
            "top_level_keys" => String[],
        )
    elseif path isa AbstractString && stat(path).size > max_text_bytes
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "skipped_size_limit",
            "message" =>
                "YAML parse skipped because file exceeds max_text_bytes=$(max_text_bytes).",
            "top_level_keys" => String[],
        )
    end

    try
        parsed = YAML.load_file(path; dicttype = OrderedCollections.OrderedDict{String,Any})
        top_keys = parsed isa AbstractDict ? String.(collect(Base.keys(parsed))) : String[]
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "ok",
            "message" => "",
            "top_level_keys" => top_keys,
        )
    catch err
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "error",
            "message" => sprint(showerror, err),
            "top_level_keys" => String[],
        )
    end
end

function _studio_input_validation_summary(
    role,
    path,
    format::AbstractString,
    available::Bool,
    max_text_bytes::Integer,
)
    empty_summary = OrderedCollections.OrderedDict{String,Any}(
        "status" => "not_validated",
        "blocking" => false,
        "issues" => OrderedCollections.OrderedDict{String,Any}[],
    )
    available || return empty_summary
    role isa AbstractString || return empty_summary
    role in STUDIO_EDITABLE_INPUT_ROLES || return empty_summary
    format == "yaml" || return empty_summary
    path isa AbstractString || return empty_summary
    if stat(path).size > max_text_bytes
        issues = OrderedCollections.OrderedDict{String,Any}[]
        _push_studio_validation_issue!(
            issues,
            "info",
            "Studio input validation skipped because file exceeds max_text_bytes=$(max_text_bytes).",
            role = role,
            path = path,
            field = "document",
            yaml_path = "document",
            physical_implication = "The GUI did not load or validate this large input file, so field-level solver readiness is unknown.",
            suggested_fix = "Increase max_text_bytes for trusted large inputs, or inspect and validate this file outside the lightweight Studio summary.",
            remediation_action = "increase_max_text_bytes_or_validate_externally",
        )
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "not_validated",
            "blocking" => false,
            "issues" => issues,
        )
    end

    text = read(path, String)
    return validate_studio_project_input_text(role, path, text)
end

function _studio_input_line_count(path::AbstractString)
    count = 0
    open(path, "r") do io
        for _ in eachline(io)
            count += 1
        end
    end
    return count
end

function _refresh_studio_records(
    project::AbstractDict,
    section::AbstractString,
    root::AbstractString,
)
    records = get(project, section, OrderedCollections.OrderedDict{String,Any}[])
    records isa AbstractVector ||
        throw(ArgumentError("Studio project section must be a vector: $section"))
    refreshed = OrderedCollections.OrderedDict{String,Any}[]
    for record in records
        issues = _studio_file_record_issues(record)
        isempty(issues) || throw(
            ArgumentError(
                "Cannot refresh invalid Studio project record: " * join(issues, "; "),
            ),
        )
        resolved_path = _resolve_manifest_file_path(record["path"], root)
        push!(
            refreshed,
            file_provenance(resolved_path; root, role = get(record, "role", nothing)),
        )
    end
    return refreshed
end

function _studio_project_record_for_role(
    project::AbstractDict,
    section::AbstractString,
    role::AbstractString,
)
    records = get(project, section, OrderedCollections.OrderedDict{String,Any}[])
    records isa AbstractVector || return nothing
    for record in records
        if record isa AbstractDict && get(record, "role", nothing) == role
            return record
        end
    end
    return nothing
end

function _studio_path_within_root(path::AbstractString, root::AbstractString)
    path_parts = splitpath(_canonical_abs_path(path))
    root_parts = splitpath(_canonical_abs_path(root))
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _studio_input_summary_input(project_or_inputs)
    if project_or_inputs isa AbstractString
        return studio_project_input_summary(project_or_inputs; include_text = true)
    elseif project_or_inputs isa AbstractDict &&
           get(project_or_inputs, "schema_version", nothing) ==
           STUDIO_INPUT_SUMMARY_SCHEMA_VERSION
        return project_or_inputs
    elseif project_or_inputs isa AbstractDict
        return studio_project_input_summary(project_or_inputs; include_text = true)
    end

    throw(
        ArgumentError("Expected a Studio project path, project manifest, or input summary"),
    )
end

function _studio_input_editor_forms_html(inputs::AbstractDict, save_action::AbstractString)
    rows = get(inputs, "files", OrderedCollections.OrderedDict{String,Any}[])
    rows isa AbstractVector || return "<p class=\"empty\">No editable input files.</p>"
    isempty(rows) && return "<p class=\"empty\">No editable input files.</p>"
    project_path = string(get(inputs, "project_path", ""))
    return join(
        [_studio_input_editor_form_html(row, project_path, save_action) for row in rows],
        "\n      ",
    )
end

function _studio_input_editor_form_html(
    row::AbstractDict,
    project_path::AbstractString,
    save_action::AbstractString,
)
    role = string(get(row, "role", "input"))
    editable = get(row, "editable", false) === true
    status = string(get(row, "status", "unknown"))
    parse_status = string(get(row, "parse_status", "unknown"))
    validation_status = string(get(row, "validation_status", "not_validated"))
    validation_blocking = get(row, "validation_blocking", false) === true
    text = string(get(row, "text", ""))
    sha = get(row, "actual_sha256", nothing)
    key_text = join(string.(get(row, "top_level_keys", String[])), ", ")
    field_id = "studio-input-$(_studio_safe_filename_label(role))"
    help_id = "$field_id-help"
    issues_id = "$field_id-issues"
    validation_issues =
        get(row, "validation_issues", OrderedCollections.OrderedDict{String,Any}[])
    describedby = isempty(validation_issues) ? help_id : "$help_id $issues_id"
    readonly_attr = editable ? "" : " readonly"
    disabled_attr = editable ? "" : " disabled"
    status_class = status == "ok" ? " ok" : " attention"
    parse_class =
        parse_status == "ok" || parse_status == "not_parsed" ? " ok" :
        parse_status == "skipped_size_limit" ? " attention" : " bad"
    validation_class =
        validation_status == "ok" || validation_status == "not_validated" ? " ok" :
        validation_blocking ? " bad" : " attention"
    button_label = editable ? "Save" : "Unavailable"

    return """
<section class="editor">
  <header>
    <div>
      <h3>$(_html_escape(role))</h3>
      <p>$(_html_escape(string(get(row, "path", ""))))</p>
      <div class="meta">
        <span class="pill$status_class">$(_html_escape(status))</span>
        <span class="pill$parse_class">$(_html_escape(parse_status))</span>
        <span class="pill$validation_class">$(_html_escape(validation_status)) validation</span>
        <span class="pill">$(_html_escape(string(get(row, "format", ""))))</span>
        <span class="pill">$(_html_escape(string(get(row, "line_count", "")))) lines</span>
      </div>
    </div>
    <div class="meta">
      <span class="pill">$(_html_escape(string(get(row, "actual_bytes", "")))) bytes</span>
    </div>
  </header>
  <form method="post" action="$(_html_escape(save_action))">
    <input type="hidden" name="project_path" value="$(_html_escape(project_path))">
    <input type="hidden" name="role" value="$(_html_escape(role))">
    <input type="hidden" name="expected_sha256" value="$(_html_escape(string(sha)))">
    <label for="$(_html_escape(field_id))">$(_html_escape(role)) input text</label>
    <p class="field-help" id="$(_html_escape(help_id))">Edit $(_html_escape(string(get(row, "path", "")))) with optimistic-lock SHA-256 save protection.</p>
    <textarea id="$(_html_escape(field_id))" name="text" aria-describedby="$(_html_escape(describedby))" aria-invalid="$(_html_escape(string(validation_blocking)))"$readonly_attr>$(_html_escape(text))</textarea>
    $(_studio_validation_issues_html(validation_issues; id = issues_id))
    <button type="submit"$disabled_attr>$(_html_escape(button_label))</button>
    <div class="meta">
      <span class="pill">$(_html_escape(string(sha)))</span>
      <span class="pill">$(_html_escape(key_text))</span>
    </div>
  </form>
</section>"""
end

function _studio_validation_issues_html(issues; id = nothing)
    issues isa AbstractVector || return ""
    isempty(issues) && return ""

    cards = [_studio_validation_issue_html(issue) for issue in issues]
    id_attr = isnothing(id) ? "" : " id=\"$(_html_escape(string(id)))\""
    return """
<div class="validation-issues"$id_attr aria-label="Validation issues">
      $(join(cards, "\n      "))
    </div>"""
end

function _studio_validation_issue_html(issue)
    issue isa AbstractDict ||
        return "<article class=\"validation-issue\"><p class=\"issue-title\">$(_html_escape(string(issue)))</p></article>"
    severity = string(get(issue, "severity", "unknown"))
    blocking = get(issue, "blocking", false) === true
    field = get(issue, "field", nothing)
    title = "$(severity)$(blocking ? " (blocking)" : "")"
    details = String[
        "<p class=\"issue-title\">$(_html_escape(title)): $(_html_escape(string(get(issue, "message", ""))))</p>",
    ]
    for (label, key) in (
        ("Field", "field"),
        ("YAML path", "yaml_path"),
        ("Physical implication", "physical_implication"),
        ("Suggested fix", "suggested_fix"),
        ("Documentation", "documentation"),
        ("Remediation action", "remediation_action"),
    )
        value = get(issue, key, nothing)
        value isa AbstractString && !isempty(value) || continue
        push!(
            details,
            "<p><strong>$label:</strong> $(_html_escape(value))</p>",
        )
    end
    class = "validation-issue$(blocking ? " blocking" : "")"
    return """
<article class="$class" data-validation-severity="$(_html_escape(severity))" data-validation-field="$(_html_escape(string(field)))">
        $(join(details, "\n        "))
      </article>"""
end

function _studio_input_summary_cards_html(summary::AbstractDict)
    keys = ["files", "editable", "parse_errors", "text_included"]
    labels = ["Files", "Editable", "Parse Errors", "Text Included"]
    cards = [
        "<p><strong>$(_html_escape(label)):</strong> $(_html_escape(string(get(summary, key, 0))))</p>"
        for (key, label) in zip(keys, labels)
    ]
    return join(cards, "\n      ")
end

function _studio_editor_artifact_links_html(workbench_href, inputs_href)
    links = String[]
    if workbench_href isa AbstractString && !isempty(workbench_href)
        push!(links, "<a href=\"$(_html_escape(workbench_href))\">Workbench</a>")
    end
    if inputs_href isa AbstractString && !isempty(inputs_href)
        push!(links, "<a href=\"$(_html_escape(inputs_href))\">Inputs YAML</a>")
    end
    isempty(links) && return ""
    return "<h3>Artifacts</h3>\n      <p>" * join(links, " &middot; ") * "</p>"
end

function _html_escape(value::AbstractString)
    escaped = replace(value, "&" => "&amp;")
    escaped = replace(escaped, "<" => "&lt;")
    escaped = replace(escaped, ">" => "&gt;")
    escaped = replace(escaped, "\"" => "&quot;")
    return replace(escaped, "'" => "&#39;")
end

function _studio_project_value(value)
    if value === missing
        return nothing
    elseif value isa AbstractDict
        normalized = OrderedCollections.OrderedDict{String,Any}()
        keys_iter =
            value isa OrderedCollections.OrderedDict ? collect(keys(value)) :
            sort(collect(keys(value)); by = string)
        for key in keys_iter
            normalized[string(key)] = _studio_project_value(value[key])
        end
        return normalized
    elseif value isa AbstractVector || value isa Tuple
        return [_studio_project_value(item) for item in value]
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
