export build_studio_project,
    write_studio_project,
    read_studio_project,
    studio_project_issues,
    validate_studio_project,
    studio_project_health,
    studio_project_generated_script_path,
    read_studio_project_generated_script,
    render_studio_home_html,
    write_studio_home_html,
    render_studio_workbench_html,
    write_studio_workbench_html,
    write_studio_workbench_bundle

const STUDIO_PROJECT_SCHEMA_VERSION = "owens-studio-project/v1"
const STUDIO_WORKBENCH_SCHEMA_VERSION = "owens-studio-workbench/v1"
const STUDIO_WORKBENCH_BUNDLE_SCHEMA_VERSION = "owens-studio-bundle/v1"
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
    status =
        isempty(project_issues) && all(row["status"] == "ok" for row in rows) ? "ok" :
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
    render_studio_workbench_html(project_or_health; health_href=nothing, script_href=nothing, open_href=nothing)

Render a dependency-light OWENS Studio workbench shell as static HTML. This is
the first GUI slice: a project health view that later Genie routes can serve
without changing the underlying health API.
"""
function render_studio_workbench_html(
    project_or_health;
    health_href = nothing,
    script_href = nothing,
    open_href = nothing,
)
    health = _studio_health_input(project_or_health)
    title = _html_escape(string(get(health, "name", "OWENS Studio")))
    status = _html_escape(string(health["status"]))

    return """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>$title - OWENS Studio</title>
  <style>
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
      grid-template-columns: 220px minmax(0, 1fr) 320px;
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
    nav a {
      display: block;
      padding: 9px 10px;
      color: var(--ink);
      text-decoration: none;
      border-radius: 6px;
      font-size: 14px;
    }
    nav a.active {
      background: #e7f0f8;
      color: var(--blue);
      font-weight: 650;
    }
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
    @media (max-width: 980px) {
      .studio { grid-template-columns: 1fr; }
      nav, aside { border: 0; border-bottom: 1px solid var(--line); }
      .cards { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
    }
  </style>
</head>
<body>
  <div class="studio">
    <nav>
      <h1>OWENS Studio</h1>
      $(_studio_nav_html())
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
      $(_studio_artifact_links_html(health_href, script_href, open_href))
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
    return """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>OWENS Studio</title>
  <style>
    :root {
      color-scheme: light;
      --ink: #17202a;
      --muted: #5c6670;
      --line: #d9dee4;
      --panel: #f6f8fa;
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
      grid-template-columns: 220px minmax(0, 1fr);
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
    nav a {
      display: block;
      padding: 9px 10px;
      color: var(--ink);
      text-decoration: none;
      border-radius: 6px;
      font-size: 14px;
    }
    nav a.active {
      background: #e7f0f8;
      color: var(--blue);
      font-weight: 650;
    }
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
    @media (max-width: 840px) {
      .studio { grid-template-columns: 1fr; }
      nav { border-right: 0; border-bottom: 1px solid var(--line); }
    }
  </style>
</head>
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
    health_file = joinpath(bundle_dir, "health.yml")
    YAML.write_file(health_file, health)

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
    )
    open(index_file, "w") do io
        write(io, html)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_WORKBENCH_BUNDLE_SCHEMA_VERSION,
        "bundle_dir" => bundle_dir,
        "project_file" => abspath(project_path),
        "index_html" => index_file,
        "health_file" => health_file,
        "script_file" => script_file,
        "project_status" => health["status"],
        "bytes" => OrderedCollections.OrderedDict{String,Any}(
            "index_html" => stat(index_file).size,
            "health_file" => stat(health_file).size,
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

function _studio_nav_html(active_item::AbstractString = "Project")
    items = [
        "Home",
        "Project",
        "Geometry",
        "Airfoils",
        "Structure",
        "Environment",
        "Controls",
        "Simulation",
        "Validation",
        "Results",
        "Reports",
    ]
    return join(
        [
            "<a class=\"$(label == active_item ? "active" : "")\" href=\"#\">$(_html_escape(label))</a>"
            for label in items
        ],
        "\n      ",
    )
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
    issues = get(row, "issues", String[])
    issue_text = isempty(issues) ? "" : join(string.(issues), "; ")
    return "<tr><td>$(_html_escape(string(row["status"])))</td><td>$(_html_escape(string(get(row, "role", ""))))</td><td>$(_html_escape(string(get(row, "path", ""))))</td><td>$(_html_escape(issue_text))</td></tr>"
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

function _studio_artifact_links_html(health_href, script_href, open_href)
    links = String[]
    if open_href isa AbstractString && !isempty(open_href)
        push!(links, "<a href=\"$(_html_escape(open_href))\">Open Payload</a>")
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
