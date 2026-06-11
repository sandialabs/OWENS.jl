export studio_gui_capability_catalog,
    studio_gui_capability_ids,
    studio_gui_quality_gate_catalog,
    studio_gui_quality_gate_ids

const STUDIO_GUI_CAPABILITY_SCHEMA_VERSION = "owens-studio-capability-catalog/v1"
const STUDIO_GUI_QUALITY_GATE_SCHEMA_VERSION = "owens-studio-quality-gates/v1"

const STUDIO_GUI_CAPABILITIES = OrderedCollections.OrderedDict{String,Any}(
    "project_manifest_health" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Project manifests and health",
        "tier" => "mvp",
        "status" => "implemented",
        "description" => "Open OWENS Studio projects, verify tracked file provenance, and summarize run-manifest health.",
        "surfaces" => ["OWENS", "OWENS_APP", "static_workbench"],
    ),
    "template_catalog" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Project templates and examples",
        "tier" => "mvp",
        "status" => "implemented",
        "description" => "List built-in project templates and committed example projects for the project gallery.",
        "surfaces" => ["OWENS", "OWENS_APP", "static_home"],
    ),
    "generated_script_export" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Generated Julia script export",
        "tier" => "mvp",
        "status" => "implemented",
        "description" => "Expose the reproducible Julia driver attached to a Studio project or bundle.",
        "surfaces" => ["OWENS", "OWENS_APP", "static_bundle"],
    ),
    "input_file_summary" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Editable input-file summary",
        "tier" => "mvp",
        "status" => "implemented",
        "description" => "Summarize project input files, YAML top-level keys, provenance, and optional editor text.",
        "surfaces" => ["OWENS", "OWENS_APP"],
    ),
    "route_contracts" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Versioned app route contracts",
        "tier" => "mvp",
        "status" => "implemented",
        "description" => "Return dependency-light route metadata that a local web shell can wrap without duplicating OWENS semantics.",
        "surfaces" => ["OWENS_APP"],
    ),
    "local_backend_server" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Launchable local backend server",
        "tier" => "mvp",
        "status" => "planned",
        "description" => "Run a persistent local HTTP/WebSocket backend over the tested service contracts.",
        "surfaces" => ["OWENS_APP"],
    ),
    "model_editors" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Modeling-options and WindIO editors",
        "tier" => "mvp",
        "status" => "planned",
        "description" => "Provide schema-aware edit/apply/revert flows for OWENS modeling options and WindIO files.",
        "surfaces" => ["frontend", "OWENS_APP"],
    ),
    "geometry_editor" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Interactive turbine geometry editor",
        "tier" => "research",
        "status" => "planned",
        "description" => "Edit blades, towers, struts, hubs, nacelles, platforms, local frames, and mesh connectivity.",
        "surfaces" => ["frontend", "OWENS"],
    ),
    "airfoil_polar_manager" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Airfoil and polar manager",
        "tier" => "research",
        "status" => "planned",
        "description" => "Import, plot, validate, extrapolate, and assign airfoil coordinates and polar families.",
        "surfaces" => ["frontend", "OWENSAero", "OWENS_APP"],
    ),
    "structural_editor" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Structural and composites editor",
        "tier" => "research",
        "status" => "planned",
        "description" => "Edit materials, laminates, sectional properties, offsets, and structural component assignments.",
        "surfaces" => ["frontend", "OWENSPreComp", "OWENSFEA"],
    ),
    "simulation_setup" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Simulation setup panels",
        "tier" => "mvp",
        "status" => "planned",
        "description" => "Expose backend-specific aerodynamic, structural, controls, wind, hydrodynamic, and restart settings.",
        "surfaces" => ["frontend", "OWENS_APP", "OWENS"],
    ),
    "hawt_aeroelastic_workflow" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "HAWT aeroelastic workflow",
        "tier" => "research",
        "status" => "experimental",
        "description" => "Label HAWT setup as an experimental, validation-gated workflow until setupOWENS, frames, controls, and backend parity are complete.",
        "surfaces" => [
            "frontend",
            "OWENS_APP",
            "OWENS",
            "OWENSAero",
            "OWENSOpenFASTWrappers",
        ],
    ),
    "run_queue" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Run queue and log streaming",
        "tier" => "mvp",
        "status" => "planned",
        "description" => "Launch, cancel, restart, clone, compare, and inspect runs while streaming progress and logs.",
        "surfaces" => ["frontend", "OWENS_APP"],
    ),
    "visualization" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "3D and plot visualization",
        "tier" => "research",
        "status" => "planned",
        "description" => "Render geometry, meshes, deformation, loads, wake data, time histories, polars, and result comparisons.",
        "surfaces" => ["frontend", "OWENS"],
    ),
    "postprocessing_reports" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Postprocessing and reports",
        "tier" => "research",
        "status" => "planned",
        "description" => "Browse result channels, compare runs, apply validation tolerances, and generate reproducible reports.",
        "surfaces" => ["frontend", "OWENS"],
    ),
    "packaging" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Distribution and dependency checks",
        "tier" => "mvp",
        "status" => "planned",
        "description" => "Package or launch OWENS Studio reliably across supported platforms with solver dependency diagnostics.",
        "surfaces" => ["OWENS_APP", "ci"],
    ),
)

const STUDIO_GUI_QUALITY_GATES = OrderedCollections.OrderedDict{String,Any}(
    "service_contract_tests" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Service and route contracts",
        "status" => "implemented",
        "required_for_done" => true,
        "description" => "Every Studio feature exposes a versioned service or route payload and has Julia-level tests for success and failure cases.",
        "evidence_required" => [
            "schema version",
            "request/response tests",
            "structured error tests",
        ],
        "current_evidence" => [
            "owens-studio-route-catalog/v1",
            "OWENS Studio app service tests",
            "OWENS Studio app route handler tests",
        ],
    ),
    "file_provenance" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "File provenance and stale-state handling",
        "status" => "implemented",
        "required_for_done" => true,
        "description" => "File-writing and run-artifact workflows must record hashes, detect stale inputs/outputs, and report actionable remediation.",
        "evidence_required" => [
            "before/after SHA-256",
            "atomic write or explicit nonwrite path",
            "stale-file tests",
        ],
        "current_evidence" => [
            "owens-studio-input-save-provenance/v1",
            "owens-run-artifact-remediation/v1",
            "Studio stale input and missing output tests",
        ],
    ),
    "generated_script_equivalence" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Generated script equivalence",
        "status" => "partial",
        "required_for_done" => true,
        "description" => "GUI-driven project changes must regenerate the same solver inputs as the corresponding headless fixture.",
        "evidence_required" => [
            "headless fixture",
            "generated script diff",
            "run manifest comparison",
        ],
        "current_evidence" => [
            "RM2 static fixture generated script export",
            "project bundle script hash checks",
        ],
    ),
    "physical_validation" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Physical and numerical validation",
        "status" => "planned",
        "required_for_done" => true,
        "description" => "Feature tests must check units, sign conventions, physical bounds, and solver-specific invariants instead of string presence only.",
        "evidence_required" => [
            "unit-aware validation",
            "physical invariant tests",
            "backend-specific reduced cases",
        ],
        "current_evidence" => String[],
    ),
    "browser_interaction_tests" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Browser interaction tests",
        "status" => "planned",
        "required_for_done" => true,
        "description" => "Interactive Studio workflows must be tested in a browser or equivalent UI harness, including form submission and state transitions.",
        "evidence_required" => [
            "browser launch",
            "form submission",
            "state transition checks",
        ],
        "current_evidence" => String[],
    ),
    "visual_regression" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Visual and responsive regression",
        "status" => "planned",
        "required_for_done" => true,
        "description" => "Workbench, editor, plots, and future 3D views need screenshot or visual regression checks across desktop and mobile viewports.",
        "evidence_required" => [
            "desktop screenshot",
            "mobile screenshot",
            "nonblank plot or canvas checks",
        ],
        "current_evidence" => String[],
    ),
    "accessibility_keyboard" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Accessibility and keyboard usability",
        "status" => "partial",
        "required_for_done" => true,
        "description" => "Controls must have labels, focus states, keyboard paths, ARIA where appropriate, and color-independent status cues.",
        "evidence_required" => [
            "label/ARIA tests",
            "keyboard navigation tests",
            "contrast checks",
        ],
        "current_evidence" => [
            "editor labels",
            "aria-describedby/aria-invalid tests",
            "focus CSS tests",
        ],
    ),
    "documentation_update" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Documentation and known limitations",
        "status" => "partial",
        "required_for_done" => true,
        "description" => "Every user-facing Studio feature must document launch/use steps, limitations, troubleshooting, and capability status.",
        "evidence_required" => [
            "quick start",
            "capability status",
            "known limitations",
        ],
        "current_evidence" => [
            "Studio diagnostics quick start",
            "capability catalog",
            "todo.md limitation tracking",
        ],
    ),
    "performance_budget" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Performance and scale budget",
        "status" => "partial",
        "required_for_done" => true,
        "description" => "Large project inputs, result files, and many-run manifests need explicit size limits, lazy loading, and performance tests.",
        "evidence_required" => [
            "size limit",
            "large fixture",
            "pagination or lazy loading",
        ],
        "current_evidence" => [
            "max_text_bytes parse and validation guard",
        ],
    ),
)

"""
    studio_gui_capability_catalog()

Return the versioned OWENS Studio capability catalog used by the app shell to
show what is implemented, planned, or intentionally gated behind solver work.
"""
function studio_gui_capability_catalog()
    rows = OrderedCollections.OrderedDict{String,Any}[]
    for (id, metadata) in STUDIO_GUI_CAPABILITIES
        row = OrderedCollections.OrderedDict{String,Any}("id" => id)
        for (key, value) in metadata
            row[key] = _studio_project_value(value)
        end
        push!(rows, row)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_GUI_CAPABILITY_SCHEMA_VERSION,
        "capabilities" => rows,
    )
end

"""
    studio_gui_capability_ids()

Return GUI capability IDs in deterministic display order.
"""
studio_gui_capability_ids() = collect(keys(STUDIO_GUI_CAPABILITIES))

"""
    studio_gui_quality_gate_catalog()

Return the versioned quality gates that define when an OWENS Studio feature can
be called done. The catalog is intentionally stricter than the current static
HTML shell so partial features remain visibly incomplete.
"""
function studio_gui_quality_gate_catalog()
    rows = OrderedCollections.OrderedDict{String,Any}[]
    for (id, metadata) in STUDIO_GUI_QUALITY_GATES
        row = OrderedCollections.OrderedDict{String,Any}("id" => id)
        for (key, value) in metadata
            row[key] = _studio_project_value(value)
        end
        push!(rows, row)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_GUI_QUALITY_GATE_SCHEMA_VERSION,
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "gates" => length(rows),
            "implemented" => count(row -> row["status"] == "implemented", rows),
            "partial" => count(row -> row["status"] == "partial", rows),
            "planned" => count(row -> row["status"] == "planned", rows),
            "required_for_done" => count(row -> row["required_for_done"] === true, rows),
        ),
        "gates" => rows,
    )
end

"""
    studio_gui_quality_gate_ids()

Return Studio quality gate IDs in deterministic display order.
"""
studio_gui_quality_gate_ids() = collect(keys(STUDIO_GUI_QUALITY_GATES))
