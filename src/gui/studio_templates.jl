export studio_project_template_catalog,
    studio_project_template_names,
    studio_example_project_catalog,
    studio_example_project_names,
    create_studio_project_template

const STUDIO_TEMPLATE_CATALOG_SCHEMA_VERSION = "owens-studio-template-catalog/v1"
const STUDIO_EXAMPLE_CATALOG_SCHEMA_VERSION = "owens-studio-example-catalog/v1"
const STUDIO_PROJECT_TEMPLATES = OrderedCollections.OrderedDict{String,Any}(
    "blank" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "Blank OWENS Studio Project",
        "description" => "Blank OWENS Studio project",
        "turbine_type" => "custom",
        "solver_path" => nothing,
        "creates_generated_script" => false,
        "creates_run_manifest" => false,
    ),
    "rm2" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "RM2 VAWT Template",
        "description" => "RM2 VAWT WindIO project",
        "turbine_type" => "VAWT",
        "solver_path" => "runOWENSWINDIO",
        "creates_generated_script" => true,
        "creates_run_manifest" => true,
    ),
)
const STUDIO_EXAMPLE_PROJECTS = OrderedCollections.OrderedDict{String,Any}(
    "rm2" => OrderedCollections.OrderedDict{String,Any}(
        "title" => "RM2 GUI Fixture",
        "description" => "Relocatable RM2 VAWT project for OWENS Studio smoke tests.",
        "turbine_type" => "VAWT",
        "template" => "rm2",
        "project_relative_path" =>
            joinpath("examples", "gui", "rm2", "owens_project.yml"),
    ),
)

"""
    studio_project_template_catalog()

Return a stable catalog of built-in OWENS Studio templates for GUI project
creation controls.
"""
function studio_project_template_catalog()
    templates = OrderedCollections.OrderedDict{String,Any}[]
    for (template, metadata) in STUDIO_PROJECT_TEMPLATES
        row = OrderedCollections.OrderedDict{String,Any}("template" => template)
        for (key, value) in metadata
            row[key] = value
        end
        push!(templates, row)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_TEMPLATE_CATALOG_SCHEMA_VERSION,
        "templates" => templates,
    )
end

"""
    studio_project_template_names()

Return the built-in OWENS Studio project templates in deterministic order.
"""
studio_project_template_names() = collect(keys(STUDIO_PROJECT_TEMPLATES))

"""
    studio_example_project_catalog()

Return committed OWENS Studio example projects that the GUI can show in an
example gallery or open directly for smoke tests.
"""
function studio_example_project_catalog()
    examples = OrderedCollections.OrderedDict{String,Any}[]
    repo_root = normpath(joinpath(module_path, ".."))
    for (example, metadata) in STUDIO_EXAMPLE_PROJECTS
        relative_path = metadata["project_relative_path"]
        project_file = normpath(joinpath(repo_root, relative_path))
        row = OrderedCollections.OrderedDict{String,Any}(
            "example" => example,
            "project_file" => project_file,
            "project_relative_path" => relative_path,
            "available" => isfile(project_file),
        )
        for (key, value) in metadata
            key == "project_relative_path" && continue
            row[key] = value
        end
        push!(examples, row)
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => STUDIO_EXAMPLE_CATALOG_SCHEMA_VERSION,
        "examples" => examples,
    )
end

"""
    studio_example_project_names()

Return committed OWENS Studio example project names in deterministic order.
"""
studio_example_project_names() = collect(keys(STUDIO_EXAMPLE_PROJECTS))

"""
    create_studio_project_template(target; template="rm2", overwrite=false, created_at_utc=nothing)

Create a file-backed OWENS Studio project template under `target` and return
the generated project, run manifest, and script paths. The RM2 template copies
the existing RM2 modeling-options and WindIO inputs, writes a reproducible
`runOWENSWINDIO` Julia script, records it in a run manifest, and creates an
`owens_project.yml` manifest that the GUI health panel can open immediately.
"""
function create_studio_project_template(
    target::AbstractString;
    template::AbstractString = "rm2",
    overwrite::Bool = false,
    created_at_utc = nothing,
)
    template_key = lowercase(string(template))
    haskey(STUDIO_PROJECT_TEMPLATES, template_key) ||
        throw(ArgumentError("Unknown OWENS Studio template: $template"))

    project_root = abspath(target)
    _prepare_studio_template_root(project_root, overwrite)

    if template_key == "blank"
        return _create_blank_studio_project_template(project_root; created_at_utc)
    elseif template_key == "rm2"
        return _create_rm2_studio_project_template(project_root; overwrite, created_at_utc)
    end
end

function _create_blank_studio_project_template(project_root; created_at_utc)
    mkpath(joinpath(project_root, "inputs"))
    mkpath(joinpath(project_root, "runs"))
    project_file = joinpath(project_root, "owens_project.yml")
    project = write_studio_project(
        project_file,
        project_root;
        project_id = "blank",
        name = "Blank OWENS Studio Project",
        description = "Blank project scaffold for building a new OWENS model.",
        metadata = OrderedCollections.OrderedDict{String,Any}(
            "template" => "blank",
            "template_description" => _studio_template_description("blank"),
        ),
        created_at_utc,
    )

    return OrderedCollections.OrderedDict{String,Any}(
        "template" => "blank",
        "project_file" => project_file,
        "project" => project,
        "run_manifest_file" => nothing,
        "run_manifest" => nothing,
        "script_file" => nothing,
        "script" => nothing,
    )
end

function _create_rm2_studio_project_template(project_root; overwrite, created_at_utc)
    template_root = normpath(joinpath(module_path, "..", "examples", "RM2"))
    source_model = joinpath(template_root, "modeling_options_OWENS_RM2.yml")
    source_windio = joinpath(template_root, "WINDIO_RM2.yaml")
    isfile(source_model) ||
        throw(ArgumentError("RM2 template modeling-options file is missing: $source_model"))
    isfile(source_windio) ||
        throw(ArgumentError("RM2 template WindIO file is missing: $source_windio"))

    inputs_dir = joinpath(project_root, "inputs")
    run_dir = joinpath(project_root, "runs", "rm2")
    mkpath(inputs_dir)
    mkpath(run_dir)

    model_file = joinpath(inputs_dir, "modeling_options_OWENS_RM2.yml")
    windio_file = joinpath(inputs_dir, "WINDIO_RM2.yaml")
    _copy_studio_template_file(source_model, model_file, overwrite)
    _copy_studio_template_file(source_windio, windio_file, overwrite)

    spec = windio_run_spec(model_file, windio_file, run_dir)
    script_file = joinpath(run_dir, "run_rm2_windio.jl")
    script = write_windio_run_script(script_file, spec)
    manifest_file = joinpath(run_dir, "run_manifest.yml")
    manifest = write_windio_run_manifest(
        manifest_file,
        spec;
        run_id = "rm2-template",
        run_name = "RM2 Studio Template",
        generated_files = [script_file],
        metadata = OrderedCollections.OrderedDict{String,Any}(
            "template" => "rm2",
            "template_description" => _studio_template_description("rm2"),
        ),
        status = "created",
        created_at_utc,
    )

    project_file = joinpath(project_root, "owens_project.yml")
    project = write_studio_project(
        project_file,
        project_root;
        project_id = "rm2",
        name = "RM2 VAWT Template",
        description = "RM2 VAWT WindIO template for OWENS Studio GUI workflows.",
        modeling_options_file = model_file,
        windio_file,
        run_manifests = [manifest_file],
        metadata = OrderedCollections.OrderedDict{String,Any}(
            "template" => "rm2",
            "template_description" => _studio_template_description("rm2"),
            "generated_script" => relpath(script_file, project_root),
        ),
        created_at_utc,
    )

    return OrderedCollections.OrderedDict{String,Any}(
        "template" => "rm2",
        "project_file" => project_file,
        "project" => project,
        "run_manifest_file" => manifest_file,
        "run_manifest" => manifest,
        "script_file" => script_file,
        "script" => script,
    )
end

_studio_template_description(template::AbstractString) =
    string(STUDIO_PROJECT_TEMPLATES[string(template)]["description"])

function _prepare_studio_template_root(project_root::AbstractString, overwrite::Bool)
    if isdir(project_root)
        if !overwrite && !isempty(readdir(project_root))
            throw(ArgumentError("Target project directory is not empty: $project_root"))
        end
    else
        mkpath(project_root)
    end

    return project_root
end

function _copy_studio_template_file(
    source::AbstractString,
    destination::AbstractString,
    overwrite::Bool,
)
    if isfile(destination) && !overwrite
        throw(ArgumentError("Template destination already exists: $destination"))
    end

    cp(source, destination; force = overwrite)
    return destination
end
