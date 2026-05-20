export file_sha256,
    file_provenance,
    collect_git_state,
    collect_package_versions,
    build_run_manifest,
    write_run_manifest,
    read_run_manifest

const RUN_MANIFEST_SCHEMA_VERSION = "owens-run-manifest/v1"

"""
    file_sha256(path)

Return the SHA-256 digest for `path` as a lowercase hexadecimal string.
"""
function file_sha256(path::AbstractString)
    isfile(path) || throw(ArgumentError("Cannot hash missing file: $path"))

    open(path, "r") do io
        return bytes2hex(SHA.sha256(io))
    end
end

"""
    file_provenance(path; root=pwd(), role=nothing)

Build a deterministic manifest record for a file. The returned record includes
the file path relative to `root`, byte count, and SHA-256 digest.
"""
function file_provenance(path::AbstractString; root::AbstractString = pwd(), role = nothing)
    isfile(path) || throw(ArgumentError("Cannot record missing file: $path"))

    record = OrderedCollections.OrderedDict{String,Any}()
    record["path"] = relpath(abspath(path), abspath(root))
    record["bytes"] = stat(path).size
    record["sha256"] = file_sha256(path)
    if !isnothing(role)
        record["role"] = string(role)
    end

    return record
end

"""
    collect_package_versions([modules])

Return package versions for the OWENS modules loaded in the current process.
"""
function collect_package_versions(
    modules::Vector{Module} = Module[
        @__MODULE__,
        OWENSAero,
        OWENSFEA,
        OWENSOpenFASTWrappers,
        OWENSPreComp,
    ],
)
    packages = OrderedCollections.OrderedDict{String,Any}()
    for module_ref in modules
        packages[string(nameof(module_ref))] = _module_version(module_ref)
    end

    return packages
end

"""
    collect_git_state(root=pwd())

Return best-effort git metadata for the repository containing `root`.
"""
function collect_git_state(root::AbstractString = pwd())
    repo_root = _git_repo_root(root)
    git_state = OrderedCollections.OrderedDict{String,Any}()

    if isnothing(repo_root)
        git_state["available"] = false
        git_state["root"] = abspath(root)
        git_state["branch"] = nothing
        git_state["commit"] = nothing
        git_state["dirty"] = nothing
        return git_state
    end

    git_dir = _git_dir(repo_root)
    branch, commit = _git_head(git_dir)
    git_state["available"] = !isnothing(git_dir)
    git_state["root"] = repo_root
    git_state["branch"] = branch
    git_state["commit"] = commit
    git_state["dirty"] = _git_dirty(repo_root)

    return git_state
end

"""
    build_run_manifest(; kwargs...)

Build an OWENS run manifest with hashed input and output files, Julia/package
versions, and best-effort git metadata. This is intentionally opt-in so solver
workflows can adopt it incrementally.
"""
function build_run_manifest(;
    run_id = nothing,
    run_name::AbstractString = "unnamed",
    project_root::AbstractString = pwd(),
    solver::AbstractString = "unspecified",
    project_file = nothing,
    modeling_options_file = nothing,
    windio_file = nothing,
    input_files = String[],
    output_files = String[],
    generated_files = String[],
    parameters = OrderedCollections.OrderedDict{String,Any}(),
    metadata = OrderedCollections.OrderedDict{String,Any}(),
    warnings = String[],
    status::AbstractString = "created",
    created_at_utc = nothing,
)
    root = abspath(project_root)
    created = isnothing(created_at_utc) ? _utc_timestamp() : string(created_at_utc)
    manifest_run_id = isnothing(run_id) ? _default_run_id() : string(run_id)

    inputs = OrderedCollections.OrderedDict{String,Any}[]
    _push_file_record!(inputs, project_file, root, "project")
    _push_file_record!(inputs, modeling_options_file, root, "modeling_options")
    _push_file_record!(inputs, windio_file, root, "windio")
    for input_file in input_files
        _push_file_record!(inputs, input_file, root, "input")
    end

    outputs = OrderedCollections.OrderedDict{String,Any}[]
    for output_file in output_files
        _push_file_record!(outputs, output_file, root, "output")
    end

    generated = OrderedCollections.OrderedDict{String,Any}[]
    for generated_file in generated_files
        _push_file_record!(generated, generated_file, root, "generated")
    end

    julia_state = OrderedCollections.OrderedDict{String,Any}(
        "version" => string(VERSION),
        "kernel" => string(Sys.KERNEL),
        "machine" => string(Sys.MACHINE),
        "threads" => Base.Threads.nthreads(),
    )

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => RUN_MANIFEST_SCHEMA_VERSION,
        "run_id" => manifest_run_id,
        "run_name" => string(run_name),
        "created_at_utc" => created,
        "project_root" => root,
        "solver" => string(solver),
        "status" => string(status),
        "julia" => julia_state,
        "packages" => collect_package_versions(),
        "git" => collect_git_state(root),
        "parameters" => _manifest_value(parameters),
        "metadata" => _manifest_value(metadata),
        "inputs" => inputs,
        "outputs" => outputs,
        "generated" => generated,
        "warnings" => [string(warning) for warning in warnings],
    )
end

"""
    write_run_manifest(path, manifest)
    write_run_manifest(path; kwargs...)

Write an existing manifest, or build and write one from `build_run_manifest`
keyword arguments. Returns the manifest that was written.
"""
function write_run_manifest(path::AbstractString, manifest::AbstractDict)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    normalized_manifest = _manifest_value(manifest)
    YAML.write_file(path, normalized_manifest)
    return normalized_manifest
end

function write_run_manifest(path::AbstractString; kwargs...)
    manifest = build_run_manifest(; kwargs...)
    write_run_manifest(path, manifest)
    return manifest
end

"""
    read_run_manifest(path)

Read a YAML run manifest using string-keyed ordered dictionaries.
"""
function read_run_manifest(path::AbstractString)
    isfile(path) || throw(ArgumentError("Cannot read missing manifest: $path"))
    return YAML.load_file(path; dicttype = OrderedCollections.OrderedDict{String,Any})
end

function _push_file_record!(records, path, root::AbstractString, role::AbstractString)
    if !isnothing(path)
        push!(records, file_provenance(string(path); root = root, role = role))
    end

    return records
end

function _manifest_value(value)
    if value isa AbstractDict
        normalized = OrderedCollections.OrderedDict{String,Any}()
        keys_iter =
            value isa OrderedCollections.OrderedDict ? collect(keys(value)) :
            sort(collect(keys(value)); by = string)
        for key in keys_iter
            normalized[string(key)] = _manifest_value(value[key])
        end
        return normalized
    elseif value isa AbstractVector || value isa Tuple
        return [_manifest_value(item) for item in value]
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

function _module_version(module_ref::Module)
    version = try
        Base.pkgversion(module_ref)
    catch
        nothing
    end

    return isnothing(version) ? nothing : string(version)
end

function _utc_timestamp()
    return Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS.sss") * "Z"
end

function _default_run_id()
    return "run-" * Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSSsss")
end

function _git_repo_root(path::AbstractString)
    root = abspath(isdir(path) ? path : dirname(path))
    while true
        git_marker = joinpath(root, ".git")
        if isdir(git_marker) || isfile(git_marker)
            return root
        end

        parent = dirname(root)
        parent == root && return nothing
        root = parent
    end
end

function _git_dir(repo_root::AbstractString)
    git_marker = joinpath(repo_root, ".git")
    if isdir(git_marker)
        return git_marker
    elseif isfile(git_marker)
        git_dir_line = strip(read(git_marker, String))
        startswith(git_dir_line, "gitdir:") || return nothing
        git_dir = strip(replace(git_dir_line, r"^gitdir:\s*" => ""))
        resolved_git_dir =
            isabspath(git_dir) ? git_dir : normpath(joinpath(repo_root, git_dir))
        return isdir(resolved_git_dir) ? resolved_git_dir : nothing
    end

    return nothing
end

function _git_head(git_dir)
    isnothing(git_dir) && return nothing, nothing

    head_path = joinpath(git_dir, "HEAD")
    isfile(head_path) || return nothing, nothing

    head = strip(read(head_path, String))
    if startswith(head, "ref: ")
        ref = strip(head[6:end])
        branch = replace(ref, r"^refs/heads/" => "")
        commit = _git_ref_commit(git_dir, ref)
        return branch, commit
    end

    return nothing, isempty(head) ? nothing : head
end

function _git_ref_commit(git_dir::AbstractString, ref::AbstractString)
    ref_path = joinpath(git_dir, ref)
    if isfile(ref_path)
        commit = strip(read(ref_path, String))
        return isempty(commit) ? nothing : commit
    end

    return _packed_ref_commit(git_dir, ref)
end

function _packed_ref_commit(git_dir::AbstractString, ref::AbstractString)
    packed_refs_path = joinpath(git_dir, "packed-refs")
    isfile(packed_refs_path) || return nothing

    for line in eachline(packed_refs_path)
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#") || startswith(stripped, "^")
            continue
        end

        parts = split(stripped)
        if length(parts) == 2 && parts[2] == ref
            return parts[1]
        end
    end

    return nothing
end

function _git_dirty(repo_root::AbstractString)
    git_exe = Sys.which("git")
    isnothing(git_exe) && return nothing

    try
        status = read(
            Cmd(
                [git_exe, "-C", repo_root, "status", "--porcelain", "--untracked-files=no"];
                ignorestatus = true,
            ),
            String,
        )
        return !isempty(strip(status))
    catch
        return nothing
    end
end
