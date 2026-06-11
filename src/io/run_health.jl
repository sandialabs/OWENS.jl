export verify_file_provenance, run_manifest_health

const RUN_MANIFEST_HEALTH_SCHEMA_VERSION = "owens-run-health/v1"
const RUN_ARTIFACT_REMEDIATION_SCHEMA_VERSION = "owens-run-artifact-remediation/v1"
const RUN_ARTIFACT_REMEDIATION_DOC = "docs/src/getting_started.md"

"""
    verify_file_provenance(record; root=pwd())

Compare a manifest file record against the current file on disk. Relative paths
are resolved against `root`; absolute paths are used as-is. The returned row is
intended for project-health panels and reports missing, modified, or malformed
artifacts without throwing for ordinary file drift.
"""
function verify_file_provenance(record; root::AbstractString = pwd())
    issues = _file_provenance_record_issues(record)
    if !isempty(issues)
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "invalid_record",
            "issues" => issues,
            "path" => _record_get_string(record, "path"),
            "resolved_path" => nothing,
            "role" => _record_get_string(record, "role"),
            "expected_bytes" => _record_get_value(record, "bytes"),
            "actual_bytes" => nothing,
            "expected_sha256" => _record_get_string(record, "sha256"),
            "actual_sha256" => nothing,
        )
    end

    resolved_path = _resolve_manifest_file_path(record["path"], root)
    if !isfile(resolved_path)
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "missing",
            "issues" => String["missing file"],
            "path" => string(record["path"]),
            "resolved_path" => resolved_path,
            "role" => get(record, "role", nothing),
            "expected_bytes" => record["bytes"],
            "actual_bytes" => nothing,
            "expected_sha256" => record["sha256"],
            "actual_sha256" => nothing,
        )
    end

    actual_bytes = stat(resolved_path).size
    actual_sha256 = file_sha256(resolved_path)
    issues = String[]
    if actual_bytes != record["bytes"]
        push!(issues, "bytes mismatch")
    end
    if actual_sha256 != record["sha256"]
        push!(issues, "sha256 mismatch")
    end

    return OrderedCollections.OrderedDict{String,Any}(
        "status" => isempty(issues) ? "ok" : "modified",
        "issues" => issues,
        "path" => string(record["path"]),
        "resolved_path" => resolved_path,
        "role" => get(record, "role", nothing),
        "expected_bytes" => record["bytes"],
        "actual_bytes" => actual_bytes,
        "expected_sha256" => record["sha256"],
        "actual_sha256" => actual_sha256,
    )
end

"""
    run_manifest_health(manifest_or_path; root=nothing, summarize_outputs=true)

Return a run health summary for a run manifest, including manifest schema
diagnostics, current file-provenance status for inputs/outputs/generated files,
and optional `outputData` channel summaries for healthy HDF5 outputs.
Relative `project_root` values are resolved against the manifest file location
when a manifest path is supplied.
"""
run_manifest_health(path::AbstractString; kwargs...) =
    run_manifest_health(read_run_manifest(path); manifest_path = path, kwargs...)

function run_manifest_health(
    manifest::AbstractDict;
    manifest_path = nothing,
    root = nothing,
    summarize_outputs::Bool = true,
)
    manifest_issues = run_manifest_issues(manifest)
    resolved_root = _manifest_health_root(manifest, manifest_path, root)

    input_rows =
        _manifest_section_health(manifest, "inputs", resolved_root; summarize = false)
    output_rows = _manifest_section_health(
        manifest,
        "outputs",
        resolved_root;
        summarize = summarize_outputs,
    )
    generated_rows =
        _manifest_section_health(manifest, "generated", resolved_root; summarize = false)
    all_rows = vcat(input_rows, output_rows, generated_rows)
    health_status =
        isempty(manifest_issues) && all(row["status"] == "ok" for row in all_rows) ? "ok" :
        "attention"

    health = OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => RUN_MANIFEST_HEALTH_SCHEMA_VERSION,
        "status" => health_status,
        "manifest_path" =>
            isnothing(manifest_path) ? nothing : abspath(string(manifest_path)),
        "root" => resolved_root,
        "manifest_issues" => manifest_issues,
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "records" => length(all_rows),
            "ok" => count(row -> row["status"] == "ok", all_rows),
            "modified" => count(row -> row["status"] == "modified", all_rows),
            "missing" => count(row -> row["status"] == "missing", all_rows),
            "invalid_record" =>
                count(row -> row["status"] == "invalid_record", all_rows),
        ),
        "inputs" => input_rows,
        "outputs" => output_rows,
        "generated" => generated_rows,
    )
    if haskey(manifest, "capability_gates")
        health["capability_gates"] = _manifest_value(manifest["capability_gates"])
    end
    return health
end

function _manifest_health_root(manifest::AbstractDict, manifest_path, root)
    if !isnothing(root)
        return _canonical_abs_path(string(root))
    elseif haskey(manifest, "project_root") && manifest["project_root"] isa AbstractString
        project_root = manifest["project_root"]
        if isabspath(project_root) || isnothing(manifest_path)
            return _canonical_abs_path(project_root)
        end
        return _canonical_abs_path(
            joinpath(dirname(abspath(string(manifest_path))), project_root),
        )
    elseif !isnothing(manifest_path)
        return _canonical_abs_path(dirname(abspath(string(manifest_path))))
    end

    return _canonical_abs_path(pwd())
end

function _manifest_section_health(
    manifest::AbstractDict,
    section::AbstractString,
    root::AbstractString;
    summarize::Bool,
)
    records = get(manifest, section, OrderedCollections.OrderedDict{String,Any}[])
    records isa AbstractVector || return OrderedCollections.OrderedDict{String,Any}[]

    rows = OrderedCollections.OrderedDict{String,Any}[]
    for (index, record) in enumerate(records)
        row = verify_file_provenance(record; root)
        if row["status"] != "ok"
            row["remediation"] = _run_artifact_remediation(row, section, index)
        end
        if summarize &&
           row["status"] == "ok" &&
           endswith(lowercase(row["resolved_path"]), ".h5")
            row["output_data_summary"] = _safe_output_data_summary(row["resolved_path"])
        end
        push!(rows, row)
    end

    return rows
end

function _safe_output_data_summary(path::AbstractString)
    try
        return [_validation_report_value(row) for row in output_data_summary(path)]
    catch err
        return OrderedCollections.OrderedDict{String,Any}(
            "status" => "summary_error",
            "message" => sprint(showerror, err),
        )
    end
end

function _run_artifact_remediation(row::AbstractDict, section::AbstractString, index::Integer)
    status = string(get(row, "status", "attention"))
    role = _record_get_string(row, "role")
    artifact = isnothing(role) || isempty(role) ? _run_artifact_label(section) : role
    code =
        status == "missing" ? "missing_run_artifact" :
        status == "modified" ? "modified_run_artifact" :
        status == "invalid_record" ? "invalid_run_artifact_record" :
        "run_artifact_attention"
    severity = status == "modified" ? "warning" : "error"

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => RUN_ARTIFACT_REMEDIATION_SCHEMA_VERSION,
        "code" => code,
        "severity" => severity,
        "section" => string(section),
        "index" => Int(index),
        "field" => "$(section)[$index].path",
        "role" => role,
        "path" => _record_get_string(row, "path"),
        "resolved_path" => _record_get_string(row, "resolved_path"),
        "status" => status,
        "issues" => String.(get(row, "issues", String[])),
        "physical_implication" =>
            _run_artifact_physical_implication(section, status, artifact),
        "suggested_fix" => _run_artifact_suggested_fix(section, status, artifact),
        "documentation" => RUN_ARTIFACT_REMEDIATION_DOC,
        "remediation_action" =>
            "$(status)_$(replace(string(section), r"[^A-Za-z0-9]+" => "_"))_artifact",
        "provenance" => OrderedCollections.OrderedDict{String,Any}(
            "expected_bytes" => get(row, "expected_bytes", nothing),
            "actual_bytes" => get(row, "actual_bytes", nothing),
            "expected_sha256" => get(row, "expected_sha256", nothing),
            "actual_sha256" => get(row, "actual_sha256", nothing),
        ),
    )
end

function _run_artifact_label(section::AbstractString)
    return section == "inputs" ? "input" :
           section == "outputs" ? "output" :
           section == "generated" ? "generated artifact" : "artifact"
end

function _run_artifact_physical_implication(
    section::AbstractString,
    status::AbstractString,
    artifact::AbstractString,
)
    if section == "outputs"
        return status == "missing" ?
               "Result browsers, validation reports, and postprocessing cannot load the recorded output artifact." :
               "The recorded output artifact no longer matches the run manifest, so comparisons and validation may not describe the original run."
    elseif section == "inputs"
        return status == "missing" ?
               "The run manifest cannot prove which input deck produced the recorded results." :
               "The recorded input artifact changed after the run, so the run is no longer reproducible from this manifest alone."
    elseif section == "generated"
        return status == "missing" ?
               "The generated script or derived artifact needed to audit the run is missing." :
               "The generated script or derived artifact changed after the run, so provenance is stale."
    end
    return "The $artifact provenance record is not healthy, so the run manifest needs attention before it is trusted."
end

function _run_artifact_suggested_fix(
    section::AbstractString,
    status::AbstractString,
    artifact::AbstractString,
)
    if status == "missing"
        if section == "outputs"
            return "Re-run the case or restore the missing output artifact at the recorded path, then refresh the run manifest provenance."
        elseif section == "inputs"
            return "Restore the missing input artifact or rebuild the run from available inputs, then refresh the run manifest provenance."
        elseif section == "generated"
            return "Regenerate or restore the missing generated artifact, then refresh the run manifest provenance."
        end
        return "Restore the missing $artifact and refresh the run manifest provenance."
    elseif status == "modified"
        return "Restore the original $artifact or explicitly refresh the run manifest only after confirming the change is intentional."
    elseif status == "invalid_record"
        return "Fix the malformed $artifact provenance record in the run manifest before using this run for validation or reporting."
    end
    return "Inspect the $artifact artifact and refresh the run manifest after the issue is resolved."
end

function _file_provenance_record_issues(record)
    issues = String[]
    if !(record isa AbstractDict)
        return String["record must be a dictionary"]
    end

    if !haskey(record, "path")
        push!(issues, "path is required")
    elseif !(record["path"] isa AbstractString)
        push!(issues, "path must be a string")
    end

    if !haskey(record, "bytes")
        push!(issues, "bytes is required")
    elseif !_is_nonnegative_integer(record["bytes"])
        push!(issues, "bytes must be a non-negative integer")
    end

    if !haskey(record, "sha256")
        push!(issues, "sha256 is required")
    elseif !_is_sha256_digest(record["sha256"])
        push!(issues, "sha256 must be a lowercase 64-character SHA-256 digest")
    end

    if haskey(record, "role") && !(record["role"] isa AbstractString)
        push!(issues, "role must be a string when present")
    end

    return issues
end

function _resolve_manifest_file_path(path::AbstractString, root::AbstractString)
    return isabspath(path) ? normpath(path) : normpath(joinpath(root, path))
end

function _canonical_abs_path(path::AbstractString)
    return joinpath(splitpath(normpath(abspath(path)))...)
end

function _record_get_string(record, key::AbstractString)
    record isa AbstractDict || return nothing
    value = get(record, key, nothing)
    return value isa AbstractString ? string(value) : nothing
end

function _record_get_value(record, key::AbstractString)
    record isa AbstractDict || return nothing
    return get(record, key, nothing)
end
