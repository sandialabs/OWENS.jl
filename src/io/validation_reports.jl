export build_output_data_validation_report,
    write_output_data_validation_report, read_output_data_validation_report

const OUTPUT_DATA_VALIDATION_REPORT_SCHEMA_VERSION = "owens-output-data-validation/v1"
const OUTPUT_DATA_VALIDATION_NOTE = "Elementwise comparison of selected registered outputData channels; no time alignment, interpolation, phase correction, or whole-revolution windowing is applied."

"""
    build_output_data_validation_report(reference_path, candidate_path; channels, atol, rtol, tolerances, case_id, metadata, created_at_utc)

Build a deterministic YAML-friendly report around `output_data_channel_metrics`.
The report is intended for regression dashboards and GUI validation panels; it
compares selected registered channels elementwise and records per-channel metric
rows, tolerances, file hashes, and a concise summary.
"""
function build_output_data_validation_report(
    reference_path::AbstractString,
    candidate_path::AbstractString;
    channels = output_data_channel_names(),
    atol::Real = 0.0,
    rtol::Real = 1e-6,
    tolerances = OrderedCollections.OrderedDict{String,Any}(),
    case_id::AbstractString = "output-data-comparison",
    metadata = OrderedCollections.OrderedDict{String,Any}(),
    created_at_utc = nothing,
)
    requested_channels = _output_summary_channels(channels)
    tolerance_map =
        _output_validation_tolerances(requested_channels, tolerances, atol, rtol)
    rows = OrderedCollections.OrderedDict{String,Any}[]
    for channel in requested_channels
        channel_tolerance = tolerance_map[channel]
        metrics = output_data_channel_metrics(
            reference_path,
            candidate_path;
            channels = channel,
            atol = channel_tolerance["atol"],
            rtol = channel_tolerance["rtol"],
        )
        push!(rows, _output_validation_metric_row(metrics[1], channel_tolerance))
    end

    passed = all(row["passed"] === true for row in rows)
    comparable_count = count(row -> row["comparable"] === true, rows)

    return OrderedCollections.OrderedDict{String,Any}(
        "schema_version" => OUTPUT_DATA_VALIDATION_REPORT_SCHEMA_VERSION,
        "case_id" => string(case_id),
        "created_at_utc" =>
            isnothing(created_at_utc) ? _utc_timestamp() : string(created_at_utc),
        "status" => passed ? "passed" : "failed",
        "comparison" => OrderedCollections.OrderedDict{String,Any}(
            "kind" => "outputData_channel_metrics",
            "note" => OUTPUT_DATA_VALIDATION_NOTE,
            "default_atol" => Float64(atol),
            "default_rtol" => Float64(rtol),
        ),
        "reference" => _absolute_file_provenance(reference_path, "reference_output_data"),
        "candidate" => _absolute_file_provenance(candidate_path, "candidate_output_data"),
        "summary" => OrderedCollections.OrderedDict{String,Any}(
            "channels_requested" => length(rows),
            "channels_comparable" => comparable_count,
            "channels_noncomparable" => length(rows) - comparable_count,
            "channels_passed" => count(row -> row["passed"] === true, rows),
            "channels_failed" => count(row -> row["passed"] === false, rows),
        ),
        "metadata" => _validation_report_value(metadata),
        "channels" => rows,
    )
end

"""
    write_output_data_validation_report(path, reference_path, candidate_path; kwargs...)

Build and write an output-data validation report. Returns the report that was
written.
"""
function write_output_data_validation_report(
    path::AbstractString,
    reference_path::AbstractString,
    candidate_path::AbstractString;
    kwargs...,
)
    report = build_output_data_validation_report(reference_path, candidate_path; kwargs...)
    parent = dirname(path)
    if !isempty(parent)
        mkpath(parent)
    end

    YAML.write_file(path, report)
    return report
end

"""
    read_output_data_validation_report(path)

Read a YAML output-data validation report using string-keyed ordered
dictionaries.
"""
function read_output_data_validation_report(path::AbstractString)
    isfile(path) || throw(ArgumentError("Cannot read missing validation report: $path"))
    return YAML.load_file(path; dicttype = OrderedCollections.OrderedDict{String,Any})
end

function _output_validation_tolerances(
    channels::AbstractVector{<:AbstractString},
    tolerances,
    default_atol::Real,
    default_rtol::Real,
)
    default_atol >= 0 || throw(
        ArgumentError("default validation atol must be non-negative, got $default_atol"),
    )
    default_rtol >= 0 || throw(
        ArgumentError("default validation rtol must be non-negative, got $default_rtol"),
    )

    tolerance_map = OrderedCollections.OrderedDict{String,Any}()
    for channel in channels
        channel_tolerance = _channel_validation_tolerance(tolerances, channel)
        atol = get(channel_tolerance, "atol", default_atol)
        rtol = get(channel_tolerance, "rtol", default_rtol)
        atol >= 0 || throw(
            ArgumentError("validation atol for $channel must be non-negative, got $atol"),
        )
        rtol >= 0 || throw(
            ArgumentError("validation rtol for $channel must be non-negative, got $rtol"),
        )
        tolerance_map[string(channel)] = OrderedCollections.OrderedDict{String,Any}(
            "atol" => Float64(atol),
            "rtol" => Float64(rtol),
        )
    end

    return tolerance_map
end

function _channel_validation_tolerance(tolerances, channel::AbstractString)
    if isnothing(tolerances)
        return OrderedCollections.OrderedDict{String,Any}()
    end

    if !(tolerances isa AbstractDict)
        throw(ArgumentError("tolerances must be a dictionary keyed by channel name"))
    end
    isempty(tolerances) && return OrderedCollections.OrderedDict{String,Any}()

    raw_tolerance =
        haskey(tolerances, channel) ? tolerances[channel] :
        get(tolerances, Symbol(channel), OrderedCollections.OrderedDict{String,Any}())
    if raw_tolerance isa NamedTuple
        return OrderedCollections.OrderedDict{String,Any}(
            string(key) => raw_tolerance[key] for key in keys(raw_tolerance)
        )
    elseif raw_tolerance isa AbstractDict
        return OrderedCollections.OrderedDict{String,Any}(
            string(key) => value for (key, value) in raw_tolerance
        )
    else
        throw(ArgumentError("tolerance for $channel must be a dictionary or named tuple"))
    end
end

function _output_validation_metric_row(metric, tolerance::AbstractDict)
    return OrderedCollections.OrderedDict{String,Any}(
        "name" => metric.name,
        "passed" => metric.passed,
        "comparable" => metric.comparable,
        "status" => metric.status,
        "units" => metric.units,
        "dimensions" => copy(metric.dimensions),
        "reference_shape" => _validation_report_value(metric.reference_shape),
        "candidate_shape" => _validation_report_value(metric.candidate_shape),
        "n" => metric.n,
        "bias" => _validation_report_value(metric.bias),
        "rmse" => _validation_report_value(metric.rmse),
        "max_abs_error" => _validation_report_value(metric.max_abs_error),
        "mean_abs_error" => _validation_report_value(metric.mean_abs_error),
        "reference_rms" => _validation_report_value(metric.reference_rms),
        "max_abs_reference" => _validation_report_value(metric.max_abs_reference),
        "relative_rmse" => _validation_report_value(metric.relative_rmse),
        "tolerance" => _validation_report_value(metric.tolerance),
        "atol" => tolerance["atol"],
        "rtol" => tolerance["rtol"],
    )
end

function _absolute_file_provenance(path::AbstractString, role::AbstractString)
    record = file_provenance(path; root = dirname(abspath(path)), role = role)
    record["path"] = abspath(path)
    return record
end

function _validation_report_value(value)
    if value === missing
        return nothing
    elseif value isa NamedTuple
        return OrderedCollections.OrderedDict{String,Any}(
            string(key) => _validation_report_value(getfield(value, key)) for
            key in keys(value)
        )
    elseif value isa AbstractDict
        normalized = OrderedCollections.OrderedDict{String,Any}()
        keys_iter =
            value isa OrderedCollections.OrderedDict ? collect(keys(value)) :
            sort(collect(keys(value)); by = string)
        for key in keys_iter
            normalized[string(key)] = _validation_report_value(value[key])
        end
        return normalized
    elseif value isa AbstractVector || value isa Tuple
        return [_validation_report_value(item) for item in value]
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
