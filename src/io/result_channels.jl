export ResultChannel,
    output_data_channels,
    output_data_channel,
    output_data_channel_names,
    annotate_output_data_channels!,
    output_data_summary,
    output_data_channel_metrics

struct ResultChannel
    name::String
    units::String
    dimensions::Vector{String}
    frame::String
    association::String
    source::String
    sign_convention::String
    description::String
end

function ResultChannel(
    name::AbstractString,
    units::AbstractString,
    dimensions::AbstractVector,
    frame::AbstractString,
    association::AbstractString,
    source::AbstractString,
    sign_convention::AbstractString,
    description::AbstractString,
)
    return ResultChannel(
        string(name),
        string(units),
        string.(collect(dimensions)),
        string(frame),
        string(association),
        string(source),
        string(sign_convention),
        string(description),
    )
end

const OUTPUT_DATA_CHANNELS = ResultChannel[
    ResultChannel(
        "t",
        "s",
        ["time"],
        "scalar",
        "time",
        "OWENS runtime",
        "positive elapsed simulation time",
        "Simulation time at each saved state.",
    ),
    ResultChannel(
        "aziHist",
        "rad",
        ["time"],
        "rotor azimuth",
        "rotor",
        "OWENS runtime",
        "positive in the OWENS rotor azimuth convention",
        "Rotor azimuth history.",
    ),
    ResultChannel(
        "OmegaHist",
        "Hz",
        ["time"],
        "rotor",
        "rotor",
        "OWENS runtime",
        "positive in the OWENS rotor-speed convention",
        "Rotor rotational speed history.",
    ),
    ResultChannel(
        "OmegaDotHist",
        "Hz/s",
        ["time"],
        "rotor",
        "rotor",
        "OWENS runtime",
        "positive in the OWENS rotor-acceleration convention",
        "Rotor rotational acceleration history.",
    ),
    ResultChannel(
        "gbHist",
        "rad",
        ["time"],
        "drivetrain",
        "gearbox",
        "drivetrain model",
        "positive in the gearbox azimuth convention",
        "Gearbox azimuth history when the drivetrain model is active.",
    ),
    ResultChannel(
        "gbDotHist",
        "Hz",
        ["time"],
        "drivetrain",
        "gearbox",
        "drivetrain model",
        "positive in the gearbox speed convention",
        "Gearbox rotational speed history when the drivetrain model is active.",
    ),
    ResultChannel(
        "gbDotDotHist",
        "Hz/s",
        ["time"],
        "drivetrain",
        "gearbox",
        "drivetrain model",
        "positive in the gearbox acceleration convention",
        "Gearbox rotational acceleration history when the drivetrain model is active.",
    ),
    ResultChannel(
        "FReactionHist",
        "N or N*m by DOF",
        ["time", "structural_dof"],
        "structural mesh",
        "node dof",
        "OWENSFEA reactions",
        "unknown",
        "Reaction force and moment history grouped by six structural degrees of freedom per node.",
    ),
    ResultChannel(
        "FTwrBsHist",
        "N or N*m by DOF",
        ["time", "tower_base_dof"],
        "global",
        "tower base",
        "floating-platform coupling",
        "unknown",
        "Tower-base force and moment history.",
    ),
    ResultChannel(
        "genTorque",
        "N*m",
        ["time"],
        "drivetrain",
        "generator",
        "generator model",
        "unknown",
        "Generator torque history.",
    ),
    ResultChannel(
        "genPower",
        "W",
        ["time"],
        "drivetrain",
        "generator",
        "generator model",
        "positive generator power convention",
        "Generator power history.",
    ),
    ResultChannel(
        "torqueDriveShaft",
        "N*m",
        ["time"],
        "drivetrain",
        "shaft",
        "drivetrain model",
        "unknown",
        "Drive-shaft torque history.",
    ),
    ResultChannel(
        "uHist",
        "m or rad by DOF",
        ["time", "structural_dof"],
        "structural mesh",
        "node dof",
        "OWENSFEA displacement state",
        "positive by structural DOF convention",
        "Structural displacement and rotation history grouped by six degrees of freedom per node.",
    ),
    ResultChannel(
        "uHist_prp",
        "m or rad by DOF",
        ["time", "platform_dof"],
        "global",
        "platform reference point",
        "floating-platform coupling",
        "positive by platform DOF convention",
        "Platform-reference-point displacement and rotation history.",
    ),
    ResultChannel(
        "epsilon_x_hist",
        "strain",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA strain recovery",
        "positive beam axial strain",
        "Axial beam strain history.",
    ),
    ResultChannel(
        "epsilon_y_hist",
        "strain",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA strain recovery",
        "positive local y shear/strain convention",
        "Local y strain history.",
    ),
    ResultChannel(
        "epsilon_z_hist",
        "strain",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA strain recovery",
        "positive local z shear/strain convention",
        "Local z strain history.",
    ),
    ResultChannel(
        "kappa_x_hist",
        "1/m",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA curvature recovery",
        "positive twist curvature convention",
        "Beam twist curvature history.",
    ),
    ResultChannel(
        "kappa_y_hist",
        "1/m",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA curvature recovery",
        "positive local y bending curvature convention",
        "Beam local-y bending curvature history.",
    ),
    ResultChannel(
        "kappa_z_hist",
        "1/m",
        ["quadrature_point", "element", "step"],
        "beam element",
        "element quadrature",
        "OWENSFEA curvature recovery",
        "positive local z bending curvature convention",
        "Beam local-z bending curvature history.",
    ),
    ResultChannel(
        "massOwens",
        "kg",
        String[],
        "scalar",
        "model",
        "OWENS mass integration",
        "positive mass",
        "Integrated OWENS structural mass.",
    ),
    ResultChannel(
        "stress_U",
        "Pa",
        ["time", "span_station", "chord_station", "stress_component"],
        "composite material",
        "blade upper laminate",
        "composite stress recovery",
        "unknown",
        "Blade upper-laminate stress history.",
    ),
    ResultChannel(
        "SF_ult_U",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "blade upper laminate",
        "composite failure postprocessing",
        "larger is safer",
        "Blade upper-laminate ultimate safety-factor history.",
    ),
    ResultChannel(
        "SF_buck_U",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "blade upper laminate",
        "composite buckling postprocessing",
        "larger is safer",
        "Blade upper-laminate buckling safety-factor history.",
    ),
    ResultChannel(
        "stress_L",
        "Pa",
        ["time", "span_station", "chord_station", "stress_component"],
        "composite material",
        "blade lower laminate",
        "composite stress recovery",
        "unknown",
        "Blade lower-laminate stress history.",
    ),
    ResultChannel(
        "SF_ult_L",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "blade lower laminate",
        "composite failure postprocessing",
        "larger is safer",
        "Blade lower-laminate ultimate safety-factor history.",
    ),
    ResultChannel(
        "SF_buck_L",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "blade lower laminate",
        "composite buckling postprocessing",
        "larger is safer",
        "Blade lower-laminate buckling safety-factor history.",
    ),
    ResultChannel(
        "stress_TU",
        "Pa",
        ["time", "span_station", "chord_station", "stress_component"],
        "composite material",
        "tower upper laminate",
        "composite stress recovery",
        "unknown",
        "Tower upper-laminate stress history.",
    ),
    ResultChannel(
        "SF_ult_TU",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "tower upper laminate",
        "composite failure postprocessing",
        "larger is safer",
        "Tower upper-laminate ultimate safety-factor history.",
    ),
    ResultChannel(
        "SF_buck_TU",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "tower upper laminate",
        "composite buckling postprocessing",
        "larger is safer",
        "Tower upper-laminate buckling safety-factor history.",
    ),
    ResultChannel(
        "stress_TL",
        "Pa",
        ["time", "span_station", "chord_station", "stress_component"],
        "composite material",
        "tower lower laminate",
        "composite stress recovery",
        "unknown",
        "Tower lower-laminate stress history.",
    ),
    ResultChannel(
        "SF_ult_TL",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "tower lower laminate",
        "composite failure postprocessing",
        "larger is safer",
        "Tower lower-laminate ultimate safety-factor history.",
    ),
    ResultChannel(
        "SF_buck_TL",
        "dimensionless",
        ["time", "span_station", "chord_station"],
        "composite material",
        "tower lower laminate",
        "composite buckling postprocessing",
        "larger is safer",
        "Tower lower-laminate buckling safety-factor history.",
    ),
    ResultChannel(
        "topstrainout_blade_U",
        "mixed",
        ["time", "span_station", "chord_station", "strain_component"],
        "composite material",
        "blade upper laminate",
        "composite strain postprocessing",
        "unknown",
        "Blade upper-laminate top strain and beam strain/curvature channels.",
    ),
    ResultChannel(
        "topstrainout_blade_L",
        "mixed",
        ["time", "span_station", "chord_station", "strain_component"],
        "composite material",
        "blade lower laminate",
        "composite strain postprocessing",
        "unknown",
        "Blade lower-laminate top strain and beam strain/curvature channels.",
    ),
    ResultChannel(
        "topstrainout_tower_U",
        "mixed",
        ["time", "span_station", "chord_station", "strain_component"],
        "composite material",
        "tower upper laminate",
        "composite strain postprocessing",
        "unknown",
        "Tower upper-laminate top strain and beam strain/curvature channels.",
    ),
    ResultChannel(
        "topstrainout_tower_L",
        "mixed",
        ["time", "span_station", "chord_station", "strain_component"],
        "composite material",
        "tower lower laminate",
        "composite strain postprocessing",
        "unknown",
        "Tower lower-laminate top strain and beam strain/curvature channels.",
    ),
    ResultChannel(
        "topDamage_blade_U",
        "unknown",
        ["span_station", "chord_station"],
        "composite material",
        "blade upper laminate",
        "fatigue postprocessing",
        "unknown",
        "Blade upper-laminate fatigue-damage output.",
    ),
    ResultChannel(
        "topDamage_blade_L",
        "unknown",
        ["span_station", "chord_station"],
        "composite material",
        "blade lower laminate",
        "fatigue postprocessing",
        "unknown",
        "Blade lower-laminate fatigue-damage output.",
    ),
    ResultChannel(
        "topDamage_tower_U",
        "unknown",
        ["span_station", "chord_station"],
        "composite material",
        "tower upper laminate",
        "fatigue postprocessing",
        "unknown",
        "Tower upper-laminate fatigue-damage output.",
    ),
    ResultChannel(
        "topDamage_tower_L",
        "unknown",
        ["span_station", "chord_station"],
        "composite material",
        "tower lower laminate",
        "fatigue postprocessing",
        "unknown",
        "Tower lower-laminate fatigue-damage output.",
    ),
]

const OUTPUT_DATA_CHANNEL_MAP = OrderedCollections.OrderedDict(
    channel.name => channel for channel in OUTPUT_DATA_CHANNELS
)

"""
    output_data_channels()

Return metadata for each root-level dataset written by `outputData`.
"""
output_data_channels() = copy(OUTPUT_DATA_CHANNELS)

"""
    output_data_channel(name)

Return the `ResultChannel` metadata for one `outputData` HDF5 dataset.
"""
function output_data_channel(name::AbstractString)
    key = string(name)
    haskey(OUTPUT_DATA_CHANNEL_MAP, key) ||
        throw(KeyError("No outputData result channel is registered for $key"))
    return OUTPUT_DATA_CHANNEL_MAP[key]
end

"""
    output_data_channel_names()

Return the dataset names written by `outputData` in writer order.
"""
output_data_channel_names() = [channel.name for channel in OUTPUT_DATA_CHANNELS]

"""
    annotate_output_data_channels!(file)

Attach channel metadata attributes to any registered `outputData` datasets that
are present in an open HDF5 file.
"""
function annotate_output_data_channels!(file)
    for channel in OUTPUT_DATA_CHANNELS
        if haskey(file, channel.name)
            dataset_attrs = HDF5.attrs(file[channel.name])
            dataset_attrs["owens_channel_name"] = channel.name
            dataset_attrs["units"] = channel.units
            dataset_attrs["dimensions"] = join(channel.dimensions, ",")
            dataset_attrs["frame"] = channel.frame
            dataset_attrs["association"] = channel.association
            dataset_attrs["source"] = channel.source
            dataset_attrs["sign_convention"] = channel.sign_convention
            dataset_attrs["description"] = channel.description
        end
    end

    return file
end

"""
    output_data_summary(path; channels=output_data_channel_names(), include_unregistered=false)

Inspect an OWENS `outputData` HDF5 file and return one summary row per
requested registered channel. Missing registered channels are returned with
`present=false`; root-level unregistered datasets are appended only when
`include_unregistered=true`. The helper reads dataset headers and attributes,
not full dataset contents.
"""
function output_data_summary(
    path::AbstractString;
    channels = output_data_channel_names(),
    include_unregistered::Bool = false,
)
    isfile(path) || throw(ArgumentError("Cannot summarize missing output file: $path"))

    requested_channels = _output_summary_channels(channels)
    for name in requested_channels
        haskey(OUTPUT_DATA_CHANNEL_MAP, name) ||
            throw(KeyError("No outputData result channel is registered for $name"))
    end

    HDF5.h5open(path, "r") do file
        summaries = [_registered_output_summary(file, name) for name in requested_channels]
        if include_unregistered
            registered_names = output_data_channel_names()
            extras = sort(setdiff(_root_dataset_names(file), registered_names))
            return vcat(
                summaries,
                [_unregistered_output_summary(file, name) for name in extras],
            )
        end

        return summaries
    end
end

_root_dataset_names(file) =
    String[name for name in keys(file) if file[name] isa HDF5.Dataset]

_output_summary_channels(channels::AbstractString) = [String(channels)]
_output_summary_channels(channels) = String.(collect(channels))

function _registered_output_summary(file, name::AbstractString)
    channel = output_data_channel(name)
    if !haskey(file, name) || !(file[name] isa HDF5.Dataset)
        return _output_summary_row(
            channel;
            present = false,
            shape = missing,
            ndims = missing,
            eltype = missing,
            has_channel_attrs = false,
            attr_mismatches = String[],
        )
    end

    dataset = file[name]
    dataset_attrs = HDF5.attrs(dataset)
    return _output_summary_row(
        channel;
        present = true,
        shape = collect(Int, size(dataset)),
        ndims = length(size(dataset)),
        eltype = string(eltype(dataset)),
        has_channel_attrs = _has_channel_attrs(dataset_attrs),
        attr_mismatches = _channel_attr_mismatches(dataset_attrs, channel),
    )
end

function _unregistered_output_summary(file, name::AbstractString)
    dataset = file[name]
    return (
        name = string(name),
        present = true,
        registered = false,
        shape = collect(Int, size(dataset)),
        ndims = length(size(dataset)),
        eltype = string(eltype(dataset)),
        units = missing,
        dimensions = missing,
        frame = missing,
        association = missing,
        source = missing,
        sign_convention = missing,
        description = missing,
        has_channel_attrs = false,
        attr_mismatches = String[],
    )
end

function _output_summary_row(
    channel::ResultChannel;
    present,
    shape,
    ndims,
    eltype,
    has_channel_attrs::Bool,
    attr_mismatches::Vector{String},
)
    return (
        name = channel.name,
        present,
        registered = true,
        shape,
        ndims,
        eltype,
        units = channel.units,
        dimensions = copy(channel.dimensions),
        frame = channel.frame,
        association = channel.association,
        source = channel.source,
        sign_convention = channel.sign_convention,
        description = channel.description,
        has_channel_attrs,
        attr_mismatches,
    )
end

const OUTPUT_DATA_CHANNEL_ATTRS = OrderedCollections.OrderedDict{String,Function}(
    "owens_channel_name" => channel -> channel.name,
    "units" => channel -> channel.units,
    "dimensions" => channel -> join(channel.dimensions, ","),
    "frame" => channel -> channel.frame,
    "association" => channel -> channel.association,
    "source" => channel -> channel.source,
    "sign_convention" => channel -> channel.sign_convention,
    "description" => channel -> channel.description,
)

function _has_channel_attrs(dataset_attrs)
    return all(haskey(dataset_attrs, key) for key in keys(OUTPUT_DATA_CHANNEL_ATTRS))
end

function _channel_attr_mismatches(dataset_attrs, channel::ResultChannel)
    mismatches = String[]
    for (key, expected_func) in OUTPUT_DATA_CHANNEL_ATTRS
        expected = expected_func(channel)
        if !haskey(dataset_attrs, key)
            push!(mismatches, "missing:$key")
        elseif dataset_attrs[key] != expected
            push!(
                mismatches,
                "$key: expected $(repr(expected)), got $(repr(dataset_attrs[key]))",
            )
        end
    end

    return mismatches
end

"""
    output_data_channel_metrics(reference_path, candidate_path; channels, atol=0.0, rtol=1e-6)

Compare selected numeric channels from two OWENS `outputData` HDF5 files. The
function reads only the requested datasets and returns per-channel bias, RMSE,
maximum absolute error, mean absolute error, reference RMS, and relative RMSE.
Missing, shape-mismatched, or non-numeric channels are returned as
`comparable=false` rows with an explanatory `status`.

A comparable channel passes when `max_abs_error <= atol + rtol * max_abs_reference`.
"""
function output_data_channel_metrics(
    reference_path::AbstractString,
    candidate_path::AbstractString;
    channels = output_data_channel_names(),
    atol::Real = 0.0,
    rtol::Real = 1e-6,
)
    isfile(reference_path) ||
        throw(ArgumentError("Cannot compare missing reference file: $reference_path"))
    isfile(candidate_path) ||
        throw(ArgumentError("Cannot compare missing candidate file: $candidate_path"))
    atol >= 0 || throw(ArgumentError("atol must be non-negative, got $atol"))
    rtol >= 0 || throw(ArgumentError("rtol must be non-negative, got $rtol"))

    requested_channels = _output_summary_channels(channels)
    for name in requested_channels
        haskey(OUTPUT_DATA_CHANNEL_MAP, name) ||
            throw(KeyError("No outputData result channel is registered for $name"))
    end

    HDF5.h5open(reference_path, "r") do reference_file
        HDF5.h5open(candidate_path, "r") do candidate_file
            return [
                _output_channel_metric(reference_file, candidate_file, name, atol, rtol) for
                name in requested_channels
            ]
        end
    end
end

function _output_channel_metric(
    reference_file,
    candidate_file,
    name::AbstractString,
    atol::Real,
    rtol::Real,
)
    channel = output_data_channel(name)
    reference_present =
        haskey(reference_file, name) && reference_file[name] isa HDF5.Dataset
    candidate_present =
        haskey(candidate_file, name) && candidate_file[name] isa HDF5.Dataset

    if !reference_present && !candidate_present
        return _noncomparable_metric_row(channel, "missing_reference_and_candidate")
    elseif !reference_present
        return _noncomparable_metric_row(
            channel,
            "missing_reference";
            candidate_shape = collect(Int, size(candidate_file[name])),
        )
    elseif !candidate_present
        return _noncomparable_metric_row(
            channel,
            "missing_candidate";
            reference_shape = collect(Int, size(reference_file[name])),
        )
    end

    reference_shape = collect(Int, size(reference_file[name]))
    candidate_shape = collect(Int, size(candidate_file[name]))
    if reference_shape != candidate_shape
        return _noncomparable_metric_row(
            channel,
            "shape_mismatch";
            reference_shape,
            candidate_shape,
        )
    end

    reference_data = HDF5.read(reference_file[name])
    candidate_data = HDF5.read(candidate_file[name])
    if !_is_real_output_data(reference_data) || !_is_real_output_data(candidate_data)
        return _noncomparable_metric_row(
            channel,
            "non_numeric";
            reference_shape,
            candidate_shape,
        )
    end

    reference_values = _real_output_values(reference_data)
    candidate_values = _real_output_values(candidate_data)
    if isempty(reference_values)
        return _noncomparable_metric_row(channel, "empty"; reference_shape, candidate_shape)
    end

    if !all(isfinite, reference_values) || !all(isfinite, candidate_values)
        return _noncomparable_metric_row(
            channel,
            "non_finite";
            reference_shape,
            candidate_shape,
        )
    end

    difference = candidate_values .- reference_values
    abs_difference = abs.(difference)
    n = length(difference)
    bias = sum(difference) / n
    rmse = sqrt(sum(abs2, difference) / n)
    max_abs_error = maximum(abs_difference)
    mean_abs_error = sum(abs_difference) / n
    reference_rms = sqrt(sum(abs2, reference_values) / n)
    max_abs_reference = maximum(abs.(reference_values))
    relative_rmse = _relative_error(rmse, reference_rms)
    tolerance = Float64(atol) + Float64(rtol) * max_abs_reference
    passed = max_abs_error <= tolerance

    return (
        name = channel.name,
        comparable = true,
        passed,
        status = passed ? "pass" : "fail",
        units = channel.units,
        dimensions = copy(channel.dimensions),
        reference_shape,
        candidate_shape,
        n,
        bias,
        rmse,
        max_abs_error,
        mean_abs_error,
        reference_rms,
        max_abs_reference,
        relative_rmse,
        tolerance,
    )
end

function _noncomparable_metric_row(
    channel::ResultChannel,
    status::AbstractString;
    reference_shape = missing,
    candidate_shape = missing,
)
    return (
        name = channel.name,
        comparable = false,
        passed = false,
        status = string(status),
        units = channel.units,
        dimensions = copy(channel.dimensions),
        reference_shape,
        candidate_shape,
        n = 0,
        bias = missing,
        rmse = missing,
        max_abs_error = missing,
        mean_abs_error = missing,
        reference_rms = missing,
        max_abs_reference = missing,
        relative_rmse = missing,
        tolerance = missing,
    )
end

_is_real_output_data(value) = value isa Real || value isa AbstractArray{<:Real}
_real_output_values(value::Real) = [Float64(value)]
_real_output_values(value::AbstractArray{<:Real}) = Float64.(vec(value))

function _relative_error(error::Real, scale::Real)
    if scale == 0
        return error == 0 ? 0.0 : Inf
    end

    return error / scale
end
