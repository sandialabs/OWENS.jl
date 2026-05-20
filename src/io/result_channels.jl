export ResultChannel,
    output_data_channels,
    output_data_channel,
    output_data_channel_names,
    annotate_output_data_channels!

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
