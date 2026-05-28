import .OWENS: SetupOutputs, ModelingOptions

"""
    TopData

Structure containing time history data and state variables for the top section of the wind turbine.

# Fields
* `delta_t::Float64`: Time step size (s)
* `numTS::Int`: Number of time steps
* `numDOFPerNode::Int`: Number of degrees of freedom per node
* `CN2H::Array{Float64}`: Coordinate transformation matrix
* `t::Array{Float64}`: Time array (s)
* `integrator::String`: Integration method
* `integrator_j::String`: Previous integration method
* `topDispOut::Array{Float64}`: Top displacement output
* `uHist::Array{Float64}`: Displacement history
* `udotHist::Array{Float64}`: Velocity history
* `uddotHist::Array{Float64}`: Acceleration history
* `epsilon_x_hist::Array{Float64}`: Axial strain history
* `epsilon_y_hist::Array{Float64}`: Transverse strain history (y)
* `epsilon_z_hist::Array{Float64}`: Transverse strain history (z)
* `kappa_x_hist::Array{Float64}`: Torsional curvature history
* `kappa_y_hist::Array{Float64}`: Bending curvature history (y)
* `kappa_z_hist::Array{Float64}`: Bending curvature history (z)
* `FReactionHist::Array{Float64}`: Reaction force history
* `FTwrBsHist::Array{Float64}`: Tower base force history
* `aziHist::Array{Float64}`: Azimuth angle history (rad)
* `OmegaHist::Array{Float64}`: Rotor speed history (rad/s)
* `OmegaDotHist::Array{Float64}`: Rotor acceleration history (rad/s²)
* `gbHist::Array{Float64}`: Gearbox position history
* `gbDotHist::Array{Float64}`: Gearbox velocity history
* `gbDotDotHist::Array{Float64}`: Gearbox acceleration history
* `genTorque::Array{Float64}`: Generator torque (N⋅m)
* `genPower::Array{Float64}`: Generator power (W)
* `torqueDriveShaft::Array{Float64}`: Drive shaft torque (N⋅m)
* `uHist_prp::Array{Float64}`: Prescribed displacement history
* `FPtfmHist::Array{Float64}`: Platform force history
* `FHydroHist::Array{Float64}`: Hydrodynamic force history
* `FMooringHist::Array{Float64}`: Mooring force history
* `rbData::Array{Float64}`: Rigid body data
* `rbDataHist::Array{Float64}`: Rigid body data history
* `u_s::Array{Float64}`: Current displacement state
* `udot_s::Array{Float64}`: Current velocity state
* `uddot_s::Array{Float64}`: Current acceleration state
* `u_sm1::Array{Float64}`: Previous displacement state
* `topDispData1::Array{Float64}`: Top displacement data set 1
* `topDispData2::Array{Float64}`: Top displacement data set 2
* `topElStrain::Array{Float64}`: Top element strain
* `gb_s::Float64`: Current gearbox position
* `gbDot_s::Float64`: Current gearbox velocity
* `gbDotDot_s::Float64`: Current gearbox acceleration
* `azi_s::Float64`: Current azimuth angle (rad)
* `Omega_s::Float64`: Current rotor speed (rad/s)
* `OmegaDot_s::Float64`: Current rotor acceleration (rad/s²)
* `genTorque_s::Float64`: Current generator torque (N⋅m)
* `torqueDriveShaft_s::Float64`: Current drive shaft torque (N⋅m)
* `topFexternal::Array{Float64}`: Current external forces
* `topFexternal_hist::Array{Float64}`: External force history
* `rotorSpeedForGenStart::Float64`: Rotor speed for generator start (rad/s)
* `top_rom::Array{Float64}`: Reduced order model data
* `topJointTransformTrans::Array{Float64}`: Joint transformation matrix
* `u_sRed::Array{Float64}`: Reduced displacement state
* `udot_sRed::Array{Float64}`: Reduced velocity state
* `uddot_sRed::Array{Float64}`: Reduced acceleration state
* `topBC::Array{Float64}`: Top boundary conditions
* `u_s2::Array{Float64}`: Secondary displacement state
* `udot_s2::Array{Float64}`: Secondary velocity state
* `uddot_s2::Array{Float64}`: Secondary acceleration state
* `top_invPhi::Array{Float64}`: Inverse modal matrix
* `eta_s::Array{Float64}`: Current modal coordinates
* `etadot_s::Array{Float64}`: Current modal velocities
* `etaddot_s::Array{Float64}`: Current modal accelerations
* `topsideMass::Float64`: Topside mass (kg)
* `topsideMOI::Array{Float64}`: Topside moment of inertia (kg⋅m²)
* `topsideCG::Array{Float64}`: Topside center of gravity (m)
* `u_j::Array{Float64}`: Previous displacement state
* `udot_j::Array{Float64}`: Previous velocity state
* `uddot_j::Array{Float64}`: Previous acceleration state
* `azi_j::Float64`: Previous azimuth angle (rad)
* `Omega_j::Float64`: Previous rotor speed (rad/s)
* `OmegaDot_j::Float64`: Previous rotor acceleration (rad/s²)
* `gb_j::Float64`: Previous gearbox position
* `gbDot_j::Float64`: Previous gearbox velocity
* `gbDotDot_j::Float64`: Previous gearbox acceleration
* `genTorque_j::Float64`: Previous generator torque (N⋅m)
* `torqueDriveShaft_j::Float64`: Previous drive shaft torque (N⋅m)
* `FReactionsm1::Array{Float64}`: Previous reaction forces
* `topFReaction_j::Array{Float64}`: Previous top reaction forces

# Notes
- This structure stores both time history data and current state variables
- Used for time integration and state tracking in unsteady analysis
- Contains data for structural, aerodynamic, and control system states
- Supports both full-order and reduced-order modeling
- Includes platform and mooring system data when applicable
"""
mutable struct TopData
    delta_t::Any
    numTS::Any
    numDOFPerNode::Any
    CN2H::Any
    t::Any
    integrator::Any
    integrator_j::Any
    topDispOut::Any
    uHist::Any
    udotHist::Any
    uddotHist::Any
    epsilon_x_hist::Any
    epsilon_y_hist::Any
    epsilon_z_hist::Any
    kappa_x_hist::Any
    kappa_y_hist::Any
    kappa_z_hist::Any
    FReactionHist::Any
    FTwrBsHist::Any
    aziHist::Any
    OmegaHist::Any
    OmegaDotHist::Any
    gbHist::Any
    gbDotHist::Any
    gbDotDotHist::Any
    genTorque::Any
    genPower::Any
    torqueDriveShaft::Any
    uHist_prp::Any
    FPtfmHist::Any
    FHydroHist::Any
    FMooringHist::Any
    rbData::Any
    rbDataHist::Any
    u_s::Any
    udot_s::Any
    uddot_s::Any
    u_sm1::Any
    topDispData1::Any
    topDispData2::Any
    topElStrain::Any
    gb_s::Any
    gbDot_s::Any
    gbDotDot_s::Any
    azi_s::Any
    Omega_s::Any
    OmegaDot_s::Any
    genTorque_s::Any
    torqueDriveShaft_s::Any
    topFexternal::Any
    topFexternal_hist::Any
    rotorSpeedForGenStart::Any
    top_rom::Any
    topJointTransformTrans::Any
    u_sRed::Any
    udot_sRed::Any
    uddot_sRed::Any
    topBC::Any
    u_s2::Any
    udot_s2::Any
    uddot_s2::Any
    top_invPhi::Any
    eta_s::Any
    etadot_s::Any
    etaddot_s::Any
    topsideMass::Any
    topsideMOI::Any
    topsideCG::Any
    u_j::Any
    udot_j::Any
    uddot_j::Any
    azi_j::Any
    Omega_j::Any
    OmegaDot_j::Any
    gb_j::Any
    gbDot_j::Any
    gbDotDot_j::Any
    genTorque_j::Any
    torqueDriveShaft_j::Any
    FReactionsm1::Any
    topFReaction_j::Any
end

const _UNSTEADY_LAND_RESTART_FORMAT = "OWENS.Unsteady_Land.restart"
const _UNSTEADY_LAND_RESTART_VERSION = 1

const _RESTART_CURRENT_FIELDS = (
    :CN2H,
    :integrator,
    :integrator_j,
    :topDispOut,
    :u_s,
    :udot_s,
    :uddot_s,
    :u_sm1,
    :topDispData1,
    :topDispData2,
    :topElStrain,
    :gb_s,
    :gbDot_s,
    :gbDotDot_s,
    :azi_s,
    :Omega_s,
    :OmegaDot_s,
    :genTorque_s,
    :torqueDriveShaft_s,
    :topFexternal,
    :rbData,
    :u_sRed,
    :udot_sRed,
    :uddot_sRed,
    :u_s2,
    :udot_s2,
    :uddot_s2,
    :eta_s,
    :etadot_s,
    :etaddot_s,
    :u_j,
    :udot_j,
    :uddot_j,
    :azi_j,
    :Omega_j,
    :OmegaDot_j,
    :gb_j,
    :gbDot_j,
    :gbDotDot_j,
    :genTorque_j,
    :torqueDriveShaft_j,
    :FReactionsm1,
    :topFReaction_j,
)

const _RESTART_STATE_HISTORY_FIELDS = (
    :uHist,
    :udotHist,
    :uddotHist,
    :FReactionHist,
    :FTwrBsHist,
    :aziHist,
    :OmegaHist,
    :OmegaDotHist,
    :gbHist,
    :gbDotHist,
    :gbDotDotHist,
    :genTorque,
    :genPower,
    :torqueDriveShaft,
    :uHist_prp,
    :FPtfmHist,
    :FHydroHist,
    :FMooringHist,
    :rbDataHist,
    :topFexternal_hist,
)

const _RESTART_STEP_HISTORY_FIELDS = (
    :epsilon_x_hist,
    :epsilon_y_hist,
    :epsilon_z_hist,
    :kappa_x_hist,
    :kappa_y_hist,
    :kappa_z_hist,
)

const _RESTART_INPUT_RUNTIME_FIELDS = (:generatorOn, :omegaControl)

const _OWENSAERO_TURBINE_RESTART_FIELDS =
    (:r, :z, :twist, :delta, :omega, :centerX, :centerY, :rhoA)

const _OWENSAERO_ENVIRONMENT_RESTART_FIELDS = (
    :V_x,
    :V_y,
    :V_z,
    :V_twist,
    :aw_warm,
    :steplast,
    :idx_RPI,
    :V_wake_old,
    :BV_DynamicFlagL,
    :BV_DynamicFlagD,
    :alpha_last,
    :accel_flap,
    :accel_edge,
    :gravity,
)

const _OWENSAERO_LB_STATE_RESTART_FIELDS = (
    :cl_ref_le_last,
    :cl_ref_last,
    :fstat_last,
    :cv_last,
    :dp,
    :df,
    :dcnv,
    :le_separation_state,
    :s_lev,
)

const _OWENSAERO_CACHE_RESTART_GLOBALS = (
    :dt,
    :last_step1,
    :last_azi,
    :last_stepL,
    :last_stepU,
    :timelast,
    :z3D,
    :z3Dnorm,
    :delta,
    :aziL_save,
    :aziU_save,
    :startingtwist,
    :CPL,
    :CPU,
    :RpL,
    :RpU,
    :TpL,
    :TpU,
    :ZpL,
    :ZpU,
    :XpL,
    :YpL,
    :XpU,
    :YpU,
    :M_addedmass_Np_L,
    :M_addedmass_Np_U,
    :M_addedmass_Tp_L,
    :M_addedmass_Tp_U,
    :F_addedmass_Np_L,
    :F_addedmass_Np_U,
    :F_addedmass_Tp_L,
    :F_addedmass_Tp_U,
    :F_buoy_L,
    :F_buoy_U,
    :alphaL,
    :alphaU,
    :clL,
    :clU,
    :cd_afL,
    :cd_afU,
    :cm_afL,
    :cm_afU,
    :M25L,
    :M25U,
    :VlocL,
    :VlocU,
    :ReL,
    :ReU,
    :thetavecL,
    :thetavecU,
    :deltaL,
    :deltaU,
    :Fx_baseL,
    :Fx_baseU,
    :Fy_baseL,
    :Fy_baseU,
    :Fz_baseL,
    :Fz_baseU,
    :Mx_baseL,
    :Mx_baseU,
    :My_baseL,
    :My_baseU,
    :Mz_baseL,
    :Mz_baseU,
    :power2L,
    :power2U,
    :powerL,
    :powerU,
)

function _dict_get(dict, key::String, default = nothing)
    if haskey(dict, key)
        return dict[key]
    end
    symbol_key = Symbol(key)
    if haskey(dict, symbol_key)
        return dict[symbol_key]
    end
    return default
end

function _field_dict(source, field_names)
    fields = Dict{String,Any}()
    for field_name in field_names
        if hasfield(typeof(source), field_name)
            fields[string(field_name)] = deepcopy(getfield(source, field_name))
        end
    end
    return fields
end

function _restart_history_prefix(
    value,
    filled_length::Integer,
    field_name::Symbol,
    axis::Integer,
)
    value === nothing && return nothing
    ndims(value) >= axis ||
        throw(DimensionMismatch("restart history $field_name does not have axis $axis"))
    0 <= filled_length <= size(value, axis) || throw(
        ArgumentError(
            "restart history $field_name has invalid filled length $filled_length",
        ),
    )
    indices = ntuple(dim -> dim == axis ? (1:filled_length) : Colon(), ndims(value))
    return copy(value[indices...])
end

function _history_dict(topdata, field_names, filled_length::Integer, axis::Integer)
    histories = Dict{String,Any}()
    for field_name in field_names
        if hasfield(typeof(topdata), field_name)
            histories[string(field_name)] = _restart_history_prefix(
                getfield(topdata, field_name),
                filled_length,
                field_name,
                axis,
            )
        end
    end
    return histories
end

function _copy_restart_value!(target, field_name::Symbol, saved_value)
    current_value = getfield(target, field_name)
    if current_value isa AbstractArray && saved_value isa AbstractArray
        size(current_value) == size(saved_value) || throw(
            DimensionMismatch(
                "restart field $field_name has size $(size(saved_value)); current model expects $(size(current_value))",
            ),
        )
    end
    setfield!(target, field_name, deepcopy(saved_value))
    return nothing
end

function _copy_restart_history_prefix!(
    topdata,
    field_name::Symbol,
    saved_value,
    filled_length::Integer,
    axis::Integer,
)
    saved_value === nothing && return nothing
    target_value = getfield(topdata, field_name)
    ndims(target_value) == ndims(saved_value) || throw(
        DimensionMismatch(
            "restart history $field_name has $(ndims(saved_value)) dimensions; current model expects $(ndims(target_value))",
        ),
    )
    size(saved_value, axis) == filled_length || throw(
        DimensionMismatch(
            "restart history $field_name has saved length $(size(saved_value, axis)); expected $filled_length",
        ),
    )
    size(target_value, axis) >= filled_length || throw(
        DimensionMismatch(
            "restart history $field_name needs $filled_length entries; current model has $(size(target_value, axis))",
        ),
    )

    for dim = 1:ndims(target_value)
        if dim != axis && size(target_value, dim) != size(saved_value, dim)
            throw(
                DimensionMismatch(
                    "restart history $field_name has size $(size(saved_value)); current model expects compatible size $(size(target_value))",
                ),
            )
        end
    end

    if filled_length > 0
        indices =
            ntuple(dim -> dim == axis ? (1:filled_length) : Colon(), ndims(target_value))
        target_value[indices...] .= saved_value
    end
    return nothing
end

function _copy_array_fields(source, field_names)
    fields = Dict{String,Any}()
    for field_name in field_names
        if hasfield(typeof(source), field_name)
            value = getfield(source, field_name)
            fields[string(field_name)] = deepcopy(value)
        end
    end
    return fields
end

function _copy_lb_state(lb_state)
    lb_state === nothing && return nothing
    return _copy_array_fields(lb_state, _OWENSAERO_LB_STATE_RESTART_FIELDS)
end

function _restore_array_field!(target, field_name::Symbol, saved_value)
    hasfield(typeof(target), field_name) || return nothing
    current_value = getfield(target, field_name)
    if current_value isa AbstractArray && saved_value isa AbstractArray
        size(current_value) == size(saved_value) || throw(
            DimensionMismatch(
                "OWENSAero restart field $field_name has size $(size(saved_value)); current model expects $(size(current_value))",
            ),
        )
        current_value .= saved_value
    elseif current_value != saved_value
        throw(
            ArgumentError(
                "OWENSAero restart field $field_name cannot be restored into the current model",
            ),
        )
    end
    return nothing
end

function _restore_lb_state!(lb_state, saved_state)
    (lb_state === nothing || saved_state === nothing) && return nothing
    for field_name in _OWENSAERO_LB_STATE_RESTART_FIELDS
        saved_value = _dict_get(saved_state, string(field_name), nothing)
        saved_value === nothing && continue
        _restore_array_field!(lb_state, field_name, saved_value)
    end
    return nothing
end

function _capture_owensaero_restart_state()
    if !(isdefined(OWENSAero, :turbslices) && isdefined(OWENSAero, :envslices))
        return nothing
    end

    turbine_states = [
        _copy_array_fields(turbine, _OWENSAERO_TURBINE_RESTART_FIELDS) for
        turbine in getfield(OWENSAero, :turbslices)
    ]
    environment_states = Dict{String,Any}[]
    for environment in getfield(OWENSAero, :envslices)
        environment_state =
            _copy_array_fields(environment, _OWENSAERO_ENVIRONMENT_RESTART_FIELDS)
        if hasfield(typeof(environment), :lb_state)
            environment_state["lb_state"] = _copy_lb_state(getfield(environment, :lb_state))
        end
        push!(environment_states, environment_state)
    end

    cache_state = Dict{String,Any}()
    for global_name in _OWENSAERO_CACHE_RESTART_GLOBALS
        if isdefined(OWENSAero, global_name)
            cache_state[string(global_name)] = deepcopy(getfield(OWENSAero, global_name))
        end
    end

    return Dict{String,Any}(
        "turbines" => turbine_states,
        "environments" => environment_states,
        "cache" => cache_state,
    )
end

function _set_module_global!(module_ref, global_name::Symbol, value)
    Core.eval(module_ref, Expr(:(=), global_name, value))
    return nothing
end

function _restore_owensaero_restart_state!(state)
    state === nothing && return nothing
    if !(isdefined(OWENSAero, :turbslices) && isdefined(OWENSAero, :envslices))
        throw(
            ArgumentError(
                "OWENSAero restart state was saved, but OWENSAero is not initialized in the current run",
            ),
        )
    end

    turbines = getfield(OWENSAero, :turbslices)
    environments = getfield(OWENSAero, :envslices)
    saved_turbines = _dict_get(state, "turbines", [])
    saved_environments = _dict_get(state, "environments", [])

    length(turbines) == length(saved_turbines) || throw(
        DimensionMismatch(
            "OWENSAero restart has $(length(saved_turbines)) turbine slices; current model has $(length(turbines))",
        ),
    )
    length(environments) == length(saved_environments) || throw(
        DimensionMismatch(
            "OWENSAero restart has $(length(saved_environments)) environment slices; current model has $(length(environments))",
        ),
    )

    for (turbine, saved_turbine) in zip(turbines, saved_turbines)
        for field_name in _OWENSAERO_TURBINE_RESTART_FIELDS
            saved_value = _dict_get(saved_turbine, string(field_name), nothing)
            saved_value === nothing && continue
            _restore_array_field!(turbine, field_name, saved_value)
        end
    end

    for (environment, saved_environment) in zip(environments, saved_environments)
        for field_name in _OWENSAERO_ENVIRONMENT_RESTART_FIELDS
            saved_value = _dict_get(saved_environment, string(field_name), nothing)
            saved_value === nothing && continue
            _restore_array_field!(environment, field_name, saved_value)
        end
        if hasfield(typeof(environment), :lb_state)
            _restore_lb_state!(
                getfield(environment, :lb_state),
                _dict_get(saved_environment, "lb_state", nothing),
            )
        end
    end

    cache_state = _dict_get(state, "cache", Dict{String,Any}())
    for (global_name_string, value) in cache_state
        _set_module_global!(OWENSAero, Symbol(global_name_string), deepcopy(value))
    end
    return nothing
end

function _capture_restart_backend_states(; capture_owensaero = false)
    backend_states = Dict{String,Any}()
    if capture_owensaero
        owensaero_state = _capture_owensaero_restart_state()
        if !isnothing(owensaero_state)
            backend_states["OWENSAero"] = owensaero_state
        end
    end
    return backend_states
end

function _restore_restart_backend_states!(backend_states)
    backend_states === nothing && return nothing
    owensaero_state = _dict_get(backend_states, "OWENSAero", nothing)
    _restore_owensaero_restart_state!(owensaero_state)
    return nothing
end

function _restart_input_runtime_state(inputs)
    inputs === nothing && return Dict{String,Any}()
    return _field_dict(inputs, _RESTART_INPUT_RUNTIME_FIELDS)
end

function _restore_restart_input_runtime_state!(inputs, input_state)
    (inputs === nothing || input_state === nothing) && return nothing
    for field_name in _RESTART_INPUT_RUNTIME_FIELDS
        saved_value = _dict_get(input_state, string(field_name), nothing)
        saved_value === nothing && continue
        hasfield(typeof(inputs), field_name) &&
            setfield!(inputs, field_name, deepcopy(saved_value))
    end
    return nothing
end

function _infer_legacy_restart_history_index(topdata)
    latest_index = 1
    for field_name in _RESTART_STATE_HISTORY_FIELDS
        hasfield(typeof(topdata), field_name) || continue
        history = getfield(topdata, field_name)
        history isa AbstractArray || continue
        ndims(history) >= 1 || continue
        for index = size(history, 1):-1:1
            row =
                ndims(history) == 1 ? history[index] :
                view(history, index, ntuple(_ -> Colon(), ndims(history)-1)...)
            row_has_data = ndims(history) == 1 ? !iszero(row) : any(!iszero, row)
            if row_has_data
                latest_index = max(latest_index, index)
                break
            end
        end
    end
    return latest_index
end

function _restart_state_from_topdata(
    topdata;
    completed_step = nothing,
    inputs = nothing,
    system = nothing,
    capture_backend_state = false,
)
    if isnothing(completed_step)
        history_index = _infer_legacy_restart_history_index(topdata)
        completed_step = history_index - 1
    else
        history_index = Int(completed_step) + 1
    end

    0 <= completed_step <= topdata.numTS - 1 || throw(
        ArgumentError(
            "restart completed step $completed_step is outside 0:$(topdata.numTS - 1)",
        ),
    )
    history_index == completed_step + 1 ||
        throw(ArgumentError("restart history index must be completed_step + 1"))

    return Dict{String,Any}(
        "format" => _UNSTEADY_LAND_RESTART_FORMAT,
        "version" => _UNSTEADY_LAND_RESTART_VERSION,
        "completed_step" => Int(completed_step),
        "history_index" => Int(history_index),
        "delta_t" => topdata.delta_t,
        "numTS_at_write" => topdata.numTS,
        "time_at_write" => topdata.t[history_index],
        "current_fields" => _field_dict(topdata, _RESTART_CURRENT_FIELDS),
        "state_histories" =>
            _history_dict(topdata, _RESTART_STATE_HISTORY_FIELDS, history_index, 1),
        "step_histories" =>
            _history_dict(topdata, _RESTART_STEP_HISTORY_FIELDS, completed_step, 3),
        "input_runtime_state" => _restart_input_runtime_state(inputs),
        "system" => deepcopy(system),
        "backend_states" =>
            _capture_restart_backend_states(capture_owensaero = capture_backend_state),
    )
end

function _write_restart_state(
    filename,
    topdata,
    completed_step::Integer;
    inputs = nothing,
    system = nothing,
    capture_backend_state = false,
)
    restart_state = _restart_state_from_topdata(
        topdata;
        completed_step = completed_step,
        inputs = inputs,
        system = system,
        capture_backend_state = capture_backend_state,
    )
    directory = dirname(abspath(filename))
    isdir(directory) || mkpath(directory)
    temporary_filename = tempname(directory) * ".jld2"
    try
        JLD2.jldsave(temporary_filename; restart_state)
        mv(temporary_filename, filename; force = true)
    catch
        isfile(temporary_filename) && rm(temporary_filename; force = true)
        rethrow()
    end
    return restart_state
end

function _load_restart_state(filename)
    isfile(filename) || throw(ArgumentError("restart file does not exist: $filename"))
    data = JLD2.load(filename)
    if haskey(data, "restart_state")
        restart_state = data["restart_state"]
    elseif haskey(data, "topdata")
        restart_state = _restart_state_from_topdata(data["topdata"])
    else
        throw(
            ArgumentError("restart file $filename does not contain an OWENS restart state"),
        )
    end
    _dict_get(restart_state, "format", _UNSTEADY_LAND_RESTART_FORMAT) ==
    _UNSTEADY_LAND_RESTART_FORMAT ||
        throw(ArgumentError("restart file $filename has an unsupported format"))
    version = _dict_get(restart_state, "version", 0)
    version <= _UNSTEADY_LAND_RESTART_VERSION || throw(
        ArgumentError(
            "restart file $filename version $version is newer than this OWENS version supports",
        ),
    )
    return restart_state
end

function _restore_restart!(topdata, restart_state; inputs = nothing, system = nothing)
    completed_step = Int(_dict_get(restart_state, "completed_step", -1))
    history_index = Int(_dict_get(restart_state, "history_index", completed_step + 1))
    completed_step >= 0 ||
        throw(ArgumentError("restart completed step must be nonnegative"))
    history_index == completed_step + 1 ||
        throw(ArgumentError("restart history index must be completed_step + 1"))
    history_index <= topdata.numTS || throw(
        ArgumentError(
            "restart history index $history_index exceeds current numTS $(topdata.numTS)",
        ),
    )
    completed_step <= topdata.numTS - 1 || throw(
        ArgumentError(
            "restart completed step $completed_step exceeds current final step $(topdata.numTS - 1)",
        ),
    )

    saved_delta_t = _dict_get(restart_state, "delta_t", topdata.delta_t)
    isapprox(
        saved_delta_t,
        topdata.delta_t;
        rtol = 0,
        atol = eps(Float64) * max(1, abs(topdata.delta_t)),
    ) || throw(
        ArgumentError(
            "restart delta_t $saved_delta_t does not match current delta_t $(topdata.delta_t)",
        ),
    )

    current_fields = _dict_get(restart_state, "current_fields", Dict{String,Any}())
    for field_name in _RESTART_CURRENT_FIELDS
        saved_value = _dict_get(current_fields, string(field_name), nothing)
        saved_value === nothing && continue
        hasfield(typeof(topdata), field_name) &&
            _copy_restart_value!(topdata, field_name, saved_value)
    end

    state_histories = _dict_get(restart_state, "state_histories", Dict{String,Any}())
    for field_name in _RESTART_STATE_HISTORY_FIELDS
        saved_value = _dict_get(state_histories, string(field_name), nothing)
        saved_value === nothing && continue
        hasfield(typeof(topdata), field_name) && _copy_restart_history_prefix!(
            topdata,
            field_name,
            saved_value,
            history_index,
            1,
        )
    end

    step_histories = _dict_get(restart_state, "step_histories", Dict{String,Any}())
    for field_name in _RESTART_STEP_HISTORY_FIELDS
        saved_value = _dict_get(step_histories, string(field_name), nothing)
        saved_value === nothing && continue
        hasfield(typeof(topdata), field_name) && _copy_restart_history_prefix!(
            topdata,
            field_name,
            saved_value,
            completed_step,
            3,
        )
    end

    _restore_restart_input_runtime_state!(
        inputs,
        _dict_get(restart_state, "input_runtime_state", nothing),
    )
    _restore_restart_backend_states!(_dict_get(restart_state, "backend_states", nothing))

    restored_system = _dict_get(restart_state, "system", system)
    return (
        completed_step = completed_step,
        history_index = history_index,
        system = restored_system,
    )
end
"""

Unsteady(model,topModel,mesh,el,aero;getLinearizedMatrices=false)

Executable function for transient analysis. Provides the interface of various
    external module with transient structural dynamics analysis capability.

    # Input
    * `inputs::Model`: see ?Model
    * `topModel::FEAModel`: see ?OWENSFEA.FEAModel
    * `mesh::Mesh`: see ?OWENSFEA.Mesh
    * `el::El`: see ?OWENSFEA.El
    * `bin::Bin`: see ?Bin
    * `aero::function`: Fexternal, Fdof = aero(t) where Fexternal is the force on each affected mesh dof and Fdof is the corresponding DOFs affected
    * `getLinearizedMatrices::Bool`: Flag to save the linearized matrices
    * `elStorage::ElStorage.ElStorage`: Optional object containing stored element matrices
    * `u_s::Array{<:float}`: Optional warm start of top deflections, of length Nnodes x Ndof


    # Output
    * `t`: time array
    * `aziHist`: azimuthal history array
    * `OmegaHist`: rotational speed array history
    * `OmegaDotHist`: rotational acceleration array history
    * `gbHist`: gearbox position history array
    * `gbDotHist`: gearbox velocity history array
    * `gbDotDotHist`: gearbox acceleration history array
    * `FReactionHist`: Nodal reaction 6dof forces history
    * `rigidDof`:
    * `genTorque`: generator torque history
    * `genPower`: generator power history
    * `torqueDriveShaft`: driveshaft torque history
    * `uHist`: mesh displacement history for each dof
    * `epsilon_x_hist`: strain history for eps_xx_0 for each dof
    * `epsilon_y_hist`: strain history for eps_xx_z for each dof
    * `epsilon_z_hist`: strain history for eps_xx_y for each dof
    * `kappa_x_hist`: strain history for gam_xz_0 for each dof
    * `kappa_y_hist`: strain history for gam_xz_y for each dof
    * `kappa_z_hist`: strain history for gam_xy_0 for each dof
    """
# New interface that takes SetupOutputs and ModelingOptions directly
function Unsteady_Land(
    setup_outputs::SetupOutputs,
    modeling_options::ModelingOptions;
    returnold::Bool = true,
    getLinearizedMatrices::Bool = false,
    topElStorage = nothing,
    bottomElStorage = nothing,
    u_s = nothing,
    meshcontrolfunction = nothing,
    userDefinedGenerator = nothing,
    turbsimfile = nothing,
    dataDumpFilename = nothing,
    datadumpfrequency::Int = 1000,
    restart::Bool = false,
)

    # Determine if using AD15
    AD15On = modeling_options.OWENS_Options.AeroModel == "AD"

    # Create Inputs struct
    inputs = OWENS.Inputs(;
        verbosity = modeling_options.OWENS_Options.verbosity,
        analysisType = modeling_options.OWENS_Options.structuralModel,
        tocp = modeling_options.OWENSAero_Options.tocp,
        Omegaocp = [
            modeling_options.OWENSAero_Options.RPM,
            modeling_options.OWENSAero_Options.RPM,
        ] ./ 60,
        tocp_Vinf = modeling_options.OWENSAero_Options.tocp_Vinf,
        Vinfocp = [
            modeling_options.OWENSAero_Options.Vinf,
            modeling_options.OWENSAero_Options.Vinf,
        ],
        numTS = modeling_options.OWENS_Options.numTS,
        delta_t = modeling_options.OWENS_Options.delta_t,
        AD15On = AD15On,
        aeroLoadsOn = modeling_options.OWENS_Options.aeroLoadsOn,
    )

    # Create FEAModel
    fea_options = modeling_options.OWENSFEA_Options
    feamodel = OWENS.FEAModel(;
        analysisType = modeling_options.OWENS_Options.structuralModel,
        dataOutputFilename = modeling_options.OWENS_Options.dataOutputFilename,
        joint = setup_outputs.myjoint,
        platformTurbineConnectionNodeNumber = fea_options.platformTurbineConnectionNodeNumber,
        pBC = fea_options.pBC,
        nlOn = fea_options.nlOn,
        numNodes = setup_outputs.mymesh.numNodes,
        RayleighAlpha = fea_options.RayleighAlpha,
        RayleighBeta = fea_options.RayleighBeta,
        iterationType = fea_options.iterationType,
    )

    # Call the original Unsteady_Land function with extracted parameters
    unsteady_outputs = Unsteady_Land(
        inputs;
        topModel = feamodel,
        topMesh = setup_outputs.mymesh,
        topEl = setup_outputs.myel,
        aero = setup_outputs.aeroForces,
        deformAero = setup_outputs.deformAero,
        system = setup_outputs.system,
        assembly = setup_outputs.assembly,
        returnold = returnold,
        getLinearizedMatrices = getLinearizedMatrices,
        topElStorage = topElStorage,
        bottomElStorage = bottomElStorage,
        u_s = u_s,
        meshcontrolfunction = meshcontrolfunction,
        userDefinedGenerator = userDefinedGenerator,
        turbsimfile = turbsimfile,
        dataDumpFilename = dataDumpFilename,
        datadumpfrequency = datadumpfrequency,
        restart = restart,
    )

    # Populate Components with Strain Data
    for icomp = 1:length(setup_outputs.components)
        component = setup_outputs.components[icomp]

        startE = component.elNumbers[1]
        stopE = component.elNumbers[end]

        component.e_x = unsteady_outputs.epsilon_x_hist[1, startE:stopE, :]
        component.e_y = unsteady_outputs.epsilon_y_hist[1, startE:stopE, :]
        component.e_z = unsteady_outputs.epsilon_z_hist[1, startE:stopE, :]
        component.k_x = unsteady_outputs.kappa_x_hist[1, startE:stopE, :]
        component.k_y = unsteady_outputs.kappa_y_hist[1, startE:stopE, :]
        component.k_z = unsteady_outputs.kappa_z_hist[1, startE:stopE, :]
        OWENS.populate_gxbeam_sectional_recovery!(component)
    end

    return unsteady_outputs
end

# Original interface for backward compatibility
function Unsteady_Land(
    inputs;
    topModel = nothing,
    topMesh = nothing,
    topEl = nothing,
    aeroVals = nothing,
    aeroDOFs = nothing,
    aero = nothing,
    deformAero = nothing,
    bottomModel = nothing,
    bottomMesh = nothing,
    bottomEl = nothing,
    bin = nothing,
    getLinearizedMatrices = false,
    system = nothing,
    assembly = nothing,
    returnold = true, #TODO: should we initialize them in here? Unify the interface for ease?
    topElStorage = nothing,
    bottomElStorage = nothing,
    u_s = nothing,
    meshcontrolfunction = nothing,
    userDefinedGenerator = nothing,
    turbsimfile = nothing,
    dataDumpFilename = nothing,
    datadumpfrequency = 1000,
    restart = false,
)

    #..........................................................................
    #                             INITIALIZATION
    #..........................................................................

    if (!inputs.topsideOn) && (!inputs.platformActive)
        error("No structure is being simulated!")
    end
    if restart && isnothing(dataDumpFilename)
        throw(
            ArgumentError(
                "restart=true requires dataDumpFilename to point to a restart file",
            ),
        )
    end
    if !isnothing(dataDumpFilename)
        datadumpfrequency = Int(datadumpfrequency)
        datadumpfrequency > 0 || throw(
            ArgumentError(
                "datadumpfrequency must be positive when dataDumpFilename is set",
            ),
        )
    end

    ## General
    delta_t = inputs.delta_t
    numTS = Int(inputs.numTS)
    numDOFPerNode = 6
    CN2H = LinearAlgebra.I(3) # hub and inertial frames initialize as copies
    t = range(0, length = numTS, step = delta_t)
    integrator = 0.0 #for generator control algorithm
    integrator_j = 0.0
    topDispOut = [] #TODO: better way to control scope

    # g = [0.0, 0.0, -9.80665]
    uHist,
    epsilon_x_hist,
    epsilon_y_hist,
    epsilon_z_hist,
    kappa_x_hist,
    kappa_y_hist,
    kappa_z_hist,
    FReactionHist,
    FTwrBsHist,
    aziHist,
    OmegaHist,
    OmegaDotHist,
    gbHist,
    gbDotHist,
    gbDotDotHist,
    genTorque,
    genPower,
    torqueDriveShaft,
    uHist_prp,
    FPtfmHist,
    FHydroHist,
    FMooringHist,
    rbData,
    rbDataHist,
    udotHist,
    uddotHist = allocate_general(inputs, topModel, topMesh, numDOFPerNode, numTS, assembly)

    # Allocate memory for topside
    u_s,
    udot_s,
    uddot_s,
    u_sm1,
    topDispData1,
    topDispData2,
    topElStrain,
    gb_s,
    gbDot_s,
    gbDotDot_s,
    azi_s,
    Omega_s,
    OmegaDot_s,
    genTorque_s,
    torqueDriveShaft_s,
    topFexternal,
    topFexternal_hist =
        allocate_topside(inputs, topMesh, topEl, topModel, numDOFPerNode, u_s, assembly)

    ## Rotor mode initialization
    rotorSpeedForGenStart = initialize_generator!(inputs)

    ## Structural dynamics initialization
    if isnothing(topElStorage)
        topElStorage = OWENSFEA.initialElementCalculations(topModel, topEl, topMesh) #perform initial element calculations for conventional structural dynamics analysis
    end

    top_rom,
    topJointTransformTrans,
    u_sRed,
    udot_sRed,
    uddot_sRed,
    topBC,
    u_s2,
    udot_s2,
    uddot_s2,
    top_invPhi,
    eta_s,
    etadot_s,
    etaddot_s = initialize_ROM(
        deepcopy(topElStorage),
        deepcopy(topModel),
        deepcopy(topMesh),
        deepcopy(topEl),
        deepcopy(u_s),
        deepcopy(udot_s),
        deepcopy(uddot_s),
    )

    topDispData1.eta_s = eta_s
    topDispData1.etadot_s = etadot_s
    topDispData1.etaddot_s = etaddot_s
    topDispData2.eta_s = eta_s
    topDispData2.etadot_s = etadot_s
    topDispData2.etaddot_s = etaddot_s

    topsideMass, topsideMOI, topsideCG = OWENSFEA.calculateStructureMassProps(topElStorage)

    # TODO: clean this up
    u_j = 0.0
    udot_j = 0.0
    uddot_j = 0.0
    azi_j = 0.0
    Omega_j = 0.0
    OmegaDot_j = 0.0
    gb_j = 0.0
    gbDot_j = 0.0
    gbDotDot_j = 0.0
    genTorque_j = 0.0
    torqueDriveShaft_j = 0.0
    FReactionsm1 = 0.0
    topFReaction_j = 0.0
    # Package up into data struct
    topdata = TopData(
        delta_t,
        numTS,
        numDOFPerNode,
        CN2H,
        t,
        integrator,
        integrator_j,
        topDispOut,
        uHist,
        udotHist,
        uddotHist,
        epsilon_x_hist,
        epsilon_y_hist,
        epsilon_z_hist,
        kappa_x_hist,
        kappa_y_hist,
        kappa_z_hist,
        FReactionHist,
        FTwrBsHist,
        aziHist,
        OmegaHist,
        OmegaDotHist,
        gbHist,
        gbDotHist,
        gbDotDotHist,
        genTorque,
        genPower,
        torqueDriveShaft,
        uHist_prp,
        FPtfmHist,
        FHydroHist,
        FMooringHist,
        rbData,
        rbDataHist,
        u_s,
        udot_s,
        uddot_s,
        u_sm1,
        topDispData1,
        topDispData2,
        topElStrain,
        gb_s,
        gbDot_s,
        gbDotDot_s,
        azi_s,
        Omega_s,
        OmegaDot_s,
        genTorque_s,
        torqueDriveShaft_s,
        topFexternal,
        topFexternal_hist,
        rotorSpeedForGenStart,
        top_rom,
        topJointTransformTrans,
        u_sRed,
        udot_sRed,
        uddot_sRed,
        topBC,
        u_s2,
        udot_s2,
        uddot_s2,
        top_invPhi,
        eta_s,
        etadot_s,
        etaddot_s,
        topsideMass,
        topsideMOI,
        topsideCG,
        u_j,
        udot_j,
        uddot_j,
        azi_j,
        Omega_j,
        OmegaDot_j,
        gb_j,
        gbDot_j,
        gbDotDot_j,
        genTorque_j,
        torqueDriveShaft_j,
        FReactionsm1,
        topFReaction_j,
    )

    topModel.jointTransform, topModel.reducedDOFList =
        OWENSFEA.createJointTransform(topModel.joint, topMesh.numNodes, 6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints

    topdata.uHist[1, :] = topdata.u_s          #store initial condition
    topdata.aziHist[1] = topdata.azi_s
    topdata.OmegaHist[1] = topdata.Omega_s
    topdata.OmegaDotHist[1] = topdata.OmegaDot_s
    topdata.FReactionsm1 = zeros(length(topdata.u_s))
    topdata.FReactionHist[1, :] = topdata.FReactionsm1
    topdata.topFReaction_j = topdata.FReactionsm1
    # topWeight = [0.0, 0.0, topsideMass*-9.80665, 0.0, 0.0, 0.0] #TODO: propogate gravity, or remove since this isn't used
    topdata.gbHist[1] = gb_s
    topdata.gbDotHist[1] = gbDot_s
    topdata.gbDotDotHist[1] = gbDotDot_s
    topdata.rbDataHist[1, :] = zeros(9)
    topdata.genTorque[1] = genTorque_s
    topdata.torqueDriveShaft[1] = torqueDriveShaft_s


    #..........................................................................
    #                                MAIN LOOP
    #..........................................................................

    ### Iterate for a solution at t+dt
    i=0
    timeconverged = false
    pbar = ProgressBars.ProgressBar(total = numTS-1)

    if restart
        println(
            "\n Restarting from intermediate results from the following file: $dataDumpFilename \n",
        )
        restart_state = _load_restart_state(dataDumpFilename)
        restart_restore =
            _restore_restart!(topdata, restart_state; inputs = inputs, system = system)
        i = restart_restore.completed_step
        system = restart_restore.system
        if inputs.AD15On && inputs.aeroLoadsOn > 0
            @warn "AeroDyn/OpenFAST native module state is not included in OWENS restart files; AD15 restarts may need backend-specific support for bitwise-continuous aerodynamic transients."
        end
        for ii = 1:i
            ProgressBars.update(pbar)
        end
    end

    while (i<numTS-1) && timeconverged == false # we compute for the next time step, so the last step of our desired time series is computed in the second to last numTS value
        i += 1

        ProgressBars.update(pbar)

        ## Check for specified rotor speed at t+dt #TODO: fix this so that it can be probably accounted for in RK4
        if (inputs.turbineStartup == 0)
            inputs.omegaControl = true #TODO: are we setting this back?
            if (inputs.usingRotorSpeedFunction) #use user specified rotor speed profile function
                _, omegaCurrent, _ = getRotorPosSpeedAccelAtTime(t[i], t[i+1], 0.0, delta_t)
                topdata.Omega_s = omegaCurrent
            else #use discreteized rotor speed profile function
                omegaCurrent, OmegaDotCurrent, terminateSimulation =
                    omegaSpecCheck(t[i+1], inputs.tocp, inputs.Omegaocp, delta_t)
                if (terminateSimulation)
                    break
                end
                topdata.Omega_s = omegaCurrent
                topdata.OmegaDot_s = OmegaDotCurrent
            end
        else
            omegaCurrent = 0.0
        end

        ## Initialize "j" Gauss-Seidel iteration variables
        topdata.u_j = topdata.u_s
        topdata.udot_j = topdata.udot_s
        topdata.uddot_j = topdata.uddot_s

        topdata.azi_j = topdata.azi_s
        topdata.Omega_j = topdata.Omega_s
        topdata.OmegaDot_j = topdata.OmegaDot_s
        topdata.gb_j = topdata.gb_s
        topdata.gbDot_j = topdata.gbDot_s
        topdata.gbDotDot_j = topdata.gbDotDot_s
        topdata.genTorque_j = topdata.genTorque_s
        topdata.torqueDriveShaft_j = topdata.torqueDriveShaft_s

        #TODO: put these in the model
        TOL = inputs.iteration_parameters.TOL#1e-4  #gauss-seidel iteration tolerance for various modules
        MAXITER = inputs.iteration_parameters.MAXITER#2 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        aziNorm = 1e5
        platNorm = 0.0 #TODO: remove this?
        gbNorm = 0.0 #initialize norms for various module states
        topdata.integrator = topdata.integrator_j #don't compound integrator within the convergence loop.

        if inputs.analysisType=="GX"
            # systemout = deepcopy(system)
            strainGX = zeros(3, length(assembly.elements))
            curvGX = zeros(3, length(assembly.elements))
        end

        newVinf = 0.0 #TODO

        ## Gauss-Seidel predictor-corrector loop
        while ((uNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))
            # println("Iteration $numIterations")

            #------------------
            # GENERATOR MODULE
            #------------------
            topdata.genTorque_j = 0

            if inputs.generatorOn
                if inputs.useGeneratorFunction
                    specifiedOmega, _, _ = omegaSpecCheck(
                        t[i]+topdata.delta_t,
                        inputs.tocp,
                        inputs.Omegaocp,
                        topdata.delta_t,
                    )
                    newVinf = safeakima(inputs.tocp_Vinf, inputs.Vinfocp, t[i]) #TODO: ifw sampling of same file as aerodyn
                    if isnothing(userDefinedGenerator)
                        genTorqueHSS0, topdata.integrator_j, controlnamecurrent =
                            internaluserDefinedGenerator(
                                newVinf,
                                t[i],
                                topdata.azi_j,
                                topdata.Omega_j,
                                topdata.OmegaHist[i],
                                topdata.OmegaDot_j,
                                topdata.OmegaDotHist[i],
                                topdata.delta_t,
                                topdata.integrator,
                                specifiedOmega,
                            ) #;operPhase
                    else
                        if !isnothing(turbsimfile) #&& inputs.AD15On
                            velocity = OWENSOpenFASTWrappers.ifwcalcoutput(
                                [0.0, 0.0, maximum(topMesh.z)],
                                t[i],
                            )
                            newVinf = velocity[1]
                        end
                        genTorqueHSS0 =
                            userDefinedGenerator(t[i], topdata.Omega_j*60, newVinf) #;operPhase
                    end
                else
                    genTorqueHSS0 = simpleGenerator(inputs, topdata.Omega_j)
                end

                #should eventually account for Omega = gbDot*gearRatio here...
                topdata.genTorque_j =
                    genTorqueHSS0*inputs.gearRatio*inputs.gearBoxEfficiency #calculate generator torque on LSS side

            #         genTorqueAppliedToTurbineRotor0 = -genTorque0
            #         genTorqueAppliedToPlatform0 = genTorqueHSS0
            else
                if !isnothing(turbsimfile) && inputs.verbosity > 0#&& inputs.AD15On
                    try
                        velocity = OWENSOpenFASTWrappers.ifwcalcoutput(
                            [0.0, 0.0, maximum(topMesh.z)],
                            t[i],
                        )
                        newVinf = velocity[1]
                    catch
                        newVinf = safeakima(inputs.tocp_Vinf, inputs.Vinfocp, t[i]) #TODO: ifw sampling of same file as aerodyn
                    end
                else
                    newVinf = safeakima(inputs.tocp_Vinf, inputs.Vinfocp, t[i]) #TODO: ifw sampling of same file as aerodyn
                end
            end


            #-------------------
            # DRIVETRAIN MODULE
            #-------------------
            topdata.torqueDriveShaft_j = topdata.genTorque_j
            gb_jLast = topdata.gb_j
            if (!inputs.omegaControl)
                if (inputs.driveTrainOn)
                    topdata.torqueDriveShaft_j = calculateDriveShaftReactionTorque(
                        inputs.driveShaftProps,
                        topdata.azi_j,
                        topdata.gb_j,
                        topdata.Omega_j*2*pi,
                        topdata.gbDot_j*2*pi,
                    )

                    topdata.gb_j, topdata.gbDot_j, topdata.gbDotDot_j = updateRotorRotation(
                        inputs.JgearBox,
                        0,
                        0,
                        -topdata.genTorque_j,
                        topdata.torqueDriveShaft_j,
                        topdata.gb_s,
                        topdata.gbDot_s,
                        topdata.gbDotDot_s,
                        topdata.delta_t,
                    )
                else
                    topdata.gb_j = topdata.azi_j
                    topdata.gbDot_j = topdata.Omega_j
                    topdata.gbDotDot_j = topdata.OmegaDot_j
                end
            else
                topdata.gb_j = topdata.azi_j
                topdata.gbDot_j = omegaCurrent*2*pi
                topdata.gbDotDot_j = 0
            end

            # Update rotor speed
            azi_jLast = topdata.azi_j
            if inputs.omegaControl
                if (inputs.usingRotorSpeedFunction)
                    topdata.azi_j, topdata.Omega_j, topdata.OmegaDot_j =
                        getRotorPosSpeedAccelAtTime(t[i], t[i+1], azi_s, delta_t)
                else
                    topdata.Omega_j = topdata.Omega_s
                    topdata.OmegaDot_j = topdata.OmegaDot_s
                    topdata.azi_j = topdata.azi_s + topdata.Omega_j*topdata.delta_t*2*pi
                end
            elseif !inputs.omegaControl
                Crotor = 0
                Krotor = 0
                topdata.azi_j, topdata.Omega_j, topdata.OmegaDot_j = updateRotorRotation(
                    topdata.topsideMOI[3, 3],
                    Crotor,
                    Krotor,
                    -topdata.topFReaction_j[6],
                    -topdata.torqueDriveShaft_j,
                    topdata.azi_s,
                    topdata.Omega_s,
                    topdata.OmegaDot_s,
                    topdata.delta_t,
                )
            else
                error("omega control option not correctly specified")
            end

            #---------------------
            # AERODYNAMICS MODULE
            #---------------------
            # Calculate new aerodynamic loading

            # Update reference frame transformation and convert aerodynamic loads to hub reference frame

            topdata.CN2H = calcHubRotMat(zeros(3), azi_j)
            CN2H_no_azi = calcHubRotMat(zeros(3), 0.0)
            rbData = zeros(9)

            CH2N = LinearAlgebra.transpose(CN2H)

            runaero = false
            #################################################################
            if !isnothing(aero)
                if inputs.aeroLoadsOn > 0 #0 off, 1 one way, 1.5 one way with deformation from last timestep, 2 two way
                    runaero = true
                    if (inputs.aeroLoadsOn==1 || inputs.aeroLoadsOn==1.5) &&
                       numIterations!=1
                        runaero = false
                    end


                    if runaero
                        if inputs.AD15On
                            aeroVals, aeroDOFs = run_aero_with_deformAD15(
                                aero,
                                deformAero,
                                [topMesh],
                                [topEl],
                                [topdata],
                                [inputs],
                                t[i],
                            )
                            aeroVals = aeroVals[1]
                            aeroDOFs = aeroDOFs[1]
                        else
                            aeroVals, aeroDOFs = run_aero_with_deform(
                                aero,
                                deformAero,
                                topMesh,
                                topEl,
                                topdata.u_j,
                                topdata.uddot_j,
                                inputs,
                                numIterations,
                                t[i],
                                topdata.azi_j,
                                topdata.Omega_j,
                                topModel.gravityOn,
                            )
                        end
                    end
                end
            end

            #################################################################
            if inputs.aeroLoadsOn > 0 && isnan(maximum(aeroVals))
                @warn "Nan detected in aero forces"
            end
            if runaero || !isnothing(aeroVals)
                if inputs.aeroLoadsOn > 0
                    if isnothing(aeroVals)
                        error("aeroVals must be specified if OWENS.Inputs.aeroLoadsOn")
                    elseif isnothing(aeroDOFs)
                        error("aeroDOFs must be specified if OWENS.Inputs.aeroLoadsOn")
                    end

                    if inputs.AD15On
                        # AD15 is in global frame, so no frame conversion???
                        topdata.topFexternal = aeroVals
                        full_aeroDOFs = aeroDOFs
                    else
                        if length(size(aeroVals))==1 || size(aeroVals)[2]==1 #i.e. the standard aero force input as a long array
                            # Fill in forces and dofs if they were specified not in full arrays TODO: make this more efficient
                            full_aeroVals = zeros(eltype(aeroVals), topMesh.numNodes * 6)
                            for i_idx = 1:length(aeroDOFs)
                                full_aeroVals[Int(aeroDOFs[i_idx])] = aeroVals[i_idx]
                            end
                            full_aeroDOFs = collect(1:(topMesh.numNodes*6))
                            for iter_i = 1:floor(Int, length(full_aeroVals)/6)
                                topdata.topFexternal[(6*(iter_i-1)+1):(6*(iter_i-1)+6)] =
                                    frame_convert(
                                        full_aeroVals[(6*(iter_i-1)+1):(6*(iter_i-1)+6)],
                                        CN2H_no_azi,
                                    )
                            end
                        else # the other aero input as a 2D array
                            topdata.topFexternal = frame_convert(aeroVals[i+1, :], CN2H)
                        end
                    end
                else
                    topdata.topFexternal = zeros(numDOFPerNode)
                    full_aeroDOFs = copy(topdata.topFexternal) .* 0.0
                end
            else
                full_aeroDOFs = collect(1:(topMesh.numNodes*6))
                topdata.topFexternal = zero(full_aeroDOFs)
            end

            if meshcontrolfunction !== nothing
                # add to the loads based on the inputs, TODO: CN2H
                meshforces, meshdofs, timeconverged =
                    meshcontrolfunction(topMesh, u_j, t[i])
                for idx_main in full_aeroDOFs
                    for (idx, meshdof_idx) in enumerate(meshdofs)
                        if idx_main == meshdof_idx
                            topdata.topFexternal[idx_main] += meshforces[idx]
                        end
                    end
                end
            end

            #------------------------------------
            # TOPSIDE STRUCTURAL MODULE
            #------------------------------------

            if inputs.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j =
                    OWENSFEA.structuralDynamicsTransientROM(
                        topModel,
                        topMesh,
                        topEl,
                        topdata.topDispData1,
                        topdata.Omega_s,
                        topdata.OmegaDot_s,
                        t[i+1],
                        topdata.delta_t,
                        topElStorage,
                        topdata.top_rom,
                        topdata.topFexternal,
                        Int.(full_aeroDOFs),
                        topdata.CN2H,
                        topdata.rbData,
                    )
            elseif inputs.analysisType=="GX"
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j, systemout =
                    structuralDynamicsTransientGX(
                        topModel,
                        topMesh,
                        topdata.topFexternal,
                        Int.(full_aeroDOFs),
                        system,
                        assembly,
                        t,
                        topdata.Omega_j,
                        topdata.OmegaDot_j,
                        topdata.delta_t,
                        numIterations,
                        i,
                        strainGX,
                        curvGX,
                    )
            else # evalulate structural dynamics using conventional representation
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j =
                    OWENSFEA.structuralDynamicsTransient(
                        topModel,
                        topMesh,
                        topEl,
                        topdata.topDispData1,
                        topdata.Omega_s,
                        topdata.OmegaDot_s,
                        t[i+1],
                        topdata.delta_t,
                        topElStorage,
                        topdata.topFexternal,
                        Int.(full_aeroDOFs),
                        topdata.CN2H,
                        topdata.rbData;
                        predef = topModel.nlParams.predef,
                    )
            end

            u_jLast = copy(topdata.u_j)
            uddot_jLast = copy(topdata.u_j)
            topdata.u_j = topdata.topDispOut.displ_sp1
            topdata.udot_j = topdata.topDispOut.displdot_sp1
            topdata.uddot_j = topdata.topDispOut.displddot_sp1

            ## calculate norms
            uNorm =
                LinearAlgebra.norm(topdata.u_j .- u_jLast)/LinearAlgebra.norm(topdata.u_j)            #structural dynamics displacement iteration norm
            uddotNorm =
                LinearAlgebra.norm(
                    topdata.uddot_j .- uddot_jLast,
                )/LinearAlgebra.norm(topdata.uddot_j)            #structural dynamics displacement iteration norm
            aziNorm =
                LinearAlgebra.norm(
                    topdata.azi_j .- azi_jLast,
                )/LinearAlgebra.norm(topdata.azi_j)  #rotor azimuth iteration norm
            if inputs.generatorOn
                gbNorm =
                    LinearAlgebra.norm(
                        topdata.gb_j .- gb_jLast,
                    )/LinearAlgebra.norm(topdata.gb_j) #gearbox states iteration norm if it is off, the norm will be zero
            end

            numIterations = numIterations + 1

            if inputs.analysisType=="GX" && (
                !(uNorm > TOL || aziNorm > TOL || gbNorm > TOL) ||
                (numIterations >= MAXITER)
            )
                system = deepcopy(systemout)
            end

            # Strain stiffening, save at the end of the simulation, at the last while loop iteration, mutates elStorage
            if (i==numTS-1 || timeconverged == true) &&
               inputs.analysisType=="TNB" &&
               topModel.nlParams.predef=="update" &&
               (
                   !(uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) ||
                   (numIterations >= MAXITER)
               )
                OWENSFEA.structuralDynamicsTransient(
                    topModel,
                    topMesh,
                    topEl,
                    topdata.topDispData2,
                    topdata.Omega_s,
                    topdata.OmegaDot_s,
                    t[i+1],
                    topdata.delta_t,
                    topElStorage,
                    topdata.topFexternal,
                    Int.(full_aeroDOFs),
                    topdata.CN2H,
                    topdata.rbData;
                    predef = topModel.nlParams.predef,
                )
            end

            if inputs.verbosity>3
                println(
                    "$(numIterations) uNorm: $(uNorm) aziNorm: $(aziNorm) gbNorm: $(gbNorm) \n",
                )
            end

            if numIterations==MAXITER
                if inputs.verbosity>0
                    @warn "Maximum Iterations Met Breaking Iteration Loop"
                end
                break
            end

        end #end iteration while loop

        if inputs.verbosity >= 3
            avePower = mean(topdata.FReactionHist[:, 6] .* topdata.OmegaHist*(2*pi))
            instPower = mean(topdata.FReactionHist[i, 6] .* topdata.OmegaHist[i]*(2*pi))
            println("Gen Torque: $(topdata.genTorque_j)\n")
            println("Base Torque: $(topdata.FReactionHist[i,6])\n")
            println("RPM: $(topdata.Omega_j*60)\n")
            println("Vinf: $(newVinf)\n")
            println("Average Power: $(avePower)")
            println("Instant Power: $(instPower)")
            println("")
            # velocitymid = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)/2],t[i])
            # velocityquarter = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)/4],t[i])
            # println("Velocity mid: $(velocitymid[1])")
            # println("Velocity quarter: $(velocityquarter[1])")
        end

        if inputs.analysisType=="ROM"
            eta_j = topdata.topDispOut.eta_sp1 #eta_j
            etadot_j = topdata.topDispOut.etadot_sp1 #etadot_j
            etaddot_j = topdata.topDispOut.etaddot_sp1 #etaddot_j
            topdata.topDispData1 = OWENSFEA.DispData(
                topdata.u_j,
                topdata.udot_j,
                topdata.uddot_j,
                topdata.u_sm1,
                eta_j,
                etadot_j,
                etaddot_j,
            )
        else
            topdata.topDispData1 = OWENSFEA.DispData(
                topdata.u_j,
                topdata.udot_j,
                topdata.uddot_j,
                topdata.u_sm1,
            )
        end

        if isnan(maximum(u_j))
            @warn "Nan detected in displacements"
            break
        end

        ## update timestepping variables and other states, store in history arrays
        ## calculate converged generator torque/power

        genPowerPlot = topdata.genTorque_s*(gbDot_j*2*pi)*inputs.gearRatio

        topdata.u_sm1 = copy(topdata.u_s)
        topdata.u_s = topdata.u_j
        topdata.udot_s = topdata.udot_j
        topdata.uddot_s = topdata.uddot_j

        topdata.azi_s = topdata.azi_j
        topdata.Omega_s = topdata.Omega_j
        topdata.OmegaDot_s = topdata.OmegaDot_j

        topdata.genTorque_s = topdata.genTorque_j
        topdata.torqueDriveShaft_s = topdata.torqueDriveShaft_j

        topdata.gb_s = topdata.gb_j
        topdata.gbDot_s = topdata.gbDot_j
        topdata.gbDotDot_s = topdata.gbDotDot_j

        topdata.uHist[i+1, :] = topdata.u_s
        topdata.udotHist[i+1, :] = topdata.udot_s
        topdata.uddotHist[i+1, :] = topdata.uddot_s
        topdata.FReactionHist[i+1, :] = topdata.topFReaction_j
        topdata.topFexternal_hist[i+1, 1:length(topdata.topFexternal)] =
            topdata.topFexternal
        # FTwrBsHist[i+1,:] = -topFReaction_j  + topFWeight_j

        topdata.aziHist[i+1] = topdata.azi_s
        topdata.OmegaHist[i+1] = topdata.Omega_s
        topdata.OmegaDotHist[i+1] = topdata.OmegaDot_s

        topdata.gbHist[i+1] = topdata.gb_s
        topdata.gbDotHist[i+1] = topdata.gbDot_s
        topdata.gbDotDotHist[i+1] = topdata.gbDotDot_s

        #genTorque[i+1] = genTorque_s
        topdata.genTorque[i+1] = topdata.genTorque_s
        topdata.genPower[i+1] = genPowerPlot
        topdata.torqueDriveShaft[i+1] = topdata.torqueDriveShaft_s

        if inputs.analysisType=="GX"
            strainGX = topdata.topElStrain[1]
            curvGX = topdata.topElStrain[2]
            for ii = 1:length(strainGX[1, :])
                topdata.epsilon_x_hist[:, ii, i] .= strainGX[1, ii]
                topdata.kappa_y_hist[:, ii, i] .= curvGX[2, ii]
                topdata.kappa_z_hist[:, ii, i] .= curvGX[3, ii]
                topdata.epsilon_z_hist[:, ii, i] .= strainGX[3, ii]
                topdata.kappa_x_hist[:, ii, i] .= curvGX[1, ii]
                topdata.epsilon_y_hist[:, ii, i] .= strainGX[2, ii]
            end
        else
            for ii = 1:length(topdata.topElStrain)
                topdata.epsilon_x_hist[:, ii, i] = topdata.topElStrain[ii].epsilon_x
                topdata.kappa_y_hist[:, ii, i] = topdata.topElStrain[ii].kappa_y
                topdata.kappa_z_hist[:, ii, i] = topdata.topElStrain[ii].kappa_z
                topdata.epsilon_z_hist[:, ii, i] = topdata.topElStrain[ii].epsilon_z
                topdata.kappa_x_hist[:, ii, i] = topdata.topElStrain[ii].kappa_x
                topdata.epsilon_y_hist[:, ii, i] = topdata.topElStrain[ii].epsilon_y
            end
        end
        ## check rotor speed for generator operation
        if topdata.Omega_s >= topdata.rotorSpeedForGenStart
            inputs.generatorOn = true
        else
            inputs.generatorOn = false
        end

        if !isnothing(dataDumpFilename) && ((i-1)%datadumpfrequency==0 || i == numTS-1)
            println("\n Saving intermediate results to $dataDumpFilename \n")
            _write_restart_state(
                dataDumpFilename,
                topdata,
                i;
                inputs = inputs,
                system = system,
                capture_backend_state = !inputs.AD15On &&
                                        inputs.aeroLoadsOn > 0 &&
                                        !isnothing(aero),
            )
        end

    end #end timestep loop

    println("Simulation Complete.")


    if inputs.AD15On
        OWENSOpenFASTWrappers.endTurb()
    end

    if returnold
        return t[1:i],
        topdata.aziHist[1:i],
        topdata.OmegaHist[1:i],
        topdata.OmegaDotHist[1:i],
        topdata.gbHist[1:i],
        topdata.gbDotHist[1:i],
        topdata.gbDotDotHist[1:i],
        topdata.FReactionHist[1:i, :],
        topdata.FTwrBsHist[1:i, :],
        topdata.genTorque[1:i],
        topdata.genPower[1:i],
        topdata.torqueDriveShaft[1:i],
        topdata.uHist[1:i, :],
        topdata.uHist_prp[1:i, :],
        topdata.epsilon_x_hist[:, :, 1:i],
        topdata.epsilon_y_hist[:, :, 1:i],
        topdata.epsilon_z_hist[:, :, 1:i],
        topdata.kappa_x_hist[:, :, 1:i],
        topdata.kappa_y_hist[:, :, 1:i],
        topdata.kappa_z_hist[:, :, 1:i],
        topdata.FPtfmHist[1:i, :],
        topdata.FHydroHist[1:i, :],
        topdata.FMooringHist[1:i, :],
        topdata.topFexternal_hist[1:i, :],
        topdata.rbDataHist[1:i, :]
    else
        return topdata
    end
end

"""

run34m(inputs,feamodel,mymesh,myel,aeroForces,deformAero;steady=true,system=nothing,assembly=nothing,VTKFilename="./outvtk")

helper function that rearranges the outputs into the expected 34m output
    # Input
    * `inputs::Model`: see ?Model
    * `topModel::FEAModel`: see ?OWENSFEA.FEAModel
    * `mesh::Mesh`: see ?OWENSFEA.Mesh
    * `el::El`: see ?OWENSFEA.El
    * `aeroForces::function`: Fexternal, Fdof = aero(t) where Fexternal is the force on each affected mesh dof and Fdof is the corresponding DOFs affected
    * `deformAero::function`: see deformTurb(azi;newOmega=-1,newVinf=-1,bld_x=-1,bld_z=-1,bld_twist=-1,steady=false)
    * `steady::bool`: run steadystate with no aero or not
    * `system`: see ?GXBeam.System
    * `assembly`: see ?GXBeam.Assembly
    * `VTKFilename::string`: Unused: path and name of VTK output
    
    # Output
    * `eps_x`: strain history for eps_xx_0 for (Nbld,N_ts,mymesh.meshSeg[2]+1)
    * `eps_y`: strain history for eps_xx_z like above
    * `eps_z`: strain history for eps_xx_y like above
    * `kappa_x`: strain history for gam_xz_0 like above
    * `kappa_y`: strain history for gam_xz_y like above
    * `kappa_z`: strain history for gam_xy_0 like above
    * `t`: time array
    * `FReactionHist`: Nodal reaction 6dof forces history
    * `OmegaHist`: rotational speed array history
    * `genTorque`: generator torque history
    * `torqueDriveShaft`: driveshaft torque history
    * `aziHist`: azimuthal history array
    * `uHist`: mesh displacement history for each dof
    * `epsilon_x_hist`: strain history for eps_xx_0 for each dof
    * `epsilon_y_hist`: strain history for eps_xx_z for each dof
    * `epsilon_z_hist`: strain history for eps_xx_y for each dof
    * `kappa_x_hist`: strain history for gam_xz_0 for each dof
    * `kappa_y_hist`: strain history for gam_xz_y for each dof
    * `kappa_z_hist`: strain history for gam_xy_0 for each dof
    """
function run34m(
    inputs,
    feamodel,
    mymesh,
    myel,
    aeroForces,
    deformAero;
    steady = true,
    system = nothing,
    assembly = nothing,
    VTKFilename = "./outvtk",
)

    if !steady
        println("running unsteady")

        t,
        aziHist,
        OmegaHist,
        OmegaDotHist,
        gbHist,
        gbDotHist,
        gbDotDotHist,
        FReactionHist,
        FTwrBsHist,
        genTorque,
        genPower,
        torqueDriveShaft,
        uHist,
        uHist_prp,
        epsilon_x_hist,
        epsilon_y_hist,
        epsilon_z_hist,
        kappa_x_hist,
        kappa_y_hist,
        kappa_z_hist,
        FPtfmHist,
        FHydroHist,
        FMooringHist = OWENS.Unsteady_Land(
            inputs;
            topModel = feamodel,
            topMesh = mymesh,
            topEl = myel,
            aero = aeroForces,
            deformAero,
            system,
            assembly,
        )

        meanepsilon_z_hist = Statistics.mean(epsilon_z_hist, dims = 1)
        meanepsilon_y_hist = Statistics.mean(epsilon_y_hist, dims = 1)

    else
        println("here")
        println("running steady")

        feamodel.analysisType = "S"

        displ=zeros(mymesh.numNodes*6)
        elStorage = OWENS.OWENSFEA.initialElementCalculations(feamodel, myel, mymesh)
        displ, elStrain, staticAnalysisSuccessful, FReaction =
            OWENS.OWENSFEA.staticAnalysis(
                feamodel,
                mymesh,
                myel,
                displ,
                inputs.OmegaInit,
                inputs.OmegaInit,
                elStorage,
            )

        # format to match the unsteady method
        eps_x = [elStrain[i].epsilon_x[1] for i = 1:length(elStrain)]
        epsilon_x_hist = zeros(1, length(eps_x), 2)
        epsilon_x_hist[1, :, 1] = eps_x
        epsilon_x_hist[1, :, 2] = eps_x

        eps_y1_OW = [elStrain[i].epsilon_y[1] for i = 1:length(elStrain)]
        eps_y2_OW = [elStrain[i].epsilon_y[2] for i = 1:length(elStrain)]
        eps_y3_OW = [elStrain[i].epsilon_y[3] for i = 1:length(elStrain)]
        eps_y4_OW = [elStrain[i].epsilon_y[4] for i = 1:length(elStrain)]
        eps_y = (eps_y1_OW .+ eps_y2_OW .+ eps_y3_OW .+ eps_y4_OW) .* 0.25#0.34785484513745385
        meanepsilon_y_hist = zeros(1, length(eps_x), 2)
        meanepsilon_y_hist[1, :, 1] = eps_y
        meanepsilon_y_hist[1, :, 2] = eps_y

        eps_z1_OW = [elStrain[i].epsilon_z[1] for i = 1:length(elStrain)]
        eps_z2_OW = [elStrain[i].epsilon_z[2] for i = 1:length(elStrain)]
        eps_z3_OW = [elStrain[i].epsilon_z[3] for i = 1:length(elStrain)]
        eps_z4_OW = [elStrain[i].epsilon_z[4] for i = 1:length(elStrain)]
        eps_z = (eps_z1_OW .+ eps_z2_OW .+ eps_z3_OW .+ eps_z4_OW) .* 0.25#0.34785484513745385
        meanepsilon_z_hist = zeros(1, length(eps_x), 2)
        meanepsilon_z_hist[1, :, 1] = eps_z
        meanepsilon_z_hist[1, :, 2] = eps_z

        kappa_x = [elStrain[i].kappa_x[1] for i = 1:length(elStrain)]
        kappa_x_hist = zeros(1, length(eps_x), 2)
        kappa_x_hist[1, :, 1] = kappa_x
        kappa_x_hist[1, :, 2] = kappa_x

        kappa_y = [elStrain[i].kappa_y[1] for i = 1:length(elStrain)]
        kappa_y_hist = zeros(1, length(eps_x), 2)
        kappa_y_hist[1, :, 1] = kappa_y
        kappa_y_hist[1, :, 2] = kappa_y

        kappa_z = [elStrain[i].kappa_z[1] for i = 1:length(elStrain)]
        kappa_z_hist = zeros(1, length(eps_x), 2)
        kappa_z_hist[1, :, 1] = kappa_z
        kappa_z_hist[1, :, 2] = kappa_z

        FReactionHist = zeros(2, 6)
        FReactionHist[1, :] = FReaction[1:6]
        FReactionHist[2, :] = FReaction[1:6]

        OmegaHist = [inputs.OmegaInit, inputs.OmegaInit]
        genTorque = FReactionHist[:, 6]
        t = [0.0, 1.0]
        torqueDriveShaft = [0.0]
        aziHist = [0.0]
        uHist = [0.0]
    end


    # Interpolate the mesh strains onto the composite layup
    # TODO: or should we interpolate the composite stations onto the mesh?  It would be much more challenging
    Nbld = size(mymesh.structuralNodeNumbers)[1]
    N_ts = length(epsilon_x_hist[1, 1, :])
    eps_x = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)
    eps_z = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)
    eps_y = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)
    kappa_x = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)
    kappa_y = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)
    kappa_z = zeros(Nbld, N_ts, mymesh.meshSeg[2]+1)

    for ibld = 1:Nbld
        start = Int(mymesh.structuralElNumbers[ibld, 1])
        stop = Int(mymesh.structuralElNumbers[ibld, end-1])+1
        x = mymesh.z[start:stop]
        x = x .- x[1] #zero
        x = x ./ x[end] #normalize
        # samplepts = numadIn_bld.span./maximum(numadIn_bld.span) #normalize #TODO: this is spanwise, while everything else is vertical-wise
        for its = 1:N_ts
            #TODO: there are strain values at each quad point, should be better than just choosing one
            eps_x[ibld, its, :] = epsilon_x_hist[1, start:stop, its]#safeakima(x,epsilon_x_hist[1,start:stop,its],samplepts)
            eps_z[ibld, its, :] = meanepsilon_z_hist[1, start:stop, its]#safeakima(x,meanepsilon_z_hist[1,start:stop,its],samplepts)
            eps_y[ibld, its, :] = meanepsilon_y_hist[1, start:stop, its]#safeakima(x,meanepsilon_y_hist[1,start:stop,its],samplepts)
            kappa_x[ibld, its, :] = kappa_x_hist[1, start:stop, its]#safeakima(x,kappa_x_hist[1,start:stop,its],samplepts)
            kappa_y[ibld, its, :] = kappa_y_hist[1, start:stop, its]#safeakima(x,kappa_y_hist[1,start:stop,its],samplepts)
            kappa_z[ibld, its, :] = kappa_z_hist[1, start:stop, its]#safeakima(x,kappa_z_hist[1,start:stop,its],samplepts)
        end
    end

    # PyPlot.figure()
    # PyPlot.plot(t[1:end-1],eps_x[1,:,15],label="eps_x")
    # PyPlot.plot(t[1:end-1],eps_z[1,:,15],label="eps_z")
    # PyPlot.plot(t[1:end-1],eps_y[1,:,15],label="eps_y")
    # PyPlot.plot(t[1:end-1],kappa_x[1,:,15],label="kappa_x")
    # PyPlot.plot(t[1:end-1],kappa_y[1,:,15],label="kappa_y")
    # PyPlot.plot(t[1:end-1],kappa_z[1,:,15],label="kappa_z")
    #
    # PyPlot.plot(t[1:end-1],eps_x[2,:,15],":",label="eps_x2")
    # PyPlot.plot(t[1:end-1],eps_z[2,:,15],":",label="eps_z2")
    # PyPlot.plot(t[1:end-1],eps_y[2,:,15],":",label="eps_y2")
    # PyPlot.plot(t[1:end-1],kappa_x[2,:,15],":",label="kappa_x2")
    # PyPlot.plot(t[1:end-1],kappa_y[2,:,15],":",label="kappa_y2")
    # PyPlot.plot(t[1:end-1],kappa_z[2,:,15],":",label="kappa_z2")
    # PyPlot.legend()

    return eps_x,
    eps_z,
    eps_y,
    kappa_x,
    kappa_y,
    kappa_z,
    t,
    FReactionHist,
    OmegaHist,
    genTorque,
    torqueDriveShaft,
    aziHist,
    uHist,
    epsilon_x_hist,
    meanepsilon_y_hist,
    meanepsilon_z_hist,
    kappa_x_hist,
    kappa_y_hist,
    kappa_z_hist
end
