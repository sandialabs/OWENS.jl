export HAWTOyeState,
    initializeHAWTOyeState, runHAWTCCBladeOyeStep, runHAWTUnsteadyAeroelastic

"""
    HAWTOyeState

Oye dynamic-inflow state for the native OWENS HAWT/CCBlade coupling. The
state stores reduced and dynamic axial induction plus the same pair for
tangential induction at the CCBlade radial stations.
"""
mutable struct HAWTOyeState
    axial_reduced::Vector{Float64}
    axial_dynamic::Vector{Float64}
    tangential_reduced::Vector{Float64}
    tangential_dynamic::Vector{Float64}
end

function _hawt_station_vector(values, name; positive = false)
    values isa AbstractVector || throw(ArgumentError("$name must be a vector"))
    isempty(values) && throw(ArgumentError("$name must not be empty"))
    all(x -> x isa Real && isfinite(x), values) ||
        throw(ArgumentError("$name must contain only finite real values"))
    if positive && any(x -> x <= 0.0, values)
        throw(ArgumentError("$name must contain only finite positive values"))
    end
    return Float64.(collect(values))
end

function _hawt_matching_station_vectors(radial_positions, chord, twist)
    r = _hawt_station_vector(radial_positions, "radial_positions"; positive = true)
    c = _hawt_station_vector(chord, "chord"; positive = true)
    theta = _hawt_station_vector(twist, "twist")
    length(r) == length(c) == length(theta) ||
        throw(ArgumentError("radial_positions, chord, and twist must have the same length"))
    all(diff(r) .> 0.0) ||
        throw(ArgumentError("radial_positions must be strictly increasing"))
    return r, c, theta
end

function _hawt_checked_state(state::HAWTOyeState, nstations)
    for field in (:axial_reduced, :axial_dynamic, :tangential_reduced, :tangential_dynamic)
        values = getfield(state, field)
        length(values) == nstations || throw(
            ArgumentError("HAWTOyeState.$field must have one value per radial station"),
        )
        all(isfinite, values) ||
            throw(ArgumentError("HAWTOyeState.$field must contain only finite values"))
    end
    return state
end

function _copy_hawt_oye_state(state::HAWTOyeState)
    return HAWTOyeState(
        copy(state.axial_reduced),
        copy(state.axial_dynamic),
        copy(state.tangential_reduced),
        copy(state.tangential_dynamic),
    )
end

"""
    initializeHAWTOyeState(radial_positions; axial_induction=0.0,
                           tangential_induction=0.0)

Create an Oye dynamic-inflow state sized to the supplied HAWT radial stations.
The default starts from zero induction, which is useful for startup/transient
smoke tests. Pass nonzero `axial_induction` or `tangential_induction` to start
from a known operating point.
"""
function initializeHAWTOyeState(
    radial_positions;
    axial_induction = 0.0,
    tangential_induction = 0.0,
)
    r = _hawt_station_vector(radial_positions, "radial_positions"; positive = true)
    axial = _hawt_initial_induction_vector(axial_induction, length(r), "axial_induction")
    tangential = _hawt_initial_induction_vector(
        tangential_induction,
        length(r),
        "tangential_induction",
    )
    return HAWTOyeState(copy(axial), copy(axial), copy(tangential), copy(tangential))
end

function _hawt_initial_induction_vector(value, nstations, name)
    if value isa AbstractVector
        values = _hawt_station_vector(value, name)
        length(values) == nstations ||
            throw(ArgumentError("$name must have one value per radial station"))
        return values
    end
    value isa Real && isfinite(value) ||
        throw(ArgumentError("$name must be a finite scalar or vector"))
    return fill(Float64(value), nstations)
end

function _hawt_value_at(value, index, time, ntimes, name)
    selected = if value isa Function
        value(time)
    elseif value isa AbstractVector
        length(value) == ntimes ||
            throw(ArgumentError("$name vector must have one value per time step"))
        value[index]
    else
        value
    end
    selected isa Real && isfinite(selected) ||
        throw(ArgumentError("$name must evaluate to a finite real value"))
    return Float64(selected)
end

function _hawt_dof_vector(values, total_dof, name; default = 0.0)
    if isnothing(values)
        return fill(Float64(default), total_dof)
    end
    length(values) == total_dof ||
        throw(ArgumentError("$name must have length mesh.numNodes * 6"))
    all(isfinite, values) || throw(ArgumentError("$name must contain only finite values"))
    return Float64.(collect(values))
end

function _hawt_time_vector(times)
    t = _hawt_station_vector(times, "times")
    all(diff(t) .> 0.0) || throw(ArgumentError("times must be strictly increasing"))
    return t
end

function _hawt_deformed_station_positions(
    mesh,
    radial_positions;
    displacements,
    hub_position,
    rotor_axis,
)
    nominal_node_radii = vec(mean(hawtStructuralRadialStations(mesh; rotor_axis), dims = 1))
    current_node_radii = vec(
        mean(
            hawtStructuralRadialStations(mesh; displacements, hub_position, rotor_axis),
            dims = 1,
        ),
    )
    return _hawt_linear_interpolate(
        nominal_node_radii,
        current_node_radii,
        radial_positions,
        "radial_positions",
    )
end

function _hawt_relative_speed_load_scale(result, dynamic_axial, dynamic_tangential)
    quasi_axial = Float64.(collect(result.a))
    quasi_tangential = Float64.(collect(result.ap))
    vx = [Float64(op.Vx) for op in result.operating_points]
    vy = [Float64(op.Vy) for op in result.operating_points]

    quasi_speed2 = @. ((1.0 - quasi_axial) * vx)^2 + ((1.0 + quasi_tangential) * vy)^2
    dynamic_speed2 = @. ((1.0 - dynamic_axial) * vx)^2 + ((1.0 + dynamic_tangential) * vy)^2

    return dynamic_speed2 ./ max.(quasi_speed2, eps(Float64))
end

"""
    runHAWTCCBladeOyeStep(mesh, radial_positions, chord, twist, airfoils, state,
                          dt, rotor_speed, inflow_speed, rho; kwargs...)

Run one native HAWT aerodynamic coupling step. The step solves a rigid CCBlade
HAWT operating point, advances Oye axial/tangential induction states, applies a
first-order relative-speed load lag, and maps the resulting normal/tangential
loads onto OWENS structural DOFs with `mapHAWTCCBladeLoads`.

This is the OWENS-side native unsteady HAWT smoke path. It is intentionally
separate from the AeroDyn wrapper path and is not a replacement for validated
AeroDyn DBEMT parity studies.
"""
function runHAWTCCBladeOyeStep(
    mesh,
    radial_positions,
    chord,
    twist,
    airfoils,
    state::HAWTOyeState,
    dt,
    rotor_speed,
    inflow_speed,
    rho;
    num_blades = size(mesh.structuralNodeNumbers, 1),
    hub_radius = 0.0,
    tip_radius = maximum(radial_positions),
    pitch = 0.0,
    precone = 0.0,
    mu = 1.7894e-5,
    asound = 340.0,
    npts = 10,
    tip_correction = OWENSAero.CCBlade.PrandtlTipHub(),
    displacements = nothing,
    hub_position = nothing,
    rotor_axis = :x,
    dynamic_inflow = true,
)
    r, c, theta = _hawt_matching_station_vectors(radial_positions, chord, twist)
    _hawt_checked_state(state, length(r))
    dt isa Real && isfinite(dt) && dt >= 0.0 ||
        throw(ArgumentError("dt must be a finite nonnegative real value"))

    station_positions =
        _hawt_deformed_station_positions(mesh, r; displacements, hub_position, rotor_axis)
    all(station_positions .> 0.0) ||
        throw(ArgumentError("deformed radial positions must remain positive"))
    all(diff(station_positions) .> 0.0) ||
        throw(ArgumentError("deformed radial positions must remain strictly increasing"))
    current_tip_radius =
        maximum(hawtStructuralRadialStations(mesh; displacements, hub_position, rotor_axis))
    effective_tip_radius = max(tip_radius, last(station_positions), current_tip_radius)

    steady = OWENSAero.ccbladeHAWTSolve(
        station_positions,
        c,
        theta,
        airfoils,
        rotor_speed,
        inflow_speed,
        rho;
        num_blades,
        hub_radius,
        tip_radius = effective_tip_radius,
        pitch,
        precone,
        mu,
        asound,
        npts,
        tip_correction,
    )

    if dynamic_inflow
        tau1, tau2 = OWENSAero.oyeDynamicInflowTimeConstants(
            effective_tip_radius,
            inflow_speed,
            mean(steady.a),
            station_positions,
        )
        state.axial_reduced, state.axial_dynamic = OWENSAero.oyeDynamicInflowStep(
            state.axial_reduced,
            state.axial_dynamic,
            steady.a,
            dt,
            tau1,
            tau2,
        )
        state.tangential_reduced, state.tangential_dynamic = OWENSAero.oyeDynamicInflowStep(
            state.tangential_reduced,
            state.tangential_dynamic,
            steady.ap,
            dt,
            tau1,
            tau2,
        )
        load_scale = _hawt_relative_speed_load_scale(
            steady,
            state.axial_dynamic,
            state.tangential_dynamic,
        )
    else
        tau1, tau2 = NaN, fill(NaN, length(r))
        state.axial_reduced .= steady.a
        state.axial_dynamic .= steady.a
        state.tangential_reduced .= steady.ap
        state.tangential_dynamic .= steady.ap
        load_scale = ones(length(r))
    end

    normal_loads = Float64.(steady.Np) .* load_scale
    tangential_loads = Float64.(steady.Tp) .* load_scale
    force_values, force_dofs = mapHAWTCCBladeLoads(
        mesh,
        station_positions,
        normal_loads,
        tangential_loads;
        hub_radius,
        tip_radius = effective_tip_radius,
        displacements,
        hub_position,
        rotor_axis,
    )
    resultants =
        hawtNodalLoadResultants(mesh, force_values; displacements, hub_position, rotor_axis)

    return (
        steady = steady,
        state = _copy_hawt_oye_state(state),
        station_positions = station_positions,
        tau1 = tau1,
        tau2 = tau2,
        load_scale = load_scale,
        normal_loads = normal_loads,
        tangential_loads = tangential_loads,
        force_values = force_values,
        force_dofs = force_dofs,
        resultants = resultants,
    )
end

"""
    runHAWTUnsteadyAeroelastic(mesh, el, feamodel, radial_positions, chord,
                               twist, airfoils, times; kwargs...)

Run a compact native HAWT aeroelastic time history. Each step calls
`runHAWTCCBladeOyeStep`; when `structural_coupling=:static`, the mapped
aerodynamic load vector is passed to `OWENSFEA.staticAnalysis`. When
`structural_coupling=:transient`, the load vector advances
`OWENSFEA.structuralDynamicsTransient` for each positive time increment. In
both cases the next aero step sees the deformed mesh through the displacement
vector.

`rotor_speed`, `inflow_speed`, and `pitch` may be scalars, vectors with one
entry per time step, or functions of time. Rotor speed is in rad/s.
"""
function runHAWTUnsteadyAeroelastic(
    mesh,
    el,
    feamodel,
    radial_positions,
    chord,
    twist,
    airfoils,
    times;
    rotor_speed,
    inflow_speed,
    rho,
    num_blades = size(mesh.structuralNodeNumbers, 1),
    hub_radius = 0.0,
    tip_radius = maximum(radial_positions),
    pitch = 0.0,
    precone = 0.0,
    mu = 1.7894e-5,
    asound = 340.0,
    npts = 10,
    tip_correction = OWENSAero.CCBlade.PrandtlTipHub(),
    rotor_axis = :x,
    initial_state = nothing,
    initial_displacements = nothing,
    initial_velocities = nothing,
    initial_accelerations = nothing,
    initial_previous_displacements = nothing,
    structural_coupling = :static,
    structural_relaxation = 1.0,
    dynamic_inflow = true,
    CN2H = LinearAlgebra.I(3),
    rbData = zeros(9),
)
    t = _hawt_time_vector(times)
    r, c, theta = _hawt_matching_station_vectors(radial_positions, chord, twist)
    structural_relaxation isa Real && 0.0 < structural_relaxation <= 1.0 ||
        throw(ArgumentError("structural_relaxation must be in (0, 1]"))

    state =
        isnothing(initial_state) ? initializeHAWTOyeState(r) :
        _copy_hawt_oye_state(_hawt_checked_state(initial_state, length(r)))

    total_dof = mesh.numNodes * 6
    displacements =
        _hawt_dof_vector(initial_displacements, total_dof, "initial_displacements")
    velocities = _hawt_dof_vector(initial_velocities, total_dof, "initial_velocities")
    accelerations =
        _hawt_dof_vector(initial_accelerations, total_dof, "initial_accelerations")
    previous_displacements =
        isnothing(initial_previous_displacements) ? copy(displacements) :
        _hawt_dof_vector(
            initial_previous_displacements,
            total_dof,
            "initial_previous_displacements",
        )

    if structural_coupling == :static || structural_coupling == :transient
        isnothing(el) && throw(ArgumentError("el is required for HAWT structural coupling"))
        isnothing(feamodel) &&
            throw(ArgumentError("feamodel is required for HAWT structural coupling"))
        el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
    elseif structural_coupling == :none
        el_storage = nothing
    else
        throw(ArgumentError("structural_coupling must be :static, :transient, or :none"))
    end
    if structural_coupling == :transient &&
       !(feamodel.analysisType in ("TNB", "TD", "stiff"))
        throw(
            ArgumentError(
                "transient HAWT coupling requires feamodel.analysisType to be TNB, TD, or stiff",
            ),
        )
    end
    CN2H_matrix = Matrix{Float64}(CN2H)
    size(CN2H_matrix) == (3, 3) || throw(ArgumentError("CN2H must be a 3x3 matrix"))
    rb_data = _hawt_station_vector(rbData, "rbData")
    length(rb_data) == 9 || throw(ArgumentError("rbData must have length 9"))

    ntimes = length(t)
    load_history = zeros(Float64, total_dof, ntimes)
    displacement_history = zeros(Float64, total_dof, ntimes)
    velocity_history = zeros(Float64, total_dof, ntimes)
    acceleration_history = zeros(Float64, total_dof, ntimes)
    reaction_history = zeros(Float64, total_dof, ntimes)
    resultant_force_history = zeros(Float64, 3, ntimes)
    resultant_moment_history = zeros(Float64, 3, ntimes)
    station_history = zeros(Float64, length(r), ntimes)
    normal_load_history = zeros(Float64, length(r), ntimes)
    tangential_load_history = zeros(Float64, length(r), ntimes)
    axial_dynamic_history = zeros(Float64, length(r), ntimes)
    tangential_dynamic_history = zeros(Float64, length(r), ntimes)
    structural_success = trues(ntimes)
    step_results = Vector{Any}(undef, ntimes)
    strain_history = Vector{Any}(undef, ntimes)

    for itime = 1:ntimes
        time = t[itime]
        dt = itime == 1 ? 0.0 : time - t[itime-1]
        omega = _hawt_value_at(rotor_speed, itime, time, ntimes, "rotor_speed")
        wind = _hawt_value_at(inflow_speed, itime, time, ntimes, "inflow_speed")
        pitch_i = _hawt_value_at(pitch, itime, time, ntimes, "pitch")

        step = runHAWTCCBladeOyeStep(
            mesh,
            r,
            c,
            theta,
            airfoils,
            state,
            dt,
            omega,
            wind,
            rho;
            num_blades,
            hub_radius,
            tip_radius,
            pitch = pitch_i,
            precone,
            mu,
            asound,
            npts,
            tip_correction,
            displacements,
            rotor_axis,
            dynamic_inflow,
        )

        if structural_coupling == :static
            omega_hz = omega / (2.0 * pi)
            new_displacements, _, success, _ = redirect_stdout(devnull) do
                OWENSFEA.staticAnalysis(
                    feamodel,
                    mesh,
                    el,
                    displacements,
                    omega_hz,
                    omega_hz,
                    el_storage;
                    Fdof = step.force_dofs,
                    Fexternal = step.force_values,
                )
            end
            structural_success[itime] = success
            displacements .=
                (1.0 - structural_relaxation) .* displacements .+
                structural_relaxation .* new_displacements
            velocities .= 0.0
            accelerations .= 0.0
        elseif structural_coupling == :transient && itime > 1
            previous_omega =
                _hawt_value_at(rotor_speed, itime - 1, t[itime-1], ntimes, "rotor_speed")
            omega_hz = omega / (2.0 * pi)
            omega_dot_hz = (omega - previous_omega) / (2.0 * pi * dt)
            disp_data = OWENSFEA.DispData(
                displacements,
                velocities,
                accelerations,
                previous_displacements,
            )
            el_strain, disp_out, reaction = redirect_stdout(devnull) do
                OWENSFEA.structuralDynamicsTransient(
                    feamodel,
                    mesh,
                    el,
                    disp_data,
                    omega_hz,
                    omega_dot_hz,
                    time,
                    dt,
                    el_storage,
                    step.force_values,
                    Int.(step.force_dofs),
                    CN2H_matrix,
                    rb_data;
                    predef = feamodel.nlParams.predef,
                )
            end
            previous_displacements .= displacements
            displacements .= disp_out.displ_sp1
            velocities .= disp_out.displdot_sp1
            accelerations .= disp_out.displddot_sp1
            reaction_history[:, itime] .= reaction
            strain_history[itime] = el_strain
        end

        load_history[:, itime] .= step.force_values
        displacement_history[:, itime] .= displacements
        velocity_history[:, itime] .= velocities
        acceleration_history[:, itime] .= accelerations
        resultant_force_history[:, itime] .= step.resultants.force
        resultant_moment_history[:, itime] .= step.resultants.moment
        station_history[:, itime] .= step.station_positions
        normal_load_history[:, itime] .= step.normal_loads
        tangential_load_history[:, itime] .= step.tangential_loads
        axial_dynamic_history[:, itime] .= step.state.axial_dynamic
        tangential_dynamic_history[:, itime] .= step.state.tangential_dynamic
        step_results[itime] = step
    end

    return (
        times = t,
        final_state = _copy_hawt_oye_state(state),
        load_history = load_history,
        displacement_history = displacement_history,
        velocity_history = velocity_history,
        acceleration_history = acceleration_history,
        reaction_history = reaction_history,
        resultant_force_history = resultant_force_history,
        resultant_moment_history = resultant_moment_history,
        station_history = station_history,
        normal_load_history = normal_load_history,
        tangential_load_history = tangential_load_history,
        axial_dynamic_history = axial_dynamic_history,
        tangential_dynamic_history = tangential_dynamic_history,
        structural_success = structural_success,
        strain_history = strain_history,
        step_results = step_results,
    )
end
