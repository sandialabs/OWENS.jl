function _hawt_axis_index(rotor_axis)
    if rotor_axis == :x || rotor_axis == 1
        return 1
    elseif rotor_axis == :y || rotor_axis == 2
        return 2
    elseif rotor_axis == :z || rotor_axis == 3
        return 3
    end
    throw(ArgumentError("rotor_axis must be :x, :y, :z, 1, 2, or 3"))
end

function _hawt_radial_indices(rotor_axis)
    axis_index = _hawt_axis_index(rotor_axis)
    return Tuple(i for i = 1:3 if i != axis_index)
end

function _hawt_axis_unit(rotor_axis)
    axis_unit = zeros(Float64, 3)
    axis_unit[_hawt_axis_index(rotor_axis)] = 1.0
    return axis_unit
end

function _hawt_displacement_vector(displacements, mesh)
    isnothing(displacements) && return nothing
    length(displacements) == mesh.numNodes * 6 ||
        throw(ArgumentError("displacements must have length mesh.numNodes * 6"))
    return displacements
end

function _hawt_node_position(mesh, node_number, displacements)
    position = [mesh.x[node_number], mesh.y[node_number], mesh.z[node_number]]
    if !isnothing(displacements)
        dof0 = 6 * (node_number - 1)
        position .+= displacements[(dof0+1):(dof0+3)]
    end
    return position
end

function _hawt_default_hub_position(mesh, displacements)
    root_nodes = Int.(mesh.structuralNodeNumbers[:, 1])
    hub_position = zeros(Float64, 3)
    for node_number in root_nodes
        hub_position .+= _hawt_node_position(mesh, node_number, displacements)
    end
    return hub_position ./ length(root_nodes)
end

function _hawt_hub_position(mesh, displacements, hub_position)
    isnothing(hub_position) && return _hawt_default_hub_position(mesh, displacements)
    length(hub_position) == 3 || throw(ArgumentError("hub_position must have length 3"))
    all(isfinite, hub_position) || throw(ArgumentError("hub_position must be finite"))
    return Float64.(collect(hub_position))
end

"""
    hawtStructuralRadialStations(mesh; displacements=nothing, hub_position=nothing,
                                 rotor_axis=:x)

Return the HAWT blade-node radial stations implied by an OWENS structural mesh.
The default convention uses the hub `x` axis as the shaft axis, so radial
distance is measured in the `y-z` rotor plane. Legacy meshes that place the
rotor in another plane can pass `rotor_axis=:y` or `rotor_axis=:z`.
Translational structural displacements may be supplied as a full OWENSFEA
`6 * mesh.numNodes` vector.
"""
function hawtStructuralRadialStations(
    mesh;
    displacements = nothing,
    hub_position = nothing,
    rotor_axis = :x,
)
    radial_indices = _hawt_radial_indices(rotor_axis)
    displacements = _hawt_displacement_vector(displacements, mesh)
    hub_position = _hawt_hub_position(mesh, displacements, hub_position)

    nblade, nnodes_per_blade = size(mesh.structuralNodeNumbers)
    radial_stations = zeros(Float64, nblade, nnodes_per_blade)
    for iblade = 1:nblade
        for inode = 1:nnodes_per_blade
            node_number = Int(mesh.structuralNodeNumbers[iblade, inode])
            position = _hawt_node_position(mesh, node_number, displacements)
            radial_stations[iblade, inode] = hypot(
                position[radial_indices[1]] - hub_position[radial_indices[1]],
                position[radial_indices[2]] - hub_position[radial_indices[2]],
            )
        end
    end
    return radial_stations
end

function _validate_hawt_station_loads(
    radial_positions,
    normal_loads,
    tangential_loads,
    nblade,
)
    radial_positions isa AbstractVector ||
        throw(ArgumentError("radial_positions must be a vector"))
    isempty(radial_positions) && throw(ArgumentError("radial_positions must not be empty"))
    r = Float64.(collect(radial_positions))
    all(isfinite, r) || throw(ArgumentError("radial_positions must be finite"))
    minimum(r) >= 0.0 || throw(ArgumentError("radial_positions must be nonnegative"))
    all(diff(r) .> 0.0) ||
        throw(ArgumentError("radial_positions must be strictly increasing"))

    normal = _hawt_load_matrix(normal_loads, nblade, length(r), "normal_loads")
    tangential = _hawt_load_matrix(tangential_loads, nblade, length(r), "tangential_loads")
    return r, normal, tangential
end

function _hawt_load_matrix(loads, nblade, nstations, name)
    all(isfinite, loads) || throw(ArgumentError("$name must be finite"))
    if loads isa AbstractVector
        length(loads) == nstations ||
            throw(ArgumentError("$name must have one value per radial station"))
        return repeat(reshape(Float64.(collect(loads)), 1, nstations), nblade, 1)
    elseif loads isa AbstractMatrix
        size(loads) == (nblade, nstations) || throw(
            ArgumentError("$name must have size (number of blades, number of stations)"),
        )
        return Float64.(loads)
    else
        throw(ArgumentError("$name must be a vector or matrix"))
    end
end

function _hawt_support_points(radial_positions, loads, hub_radius, tip_radius)
    hub_radius >= 0.0 || throw(ArgumentError("hub_radius must be nonnegative"))
    tip_radius >= last(radial_positions) ||
        throw(ArgumentError("tip_radius must be at least the last radial station"))
    first(radial_positions) >= hub_radius ||
        throw(ArgumentError("radial_positions must not be inside hub_radius"))

    support_r = Float64[]
    support_loads = Float64[]
    if hub_radius > 0.0
        push!(support_r, 0.0)
        push!(support_loads, 0.0)
        if hub_radius < first(radial_positions)
            push!(support_r, hub_radius)
            push!(support_loads, 0.0)
        end
    elseif first(radial_positions) > 0.0
        push!(support_r, 0.0)
        push!(support_loads, 0.0)
    end

    append!(support_r, radial_positions)
    append!(support_loads, loads)

    if tip_radius > last(radial_positions)
        push!(support_r, tip_radius)
        push!(support_loads, 0.0)
    end
    return support_r, support_loads
end

function _hawt_linear_interpolate(x, y, targets, name)
    all(diff(x) .> 0.0) || throw(
        ArgumentError("$name interpolation support points must be strictly increasing"),
    )
    values = zeros(Float64, length(targets))
    tolerance = 100 * eps(Float64) * max(1.0, maximum(abs.(x)))
    for (itarget, target) in enumerate(targets)
        if target < first(x) - tolerance || target > last(x) + tolerance
            throw(
                ArgumentError(
                    "$name interpolation target $target is outside $(first(x)) to $(last(x))",
                ),
            )
        elseif target <= first(x) + tolerance
            values[itarget] = first(y)
        elseif target >= last(x) - tolerance
            values[itarget] = last(y)
        else
            ilower = searchsortedlast(x, target)
            if x[ilower] == target
                values[itarget] = y[ilower]
            else
                fraction = (target - x[ilower]) / (x[ilower+1] - x[ilower])
                values[itarget] = y[ilower] + fraction * (y[ilower+1] - y[ilower])
            end
        end
    end
    return values
end

function _hawt_node_unit_vectors(node_positions, hub_position, rotor_axis)
    nnodes = size(node_positions, 2)
    axis_index = _hawt_axis_index(rotor_axis)
    axis_unit = _hawt_axis_unit(rotor_axis)
    radial_indices = _hawt_radial_indices(rotor_axis)
    radial_units = zeros(Float64, 3, nnodes)
    tangential_units = zeros(Float64, 3, nnodes)

    fallback = nothing
    for inode = 1:nnodes
        radial_components = [
            node_positions[radial_indices[1], inode] - hub_position[radial_indices[1]],
            node_positions[radial_indices[2], inode] - hub_position[radial_indices[2]],
        ]
        radius = hypot(radial_components[1], radial_components[2])
        if radius > 100 * eps(Float64)
            fallback = zeros(Float64, 3)
            fallback[radial_indices[1]] = radial_components[1] / radius
            fallback[radial_indices[2]] = radial_components[2] / radius
            break
        end
    end
    if isnothing(fallback)
        fallback = zeros(Float64, 3)
        fallback[radial_indices[1]] = 1.0
    end

    for inode = 1:nnodes
        radial_components = [
            node_positions[radial_indices[1], inode] - hub_position[radial_indices[1]],
            node_positions[radial_indices[2], inode] - hub_position[radial_indices[2]],
        ]
        radius = hypot(radial_components[1], radial_components[2])
        radial_unit = copy(fallback)
        if radius > 100 * eps(Float64)
            radial_unit .= 0.0
            radial_unit[radial_indices[1]] = radial_components[1] / radius
            radial_unit[radial_indices[2]] = radial_components[2] / radius
        end
        radial_unit[axis_index] = 0.0
        radial_units[:, inode] .= radial_unit
        tangential_units[:, inode] .= LinearAlgebra.cross(axis_unit, radial_unit)
    end
    return radial_units, tangential_units
end

"""
    mapHAWTCCBladeLoads(mesh, radial_positions, normal_loads, tangential_loads;
                        hub_radius=0.0, tip_radius=maximum(radial_positions),
                        displacements=nothing, hub_position=nothing,
                        rotor_axis=:x)

Map CCBlade-style HAWT blade distributed loads onto OWENSFEA nodal degrees of
freedom. `normal_loads` are force per span along the positive shaft axis.
`tangential_loads` are force per span in the positive rotor-rotation direction
about that axis. Loads may be one vector shared by all blades or a
`nblade x nstation` matrix.

The return value is `(ForceValHist, ForceDof)`, matching the existing OWENS
aero-load mapping convention. `ForceValHist` is a full `6 * mesh.numNodes`
vector so it can be passed directly to OWENSFEA transient/static solves.
The default `rotor_axis=:x` matches the intended HAWT hub convention; legacy
HAWT meshes generated in the `x-y` plane should pass `rotor_axis=:z`.
"""
function mapHAWTCCBladeLoads(
    mesh,
    radial_positions,
    normal_loads,
    tangential_loads;
    hub_radius = 0.0,
    tip_radius = maximum(radial_positions),
    displacements = nothing,
    hub_position = nothing,
    rotor_axis = :x,
)
    _hawt_axis_index(rotor_axis)
    displacements = _hawt_displacement_vector(displacements, mesh)
    hub_position = _hawt_hub_position(mesh, displacements, hub_position)

    nblade, nnodes_per_blade = size(mesh.structuralNodeNumbers)
    r, normal, tangential = _validate_hawt_station_loads(
        radial_positions,
        normal_loads,
        tangential_loads,
        nblade,
    )

    ForceValHist = zeros(Float64, mesh.numNodes * 6)
    ForceDof = collect(1:(mesh.numNodes*6))
    axial_unit = _hawt_axis_unit(rotor_axis)

    for iblade = 1:nblade
        node_numbers = Int.(mesh.structuralNodeNumbers[iblade, :])
        node_positions = zeros(Float64, 3, nnodes_per_blade)
        for (inode, node_number) in enumerate(node_numbers)
            node_positions[:, inode] .=
                _hawt_node_position(mesh, node_number, displacements)
        end

        node_radii = vec(
            hawtStructuralRadialStations(mesh; displacements, hub_position, rotor_axis)[
                iblade,
                :,
            ],
        )
        normal_r, normal_support =
            _hawt_support_points(r, normal[iblade, :], hub_radius, tip_radius)
        tangential_r, tangential_support =
            _hawt_support_points(r, tangential[iblade, :], hub_radius, tip_radius)
        normal_at_nodes =
            _hawt_linear_interpolate(normal_r, normal_support, node_radii, "normal_loads")
        tangential_at_nodes = _hawt_linear_interpolate(
            tangential_r,
            tangential_support,
            node_radii,
            "tangential_loads",
        )
        _, tangential_units =
            _hawt_node_unit_vectors(node_positions, hub_position, rotor_axis)

        for isegment = 1:(nnodes_per_blade-1)
            node1 = node_numbers[isegment]
            node2 = node_numbers[isegment+1]
            segment_vector = node_positions[:, isegment+1] - node_positions[:, isegment]
            segment_length = LinearAlgebra.norm(segment_vector)
            segment_length > 0.0 ||
                throw(ArgumentError("HAWT blade segment length must be positive"))

            q1 =
                normal_at_nodes[isegment] .* axial_unit .+
                tangential_at_nodes[isegment] .* tangential_units[:, isegment]
            q2 =
                normal_at_nodes[isegment+1] .* axial_unit .+
                tangential_at_nodes[isegment+1] .* tangential_units[:, isegment+1]

            f1 = segment_length .* (2.0 .* q1 .+ q2) ./ 6.0
            f2 = segment_length .* (q1 .+ 2.0 .* q2) ./ 6.0

            dof1 = 6 * (node1 - 1)
            dof2 = 6 * (node2 - 1)
            ForceValHist[(dof1+1):(dof1+3)] .+= f1
            ForceValHist[(dof2+1):(dof2+3)] .+= f2
        end
    end

    return ForceValHist, ForceDof
end

"""
    hawtNodalLoadResultants(mesh, force_values; displacements=nothing,
                            hub_position=nothing, rotor_axis=:x)

Return total force and moment about the HAWT hub for a full OWENS nodal load
vector. This is intended for coupling verification and sign-convention checks.
"""
function hawtNodalLoadResultants(
    mesh,
    force_values;
    displacements = nothing,
    hub_position = nothing,
    rotor_axis = :x,
)
    _hawt_axis_index(rotor_axis)
    length(force_values) == mesh.numNodes * 6 ||
        throw(ArgumentError("force_values must have length mesh.numNodes * 6"))
    displacements = _hawt_displacement_vector(displacements, mesh)
    hub_position = _hawt_hub_position(mesh, displacements, hub_position)

    total_force = zeros(Float64, 3)
    total_moment = zeros(Float64, 3)
    for node_number = 1:mesh.numNodes
        dof0 = 6 * (node_number - 1)
        force = Float64.(collect(force_values[(dof0+1):(dof0+3)]))
        moment = Float64.(collect(force_values[(dof0+4):(dof0+6)]))
        position = _hawt_node_position(mesh, node_number, displacements)
        total_force .+= force
        total_moment .+= moment .+ LinearAlgebra.cross(position .- hub_position, force)
    end
    return (force = total_force, moment = total_moment)
end

"""
    mapAD15(t,mesh)

map AD15 forces to OWENS mesh dofs

# Inputs
* `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
* `mesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `mesh::OWENSFEA.El`: see ?OWENSFEA.El

# Outputs:
* `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
* `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)

"""
function mapAD15(
    t,
    azi_j,
    mesh,
    advanceAD15;
    numAeroTS = 1,
    alwaysrecalc = true,
    verbosity = 0,
)
    Nturb = length(mesh)
    n_steps, Fx, Fy, Fz, Mx, My, Mz = advanceAD15(t, mesh, azi_j)

    # NOTE on AD15 advanceTurb values (Fx,Fy,Fz,Mx,My,Mz)
    #       - forces/moments are in hub coordinates (converted in advanceAD15)
    #       - array length is the number of OWENS mesh points
    #       - This includes the struts (and could include tower when we add that to the AD15 interface)


    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    T = Float64 # TODO: This is a placeholer, likely needs to be something else for full AD support.
    ForceValHist = [zeros(T, Int(mesh[iturb].numNodes*6), numAeroTS) for iturb = 1:Nturb]
    # DOFs are sequential through all nodes
    ForceDof=[collect(1:1:(mesh[iturb].numNodes*6)) for iturb = 1:Nturb]

    for iturb = 1:Nturb
        # Map loads over from advanceTurb
        Fx_base = zeros(T, numAeroTS)
        Fy_base = zeros(T, numAeroTS)
        Fz_base = zeros(T, numAeroTS)
        Mz_base = zeros(T, numAeroTS)
        for i = 1:mesh[iturb].numNodes
            ForceValHist[iturb][(i-1)*6+1, :] = Fx[iturb][i, 1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+2, :] = Fy[iturb][i, 1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+3, :] = Fz[iturb][i, 1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+4, :] = Mx[iturb][i, 1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+5, :] = My[iturb][i, 1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+6, :] = Mz[iturb][i, 1:numAeroTS]

            Fx_base .+= Fx[iturb][i, 1:numAeroTS]
            Fy_base .+= Fy[iturb][i, 1:numAeroTS]
            Fz_base .+= Fz[iturb][i, 1:numAeroTS]

            Mz_base .+=
                Fy[iturb][i, 1:numAeroTS] .*
                mesh[iturb].x[i]-Fx[iturb][i, 1:numAeroTS] .* mesh[iturb].y[i]
        end

        # TODO: This assumes that the 1st node is the tower connection node, which is constrained so it enables passing of values without using them on the structural side
        ForceValHist[iturb][1, 1:numAeroTS] .= Fx_base
        ForceValHist[iturb][2, 1:numAeroTS] .= Fy_base
        ForceValHist[iturb][3, 1:numAeroTS] .= Fz_base
        ForceValHist[iturb][6, 1:numAeroTS] .= Mz_base
    end

    return ForceValHist, ForceDof
end


"""
    mapACDMS(t,mesh,el)

map OWENSAero forces to OWENS mesh dofs

# Inputs
* `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
* `mesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `mesh::OWENSFEA.El`: see ?OWENSFEA.El

# Outputs:
* `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
* `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)

"""
function mapACDMS(
    t,
    azi_j,
    mesh,
    el,
    advanceTurb;
    numAeroTS = 1,
    alwaysrecalc = true,
    outputfile = nothing,
    offsetmomentarm = 0.0,
)
    aero_result = advanceTurb(t; azi = azi_j, alwaysrecalc) #add 3pi/2 to align aero with structural azimuth
    CP,
    Rp,
    Tp,
    Zp,
    alpha,
    cl,
    cd_af,
    Vloc,
    Re,
    thetavec,
    n_steps,
    Fx_base,
    Fy_base,
    Fz_base,
    Mx_base,
    My_base,
    Mz_base,
    power,
    power2,
    rev_step,
    z3Dnorm,
    delta,
    Xp,
    Yp,
    M_addedmass_Np,
    M_addedmass_Tp,
    F_addedmass_Np,
    F_addedmass_Tp = aero_result
    aero_M25 = _optional_acdms_m25(aero_result, size(Rp))

    NBlade = length(Rp[:, 1, 1])
    Nslices = length(Rp[1, :, 1])

    # Initialize bladeForces
    # TODO: This should be `eltype(Rp)` or similar but currently that return `Real` instead of a concrete type
    TT = typeof(first(Rp))
    N = zeros(TT, NBlade, numAeroTS, Nslices)
    T = zeros(TT, NBlade, numAeroTS, Nslices)
    X = zeros(TT, NBlade, numAeroTS, Nslices)
    Y = zeros(TT, NBlade, numAeroTS, Nslices)
    Z = zeros(TT, NBlade, numAeroTS, Nslices)
    M25 = zeros(TT, NBlade, numAeroTS, Nslices)

    for iTS = 1:numAeroTS
        if numAeroTS == 1
            t_idx = length(Rp[1, 1, :])
        else
            t_idx = iTS
        end
        for jbld = 1:NBlade
            for islice = 1:Nslices
                N[jbld, iTS, islice] = Rp[jbld, islice, t_idx] #Normal force on the structure is inward, OWENSAero normal is inward positive #we multiply by cos(delta) to go from force per height to force per span, and then divide by cos(delta) to go from radial to normal, so they cancel
                T[jbld, iTS, islice] = -Tp[jbld, islice, t_idx]*cos(delta[jbld, islice]) ##Tangential force on the structure is against turbine rotation, OWENSAero tangential is with rotation positive # multiply by delta to convert from force per height to force per span
                M25[jbld, iTS, islice] =
                    (
                        aero_M25 === nothing ? zero(Rp[jbld, islice, t_idx]) :
                        aero_M25[jbld, islice, t_idx]
                    ) + Rp[jbld, islice, t_idx]*offsetmomentarm
                X[jbld, iTS, islice] = Xp[jbld, islice, t_idx]#*cos(-azi_j) + Yp[jbld,islice,t_idx]*sin(-azi_j) #*cos(delta[jbld,islice]) 
                Y[jbld, iTS, islice] = Yp[jbld, islice, t_idx]#*sin(-azi_j) + Yp[jbld,islice,t_idx]*cos(-azi_j)#*cos(delta[jbld,islice]) 
                Z[jbld, iTS, islice] = Zp[jbld, islice, t_idx]#*cos(delta[jbld,islice])
            end
        end
    end

    spanLocNorm = zeros(eltype(z3Dnorm), NBlade, Nslices)

    for i = 1:NBlade
        spanLocNorm[i, :] = z3Dnorm #Note that the lookup for this is not the span position, but the vertical position
    end

    structuralSpanLocNorm = mesh.structuralSpanLocNorm # this is also just the blade z position
    structuralNodeNumbers = mesh.structuralNodeNumbers
    structuralElNumbers = mesh.structuralElNumbers

    #Initialize structuralLoad
    struct_N = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    struct_T = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    struct_M25 = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    struct_X = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    struct_Y = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    struct_Z = zeros(TT, NBlade, numAeroTS, length(structuralElNumbers[1, :]))
    if maximum(structuralSpanLocNorm)>1.0001 || minimum(structuralSpanLocNorm)<-1e-4
        @warn "extrapolating on akima spline, unexpected behavior may occur (very large numbers)."
    end
    for i = 1:NBlade
        for j = 1:numAeroTS
            # AC and DMS calculate inbetween aero slices, so we add the 0 and 1 normed values here to ensure we don't extrapolate
            struct_N[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; N[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
            struct_T[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; T[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
            struct_M25[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; M25[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
            struct_X[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; X[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
            struct_Y[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; Y[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
            struct_Z[i, j, :] = safeakima(
                [0.0; spanLocNorm[i, :]; 1.0],
                [0.0; Z[i, j, :]; 0.0],
                structuralSpanLocNorm[i, :],
            )
        end
    end

    _, numNodesPerBlade = size(structuralNodeNumbers)

    #integrate over elements

    #read element aero_data in
    numDOFPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(TT, mesh.numNodes * 6, numAeroTS)
    Fg_global = zeros(TT, mesh.numNodes * 6, numAeroTS)
    for i = 1:numAeroTS
        for j = 1:NBlade
            for k = 1:(numNodesPerBlade-1)
                #get element aero_data
                # orientation angle,xloc,sectionProps,element order]
                elNum = Int(structuralElNumbers[j, k])
                #get dof map
                node1 = Int(structuralNodeNumbers[j, k])
                node2 = Int(structuralNodeNumbers[j, k+1])
                dofList =
                    [(node1-1)*numDOFPerNode .+ (1:6) (node2-1)*numDOFPerNode .+ (1:6)]

                elementOrder = 1
                x = [mesh.x[node1], mesh.x[node2]]
                elLength = sqrt(
                    (mesh.x[node2]-mesh.x[node1])^2 +
                    (mesh.y[node2]-mesh.y[node1])^2 +
                    (mesh.z[node2]-mesh.z[node1])^2,
                )
                elHeight = abs(mesh.z[node2]-mesh.z[node1])
                xloc = [0 elLength]
                twist = el.props[elNum].twist
                sweepAngle = el.psi[elNum]
                coneAngle = el.theta[elNum]
                rollAngle = el.roll[elNum]

                extDistF2Node = [struct_T[j, i, k] struct_T[j, i, k+1]]
                extDistF3Node = [struct_N[j, i, k] struct_N[j, i, k+1]]
                extDistF4Node = [struct_M25[j, i, k] struct_M25[j, i, k+1]]

                Fe = OWENSFEA.calculateLoadVecFromDistForce(
                    elementOrder,
                    x,
                    xloc,
                    twist,
                    sweepAngle,
                    coneAngle,
                    rollAngle,
                    extDistF2Node,
                    extDistF3Node,
                    extDistF4Node,
                )
                Fe_global = [
                    struct_X[j, i, k]*elHeight,
                    struct_Y[j, i, k]*elHeight,
                    struct_Z[j, i, k]*elHeight,
                    0.0,
                    0.0,
                    0.0,
                    struct_X[j, i, k+1]*elHeight,
                    struct_Y[j, i, k+1]*elHeight,
                    struct_Z[j, i, k+1]*elHeight,
                    0.0,
                    0.0,
                    0.0,
                ]

                #assembly
                for m = 1:length(dofList)
                    Fg[dofList[m], i] = Fg[dofList[m], i]+Fe[m]
                    Fg_global[dofList[m], i] = Fg_global[dofList[m], i]+Fe_global[m]
                end

            end
        end
    end

    #reduce Fg to nonzero components
    #assumes any loaded DOF will never be identically zero throughout time
    #history
    # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
    # ForceDof = zeros(sum(Fg[:,1].!=0),1)
    ForceValHist = zeros(TT, size(Fg))
    ForceDof = zeros(Int, size(Fg, 1))
    index = 1
    for i = 1:Int(mesh.numNodes*6)
        # if !isempty(findall(x->x!=0,Fg[i,:]))

        ForceValHist[index, :] = Fg[i, :]
        ForceDof[index] = i
        index = index + 1
        # end
    end

    if outputfile!=nothing
        DelimitedFiles.open(string("$(outputfile)_fullmesh.txt"), "a") do io
            if t==0
                header1 = ["t" "azi" "nodenum" "Fx" "Fy" "Fz" "Mx" "My" "Mz"]
                header2 =
                    ["(s)" "(rad)" "(#)" "(m/s)" "(N)" "(N)" "(N)" "(N-m)" "(N-m)" "(N-m)"]

                DelimitedFiles.writedlm(io, header1, '\t')
                DelimitedFiles.writedlm(io, header2, '\t')
            end
            Fx = ForceValHist[1:6:end, end]
            Fy = ForceValHist[2:6:end, end]
            Fz = ForceValHist[3:6:end, end]
            Mx = ForceValHist[4:6:end, end]
            My = ForceValHist[5:6:end, end]
            Mz = ForceValHist[6:6:end, end]
            for inode = 1:Int(mesh.numNodes)

                data =
                    [t azi_j inode Fx[inode] Fy[inode] Fz[inode] Mx[inode] My[inode] Mz[inode]]

                DelimitedFiles.writedlm(io, data, '\t')

            end

        end

        for j = 1:NBlade
            DelimitedFiles.open(string("$(outputfile)_blade$j.txt"), "a") do io
                if t==0
                    header1 = ["t" "azi" "nodenum" "Fx" "Fy" "Fz" "Mx" "My" "Mz"]
                    header2 =
                        ["(s)" "(rad)" "(#)" "(m/s)" "(N)" "(N)" "(N)" "(N-m)" "(N-m)" "(N-m)"]

                    DelimitedFiles.writedlm(io, header1, '\t')
                    DelimitedFiles.writedlm(io, header2, '\t')
                end
                Fx = ForceValHist[1:6:end, end]
                Fy = ForceValHist[2:6:end, end]
                Fz = ForceValHist[3:6:end, end]
                Mx = ForceValHist[4:6:end, end]
                My = ForceValHist[5:6:end, end]
                Mz = ForceValHist[6:6:end, end]
                for k = 1:numNodesPerBlade
                    #get element aero_data
                    # orientation angle,xloc,sectionProps,element order]
                    elNum = Int(structuralElNumbers[j, k])
                    #get dof map
                    inode = Int(structuralNodeNumbers[j, k])

                    data =
                        [t azi_j inode Fx[inode] Fy[inode] Fz[inode] Mx[inode] My[inode] Mz[inode]]

                    DelimitedFiles.writedlm(io, data, '\t')

                end

            end
        end
    end

    # This assumes that the 1st node is the tower connection node, which is constrained so it enables passing of values without using them on the structural side
    ForceValHist[1, 1:numAeroTS] .= Fx_base
    ForceValHist[2, 1:numAeroTS] .= Fy_base
    ForceValHist[3, 1:numAeroTS] .= Fz_base
    ForceValHist[4, 1:numAeroTS] .= Mx_base
    ForceValHist[5, 1:numAeroTS] .= My_base
    ForceValHist[6, 1:numAeroTS] .= Mz_base

    # return Fexternal, Fdof
    return ForceValHist[:, 1:numAeroTS],
    ForceDof,
    Fg_global,
    ForceDof,
    ForceValHist[:, 1:numAeroTS],
    z3Dnorm
end

function _optional_acdms_m25(aero_result, reference_size)
    for index = length(aero_result):-1:29
        candidate = aero_result[index]
        if candidate isa AbstractArray &&
           ndims(candidate) == 3 &&
           size(candidate) == reference_size
            return candidate
        end
    end
    return nothing
end

# """
#     mapACDMS(t,mesh,el,loadsFn)
#
# map OWENSAero forces to OWENS mesh dofs using a file of loads TODO: merge the two functions together
#
# # Inputs
# * `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
# * `mesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
# * `mesh::OWENSFEA.El`: see ?OWENSFEA.El
# * `loadsFn::string`: path/name to loads filename, in cactus Element_Data format
#
#
# # Outputs:
# * `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
# * `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)
#
# """
# function mapCACTUSFILE_minimalio(t,mesh,el,loadsFn)
#     # TODO: if this function is used, either export the required data from OWENSAero to get rid of these globals, or put the function back into OWENSAero scope
#     global turbslices
#     global envslices
#     NBlade = turbslices[1].B
#
#     aero_data = DelimitedFiles.readdlm(loadsFn,',',skipstart = 1)
#
#     #define these from params file
#     rho = envslices[1].rho
#
#     RefR = turbslices[1].R
#     NElem = length(turbslices)
#     #TODO: get from OWENSAero structs
#     # chord = RefR.*[1.11093e-01,1.11093e-01,8.20390e-02,6.35940e-02,5.68145e-02,5.55467e-02,5.88008e-02,6.82860e-02,9.11773e-02,1.11093e-01,1.11093e-01]
#     # chord = (chord[1:end-1]+chord[2:end])./2
#     chord = [turb.chord for turb in turbslices]
#     V = envslices[1].V_x[1] #m/s #TODO: get nominal vinf
#     global z3Dnorm
#     # z3Dnorm = 1/2.44680.*[1.22340e-01,3.67020e-01,6.11700e-01,8.56380e-01,1.10106e+00,1.34574e+00,1.59042e+00,1.83510e+00,2.07978e+00,2.32446e+00]
#
#     normTime = aero_data[:,1]
#
#     numAeroEl = 0
#     for i=1:NBlade
#         numAeroEl = numAeroEl + NElem
#     end
#
#     len,_ = size(aero_data)
#
#     numAeroTS = Int(len/numAeroEl)
#
#     time = normTime[1:Int(numAeroEl):end,1].*RefR[1]./V[1]
#
#     urel = aero_data[:,15]
#     uloc = urel.*V
#
#     cn = aero_data[:,24]
#     ct = aero_data[:,25]
#     cm25 = aero_data[:,22]
#
#     NperSpan = zeros(len)
#     TperSpan = zeros(len)
#     M25perSpan = zeros(len)
#
#     for i=1:len
#         NperSpan[i] =  cn[i]  * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
#         TperSpan[i] =  ct[i]  * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
#         M25perSpan[i] = cm25[i] * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)*blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR
#     end
#
#     # Initialize bladeForces
#     N = zeros(NBlade,numAeroTS,NElem)
#     T = zeros(NBlade,numAeroTS,NElem)
#     M25 = zeros(NBlade,numAeroTS,NElem)
#
#     index = 1
#     for i=1:numAeroTS
#         for j=1:NBlade
#             for k=1:NElem
#                 N[j,i,k] = NperSpan[index]
#                 T[j,i,k] = TperSpan[index]
#                 M25[j,i,k] = M25perSpan[index]
#                 index = index + 1
#             end
#         end
#     end
#
#     #Apply chord since it was pulled out above
#     for i=1:numAeroTS
#         for j=1:NBlade
#             N[j,i,:] = N[j,i,:].*chord
#             T[j,i,:] = T[j,i,:].*chord
#             M25[j,i,:] = M25[j,i,:].*chord
#         end
#     end
#
#     spanLocNorm = zeros(NBlade,NElem)
#
#     for i=1:NBlade
#         spanLocNorm[i,:] = z3Dnorm
#     end
#
#     structuralSpanLocNorm = mesh.structuralSpanLocNorm
#     structuralNodeNumbers = mesh.structuralNodeNumbers
#     structuralElNumbers = mesh.structuralElNumbers
#
#     #Initialize structuralLoad
#     struct_N = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#     struct_T = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#     struct_M25 = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#
#     for i=1:NBlade
#         for j=1:numAeroTS
#             struct_N[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],N[i,j,:],structuralSpanLocNorm[i,:])
#             struct_T[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],T[i,j,:],structuralSpanLocNorm[i,:])
#             struct_M25[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],M25[i,j,:],structuralSpanLocNorm[i,:])
#         end
#     end
#
#     _,numNodesPerBlade = size(structuralNodeNumbers)
#
#     #integrate over elements
#
#     #read element aero_data in
#     numDOFPerNode = 6
#     #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
#     Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6),numAeroTS)
#     for i=1:numAeroTS
#         for j = 1:NBlade
#             for k = 1:numNodesPerBlade-1
#                 #get element aero_data
#                 # orientation angle,xloc,sectionProps,element order]
#                 elNum = Int(structuralElNumbers[j,k])
#                 #get dof map
#                 node1 = Int(structuralNodeNumbers[j,k])
#                 node2 = Int(structuralNodeNumbers[j,k+1])
#                 dofList = [(node1-1)*numDOFPerNode.+(1:6) (node2-1)*numDOFPerNode.+(1:6)]
#
#                 elementOrder = 1
#                 x = [mesh.x[node1], mesh.x[node2]]
#                 elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
#                 xloc = [0 elLength]
#                 twist = el.props[elNum].twist
#                 sweepAngle = el.psi[elNum]
#                 coneAngle = el.theta[elNum]
#                 rollAngle = el.roll[elNum]
#
#                 extDistF2Node =  [struct_T[j,i,k]    struct_T[j,i,k+1]]
#                 extDistF3Node = -[struct_N[j,i,k]    struct_N[j,i,k+1]]
#                 extDistF4Node = -[struct_M25[j,i,k]  struct_M25[j,i,k+1]]
#
#                 Fe = OWENSFEA.calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)
#
#                 #asssembly
#                 for m = 1:length(dofList)
#                     Fg[dofList[m],i] =  Fg[dofList[m],i]+Fe[m]
#                 end
#
#             end
#         end
#     end
#
#     #reduce Fg to nonzero components
#     #assumes any loaded DOF will never be identically zero throughout time
#     #history
#     # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
#     # ForceDof = zeros(sum(Fg[:,1].!=0),1)
#     ForceValHist = zeros(length(Fg[:,1]),length(Fg[1,:]))
#     ForceDof = zeros(length(Fg[:,1]),1)
#     index = 1
#     for i=1:Int(maximum(maximum(structuralNodeNumbers))*6)
#         # if !isempty(findall(x->x!=0,Fg[i,:]))
#
#             ForceValHist[index,:] = Fg[i,:]
#             ForceDof[index] = i
#             index = index + 1
#         # end
#     end
#
#     #TODO: wrap the function at this level for time so you don't read in the file each time
#     Fexternal = zeros(length(ForceDof))
#     for i = 1:length(ForceDof)
#         Fexternal[i] = FLOWMath.linear(time,ForceValHist[i,:],t)
#     end
#
#     return Fexternal, ForceDof
# end
