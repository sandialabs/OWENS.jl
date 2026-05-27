"""
    GXBeamSectionalProperties

Container for a GXBeam cross-section mesh, its reusable section cache, and the
section matrices computed from that mesh.
"""
struct GXBeamSectionalProperties
    nodes::Any
    elements::Any
    cache::Any
    compliance::Any
    stiffness::Any
    mass::Any
    shear_center::Any
    tension_center::Any
    mass_center::Any
end

"""
    GXBeamSectionalRecovery

Recovered sectional strain/stress arrays for one structural element at one
time step. Arrays follow GXBeam's section-mesh element ordering.
"""
struct GXBeamSectionalRecovery
    beam_strain::Any
    beam_stress::Any
    ply_strain::Any
    ply_stress::Any
end

_section_source_symbol(source) = source isa Symbol ? source : Symbol(source)

function _precomp_section_real_type(pc_input)
    return promote_type(
        typeof(pc_input.chord),
        typeof(pc_input.le_loc),
        eltype(pc_input.xnode),
        eltype(pc_input.ynode),
        eltype(pc_input.e1),
        eltype(pc_input.e2),
        eltype(pc_input.g12),
        eltype(pc_input.anu12),
        eltype(pc_input.density),
        eltype(pc_input.n_pliesU),
        eltype(pc_input.t_lamU),
        eltype(pc_input.tht_lamU),
        eltype(pc_input.n_pliesL),
        eltype(pc_input.t_lamL),
        eltype(pc_input.tht_lamL),
        eltype(pc_input.n_pliesW),
        eltype(pc_input.t_lamW),
        eltype(pc_input.tht_lamW),
    )
end

function _gxbeam_materials_from_precomp(pc_input)
    T = _precomp_section_real_type(pc_input)
    materials = GXBeam.Material{T}[]

    for i in eachindex(pc_input.e1)
        E1 = T(pc_input.e1[i])
        E2 = T(pc_input.e2[i])
        G12 = T(pc_input.g12[i])
        nu12 = T(pc_input.anu12[i])

        # OWENSPreComp only carries the plane-stress lamina data.  GXBeam's
        # section solver needs the 3-D orthotropic form, so use the available
        # transverse values consistently for the missing through-thickness terms.
        push!(
            materials,
            GXBeam.Material{T}(
                E1,
                E2,
                E2,
                G12,
                G12,
                G12,
                nu12,
                nu12,
                nu12,
                T(pc_input.density[i]),
            ),
        )
    end

    return materials
end

function _gxbeam_layer_sets(materials, n_lamina, n_plies, t_lam, tht_lam, mat_lam, T)
    layer_sets = Vector{Vector{GXBeam.Layer{T}}}(undef, length(n_lamina))
    angle_scale = one(T) * (pi / 180)
    idx = 1

    for isegment in eachindex(n_lamina)
        layers = GXBeam.Layer{T}[]
        for _ = 1:n_lamina[isegment]
            thickness = T(n_plies[idx]) * T(t_lam[idx])
            if thickness > zero(T)
                push!(
                    layers,
                    GXBeam.Layer(
                        materials[mat_lam[idx]],
                        thickness,
                        T(tht_lam[idx]) * angle_scale,
                    ),
                )
            end
            idx += 1
        end
        isempty(layers) && throw(
            ArgumentError(
                "GXBeam sectional meshing requires at least one positive-thickness layer per PreComp sector",
            ),
        )
        layer_sets[isegment] = layers
    end

    return layer_sets
end

function _unique_section_breakpoints(values)
    T = eltype(values)
    raw = sort!(clamp.(collect(values), zero(T), one(T)))
    out = T[]
    for value in raw
        if isempty(out) || abs(value - out[end]) > T(1e-8)
            push!(out, value)
        end
    end
    isempty(out) && append!(out, T[zero(T), one(T)])
    out[1] = zero(T)
    out[end] = one(T)
    return out
end

function _sector_layers_on_breakpoints(source_breaks, source_layers, target_breaks)
    layers = Vector{eltype(source_layers)}(undef, length(target_breaks) - 1)

    for i in eachindex(layers)
        midpoint = (target_breaks[i] + target_breaks[i+1]) / 2
        idx = searchsortedlast(source_breaks, midpoint)
        idx = clamp(idx, 1, length(source_layers))
        layers[i] = source_layers[idx]
    end

    return layers
end

function _precomp_airfoil_to_gxbeam(pc_input)
    xnode = collect(pc_input.xnode)
    ynode = collect(pc_input.ynode)
    T = promote_type(eltype(xnode), eltype(ynode))

    if isapprox(xnode[1], one(T); atol = T(1e-8)) &&
       isapprox(xnode[end], one(T); atol = T(1e-8))
        xaf = copy(xnode)
        yaf = copy(ynode)
        if length(yaf) > 2 && yaf[3] < zero(T)
            reverse!(xaf)
            reverse!(yaf)
        end
        return xaf, yaf
    end

    te_idx = argmax(xnode)
    te_idx == 1 && throw(
        ArgumentError(
            "PreComp airfoil input must include an upper surface before the trailing edge",
        ),
    )

    upper_x = xnode[1:te_idx]
    upper_y = ynode[1:te_idx]
    lower_x = xnode[te_idx:end]
    lower_y = ynode[te_idx:end]

    if mean(upper_y) < mean(lower_y)
        upper_x, lower_x = lower_x, upper_x
        upper_y, lower_y = lower_y, upper_y
    end

    xaf = [reverse(upper_x); reverse(lower_x[1:(end-1)])]
    yaf = [reverse(upper_y); reverse(lower_y[1:(end-1)])]
    xaf[1] = one(T)
    xaf[end] = one(T)

    return xaf, yaf
end

function _asymmetric_gxbeam_airfoil_mesh(
    xaf,
    yaf,
    chord,
    xbreak,
    webloc,
    segments_upper,
    segments_lower,
    webs;
    ds = nothing,
    dt = nothing,
    ns = nothing,
    nt = nothing,
    wns = 4,
    wnt = nothing,
)
    skin_segments, webs =
        GXBeam.preprocess_layers([segments_upper; segments_lower], webs, dt, nt, wnt)
    n_upper = length(segments_upper)
    segments_upper = skin_segments[1:n_upper]
    segments_lower = skin_segments[(n_upper+1):end]

    xu, yu, xl, yl = GXBeam.parseairfoil(xaf, yaf, xbreak, ds, ns)

    xu *= chord
    yu *= chord
    xl *= chord
    yl *= chord
    xbreak *= chord

    txu, tyu = GXBeam.tangential(xu, yu; upper = true)
    txl, tyl = GXBeam.tangential(xl, yl; upper = false)

    xiu, yiu = GXBeam.find_inner_surface(xu, yu, txu, tyu, segments_upper, xbreak)
    xil, yil = GXBeam.find_inner_surface(xl, yl, txl, tyl, segments_lower, xbreak)

    nw = length(webs)
    idx_webu = Vector{Int}(undef, nw)
    idx_webl = Vector{Int}(undef, nw)
    nx_web = Vector{Int}(undef, nw)

    for i = 1:nw
        idx_webu[i], xiu, yiu, xu, yu, txu, tyu =
            GXBeam.web_intersections(xiu, yiu, xu, yu, txu, tyu, chord, webloc[i], webs[i])
        idx_webl[i], xil, yil, xl, yl, txl, tyl =
            GXBeam.web_intersections(xil, yil, xl, yl, txl, tyl, chord, webloc[i], webs[i])
        nx_web[i] = length(webs[i]) + 1
    end

    x_te, y_te, xu, yu, xl, yl =
        GXBeam.te_inner_intersection(xiu, yiu, xil, yil, xu, yu, xl, yl)
    nodesu, elementsu =
        GXBeam.nodes_half(xu, yu, txu, tyu, xbreak, segments_upper, chord, x_te, y_te)
    nodesl, elementsl =
        GXBeam.nodes_half(xl, yl, txl, tyl, xbreak, segments_lower, chord, x_te, y_te)

    nlayer = length(segments_upper[1])
    nodes, elements =
        GXBeam.combine_halfs(nodesu, elementsu, nodesl, elementsl, nlayer, x_te)

    if nw > 0
        nodes, elements = GXBeam.addwebs(
            idx_webu,
            idx_webl,
            nx_web,
            nodes,
            elements,
            webs,
            length(nodesu),
            nlayer,
            wns,
        )
    end

    return nodes, elements
end

"""
    gxbeam_sectional_mesh_from_precomp(precompinput; kwargs...)

Build a GXBeam shell mesh from one `OWENSPreComp.Input`, preserving distinct
upper-surface, lower-surface, and web layups.
"""
function gxbeam_sectional_mesh_from_precomp(
    pc_input;
    ds = nothing,
    dt = nothing,
    ns = nothing,
    nt = nothing,
    wns = 4,
    wnt = nothing,
)
    T = _precomp_section_real_type(pc_input)
    materials = _gxbeam_materials_from_precomp(pc_input)

    upper_layers = _gxbeam_layer_sets(
        materials,
        pc_input.n_laminaU,
        pc_input.n_pliesU,
        pc_input.t_lamU,
        pc_input.tht_lamU,
        pc_input.mat_lamU,
        T,
    )
    lower_layers = _gxbeam_layer_sets(
        materials,
        pc_input.n_laminaL,
        pc_input.n_pliesL,
        pc_input.t_lamL,
        pc_input.tht_lamL,
        pc_input.mat_lamL,
        T,
    )
    web_layers = _gxbeam_layer_sets(
        materials,
        pc_input.n_laminaW,
        pc_input.n_pliesW,
        pc_input.t_lamW,
        pc_input.tht_lamW,
        pc_input.mat_lamW,
        T,
    )

    web_pairs = [
        (T(pc_input.loc_web[iweb]), web_layers[iweb]) for
        iweb in eachindex(web_layers) if
        T(1e-8) < T(pc_input.loc_web[iweb]) < one(T) - T(1e-8)
    ]
    sort!(web_pairs; by = first)
    web_breaks = first.(web_pairs)
    web_layers = last.(web_pairs)
    xbreak = _unique_section_breakpoints(
        [T.(pc_input.xsec_nodeU); T.(pc_input.xsec_nodeL); web_breaks],
    )

    upper_layers =
        _sector_layers_on_breakpoints(T.(pc_input.xsec_nodeU), upper_layers, xbreak)
    lower_layers =
        _sector_layers_on_breakpoints(T.(pc_input.xsec_nodeL), lower_layers, xbreak)

    xaf, yaf = _precomp_airfoil_to_gxbeam(pc_input)
    nodes, elements = _asymmetric_gxbeam_airfoil_mesh(
        T.(xaf),
        T.(yaf),
        T(pc_input.chord),
        xbreak,
        web_breaks,
        upper_layers,
        lower_layers,
        web_layers;
        ds,
        dt,
        ns,
        nt,
        wns,
        wnt,
    )

    reference_x = T(pc_input.le_loc) * T(pc_input.chord)
    nodes = [GXBeam.Node(node.x - reference_x, node.y) for node in nodes]

    return nodes, elements
end

"""
    gxbeam_sectional_properties_from_precomp(precompinput; kwargs...)

Build a GXBeam sectional mesh from an `OWENSPreComp.Input` and compute the
compliance, stiffness, and mass matrices needed by GXBeam.
"""
function gxbeam_sectional_properties_from_precomp(
    pc_input;
    shear_center = false,
    ds = nothing,
    dt = nothing,
    ns = nothing,
    nt = nothing,
    wns = 4,
    wnt = nothing,
)
    nodes, elements = gxbeam_sectional_mesh_from_precomp(pc_input; ds, dt, ns, nt, wns, wnt)
    cache = GXBeam.initialize_cache(nodes, elements)
    compliance, shear_center_xy, tension_center_xy =
        GXBeam.compliance_matrix(nodes, elements; cache, gxbeam_order = true, shear_center)
    mass, mass_center_xy = GXBeam.mass_matrix(nodes, elements)

    stiffness = LinearAlgebra.inv(Matrix(compliance))
    stiffness = Matrix(LinearAlgebra.Symmetric((stiffness + stiffness') / 2))

    return GXBeamSectionalProperties(
        nodes,
        elements,
        cache,
        Matrix(compliance),
        stiffness,
        Matrix(mass),
        shear_center_xy,
        tension_center_xy,
        mass_center_xy,
    )
end

function _interpolate_gxbeam_section_matrices(
    orig_unit_span,
    used_unit_span,
    section_properties,
    added_M22,
    added_M33,
)
    element_span = used_unit_span[1:(end-1)]
    n_element = length(element_span)
    stiff = [zeros(eltype(section_properties[1].stiffness), 6, 6) for _ = 1:n_element]
    mass = [zeros(eltype(section_properties[1].mass), 6, 6) for _ = 1:n_element]

    for irow = 1:6, icol = 1:6
        stiffness_values = [section.stiffness[irow, icol] for section in section_properties]
        mass_values = [section.mass[irow, icol] for section in section_properties]
        stiffness_used = safeakima(orig_unit_span, stiffness_values, element_span)
        mass_used = safeakima(orig_unit_span, mass_values, element_span)

        for ielem = 1:n_element
            stiff[ielem][irow, icol] = stiffness_used[ielem]
            mass[ielem][irow, icol] = mass_used[ielem]
        end
    end

    for ielem = 1:n_element
        stiff[ielem] = Matrix(LinearAlgebra.Symmetric((stiff[ielem] + stiff[ielem]') / 2))
        mass[ielem][2, 2] += added_M22[ielem]
        mass[ielem][3, 3] += added_M33[ielem]
        mass[ielem] = Matrix(LinearAlgebra.Symmetric((mass[ielem] + mass[ielem]') / 2))
    end

    return stiff, mass
end

function _nearest_gxbeam_sections(orig_unit_span, used_unit_span, section_properties)
    element_span = used_unit_span[1:(end-1)]
    return map(element_span) do span
        idx = argmin(abs.(orig_unit_span .- span))
        section_properties[idx]
    end
end

"""
    recover_gxbeam_sectional_load(section, generalized_load)

Recover beam and ply strain/stress in a GXBeam sectional mesh from a 6-vector
of beam forces and moments in GXBeam order.
"""
function recover_gxbeam_sectional_load(section::GXBeamSectionalProperties, generalized_load)
    length(generalized_load) == 6 ||
        throw(ArgumentError("generalized_load must contain 6 force/moment components"))

    return GXBeam.strain_recovery(
        generalized_load[1:3],
        generalized_load[4:6],
        section.nodes,
        section.elements,
        section.cache;
        gxbeam_order = true,
    )
end

"""
    recover_gxbeam_sectional_strain(section, generalized_strain[, stiffness])

Map a 6-vector of beam strains/curvatures in GXBeam order through the section
stiffness and recover element-level beam/ply strains and stresses.
"""
function recover_gxbeam_sectional_strain(
    section::GXBeamSectionalProperties,
    generalized_strain,
    stiffness = section.stiffness,
)
    length(generalized_strain) == 6 || throw(
        ArgumentError("generalized_strain must contain 6 strain/curvature components"),
    )

    generalized_load = stiffness * generalized_strain
    return recover_gxbeam_sectional_load(section, generalized_load)
end

"""
    recover_gxbeam_component_sectional_strain(component)

Recover GXBeam section-mesh beam/ply strains and stresses for each structural
element and time step in a component. The component must have element-aligned
`gxBeamSectionalProperties` plus `e_x`, `e_y`, `e_z`, `k_x`, `k_y`, and `k_z`
histories with dimensions `(nelem, ntime)`.
"""
function recover_gxbeam_component_sectional_strain(component)
    sections = component.gxBeamSectionalProperties
    isnothing(sections) && return nothing

    histories = (
        component.e_x,
        component.e_y,
        component.e_z,
        component.k_x,
        component.k_y,
        component.k_z,
    )
    if any(isnothing, histories)
        throw(
            ArgumentError(
                "GXBeam sectional recovery requires component beam strain and curvature histories",
            ),
        )
    end

    reference_size = size(component.e_x)
    length(reference_size) == 2 || throw(
        ArgumentError("component strain histories must have dimensions (nelem, ntime)"),
    )
    all(size(history) == reference_size for history in histories) || throw(
        DimensionMismatch(
            "component strain and curvature histories must have matching sizes",
        ),
    )

    nelem, ntime = reference_size
    length(sections) == nelem || throw(
        DimensionMismatch(
            "component has $(length(sections)) GXBeam sectional meshes but $nelem strain-history elements",
        ),
    )

    recovery = Matrix{GXBeamSectionalRecovery}(undef, nelem, ntime)
    for it = 1:ntime, ielem = 1:nelem
        generalized_strain = [
            component.e_x[ielem, it],
            component.e_y[ielem, it],
            component.e_z[ielem, it],
            component.k_x[ielem, it],
            component.k_y[ielem, it],
            component.k_z[ielem, it],
        ]
        beam_strain, beam_stress, ply_strain, ply_stress =
            recover_gxbeam_sectional_strain(sections[ielem], generalized_strain)
        recovery[ielem, it] =
            GXBeamSectionalRecovery(beam_strain, beam_stress, ply_strain, ply_stress)
    end

    return recovery
end

function populate_gxbeam_sectional_recovery!(component)
    component.gxBeamSectionalRecovery = recover_gxbeam_component_sectional_strain(component)
    return component
end
