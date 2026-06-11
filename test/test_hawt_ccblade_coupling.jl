using OWENS
using Test
using LinearAlgebra

function _hawt_test_mesh()
    mesh, _, _, _, _ = OWENS.create_hawt_mesh(;
        hub_depth = 10.0,
        tip_precone = 0.0,
        R = 5.0,
        AD15hubR = 0.0,
        nblade = 3,
        ntelem = 1,
        nbelem = 4,
        bshapex = LinRange(0.0, 1.0, 5),
        bshapez = zeros(5),
    )
    return mesh
end

function _shaft_x_hawt_test_mesh()
    mesh = deepcopy(_hawt_test_mesh())
    old_x = copy(mesh.x)
    old_y = copy(mesh.y)
    old_z = copy(mesh.z)
    mesh.x .= 0.0
    mesh.y .= old_x
    mesh.z .= old_z .+ old_y
    return mesh
end

function _shaft_y_hawt_test_mesh()
    mesh = deepcopy(_hawt_test_mesh())
    old_x = copy(mesh.x)
    old_y = copy(mesh.y)
    old_z = copy(mesh.z)
    mesh.x .= old_x
    mesh.y .= 0.0
    mesh.z .= old_z .+ old_y
    return mesh
end

function _shaft_x_hawt_test_structure(; analysisType = "S")
    mesh, _, joints, _, _ = OWENS.create_hawt_mesh(;
        hub_depth = 10.0,
        tip_precone = 0.0,
        R = 5.0,
        AD15hubR = 0.0,
        nblade = 3,
        ntelem = 2,
        nbelem = 4,
        bshapex = LinRange(0.0, 1.0, 5),
        bshapez = zeros(5),
    )

    old_x = copy(mesh.x)
    old_y = copy(mesh.y)
    old_z = copy(mesh.z)
    mesh.x .= 0.0
    mesh.y .= old_x
    mesh.z .= old_z .+ old_y
    ort = OWENS.calculateElementOrientation(mesh)

    area = 0.08
    chordwise_depth = 0.4
    flapwise_depth = 0.2
    elastic_modulus = 5.0e10
    shear_modulus = elastic_modulus / (2 * (1 + 0.3))
    rho = 1600.0
    iyy = chordwise_depth * flapwise_depth^3 / 12
    izz = flapwise_depth * chordwise_depth^3 / 12
    polar = iyy + izz
    pair(value) = fill(value, 2)

    section_props = Array{OWENS.OWENSFEA.SectionPropsArray,1}(undef, mesh.numEl)
    for iel = 1:mesh.numEl
        section_props[iel] = OWENS.OWENSFEA.SectionPropsArray(
            pair(0.0),
            pair(0.0),
            fill(rho * area, 2),
            fill(elastic_modulus * iyy, 2),
            fill(elastic_modulus * izz, 2),
            fill(shear_modulus * polar, 2),
            fill(elastic_modulus * area, 2),
            fill(rho * iyy, 2),
            fill(rho * izz, 2),
            fill(rho * polar, 2),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            pair(0.0),
            nothing,
            nothing,
            pair(0.0),
            pair(0.0),
            fill(shear_modulus * area, 2),
            fill(shear_modulus * area, 2),
        )
    end

    el = OWENS.OWENSFEA.El(
        section_props,
        ort.Length,
        ort.Psi_d,
        ort.Theta_d,
        ort.Twist_d,
        ones(mesh.numEl),
    )
    root_fixed = [
        1 1 0
        1 2 0
        1 3 0
        1 4 0
        1 5 0
        1 6 0
    ]
    feamodel = OWENS.OWENSFEA.FEAModel(;
        analysisType,
        dataOutputFilename = "none",
        joint = joints,
        pBC = root_fixed,
        numNodes = mesh.numNodes,
        nlOn = false,
        gravityOn = false,
        iterationType = "LINEAR",
        maxNumLoadSteps = 3,
    )
    return mesh, el, feamodel
end

_hawt_linear_airfoil(alpha, reynolds_number, mach) =
    (2 * pi * alpha, 0.01 + 0.02 * alpha^2, -0.05)

@testset "HAWT mesh topology uses integer indices and correct element types" begin
    mesh, _, joints, blade_node_ranges, blade_element_ranges = OWENS.create_hawt_mesh(;
        hub_depth = 10.0,
        tip_precone = 0.0,
        R = 5.0,
        AD15hubR = 0.0,
        nblade = 2,
        ntelem = 2,
        nbelem = 3,
        bshapex = LinRange(0.0, 1.0, 4),
        bshapez = zeros(4),
    )

    @test mesh.numNodes == 11
    @test mesh.numEl == 8
    @test mesh.nodeNum == collect(1:mesh.numNodes)
    @test mesh.elNum == collect(1:mesh.numEl)
    @test eltype(mesh.x) == Float64
    @test eltype(mesh.y) == Float64
    @test eltype(mesh.z) == Float64
    @test mesh.meshSeg == [2, 3, 3]
    @test mesh.type == [1, 1, 0, 0, 0, 0, 0, 0]
    @test mesh.conn == [1 2; 2 3; 4 5; 5 6; 6 7; 8 9; 9 10; 10 11]
    @test Int.(joints[:, 2]) == [3, 3]
    @test Int.(joints[:, 3]) == [4, 8]
    @test blade_node_ranges == [7 4; 11 8]
    @test blade_element_ranges == [5 3; 8 6]
end

@testset "HAWT mesh hub-radius adjustment is silent" begin
    result = nothing
    captured_stdout = mktemp() do path, io
        result = redirect_stdout(io) do
            OWENS.create_hawt_mesh(;
                hub_depth = 10.0,
                tip_precone = 0.0,
                R = 5.0,
                AD15hubR = 1.0,
                nblade = 2,
                ntelem = 2,
                nbelem = 4,
                bshapex = LinRange(0.0, 1.0, 5),
                bshapez = zeros(5),
            )
        end
        flush(io)
        read(path, String)
    end
    mesh, _, _, blade_node_ranges, _ = result

    @test captured_stdout == ""
    @test mesh.numNodes == 13
    @test blade_node_ranges == [8 5; 13 10]
end

@testset "HAWT mesh input validation" begin
    @test_throws ArgumentError OWENS.create_hawt_mesh(; nblade = 0)
    @test_throws ArgumentError OWENS.create_hawt_mesh(; ntelem = 0)
    @test_throws ArgumentError OWENS.create_hawt_mesh(; nbelem = 0)
    @test_throws DimensionMismatch OWENS.create_hawt_mesh(;
        nbelem = 3,
        bshapex = LinRange(0.0, 1.0, 3),
        bshapez = zeros(4),
    )
    @test_throws DimensionMismatch OWENS.create_hawt_mesh(;
        nbelem = 3,
        bshapex = LinRange(0.0, 1.0, 4),
        bshapez = zeros(3),
    )
    @test_throws ArgumentError OWENS.create_hawt_mesh(;
        nbelem = 3,
        bshapex = zeros(4),
        bshapez = zeros(4),
    )
end

@testset "HAWT CCBlade load mapping defaults to shaft-x" begin
    mesh = _shaft_x_hawt_test_mesh()
    radial_positions = collect(LinRange(0.0, 5.0, 5))
    normal_loads = fill(2.0, length(radial_positions))
    tangential_loads = fill(3.0, length(radial_positions))

    force_values, _ = OWENS.mapHAWTCCBladeLoads(
        mesh,
        radial_positions,
        normal_loads,
        tangential_loads;
        tip_radius = 5.0,
    )
    resultants = OWENS.hawtNodalLoadResultants(mesh, force_values)

    @test resultants.force[1] ≈ 3 * 2.0 * 5.0 atol = 1e-12
    @test resultants.force[2] ≈ 0.0 atol = 1e-12
    @test resultants.force[3] ≈ 0.0 atol = 1e-12
    @test resultants.moment[1] ≈ 3 * 3.0 * 5.0^2 / 2 atol = 1e-12
    @test resultants.moment[2] ≈ 0.0 atol = 1e-12
    @test resultants.moment[3] ≈ 0.0 atol = 1e-12
end

@testset "HAWT CCBlade load mapping conserves resultants" begin
    mesh = _hawt_test_mesh()
    radial_positions = collect(LinRange(0.0, 5.0, 5))
    normal_loads = fill(2.0, length(radial_positions))
    tangential_loads = fill(3.0, length(radial_positions))

    force_values, force_dofs = OWENS.mapHAWTCCBladeLoads(
        mesh,
        radial_positions,
        normal_loads,
        tangential_loads;
        hub_radius = 0.0,
        tip_radius = 5.0,
        rotor_axis = :z,
    )
    resultants = OWENS.hawtNodalLoadResultants(mesh, force_values; rotor_axis = :z)

    @test force_dofs == collect(1:(mesh.numNodes*6))
    @test length(force_values) == mesh.numNodes * 6
    @test resultants.force[1] ≈ 0.0 atol = 1e-12
    @test resultants.force[2] ≈ 0.0 atol = 1e-12
    @test resultants.force[3] ≈ 3 * 2.0 * 5.0 atol = 1e-12
    @test resultants.moment[1] ≈ 0.0 atol = 1e-12
    @test resultants.moment[2] ≈ 0.0 atol = 1e-12
    @test resultants.moment[3] ≈ 3 * 3.0 * 5.0^2 / 2 atol = 1e-12
end

@testset "HAWT CCBlade load mapping signs are consistent across rotor axes" begin
    radial_positions = collect(LinRange(0.0, 5.0, 5))
    normal_loads = fill(2.0, length(radial_positions))
    positive_tangential_loads = fill(3.0, length(radial_positions))
    negative_tangential_loads = fill(-3.0, length(radial_positions))
    expected_force = 3 * 2.0 * 5.0
    expected_moment = 3 * 3.0 * 5.0^2 / 2
    axis_cases = (
        (:x, _shaft_x_hawt_test_mesh()),
        (:y, _shaft_y_hawt_test_mesh()),
        (:z, _hawt_test_mesh()),
    )

    for (rotor_axis, mesh) in axis_cases
        axis_index = findfirst(==(rotor_axis), (:x, :y, :z))
        transverse_indices = setdiff(1:3, [axis_index])

        force_values, _ = OWENS.mapHAWTCCBladeLoads(
            mesh,
            radial_positions,
            normal_loads,
            positive_tangential_loads;
            tip_radius = 5.0,
            rotor_axis,
        )
        resultants = OWENS.hawtNodalLoadResultants(mesh, force_values; rotor_axis)

        @test resultants.force[axis_index] ≈ expected_force atol = 1e-12
        @test resultants.moment[axis_index] ≈ expected_moment atol = 1e-12
        @test all(isapprox.(resultants.force[transverse_indices], 0.0; atol = 1e-12))
        @test all(isapprox.(resultants.moment[transverse_indices], 0.0; atol = 1e-12))

        reverse_force_values, _ = OWENS.mapHAWTCCBladeLoads(
            mesh,
            radial_positions,
            normal_loads,
            negative_tangential_loads;
            tip_radius = 5.0,
            rotor_axis,
        )
        reverse_resultants =
            OWENS.hawtNodalLoadResultants(mesh, reverse_force_values; rotor_axis)

        @test reverse_resultants.force[axis_index] ≈ expected_force atol = 1e-12
        @test reverse_resultants.moment[axis_index] ≈ -expected_moment atol = 1e-12
        @test all(isapprox.(reverse_resultants.force[transverse_indices], 0.0; atol = 1e-12))
        @test all(isapprox.(reverse_resultants.moment[transverse_indices], 0.0; atol = 1e-12))
    end
end

@testset "HAWT CCBlade load mapping supports per-blade loads" begin
    mesh = _hawt_test_mesh()
    radial_positions = collect(LinRange(0.0, 5.0, 5))
    normal_loads = [
        fill(1.0, length(radial_positions))';
        fill(2.0, length(radial_positions))';
        fill(3.0, length(radial_positions))'
    ]
    tangential_loads = zeros(size(normal_loads))

    force_values, _ = OWENS.mapHAWTCCBladeLoads(
        mesh,
        radial_positions,
        normal_loads,
        tangential_loads;
        tip_radius = 5.0,
        rotor_axis = :z,
    )
    resultants = OWENS.hawtNodalLoadResultants(mesh, force_values; rotor_axis = :z)

    @test resultants.force[3] ≈ (1.0 + 2.0 + 3.0) * 5.0 atol = 1e-12
    @test resultants.moment[3] ≈ 0.0 atol = 1e-12
end

@testset "HAWT structural stations include deflection" begin
    mesh = _hawt_test_mesh()
    stations = OWENS.hawtStructuralRadialStations(mesh; rotor_axis = :z)
    displacements = zeros(mesh.numNodes * 6)
    tip_node = Int(mesh.structuralNodeNumbers[1, end])
    displacements[6 * (tip_node - 1) + 1] = 4.0
    displaced_stations =
        OWENS.hawtStructuralRadialStations(mesh; displacements, rotor_axis = :z)

    @test stations[1, end] ≈ 5.0 atol = 1e-12
    @test displaced_stations[1, end] ≈ hypot(5.0, 4.0) atol = 1e-12
end

@testset "HAWT CCBlade load mapping validates inputs" begin
    mesh = _hawt_test_mesh()
    @test_throws ArgumentError OWENS.mapHAWTCCBladeLoads(
        mesh,
        [0.0, 1.0],
        [1.0],
        [1.0, 1.0],
    )
    @test_throws ArgumentError OWENS.mapHAWTCCBladeLoads(
        mesh,
        [0.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0];
        hub_radius = 0.5,
        rotor_axis = :z,
    )
    @test_throws ArgumentError OWENS.mapHAWTCCBladeLoads(
        mesh,
        [1.0, 0.0],
        [1.0, 1.0],
        [1.0, 1.0],
    )
    @test_throws ArgumentError OWENS.hawtStructuralRadialStations(
        mesh;
        displacements = zeros(mesh.numNodes * 6 - 1),
    )
end

@testset "Native HAWT CCBlade Oye step maps unsteady loads" begin
    mesh = _shaft_x_hawt_test_mesh()
    radial_positions = [1.0, 2.5, 4.5]
    chord = [0.5, 0.42, 0.32]
    twist = deg2rad.([10.0, 6.0, 2.0])
    state = OWENS.initializeHAWTOyeState(radial_positions)

    step = OWENS.runHAWTCCBladeOyeStep(
        mesh,
        radial_positions,
        chord,
        twist,
        _hawt_linear_airfoil,
        state,
        0.2,
        2.0,
        8.0,
        1.225;
        hub_radius = 0.5,
        tip_radius = 5.0,
        pitch = deg2rad(1.0),
        npts = 8,
    )

    @test length(step.force_values) == mesh.numNodes * 6
    @test step.force_dofs == collect(1:(mesh.numNodes*6))
    @test all(isfinite, step.force_values)
    @test all(isfinite, step.normal_loads)
    @test all(isfinite, step.tangential_loads)
    @test step.resultants.force[1] > 0.0
    @test step.resultants.moment[1] > 0.0
    @test any(abs.(step.state.axial_dynamic) .> 0.0)
    @test any(abs.(step.state.tangential_dynamic) .> 0.0)
    @test state.axial_dynamic ≈ step.state.axial_dynamic
    @test state.tangential_dynamic ≈ step.state.tangential_dynamic
end

@testset "Native HAWT unsteady aeroelastic smoke simulation runs" begin
    mesh, el, feamodel = _shaft_x_hawt_test_structure()
    radial_positions = [1.0, 2.5, 4.5]
    chord = [0.5, 0.42, 0.32]
    twist = deg2rad.([10.0, 6.0, 2.0])
    times = [0.0, 0.1, 0.2]

    result = OWENS.runHAWTUnsteadyAeroelastic(
        mesh,
        el,
        feamodel,
        radial_positions,
        chord,
        twist,
        _hawt_linear_airfoil,
        times;
        rotor_speed = 2.0,
        inflow_speed = t -> 8.0 + t,
        rho = 1.225,
        hub_radius = 0.5,
        tip_radius = 5.0,
        pitch = [deg2rad(1.0), deg2rad(1.0), deg2rad(1.5)],
        npts = 8,
        structural_relaxation = 0.5,
    )

    @test result.times == times
    @test size(result.load_history) == (mesh.numNodes * 6, length(times))
    @test size(result.displacement_history) == (mesh.numNodes * 6, length(times))
    @test all(result.structural_success)
    @test all(isfinite, result.load_history)
    @test all(isfinite, result.displacement_history)
    @test maximum(abs.(result.load_history)) > 0.0
    @test maximum(abs.(result.displacement_history)) > 0.0
    @test result.resultant_force_history[1, end] > 0.0
    @test result.resultant_moment_history[1, end] > 0.0
    @test any(abs.(result.axial_dynamic_history[:, end]) .> 0.0)
    @test result.station_history[:, end] != result.station_history[:, 1]
end

@testset "Native HAWT transient aeroelastic simulation runs" begin
    mesh, el, feamodel = _shaft_x_hawt_test_structure(; analysisType = "TNB")
    radial_positions = [1.0, 2.5, 4.5]
    chord = [0.5, 0.42, 0.32]
    twist = deg2rad.([10.0, 6.0, 2.0])
    times = [0.0, 0.05, 0.10, 0.15]

    result = OWENS.runHAWTUnsteadyAeroelastic(
        mesh,
        el,
        feamodel,
        radial_positions,
        chord,
        twist,
        _hawt_linear_airfoil,
        times;
        rotor_speed = [1.8, 1.9, 2.0, 2.1],
        inflow_speed = 8.0,
        rho = 1.225,
        hub_radius = 0.5,
        tip_radius = 5.0,
        pitch = deg2rad(1.0),
        npts = 8,
        structural_coupling = :transient,
    )

    @test result.times == times
    @test size(result.velocity_history) == (mesh.numNodes * 6, length(times))
    @test size(result.acceleration_history) == (mesh.numNodes * 6, length(times))
    @test size(result.reaction_history) == (mesh.numNodes * 6, length(times))
    @test all(result.structural_success)
    @test all(isfinite, result.load_history)
    @test all(isfinite, result.displacement_history)
    @test all(isfinite, result.velocity_history)
    @test all(isfinite, result.acceleration_history)
    @test maximum(abs.(result.displacement_history[:, 2:end])) > 0.0
    @test maximum(abs.(result.velocity_history[:, 2:end])) > 0.0
    @test maximum(abs.(result.acceleration_history[:, 2:end])) > 0.0
    @test maximum(abs.(result.reaction_history[:, 2:end])) > 0.0
    @test result.station_history[:, end] != result.station_history[:, 1]
    @test !isnothing(result.strain_history[end])
end
