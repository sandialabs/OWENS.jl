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
