using Test
import OWENS

@testset "VTK stress history selection" begin
    full_stress = reshape(collect(1.0:24.0), 4, 3, 2)
    tsave_idx = [1, 3]

    saved_stress = OWENS._saved_vtk_stress_history(full_stress, tsave_idx, 4)
    @test size(saved_stress) == (2, 3, 2)
    @test saved_stress[1, 1, 1] == 1.0
    @test saved_stress[1, 3, 2] == full_stress[1, 3, 2]
    @test saved_stress[2, 2, 1] == full_stress[3, 2, 1]
    @test saved_stress[2, 3, 2] == 23.0

    presliced_stress = full_stress[tsave_idx, :, :]
    @test OWENS._saved_vtk_stress_history(presliced_stress, tsave_idx, 4) ===
          presliced_stress
    @test isnothing(OWENS._saved_vtk_stress_history(nothing, tsave_idx, 4))

    malformed_stress = reshape(collect(1.0:18.0), 3, 3, 2)
    @test_throws DimensionMismatch OWENS._saved_vtk_stress_history(
        malformed_stress,
        tsave_idx,
        4,
    )
    @test_throws ArgumentError OWENS._saved_vtk_stress_history(ones(2, 3), tsave_idx, 4)
    @test_throws ArgumentError OWENS._saved_vtk_stress_history(full_stress, [1.0, 3.0], 4)
    @test_throws ArgumentError OWENS._saved_vtk_stress_history(full_stress, [true, false], 4)
    @test_throws ArgumentError OWENS._saved_vtk_stress_history(full_stress, [1, 5], 4)
end

@testset "VTK history input validation" begin
    assembly = (; points = [1, 2, 3], elements = [1, 2])
    history = [:state1, :state2]
    time = [0.0, 0.1]
    sections = zeros(3, 2, 3)
    user_point_names = ["load", "Principal_Surface_Layer_Stress"]
    user_point_data = zeros(2, 2, 3)
    stress = zeros(2, 3, 2)

    @test OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        theta_z = [0.0, pi / 2],
        userPointNames = user_point_names,
        userPointData = user_point_data,
        stress,
    ) === nothing
    @test OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = user_point_names,
        userPointData = user_point_data,
        stress = nothing,
    ) === nothing

    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        [:state1],
        time,
        assembly,
        sections,
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        theta_z = [0.0],
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        zeros(2, 2, 3),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        zeros(3, 2, 2),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = user_point_names,
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointData = user_point_data,
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = ["load"],
        userPointData = user_point_data,
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = user_point_names,
        userPointData = zeros(2, 2),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = user_point_names,
        userPointData = zeros(2, 1, 3),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        userPointNames = user_point_names,
        userPointData = zeros(2, 2, 2),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        stress = zeros(2, 3),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        stress = zeros(1, 3, 2),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        stress = zeros(2, 2, 2),
    )
    @test_throws ArgumentError OWENS._validate_vtk_history_inputs(
        history,
        time,
        assembly,
        sections;
        stress = zeros(2, 3, 1),
    )
end
