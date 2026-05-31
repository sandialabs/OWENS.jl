using Test
import OWENS

mutable struct RecordingAxis
    plots::Vector{Tuple}
    labels::Vector{Pair{Symbol,String}}
    views::Vector{Tuple{Any,Any}}
    grids::Vector{Any}
end

RecordingAxis() = RecordingAxis(Tuple[], Pair{Symbol,String}[], Tuple{Any,Any}[], Any[])

function Base.getproperty(axis::RecordingAxis, name::Symbol)
    if name === :plot
        return (args...) -> (push!(getfield(axis, :plots), Tuple(args)); nothing)
    elseif name === :set_xlabel
        return label -> (push!(getfield(axis, :labels), :x => string(label)); nothing)
    elseif name === :set_ylabel
        return label -> (push!(getfield(axis, :labels), :y => string(label)); nothing)
    elseif name === :set_zlabel
        return label -> (push!(getfield(axis, :labels), :z => string(label)); nothing)
    elseif name === :view_init
        return (azimuth, elevation) ->
            (push!(getfield(axis, :views), (azimuth, elevation)); nothing)
    elseif name === :grid
        return setting -> (push!(getfield(axis, :grids), setting); nothing)
    end
    return getfield(axis, name)
end

mutable struct RecordingFigure
    axes::Vector{RecordingAxis}
end

RecordingFigure() = RecordingFigure(RecordingAxis[])

function Base.getproperty(figure::RecordingFigure, name::Symbol)
    if name === :add_subplot
        return (args...; kwargs...) -> begin
            axis = RecordingAxis()
            push!(getfield(figure, :axes), axis)
            return axis
        end
    end
    return getfield(figure, name)
end

mutable struct RecordingPyPlot
    figure_handle::RecordingFigure
    titles::Vector{String}
    saved::Vector{Pair{String,Bool}}
end

RecordingPyPlot() = RecordingPyPlot(RecordingFigure(), String[], Pair{String,Bool}[])

function Base.getproperty(pyplot::RecordingPyPlot, name::Symbol)
    if name === :figure
        return _ -> getfield(pyplot, :figure_handle)
    elseif name === :title
        return title -> (push!(getfield(pyplot, :titles), string(title)); nothing)
    elseif name === :savefig
        return (path; transparent = false) ->
            (push!(getfield(pyplot, :saved), string(path) => transparent); nothing)
    end
    return getfield(pyplot, name)
end

mutable struct SimpleVizMesh
    conn::Matrix{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    meshSeg::Vector{Int}
end

@testset "plotMesh records segment coordinates and labels" begin
    mesh = SimpleVizMesh(
        [1 2; 2 3; 3 4],
        [0.0, 1.0, 2.0, 3.0],
        [0.0, 0.5, 1.0, 1.5],
        [0.0, 0.0, 0.0, 1.0],
        [2, 5],
    )
    ax1, ax2, ax3, ax4 = (RecordingAxis(), RecordingAxis(), RecordingAxis(), RecordingAxis())

    OWENS.plotMesh(mesh, ".k", mesh.meshSeg, ax1, ax2, ax3, ax4)

    @test length(ax1.plots) == 2
    @test ax1.plots[1][1] == [0.0, 1.0, 1.0, 2.0]
    @test ax1.plots[1][2] == [0.0, 0.0, 0.0, 0.0]
    @test ax2.plots[1][1] == [0.0, 1.0, 1.0, 2.0]
    @test ax2.plots[1][2] == [0.0, 0.5, 0.5, 1.0]
    @test ax3.plots[1][1] == [0.0, 0.5, 0.5, 1.0]
    @test ax3.plots[1][2] == [0.0, 0.0, 0.0, 0.0]
    @test ax4.plots[1][1] == [0.0, 1.0, 1.0, 2.0]
    @test ax4.plots[1][2] == [0.0, 0.5, 0.5, 1.0]
    @test ax4.plots[1][3] == [0.0, 0.0, 0.0, 0.0]
    @test ax4.plots[1][4] == ".k"
    @test (:x => "h_1") in ax1.labels
    @test (:y => "h_3") in ax1.labels
    @test (:z => "h_3") in ax4.labels
end

@testset "viz deforms modal meshes before plotting" begin
    mesh = SimpleVizMesh(
        [1 2; 2 3],
        [0.0, 1.0, 2.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],
        [2],
    )
    nmodes = 2
    nnodes = length(mesh.x)
    zeros_modes = zeros(nnodes, nmodes)
    U_x_0 = copy(zeros_modes)
    U_y_90 = copy(zeros_modes)
    U_x_0[2, 2] = 0.25
    U_y_90[3, 2] = -0.5
    pyplot = RecordingPyPlot()

    OWENS.viz(
        pyplot;
        mesh,
        selectedMode = 2,
        sf = 4.0,
        freq = [0.0, 1.5],
        damp = zeros(nmodes),
        U_x_0,
        U_y_0 = zeros_modes,
        U_z_0 = zeros_modes,
        theta_x_0 = zeros_modes,
        theta_y_0 = zeros_modes,
        theta_z_0 = zeros_modes,
        U_x_90 = zeros_modes,
        U_y_90,
        U_z_90 = zeros_modes,
        theta_x_90 = zeros_modes,
        theta_y_90 = zeros_modes,
        theta_z_90 = zeros_modes,
        savename = "mode-shape.pdf",
    )

    @test pyplot.titles == ["Frequency 1.5"]
    @test pyplot.saved == ["mode-shape.pdf" => true]
    @test length(pyplot.figure_handle.axes) == 4
    @test length(pyplot.figure_handle.axes[1].plots) == 3
    @test pyplot.figure_handle.axes[1].plots[1][1] == [0.0, 2.0, 2.0, 2.0]
    @test pyplot.figure_handle.axes[2].plots[2][2] == [0.0, 0.0, 0.0, -1.0]
    @test pyplot.figure_handle.axes[4].views == [(45, 45)]
    @test pyplot.figure_handle.axes[4].grids == ["off"]
end

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
