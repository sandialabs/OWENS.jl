path = splitdir(@__FILE__)[1]
using Test

@testset "DLC Parameter Helpers" begin
    include("$path/test_dlc.jl")
end

@testset "TurbSim Template Rendering" begin
    include("$path/test_turbsim_template.jl")
end

@testset "Uniform Wind Rendering" begin
    include("$path/test_uniform_wind.jl")
end

@testset "Initial Conditions" begin
    include("$path/test_initial_conditions.jl")
end

@testset "MoorDyn Buffer Helpers" begin
    include("$path/test_moordyn_buffers.jl")
end

@testset "HydroDyn PotFile Resolver" begin
    include("$path/test_hydrodyn_potfile.jl")
end

@testset "Unsteady Completed History Ranges" begin
    include("$path/test_unsteady_history_ranges.jl")
end

@testset "OpenFAST Module Cleanup" begin
    include("$path/test_openfast_cleanup.jl")
end

@testset "Run Manifest Provenance" begin
    include("$path/test_run_manifest.jl")
end

@testset "VTK History Output" begin
    include("$path/test_visualization_history.jl")
end

@testset "SNL5MW With CACTUS One-Way Coupling, Preprepared Input Files" begin
    include("$path/SNL5MW_CACT_oneway_unit.jl")
end

@testset "SNL34m With impulse load ROM test, Preprepared Input Files" begin
    include("$path/34mROMtest.jl")
end

@testset "SNL34m Campbell Diagram Test Against GXBeam" begin
    include("$path/Fig4_5_campbell2.jl")
end

@testset "Rainflow Check" begin
    include("$path/testrainflow.jl")
end

@testset "Mesh Check" begin
    include("$path/check_mesh.jl")
end

@testset "Floating Platform Check" begin
    include("$path/testfloating.jl")
end

@testset "Added Mass" begin
    include("$path/../examples/AddedMass_Buoyancy/Flapping_Added_Mass.jl")
end

@testset "Buoyancy" begin
    path = splitdir(@__FILE__)[1]
    include("$path/../examples/AddedMass_Buoyancy/Buoyancy.jl")
end

@testset "WindIO" begin
    path = splitdir(@__FILE__)[1]
    include("$path/../examples/Optimization/windio_example.jl")
end
