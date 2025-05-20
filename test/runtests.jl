path = splitdir(@__FILE__)[1]
using Test

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
