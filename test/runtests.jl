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

@testset "Mesh Check" begin
    include("$path/check_mesh.jl")
end
