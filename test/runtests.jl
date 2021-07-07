path = splitdir(@__FILE__)[1]
using Test

@testset "SNL5MW With CACTUS One-Way Coupling, Preprepared Input Files" begin
    include("$path/SNL5MW_CACT_oneway_unit.jl")
end
