path = splitdir(@__FILE__)[1]
# using OWENS
# include("$(path)/../src/OWENS.jl")
using Test

@testset "Mesh Check" begin
    include("$path/check_mesh.jl")
end

@testset "_15mTower_transient_dvawt_c_2_lcdt" begin
    include("$path/test_owens.jl")
end

@testset "Run Hydro" begin
    include("$path/simple_tlp_unit_test.jl")
end
