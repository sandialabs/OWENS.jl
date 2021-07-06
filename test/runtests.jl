path = splitdir(@__FILE__)[1]
# using OWENS
# include("$(path)/../src/OWENS.jl")
using Test

@testset "_15mTower_transient_dvawt_c_2_lcdt" begin
    include("$path/test_owens_with_fileio.jl")
end
