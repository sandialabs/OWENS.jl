path,_ = splitdir(@__FILE__)
using OWENS
# include("$(path)/../src/OWENS.jl")
using Test

@testset "OWENS.jl" begin
    atol = 1e-6
    @test isapprox(1.0,1.0;atol)
end

@testset "Mesh Check" begin
    include("$path/check_mesh.jl")
end

@testset "_15mTower_transient_dvawt_c_2_lcdt" begin
    include("$path/test_owens.jl")
end

@testset "Run Hydro" begin
    p = OWENS.hydro.tlp_platform(r_spar=2,draft=30, height=5, width=2, length=10, num=3, ofst=1)
    p.make_mesh(
        mshRefFactor=1,
        # clcurv=360/200,
        write_file=true)
    p.run_hydro()
    println("We ran, and so for now we pass!")
    @test true
end
