using OWENS
using Test

@testset "OWENS.jl" begin
    atol = 1e-6
    @test isapprox(1.0,1.0;atol)
end
