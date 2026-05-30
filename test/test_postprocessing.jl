using OWENS
using Test

@testset "Composite ply strain includes shear" begin
    nply = [1]
    tply = [0.2]
    theta = [0.0]
    resultantstrain = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    lowerplystrain, upperplystrain =
        OWENS.my_getplystrain(nply, tply, theta, resultantstrain)

    @test lowerplystrain isa Vector{Vector{Float64}}
    @test upperplystrain isa Vector{Vector{Float64}}
    @test length(lowerplystrain) == 1
    @test length(upperplystrain) == 1
    @test lowerplystrain[1] == [0.5, 1.6, 2.4]
    @test upperplystrain[1] == [1.5, 2.4, 3.6]

    offset = [0.0, 0.25, 0.0]
    lower_offset, upper_offset =
        OWENS.my_getplystrain(nply, tply, theta, resultantstrain, offset)

    @test lower_offset[1] == [2.0, 1.6, 2.4]
    @test upper_offset[1] == [3.0, 2.4, 3.6]
end
