using Test
import OWENS

@testset "MoorDyn tension buffers" begin
    buffer = OWENS.mooringTensionBuffer(3)
    @test buffer isa Vector{Float32}
    @test length(buffer) == 6

    empty_buffer = OWENS.mooringTensionBuffer(0)
    @test empty_buffer isa Vector{Float32}
    @test isempty(empty_buffer)

    existing_buffer = zeros(Float32, 8)
    @test OWENS.validateMooringTensionBuffer(existing_buffer) === existing_buffer

    @test_throws ArgumentError OWENS.mooringTensionBuffer(-1)
    @test_throws MethodError OWENS.mooringTensionBuffer(1.5)
    @test_throws ArgumentError OWENS.validateMooringTensionBuffer(zeros(Float32, 5))
    @test_throws ArgumentError OWENS.validateMooringTensionBuffer((1.0, 2.0))
end
