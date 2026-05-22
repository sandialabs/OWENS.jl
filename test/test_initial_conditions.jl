using Test
import OWENS

@testset "Initial condition arrays" begin
    legacy_conditions = OWENS.createInitCondArray([0.0, 2.5, 0.0, -1.0], 3, 4)
    @test legacy_conditions isa Matrix{Float64}
    @test legacy_conditions == [
        1.0 2.0 2.5
        2.0 2.0 2.5
        3.0 2.0 2.5
        1.0 4.0 -1.0
        2.0 4.0 -1.0
        3.0 4.0 -1.0
    ]

    explicit_node_conditions = OWENS.createInitCondArray([1.0, 0.0, -0.5], [10, 20], 6)
    @test explicit_node_conditions isa Matrix{Float64}
    @test explicit_node_conditions == [
        10.0 1.0 1.0
        20.0 1.0 1.0
        10.0 3.0 -0.5
        20.0 3.0 -0.5
    ]

    @test OWENS.createInitCondArray(zeros(3), 2, 3) == []
    @test OWENS.createInitCondArray(zeros(3), [10, 20], 3) == []

    @test_throws ArgumentError OWENS.createInitCondArray([1.0], -1, 1)
    @test_throws ArgumentError OWENS.createInitCondArray([1.0, 2.0], [1, 2], 1)
    @test_throws ArgumentError OWENS.createInitCondArray([1.0], [1.5], 1)
end
