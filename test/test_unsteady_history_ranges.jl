using Test
import OWENS

@testset "Unsteady completed history ranges" begin
    ranges = OWENS.completedHistoryRanges(1, 5)
    @test ranges.state == 1:1
    @test isempty(ranges.step)

    ranges = OWENS.completedHistoryRanges(4, 5)
    @test ranges.state == 1:4
    @test ranges.step == 1:3

    ranges = OWENS.completedHistoryRanges(5, 5)
    @test ranges.state == 1:5
    @test ranges.step == 1:4

    @test_throws ArgumentError OWENS.completedHistoryRanges(0, 5)
    @test_throws ArgumentError OWENS.completedHistoryRanges(6, 5)
    @test_throws ArgumentError OWENS.completedHistoryRanges(1, 0)
end
