using Test
import OWENS

const test_path = splitdir(@__FILE__)[1]

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

@testset "MoorDyn input summary" begin
    summary =
        OWENS.moorDynInputSummary(joinpath(test_path, "data", "MoorDyn_CCT2_test.dat"))
    @test summary isa NamedTuple{(:numMooringLines, :waterDepth),Tuple{Int64,Float64}}
    @test summary.numMooringLines == 3
    @test summary.waterDepth == 200.0

    mktempdir() do dir
        input_file = joinpath(dir, "MoorDyn_two_line.dat")
        write(
            input_file,
            """
            ---------------------- POINTS --------------------------------
            ID Attachment X Y Z M V CdA CA
            (-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
            1 Fixed 0.0 0.0 -50.0 0 0 0 0
            2 Vessel 0.0 0.0 -5.0 0 0 0 0
            3 Fixed 10.0 0.0 -75.0 0 0 0 0
            4 Vessel 10.0 0.0 -5.0 0 0 0 0
            ---------------------- LINES ---------------------------------
            ID LineType AttachA AttachB UnstrLen NumSegs Outputs
            (-) (-) (-) (-) (m) (-) (-)
            1 main 1 2 80.0 10 -
            2 main 3 4 90.0 10 -
            ---------------------- SOLVER OPTIONS ------------------------
            """,
        )

        summary = OWENS.moorDynInputSummary(input_file)
        @test summary.numMooringLines == 2
        @test summary.waterDepth == 75.0
    end

    mktempdir() do dir
        input_file = joinpath(dir, "MoorDyn_no_fixed.dat")
        write(
            input_file,
            """
            ---------------------- POINTS --------------------------------
            ID Attachment X Y Z M V CdA CA
            1 Vessel 0.0 0.0 -5.0 0 0 0 0
            ---------------------- LINES ---------------------------------
            ID LineType AttachA AttachB UnstrLen NumSegs Outputs
            1 main 1 1 80.0 10 -
            """,
        )
        @test_throws ArgumentError OWENS.moorDynInputSummary(input_file)
    end

    mktempdir() do dir
        input_file = joinpath(dir, "MoorDyn_no_lines.dat")
        write(
            input_file,
            """
            ---------------------- POINTS --------------------------------
            ID Attachment X Y Z M V CdA CA
            1 Fixed 0.0 0.0 -50.0 0 0 0 0
            ---------------------- LINES ---------------------------------
            ID LineType AttachA AttachB UnstrLen NumSegs Outputs
            """,
        )
        @test_throws ArgumentError OWENS.moorDynInputSummary(input_file)
    end
end
