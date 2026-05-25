using Test
using OWENS
using OrderedCollections: OrderedDict

@testset "legacy Time Dean option is rejected" begin
    @test_throws ArgumentError OWENS.Inputs(; analysisType = "TD")
    @test_throws ArgumentError OWENS.Inputs(; analysisType = "time_dean")
    @test_throws ArgumentError OWENS.MasterInput(; structuralModel = "TD")
    @test_throws ArgumentError OWENS.OWENS_Options(
        OrderedDict{Symbol,Any}(:structuralModel=>"TD"),
    )

    mktempdir() do dir
        options_file = joinpath(dir, "modeling_options.yml")
        write(
            options_file,
            """
            OWENS_Options:
              structuralModel: TD
            """,
        )

        @test_throws ArgumentError OWENS.ModelingOptions(options_file)
    end
end
