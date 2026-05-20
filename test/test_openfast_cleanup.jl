using Test
import OWENS

@testset "OpenFAST module cleanup helper" begin
    calls = Symbol[]
    fake_openfast = (
        HD_End = () -> push!(calls, :HD_End),
        MD_End = () -> push!(calls, :MD_End),
        endTurb = () -> push!(calls, :endTurb),
    )

    inputs = (platformActive = true, AD15On = true)
    errors = OWENS.endOpenFASTModules(inputs; openfast = fake_openfast)
    @test calls == [:HD_End, :MD_End, :endTurb]
    @test isempty(errors)

    empty_calls = Symbol[]
    inactive_inputs = (platformActive = false, AD15On = false)
    inactive_openfast = (
        HD_End = () -> push!(empty_calls, :HD_End),
        MD_End = () -> push!(empty_calls, :MD_End),
        endTurb = () -> push!(empty_calls, :endTurb),
    )
    errors = OWENS.endOpenFASTModules(inactive_inputs; openfast = inactive_openfast)
    @test isempty(empty_calls)
    @test isempty(errors)

    calls_with_failure = Symbol[]
    failing_openfast = (
        HD_End = () -> error("hydrodyn close failed"),
        MD_End = () -> push!(calls_with_failure, :MD_End),
        endTurb = () -> push!(calls_with_failure, :endTurb),
    )
    errors =
        OWENS.endOpenFASTModules(inputs; openfast = failing_openfast, warn_on_error = false)
    @test calls_with_failure == [:MD_End, :endTurb]
    @test length(errors) == 1
    @test errors[1].label == "HydroDyn"
    @test errors[1].exception isa ErrorException
    @test occursin("hydrodyn close failed", sprint(showerror, errors[1].exception))

    calls_with_flags = Symbol[]
    errors = OWENS.endOpenFASTModules(
        inputs;
        openfast = (
            HD_End = () -> push!(calls_with_flags, :HD_End),
            MD_End = () -> push!(calls_with_flags, :MD_End),
            endTurb = () -> push!(calls_with_flags, :endTurb),
        ),
        platform_initialized = false,
        ad_initialized = true,
    )
    @test calls_with_flags == [:endTurb]
    @test isempty(errors)

    partial_calls = Symbol[]
    errors = OWENS.endOpenFASTModules(
        inputs;
        openfast = (
            HD_End = () -> push!(partial_calls, :HD_End),
            MD_End = () -> push!(partial_calls, :MD_End),
            endTurb = () -> push!(partial_calls, :endTurb),
        ),
        hd_initialized = true,
        md_initialized = false,
        ad_initialized = false,
    )
    @test partial_calls == [:HD_End]
    @test isempty(errors)
end

@testset "OpenFAST allocation cleanup" begin
    mktempdir() do dir
        md_input_file = joinpath(dir, "MoorDyn_minimal.dat")
        write(
            md_input_file,
            """
            ---------------------- POINTS --------------------------------
            ID Attachment X Y Z M V CdA CA
            1 Fixed 0.0 0.0 -75.0 0 0 0 0
            2 Vessel 0.0 0.0 -5.0 0 0 0 0
            ---------------------- LINES ---------------------------------
            ID LineType AttachA AttachB UnstrLen NumSegs Outputs
            1 main 1 2 80.0 10 -
            ---------------------- SOLVER OPTIONS ------------------------
            """,
        )

        inputs = (
            dataOutputFilename = "none",
            hd_input_file = "none",
            ss_input_file = "none",
            md_input_file = md_input_file,
            potflowfile = nothing,
            interpOrder = 2,
            platformActive = true,
            AD15On = false,
        )
        bottom_mesh = (numNodes = 1,)
        bottom_el = (;)
        bottom_model = (initCond = zeros(0, 3),)
        bin = (hydrodynLibPath = "libhydrodyn", moordynLibPath = "libmoordyn")
        times = range(0.0, length = 2, step = 0.1)

        init_calls = Symbol[]
        init_kwargs = Dict{Symbol,Any}()
        openfast = (
            HD_Init = (; kwargs...) -> begin
                push!(init_calls, :HD_Init)
                init_kwargs[:hd_input_file] = kwargs[:hd_input_file]
                nothing
            end,
            MD_Init = (; kwargs...) -> begin
                push!(init_calls, :MD_Init)
                init_kwargs[:WtrDpth] = kwargs[:WtrDpth]
                nothing
            end,
            HD_End = () -> push!(init_calls, :HD_End),
            MD_End = () -> push!(init_calls, :MD_End),
            endTurb = () -> push!(init_calls, :endTurb),
        )

        result = OWENS.allocate_bottom(
            times,
            2,
            0.1,
            inputs,
            bottom_mesh,
            bottom_el,
            bottom_model,
            bin,
            6;
            openfast = openfast,
        )
        @test result[1] == 6
        @test result[12] == 1
        @test init_calls == [:HD_Init, :MD_Init]
        @test init_kwargs[:hd_input_file] == "none"
        @test init_kwargs[:WtrDpth] == 75.0

        partial_calls = Symbol[]
        md_fails = (
            HD_Init = (; kwargs...) -> push!(partial_calls, :HD_Init),
            MD_Init = (; kwargs...) -> begin
                push!(partial_calls, :MD_Init)
                error("moordyn init failed")
            end,
            HD_End = () -> push!(partial_calls, :HD_End),
            MD_End = () -> push!(partial_calls, :MD_End),
            endTurb = () -> push!(partial_calls, :endTurb),
        )
        @test_throws ErrorException OWENS.allocate_bottom(
            times,
            2,
            0.1,
            inputs,
            bottom_mesh,
            bottom_el,
            bottom_model,
            bin,
            6;
            openfast = md_fails,
        )
        @test partial_calls == [:HD_Init, :MD_Init, :HD_End]

        hd_calls = Symbol[]
        hd_fails = (
            HD_Init = (; kwargs...) -> error("hydrodyn init failed"),
            MD_Init = (; kwargs...) -> push!(hd_calls, :MD_Init),
            HD_End = () -> push!(hd_calls, :HD_End),
            MD_End = () -> push!(hd_calls, :MD_End),
            endTurb = () -> push!(hd_calls, :endTurb),
        )
        @test_throws ErrorException OWENS.allocate_bottom(
            times,
            2,
            0.1,
            inputs,
            bottom_mesh,
            bottom_el,
            bottom_model,
            bin,
            6;
            openfast = hd_fails,
        )
        @test isempty(hd_calls)
    end
end
