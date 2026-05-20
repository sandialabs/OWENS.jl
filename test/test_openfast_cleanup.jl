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
end
