using Test
import OWENS

const test_path = splitdir(@__FILE__)[1]

function hydrodyn_potfile_line(filename)
    return only(filter(line -> occursin(r"\bPotFile\b", line), readlines(filename)))
end

@testset "HydroDyn PotFile resolver" begin
    hd_input_file = joinpath(test_path, "data", "HydroDyn_CCT2_test.dat")
    staged_file = OWENS.hydrodynInputWithResolvedPotFile(hd_input_file, nothing)
    expected_root =
        abspath(joinpath(test_path, "data", "potential_flow_data", "marin_semi"))

    @test staged_file isa String
    @test staged_file != hd_input_file
    @test isfile(staged_file)
    @test occursin("\"$(expected_root)\"", hydrodyn_potfile_line(staged_file))
    @test OWENS._hydrodynPotentialFlowRootExists(expected_root) === true

    mktempdir() do dir
        input_file = joinpath(dir, "HydroDyn_override.dat")
        override_root = joinpath(dir, "explicit", "override_case")
        write(
            input_file,
            """
            ------- POTENTIAL FLOW INPUTS --------------------------------
            "relative/default"    PotFile       - Root name of potential-flow model data
            """,
        )

        staged_file = OWENS.hydrodynInputWithResolvedPotFile(input_file, override_root)
        @test occursin("\"$(abspath(override_root))\"", hydrodyn_potfile_line(staged_file))
    end

    mktempdir() do dir
        no_potfile_input = joinpath(dir, "HydroDyn_no_potfile.dat")
        write(no_potfile_input, "False Echo - Echo input data to <RootName>.ech\n")
        @test OWENS.hydrodynInputWithResolvedPotFile(no_potfile_input, nothing) ==
              no_potfile_input

        unused_potfile_input = joinpath(dir, "HydroDyn_unused_potfile.dat")
        write(
            unused_potfile_input,
            "\"unused\"    PotFile       - Root name of potential-flow model data\n",
        )
        @test OWENS.hydrodynInputWithResolvedPotFile(unused_potfile_input, nothing) ==
              unused_potfile_input

        malformed_potfile_input = joinpath(dir, "HydroDyn_malformed_potfile.dat")
        write(
            malformed_potfile_input,
            "      PotFile       - missing potential-flow root\n",
        )
        @test_throws ArgumentError OWENS.hydrodynInputWithResolvedPotFile(
            malformed_potfile_input,
            nothing,
        )

        unresolved_relative =
            OWENS._resolvedHydroDynPotFileRoot(
                "missing_case    PotFile\n",
                no_potfile_input,
                nothing,
            )
        @test unresolved_relative == abspath(joinpath(dir, "missing_case"))
        @test_throws ArgumentError OWENS._hydrodynPotFileLineWithRoot(
            "      PotFile       - missing potential-flow root\n",
            "unused",
        )
    end

    @test OWENS.hydrodynInputWithResolvedPotFile("none", nothing) == "none"
end
