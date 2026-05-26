using Test
import OWENS

@testset "Element reader section-property conventions" begin
    mktempdir() do dir
        elfile = joinpath(dir, "section.el")
        ortfile = joinpath(dir, "section.ort")

        write(
            elfile,
            join(
                [
                    "0.0 0.4 1.0 2.0 3.0 4.0 5.0 6.0 0.7 8.0 9.0 0.0 0.0 0.11 0.12 9.0 0.25",
                    "1.0 0.6 1.5 2.5 3.5 4.5 5.5 6.5 0.8 8.5 9.5 0.0 0.0 0.21 0.22 11.0 0.50",
                ],
                "\n",
            ) * "\n",
        )
        write(ortfile, "1 10.0 20.0 30.0 4.0\n")

        blade_remaining = zeros(2, 12)
        blade_remaining[:, 10] .= [4.0, 6.0]
        blade_remaining[:, 12] .= [5.0, 7.0]
        blade_data =
            OWENS.BladeData(1, [1, 1], [0.0, 1.0], [1, 2], [1, -1], blade_remaining)

        el = OWENS.readElementData(1, elfile, ortfile, blade_data)

        @test el.elLen == [4.0]
        @test el.psi == [10.0]
        @test el.theta == [20.0]
        @test el.roll == [30.0]

        props = el.props[1]
        @test props.zcm == [0.11, 0.21]
        @test props.ycm == [0.12, 0.22]
        @test props.flapwiseEAOffset == [9.0, 11.0]
        @test props.edgewiseEAOffset == [0.25, 0.50]
        @test props.b == [2.0, 3.0]
        @test props.a0 == [5.0, 7.0]
        @test props.ac ≈ [0.2, -0.2]
        @test props.a ≈ [0.25 / 2.0 - 0.5, 0.50 / 3.0 - 0.5]
    end
end
