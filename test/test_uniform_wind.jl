using Test
import OWENS

struct NonIndexableSeries end

function sample_uniform_wind_dlc()
    return OWENS.DLC_internal(
        [12.5],
        "U",
        "normal",
        40071,
        3,
        3,
        0.5,
        30.0,
        10.0,
        20.0,
        20.0,
        0.0,
        0.0,
        "\"IECKAI\"",
        "\"1-ED3\"",
        "\"A\"",
        "\"ECD\"",
        30.0,
        12.5,
        [0.0, 5.0, 10.0],
        nothing,
        [0.0, 10.0, 20.0],
        [0.1, 0.2, 0.3],
        [0.01, 0.02, 0.03],
        [0.2, 0.25, 0.3],
        [0.0, 0.1, 0.2],
        [0.0, 1.0, 0.0],
        [1.0, 1.5, 2.0],
    )
end

@testset "Uniform wind rendering" begin
    expected_lines = [
        "! OpenFAST Deterministic Wind File",
        "#",
        "# Comment lines begin with \"!\" or \"#\" or \"%\", then the data lines must contain the following columns:",
        "#",
        "# If there are only 8 columns, upflow is assumed to be 0.",
        "#",
        "# Parameters are interpolated linearly between time steps; using nearest neighbor before the first time ",
        "# listed in this file and after the last time listed in the file. ",
        "#",
        "! Time     Wind    Wind    Vertical    Horiz.      Pwr.Law     Lin.Vert.   Gust     Upflow",
        "!          Speed   Dir     Speed       Shear       Vert.Shr    Shear       Speed    Angle ",
        "! (sec)    (m/s)   (Deg)   (m/s)                                            (m/s)   (deg)",
        "0.0 12.5 0.0 0.1 0.01 0.2 0.0 0.0 1.0",
        "5.0 12.5 10.0 0.2 0.02 0.25 0.1 1.0 1.5",
        "10.0 12.5 20.0 0.3 0.03 0.3 0.2 0.0 2.0",
    ]

    params = sample_uniform_wind_dlc()
    @test OWENS.renderUniformWindLines(params) == expected_lines

    mktempdir() do dir
        filename = joinpath(dir, "uniform.wnd")
        @test OWENS.generateUniformwind(params, filename) === nothing
        @test readlines(filename) == expected_lines
    end
end

@testset "Uniform wind validation" begin
    params = sample_uniform_wind_dlc()
    params.winddir = nothing
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.gustvel = [0.0, 1.0]
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.windvertvel = [0.1, Inf, 0.3]
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.winddir = 0.0
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.time = Float64[]
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.URef = Inf
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.time = NonIndexableSeries()
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)

    params = sample_uniform_wind_dlc()
    params.winddir = NonIndexableSeries()
    @test_throws ArgumentError OWENS.renderUniformWindLines(params)
end
