using Test
import OWENS

function minimal_dlc_model_options()
    return (
        OWENS_Options = (numTS = 200, delta_t = 0.1),
        Mesh_Options = (ntelem = 4, nbelem = 5),
    )
end

function minimal_dlc_design_parameters()
    return Dict(
        :assembly => Dict(:hub_height => 30.0),
        :components => Dict(
            :blade => Dict(
                :outer_shape_bem => Dict(
                    :reference_axis => Dict(
                        :x => Dict(:values => [2.0, 3.0, 5.0]),
                        :y => Dict(:values => [0.0, 4.0, 0.0]),
                        :z => Dict(:values => [0.0, 10.0, 20.0]),
                    ),
                ),
            ),
        ),
    )
end

function dlc_params(case; Vref = 42.0, Vinf_range = [5.0, 10.0, 15.0])
    return OWENS.getDLCparams(
        case,
        minimal_dlc_model_options(),
        minimal_dlc_design_parameters(),
        Vinf_range,
        11.0,
        Vref,
        "\"A\"",
        1,
        "\"1-ED3\"";
        grid_oversize = 1.1,
        simtime_turbsim = nothing,
        delta_t_turbsim = nothing,
        NumGrid_Z = nothing,
        NumGrid_Y = nothing,
        RandSeed1 = 40071,
    )
end

@testset "IEC extreme wind speeds from Vref" begin
    speeds = OWENS.iecExtremeWindSpeeds(42.0)

    @test speeds isa NamedTuple{(:Ve50, :Ve1),Tuple{Float64,Float64}}
    @test speeds.Ve50 == 58.8
    @test speeds.Ve1 == 47.04

    custom_speeds = OWENS.iecExtremeWindSpeeds(10.0; ve50_factor = 1.2, ve1_factor = 0.75)
    @test custom_speeds isa NamedTuple{(:Ve50, :Ve1),Tuple{Float64,Float64}}
    @test custom_speeds.Ve50 == 12.0
    @test custom_speeds.Ve1 == 9.0

    @test_throws ArgumentError OWENS.iecExtremeWindSpeeds(0.0)
    @test_throws ArgumentError OWENS.iecExtremeWindSpeeds(-1.0)
    @test_throws ArgumentError OWENS.iecExtremeWindSpeeds(Inf)
    @test_throws ArgumentError OWENS.iecExtremeWindSpeeds(10.0; ve50_factor = 0.0)
    @test_throws ArgumentError OWENS.iecExtremeWindSpeeds(10.0; ve1_factor = -0.5)
end

@testset "DLC parked extreme winds use Vref" begin
    dlc_6_1 = dlc_params("6_1")
    @test dlc_6_1.Vinf_range_used isa Vector{Float64}
    @test dlc_6_1.Vinf_range_used == [58.8]
    @test dlc_6_1.analysis_type == "U"
    @test dlc_6_1.controlStrategy == "parked"
    @test dlc_6_1.IEC_WindType == "\"1EWM50\""
    @test dlc_6_1.GridHeight == 22.0
    @test dlc_6_1.GridWidth == 11.0
    @test dlc_6_1.HubHt == 30.0
    @test dlc_6_1.TimeStep == 0.1
    @test dlc_6_1.AnalysisTime == 20.0

    dlc_6_2 = dlc_params("6_2")
    @test dlc_6_2.Vinf_range_used isa Vector{Float64}
    @test dlc_6_2.Vinf_range_used == [58.8]
    @test dlc_6_2.analysis_type == "U"
    @test dlc_6_2.controlStrategy == "parked_idle"
    @test dlc_6_2.IEC_WindType == "\"1EWM50\""

    dlc_6_3 = dlc_params("6_3")
    @test dlc_6_3.Vinf_range_used isa Vector{Float64}
    @test dlc_6_3.Vinf_range_used == [47.04]
    @test dlc_6_3.analysis_type == "U"
    @test dlc_6_3.controlStrategy == "parked_yaw"
    @test dlc_6_3.IEC_WindType == "\"1EWM1\""

    dlc_6_4 = dlc_params("6_4")
    @test dlc_6_4.Vinf_range_used isa Vector{Float64}
    @test dlc_6_4.Vinf_range_used == [41.16]
    @test dlc_6_4.analysis_type == "F"
    @test dlc_6_4.controlStrategy == "parked"
    @test dlc_6_4.IEC_WindType == "\"1NTM\""

    dlc_7_1 = dlc_params("7_1")
    @test dlc_7_1.Vinf_range_used isa Vector{Float64}
    @test dlc_7_1.Vinf_range_used == [47.04]
    @test dlc_7_1.analysis_type == "U"
    @test dlc_7_1.controlStrategy == "parked"
    @test dlc_7_1.IEC_WindType == "\"1EWM1\""
end

@testset "DLC normal-operation wind range remains explicit" begin
    requested_range = [4.0, 8.0, 12.0]
    dlc_1_1 = dlc_params("1_1"; Vref = 99.0, Vinf_range = requested_range)

    @test dlc_1_1.Vinf_range_used === requested_range
    @test dlc_1_1.Vinf_range_used isa Vector{Float64}
    @test dlc_1_1.analysis_type == "UF"
    @test dlc_1_1.controlStrategy == "normal"
    @test dlc_1_1.IEC_WindType == "\"1NTM\""
end

@testset "IEC 61400-1 DLC branch table" begin
    expected_cases = [
        ("1_2", "normal", [5.0, 10.0, 15.0], "UF", "\"1NTM\""),
        ("1_3", "normal", [5.0, 10.0, 15.0], "U", "\"1ETM\""),
        ("1_4", "normal", [9.0, 13.0], "U", "\"1ECD\""),
        ("1_5", "normal", [5.0, 10.0, 15.0], "U", "\"1EWS\""),
        ("2_1", "freewheelatNormalOperatingRPM", [5.0, 10.0, 15.0], "U", "\"1NTM\""),
        ("2_2", "freewheelatNormalOperatingRPM", [5.0, 10.0, 15.0], "U", "\"1NTM\""),
        ("2_3", "freewheelatNormalOperatingRPM", [9.0, 13.0, 15.0], "U", "\"1EOG\""),
        ("2_4", "freewheelatNormalOperatingRPM", [5.0, 10.0, 15.0], "U", "\"1NTM\""),
        ("3_1", "startup", [5.0, 10.0, 15.0], "F", "\"1NWP\""),
        ("3_2", "startup", [5.0, 9.0, 13.0, 15.0], "U", "\"1EOG\""),
        ("3_3", "startup", [5.0, 9.0, 13.0, 15.0], "U", "\"1ECD\""),
        ("4_1", "shutdown", [5.0, 10.0, 15.0], "F", "\"1NWP\""),
        ("4_2", "shutdown", [9.0, 13.0, 15.0], "U", "\"1EOG\""),
        ("5_1", "emergencyshutdown", [9.0, 13.0, 15.0], "U", "\"1NTM\""),
        ("8_1", "transport", [11.0], "U", "\"1EWM1\""),
        ("CPCurve", "normal", [5.0, 10.0, 15.0], "F", "\"1NWP\""),
    ]

    for (case, strategy, wind_range, analysis_type, wind_type) in expected_cases
        params = dlc_params(case)
        @test params.controlStrategy == strategy
        @test collect(params.Vinf_range_used) == wind_range
        @test params.analysis_type == analysis_type
        @test params.IEC_WindType == wind_type
        @test params.RandSeed1 == 40071
        @test params.NumGrid_Z == 9
        @test params.NumGrid_Y == 5
        @test params.RefHt == 30.0
    end

    ecd = dlc_params("1_4")
    @test ecd.time == [0.0, 10.0, 15.0, 20.0, 25.0, 30.0, 10000.0]
    @test ecd.winddir == [0.0, 0.0, 45.0, 90.0, 45.0, 0.0, 0.0]
    @test ecd.gustvel == [0.0, 0.0, 7.0, 15.0, 7.5, 0.0, 0.0]
    @test ecd.UpflowAngle == zeros(7)

    ews = dlc_params("1_5")
    @test ews.horizshear == [0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0]
    @test ews.LinVertShear == [0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0]
    @test ews.gustvel == zeros(7)

    eog = dlc_params("2_3")
    @test eog.time isa LinRange{Float64,Int64}
    @test length(eog.time) == 100
    @test first(eog.time) == 0.0
    @test last(eog.time) == 30.0
    @test eog.gustvel[1] == 0.0
    @test eog.gustvel[50] ≈ 10.962142970769854 atol=1e-14
    @test eog.gustvel[67] == 0.0
    @test eog.gustvel[83] ≈ 10.962142970769854 atol=1e-14

    nwp = dlc_params("3_1")
    @test collect(nwp.time) == collect(LinRange(0.0, 30.0, 10))
    @test nwp.pwrLawVertShear == fill(0.2, 10)
    @test nwp.winddir == zeros(10)

    @test_throws ErrorException dlc_params("9_9")
    @test_throws ErrorException OWENS.getDLCparams(
        "1_1",
        minimal_dlc_model_options(),
        minimal_dlc_design_parameters(),
        [5.0],
        11.0,
        42.0,
        "\"A\"",
        1,
        "\"3-ED1\"";
    )
    @test_throws ErrorException OWENS.getDLCparams(
        "1_1",
        minimal_dlc_model_options(),
        minimal_dlc_design_parameters(),
        [5.0],
        11.0,
        42.0,
        "\"A\"",
        1,
        "\"2\"";
    )
end

@testset "DLC gust and WindIO control helpers" begin
    @test OWENS.simpleGustVel(5.0, 10.0, 15.0, 10.0) == 0.0
    @test OWENS.simpleGustVel(10.0, 10.0, 15.0, 10.0) == 0.0
    @test OWENS.simpleGustVel(15.0, 10.0, 15.0, 10.0) ≈ 11.1 atol=1e-14
    @test OWENS.simpleGustVel(25.0, 10.0, 15.0, 10.0) == 0.0

    @test OWENS.getGustVel(0.0, 12.0, 4.0, 3.0, 10.0, 8.0) == 12.0
    @test OWENS.getGustVel(17.0, 12.0, 4.0, 3.0, 10.0, 8.0) ≈ 14.22 atol=1e-14

    prescribed_time = [0.0, 5.0, 10.0]
    prescribed_rpm = [60.0, 90.0, 120.0]
    strategy_expectations = Dict(
        "normal" => ([0.0, 10.0, 30.0, 1000.0], fill(2.0, 4), 0),
        "freewheelatNormalOperatingRPM" => ([0.0, 10.0, 30.0, 1000.0], fill(2.0, 4), 2),
        "startup" => ([0.0, 30.0, 1000.0], [0.5, 2.0, 2.0], 0),
        "shutdown" => ([0.0, 10.0, 15.0, 30.0, 1000.0], [2.0, 2.0, 1.0, 0.5, 0.0], 0),
        "emergencyshutdown" => ([0.0, 10.0, 15.0, 30.0, 1000.0], [2.0, 2.0, 0.5, 0.0, 0.0], 0),
        "parked" => ([0.0, 10.0, 30.0, 1000.0], fill(0.5, 4), 0),
        "parked_idle" => ([0.0, 10.0, 30.0, 1000.0], fill(0.5, 4), 2),
        "parked_yaw" => ([0.0, 10.0, 30.0, 1000.0], fill(0.5, 4), 0),
        "transport" => ([0.0, 10.0, 30.0, 1000.0], fill(0.5, 4), 0),
    )

    for (strategy, (tocp_expected, omega_expected, startup_expected)) in strategy_expectations
        tocp, Omegaocp, turbineStartup, generatorOn, useGeneratorFunction =
            OWENS.handle_control_strategy(strategy, 120.0, 30.0, prescribed_time, prescribed_rpm)
        @test tocp == tocp_expected
        @test Omegaocp == omega_expected
        @test turbineStartup == startup_expected
        @test generatorOn === false
        @test useGeneratorFunction === false
    end

    tocp, Omegaocp, turbineStartup, generatorOn, useGeneratorFunction =
        OWENS.handle_control_strategy(
            "prescribedRPM",
            120.0,
            30.0,
            prescribed_time,
            prescribed_rpm,
        )
    @test tocp === prescribed_time
    @test Omegaocp == [1.0, 1.5, 2.0]
    @test turbineStartup == 0
    @test generatorOn === false
    @test useGeneratorFunction === false

    unknown = @test_logs (:warn, "ControlStrategy unsupported not recognized, using prescribed RPM spline control points") OWENS.handle_control_strategy(
        "unsupported",
        120.0,
        30.0,
        prescribed_time,
        prescribed_rpm,
    )
    @test unknown == (prescribed_time, [1.0, 1.5, 2.0], 0, false, false)
end
