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
