using Test
import OWENS

function sample_dlc_internal()
    return OWENS.DLC_internal(
        [12.5],
        "U",
        "parked",
        12345,
        17,
        19,
        0.25,
        75.0,
        60.0,
        110.0,
        120.0,
        1.5,
        -2.5,
        "\"IECKAI\"",
        "\"1-ED3\"",
        "\"B\"",
        "\"1EWM50\"",
        75.0,
        12.5,
        [0.0, 1.0],
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

@testset "TurbSim template renderer" begin
    template_lines = [
        "False     Echo        - Echo input data to <RootName>.ech (flag)",
        "40071     RandSeed1   - First random seed (-2147483648 to 2147483647)",
        "30        NumGrid_Z   - Vertical grid-point matrix dimension",
        "0.05      TimeStep    - Time step [s]",
        "58.0      HubHt       - Hub height [m] (should be > 0.5*GridHeight)",
        "0         VFlowAng    - Vertical mean flow (uptilt) angle [degrees]",
        "\"A\"       IECturbc    - IEC turbulence characteristic",
        "\"1ETM\"    IEC_WindType - IEC turbulence type",
        "5.0       URef        - Mean wind speed at the reference height [m/s]",
        "placeholder URefExtra - Similar name that must not match URef",
        "normal    controlStrategy - OWENS-only field that is not a TurbSim key",
        "UF        analysis_type - OWENS-only field that is not a TurbSim key",
        "line without comment separator",
    ]

    rendered_lines = OWENS.renderTurbsimInputLines(template_lines, sample_dlc_internal())

    @test rendered_lines isa Vector{String}
    @test rendered_lines == [
        "False     Echo        - Echo input data to <RootName>.ech (flag)",
        "12345 RandSeed1 - First random seed (-2147483648 to 2147483647)",
        "17 NumGrid_Z - Vertical grid-point matrix dimension",
        "0.25 TimeStep - Time step [s]",
        "75.0 HubHt - Hub height [m] (should be > 0.5*GridHeight)",
        "1.5 VFlowAng - Vertical mean flow (uptilt) angle [degrees]",
        "\"B\" IECturbc - IEC turbulence characteristic",
        "\"1EWM50\" IEC_WindType - IEC turbulence type",
        "12.5 URef - Mean wind speed at the reference height [m/s]",
        "placeholder URefExtra - Similar name that must not match URef",
        "normal    controlStrategy - OWENS-only field that is not a TurbSim key",
        "UF        analysis_type - OWENS-only field that is not a TurbSim key",
        "line without comment separator",
    ]
    @test template_lines[2] ==
          "40071     RandSeed1   - First random seed (-2147483648 to 2147483647)"
end
