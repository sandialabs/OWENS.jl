using OrderedCollections: OrderedDict
using YAML

export DLC_internal,
    runDLC, getDLCparams, generateUniformwind, generateTurbsimBTS, getGustVel, simpleGustVel

"""
    DLC_internal

Internal structure for Design Load Case (DLC) analysis parameters and results.

# Fields
* `Vinf_range_used::Array{Float64}`: Range of wind speeds used in the analysis
* `analysis_type::String`: Type of analysis ("U"=unsteady, "F"=fatigue, "UF"=unsteady fatigue)
* `controlStrategy::String`: Control strategy type
* `RandSeed1::Int`: Random seed for turbulent wind generation
* `NumGrid_Z::Int`: Number of vertical grid points
* `NumGrid_Y::Int`: Number of horizontal grid points
* `TimeStep::Float64`: Time step size in seconds
* `HubHt::Float64`: Hub height in meters
* `AnalysisTime::Float64`: Total analysis time in seconds
* `GridHeight::Float64`: Grid height in meters
* `GridWidth::Float64`: Grid width in meters
* `VFlowAng::Float64`: Vertical mean flow angle in degrees
* `HFlowAng::Float64`: Horizontal mean flow angle in degrees
* `TurbModel::String`: Turbulence model specification
* `IECstandard::String`: IEC standard reference (e.g., "1-ED3")
* `IECturbc::String`: IEC turbulence characteristic ("A", "B", "C" or TI%)
* `IEC_WindType::String`: IEC turbulence type ("NTM", "xETM", "xEWM1", "xEWM50")
* `RefHt::Float64`: Reference height for wind speed in meters
* `URef::Float64`: Reference wind speed in m/s
* `time::Array{Float64}`: Time series array
* `windvel::Array{Float64}`: Wind velocity time series
* `winddir::Array{Float64}`: Wind direction time series
* `windvertvel::Array{Float64}`: Vertical wind velocity time series
* `horizshear::Array{Float64}`: Horizontal wind shear profile
* `pwrLawVertShear::Float64`: Power law vertical shear coefficient
* `LinVertShear::Float64`: Linear vertical shear coefficient
* `gustvel::Array{Float64}`: Gust velocity time series
* `UpflowAngle::Float64`: Upflow angle in degrees
"""
mutable struct DLC_internal
    Vinf_range_used::Any
    analysis_type::Any # "U", "F", "UF"
    controlStrategy::Any # "constRPM", function handle
    RandSeed1::Any # Turbulent Random Seed Number
    NumGrid_Z::Any # Vertical grid-point matrix dimension
    NumGrid_Y::Any # Horizontal grid-point matrix dimension
    TimeStep::Any # Time step [s]
    HubHt::Any # Hub height [m] (should be > 0.5*GridHeight)
    AnalysisTime::Any # Length of analysis time series [s] (program will add time if necessary)
    GridHeight::Any # Grid height [m]
    GridWidth::Any # Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))
    VFlowAng::Any # Vertical mean flow (uptilt) angle [degrees]
    HFlowAng::Any # Horizontal mean flow (skew) angle [degrees]
    TurbModel::Any # Turbulence model (see Table 4 for valid codes)
    IECstandard::Any # Number of the IEC standard (61400-x, x=1,2,3) with optional 61400-1 ed. number
    IECturbc::Any # IEC turbulence characteristic ("A", "B", "C" or TI in %) or KHTEST
    IEC_WindType::Any # IEC turbulence type ("NTM", "xETM", "xEWM1", or "xEWM50" for x=class 1, 2, or 3)
    RefHt::Any # Height of the reference wind speed [m]
    URef::Any # Mean wind speed at the reference height [m/s] [must be 1-hr mean for API model]
    time::Any
    windvel::Any
    winddir::Any
    windvertvel::Any
    horizshear::Any
    pwrLawVertShear::Any
    LinVertShear::Any
    gustvel::Any
    UpflowAngle::Any
end

"""

runDLC(modelopt,designparams,path;runScript = OWENS.runOWENSWINDIO)
   
    # Input
    * `modelopt::OWENS.ModelingOptions`: see ?OWENS.ModelingOption
    * `designparams::OWENS.Design_Data`: see ?OWENS.Design_Data
    * `path`: desired path to run everything
    * `runScript`: function handle to run script, defaults to OWENS.runOWENSWINDIO


    # Output
    * `nothing`: 
"""
function runDLC(modelopt, designparams, path; runScript = OWENS.runOWENSWINDIO)

    if isa(designparams, String)
        designparams = Design_Data(designparams)
    end

    if isa(modelopt, String)
        modelopt = ModelingOptions(modelopt)
    end

    DLCs = modelopt.DLC_Options.DLCs
    Vinf_range = modelopt.DLC_Options.Vinf_range # = LinRange(5,20,16),
    IEC_std = modelopt.DLC_Options.IEC_std # = "\"1-ED3\"",
    WindChar = modelopt.DLC_Options.WindChar # = "\"A\"",
    WindClass = modelopt.DLC_Options.WindClass # = 1,
    turbsimsavepath = modelopt.DLC_Options.turbsimsavepath # = "./turbsimfiles",
    pathtoturbsim = modelopt.DLC_Options.pathtoturbsim # = nothing,
    NumGrid_Z = modelopt.DLC_Options.NumGrid_Z # = 38,
    NumGrid_Y = modelopt.DLC_Options.NumGrid_Y # = 26,
    Vref = modelopt.DLC_Options.Vref # = 10.0,
    Vdesign = modelopt.DLC_Options.Vdesign # = 11.0, # Design or rated speed
    grid_oversize = modelopt.DLC_Options.grid_oversize # = 1.1,
    regenWindFiles = modelopt.DLC_Options.regenWindFiles # = false,
    delta_t_turbsim = modelopt.DLC_Options.delta_t_turbsim # = 0.05,
    simtime_turbsim = modelopt.DLC_Options.simtime_turbsim # = 600.0,
    RandSeed1 = modelopt.DLC_Options.RandSeed1

    if !isdir(turbsimsavepath)
        mkdir(turbsimsavepath)
    end

    # Fill in DLC parameters based on model inputs
    DLCParams = Array{DLC_internal,1}(undef, length(DLCs))

    for (iDLC, DLC) in enumerate(DLCs) #TODO parallelize this

        DLCParams[iDLC] = getDLCparams(
            DLC,
            modelopt,
            designparams,
            Vinf_range,
            Vdesign,
            Vref,
            WindChar,
            WindClass,
            IEC_std;
            grid_oversize,
            simtime_turbsim,
            delta_t_turbsim,
            NumGrid_Z,
            NumGrid_Y,
            RandSeed1,
        )

        # Run Simulation at each Wind Speed
        for windspeed in DLCParams[iDLC].Vinf_range_used #TODO: parallelize this

            DLCParams[iDLC].URef = windspeed
            # Check if turbulent inflow file exists, if not create it
            windspeedStr = round(windspeed; digits = 2)
            windspeedStr = lpad(windspeedStr, 4, "0")
            println("Running DLC $DLC at Vinf $windspeedStr m/s")
            windINPfilename = "$turbsimsavepath/DLC$(DLC)Vinf$(windspeedStr).inp"

            if contains(DLCParams[iDLC].IEC_WindType, "NTM") ||
               contains(DLCParams[iDLC].IEC_WindType, "ETM") ||
               contains(DLCParams[iDLC].IEC_WindType, "EWM")
                if !isfile(windINPfilename) || regenWindFiles
                    generateTurbsimBTS(DLCParams[iDLC], windINPfilename, pathtoturbsim)
                end
                modelopt.OWENSOpenFASTWrappers_Options.WindType = 3
                modelopt.OWENSOpenFASTWrappers_Options.windINPfilename = "$(windINPfilename[1:end-4]).bts"
            else
                if !isfile(windINPfilename) || regenWindFiles
                    generateUniformwind(DLCParams[iDLC], windINPfilename)
                end
                modelopt.OWENSOpenFASTWrappers_Options.windINPfilename = windINPfilename
                modelopt.OWENSOpenFASTWrappers_Options.WindType = 2
            end

            modelopt.OWENSAero_Options.ifw = true
            modelopt.OWENS_Options.controlStrategy = DLCParams[iDLC].controlStrategy
            modelopt.DLC_Options.DLCParams = DLCParams[iDLC]
            # run owens simulation
            runScript(modelopt, designparams, path)
        end
    end
end

function getDLCparams(
    DLC_case,
    modelopt,
    designparams,
    Vinf_range,
    Vdesign,
    Vref,
    WindChar,
    WindClass,
    IEC_std;
    RandSeed1 = 40071,
    grid_oversize = 1.2,
    simtime_turbsim = nothing,
    delta_t_turbsim = nothing,
    NumGrid_Z = nothing,
    NumGrid_Y = nothing,
)

    hub_height = designparams[:assembly][:hub_height]
    blade_x =
        designparams[:components][:blade][:outer_shape_bem][:reference_axis][:x][:values] #Used
    blade_y =
        designparams[:components][:blade][:outer_shape_bem][:reference_axis][:y][:values] #Used
    blade_z =
        designparams[:components][:blade][:outer_shape_bem][:reference_axis][:z][:values] #Used
    Blade_Height = maximum(blade_z) #TODO: resolve DLC dependence
    Blade_Radius = maximum(sqrt.(blade_x .^ 2 .+ blade_y .^ 2))

    Htwr_base = hub_height-Blade_Height/2

    Ve50 = 50.0 #TODO change by class etc
    Ve1 = 30.0 #TODO

    numTS = modelopt.OWENS_Options.numTS
    delta_t = modelopt.OWENS_Options.delta_t
    simtime = numTS*delta_t

    GridHeight = Blade_Height * grid_oversize
    GridWidth = Blade_Radius * 2.0 * grid_oversize
    HubHt = Htwr_base+Blade_Height/2

    if !isnothing(NumGrid_Z)
        NumGrid_Z = NumGrid_Z
        NumGrid_Y = NumGrid_Y
    else
        NumGrid_Z = modelopt.Mesh_Options.ntelem+modelopt.Mesh_Options.nbelem
        NumGrid_Y = modelopt.Mesh_Options.nbelem
    end

    if !isnothing(simtime_turbsim)
        AnalysisTime = simtime_turbsim
    else
        AnalysisTime = simtime
    end

    VFlowAng = 0.0
    HFlowAng = 0.0

    IECstandard = IEC_std
    IECturbc = WindChar
    TurbModel = "\"IECKAI\""

    RefHt = round(HubHt)
    URef = 0.0 #gets filled in later from the Vinf_range when the .bst is generated

    if !isnothing(delta_t_turbsim)
        TimeStep = delta_t_turbsim
    else
        TimeStep = delta_t
    end

    time = LinRange(0, 10, 10)
    windvel = nothing # gets supersceded   
    winddir = nothing
    windvertvel = nothing
    horizshear = nothing
    pwrLawVertShear = nothing
    LinVertShear = nothing
    gustvel = nothing
    UpflowAngle = nothing

    if contains(IEC_std, "1-")
        if DLC_case == "1_1" || DLC_case == "1_2"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""

        elseif DLC_case == "1_3"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ETM\""

        elseif DLC_case == "1_4"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign-2.0, Vdesign+2.0]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0, 10, 15, 20, 25, 30, 10000.0]#LinRange(0,10,10)
            winddir = [0, 0, 45, 90, 45, 0, 0.0]
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = zeros(length(time))
            LinVertShear = zeros(length(time))
            gustvel = [0, 0, 7.0, 15, 7.5, 0.0, 0.0]
            UpflowAngle = zeros(length(time))

        elseif DLC_case == "1_5"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWS\""

            time = [0, 10, 15, 20, 25, 30, 10000.0]#LinRange(0,10,10)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = [0, 0, 5.0, 0, 0, 0, 0]#ones(length(time)).*10.0   
            pwrLawVertShear = zeros(length(time))
            LinVertShear = [0, 0, 0, 0, 5.0, 0, 0]#ones(length(time)).*10.0 
            gustvel = zeros(length(time))
            UpflowAngle = zeros(length(time))


        elseif DLC_case == "2_1" || DLC_case == "2_2" || DLC_case == "2_4"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "2_3"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used =
                [collect(LinRange(Vdesign-2.0, Vdesign+2.0, 2)); Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0, 30, 100)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = zeros(length(time))
            LinVertShear = zeros(length(time))
            UpflowAngle = zeros(length(time))

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel =
                simpleGustVel.(time, time_delay, G_amp, gustT) .+
                simpleGustVel.(time, time_delay2, G_amp, gustT)

        elseif DLC_case == "3_1"
            ControlStrategy = "startup"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0, 30, 10)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = ones(length(time)) .* 0.2
            LinVertShear = zeros(length(time))
            gustvel = zeros(length(time))
            UpflowAngle = zeros(length(time))

        elseif DLC_case == "3_2"
            ControlStrategy = "startup"
            Vinf_range_used = [
                Vinf_range[1];
                collect(LinRange(Vdesign-2.0, Vdesign+2.0, 2));
                Vinf_range[end]
            ]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0, 30, 100)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = zeros(length(time))
            LinVertShear = zeros(length(time))
            UpflowAngle = zeros(length(time))

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel =
                simpleGustVel.(time, time_delay, G_amp, gustT) .+
                simpleGustVel.(time, time_delay2, G_amp, gustT)

        elseif DLC_case == "3_3"
            ControlStrategy = "startup"
            Vinf_range_used = [
                Vinf_range[1];
                collect(LinRange(Vdesign-2.0, Vdesign+2.0, 2));
                Vinf_range[end]
            ]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0, 10, 15, 20, 25, 30, 10000.0]#LinRange(0,10,10)
            winddir = [0, 0, 45, 90, 45, 0, 0.0]
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = zeros(length(time))
            LinVertShear = zeros(length(time))
            gustvel = [0, 0, 7.0, 15, 7.5, 0.0, 0.0]
            UpflowAngle = zeros(length(time))

        elseif DLC_case == "4_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0, 30, 10)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = ones(length(time)) .* 0.2
            LinVertShear = zeros(length(time))
            gustvel = zeros(length(time))
            UpflowAngle = zeros(length(time))

        elseif DLC_case == "4_2"
            ControlStrategy = "shutdown"
            Vinf_range_used =
                [collect(LinRange(Vdesign-2.0, Vdesign+2.0, 2)); Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0, 30, 100)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = zeros(length(time))
            LinVertShear = zeros(length(time))
            UpflowAngle = zeros(length(time))

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel =
                simpleGustVel.(time, time_delay, G_amp, gustT) .+
                simpleGustVel.(time, time_delay2, G_amp, gustT)


        elseif DLC_case == "5_1"
            ControlStrategy = "emergencyshutdown"
            Vinf_range_used =
                [collect(LinRange(Vdesign-2.0, Vdesign+2.0, 2)); Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""

        elseif DLC_case == "6_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM50\""

        elseif DLC_case == "6_2"
            ControlStrategy = "parked_idle"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM50\""

        elseif DLC_case == "6_3"
            ControlStrategy = "parked_yaw"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""

        elseif DLC_case == "6_4"
            ControlStrategy = "parked"
            Vinf_range_used = [0.7*Ve50]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""

        elseif DLC_case == "7_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""

        elseif DLC_case == "8_1" #Startup
            ControlStrategy = "transport"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""

        elseif DLC_case == "CPCurve"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0, 30, 10)
            winddir = zeros(length(time))
            windvertvel = zeros(length(time))
            horizshear = zeros(length(time))
            pwrLawVertShear = ones(length(time)) .* 0.2
            LinVertShear = zeros(length(time))
            gustvel = zeros(length(time))
            UpflowAngle = zeros(length(time))

        else
            error("IEC61400_1 DLC_cases such as 1_1, 1_2 defined, you requested $DLC_case")
        end

    elseif contains(IEC_std, "2")
        error("IEC61400_2 DLC_cases are not fully defined")
        if DLC_case == "1_1"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "1_2"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)ECD\""


        elseif DLC_case == "1_3"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG50\""


        elseif DLC_case == "1_4"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD50\""


        elseif DLC_case == "1_5"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECG\""


        elseif DLC_case == "2_1"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NWP\""


        elseif DLC_case == "2_2"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "2_3"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG1\""


        elseif DLC_case == "3_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "3_2"
            ControlStrategy = "shutdown"
            Vinf_range_used = [Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG1\""


        elseif DLC_case == "4_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "5_1"
            ControlStrategy = "freewheelatIdle"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""


        elseif DLC_case == "5_2"
            ControlStrategy = "idle"
            Vinf_range_used = [Vdesign]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""


        elseif DLC_case == "6_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""

        elseif DLC_case == "8_1" #Startup
            ControlStrategy = "startup"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""

        else
            error(
                "IEC61400_2 DLC_cases [1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,3.1,3.2,4.1,5.1,5.2,6.1] defined, you requested $DLC",
            )
        end
    else
        error("IEC_std 61400 1-ED3 and 2 defined, you requested $IEC_std")
    end

    return DLC_internal(
        Vinf_range_used,
        analysis_type, # array of windspeeds m/s
        ControlStrategy, # "constRPM", function handle
        RandSeed1, # Turbulent Random Seed Number
        NumGrid_Z, # Vertical grid-point matrix dimension
        NumGrid_Y, # Horizontal grid-point matrix dimension
        TimeStep, # Turbsim time step [s]
        HubHt, # Hub height [m] (should be > 0.5*GridHeight)
        AnalysisTime, # Length of analysis time series [s] (program will add time if necessary)
        GridHeight, # Grid height [m]
        GridWidth, # Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))
        VFlowAng, # Vertical mean flow (uptilt) angle [degrees]
        HFlowAng, # Horizontal mean flow (skew) angle [degrees]
        TurbModel, # Turbulence model (see Table 4 for valid codes)
        IECstandard, # Number of the IEC standard (61400-x, x=1,2,3) with optional 61400-1 ed. number
        IECturbc, # IEC turbulence characteristic ("A", "B", "C" or TI in %) or KHTEST
        IEC_WindType, # IEC turbulence type ("NTM", "xETM", "xEWM1", or "xEWM50" for x=class 1, 2, or 3)
        RefHt, # Height of the reference wind speed [m]
        URef, # Mean wind speed at the reference height [m/s] [must be 1-hr mean for API model]
        time,
        windvel,
        winddir,
        windvertvel,
        horizshear,
        pwrLawVertShear,
        LinVertShear,
        gustvel,
        UpflowAngle,
    )
end

function generateUniformwind(DLCParams, windINPfilename)

    time = DLCParams.time
    windvel = ones(length(DLCParams.time)) .* DLCParams.URef
    winddir = DLCParams.winddir
    windvertvel = DLCParams.windvertvel
    horizShear = DLCParams.horizshear
    pwrLawVertShear = DLCParams.pwrLawVertShear
    LinVertShear = DLCParams.LinVertShear
    gustvel = DLCParams.gustvel
    UpflowAngle = DLCParams.UpflowAngle

    lines = [
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
    ]

    for itime = 1:length(time)
        lines = [
            lines;
            "$(time[itime]) $(windvel[itime]) $(winddir[itime]) $(windvertvel[itime]) $(horizShear[itime]) $(pwrLawVertShear[itime]) $(LinVertShear[itime]) $(gustvel[itime]) $(UpflowAngle[itime])"
        ]
    end

    # Write the new file
    open(windINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end
end

function generateTurbsimBTS(
    DLCParams,
    windINPfilename,
    pathtoturbsim;
    templatefile = "$module_path/template_files/templateTurbSim.inp",
)

    lines = readlines(templatefile)

    for fieldname in fieldnames(typeof(DLCParams))
        turbsimKeyName = String(fieldname)
        myvalue = getfield(DLCParams, fieldname)
        if turbsimKeyName != "Vinf_range_used" ||
           turbsimKeyName != "analysis_type" ||
           turbsimKeyName != "ControlStrategy"# || other strings
            for (iline, line) in enumerate(lines)
                if contains(line, " - ") #TODO: this assumes that the keys aren't in the comments
                    linenocomments, comments = split(line, " - ")
                    if contains(linenocomments, turbsimKeyName)
                        value, descriptor = split(linenocomments)
                        newline = "$myvalue $turbsimKeyName - $comments"
                        lines[iline] = newline
                        break
                    end
                end
            end
        end
    end

    # Write the new file
    open(windINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end

    if isnothing(pathtoturbsim) || pathtoturbsim=="nothing"
        run(`$(OWENSOpenFASTWrappers.turbsim()) $windINPfilename`)
    else
        run(`$pathtoturbsim $windINPfilename`)
    end
end

"""
* `time::TF`: in seconds
* `nominalVinf::TF`: Nominal velocity used to calculate the IEC gust size (m/s)
* `R::TF`: Turbine Radius (m)
* `G_amp::TF`: IEC gust amplitude (m/s)
* `gustT::TF`: IEC gust duration (s)
* `gustDelayT::TF`: IEC gust delay time
"""
function getGustVel(time, nominalVinf, R, G_amp, gustT, gustDelayT)
    ele_x = 0.0 #TODO: I don't think inflowwind takes in account the 3D nature of a vawt

    gustT = gustT * nominalVinf / R
    tr = time .- ele_x .- gustDelayT / R
    if (tr >= 0) && (tr<=gustT)
        IECGustFactor =
            1.0 - 0.37 * G_amp/nominalVinf * sin(3*pi*tr/gustT) * (1.0 - cos(2*pi*tr/gustT))
        return nominalVinf*IECGustFactor
    else
        return nominalVinf
    end

end

function simpleGustVel(time, time_delay, G_amp, gustT)
    timeused = time - time_delay
    if (timeused >= 0) && (timeused<=gustT)
        gustV = -0.37 * G_amp * sin(3*pi*timeused/gustT) * (1.0 - cos(2*pi*timeused/gustT))
    else
        gustV = 0.0
    end
    return gustV
end
