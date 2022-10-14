struct Bin
"""
Inputs pointing to the file paths of compiled binaries of external libraries
"""
    hydrodynLibPath
    moordynLibPath
end

mutable struct Inputs
    analysisType
    turbineStartup
    usingRotorSpeedFunction
    tocp
    tocp_Vinf
    numTS
    delta_t
    Omegaocp
    Vinfocp
    driveTrainOn
    generatorOn
    aeroLoadsOn
    hydroOn
    topsideOn
    interpOrder
    hd_input_file
    md_input_file
    JgearBox
    gearRatio
    gearBoxEfficiency
    useGeneratorFunction
    generatorProps
    ratedTorque
    zeroTorqueGenSpeed
    pulloutRatio
    ratedGenSlipPerc
    OmegaGenStart
    omegaControl
    OmegaInit
    rigid
    aeroloadfile
    owensfile
    potflowfile
    outFilename
    bladeData
    driveShaftProps
end

# this way you can use defaults and pass in what is different, and it's mapped
# by keyword so it doesn't have to be in order.
"""
Inputs(;analysisType = "TNB",
    turbineStartup = 0,
    usingRotorSpeedFunction = false,
    tocp = [0.0,1.1],
    tocp_Vinf = [0.0,1.1],
    numTS = 50.0,
    delta_t = 2e-3,
    Omegaocp = [7.2,7.2] ./ 60,
    Vinfocp = [12.0,12.0],
    aeroLoadsOn = 1,
    driveTrainOn = false,
    generatorOn = false,
    hydroOn = false,
    topsideOn = true,
    interpOrder = 2,
    hd_input_file = "none",
    md_input_file = "none",
    JgearBox = 0.0,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    ratedTorque = 0.0,
    zeroTorqueGenSpeed = 0.0,
    pulloutRatio = 0.0,
    ratedGenSlipPerc = 0.0,
    OmegaGenStart = 0.0,
    omegaControl = false,
    OmegaInit = 7.2/60, #TODO: simplify this in the code since it is redundant
    rigid = false, #turn off structural dynamics
    aeroloadfile = "module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv",
    owensfile = "module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens",
    outFilename = "none",
    numDofPerNode = 6,
    bladeData = [],
    driveShaftProps = DriveShaftProps(0.0,0.0)
    )

Model inputs for OWENS coupled analysis, struct

# Inputs
* `analysisType::string`: Newmark Beta time stepping "TNB", Dean time stepping "TD", modal "M"
* `turbineStartup::int`: 1 forced start-up using generator as motor, 2 self-starting mode, 0 specified rotor speed mode")
* `usingRotorSpeedFunction::bool`: use user specified rotor speed profile function
* `tocp::Array{<:float}`: = time points for rotor speed profile (s)
* `tocp_Vinf::Array{<:float}`: = time points for specified Vinf profile (s)
* `numTS::int`: total number of timesteps to run
* `delta_t::float`: timestep interval (s)
* `Omegaocp::Array{<:float}`: = rotor speed points for rotor speed profile (Hz)
* `Vinfocp::Array{<:float}`: = rotor speed points for specified Vinf profile (Hz)
* `aeroLoadsOn::bool`: flag to trigger aero loads being applied
* `driveTrainOn::bool`: flag to include drivetrain effects
* `generatorOn::bool`: flag to include generator effects
* `hydroOn::bool`: flag to include platform coupling
* `interpOrder::int`: order used for extrapolating inputs and states, 0 flat, 1 linear, 2 quadratic
* `hd_input_file::string`: file path to the HydroDyn .dat input file
* `md_input_file::string`: file path to the MoorDyn .dat input file
* `JgearBox::float`: gearbox intertia, standard SI units
* `gearRatio::float`: gearbox gear ratio
* `gearBoxEfficiency::float`: gearbox efficiency (typically 0-1)
* `useGeneratorFunction::bool`: = flag to use user specified generator profile
* `generatorProps::float`: not used, should clean up
* `ratedTorque::float`: Generator rated max torque
* `zeroTorqueGenSpeed::float`: rated generator speed (minus slippage)
* `pulloutRatio::float`: Fraction of the min/max torque that the generator engages/disengages
* `ratedGenSlipPerc::float`: extra speed from slipping?
* `OmegaGenStart::float`: speed (Hz) at which generator would kick in
* `omegaControl::bool`: false for fixed speed, true for dynamic
* `OmegaInit::float`: initial rotor speed (Hz)
* `aeroloadfile::string`: string of the name and path for the cactus aeroloads if using the old serial owens call
* `owensfile::string`: string of the name and path for the owens input file if using the old serial owens call
* `potflowfile::string`: string of the prefix and path for the directory containing the potential flow files from WAMIT (required by HydroDyn)
* `outFilename::string`: path and name of output file, will be overwritten if already exists
* `numDofPerNode::int`: number of degrees of freedom per node
* `bladeData::BladeData`: see ?BladeData, only used if calling the old serial owens function
* `driveShaftProps::DriveShaftProps`: see ?DriveShaftProps

# Outputs:
* `none`:
"""

function Inputs(;analysisType = "TNB",
    turbineStartup = 0,
    usingRotorSpeedFunction = false,
    tocp = [0.0,1.1],
    tocp_Vinf = [0.0,1e6],
    numTS = 50.0,
    delta_t = 2e-3,
    Omegaocp = [7.2,7.2] ./ 60,
    Vinfocp = [10.0,10.0],
    driveTrainOn = false,
    generatorOn = false,
    aeroLoadsOn = false, #this need to get cleaned up in the code
    hydroOn = false,
    topsideOn = true,
    interpOrder = 2,
    hd_input_file = "none",
    md_input_file = "none",
    JgearBox = 0.0,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    ratedTorque = 0.0,
    zeroTorqueGenSpeed = 0.0,
    pulloutRatio = 0.0,
    ratedGenSlipPerc = 0.0,
    OmegaGenStart = 0.0,
    omegaControl = false,
    OmegaInit = 7.2/60, #TODO: simplify this in the code since it is redundant
    aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv",
    owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens",
    potflowfile = "$module_path/../test/data/potential_flow_data",
    outFilename = "none",
    numDofPerNode = 6,
    bladeData = [],
    rigid = false,
    driveShaftProps = DriveShaftProps(0.0,0.0)
    )

    return Inputs(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,tocp_Vinf,numTS,delta_t,Omegaocp,Vinfocp,
    driveTrainOn,generatorOn,aeroLoadsOn,hydroOn,topsideOn,interpOrder,hd_input_file,md_input_file,
    JgearBox,gearRatio,gearBoxEfficiency,useGeneratorFunction,generatorProps,ratedTorque,
    zeroTorqueGenSpeed,pulloutRatio,ratedGenSlipPerc,OmegaGenStart,omegaControl,OmegaInit,rigid,
    aeroloadfile,owensfile,potflowfile,outFilename,bladeData,driveShaftProps)
end

"""
Internal, driveshaft stiffness k and damping c
"""
mutable struct DriveShaftProps
    k
    c
end

# Cactus Related Structs
"""
Internal, struct containing the CACTUS geometry file data
"""
mutable struct CactusGeom
    NBlade
    NStrut
    RotN
    RotP
    RefAR
    RefR
    blade
    strut
end

"""
Internal, struct containing the CACTUS geometry file data for a blade
"""
mutable struct Blade
    NElem
    FlipN
    QCx
    QCy
    QCz
    tx
    ty
    tz
    CtoR
    PEx
    PEy
    PEz
    tEx
    tEy
    tEz
    nEx
    nEy
    nEz
    sEx
    sEy
    sEz
    ECtoR
    EAreaR
    iSect
end

"""
Internal, struct containing the CACTUS geometry file data for a strut
"""
mutable struct Strut
    NElem
    TtoC
    MCx
    MCy
    MCz
    CtoR
    PEx
    PEy
    PEz
    sEx
    sEy
    sEz
    ECtoR
    EAreaR
    BIndS
    EIndS
    BIndE
    EIndE
end

"""
Internal, struct containing blade specific data and location within the mesh
"""
mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
end
