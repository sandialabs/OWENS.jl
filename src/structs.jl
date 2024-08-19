"""
Inputs pointing to the file paths of compiled binaries of external libraries
"""
struct Bin
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
    AD15On
    hydroOn
    topsideOn
    interpOrder
    hd_input_file
    ss_input_file
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
    iteration_parameters
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
    AD15On = false,
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
    TOl = 1e-4,
    MAXITER = 300,
    iterwarnings = true,
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
* `aeroLoadsOn::bool`: #0 off, 1 one way, 1.5 one way with deformation from last timestep, 2 two way
* `AD15On::bool`: flag to use AD15 for aero
* `driveTrainOn::bool`: flag to include drivetrain effects
* `generatorOn::bool`: flag to include generator effects
* `hydroOn::bool`: flag to include platform coupling
* `interpOrder::int`: order used for extrapolating inputs and states, 0 flat, 1 linear, 2 quadratic
* `hd_input_file::string`: file path to the HydroDyn .dat input file
* `ss_input_file::string`: file path to the HydroDyn sea states input file
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
* `TOl::float`: gauss-seidel iteration tolerance
* `MAXITER::int`: gauss-seidel maximum iterations
* `iterwarnings::bool`: iteration warnings flag


# Outputs:
* `OWENS.Inputs`:
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
    AD15On = false,
    hydroOn = false,
    topsideOn = true,
    interpOrder = 2,
    hd_input_file = "none",
    ss_input_file = "none",
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
    driveShaftProps = DriveShaftProps(0.0,0.0),
    TOl = 1e-4,
    MAXITER = 300,
    iterwarnings = true,
    )

    return Inputs(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,tocp_Vinf,numTS,delta_t,Omegaocp,Vinfocp,
    driveTrainOn,generatorOn,aeroLoadsOn,AD15On,hydroOn,topsideOn,interpOrder,hd_input_file,ss_input_file,md_input_file,
    JgearBox,gearRatio,gearBoxEfficiency,useGeneratorFunction,generatorProps,ratedTorque,
    zeroTorqueGenSpeed,pulloutRatio,ratedGenSlipPerc,OmegaGenStart,omegaControl,OmegaInit,rigid,
    aeroloadfile,owensfile,potflowfile,outFilename,bladeData,driveShaftProps,Iteration_Parameters(TOl,MAXITER,iterwarnings))
end

"""
Internal, driveshaft stiffness k and damping c
"""
mutable struct DriveShaftProps
    k
    c
end

"""
Internal, gauss-seidel iteration parameters 
"""
mutable struct Iteration_Parameters            
    TOL # = 1e-4  #gauss-seidel iteration tolerance for various modules
    MAXITER # = 2 #max iteration for various modules
    iterwarnings
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


"""
NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin,web_stack,web_dp)

Parameters defining the blade composite layup. See NuMad user guide SAND2012_7028 appendix B for more details

**Arguments**
- `n_web::Int64`: number of shear webs
- `n_stack::Int64`: number of predefined composite stacks
- `n_segments::Int64`: number of segments around the airfoil
- `span::Vector{Float64}`: span-wise position
- `airfoil::Vector{String}`: airfoil name
- `te_type::Vector{String}`: trailing edge type
- `twist_d::Vector{Float64}`: twist_d in degrees
- `chord::Vector{Float64}`: chord length
- `xoffset::Vector{Float64}`: The distance from the “nose” of a station to the blade reference axis.
- `aerocenter::Vector{Float64}`: This is an aerodynamic parameter that is an output from aerodynamic performance analysis of a two-dimensional airfoil section. The aerodynamic center is the point along the chord where the aerodynamic pitching moment does not vary with changes in angle of attack.
- `stack_mat_types::Vector{Int64}`: Material numbers used that correspond to each stack number
- `stack_layers::Array{Int64,2}`: number of layers at each span used corresponding to each material type (first index corresponds to spanwise position, second index corresponds to the stack number)
- `segments::Array{Float64,2}`: normalized starting and stopping points of each section (i.e. leading edge, sparcap, etc).  First index corresponds to  spanwise position, second index corresponds to the section, except there is an extra first column starting at -1 for the trailing edge. There must be a leading edge position at 0, and the last column must be 1 corresponding to the trailing edge again.  Positions are fractions of the chord, lower (HP) is negative, upper (LP) is positive
- `DPtypes::Array{Int64,2}`: division point types (NOTE THAT THIS ISN'T IMPLEMENTED AND DOES NOTHING CURRENTLY, i.e. only SINGLE is being used). First index corresponds to  spanwise position, second corresponds to section number
- `skin_seq::Array{Seq,2}`: stack sequence, is an array of structures, each containing a Vector{Int64} of the sequence (i.e. skin[2,5].seq). First index corresponds to spanwise positoin, second index the section
- `web_seq::Array{Seq,2}`: same format and meaning as skin sequence, but for the webs with the second index corresponding to the web number
- `web_dp::Array{Seq,2}`: same format as skin sequence, but this corresponds to the section numbers the web connects to at the top and bottom at both edges. There are always four entries in the CSV list and the order goes as follows: inboard LP, inboard HP, outboard HP, outboard LP.
"""
mutable struct NuMad
    n_web
    n_stack
    n_segments
    span
    airfoil
    te_type
    twist_d
    chord
    xoffset
    aerocenter
    stack_mat_types
    stack_layers
    segments
    DPtypes
    skin_seq
    web_seq
    web_dp
end

mutable struct Seq
    seq
end

struct plyproperties
    names#::Array{String,1}
    plies#::Array{Composites.Material,1}
    costs#::Array{Float,1}
    SN_stressMpa # control points for the SN curve, matrix material x 6 where 6 is the current number of control points in an akima spline
    Log_SN_cycles2Fail # control points for the SN curve, matrix material x 6 where 6 is the current number of control points in an akima spline
end
plyproperties(names,plies) = plyproperties(names,plies,zeros(length(names)),collect(cat(fill(collect(LinRange(1e12,0,6)),length(names))[:,:]...,dims=2)'),collect(cat(fill(collect(LinRange(0,7,6)),length(names))[:,:]...,dims=2)')) #Backwards compatible convenience function

"""
Struct containing

material names

Composites.Material structs for each material name - see ?Composites.Material
"""
function plyproperties()
    #TODO: read in numad xls materials file
    names = ["highmodulus_uni",
    "standard_uni",
    "MR60H",
    "T3900_uni",
    "T700_uni",
    "ELT5500",
    "UDCarbon",
    "highmodulus_weave",
    "standard_weave",
    "T3900_weave",
    "T700_weave",
    "Gelcoat",
    "Triax",
    "Saertex",
    "taylor_foam",
    "SNL_foam",
    "aluminum6061",
    "aluminum6063",
    "CLA_5500",
    "CBX_2400",
    "ETLX_2400",
    "Airex_C70_55",
    "EBX_2400_x10",
    "ETLX_2400_x10",
    "Airex_C70_55_x10"]

    n_materials = 25
    e1  =zeros(n_materials) #pa
    e2  =zeros(n_materials) #pa
    g12 =zeros(n_materials) #pa
    anu =zeros(n_materials) #ratio
    rho =zeros(n_materials) #kg/m3
    xt  =zeros(n_materials) #pa
    xc  =zeros(n_materials) #pa
    yt  =zeros(n_materials) #pa
    yc  =zeros(n_materials) #pa
    s   =zeros(n_materials) #pa
    plythickness =zeros(n_materials) #meters

    # "highmodulus_uni"
    e1[1]  = 175.0e9
    e2[1]  = 8.0e9
    g12[1] = 5.0e9
    anu[1] = 0.30
    rho[1] = 1600.0
    xt[1]  = min(1.0,1355.656/1682.011)*1000e6 #mean -> A-basis
    xc[1]  = min(1.0,1103.943/1396.504)*850e6 #mean -> A-basis
    yt[1]  = min(1.0,39.226/52.975)*40e6 #mean -> A-basis
    yc[1]  = min(1.0,235.434/282.439)*200e6 #mean -> A-basis
    s[1]   = min(1.0,142.411/159.516)*60e6 #mean -> A-basis
    plythickness[1] = 0.152e-3
    # "standard_uni"
    e1[2]  = 135.0e9
    e2[2]  = 10.0e9
    g12[2] = 5.0e9
    anu[2] = 0.30
    rho[2] = 1600.0
    xt[2]  = min(1.0,1355.656/1682.011)*1500e6 #mean -> A-basis
    xc[2]  = min(1.0,1103.943/1396.504)*1200e6 #mean -> A-basis
    yt[2]  = min(1.0,39.226/52.975)*50e6 #mean -> A-basis
    yc[2]  = min(1.0,235.434/282.439)*250e6 #mean -> A-basis
    s[2]   = min(1.0,142.411/159.516)*70e6 #mean -> A-basis
    plythickness[2] = 0.152e-3
    # "MR60H"
    e1[3]  = (165e9+150e9)/2.0
    e2[3]  = 8.56e9
    g12[3] = 4.39e9
    anu[3] = 0.326
    rho[3] = 1810.0
    xt[3]  = min(1.0,1355.656/1682.011)*3190e6 #mean -> A-basis
    xc[3]  = min(1.0,1103.943/1396.504)*1440e6 #mean -> A-basis
    yt[3]  = min(1.0,39.226/52.975)*82.0e6 #mean -> A-basis
    yc[3]  = min(1.0,235.434/282.439)*200.0e6 #mean -> A-basis
    s[3]   = min(1.0,142.411/159.516)*141e6 #mean -> A-basis
    plythickness[3] = 0.152e-3
    # "T3900_uni"
    e1[4]  = (148e9+131e9)/2.0
    e2[4]  = (9.7e9+9.7e9)/2.0
    g12[4] = 4.83e9
    anu[4] = 0.33
    rho[4] = 1573.0
    xt[4]  = min(1.0,1355.656/1682.011)*2830e6 #CTD mean -> A-basis
    xc[4]  = min(1.0,1103.943/1396.504)*1772e6 #CTD mean -> A-basis
    yt[4]  = min(1.0,39.226/52.975)*56.9e6 #CTD mean -> A-basis
    yc[4]  = min(1.0,235.434/282.439)*303e6 #CTD mean -> A-basis
    s[4]   = min(1.0,142.411/159.516)*89.6e6 #CTD mean -> A-basis
    plythickness[4] = 0.191e-3
    # "T700_uni"
    e1[5]=120.8e9
    e2[5]=11.57e9
    g12[5]=5.219e9
    anu[5]=0.350
    rho[5]=1525.0
    xt[5]=1356e6
    xc[5]=1104e6
    yt[5]=39.23e6
    yc[5]=235.4e6
    s[5]=142.4e6
    plythickness[5] = 0.152e-3
    # "ELT5500"
    e1[6]=41.8e9
    e2[6]=14.00e9
    g12[6]=2.630e9
    anu[6]=0.280
    rho[6]=1920.0 #g/cc * 1000
    xt[6]=972.0e6
    xc[6]=702.0e6
    yt[6]=100.0e6 #made up
    yc[6]=100.0e6 #made up
    s[6]=100.0e6 #made up
    plythickness[6] = 0.47e-3
    # "UDCarbon"
    e1[7]=114.5e9
    e2[7]=8.39e9
    g12[7]=5.990e9
    anu[7]=0.270
    rho[7]=1220.0
    xt[7]=1546.0e6
    xc[7]=1047.0e6
    yt[7]=100.0e6 #made up
    yc[7]=100.0e6 #made up
    s[7]=100.0e6 #made up
    plythickness[7] = 0.47e-3

    # FABRICS
    # "highmodulus_weave"
    e1[8]  = 85.0e9
    e2[8]  = 85.0e9
    g12[8] = 5.0e9
    anu[8] = 0.10
    rho[8] = 1600.0
    xt[8]  = min(1.0,1355.656/1682.011)*350e6 #mean -> A-basis
    xc[8]  = min(1.0,1103.943/1396.504)*150e6 #mean -> A-basis
    yt[8]  = min(1.0,39.226/52.975)*350e6 #mean -> A-basis
    yc[8]  = min(1.0,235.434/282.439)*150e6 #mean -> A-basis
    s[8]   = min(1.0,142.411/159.516)*35e6 #mean -> A-basis
    plythickness[8] = 0.218e-3
    # "standard_weave"
    e1[9]  = 70.0e9
    e2[9]  = 70.0e9
    g12[9] = 5.0e9
    anu[9] = 0.10
    rho[9] = 1600.0
    xt[9]  = min(1.0,1355.656/1682.011)*600e6 #mean -> A-basis
    xc[9]  = min(1.0,1103.943/1396.504)*570e6 #mean -> A-basis
    yt[9]  = min(1.0,39.226/52.975)*600e6 #mean -> A-basis
    yc[9]  = min(1.0,235.434/282.439)*570e6 #mean -> A-basis
    s[9]   = min(1.0,142.411/159.516)*90e6 #mean -> A-basis
    plythickness[9] = 0.218e-3
    # "T3900_weave"
    e1[10]  = (70.3e9+71.0e9)/2.0
    e2[10]  = (68.9e9+67.6e9)/2.0
    g12[10] = 4.6e9
    anu[10] = 0.032
    rho[10] = 1551.0
    xt[10]  = min(1.0,701.302/803.236)*1055e6 #CTD mean -> A-basis
    xc[10]  = min(1.0,549.748/749.955)*676e6 #CTD mean -> A-basis
    yt[10]  = min(1.0,557.575/722.602)*945e6 #CTD mean -> A-basis
    yc[10]  = min(1.0,604.067/741.866)*614e6 #CTD mean -> A-basis
    s[10]   = min(1.0,138.440/154.888)*79.3e6 #CTD mean -> A-basis
    plythickness[10] = 0.218e-3
    # "T700_weave"
    e1[11]=55.82e9
    e2[11]=52.10e9
    g12[11]=4.295e9
    anu[11]=0.085
    rho[11]=1501.0
    xt[11]=701.32e6
    yt[11]=557.59e6
    xc[11]=549.77e6
    yc[11]=604.08e6
    s[11]=138.44e6
    plythickness[11] = 0.218e-3
    # "Gelcoat"
    e1[12]=3.44e9
    e2[12]=3.44e9
    g12[12]=1.38e9
    anu[12]=0.3
    rho[12]=1235.0
    xt[12]=100.0e6 #made up
    xc[12]=100.0e6 #made up
    yt[12]=100.0e6 #made up
    yc[12]=100.0e6 #made up
    s[12]=100.0e6 #made up
    plythickness[12] = 0.05e-3
    # "Triax"
    e1[13]=27.7e9
    e2[13]=13.65e9
    g12[13]=7.2e9
    anu[13]=0.39
    rho[13]=1850.0
    xt[13]=700.0e6
    xc[13]=100.0e6 #made up
    yt[13]=700.0e6
    yc[13]=100.0e6 #made up
    s[13]=100.0e6 #made up
    plythickness[13] = 0.94e-3
    # "Saertex"
    e1[14]=13.6e9
    e2[14]=13.3e9
    g12[14]=11.8e9
    anu[14]=0.49
    rho[14]=1780.0
    xt[14]=144.0e6
    xc[14]=213.0e6
    yt[14]=144.0e6
    yc[14]=213.0e6
    s[14]=100.0e6 #made up
    plythickness[14] = 1.0e-3

    #FOAMS
    # "taylor_foam"
    e1[15]=48.0e6
    e2[15]=48.0e6
    g12[15]=28.0e6
    anu[15]=0.3
    rho[15]=75.0
    xt[15]=100.0e6 #made up
    yt[15]=100.0e6 #made up
    xc[15]=100.0e6 #made up
    yc[15]=100.0e6 #made up
    s[15]=100.0e6 #made up
    plythickness[15]=1.0E-3
    # "SNL_foam"
    e1[16]=256.0e6
    e2[16]=256.0e6
    g12[16]=22.0e6
    anu[16]=0.3
    rho[16]=200.0
    xt[16]=100.0e6 #made up
    yt[16]=100.0e6 #made up
    xc[16]=100.0e6 #made up
    yc[16]=100.0e6 #made up
    s[16]=100.0e6 #made up
    plythickness[16]=1.0E-3

    # Metals
    # "aluminum6061" http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061T6
    e1[17]=68.9e9
    e2[17]=68.9e9
    g12[17]=26.0e9
    anu[17]=0.33
    rho[17]=2700.0 #g/cc * 1000
    xt[17]=276.0e6 #use yield for metal?
    yt[17]=276.0e6 #use yield for metal?
    xc[17]=386.0e6 #use yield for metal?
    yc[17]=386.0e6 #use yield for metal?
    s[17]=207.0e6
    plythickness[17]=1.0E-3

    # "aluminum6063" http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T6
    e1[18]=68.9e9
    e2[18]=68.9e9
    g12[18]=25.8e9
    anu[18]=0.33
    rho[18]=2700.0 #g/cc * 1000
    xt[18]=214.0e6 #use yield for metal?
    yt[18]=214.0e6 #use yield for metal?
    xc[18]=276.0e6 #use yield for metal?
    yc[18]=276.0e6 #use yield for metal?
    s[18]=152.0e6
    plythickness[18]=1.0E-3

    # "CLA_5500"
    e1[19]=98.240e9 #pa
    e2[19]=5.102e9 #pa
    g12[19]=4.274e9 #pa
    anu[19]=0.3 #ratio
    rho[19]=1540.0 #g/cc * 1000 #kg/m3
    xt[19]=875.634139e6 #pa
    xc[19]=592.949102e6 #pa
    yt[19]=100.0e6 #made up
    yc[19]=100.0e6 #made up
    s[19]=100.0e6 #made up
    plythickness[19]=0.66E-3 #meters

    # "CBX_2400"
    e1[20]=14.931e9 #pa
    e2[20]=14.931e9 #pa
    g12[20]=23.890e9 #pa
    anu[20]=0.3 #ratio
    rho[20]=1530.0 #g/cc * 1000 #kg/m3
    xt[20]=455.053962e6 #pa
    xc[20]=455.053962e6 #pa
    yt[20]=100.0e6 #made up
    yc[20]=100.0e6 #made up
    s[20]=100.0e6 #made up
    plythickness[20]=0.81E-3 #meters

    # "ETLX_2400"
    e1[21]=20.333e9 #pa
    e2[21]=9.305e9 #pa
    g12[21]=4.756e9 #pa
    anu[21]=0.3 #ratio
    rho[21]=1900.0 #g/cc * 1000 #kg/m3
    xt[21]=530.896289e6 #pa
    xc[21]=530.896289e6 #pa
    yt[21]=100.0e6 #made up
    yc[21]=100.0e6 #made up
    s[21]=100.0e6 #made up
    plythickness[21]=0.66E-3 #meters

    # "Airex_C70_55"
    e1[22]=0.045e9 #pa
    e2[22]=0.045e9 #pa
    g12[22]=0.022e9 #pa
    anu[22]=0.2 #ratio
    rho[22]=59.0 #g/cc * 1000 #kg/m3
    xt[22]=100.0e6 #pa #made up
    xc[22]=100.0e6 #pa
    yt[22]=100.0e6 #made up
    yc[22]=100.0e6 #made up
    s[22]=100.0e6 #made up
    plythickness[22]=1.0E-3 #meters

    # "EBX_2400_x10"
    e1[23]=982.400e9 #pa
    e2[23]=51.020e9 #pa
    g12[23]=42.740e9 #pa
    anu[23]=0.3 #ratio
    rho[23]=15300.0 #g/cc * 1000 #kg/m3
    xt[23]=4550.53962e6 #pa
    xc[23]=4550.53962e6 #pa
    yt[23]=100.0e6 #made up
    yc[23]=100.0e6 #made up
    s[23]=100.0e6 #made up
    plythickness[23]=0.07E-3 #meters

    # "ETLX_2400_x10"
    e1[24]=149.310e9 #pa
    e2[24]=149.310e9 #pa
    g12[24]=238.900e9 #pa
    anu[24]=0.3 #ratio
    rho[24]=19000.0 #g/cc * 1000 #kg/m3
    xt[24]=5308.96289e6 #pa
    xc[24]=5308.96289e6 #pa
    yt[24]=100.0e6 #made up
    yc[24]=100.0e6 #made up
    s[24]=100.0e6 #made up
    plythickness[24]=0.08E-3 #meters

    # "Airex_C70_55_x10"
    e1[25]=203.335e9 #pa
    e2[25]=93.051e9 #pa
    g12[25]=47.560e9 #pa
    anu[25]=0.2 #ratio
    rho[25]=590.0 #g/cc * 1000 #kg/m3
    xt[25]=100.0e6 #pa #made up
    xc[25]=100.0e6 #pa
    yt[25]=100.0e6 #made up
    yc[25]=100.0e6 #made up
    s[25]=100.0e6 #made up
    plythickness[25]=0.07E-3 #meters

    return plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness))
end