using OrderedCollections: OrderedDict
using YAML

export MasterInput, OWENS_Options, DLC_Options, OWENSAero_Options, OWENSFEA_Options, 
       OWENSOpenFASTWrappers_Options, Mesh_Options, Drivetrain_Options, Unified_Options, 
       ModelingOptions, Design_Data

"""
    MasterInput

Configuration structure containing all input parameters for OWENS analysis.

# Fields
* `analysisType::String`: Type of analysis to perform ("unsteady", "steady", "modal")
* `turbineType::String`: Type of turbine ("Darrieus", "H-VAWT", "ARCUS")
* `eta::Float64`: Blade mount point ratio (0.5 = blade half chord perpendicular to axis of rotation)
* `Nbld::Int`: Number of blades
* `towerHeight::Float64`: Tower extension height below blades in meters
* `rho::Float64`: Air density in kg/m³
* `Vinf::Float64`: Inflow wind speed in m/s
* `controlStrategy::String`: Control strategy type
* `RPM::Float64`: Rotor speed in RPM
* `Nslices::Int`: Number of VAWTAero discretizations
* `ntheta::Int`: Number of VAWTAero azimuthal discretizations
* `structuralModel::String`: Structural model type ("GX", "TNB", "ROM")
* `ntelem::Int`: Number of tower elements
* `nbelem::Int`: Number of blade elements
* `ncelem::Int`: Number of central cable elements (for ARCUS)
* `nselem::Int`: Number of strut elements
* `AeroModel::String`: Aerodynamic model type
* `ifw::Bool`: Inflow wind flag
* `WindType::Int`: Wind type specification
* `windINPfilename::String`: Path to wind input file
* `ifw_libfile::String`: Path to inflow wind library
* `adi_lib::String`: Path to aerodyn library
* `adi_rootname::String`: Root name for aerodyn files
* `Blade_Height::Float64`: Blade height in meters
* `Blade_Radius::Float64`: Blade radius in meters
* `numTS::Int`: Number of timesteps
* `delta_t::Float64`: Timestep size in seconds
* `NuMad_geom_xlscsv_file_twr::String`: Path to tower geometry file
* `NuMad_mat_xlscsv_file_twr::String`: Path to tower material file
* `NuMad_geom_xlscsv_file_bld::String`: Path to blade geometry file
* `NuMad_mat_xlscsv_file_bld::String`: Path to blade material file
* `NuMad_geom_xlscsv_file_strut::String`: Path to strut geometry file
* `NuMad_mat_xlscsv_file_strut::String`: Path to strut material file

# Constructor
```julia
MasterInput(;
    analysisType = "unsteady",
    turbineType = "Darrieus",
    eta = 0.5,
    Nbld = 3,
    towerHeight = 3.0,
    rho = 1.225,
    Vinf = 17.2,
    controlStrategy = "constantRPM",
    RPM = 17.2,
    Nslices = 30,
    ntheta = 30,
    structuralModel = "GX",
    ntelem = 10,
    nbelem = 60,
    ncelem = 10,
    nselem = 5,
    AeroModel = "AD",
    ifw = false,
    WindType = 1,
    ifw_libfile = "./../openfast/build/modules/inflowwind/libifw_c_binding",
    adi_lib = "./../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    numTS = 100,
    delta_t = 0.01,
    windINPfilename = "path/to/wind/file.bts",
    NuMad_geom_xlscsv_file_twr = "none",
    NuMad_mat_xlscsv_file_twr = "none",
    NuMad_geom_xlscsv_file_bld = "none",
    NuMad_mat_xlscsv_file_bld = "none",
    NuMad_geom_xlscsv_file_strut = "none",
    NuMad_mat_xlscsv_file_strut = "none"
)
```
"""
mutable struct MasterInput
    analysisType
    turbineType
    eta
    Nbld #WindIO
    towerHeight
    rho
    Vinf
    controlStrategy
    RPM
    Nslices
    ntheta
    structuralModel
    ntelem
    nbelem
    ncelem
    nselem
    AeroModel
    ifw
    WindType
    windINPfilename
    ifw_libfile
    adi_lib
    adi_rootname
    Blade_Height
    Blade_Radius
    numTS
    delta_t
    NuMad_geom_xlscsv_file_twr
    NuMad_mat_xlscsv_file_twr
    NuMad_geom_xlscsv_file_bld
    NuMad_mat_xlscsv_file_bld
    NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut
end

function MasterInput(;
    analysisType =  "unsteady", # unsteady, steady, modal
    turbineType =  "Darrieus", #Darrieus, H-VAWT, ARCUS
    eta =  0.5, # blade mount point ratio, 0.5 is the blade half chord is perpendicular with the axis of rotation, 0.25 is the quarter chord, etc
    Nbld =  3, # number of blade
    Blade_Height = 54.01123056,
    Blade_Radius = 110.1829092,
    towerHeight =  3.0, # m tower extension height below blades
    rho =  1.225, # air density
    Vinf =  17.2, # m/s
    controlStrategy = "constantRPM", # TODO: incorporate the others
    RPM =  17.2, #RPM
    Nslices =  30, # number of VAWTAero discritizations 
    ntheta =  30, # number of VAWTAero azimuthal discretizations
    structuralModel = "GX", #GX, TNB, ROM
    ntelem =  10, #tower elements in each 
    nbelem =  60, #blade elements in each 
    ncelem =  10, #central cable elements in each if turbineType is ARCUS
    nselem =  5, #strut elements in each if turbineType has struts
    AeroModel = "AD",
    ifw = false,
    WindType = 1,
    ifw_libfile = "./../openfast/build/modules/inflowwind/libifw_c_binding",
    adi_lib = "./../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    numTS = 100,
    delta_t = 0.01,
    windINPfilename ="$module_path/../test/data/turbsim/115mx115m_30x30_20.0msETM.bts",
    NuMad_geom_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv",
    NuMad_mat_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
    NuMad_mat_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_Struts.csv",
    NuMad_mat_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv"
    )

    return MasterInput(analysisType,turbineType,eta,Nbld,towerHeight,rho,Vinf,controlStrategy,
    RPM,Nslices,ntheta,structuralModel,ntelem,nbelem,ncelem,nselem,AeroModel,ifw,WindType,windINPfilename,ifw_libfile,adi_lib,adi_rootname,
    Blade_Height,Blade_Radius,numTS,delta_t,NuMad_geom_xlscsv_file_twr,NuMad_mat_xlscsv_file_twr,
    NuMad_geom_xlscsv_file_bld,NuMad_mat_xlscsv_file_bld,NuMad_geom_xlscsv_file_strut,NuMad_mat_xlscsv_file_strut)
end

"""

OWENS_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `analysisType`: default "Unsteady",  Unsteady, DLC, Campbell, todo: steady, flutter may be re-activated in the future.
    * `AeroModel`: default "DMS",  OWENSAero model "DMS" for double multiple streamtube or "AC" for actuator cylinder, or "AD" for aerodyn
    * `structuralModel`: default "TNB",  Structural models available: TNB full timoshenko beam elements with time newmark beta time stepping, ROM reduced order modal model of the timoshenko elements, GX with GXBeam's methods for geometrically exact beam theory and more efficient methods and time stepping
    * `controlStrategy`: default "normal",  should be in WindIO?- yes, 
    * `numTS`: default 10,  number of time steps TODO: change to sim time and make this derived
    * `delta_t`: default 0.05,  time step in seconds
    * `platformActive`: default false,  flag to indicate if the floating platform model is active.  
    * `topsideOn`: default true,  flat to be able to turn off the rotor and just run the floating portions
    * `interpOrder`: default 2,  if platformActive, order used for extrapolating inputs and states, 0 flat, 1 linear, 2 quadratic
    * `dataOutputFilename`: default nothing,  data output filename with path, set to nothing or don't specify to not output anything
    * `rigid`: default false,  this bypasses the structural solve and just mapps the applied loads as the reaction loads, and the deflections remain 0
    * `TOL`: default 1e-4,  gauss-seidel iteration tolerance - i.e. the two-way iteration tolerance
    * `MAXITER`: default 300,  gauss-seidel max iterations - i.e. the two-way iterations
    * `verbosity`: default 2,  verbosity where 0 is nothing, 1 is warnings, 2 is summary outputs, 3 is detailed outputs, and 4 is everything
    * `VTKsaveName`: default "./vtk/windio",  Path and name of the VTK outputs, recommended to put it in its own folder (which it will automatically create if needed)
    * `aeroLoadsOn`: default 2,  Level of aero coupling 0 structures only, 1 no deformation passed to the aero, 2 two-way coupling, 1.5 last time step's deformations passed to this timesteps aero and no internal iteration.
    * `Prescribed_RPM_time_controlpoints`: default [0.0,100000.1],  If controlStrategy is "fixedRPM", array of time control points for the internal spline
    * `Prescribed_RPM_RPM_controlpoints`: default [17.2,17.2],  If controlStrategy is "fixedRPM", array of RPM control points for the internal spline
    * `Prescribed_Vinf_time_controlpoints`: default [0.0,100000.1],  If AeroModel is "DMS" or "AC, and ifw is false, array of time control points for the internal spline
    * `Prescribed_Vinf_Vinf_controlpoints`: default [17.2,17.2],  If AeroModel is "DMS" or "AC, and ifw is false, array of Vinf control points for the internal spline 


    # Output
    * `OWENS_Options`: 
"""
mutable struct OWENS_Options
    analysisType
    AeroModel
    structuralModel
    controlStrategy
    numTS
    delta_t
    platformActive
    topsideOn
    interpOrder
    dataOutputFilename
    rigid
    TOL
    MAXITER
    verbosity
    VTKsaveName
    aeroLoadsOn
    Prescribed_RPM_time_controlpoints
    Prescribed_RPM_RPM_controlpoints
    Prescribed_Vinf_time_controlpoints
    Prescribed_Vinf_Vinf_controlpoints

    # Constructor that takes a dictionary
    function OWENS_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(
            get(dict_in,:analysisType, "Unsteady"), # Unsteady, DLC, Campbell, todo: steady, flutter may be re-activated in the future.
            get(dict_in,:AeroModel, "DMS"), # OWENSAero model "DMS" for double multiple streamtube or "AC" for actuator cylinder, or "AD" for aerodyn
            get(dict_in,:structuralModel, "TNB"), # Structural models available: TNB full timoshenko beam elements with time newmark beta time stepping, ROM reduced order modal model of the timoshenko elements, GX with GXBeam's methods for geometrically exact beam theory and more efficient methods and time stepping
            get(dict_in,:controlStrategy, "normal"), # should be in WindIO?- yes, 
            get(dict_in,:numTS, 10), # number of time steps TODO: change to sim time and make this derived
            get(dict_in,:delta_t, 0.05), # time step in seconds
            get(dict_in,:platformActive, false), # flag to indicate if the floating platform model is active.  
            get(dict_in,:topsideOn, true), # flat to be able to turn off the rotor and just run the floating portions
            get(dict_in,:interpOrder, 2), # if platformActive, order used for extrapolating inputs and states, 0 flat, 1 linear, 2 quadratic
            get(dict_in,:dataOutputFilename, "./default_savename"), # data output filename with path, set to nothing or don't specify to not output anything
            get(dict_in,:rigid, false), # this bypasses the structural solve and just mapps the applied loads as the reaction loads, and the deflections remain 0
            get(dict_in,:TOL, 1e-4), # gauss-seidel iteration tolerance - i.e. the two-way iteration tolerance
            get(dict_in,:MAXITER, 300), # gauss-seidel max iterations - i.e. the two-way iterations
            get(dict_in,:verbosity, 2), # verbosity where 0 is nothing, 1 is warnings, 2 is summary outputs, 3 is detailed outputs, and 4 is everything
            get(dict_in,:VTKsaveName, "./vtk/windio"), # Path and name of the VTK outputs, recommended to put it in its own folder (which it will automatically create if needed)
            get(dict_in,:aeroLoadsOn, 2), # Level of aero coupling 0 structures only, 1 no deformation passed to the aero, 2 two-way coupling, 1.5 last time step's deformations passed to this timesteps aero and no internal iteration.
            get(dict_in,:Prescribed_RPM_time_controlpoints, [0.0,100000.1]), # If controlStrategy is "fixedRPM", array of time control points for the internal spline
            get(dict_in,:Prescribed_RPM_RPM_controlpoints, [17.2,17.2]), # If controlStrategy is "fixedRPM", array of RPM control points for the internal spline
            get(dict_in,:Prescribed_Vinf_time_controlpoints, [0.0,100000.1]), # If AeroModel is "DMS" or "AC, and ifw is false, array of time control points for the internal spline
            get(dict_in,:Prescribed_Vinf_Vinf_controlpoints, [17.2,17.2]), # If AeroModel is "DMS" or "AC, and ifw is false, array of Vinf control points for the internal spline 
        )
    end
end

"""

DLC_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `DLCs`: default ["none"], name of DLC
    * `Vinf_range`: default inRange(5,20,16), inflow Cutin to cutout and discretization
    * `IEC_std`: default "\"1-ED3\"", turbsim input file IEC standard
    * `WindChar`: default "\"A\"", turbsim wind charasteric 
    * `WindClass`: default , DLC turbsim wind class
    * `turbsimsavepath`: default "./turbsimfiles", path where the turbsim files are saved
    * `pathtoturbsim`: default othing, path to the turbsim executable
    * `NumGrid_Z`: default 8, turbsim vertical discretizations 
    * `NumGrid_Y`: default 6, turbsim horizontal discretizations
    * `Vref`: default 0.0, reference/nominal wind speed m/s for turbsim or other inflow wind input file (depending on which DLC is selected)
    * `Vdesign`: default 1.0, Design or rated speed of turbine, used for certain DLC cases
    * `grid_oversize`: default .1, amount that the turbsim inflow is oversized compared to the turbine to allow for deflection
    * `regenWindFiles`: default false, force regeneration of turbsim files even if they already exist
    * `delta_t_turbsim`: default .05, turbsim timestep
    * `simtime_turbsim`: default 00.0, turbsim total time, which loops if simtime exceeds turbsim time
    * `RandSeed1`: default 0071, turbsim random seed number
    * `DLCParams`: see ?OWENS.DLCParams, the current DLC parameters for the run, used as internal state information


    # Output
    * `DLC_Options`: 
"""
mutable struct DLC_Options
    DLCs
    Vinf_range
    IEC_std
    WindChar
    WindClass
    turbsimsavepath
    pathtoturbsim
    NumGrid_Z
    NumGrid_Y
    Vref
    Vdesign
    grid_oversize
    regenWindFiles
    delta_t_turbsim
    simtime_turbsim
    RandSeed1
    DLCParams

    # Constructor that takes a dictionary
    function DLC_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(
            get(dict_in,:DLCs,["none"]), # name of DLC
            get(dict_in,:Vinf_range,LinRange(5,20,16)), # inflow Cutin to cutout and discretization
            get(dict_in,:IEC_std,"\"1-ED3\""), # turbsim input file IEC standard
            get(dict_in,:WindChar,"\"A\""), # turbsim wind charasteric 
            get(dict_in,:WindClass,1), # DLC turbsim wind class
            get(dict_in,:turbsimsavepath,"./turbsimfiles"), # path where the turbsim files are saved
            get(dict_in,:pathtoturbsim,nothing), # path to the turbsim executable
            get(dict_in,:NumGrid_Z,38), # turbsim vertical discretizations 
            get(dict_in,:NumGrid_Y,26), # turbsim horizontal discretizations
            get(dict_in,:Vref,10.0), # reference/nominal wind speed m/s for turbsim or other inflow wind input file (depending on which DLC is selected)
            get(dict_in,:Vdesign,11.0), # Design or rated speed of turbine, used for certain DLC cases
            get(dict_in,:grid_oversize,1.1), # amount that the turbsim inflow is oversized compared to the turbine to allow for deflection
            get(dict_in,:regenWindFiles,false), #, force regeneration of turbsim files even if they already exist
            get(dict_in,:delta_t_turbsim,0.05), # turbsim timestep
            get(dict_in,:simtime_turbsim,600.0), # turbsim total time, which loops if simtime exceeds turbsim time
            get(dict_in,:RandSeed1,40071), # turbsim random seed number
            get(dict_in,:DLCParams,nothing), # must be filled in with the DLC generator
            
        )
    end
end
    
"""

OWENSAero_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `Nslices`: default 20, number of 3-D slices for the strip method to go from 2D to 3D considering curved deforming blades
    * `ntheta`: default 30, number of azimuthal discretizations
    * `ifw`: default false, use the inflow wind coupling to get inflow velocities TODO: change ifw to inflowwind inflowwind_active etc everywhere
    * `DynamicStallModel`: default "BV", dynamic stall model, should be under an OWENSAero options
    * `RPI`: default true, rotating point iterative method (i.e. it just calculates at the blade positions and is much faster)
    * `Aero_Buoyancy_Active`: default false, flag to turn buoyancy on for the blades.  This is likely to be replaced by a different model
    * `Aero_AddedMass_Active`: default false, flag to turn added mass forces on, don't turn on if the added mass in the structures are on
    * `Aero_RotAccel_Active`: default false, flag to turn added mass forces on, don't turn on if the added mass in the structures are on
    * `eta`: default 0.5, blade mount point ratio, 0.5 is the blade half chord is perpendicular with the axis of rotation, 0.25 is the quarter chord, etc
    * `rho`: default 1.225, air density in kg/m³
    * `Vinf`: default 17.2, inflow wind speed in m/s
    * `RPM`: default 17.2, rotor speed in RPM

    # Output
    * `OWENSAero_Options`: 
"""
mutable struct OWENSAero_Options #TODO: move these downstream to their respective packages and unify the options with those
    Nslices
    ntheta
    ifw
    DynamicStallModel
    RPI
    Aero_Buoyancy_Active
    Aero_AddedMass_Active
    Aero_RotAccel_Active
    eta
    rho
    Vinf
    RPM
    
    # Constructor that takes a dictionary
    function OWENSAero_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(
            get(dict_in,:Nslices, 20), # number of 3-D slices for the strip method to go from 2D to 3D considering curved deforming blades
            get(dict_in,:ntheta, 30), # number of azimuthal discretizations
            get(dict_in,:ifw, false), # use the inflow wind coupling to get inflow velocities TODO: change ifw to inflowwind inflowwind_active etc everywhere
            get(dict_in,:DynamicStallModel,"BV"), # dynamic stall model, should be under an OWENSAero options
            get(dict_in,:RPI, true), # rotating point iterative method (i.e. it just calculates at the blade positions and is much faster)
            get(dict_in,:Aero_Buoyancy_Active, false), # flag to turn buoyancy on for the blades.  This is likely to be replaced by a different model
            get(dict_in,:Aero_AddedMass_Active, false), # flag to turn added mass forces on, don't turn on if the added mass in the structures are on
            get(dict_in,:Aero_RotAccel_Active, false), # flag to turn added mass forces on, don't turn on if the added mass in the structures are on
            get(dict_in,:eta, 0.5), # blade mount point ratio, 0.5 is the blade half chord is perpendicular with the axis of rotation, 0.25 is the quarter chord, etc
            get(dict_in,:rho, 1.225), # air density in kg/m³
            get(dict_in,:Vinf, 17.2), # inflow wind speed in m/s
            get(dict_in,:RPM, 17.2), # rotor speed in RPM
        )
    end
end

"""

OWENSFEA_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `nlOn`: default true, nonlinear effects
    * `RayleighAlpha`: default 0.05, damping coefficient scalar on the stiffness matrix
    * `RayleighBeta`: default 0.05, damping coefficient scalar on the mass matrix
    * `iterationType`: default "DI", internal iteration type DI direct iteration, NR newton rhapson (which is less stable than DI)
    * `guessFreq`: default 0.0, for the built in flutter model frequency guessed for the flutter frequency 
    * `numModes`: default 20, ROM model, number of modes used in the analysis type.  Less is faster but less accurate
    * `adaptiveLoadSteppingFlag`: default true, for steady analysis if convergence fails, it will reduce the load and retry then increase the load
    * `minLoadStepDelta`: default 0.0500, minimum change in load step
    * `minLoadStep`: default 0.0500, minimum value of reduced load
    * `prescribedLoadStep`: default 0.0, optional prescribed load fraction
    * `maxNumLoadSteps`: default 20, used in static (steady state) analysis, max load steps for adaptive load stepping
    * `tolerance`: default 1.0000e-06, total mesh unsteady analysis convergence tolerance for a timestep within the structural model
    * `maxIterations`: default 50, total mesh unsteady analysis convergence max iterations for a timestep
    * `elementOrder`: default 1, Element order, 1st order, 2nd order etc; determines the number of nodes per element (order +1).  Orders above 1 have not been tested in a long time  
    * `alpha`: default 0.5, newmark time integration alpha parameter
    * `gamma`: default 0.5, newmark time integration gamma parameter
    * `AddedMass_Coeff_Ca`: default 0.0, added mass coefficient, scaling factor (typically 0-1) on the cones of water mass applied to each structural element in the 22 and 33 diagonal terms. 0 turns this off
    * `platformTurbineConnectionNodeNumber`: default 1, TODO: reconnect this
    * `aeroElasticOn`: default false, OWENSFEA for the built in flutter model
    * `spinUpOn`: default true, TODO: remove this since it should always be true since that is how its used. To turn it off, just set RPM and gravity to 0.  OWENSFEA modal analysis, calculates steady centrifugal strain stiffening and then passes that model to the modal analysis
    * `predef`: default false, Predeformation flag for two state analysis where a portion of the blade is deformed and the nonlinear strain stiffening terms are "update"-d, then "use"-d in two different analysis

    # Output
    * `OWENSFEA_Options`: 
"""
mutable struct OWENSFEA_Options
    nlOn
    RayleighAlpha
    RayleighBeta
    iterationType
    guessFreq
    numModes
    adaptiveLoadSteppingFlag
    minLoadStepDelta
    minLoadStep
    prescribedLoadStep
    maxNumLoadSteps
    tolerance
    maxIterations
    elementOrder
    alpha
    gamma
    AddedMass_Coeff_Ca
    platformTurbineConnectionNodeNumber
    aeroElasticOn
    spinUpOn
    predef

    # Constructor that takes a dictionary
    function OWENSFEA_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(       
            get(dict_in,:nlOn, true), # nonlinear effects
            get(dict_in,:RayleighAlpha, 0.05), # damping coefficient scalar on the stiffness matrix
            get(dict_in,:RayleighBeta, 0.05), # damping coefficient scalar on the mass matrix
            get(dict_in,:iterationType, "DI"), # internal iteration type DI direct iteration, NR newton rhapson (which is less stable than DI)
            get(dict_in,:guessFreq, 0.0), # for the built in flutter model frequency guessed for the flutter frequency 
            get(dict_in,:numModes, 20), # ROM model, number of modes used in the analysis type.  Less is faster but less accurate
            get(dict_in,:adaptiveLoadSteppingFlag, true), # for steady analysis if convergence fails, it will reduce the load and retry then increase the load
            get(dict_in,:minLoadStepDelta, 0.0500), # minimum change in load step
            get(dict_in,:minLoadStep, 0.0500), # minimum value of reduced load
            get(dict_in,:prescribedLoadStep, 0.0), # optional prescribed load fraction
            get(dict_in,:maxNumLoadSteps, 20), # used in static (steady state) analysis, max load steps for adaptive load stepping
            get(dict_in,:tolerance, 1.0000e-06), # total mesh unsteady analysis convergence tolerance for a timestep within the structural model
            get(dict_in,:maxIterations, 50), # total mesh unsteady analysis convergence max iterations for a timestep
            get(dict_in,:elementOrder, 1), # Element order, 1st order, 2nd order etc; determines the number of nodes per element (order +1).  Orders above 1 have not been tested in a long time  
            get(dict_in,:alpha, 0.5), # newmark time integration alpha parameter
            get(dict_in,:gamma, 0.5), # newmark time integration gamma parameter
            get(dict_in,:AddedMass_Coeff_Ca, 0.0), # added mass coefficient, scaling factor (typically 0-1) on the cones of water mass applied to each structural element in the 22 and 33 diagonal terms. 0 turns this off
            get(dict_in,:platformTurbineConnectionNodeNumber, 1), # TODO: reconnect this
            get(dict_in,:aeroElasticOn, false), # OWENSFEA for the built in flutter model
            get(dict_in,:spinUpOn, true), # TODO: remove this since it should always be true since that is how its used. To turn it off, just set RPM and gravity to 0.  OWENSFEA modal analysis, calculates steady centrifugal strain stiffening and then passes that model to the modal analysis
            get(dict_in,:predef, false), # Predeformation flag for two state analysis where a portion of the blade is deformed and the nonlinear strain stiffening terms are "update"-d, then "use"-d in two different analysis
        )
    end
end

"""

OWENSOpenFASTWrappers_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `windINPfilename`: default nothing, If ifw or AeroDyn is being used, gets overwritten if using the DLC analysis type, the moordyn file location, like in the unit test
    * `ifw_libfile`: default nothing, location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
    * `hd_lib`: default nothing, location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
    * `md_lib`: default nothing, location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
    * `adi_lib`: default nothing, location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
    * `adi_rootname`: default "/aerodyn", location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
    * `hd_input_file`: default "none", If platformActive, the hydrodyn file location, like in the unit test
    * `ss_input_file`: default "none", If platformActive, the sea state file location, like in the unit test
    * `md_input_file`: default "none", If platformActive, the moordyn file location, like in the unit test
    * `potflowfile`: default nothing, If platformActive, the potential flow files location, like in the unit test
    * `WindType`: default 3, Derived parameter, inflowwind wind file type when DLC generator is active, matches inflowwind WindType 
            
    # Output
    * `OWENSOpenFASTWrappers_Options`: 
"""
mutable struct OWENSOpenFASTWrappers_Options
    windINPfilename
    ifw_libfile
    hd_lib
    md_lib
    adi_lib
    adi_rootname
    hd_input_file
    ss_input_file
    md_input_file
    potflowfile
    WindType

    # Constructor that takes a dictionary
    function OWENSOpenFASTWrappers_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(       
            get(dict_in,:windINPfilename, nothing), # If ifw or AeroDyn is being used, gets overwritten if using the DLC analysis type, the moordyn file location, like in the unit test
            get(dict_in,:ifw_libfile, nothing), # location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
            get(dict_in,:hd_lib, nothing),# location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
            get(dict_in,:md_lib, nothing),# location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
            get(dict_in,:adi_lib, nothing),# location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
            get(dict_in,:adi_rootname, "/aerodyn"),# location of the respective OpenFAST library, if nothing it will use the internal OWENS installation
            get(dict_in,:hd_input_file, "none"), # If platformActive, the hydrodyn file location, like in the unit test
            get(dict_in,:ss_input_file, "none"), # If platformActive, the sea state file location, like in the unit test
            get(dict_in,:md_input_file, "none"), # If platformActive, the moordyn file location, like in the unit test
            get(dict_in,:potflowfile, nothing),# If platformActive, the potential flow files location, like in the unit test
            get(dict_in,:WindType, 3),#Derived parameter, inflowwind wind file type when DLC generator is active, matches inflowwind WindType 
        )
    end
end


"""
Mesh_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `ntelem`: default 20, number of tower elements in each blade, plus nodes wherever there is a component overlap
    * `nbelem`: default 30, number of blade elements in each blade, plus nodes wherever there is a component overlap
    * `ncelem`: default 30, number of cable elements in each cable if ARCUS
    * `nselem`: default 10, number of elements in each strut
    * `angularOffset`: default 0.0, moves the structure to align with the aero model
    * `joint_type`: default 0, optionally can specify the strut to blade joints to be pinned about different axes, or 0 for welded
    * `c_mount_ratio`: default 0.05, for ARCUS, where the cable mounts on the lower side of the blade
    * `AD15hubR`: default 0.1, parameter, used in aerodyn coupling for the hub radius so that the vortex sheets don't go within the hub
    * `cables_connected_to_blade_base`: default true, for ARCUS, for the two part simulation of the blade bending
    * `turbineType`: default "Darrieus", mesh Darrieus, H-VAWT, controls if the tips of the blades are joined to the tower in the mesh or not.
    * `Nbld`: default 3, number of blades
    * `Blade_Height`: default 54.01123056, blade height in meters
    * `Blade_Radius`: default 110.1829092, blade radius in meters
    * `towerHeight`: default 3.0, tower extension height below blades in meters
    * `NuMad_geom_xlscsv_file_twr`: default "none", path to tower geometry file
    * `NuMad_mat_xlscsv_file_twr`: default "none", path to tower material file
    * `NuMad_geom_xlscsv_file_bld`: default "none", path to blade geometry file
    * `NuMad_mat_xlscsv_file_bld`: default "none", path to blade material file
    * `NuMad_geom_xlscsv_file_strut`: default "none", path to strut geometry file
    * `NuMad_mat_xlscsv_file_strut`: default "none", path to strut material file
            
    # Output
    * `Mesh_Options`: 
"""
mutable struct Mesh_Options
    ntelem
    nbelem
    ncelem
    nselem
    angularOffset
    joint_type
    c_mount_ratio
    AD15hubR
    cables_connected_to_blade_base
    turbineType
    Nbld
    Blade_Height
    Blade_Radius
    towerHeight
    NuMad_geom_xlscsv_file_twr
    NuMad_mat_xlscsv_file_twr
    NuMad_geom_xlscsv_file_bld
    NuMad_mat_xlscsv_file_bld
    NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut
    
    # Constructor that takes a dictionary
    function Mesh_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(       
            get(dict_in,:ntelem, 20), # number of tower elements in each blade, plus nodes wherever there is a component overlap
            get(dict_in,:nbelem, 30), # number of blade elements in each blade, plus nodes wherever there is a component overlap
            get(dict_in,:ncelem, 30), # number of cable elements in each cable if ARCUS
            get(dict_in,:nselem, 10), # number of elements in each strut
            get(dict_in,:angularOffset, 0.0), # moves the structure to align with the aero model
            get(dict_in,:joint_type, 0), # optionally can specify the strut to blade joints to be pinned about different axes, or 0 for welded
            get(dict_in,:c_mount_ratio, 0.05), # for ARCUS, where the cable mounts on the lower side of the blade
            get(dict_in,:AD15hubR, 0.1), # parameter, used in aerodyn coupling for the hub radius so that the vortex sheets don't go within the hub
            get(dict_in,:cables_connected_to_blade_base, true), # for ARCUS, for the two part simulation of the blade bending
            get(dict_in,:turbineType, "Darrieus"), #mesh Darrieus, H-VAWT, controls if the tips of the blades are joined to the tower in the mesh or not.
            get(dict_in,:Nbld, 3), # number of blades
            get(dict_in,:Blade_Height, 54.01123056), # blade height in meters
            get(dict_in,:Blade_Radius, 110.1829092), # blade radius in meters
            get(dict_in,:towerHeight, 3.0), # tower extension height below blades in meters
            get(dict_in,:NuMad_geom_xlscsv_file_twr, "none"), # path to tower geometry file
            get(dict_in,:NuMad_mat_xlscsv_file_twr, "none"), # path to tower material file
            get(dict_in,:NuMad_geom_xlscsv_file_bld, "none"), # path to blade geometry file
            get(dict_in,:NuMad_mat_xlscsv_file_bld, "none"), # path to blade material file
            get(dict_in,:NuMad_geom_xlscsv_file_strut, "none"), # path to strut geometry file
            get(dict_in,:NuMad_mat_xlscsv_file_strut, "none"), # path to strut material file
        )
    end
end


"""

Drivetrain_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
   
    # Input
    * `turbineStartup`: default 0, TODO: clean up since it should be derived from control strategy
    * `usingRotorSpeedFunction`: default false, TODO: clean up the speed function since the omegaocp RPM gets splined already
    * `driveTrainOn`: default false, flag to turn on the drivetrain model TODO: clean this up to make it always use the drivetrain model, with default 100% efficiency and ratio of 1 so it outputs the values
    * `JgearBox`: default 0.0, torsional stiffness of the gearbox TODO: resolve units
    * `gearRatio`: default 1.0, ratio between the turbine driveshaft and generator shaft
    * `gearBoxEfficiency`: default 1.0, efficiency of the gearbox, just decreases the torque that the generator model sees
    * `generatorOn`: default false, TODO: clean up the generator options
    * `useGeneratorFunction`: default false, TODO: clean up the generator options
    * `generatorProps`: default 0.0, TODO: clean up the generator options
    * `ratedTorque`: default 0.0, TODO: clean up the generator options
    * `zeroTorqueGenSpeed`: default 0.0, TODO: clean up the generator options
    * `pulloutRatio`: default 0.0, TODO: clean up the generator options
    * `ratedGenSlipPerc`: default 0.0, TODO: clean up the generator options
    * `OmegaGenStart`: default 0.0, TODO: clean up the generator options
    * `driveShaftProps_K`: default 0.0, TODO: break this out, driveshaft stiffness and damping
    * `driveShaftProps_C`: default 0.0, TODO: break this out, driveshaft stiffness and damping

    # Output
    * `Drivetrain_Options`: 
"""
mutable struct Drivetrain_Options
    turbineStartup
    usingRotorSpeedFunction
    driveTrainOn
    JgearBox
    gearRatio
    gearBoxEfficiency
    generatorOn
    useGeneratorFunction
    generatorProps
    ratedTorque
    zeroTorqueGenSpeed
    pulloutRatio
    ratedGenSlipPerc
    OmegaGenStart
    driveShaft_K
    driveShaft_C

    # Constructor that takes a dictionary
    function Drivetrain_Options(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Use get to provide default values for missing fields
        new(       
            get(dict_in,:turbineStartup, 0), # TODO: clean up since it should be derived from control strategy
            get(dict_in,:usingRotorSpeedFunction, false), #TODO: clean up the speed function since the omegaocp RPM gets splined already
            get(dict_in,:driveTrainOn, false), #flag to turn on the drivetrain model #TODO: clean this up to make it always use the drivetrain model, with default 100% efficiency and ratio of 1 so it outputs the values
            get(dict_in,:JgearBox, 0.0), # torsional stiffness of the gearbox TODO: resolve units
            get(dict_in,:gearRatio, 1.0), # ratio between the turbine driveshaft and generator shaft
            get(dict_in,:gearBoxEfficiency, 1.0), # efficiency of the gearbox, just decreases the torque that the generator model sees
            get(dict_in,:generatorOn, false), #TODO: clean up the generator options
            get(dict_in,:useGeneratorFunction, false), #TODO: clean up the generator options
            get(dict_in,:generatorProps, 0.0), #TODO: clean up the generator options
            get(dict_in,:ratedTorque, 0.0), #TODO: clean up the generator options
            get(dict_in,:zeroTorqueGenSpeed, 0.0), #TODO: clean up the generator options
            get(dict_in,:pulloutRatio, 0.0), #TODO: clean up the generator options
            get(dict_in,:ratedGenSlipPerc, 0.0), #TODO: clean up the generator options
            get(dict_in,:OmegaGenStart, 0.0), #TODO: clean up the generator options
            get(dict_in,:driveShaftProps_K, 0.0), #TODO: break this out, driveshaft stiffness and damping
            get(dict_in,:driveShaftProps_C, 0.0), #TODO: break this out, driveshaft stiffness and damping
        )
    end
end

"""

Unified_Options
   
    # Input
    * `OWENS_Options::OWENS_Options`:
    * `DLC_Options::DLC_Options`:
    * `OWENSAero_Options::OWENSAero_Options`:
    * `OWENSFEA_Options::OWENSFEA_Options`:
    * `OWENSOpenFASTWrappers_Options::OWENSOpenFASTWrappers_Options`:
    * `Mesh_Options::Mesh_Options`:
    * `Drivetrain_Options::Drivetrain_Options`:

    # Output
    * `Unified_Options`: 
"""
mutable struct Unified_Options
    OWENS_Options::OWENS_Options
    DLC_Options::DLC_Options
    OWENSAero_Options::OWENSAero_Options
    OWENSFEA_Options::OWENSFEA_Options
    OWENSOpenFASTWrappers_Options::OWENSOpenFASTWrappers_Options
    Mesh_Options::Mesh_Options
    Drivetrain_Options::Drivetrain_Options
end

"""

ModelingOptions(yamlInputfile)
   
    # Input
    * `yamlInputfile::string`: yaml file containing ordered inputs matching the default keys


    # Output
    * `Unified_Options::Unified_Options`: Struct of structs containing all of the OWENS Options
"""
function ModelingOptions(yamlInputfile=nothing)

    if !isnothing(yamlInputfile)
        yamlInput = YAML.load_file(yamlInputfile;dicttype=OrderedCollections.OrderedDict{Symbol,Any})
    else #just use defaults by supplying a dummy dictionary up front
        yamlInput = OrderedCollections.OrderedDict(:nothing=>0.0,:nothing2=>"string")
    end
    
    # Unpack YAML
    dummy_dict = OrderedCollections.OrderedDict(:nothing=>0.0,:nothing2=>"string")

    if haskey(yamlInput,:DLC_Options)
        dlc_options = DLC_Options(yamlInput[:DLC_Options])
    else
        dlc_options = DLC_Options(dummy_dict)
    end

    if haskey(yamlInput,:OWENSAero_Options)
        owensaero_options = OWENSAero_Options(yamlInput[:OWENSAero_Options])
    else
        owensaero_options = OWENSAero_Options(dummy_dict)
    end

    if haskey(yamlInput,:OWENS_Options)
        owens_options = OWENS_Options(yamlInput[:OWENS_Options])
    else
        owens_options = OWENS_Options(dummy_dict)
    end
    
    if haskey(yamlInput,:OWENSFEA_Options)
        owensfea_options = OWENSFEA_Options(yamlInput[:OWENSFEA_Options])
    else
        owensfea_options = OWENSFEA_Options(dummy_dict)
    end

    if haskey(yamlInput,:OWENSOpenFASTWrappers_Options)
        owensopenfastwrappers_options = OWENSOpenFASTWrappers_Options(yamlInput[:OWENSOpenFASTWrappers_Options])
    else
        owensopenfastwrappers_options = OWENSOpenFASTWrappers_Options(dummy_dict)
    end

    if haskey(yamlInput,:Mesh_Options)
        mesh_options = Mesh_Options(yamlInput[:Mesh_Options])
    else
        mesh_options = Mesh_Options(dummy_dict)
    end

    if haskey(yamlInput,:Drivetrain_Options)
        drivetrain_options = Drivetrain_Options(yamlInput[:Drivetrain_Options])
    else
        drivetrain_options = Drivetrain_Options(dummy_dict)
    end

    return Unified_Options(owens_options,dlc_options,owensaero_options,owensfea_options,owensopenfastwrappers_options,mesh_options,drivetrain_options)
end

"""
    convertMasterInputToUnifiedOptions(masterInput::MasterInput)

Converts a MasterInput struct to Unified_Options by mapping the fields appropriately.

# Arguments
- `masterInput::MasterInput`: The legacy MasterInput struct

# Returns
- `Unified_Options`: The new unified options structure
"""
function convertMasterInputToUnifiedOptions(masterInput::MasterInput)
    # Create default option structs
    dummy_dict = OrderedCollections.OrderedDict(:nothing=>0.0,:nothing2=>"string")
    
    # Create OWENS_Options with values from MasterInput
    owens_options_dict = OrderedCollections.OrderedDict{Symbol,Any}(
        :analysisType => masterInput.analysisType,
        :AeroModel => masterInput.AeroModel,
        :structuralModel => masterInput.structuralModel,
        :controlStrategy => masterInput.controlStrategy,
        :numTS => masterInput.numTS,
        :delta_t => masterInput.delta_t,
        :Prescribed_RPM_RPM_controlpoints => [masterInput.RPM],
        :Prescribed_Vinf_Vinf_controlpoints => [masterInput.Vinf]
    )
    owens_options = OWENS_Options(owens_options_dict)
    
    # Create OWENSAero_Options with values from MasterInput
    aero_options_dict = OrderedCollections.OrderedDict{Symbol,Any}(
        :Nslices => masterInput.Nslices,
        :ntheta => masterInput.ntheta,
        :ifw => masterInput.ifw,
        :eta => masterInput.eta,
        :rho => masterInput.rho,
        :Vinf => masterInput.Vinf,
        :RPM => masterInput.RPM
    )
    owensaero_options = OWENSAero_Options(aero_options_dict)
    
    # Create Mesh_Options with values from MasterInput
    mesh_options_dict = OrderedCollections.OrderedDict{Symbol,Any}(
        :ntelem => masterInput.ntelem,
        :nbelem => masterInput.nbelem,
        :ncelem => masterInput.ncelem,
        :nselem => masterInput.nselem,
        :turbineType => masterInput.turbineType,
        :Nbld => masterInput.Nbld,
        :Blade_Height => masterInput.Blade_Height,
        :Blade_Radius => masterInput.Blade_Radius,
        :towerHeight => masterInput.towerHeight,
        :NuMad_geom_xlscsv_file_twr => masterInput.NuMad_geom_xlscsv_file_twr,
        :NuMad_mat_xlscsv_file_twr => masterInput.NuMad_mat_xlscsv_file_twr,
        :NuMad_geom_xlscsv_file_bld => masterInput.NuMad_geom_xlscsv_file_bld,
        :NuMad_mat_xlscsv_file_bld => masterInput.NuMad_mat_xlscsv_file_bld,
        :NuMad_geom_xlscsv_file_strut => masterInput.NuMad_geom_xlscsv_file_strut,
        :NuMad_mat_xlscsv_file_strut => masterInput.NuMad_mat_xlscsv_file_strut
    )
    mesh_options = Mesh_Options(mesh_options_dict)
    
    # Create OpenFAST wrapper options
    openfast_options_dict = OrderedCollections.OrderedDict{Symbol,Any}(
        :windINPfilename => masterInput.windINPfilename,
        :ifw_libfile => masterInput.ifw_libfile,
        :adi_lib => masterInput.adi_lib,
        :adi_rootname => masterInput.adi_rootname,
        :WindType => masterInput.WindType
    )
    openfast_options = OWENSOpenFASTWrappers_Options(openfast_options_dict)
    
    # Create other options with defaults
    dlc_options = DLC_Options(dummy_dict)
    fea_options = OWENSFEA_Options(dummy_dict)
    drivetrain_options = Drivetrain_Options(dummy_dict)
    
    # Create and return the unified options
    unified_options = Unified_Options(owens_options, dlc_options, owensaero_options, 
                                    fea_options, openfast_options, mesh_options, drivetrain_options)
    
    return unified_options
end

function Design_Data(file_path=nothing; design_defaults_yaml="$(module_path)/template_files/design_defaults.yml")
    # Load the YAML files
    if !isnothing(file_path)
        windio = YAML.load_file(file_path; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
        println("Running: $(windio[:name])")
    else
        windio = OrderedCollections.OrderedDict(:nothing=>0.0,:nothing2=>"string")
    end

    defaults = YAML.load_file(design_defaults_yaml; dicttype=OrderedCollections.OrderedDict{Symbol,Any})

    # Create a new dictionary that merges loaded data with defaults
    Design_Data = OrderedCollections.OrderedDict{Symbol, Any}()

    # Fill in the Design_Data with defaults and loaded values
    for (key, default_value) in defaults
        Design_Data[key] = get(windio, key, default_value)
    end

    return Design_Data
end

function MasterInput(yamlInputfile)

    @warn "The old combined yaml file is being depreciated in favor of the windio yaml and the modelingoptions format"

    yamlInput = YAML.load_file(yamlInputfile)
    # Unpack YAML
    general = yamlInput["general"]
        analysisType = general["analysisType"]
        turbineType = general["turbineType"]

    designParameters = yamlInput["designParameters"]
        eta = designParameters["eta"]
        Nbld = designParameters["Nbld"]
        Blade_Height = designParameters["Blade_Height"]
        Blade_Radius = designParameters["Blade_Radius"]
        towerHeight = designParameters["towerHeight"]

    operationParameters = yamlInput["operationParameters"]
        rho = operationParameters["rho"]
        Vinf = operationParameters["Vinf"]

    controlParameters = yamlInput["controlParameters"]
        controlStrategy = controlParameters["controlStrategy"]
        RPM = controlParameters["RPM"]
        numTS = controlParameters["numTS"]
        delta_t = controlParameters["delta_t"]

    AeroParameters = yamlInput["AeroParameters"]
        Nslices = AeroParameters["Nslices"]
        ntheta = AeroParameters["ntheta"]
        AeroModel = AeroParameters["AeroModel"]
        adi_lib = AeroParameters["adi_lib"]
        adi_rootname = AeroParameters["adi_rootname"]

    turbulentInflow = yamlInput["turbulentInflow"]
        ifw = turbulentInflow["ifw"]
        WindType = turbulentInflow["WindType"]
        windINPfilename = turbulentInflow["windINPfilename"]
        ifw_libfile = turbulentInflow["ifw_libfile"]

    structuralParameters = yamlInput["structuralParameters"]
        structuralModel = structuralParameters["structuralModel"]
        ntelem = structuralParameters["ntelem"]
        nbelem = structuralParameters["nbelem"]
        ncelem = structuralParameters["ncelem"]
        nselem = structuralParameters["nselem"]
        NuMad_geom_xlscsv_file_twr = structuralParameters["NuMad_geom_xlscsv_file_twr"]
        NuMad_mat_xlscsv_file_twr = structuralParameters["NuMad_mat_xlscsv_file_twr"]
        NuMad_geom_xlscsv_file_bld = structuralParameters["NuMad_geom_xlscsv_file_bld"]
        NuMad_mat_xlscsv_file_bld = structuralParameters["NuMad_mat_xlscsv_file_bld"]
        NuMad_geom_xlscsv_file_strut = structuralParameters["NuMad_geom_xlscsv_file_strut"]
        NuMad_mat_xlscsv_file_strut = structuralParameters["NuMad_mat_xlscsv_file_strut"]

    return MasterInput(analysisType,turbineType,eta,Nbld,towerHeight,rho,Vinf,
    controlStrategy,RPM,Nslices,ntheta,structuralModel,ntelem,nbelem,ncelem,
    nselem,AeroModel,ifw,WindType,windINPfilename,ifw_libfile,adi_lib,adi_rootname,Blade_Height,Blade_Radius,numTS,
    delta_t,NuMad_geom_xlscsv_file_twr,NuMad_mat_xlscsv_file_twr,
    NuMad_geom_xlscsv_file_bld,NuMad_mat_xlscsv_file_bld,NuMad_geom_xlscsv_file_strut,NuMad_mat_xlscsv_file_strut)
end
