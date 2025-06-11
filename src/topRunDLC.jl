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
    * `DynamicStallModel`: defaultBV", dynamic stall model, should be under an OWENSAero options
    * `RPI`: default true, rotating point iterative method (i.e. it just calculates at the blade positions and is much faster)
    * `Aero_Buoyancy_Active`: default false, flag to turn buoyancy on for the blades.  This is likely to be replaced by a different model
    * `Aero_AddedMass_Active`: default false, flag to turn added mass forces on, don't turn on if the added mass in the structures are on
    * `Aero_RotAccel_Active`: default false, flag to turn added mass forces on, don't turn on if the added mass in the structures are on


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
    * `Setup_Options::SetupOptions`:

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
    Setup_Options::SetupOptions
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

    if haskey(yamlInput,:Setup_Options)
        setup_options = SetupOptions(yamlInput[:Setup_Options])
    else
        setup_options = SetupOptions()
    end

    return Unified_Options(owens_options,dlc_options,owensaero_options,owensfea_options,owensopenfastwrappers_options,mesh_options,drivetrain_options,setup_options)
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

function runOWENS(Inp,designparams,path;verbosity=2)
    analysisType = Inp.analysisType
    turbineType = Inp.turbineType
    eta = Inp.eta
    Nbld = Inp.Nbld
    towerHeight = Inp.towerHeight
    rho = Inp.rho
    Vinf = Inp.Vinf
    controlStrategy = Inp.controlStrategy
    RPM = Inp.RPM
    Nslices = Inp.Nslices
    ntheta = Inp.ntheta
    structuralModel = Inp.structuralModel
    ntelem = Inp.ntelem
    nbelem = Inp.nbelem
    ncelem = Inp.ncelem
    nselem = Inp.nselem
    ifw = Inp.ifw
    WindType = Inp.WindType
    AeroModel = Inp.AeroModel
    windINPfilename = "$(path)$(Inp.windINPfilename)"
    ifw_libfile = Inp.ifw_libfile
    if ifw_libfile == "nothing"
        ifw_libfile = nothing
    end
    Blade_Height = Inp.Blade_Height
    Blade_Radius = Inp.Blade_Radius
    numTS = Inp.numTS
    delta_t = Inp.delta_t
    NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
    NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
    NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld)"
    NuMad_mat_xlscsv_file_bld = "$(path)$(Inp.NuMad_mat_xlscsv_file_bld)"
    NuMad_geom_xlscsv_file_strut = "$(path)$(Inp.NuMad_geom_xlscsv_file_strut)"
    NuMad_mat_xlscsv_file_strut = "$(path)$(Inp.NuMad_mat_xlscsv_file_strut)"
    adi_lib = Inp.adi_lib
    if adi_lib == "nothing"
        adi_lib = nothing
    end
    adi_rootname = "$(path)$(Inp.adi_rootname)"

    B = Nbld
    R = Blade_Radius#177.2022*0.3048 #m
    H = Blade_Height#1.02*R*2 #m

    shapeZ = collect(LinRange(0,H,Nslices+1))
    shapeX = R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

    nothing

    # Call the helper function that builds the mesh, calculates the sectional properties,
    # and aligns the sectional properties to the mesh elements, 


mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    RPM,
    Vinf,
    eta,
    B,
    H,
    R,
    shapeZ,
    shapeX,
    ifw,
    WindType,
    delta_t,
    numTS,
    adi_lib,
    adi_rootname,
    windINPfilename,
    ifw_libfile,
    NuMad_geom_xlscsv_file_twr,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS_Cables.csv",
    NuMad_mat_xlscsv_file_twr,# = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv",
    NuMad_geom_xlscsv_file_bld,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS.csv",
    NuMad_mat_xlscsv_file_bld,# = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
    NuMad_geom_xlscsv_file_strut,
    NuMad_mat_xlscsv_file_strut,
    Htwr_base=towerHeight,
    ntelem, 
    nbelem, 
    ncelem,
    nselem,
    joint_type = 0,
    c_mount_ratio = 0.05,
    AeroModel, #AD, DMS, AC
    DynamicStallModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)

    nothing

    # Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
    # this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
    # includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
    # to work with OWENS inputs as well.

    nothing

    # If the sectional properties material files includes cost information, that is combined with the density 
    # to estimate the overall material cost of of materials in the blades

    if verbosity>0
        
        println("\nBlades' Mass Breakout")
        for (i,name) in enumerate(plyprops_bld.names)
            println("$name $(mass_breakout_blds[i]) kg, $(plyprops_bld.costs[i]) \$/kg: \$$(mass_breakout_blds[i]*plyprops_bld.costs[i])")
        end
        
        println("\nTower Mass Breakout")
        for (i,name) in enumerate(plyprops_twr.names)
            println("$name $(mass_breakout_twr[i]) kg, $(plyprops_twr.costs[i]) \$/kg: \$$(mass_breakout_twr[i]*plyprops_twr.costs[i])")
        end
        
        println("Total Material Cost Blades: \$$(sum(mass_breakout_blds.*plyprops_bld.costs))")
        println("Total Material Cost Tower: \$$(sum(mass_breakout_twr.*plyprops_twr.costs))")
        println("Total Material Cost: \$$(sum(mass_breakout_blds.*plyprops_bld.costs)+ sum(mass_breakout_twr.*plyprops_twr.costs))")
        
    end

    nothing

    # Here we apply the boundary conditions.  For this case, with a regular cantelever tower, the tower base node which is 
    # 1 is constrained in all 6 degrees of freedom to have a displacement of 0.  You can change this displacement to allow for things
    # like pretension, and you can apply boundary conditions to any node.

    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0]

    nothing

    # There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    # Handle the control strategy

    normalRPM = RPM
    slowRPM = 1.0 #RPM
    tocp = [0.0,10.0,30.0,1000.0]
    turbineStartup = 0
    generatorOn = false
    useGeneratorFunction = false
    throwawayTimeSteps=1#round(Int,10.0/delta_t) # 10 seconds
    if controlStrategy == "normal"
        Omegaocp = [normalRPM,normalRPM,normalRPM,normalRPM]./60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "freewheelatNormalOperatingRPM"
        Omegaocp = [normalRPM,normalRPM,normalRPM,normalRPM]./60 #hz
        turbineStartup = 2
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "startup"
        throwawayTimeSteps = 1
        tocp = [0.0,30.0,1000.0]
        Omegaocp = [slowRPM,normalRPM,normalRPM]./60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "shutdown"
        tocp = [0.0,10.0,15.0,30.0,1000.0]
        Omegaocp = [normalRPM,normalRPM,normalRPM/2,slowRPM,0.0]./60 #hz
    elseif controlStrategy == "emergencyshutdown"
        tocp = [0.0,10.0,15.0,30.0,1000.0]
        Omegaocp = [normalRPM,normalRPM,slowRPM,0.0,0.0]./60 #hz
    elseif controlStrategy == "parked"
        Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
    elseif controlStrategy == "parked_idle"
        Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
        turbineStartup = 2
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "parked_yaw"
        Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
    elseif controlStrategy == "parked"
        Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
    elseif controlStrategy == "transport"
        Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
    end

    println("controlStrategy: $controlStrategy")

    inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
    tocp,
    Omegaocp,
    tocp_Vinf = [0.0,100000.1],
    Vinfocp = [Vinf,Vinf],
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = 2,
    turbineStartup,
    generatorOn,
    useGeneratorFunction,
    OmegaInit = Omegaocp[1])

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

    feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    dataOutputFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    gravityOn = true,
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time.

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
    topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

    if AeroModel=="AD"
        OWENSOpenFASTWrappers.endTurb()
    end

    nothing

    # Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

    azi=aziHist#./aziHist*1e-6
    VTKsaveName = "$path/vtk/SNL5MW"
    OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
        epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
        FReactionHist,topFexternal_hist)

    nothing

    # This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
    # it also calculates analytical simply supported buckling safety factors

    ##########################################
    #### Ultimate Failure #####
    ##########################################

    massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
    SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
    topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
    topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,B,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

    outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)

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
function runDLC(modelopt,designparams,path;runScript = OWENS.runOWENSWINDIO)

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
    DLCParams = Array{DLC_internal, 1}(undef, length(DLCs))

    for (iDLC, DLC) in enumerate(DLCs) #TODO parallelize this

        DLCParams[iDLC] = getDLCparams(DLC, modelopt, designparams, Vinf_range, Vdesign, Vref, WindChar,WindClass, IEC_std;grid_oversize,simtime_turbsim,delta_t_turbsim,NumGrid_Z,NumGrid_Y,RandSeed1)

        # Run Simulation at each Wind Speed
        for windspeed in DLCParams[iDLC].Vinf_range_used #TODO: parallelize this

            DLCParams[iDLC].URef = windspeed
            # Check if turbulent inflow file exists, if not create it
            windspeedStr = round(windspeed;digits=2)
            windspeedStr = lpad(windspeedStr,4,"0")
            println("Running DLC $DLC at Vinf $windspeedStr m/s")
            windINPfilename = "$turbsimsavepath/DLC$(DLC)Vinf$(windspeedStr).inp"
            
            if contains(DLCParams[iDLC].IEC_WindType, "NTM") || contains(DLCParams[iDLC].IEC_WindType, "ETM") || contains(DLCParams[iDLC].IEC_WindType, "EWM")
                if !isfile(windINPfilename) || regenWindFiles
                    generateTurbsimBTS(DLCParams[iDLC],windINPfilename,pathtoturbsim)
                end
                modelopt.OWENSOpenFASTWrappers_Options.WindType = 3
                modelopt.OWENSOpenFASTWrappers_Options.windINPfilename = "$(windINPfilename[1:end-4]).bts"
            else
                if !isfile(windINPfilename) || regenWindFiles
                    generateUniformwind(DLCParams[iDLC],windINPfilename)
                end
                modelopt.OWENSOpenFASTWrappers_Options.windINPfilename = windINPfilename
                modelopt.OWENSOpenFASTWrappers_Options.WindType = 2
            end

            modelopt.OWENSAero_Options.ifw = true
            modelopt.OWENS_Options.controlStrategy = DLCParams[iDLC].controlStrategy
            modelopt.DLC_Options.DLCParams = DLCParams[iDLC]
            # run owens simulation
            runScript(modelopt,designparams,path)
        end
    end
end

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
    Vinf_range_used
    analysis_type # "U", "F", "UF"
    controlStrategy # "constRPM", function handle
    RandSeed1 # Turbulent Random Seed Number
    NumGrid_Z # Vertical grid-point matrix dimension
    NumGrid_Y # Horizontal grid-point matrix dimension
    TimeStep # Time step [s]
    HubHt # Hub height [m] (should be > 0.5*GridHeight)
    AnalysisTime # Length of analysis time series [s] (program will add time if necessary)
    GridHeight # Grid height [m]
    GridWidth # Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))
    VFlowAng # Vertical mean flow (uptilt) angle [degrees]
    HFlowAng # Horizontal mean flow (skew) angle [degrees]
    TurbModel # Turbulence model (see Table 4 for valid codes)
    IECstandard # Number of the IEC standard (61400-x, x=1,2,3) with optional 61400-1 ed. number
    IECturbc # IEC turbulence characteristic ("A", "B", "C" or TI in %) or KHTEST
    IEC_WindType # IEC turbulence type ("NTM", "xETM", "xEWM1", or "xEWM50" for x=class 1, 2, or 3)
    RefHt # Height of the reference wind speed [m]
    URef # Mean wind speed at the reference height [m/s] [must be 1-hr mean for API model]
    time
    windvel
    winddir
    windvertvel
    horizshear
    pwrLawVertShear
    LinVertShear
    gustvel
    UpflowAngle
end


function getDLCparams(DLC_case, modelopt, designparams, Vinf_range, Vdesign, Vref, WindChar, WindClass, IEC_std;
    RandSeed1 = 40071,
    grid_oversize=1.2,
    simtime_turbsim=nothing,
    delta_t_turbsim=nothing,
    NumGrid_Z=nothing,
    NumGrid_Y=nothing)
    
    hub_height = designparams[:assembly][:hub_height]
    blade_x = designparams[:components][:blade][:outer_shape_bem][:reference_axis][:x][:values] #Used
    blade_y = designparams[:components][:blade][:outer_shape_bem][:reference_axis][:y][:values] #Used
    blade_z = designparams[:components][:blade][:outer_shape_bem][:reference_axis][:z][:values] #Used
    Blade_Height = maximum(blade_z) #TODO: resolve DLC dependence
    Blade_Radius = maximum(sqrt.(blade_x.^2 .+ blade_y.^2))
    
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

    time = LinRange(0,10,10)
    windvel = nothing # gets supersceded   
    winddir = nothing  
    windvertvel = nothing  
    horizshear = nothing  
    pwrLawVertShear = nothing  
    LinVertShear = nothing  
    gustvel = nothing  
    UpflowAngle = nothing  

    if contains(IEC_std,"1-")
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
            Vinf_range_used = [Vdesign-2.0,Vdesign+2.0]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = [0,0,45,90,45,0,0.0] 
            windvertvel = zeros(length(time))   
            horizshear = zeros(length(time))  
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = zeros(length(time))  
            gustvel = [0,0,7.0,15,7.5,0.0,0.0]  
            UpflowAngle = zeros(length(time))    

        elseif DLC_case == "1_5"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWS\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = zeros(length(time))  
            windvertvel = zeros(length(time))   
            horizshear = [0,0,5.0,0,0,0,0]#ones(length(time)).*10.0   
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = [0,0,0,0,5.0,0,0]#ones(length(time)).*10.0 
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))  
            

        elseif DLC_case == "2_1" || DLC_case == "2_2" || DLC_case == "2_4"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC_case == "2_3"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
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
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
        elseif DLC_case == "3_1"
            ControlStrategy = "startup"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0,30,10)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = ones(length(time)).*0.2  
            LinVertShear = zeros(length(time))   
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))     
            
        elseif DLC_case == "3_2"
            ControlStrategy = "startup"
            Vinf_range_used = [Vinf_range[1];collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
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
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
        elseif DLC_case == "3_3"
            ControlStrategy = "startup"
            Vinf_range_used = [Vinf_range[1];collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = [0,0,45,90,45,0,0.0] 
            windvertvel = zeros(length(time))   
            horizshear = zeros(length(time))  
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = zeros(length(time))  
            gustvel = [0,0,7.0,15,7.5,0.0,0.0]  
            UpflowAngle = zeros(length(time)) 
            
        elseif DLC_case == "4_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0,30,10)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = ones(length(time)).*0.2  
            LinVertShear = zeros(length(time))   
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))     
            
        elseif DLC_case == "4_2"
            ControlStrategy = "shutdown"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
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
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
            
        elseif DLC_case == "5_1"
            ControlStrategy = "emergencyshutdown"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
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

            time = LinRange(0,30,10)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = ones(length(time)).*0.2  
            LinVertShear = zeros(length(time))   
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))  
            
        else
            error("IEC61400_1 DLC_cases such as 1_1, 1_2 defined, you requested $DLC_case")
        end

    elseif contains(IEC_std,"2")
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
            error("IEC61400_2 DLC_cases [1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,3.1,3.2,4.1,5.1,5.2,6.1] defined, you requested $DLC")
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

function generateUniformwind(DLCParams,windINPfilename)

    time = DLCParams.time
    windvel = ones(length(DLCParams.time)) .* DLCParams.URef
    winddir = DLCParams.winddir
    windvertvel = DLCParams.windvertvel
    horizShear = DLCParams.horizshear
    pwrLawVertShear = DLCParams.pwrLawVertShear
    LinVertShear = DLCParams.LinVertShear
    gustvel = DLCParams.gustvel
    UpflowAngle = DLCParams.UpflowAngle

    lines = ["! OpenFAST Deterministic Wind File",
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
    "! (sec)    (m/s)   (Deg)   (m/s)                                            (m/s)   (deg)"]
    # "0.000000   10   0.000000   0          0.000000   0.300000   0.000000   0.000000      8",
    # "10.000000   12   0.000000   0          0.000000   0.300000   0.000000   0.000000      8",

    for itime = 1:length(time)
        lines = [lines; "$(time[itime]) $(windvel[itime]) $(winddir[itime]) $(windvertvel[itime]) $(horizShear[itime]) $(pwrLawVertShear[itime]) $(LinVertShear[itime]) $(gustvel[itime]) $(UpflowAngle[itime])"]
    end

    
    # Write the new file
    open(windINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end
end

function generateTurbsimBTS(DLCParams,windINPfilename,pathtoturbsim;templatefile="$module_path/template_files/templateTurbSim.inp") 

    lines = readlines(templatefile)

    for fieldname in fieldnames(typeof(DLCParams))
        turbsimKeyName = String(fieldname)
        myvalue = getfield(DLCParams,fieldname)
        if turbsimKeyName != "Vinf_range_used" || turbsimKeyName != "analysis_type" || turbsimKeyName != "ControlStrategy"# || other strings
            for (iline,line) in enumerate(lines)
                if contains(line," - ") #TODO: this assumes that the keys aren't in the comments
                    linenocomments,comments = split(line," - ")
                    if contains(linenocomments,turbsimKeyName)
                        value,descriptor = split(linenocomments)
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
function getGustVel(time,nominalVinf,R,G_amp,gustT,gustDelayT)
    ele_x = 0.0 #TODO: I don't think inflowwind takes in account the 3D nature of a vawt

    gustT = gustT * nominalVinf / R
    tr = time .- ele_x .- gustDelayT / R
    if (tr >= 0) && (tr<=gustT)
        IECGustFactor = 1.0 - 0.37 * G_amp/nominalVinf * sin(3*pi*tr/gustT)  * (1.0 - cos(2*pi*tr/gustT))
        return nominalVinf*IECGustFactor
    else
        return nominalVinf
    end

end

function simpleGustVel(time, time_delay, G_amp,gustT)
    timeused = time - time_delay
    if (timeused >= 0) && (timeused<=gustT)
        gustV = -0.37 * G_amp * sin(3*pi*timeused/gustT)  * (1.0 - cos(2*pi*timeused/gustT))
    else
        gustV = 0.0
    end
    return gustV
end

# # Test
# Inp = OWENS.MasterInput(;numTS=3,ifw_libfile="$localpath/../../openfastandy/build/modules/inflowwind/libifw_c_binding")
# OWENS.runDLC(["1_1","1_2"],Inp,localpath;Vinf_range=LinRange(5,20,2),regenWindFiles=true,pathtoturbsim="$localpath/../../openfastandy/build/modules/turbsim/turbsim")
