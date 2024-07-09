# WindIO Yaml File Contains, for the turbine:
# Assembly: general size and turbulence conditions
# Components: general blade, tower, strut shape and joints
# Airfoils: geometry, polars across Reynolds numbers
# materials: material properties
# Control: gains
# Environment: properties such as density
# BOS: balance of system variables
# Costs: cost components

# As of the writing of this file, the following components have their own files but are planned to be switched to the WindIO file input
# Airfoils: airfoil polar files (see examples/Optimization/airfoils/ .dat files), airfoil geometry (see examples/Optimization/airfoils/ .txt files)
# Components->internal structure: Numad Files (see the examples/Optimization/data files and page 85 of https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf).  
# Materials: Numad Files in the examples/Optimization/data folder, see also page 85 of https://energy.sandia.gov/wp-content/gallery/uploads/NuMAD_UserGuide_SAND2012-7028.pdf.  
# The specifc items being used from the WindIO yaml format currently (which is subject to rapid development) can be seen in OWENS.jl/src/windio.jl

# Then the OWENS modeling options currently control 
# Controls: currently constant RPM (in work is to expose the scripting options of the discon and custom control options at the super simple file level)
# analysisType = "M" for modal, "S" for steady, "TNB" for time newmark time stepping (must be paired with TNB or GX structural models), "ROM" time stepping (must be paired with ROM structural model), 
# turbineType: "Darrieus" (where the blade ends connect to the tower), "H-VAWT" (where the blade ends do not connect to the tower), "ARCUS" where there are only blades and a cable for each blade from the top to the stubby tower base.
# eta: blade mount point ratio, or the chordwise normalized position of the rotation axis on the airfoil, extremely important as it changes the effective pitch offset for a VAWT
# towerHeight: Base tower offset to the blades in the mesh
# Vinf: Windspeed if ifw==false, or if WindType=1
# RPM: Rotational rate if constant RPM control, is also the initialized RPM
# Nslices: int, Number of vertical slices in the strip theory to get the simplified aero model to be 3D
# ntheta: int, Number of azimuthal discretizations of the simplified aero model, if active, must by wholly divisible by NBlades
# structuralModel: TNB time newmark beta time stepping linear or nonlinear, ROM reduced order modal linear or nonlinear, GX geometrically exact beam theory 
# ntelem: int, number of tower elements, in the automatically generated mesh
# nbelem: int, number of blade elements, in the automatically generated mesh
# ncelem: int, number of cable elements, if active, in the automatically generated mesh
# nselem: int, number of strut elements, in the automatically generated mesh
# ifw: true or false, tells the AC or DMS models to use constant inflow or to read the specified windINPfilename
# WindType: Inflow wind wind type (1 2 3), this is automatically handled by the DLC processor when it is used
# AModel: Aero model, "AC" actuator cylinder (slower but slightly more accurate with array capability at the planar aero only level) "DMS" (faster and nearly as accurate) or "AD" (aerodyn interface, specifically OLAF, which is a relatively slow but higher fidelity vortex line wake model.  This includes strut drag and tip losses.)
# windINPfilename: .bts or .inp files, see inflow wind documentation.  If running with the DLC script, it will automatically generate these and run with them, and not generate them if they already exist unless told to do so.
# ifw_libfile: path to the inflow wind interface library, "nothing" if you are using the one that is automatically build (on mac/linux) with OWENSOpenFASTWrappers.jl
# numTS: total number of time steps
# delta_t: time step, in seconds
# NuMad_geom_xlscsv_file_twr: numad file inputs for twr geom/composite inputs, this is also where the airfoil name is specified
# NuMad_mat_xlscsv_file_twr: numad file inputs for twr materials
# NuMad_geom_xlscsv_file_bld: numad file inputs for bld geom/composite inputs, this is also where the airfoil name is specified
# NuMad_mat_xlscsv_file_bld: numad file inputs for bld materials
# NuMad_geom_xlscsv_file_strut: numad file inputs for strut geom/composite inputs, this is also where the airfoil name is specified
# NuMad_mat_xlscsv_file_strut: numad file inputs for strut materials
# adi_lib: path to the aerodyn interface library, "nothing" if you are using the one that is automatically build (on mac/linux) with OWENSOpenFASTWrappers.jl


import OWENS # Main module for OWENS
import HDF5 # used for loading and saving binary HDF5 files that OWENS saves
import YAML # used to load YAML files and create dictionaries of the inputs
import OrderedCollections # used by the YAML reader to create a dictionary that is in the same order as the file
import FLOWMath # used to create optimization friendly splines and max calculations
import Statistics:mean # used to calculate the average of an array
using SNOW # the sparse nonlinear optimization wrapper that provides a friendly interface for these types of large evaluation problems

runpath = splitdir(@__FILE__)[1] # This path variable ensures that you can invoke this script from anywhere and the program will know where this file is and the subsequent dependent files like airfoil polars, etc.

# This is an example objective and constraint function, which passes by reference in the constraint array that gets mutated by the function, while the variables are the design inputs from the optimizer.
# The return is simply the objective 
function messyoptfun!(constraints,Vars)

    OWENS_Options = OWENS.MasterInput("$runpath/modeling_options_OWENS_windioExample.yml")
    # These options are the OWENS specific runtime options, as opposed to the design options coming from the WindIO

    WINDIO_filename = "$runpath/WINDIO_example.yaml"
    windio = YAML.load_file(WINDIO_filename; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
    # The windio filename can be passed into the runtime function below, or you can read the yaml and mutate it directly here to bypass filereading.
    # It is planned to have all aspects of the WINDIO design details incorporated as the main input method to OWENS.  
    # However, due to the complexity of the transition, the airfoil polars, geometry, and composite layup information and material properties are still in work


    # Design Variables, here are three examples of possible inputs for a constant RPM simulation.  
    quarter_span_radius = Vars[1]#24.0
    max_radius = Vars[2]#42.0
    RPM = Vars[3]#25.0

    # Additional design variables recommended for initial point design
    # chord control points
    # composite thickness control points
    # twist offset control points

    # Additional design varibles recommended for other critical DLCs (parked and emergency shutdown)

    # Since we are using constant RPM and there isn't a constant RPM option in the WindIO control block, it is within the OWENS Options
    OWENS_Options.RPM = RPM

    # Now, rather than waste dozens of design variables on the overall shape of the blade, take the two control points and spline them to a full definition
    # Note that splining at this level gives you more control, but the same thing happens internally, so you could specify the sparse control points and it will do this for you.
    x_control_pts = [0.0,quarter_span_radius,max_radius,quarter_span_radius,0.0]
    x_control_pts_grid = [0.0,0.25,0.5,0.75,1.0]

    x_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:grid]
    x_values = FLOWMath.akima(x_control_pts_grid,x_control_pts,x_grid)

    windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:values] = x_values 

    # Run with the WindIO standard, which as mentioned before, is in progress to include all of the design components.
    OWENS.runOWENSWINDIO(windio,OWENS_Options,runpath)

    # Read in all of the current OWENS outputs from the run. This is a portion of what is calculated internally.
    file = "$runpath/InitialDataOutputs.h5"
    t = HDF5.h5read(file,"t") # Array of time that the simulation was run at
    aziHist = HDF5.h5read(file,"aziHist") # Array of azimuthal positions for each time
    OmegaHist = HDF5.h5read(file,"OmegaHist") # Array of rotational rates at each time
    OmegaDotHist = HDF5.h5read(file,"OmegaDotHist") # Array of rotational acceleration at each time 
    gbHist = HDF5.h5read(file,"gbHist") # if active, gearbox model azimuthal position
    gbDotHist = HDF5.h5read(file,"gbDotHist") # if active, gearbox model rotational speed
    gbDotDotHist = HDF5.h5read(file,"gbDotDotHist") # if active, gearbox model rotational acceleration
    FReactionHist = HDF5.h5read(file,"FReactionHist") # array sized (N_timesteps,NDOF) of reaction forces where each DOF is aligned with the mesh nodes, so first node is values 1:6, second is 7:12, hub FOR x,y,z, mx,my,mz forces and moments
    FTwrBsHist = HDF5.h5read(file,"FTwrBsHist") # if running floating, gives the tower base forces and moments
    genTorque = HDF5.h5read(file,"genTorque") # if generator model active, history of generator torque
    genPower = HDF5.h5read(file,"genPower") #  if generator model active, history of generator power
    torqueDriveShaft = HDF5.h5read(file,"torqueDriveShaft") # if drivetrain model active, driveshaft torque
    uHist = HDF5.h5read(file,"uHist") # displacement history sized (N_timesteps, NDOF), maps to the mesh nodes similarly to FReactionHist
    uHist_prp = HDF5.h5read(file,"uHist_prp") # if active, platform displacement history
    epsilon_x_hist = HDF5.h5read(file,"epsilon_x_hist") # beam strain along the span for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    epsilon_y_hist = HDF5.h5read(file,"epsilon_y_hist") # beam strain in the chorwise direction for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    epsilon_z_hist = HDF5.h5read(file,"epsilon_z_hist") # beam strain in the chord thickness direction for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    kappa_x_hist = HDF5.h5read(file,"kappa_x_hist") # beam curvature about spanwise direction (twist) for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    kappa_y_hist = HDF5.h5read(file,"kappa_y_hist") # beam curvature about the chordwise direction (flap) for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    kappa_z_hist = HDF5.h5read(file,"kappa_z_hist") # beam curvature about the chord thickness direction (lag) for each element in the mesh, sized (NQuadPoints,N_Elements,N_timesteps)
    massOwens = HDF5.h5read(file,"massOwens") # integral mass of the composite structures
    stress_U = HDF5.h5read(file,"stress_U") # composite stress of the nominal low pressure side of the blade in the upper laminate, sized (N_timesteps,N_composite_input_stations,N_chordwise_input_stations, (principal,transverse, shear))
    SF_ult_U = HDF5.h5read(file,"SF_ult_U") # ultimate safety factor for the specified failure model, size(N_timesteps,N_composite_input_stations,N_chordwise_input_stations) 
    SF_buck_U = HDF5.h5read(file,"SF_buck_U") # simple buckling safety factor, size(N_timesteps,N_composite_input_stations,N_chordwise_input_stations) 
    stress_L = HDF5.h5read(file,"stress_L") # for blade lower surface, same as upper
    SF_ult_L = HDF5.h5read(file,"SF_ult_L") # for blade lower surface, same as upper
    SF_buck_L = HDF5.h5read(file,"SF_buck_L") # for blade lower surface, same as upper
    stress_TU = HDF5.h5read(file,"stress_TU") # for tower upper surface, same format as blade
    SF_ult_TU = HDF5.h5read(file,"SF_ult_TU") # for tower upper surface, same format as blade
    SF_buck_TU = HDF5.h5read(file,"SF_buck_TU") # for tower upper surface, same format as blade
    stress_TL = HDF5.h5read(file,"stress_TL") # for tower lower surface, same format as blade
    SF_ult_TL = HDF5.h5read(file,"SF_ult_TL") # for tower lower surface, same format as blade
    SF_buck_TL = HDF5.h5read(file,"SF_buck_TL") # for tower lower surface, same format as blade
    topstrainout_blade_U = HDF5.h5read(file,"topstrainout_blade_U") # strain in the composite frame of reference, for upper side of blade, sized (N_timesteps,N_composite_input_stations,N_chordwise_input_stations, (principal,transverse, shear, beam strain x,beam strain y,beam strain z,beam curvature x,beam curvature y,beam curvature z))
    topstrainout_blade_L = HDF5.h5read(file,"topstrainout_blade_L") # same as above, but for blade lower 
    topstrainout_tower_U = HDF5.h5read(file,"topstrainout_tower_U") # same as above, but for lower upper surface
    topstrainout_tower_L = HDF5.h5read(file,"topstrainout_tower_L") # same as above, but for lower upper surface
    topDamage_blade_U = HDF5.h5read(file,"topDamage_blade_U") # for the blade upper miner's damage accumulated for the simulated time, in the worst case ply, size (N_composite_input_stations,N_chordwise_input_stations)
    topDamage_blade_L = HDF5.h5read(file,"topDamage_blade_L") # for the blade lower, like above 
    topDamage_tower_U = HDF5.h5read(file,"topDamage_tower_U") # for the tower upper, like above 
    topDamage_tower_L = HDF5.h5read(file,"topDamage_tower_L") # for the tower lower, like above

    # Formulate our objective function to be a pseudo LCOE.  You should consider scaling your outputs to achieve well scaled gradients
    power = (mean(-FReactionHist[:,6])*(RPM*2*pi/60))
    pseudoLCOE = massOwens/power

    # Constraints on the minimum allowable safety factor, the minimum allowable lifetime damage, and that the power must be greater than 0.0 so the optimizer doesn't invert the LCOE equation.
    minSF = FLOWMath.ksmin(SF_ult_U)
    println(FLOWMath.ksmax(topDamage_blade_U))
    maxFatiguePer20yr = FLOWMath.ksmax(topDamage_blade_U/t[end]*60*60*20*365*24,1e10)

    constraints[1] = 1.0 - minSF # 1.0<SF
    constraints[2] = maxFatiguePer20yr - 1.0 # fatigueDamage < 1.0
    constraints[3] = 0.0 - power # 0<power i.e. power must be positive

    println("Vars: $Vars")
    println("Obj: $pseudoLCOE")
    println("Con: $constraints")

    return pseudoLCOE

end

# Let's test the function
Vars = [24.0, #quarter height radius
42.0, #mid span radius
15.0] # RPM
constraints = zeros(3)
objective = messyoptfun!(constraints,Vars)

println("Objective: $objective")
println("constraints: $constraints")

# # Now let's throw a real optimizer at it.  !!!NOTE!!!: This example may need more design freedom and variables to produce a feasible design.

# x0 = [24.0; 42.0; 15.0]  # starting point
# lx = [10.0; 20.0; 5.0]  # lower bounds on x
# ux = [20.0; 60.0; 30.0]  # upper bounds on x
# ng = 3  # number of constraints
# lg = -Inf*ones(ng)  # lower bounds on g
# ug = zeros(ng)  # upper bounds on g
# options = Options(solver=IPOPT())  # choosing IPOPT solver

# # Run the optimization
# xopt, fopt, info = minimize(messyoptfun!, x0, ng, lx, ux, lg, ug, options)

# println("xstar = ", xopt)
# println("fstar = ", fopt)
# println("info = ", info)

# Formatting for 