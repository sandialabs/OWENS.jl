import OWENS
import HDF5
import Test
import YAML
import OrderedCollections
import FLOWMath
import Statistics:mean
using SNOW

# Monday: write lots of descriptions, add more to windio 

runpath = splitdir(@__FILE__)[1]

function messyoptfun!(constraints,Vars)

    OWENS_Options = OWENS.MasterInput("$runpath/modeling_options_OWENS_windioExample.yml")

    WINDIO_filename = "$runpath/WINDIO_example.yaml"
    windio = YAML.load_file(WINDIO_filename; dicttype=OrderedCollections.OrderedDict{Symbol,Any})

    # Design Variables
    quarter_span_radius = Vars[1]#24.0
    max_radius = Vars[2]#42.0
    RPM = Vars[3]#25.0

    OWENS_Options.RPM = RPM

    x_control_pts = [0.0,quarter_span_radius,max_radius,quarter_span_radius,0.0]
    x_control_pts_grid = [0.0,0.25,0.5,0.75,1.0]

    x_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:grid]

    x_values = FLOWMath.akima(x_control_pts_grid,x_control_pts,x_grid)
    windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:values] = x_values 

    # Run
    OWENS.runOWENSWINDIO(windio,OWENS_Options,runpath)

    # Post Processing
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

    # objective
    power = (mean(FReactionHist[:6])*(RPM*2*pi/60))
    pseudoLCOE = massOwens/power

    # Constraints
    minSF = FLOWMath.ksmin(SF_ult_U)
    maxFatiguePerHour = FLOWMath.ksmax(topDamage_blade_U)/t[end]*60*60
    maxFatiguePer20yr = maxFatiguePerHour*20*365*24

    constraints[1] = 1.0 - minSF # 1.0<SF
    constraints[2] = maxFatiguePer20yr - 1.0

    return pseudoLCOE

end

# Test function
Vars = [24.0, #quarter height radius
42.0, #mid span radius
15.0] # RPM
constraints = zeros(2)
objective = messyoptfun!(constraints,Vars)

println("Objective: $objective")
println("constraints: $constraints")

x0 = [24.0; 42.0; 15.0]  # starting point
lx = [10.0; 20.0; 5.0]  # lower bounds on x
ux = [20.0; 60.0; 30.0]  # upper bounds on x
ng = 2  # number of constraints
lg = -Inf*ones(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
options = Options(solver=IPOPT())  # choosing IPOPT solver

xopt, fopt, info = minimize(messyoptfun!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)