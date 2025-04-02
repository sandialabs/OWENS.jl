
# Load in all of the modules that we'll be using, and prettify the plotting options
import OWENS
import OWENSAero
import DelimitedFiles
using Statistics:mean
using Test
import FLOWMath
import ForwardDiff
import FiniteDiff
import SNOW

## import PyPlot
## PyPlot.pygui(true)
## PyPlot.rc("figure", figsize=(4.5, 3))
## PyPlot.rc("font", size=10.0)
## PyPlot.rc("lines", linewidth=1.5)
## PyPlot.rc("lines", markersize=3.0)
## PyPlot.rc("legend", frameon=false)
## PyPlot.rc("axes.spines", right=false, top=false)
## PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
## PyPlot.rc("figure",max_open_warning=500)
## plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

path = splitdir(@__FILE__)[1]
# path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]

nothing

function objcon!(constraints,Vars;returnSNOW=true)
    turbine_height = Vars[1]
    turbine_radius = Vars[2]
    chord_scale = [Vars[3],Vars[3]]
    stack_layers_scale = [Vars[4],Vars[4]]
    RPM = Vars[5]
    # Here we set the parametric inputs to the preprocessing and run methods.  Defaults are used otherwise
    # Note that for the OpenFAST coupling, such as adi_lib, setting these to nothing will force the program 
    # to use the precompiled binaries that are part of the OWENS installation.
    # Again, please refer the to API reference for more detailed options of inputs, but of note is that we are 
    # running with AeroModel = "AD" which will automatically generate all of the AeroDyn files necessary to run the
    # structural definition generated here.  Also note that this scripting method relys on the NuMAD formatted CSV files
    # for the composite layup that gets run through PreComp.jl to calculate the sectional properties. Finally, this
    # has been shortened to enable automated deployment, update as desired.

    turbineType = "H-VAWT" # turbine type, for the automatic meshing
    Vinf = 1.2 # inflow velocity
    TSRrange = [3.0]#LinRange(1.0,5.0,2) range of tip speed ratios
    Nslices = 20 # vertical discretizations if DMS or AC aero model
    ntheta = 30 # azimuthal discretizations if DMS or AC aero model
    structuralModel = "TNB"
    ntelem = 100 # tower elements
    nbelem = 30 # blade elements
    nselem = 10 # strut elements
    ifw = false # use inflow wind, if DMS or AC aero model
    numTS = 20#321 # number of simulation time steps
    delta_t = 0.01 # simulation time step spacing
    adi_lib = nothing#"$path/../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" 
    adi_rootname = "$path/RM2" # path and name that all the aerodyn files are saved with
    VTKsaveName = "$path/vtk/RM2_medium" # path and name that all OWENS VTK files are saved with
    tsave_idx = 1:1:numTS-1 #you don't have to save every timestep in VTK 
    ifw_libfile = nothing#"$path/../../../openfast/build/modules/inflowwind/libifw_c_binding"
    fluid_density = 1000.0
    fluid_dyn_viscosity = 1.792E-3
    AddedMass_Coeff_Ca = 1.0 #For structural side added mass 
    Aero_Buoyancy_Active = true # For buoyancy forcing, handled by the OWENSAero module 
    AeroModel = "DMS"
    verbosity = 1 # verbosity level where higher is more
    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    eta = 0.5 # blade-strut mount point ratio, fraction from leading edge
    number_of_blades = 3 # number of blades
    towerHeight = 0.2165 # height of the tower past the blades
    Blade_Height = turbine_height#0.807 # height of the blades
    Blade_Radius = turbine_radius#0.5375 # radius of the turbine

    NuMad_geom_xlscsv_file_twr = "$path/data_RM2/TowerGeom.csv" # 
    NuMad_mat_xlscsv_file_twr = "$path/data_RM2/TowerMaterials.csv"
    NuMad_geom_xlscsv_file_bld = "$path/data_RM2/GeomBlades$AeroModel.csv"
    NuMad_mat_xlscsv_file_bld = "$path/data_RM2/materials_NuMAD.csv"
    NuMad_geom_xlscsv_file_strut = "$path/data_RM2/GeomStruts$AeroModel.csv"
    NuMad_mat_xlscsv_file_strut = "$path/data_RM2/materials_NuMAD.csv"

    nothing

    # For this simulation, using aerodyn ("AD"), we will use a turbulent inflow file, indicated as WindType 3
    # We can let the built in OWENS library for turbsim generate the file as so. Modify the inp file to your liking and 
    # comment out the run command to run your own with more time
    WindType = 3
    windINPfilename = "$path/data_RM2/3mx3m1pt2msNTM.bts"
    # run(`$(OWENS.OWENSOpenFASTWrappers.turbsim()) $(windINPfilename[1:end-3])inp`)

    nothing
    # Here we would set up to run multiple different inflow conditions, and we can
    # allow the preprocessing to re-run each time since it only takes a few seconds.
    # The for loop is commented out to enable this example to be parsed by the literate 
    # package that auto generates the docs.

    CP = zeros(length(TSRrange))
    ## for (iTSR,TSR) in enumerate(collect(TSRrange))
    TSR = TSRrange[1]
    iTSR = 1
    
    # omega = Vinf/Blade_Radius*TSR  
    # RPM = omega * 60 / (2*pi)

    # println(RPM)

    nothing

    # Here we run the OWENS preprocessing, which creates the beam mesh and joints, calculates the sectional properties
    # maps the sectional properties to the elements.  It also sets up the aero module being used.

    shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
    helix_angle = 0.0#-pi/4
    shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius
    shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius

    mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
    mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
        rho=fluid_density,
        mu=fluid_dyn_viscosity,
        Nslices,
        ntheta,
        RPM,
        Vinf,
        eta,
        B=number_of_blades,
        H = Blade_Height, 
        R = Blade_Radius, 
        shapeZ, 
        shapeX,
        shapeY, 
        ifw,
        WindType,
        delta_t,
        numTS,
        adi_lib,
        adi_rootname,
        windINPfilename,
        ifw_libfile,
        stack_layers_scale,
        chord_scale,
        NuMad_geom_xlscsv_file_twr,
        NuMad_mat_xlscsv_file_twr,
        NuMad_geom_xlscsv_file_bld,
        NuMad_mat_xlscsv_file_bld,
        NuMad_geom_xlscsv_file_strut,
        NuMad_mat_xlscsv_file_strut,
        Htwr_base=towerHeight,
        Htwr_blds = Blade_Height+towerHeight,
        ntelem, 
        nbelem, 
        nselem,
        joint_type = 0,
        strut_twr_mountpoint = [0.5],
        strut_bld_mountpoint = [0.5],
        AeroModel, #AD, DMS, AC
        DynamicStallModel="BV",
        Aero_AddedMass_Active = false,
        Aero_RotAccel_Active = false,
        AddedMass_Coeff_Ca,
        Aero_Buoyancy_Active,
        verbosity,
        meshtype = turbineType)

    nothing

    # PyPlot.figure()
    # for idot = 1:length(sectionPropsArray[170].xaf)
    #     PyPlot.scatter(sectionPropsArray[170].xaf[idot],sectionPropsArray[170].yaf[idot])
    #     sleep(0.001)
    # end

    ## This plots the mesh and node numbering of the resulting mesh and overlays the joint connections

    ## PyPlot.figure()
    ## for icon = 1:length(mymesh.conn[:,1])
    ##     idx1 = mymesh.conn[icon,1]
    ##     idx2 = mymesh.conn[icon,2]
    ##     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    ##     PyPlot.plot3D([1,1],[1,1],[1,1],"k.-")
    ##     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    ## end
    ## for ijoint = 1:length(myjoint[:,1])
    ##     idx2 = Int(myjoint[ijoint,2])
    ##     idx1 = Int(myjoint[ijoint,3])
    ##     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    ##     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",color="r",ha="center",va="center")
    ##     PyPlot.text3D(mymesh.x[idx2].+rand()/30,mymesh.y[idx2].+rand()/30,mymesh.z[idx2].+rand()/30,"$idx2",color="r",ha="center",va="center")
    ## end
    ## PyPlot.xlabel("x")
    ## PyPlot.ylabel("y")
    ## PyPlot.zlabel("z")
    ## PyPlot.axis("equal")

    nothing

    # Here we apply the boundary conditions.  For this case, the tower base node which is 
    # 1 is constrained in all 6 degrees of freedom to have a displacement of 0, and the tower top is also constrained.
    # You can change this displacement to allow for things like pretension, and you can apply boundary conditions to any node.

    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0
    AD15bldElIdxRng[1,2]-1 1 0
    AD15bldElIdxRng[1,2]-1 2 0
    AD15bldElIdxRng[1,2]-1 3 0
    AD15bldElIdxRng[1,2]-1 4 0
    AD15bldElIdxRng[1,2]-1 5 0]

    nothing

    # There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options: ? OWENS.Inputs()

    tocp = [0.0;10.0; 1e6]
    Omegaocp = [RPM; RPM; RPM]./60 #control inputs
    
    tocp_Vinf = [0.0;10.0; 1e6]
    Vinfocp = [Vinf;Vinf;Vinf]

    inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
    tocp,
    dataOutputFilename = "./InitialDataOutputs_scripting.out",
    Omegaocp,
    tocp_Vinf,
    Vinfocp,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = 2,
    MAXITER = 20)

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options

    FEAinputs = OWENS.FEAModel(;analysisType = structuralModel,
    dataOutputFilename = "./InitialDataOutputs_scripting.out",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = true,
    gravityOn = [0,0,9.81], #positive since the turbine is upside down
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    AddedMass_Coeff_Ca,
    iterationType = "DI")

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time. Note that for the sake of quickly deployable docs, a shortened simulation is conducted and the CP will not be correct

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
    topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=FEAinputs,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,turbsimfile = windINPfilename)

    area = Blade_Height*2*Blade_Radius
    full_rev_N_timesteps = round(Int,RPM/60/delta_t)
    if full_rev_N_timesteps>numTS
        idx_start = 1
    else
        idx_start = numTS-full_rev_N_timesteps
    end
    CP[iTSR] = mean(FReactionHist[idx_start:end,6].*OmegaHist[idx_start:end]*2*pi)/(0.5*fluid_density*mean(Vinfocp)^3*area)
    TSR = mean(OmegaHist*2*pi*Blade_Radius/mean(Vinfocp))
    ReD = fluid_density*mean(Vinfocp)*Blade_Radius*2/fluid_dyn_viscosity

    nothing

    # We can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function
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
    mymesh,myel,myort,number_of_blades,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor 
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,delta_t,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing)

    # OWENS.outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)
        
    # OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
    #     epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    #     FReactionHist,topFexternal_hist;tsave_idx)

    # Formulate our objective function to be a pseudo LCOE.  You should consider scaling your outputs to achieve well scaled gradients
    power = (mean(-FReactionHist[:,6])*(RPM*2*pi/60))
    pseudoLCOE = massOwens/power

    # Constraints on the minimum allowable safety factor, the minimum allowable lifetime damage, and that the power must be greater than 0.0 so the optimizer doesn't invert the LCOE equation.
    minSF = FLOWMath.ksmin(SF_ult_U)
    # println(FLOWMath.ksmax(topDamage_blade_U))
    maxFatiguePer20yr = FLOWMath.ksmax(topDamage_blade_U/t[end]*60*60*20*365*24,1e10)

    constraints[1] = 1.0 - minSF # 1.0<SF
    constraints[2] = maxFatiguePer20yr - 1.0 # fatigueDamage < 1.0
    constraints[3] = 0.0 - power # 0<power i.e. power must be positive

    println("Vars: $Vars")
    println("Obj: $pseudoLCOE")
    println("Con: $constraints")

    if returnSNOW
        return pseudoLCOE
    else
        return [pseudoLCOE;constraints]
    end
end


# Let's test the function
Vars = zeros(5)
Vars[1] = 0.807#turbine_height
Vars[2] = 0.5375#turbine_radius
Vars[3] = 1.0#chord_scale
Vars[4] = 1.0#stack_layers_scale
Vars[5] = 30.0#RPM

constraints = zeros(3)

# First run the function alone
objective = objcon!(constraints,Vars)
println("Objective: $objective")
println("constraints: $constraints")

# Now create a wrapper that is compatible with the diff tools
objcon0(x) = objcon!(constraints,x;returnSNOW=false)

# now create a wrapper so it is compatible with FiniteDiff
function objcon2(x)
    return Float64.(objcon0(x))
end

# Run finite diff
J2 = FiniteDiff.finite_difference_jacobian(objcon2,Vars)

# Run ForwardDiff
J = ForwardDiff.jacobian(objcon0, Vars)

#Check that they are the same
for i1 = 1:length(J[:,1])
    for i2 = 1:length(J[1,:])
        @test isapprox(J[i1,i2],J2[i1,i2];atol=1e-5)
    end
end

# # Now let's throw a real optimizer at it.

# x0 = Vars # starting point
# lx = Vars./2  # lower bounds on x
# ux = Vars.*2  # upper bounds on x
# ng = 3  # number of constraints
# lg = -Inf*ones(ng)  # lower bounds on g
# ug = zeros(ng)  # upper bounds on g
# options = SNOW.Options(solver=SNOW.IPOPT())  # choosing IPOPT solver

# # Run the optimization
# xopt, fopt, info = SNOW.minimize(objcon!, x0, ng, lx, ux, lg, ug, options)

# println("xstar = ", xopt)
# println("fstar = ", fopt)
# println("info = ", info)
