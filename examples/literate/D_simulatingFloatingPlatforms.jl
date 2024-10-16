# # [Simulating Floating Platforms](@id simple2)
# 
# In this example, we make the OWENS model more complex by simulation a floating wind
# turbine instead of a land-based one. This involves coupling two meshes together: one
# representing the wind turbine (this will be similar to meshes used in previous tutorials),
# and one representing the floating platform.
#
# The wind turbine mesh (or "topside") will be created similarly to the previous tutorial,
# though we will overwrite some of the inputs so the simulation works better for a floating system.
# This will help with revealing some more internals of OWENS and the potential inputs
# available to users, which is helpful for explaining how the floating platform mesh
# (or "bottom side") is defined and differs from the topside.


import OWENS
import OWENSFEA

path = runpath = "/home/runner/work/OWENS.jl/OWENS.jl/docs/src/literate" #splitdir(@__FILE__)[1]
verbosity = 1

# First, create the inputs for the topside mesh as done previous in tutorials A and B.
Inp = OWENS.MasterInput("$runpath/sampleOWENS.yml")

analysisType = "TNB" # this differs from previous tutorials, as it has been verified to work well with floating capabilities
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
AModel = Inp.AModel
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

# We will continue to use helper functions here to fully define the topside mesh, sectional
# properties, and their connection. However, note that our naming convention will be
# different in order to differentiate this mesh from the
# bottom side mesh. The outputs to the helper functions not being used elsewhere in the
# tutorial are returned as "_" to improve clarity. TODO finish adding _s

topMesh, topEl, topOrt, topJoint, topSectionProps, _, _,
_, _,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho, Nslices, ntheta, RPM, Vinf, eta, B, H, R, shapeZ, shapeX,
    ifw, WindType, delta_t, numTS, adi_lib, adi_rootname, indINPfilename, ifw_libfile,
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
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)

nothing

# Now, unlike the previous examples, we will **not** apply a boundary condition to this mesh.
# The boundary condition previously was essentially fixing the bottom of the turbine tower in place,
# since it is installed in the ground.
# However, for our floating platform, we want the loads from the topside to propagate to the platform,
# so our restoring hydrodynamics and mooring loads will "constaint" the combined meshes as we want.

topBC = []

nothing

# Bottom side mesh definition TODO -- be sure to mention that it's a single-element mesh that represents
# a rigid body.

# We also want to provide some additional general inputs to enable floating simulation.
# `hydroOn` being the principal and obvious one, but also the input files needed by
# HydroDyn and MoorDyn, which are the external libraries that calculate the hydrodynamic
# and mooring loads to send to the bottom side mesh. These input files include:
#  - interpOrder: the degree of interpolation used to predict the future states within
#    HydroDyn and MoorDyn, used in the `HD_UpdateStates` and `MD_UpdateStates` functions
#  - hd_input_file: the base HydroDyn input file with a `.dat` extension
#  - md_input_file: the base MoorDyn input file with a `.dat` extension
#  - ss_input_file: the sea state input file with a `.dat` extension used to define the environmental conditions of the sea.
#    This is used within HydroDyn.
#  - potfilefile: the directory containing the potential flow files.
#    At minimum, the directory must contain the .1, .3, and .hst WAMIT output files.
#    This is also used within HydroDyn.
# See the OpenFAST documentation for HydroDyn (https://openfast.readthedocs.io/en/main/source/user/hydrodyn/input_files.html)
# and MoorDyn (https://moordyn.readthedocs.io/en/latest/inputs.html) for more
# information about how to format these input files. For simplicity here, we will
# use predefined input files for the OC4 semisubmersible platform, which comes with
# OpenFAST and is copied to the `data` folder here.


if AModel=="AD"
    AD15On = true
else
    AD15On = false
end
hydroOn = true
interpOrder = 2
hd_input_file = "data/HydroDyn.dat"
md_input_file = "data/MoorDyn.dat"
ss_input_file = "data/SeaState.dat"
potflowfile = "data/potflowdata"

inputs = OWENS.Inputs(;analysisType = structuralModel,
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2,
hydroOn,
interpOrder,
hd_input_file,
md_input_file,
ss_input_file,
potflowfile)

nothing

# As in previous examples, we need to define the inputs for the finite element models, though now we need
# a definition for both the topside and the bottom side.

# The changes to the topside from example A are the lack of boundary conditions,
# non-default gamma and alpha terms (for better convergence with the platform mesh),
# and gravityOn specified TODO do I need gravityOn or is the default fine?

topFEAModel = OWENS.FEAModel(;analysisType = structuralModel,
outFilename = "none",
joint = topJoint,
platformTurbineConnectionNodeNumber = 1,
nlOn = true,
numNodes = topMesh.numNodes,
RayleighAlpha = 0.05,
RayleighBeta = 0.05,
iterationType = "DI",
gamma = 1.0,
alpha = 0.5,
gravityOn = [0.0, 0.0, 9.80665])

nothing

# The bottom finite element mesh will use the same inputs, but also
# includes a concentrated term at its bottom node to represent the rigid body mass
# of the platform.
# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

nothing

# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

nothing

# This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
# it also calculates analytical simply supported buckling safety factors

##########################################
#### Ultimate Failure #####
##########################################

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
topstrainout_tower_U,topstrainout_tower_LtopDamage_blade_U,
topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
Twr_LE_U_idx=1,Twr_LE_L_idx=1,
AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

nothing