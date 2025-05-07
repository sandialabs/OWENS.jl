# # [Simulating Floating Platforms](@id simple4)
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
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`D_simulatingFloatingPlatforms.ipynb`](@__NBVIEWER_ROOT_URL__/examples/D_simulatingFloatingPlatforms.ipynb).
#-

## import PyPlot

import OWENS
import OWENSFEA
import OWENSAero

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate"
##runpath = path = splitdir(@__FILE__)[1] # use to run locally

nothing

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
AeroModel = "DMS"#Inp.AeroModel
windINPfilename = "$(path)$(Inp.windINPfilename)"
ifw_libfile = Inp.ifw_libfile
if ifw_libfile == "nothing"
    ifw_libfile = nothing
end
Blade_Height = Inp.Blade_Height
Blade_Radius = Inp.Blade_Radius
numTS = 5 # shortened since the floating simulation can take awhile
delta_t = Inp.delta_t
NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld[1:end-4])DMS.csv"
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

# We will continue to use helper functions here to fully define the topside mesh, sectional
# properties, and their connection. However, note that our naming convention will be
# different in order to differentiate this mesh from the
# bottom side mesh. The outputs to the helper functions not being used elsewhere in the
# tutorial are returned as "_" to improve clarity.
#
# Note that, unlike the previous examples, we will **not** apply a boundary condition to this mesh.
# The boundary condition previously was essentially fixing the bottom of the turbine tower in place,
# since it is installed in the ground.
# However, for our floating platform, we want the loads from the topside to propagate to the platform,
# so our restoring hydrodynamics and mooring loads will "constain" the combined meshes as we want.

topMesh, topEl, topOrt, topJoint, topSectionProps,
_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
aeroForces, deformAero, _, _,topSystem,topAssembly, sections, _, _ = OWENS.setupOWENS(OWENSAero,path;
    rho, Nslices, ntheta, RPM, Vinf, eta, B, H, R, shapeZ, shapeX,
    ifw, WindType, delta_t, numTS, adi_lib, adi_rootname, windINPfilename, ifw_libfile,
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

# Now we need to define the properties of the bottom side mesh, which we will do manually here
# instead of reading in values from a YAML file.
# This mesh is very simple: it consists of a single element with a nodal connection on either end.
# OWENS assumes the platform is a rigid body, so the element has very high stiffness to effectively function as a point mass.

nothing

# Since the mesh is so simple, most of these inputs are trivial and affect little.
# However, a few items are of note:
#   - The z positions in `bottomMesh` start at zero for the first node, and the second node is at the
#     height of the tower base above the platform reference point (aka the platform location at the still water
#     level above the platform centroid).
#     In the case of the OC4 platform, this is 10 meters. This is also reflected in the Length argument in bottomOrt
#   - The `bottomSectionProps` structure has very large values for EIyy, EIzz, GJ, and EA, which are all the
#     relevant stiffness properties, in order to effectively model the platform as a rigid body.
#   - The `zcm` parameter in `bottomSectionProps` is set to the vertical distance of the platform with respect to the platform reference point.
#     For the OC4 platform, this is 8.6588 meters below sea level.

bottomMesh = OWENSFEA.Mesh(
    [1, 2], #nodeNum
    1, #numElx
    2, #numNodes
    [0.0, 0.0], #x
    [0.0, 0.0], #y
    [0.0, 10.0], #z
    [1], #elNum
    [1 2], #conn
    [0], #type
    [0, 1], #meshSeg
    zeros(1,1), #structuralSpanLocNorm
    zeros(Int,1,1), #structuralNodeNumbers
    zeros(Int,1,1) #structuralElNumbers
    )
numBottomNodes = length(bottomMesh.z)
bottomOrt = OWENSFEA.Ort(
    [180.0], #Psi_d
    [-90.0],  #Theta_d
    [90.0], #Twist_d
    [10.0], #Length
    ones(1,1), #elNum
    zeros(3,1)#Offset
    )
n_ptfm_elem = 1
bottomSectionProps = Array{OWENSFEA.SectionPropsArray,1}(undef, n_ptfm_elem)
for ii = 1:n_ptfm_elem
    bottomSectionProps[ii] = OWENSFEA.SectionPropsArray(
        [0.0, 0.0], #ac
        [0.0, 0.0], #twist_d
        [0.0, 0.0], #rhoA
        [1e18, 1e18], #EIyy
        [1e18, 1e18], #EIzz
        [1e18, 1e18], #GJ
        [1e18, 1e18], #EA
        [0.0, 0.0], #rhoIyy
        [0.0, 0.0], #rhoIzz
        [0.0, 0.0], #rhoJ
        [-8.6588, 0.0], #zcm
        [0.0, 0.0], #ycm,
        [0.0, 0.0], #a
        [0.0, 0.0], #EIyz
        [0.0, 0.0], #alpha1
        [0.0, 0.0], #alpha2
        [0.0, 0.0], #alpha3
        [0.0, 0.0], #alpha4
        [0.0, 0.0], #alpha5
        [0.0, 0.0], #alpha6
        [0.0, 0.0], #rhoIyz
        [0.0, 0.0], #b
        [0.0, 0.0], #a0
        [0.0, 0.0], #aeroCenterOffset
    )
end

nothing

# We also want to provide some additional general inputs to enable floating simulation.
# `platformActive` being the principal and obvious one, but also the input files needed by
# HydroDyn and MoorDyn, which are the external libraries that calculate the hydrodynamic
# and mooring loads to send to the bottom side mesh. These input files include:
#  - interpOrder: the degree of interpolation used to predict the future states within
#    HydroDyn and MoorDyn, used in the `HD_UpdateStates` and `MD_UpdateStates` functions
#  - hd_input_file: the base HydroDyn input file with a `.dat` extension
#  - md_input_file: the base MoorDyn input file with a `.dat` extension
#  - ss_input_file: the sea state input file with a `.dat` extension used to define the environmental conditions of the sea.
#    This is used within HydroDyn.
#  - potfilefile: the prefix of the potential flow files.
#    At minimum, the directory must contain the .1, .3, and .hst WAMIT output files.
#    This is also used within HydroDyn.
# See the OpenFAST documentation for HydroDyn (https://openfast.readthedocs.io/en/main/source/user/hydrodyn/input_files.html)
# and MoorDyn (https://moordyn.readthedocs.io/en/latest/inputs.html) for more
# information about how to format these input files. For simplicity here, we will
# use predefined input files for the OC4 semisubmersible platform, which comes with
# OpenFAST and is copied to the `data` folder here.

if AeroModel=="AD"
    AD15On = true
else
    AD15On = false
end
platformActive = true
interpOrder = 2
hd_input_file = "$(path)/data/HydroDyn.dat"
md_input_file = "$(path)/data/MoorDyn.dat"
ss_input_file = "$(path)/data/SeaState.dat"
potflowfile = "$(path)/data/potflowdata/marin_semi"

inputs = OWENS.Inputs(;analysisType = analysisType,
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2,
platformActive,
interpOrder,
hd_input_file,
md_input_file,
ss_input_file,
potflowfile)

nothing

# We also will need to specify where the HydroDyn and MoorDyn libraries exist on the computer.
# We will set these to `nothing` so the auto-installed libraries will be used by default.
# Change these variables as needed if the location of these libraries are different.
hd_lib = nothing
md_lib = nothing
bin = OWENS.Bin(hd_lib, md_lib)

nothing

# As in previous examples, we need to define the inputs for the finite element models, though now we need
# a definition for both the topside and the bottom side.

# The changes to the topside from example A are the lack of boundary conditions,
# non-default gamma and alpha terms (for better convergence with the platform mesh).

topFEAModel = OWENS.FEAModel(;
    analysisType = analysisType,
    dataOutputFilename = "none",
    joint = topJoint,
    platformTurbineConnectionNodeNumber = 1,
    nlOn = true,
    numNodes = topMesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI",
    gamma = 1.0,
    alpha = 0.5)

nothing

# The bottom finite element mesh will use the same inputs with two exceptions: `joints` is now empty
# (since we're considering the platform as a single rigid body), and we add an additional concentrated term
# at the bottom node to represent the rigid body mass matrix of the platform.
# The values here are pulled from the report "Definition of the Semisubmersible Floating System for Phase II of OC4"
# by Robertson et al. (2014) (https://doi.org/10.2172/1155123).
numDOFPerNode = 6

ptfm_mass = 3.85218e6
ptfm_roll_iner = 2.56193e9
ptfm_pitch_iner = 2.56193e9
ptfm_yaw_iner = 4.24265e9
ptfm_Ixx = ptfm_roll_iner
ptfm_Iyy = ptfm_pitch_iner
ptfm_Izz = ptfm_yaw_iner
ptfm_Ixy = 0.0
ptfm_Iyz = 0.0
ptfm_Ixz = 0.0
ptfm_mass = [
    ptfm_mass 0.0 0.0 0.0 0.0 0.0
    0.0 ptfm_mass 0.0 0.0 0.0 0.0
    0.0 0.0 ptfm_mass 0.0 0.0 0.0
    0.0 0.0 0.0 ptfm_Ixx ptfm_Ixy ptfm_Ixz
    0.0 0.0 0.0 ptfm_Ixy ptfm_Iyy -ptfm_Iyz
    0.0 0.0 0.0 ptfm_Ixz -ptfm_Iyz ptfm_Izz
]
bottomConcInputs = [
    1 "M6" 1 1 ptfm_mass[1, 1]
    1 "M6" 1 2 ptfm_mass[1, 2]
    1 "M6" 1 3 ptfm_mass[1, 3]
    1 "M6" 1 4 ptfm_mass[1, 4]
    1 "M6" 1 5 ptfm_mass[1, 5]
    1 "M6" 1 6 ptfm_mass[1, 6]
    1 "M6" 2 1 ptfm_mass[2, 1]
    1 "M6" 2 2 ptfm_mass[2, 2]
    1 "M6" 2 3 ptfm_mass[2, 3]
    1 "M6" 2 4 ptfm_mass[2, 4]
    1 "M6" 2 5 ptfm_mass[2, 5]
    1 "M6" 2 6 ptfm_mass[2, 6]
    1 "M6" 3 1 ptfm_mass[3, 1]
    1 "M6" 3 2 ptfm_mass[3, 2]
    1 "M6" 3 3 ptfm_mass[3, 3]
    1 "M6" 3 4 ptfm_mass[3, 4]
    1 "M6" 3 5 ptfm_mass[3, 5]
    1 "M6" 3 6 ptfm_mass[3, 6]
    1 "M6" 4 1 ptfm_mass[4, 1]
    1 "M6" 4 2 ptfm_mass[4, 2]
    1 "M6" 4 3 ptfm_mass[4, 3]
    1 "M6" 4 4 ptfm_mass[4, 4]
    1 "M6" 4 5 ptfm_mass[4, 5]
    1 "M6" 4 6 ptfm_mass[4, 6]
    1 "M6" 5 1 ptfm_mass[5, 1]
    1 "M6" 5 2 ptfm_mass[5, 2]
    1 "M6" 5 3 ptfm_mass[5, 3]
    1 "M6" 5 4 ptfm_mass[5, 4]
    1 "M6" 5 5 ptfm_mass[5, 5]
    1 "M6" 5 6 ptfm_mass[5, 6]
    1 "M6" 6 1 ptfm_mass[6, 1]
    1 "M6" 6 2 ptfm_mass[6, 2]
    1 "M6" 6 3 ptfm_mass[6, 3]
    1 "M6" 6 4 ptfm_mass[6, 4]
    1 "M6" 6 5 ptfm_mass[6, 5]
    1 "M6" 6 6 ptfm_mass[6, 6]
]
bottomConcTerms = OWENSFEA.applyConcentratedTerms(
    bottomMesh.numNodes,
    numDOFPerNode,
    data=bottomConcInputs,
    jointData=[])
bottomFEAModel = OWENS.FEAModel(;
    analysisType = analysisType,
    dataOutputFilename = "none",
    joint = [],
    platformTurbineConnectionNodeNumber = 1,
    nlOn = true,
    numNodes = bottomMesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI",
    gamma = 1.0,
    alpha = 0.5,
    nodalTerms=bottomConcTerms)

nothing

# Finally, some additional structure definitions for the bottom mesh.
# These don't contain anything new, but are required for OWENS to run properly.
bottomEl = OWENSFEA.El(bottomSectionProps,
                       bottomOrt.Length,
                       bottomOrt.Psi_d,
                       bottomOrt.Theta_d,
                       bottomOrt.Twist_d,
                       ones(bottomMesh.numEl))

nothing

# We can now run our floating simulation in OWENS.
# Note that we are using `OWENS.Unsteady` now instead of `OWENS.Unsteady_Land`.

println("Running Unsteady")
t, aziHist, OmegaHist, OmegaDotHist, gbHist, gbDotHist, gbDotDotHist, FReactionHist, FTwrBsHist,
genTorque, genPower, torqueDriveShaft, uHist, uHist_prp,
epsilon_x_hist, epsilon_y_hist, epsilon_z_hist, kappa_x_hist, kappa_y_hist, kappa_z_hist,
FPtfmHist, FHydroHist, FMooringHist = OWENS.Unsteady(inputs,
    system=topSystem,
    assembly=topAssembly,
    topModel=topFEAModel,
    topMesh=topMesh,
    topEl=topEl,
    aero=aeroForces,
    deformAero=deformAero,
    bottomModel=bottomFEAModel,
    bottomMesh=bottomMesh,
    bottomEl=bottomEl,
    bin=bin)

if AD15On #TODO: move this into the run functions
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

nothing

# Since we ran the simulation for such a short time period for sake of demonstrating
# completion of the simulation, the plots below are not too interesting.
# Try experimenting with changing the `numTS` variable above to a longer simulation to
# see how the floating VAWT responds over time.
# Note that the platform modeled here is not necessarily design for the Darrieus VAWT it
# is supporting, so your mileage may vary on the simulation results.
# PyPlot.plot(t, uHist_prp[:, 1])
# PyPlot.title("Platform Surge Response")
# PyPlot.xlabel("Time [s]")
# PyPlot.ylabel("Force [N]")
# PyPlot.plot(t, FHydroHist[:, 1])
# PyPlot.xlabel("Time [s]")
# PyPlot.ylabel("Force [N]")
# PyPlot.title("Hydrodynamic Surge Loading")

nothing