
import OWENS
import OWENSAero
import DelimitedFiles
using Statistics:mean
using Test

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# function runprofilefunction()
path = runpath = splitdir(@__FILE__)[1]

# Unpack inputs, or you could directly input them here and bypass the file 

verbosity = 1


analysisType = "unsteady"
turbineType = "H-VAWT"
eta = 0.5
Nbld = 2
towerHeight = 0.5
Vinf = 1.2
controlStrategy = "constantRPM"
Nslices = 20
ntheta = 30
structuralModel = "TNB"
ntelem = 20
nbelem = 60
ncelem = 10
nselem = 10
ifw = false
AModel = "DMS"
windINPfilename = "$path/300mx300m12msETM_Coarse.bts"
ifw_libfile = nothing#"$path/../../openfast/build/modules/inflowwind/libifw_c_binding"
Blade_Height = 20.0
Blade_Radius = 4.5
area = Blade_Height*2*Blade_Radius
numTS = 75
delta_t = 0.005
NuMad_geom_xlscsv_file_twr = "$path/TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$path/TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$path/GeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$path/Materials34m.csv"
NuMad_geom_xlscsv_file_strut = "$path/GeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$path/Materials34m.csv"
adi_lib = nothing#"$path/../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" 
adi_rootname = "$path/helical"

# TSR = 3.0

# omega = Vinf/Blade_Radius*TSR
RPM = 1e-6#omega / (2*pi) * 60

fluid_density = 998.0
fluid_dyn_viscosity = 1.792E-3
number_of_blades = Nbld
WindType = 3

AM_flag = true
rotAccel_flag = false
buoy_flag = false

##############################################
# Setup
#############################################

tocp = [0.0;10.0; 1e6]
Omegaocp = [RPM; RPM; RPM]./60 #control inputs

tocp_Vinf = [0.0;10.0; 1e6]
Vinfocp = [Vinf;Vinf;Vinf]

shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
helix_angle = 0.0#-pi/4
shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius#SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ) ones(length(shapeZ)).*Blade_Radius#
shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius # zeros(length(shapeX))#

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
    H = Blade_Height, #windio
    R = Blade_Radius, #windio
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
    strut_twr_mountpoint = [0.5],
    strut_bld_mountpoint = [0.5],
    AModel, #AD, DMS, AC
    DSModel="BV",
    AM_flag,
    rotAccel_flag,
    buoy_flag,
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)



PyPlot.figure()
for icon = 1:length(mymesh.conn[:,1])
    idx1 = mymesh.conn[icon,1]
    idx2 = mymesh.conn[icon,2]
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    # sleep(0.1)
end

for ijoint = 1:length(myjoint[:,1])
    idx2 = Int(myjoint[ijoint,2])
    idx1 = Int(myjoint[ijoint,3])
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    sleep(0.1)
end

PyPlot.scatter3D(1,1,1,"b.")

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
1 6 0
24 1 0
24 2 0
24 3 0
24 4 0
24 5 0
24 6 0
87 1 0
87 2 0
87 3 0
87 4 0
87 5 0
87 6 0
146 1 0
146 2 0
146 3 0
146 4 0
146 5 0
146 6 0
]

nothing

# There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

if AModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;analysisType = structuralModel,
tocp,
Omegaocp,
tocp_Vinf,
Vinfocp,
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2)

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
gravityOn = [0,0,0.0],
numNodes = mymesh.numNodes,
RayleighAlpha = 0.00,
RayleighBeta = 0.00,
iterationType = "DI")

nothing

# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.
forced_node = 85
function flappingForces(t,azi)
    if t<0.1
        Fexternal = [-1000000.0]
        Fdof = [(forced_node-1)*6+1]
    else
        Fexternal,Fdof = aeroForces(t,azi)
    end
    return Fexternal, Fdof
end

Keff = 6e5 #N/m
L = 9.5 #m 
# K = 3EI/L^3
EI = 2e8 #this is what was in the paper, which is not the same as the equation #Keff*L^3/3

# Set the element properties
for iprop = 1:length(myel.props)
    # do just the non-super stiff elements
    if myel.props[iprop].EIyy[2] < 1e9
        myel.props[iprop].EIyy .= EI
        myel.props[iprop].EIzz .= EI
        myel.props[iprop].GJ .= EI*5
        myel.props[iprop].EA .= EI*100
    end
end



println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=flappingForces,deformAero,verbosity=5)


nothing

OWENS_tip_displ = uHist[:,(forced_node-1)*6+1]

if AM_flag
    omega_OF = 4.5105 * 2*pi
else
    omega_OF = 13.53* 2*pi
end

ofast_tdispl = cos.(t.*omega_OF)

PyPlot.figure("important")
PyPlot.plot(t,ofast_tdispl,label="OpenFAST Added Mass $AM_flag")
PyPlot.plot(t,OWENS_tip_displ,label="OWENS Added Mass $AM_flag")
PyPlot.legend()


# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

azi=aziHist#./aziHist*1e-6
saveName = "$path/vtk/flapping_added_mass2"
tsave_idx=1:1:numTS-1
OWENS.OWENSVTK(saveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
    epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FReactionHist,topFexternal_hist;tsave_idx)

nothing

# This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
# it also calculates analytical simply supported buckling safety factors

##########################################
#### Ultimate Failure #####
##########################################

# massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
# SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
# topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
# topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
# bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
# twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
# mymesh,myel,myort,number_of_blades,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
# kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
# LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
# LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
# Twr_LE_U_idx=1,Twr_LE_L_idx=1,
# AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

# ffmpeg -i smallerhelical.%04d.png -vf palettegen=reserve_transparent=1 palette.png
# ffmpeg -framerate 30 -i smallerhelical.%04d.png -i palette.png -lavfi paletteuse=alpha_threshold=30 -gifflags -offsetting smallerhelical.gif
