import PyPlot
PyPlot.pygui(true)
PyPlot.close("all")
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

import Statistics:mean
import DelimitedFiles
import Dierckx
import QuadGK
import FLOWMath
import HDF5

import ModelGen
import GyricFEA
import OWENS
import VAWTAero
import Composites

path,_ = splitdir(@__FILE__)

# include("$(path)/../../../../OWENS.jl/src/OWENS.jl")
# include("$(path)/../../../../ModelGen.jl/src/ModelGen.jl")
println("Set up Macro Geometry/Inputs")
rho = 1.225
Nslices = 30
ntheta = 30
RPM = 7.2
Vinf = 22.26786 * 0.3048
slc1 = 3
slc2 = 5
eta = 0.5
B = Nbld = 2
R = 177.2022*0.3048 #m
H = 1.02*R*2 #m
chord = 5.0*ones(Nslices)
omega = RPM / 60 * 2 * pi
tsr = omega*R/Vinf

shapeY = collect(LinRange(0,H,Nslices+1))
shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)
shapeX_spline = FLOWMath.Akima(shapeY, shapeX)
h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
h_elem = (shapeY[2:end] - shapeY[1:end-1])
h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
RefArea = RefArea_half*2
delta_xs = shapeX[2:end] - shapeX[1:end-1]
delta_zs = shapeY[2:end] - shapeY[1:end-1]

delta3D = atan.(delta_xs./delta_zs)

#########################################
### Set up aero forces
#########################################
println("Initialize Aerodynamics")
VAWTAero.setupTurb(shapeX,shapeY,B,chord,tsr,Vinf;AModel="DMS",DSModel="BV",
afname = "$(path)/airfoils/NACA_0021.dat",
ifw=true,
ifw_libfile = joinpath("$(path)/../../openfast/build/modules/inflowwind/libifw_c_binding"),
turbsim_filename="$(path)/data/300mx300m12msETM_Coarse.bts",
ntheta,Nslices,rho,eta,RPI=true)

#########################################
### Set up mesh
#########################################
println("Create Mesh")
mymesh,myort,myjoint = ModelGen.create_mesh_struts(;Ht=15.0,
Hb = H, #blade height
R, # m bade radius
nblade = 2,
ntelem = 30, #tower elements
nbelem = 60, #blade elements
nselem = 10,
strut_mountpointbot = 0.1,
strut_mountpointtop = 0.1,
bshapex = bshapex=shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
bshapez = shapeY,
angularOffset = -pi/2) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

#########################################
### Set up Sectional Properties
#########################################
println("Calculate/Set up sectional properties")
#Tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn_twr = ModelGen.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

#Add the full path
for (i,airfoil) in enumerate(numadIn_twr.airfoil)
    numadIn_twr.airfoil[i] = "$path/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops_twr = ModelGen.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = ModelGen.getPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
nTwrElem = Int(mymesh.meshSeg[1])+1
sectionPropsArray_twr = ModelGen.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
stiff_twr, mass_twr = ModelGen.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

#Blades
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = ModelGen.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops_bld = ModelGen.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# Get blade spanwise position
bld1start = Int(mymesh.structuralNodeNumbers[1,1])
bld1end = Int(mymesh.structuralNodeNumbers[1,end])
spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = ModelGen.getPreCompOutput(numadIn_bld;plyprops = plyprops_bld)
sectionPropsArray_bld = ModelGen.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
stiff_bld, mass_bld = ModelGen.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
Nremain = mymesh.numEl-length(sectionPropsArray_twr)-length(bldssecprops) #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;bldssecprops;fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = ModelGen.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld;VTKmeshfilename="$path/vtk/SNL5MW")

#########################################
### Create Aero Functions
#########################################

# # Example aero forces file input
# d1 = "$path/data/Simplified5MW2.geom"
# d2 = "$path/data/Simplified5MW2_ElementData.csv"
# d3 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.bld"
# d4 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.el"
# d5 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.ort"
# d6 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.mesh"
#
# aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = OWENS.mapCactusLoadsFile(d1,d2,d3,d4,d5,d6)
# aeroForcesCACTUS(t,azi) = OWENS.externalForcing(t,aerotimeArray,aeroForceValHist,aeroForceDof)

# aeroForcesDMSFile(t) = VAWTAero.mapCACTUSFILE_minimalio(t,mymesh,myel,"$path/data/DMS_Simplified5MW_ElementData.csv")
# Example realtime aero calculation setup
aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.AdvanceTurbineInterpolate;alwaysrecalc=true)

######################################
#### Perform Aerostructural One Way Test
#######################################

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

model = OWENS.Inputs(;analysisType = "ROM",
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS = 40.0,
delta_t = 0.05,
aeroLoadsOn = 2)

feamodel = GyricFEA.FEAModel(;analysisType = "ROM",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
numNodes = mymesh.numNodes,
RayleighAlpha = 0.0,
RayleighBeta = 0.0,
iterationType = "DI")

# Choose which aeroforces
# aeroForces = aeroForcesCACTUS
aeroForces = aeroForcesDMS

# deformaero3(x;newOmega=-1,newVinf=-1) = 0

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist = OWENS.Unsteady(model;
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero=VAWTAero.deformTurb,system,assembly)

println("Saving VTK time domain files")
ModelGen.gyricFEA_VTK("$path/vtk/SNL5MW_timedomain",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

##########################################
#### Get strain values at the blades #####
##########################################

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL = OWENS.extractSF(bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
kappa_x_hist,epsilon_y_hist;verbosity=1, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor
# epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
Twr_LE_U_idx=1,Twr_LE_L_idx=1)

######################################
#### Plot
#######################################

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,1])
PyPlot.ylabel("FReaction Hist 1")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction1.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,2])
PyPlot.ylabel("FReaction Hist 2")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction2.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,3])
PyPlot.ylabel("FReaction Hist 3")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction3.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,4])
PyPlot.ylabel("FReaction Hist 4")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction4.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,5])
PyPlot.ylabel("FReaction Hist 5")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction5.pdf",transparent = true)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,6])
PyPlot.ylabel("FReaction Hist 6")
# PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction6.pdf",transparent = true)

# PyPlot.figure()
# for ii = 1:length(uHist[1,:])
# PyPlot.plot(1:length(old_uHist[ii,:]),old_uHist[ii,:],"k-")
#     PyPlot.plot(1:length(uHist[ii,:]),uHist[ii,:],"k--")
#     if ii%10 == 0.0
#         # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_Uhist$ii.pdf",transparent = true)
#         PyPlot.figure()
#     end
# end
