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

import OWENSFEA
import OWENS
import OWENSAero
import Composites

path,_ = splitdir(@__FILE__)

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
R = Blade_Radius = 177.2022*0.3048 #m
H = Blade_Height = 1.02*R*2 #m
towerHeight = 15.0
chord = 5.0*ones(Nslices)
omega = RPM / 60 * 2 * pi
tsr = omega*R/Vinf
ifw=false
delta_t = 0.05
simtime = 6.0
numTS = simtime/delta_t

AModel = "DMS"
turbineType = "Darrieus"

if AModel=="AD" #TODO: unify flag
    AD15On=true #AD for AeroDyn, DMS for double multiple streamtube, AC for actuator cylinder
else
    AD15On=false
end

ntelem = 30 #tower elements
nbelem = 60 #blade elements
nselem = 10

adi_lib=nothing
adi_rootname="$path/SNL5MW"
windINPfilename = "$(path)/data/300mx300m12msETM_Coarse.bts"
ifw_libfile = nothing

NuMad_geom_xlscsv_file_twr = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
NuMad_mat_xlscsv_file_twr = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
NuMad_geom_xlscsv_file_bld = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkinDMS.csv"
NuMad_mat_xlscsv_file_bld = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
NuMad_geom_xlscsv_file_strut = ["$path/data/NuMAD_Geom_SNL_5MW_strutsDMS.csv","$path/data/NuMAD_Geom_SNL_5MW_strutsDMS.csv"]
NuMad_mat_xlscsv_file_strut = NuMad_mat_xlscsv_file_bld

shapeZ = collect(LinRange(0,H,Nslices+1))
shapeX = R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)
shapeX_spline = FLOWMath.Akima(shapeZ, shapeX)
h_frac = (shapeZ[2:end] - shapeZ[1:end-1])./shapeZ[end];
h_elem = (shapeZ[2:end] - shapeZ[1:end-1])
h = (shapeZ[2:end] + shapeZ[1:end-1])/2.0;

RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
RefArea = RefArea_half*2
delta_xs = shapeX[2:end] - shapeX[1:end-1]
delta_zs = shapeZ[2:end] - shapeZ[1:end-1]

delta3D = atan.(delta_xs./delta_zs)

# #########################################
# ### Set up aero forces
# #########################################
# println("Initialize Aerodynamics")
# OWENSAero.setupTurb(shapeX,shapeZ,B,chord,tsr,Vinf;AModel="DMS",DSModel="BV",
# afname = "$(path)/airfoils/NACA_0021.dat",
# ifw=false,
# ifw_libfile = nothing,
# turbsim_filename="$(path)/data/300mx300m12msETM_Coarse.bts",
# ntheta,Nslices,rho,eta,RPI=true)

# #########################################
# ### Set up mesh
# #########################################
# println("Create Mesh")
# mymesh,myort,myjoint = OWENS.create_mesh_struts(;Htwr_base=15.0,
# Hbld = H, #blade height
# R, # m bade radius
# Htwr_blds=H,
# nblade = 2,
# ntelem = 30, #tower elements
# nbelem = 60, #blade elements
# nselem = 10,
# strut_twr_mountpoint = [0.1,0.9], # This puts struts at top and bottom
# strut_bld_mountpoint = [0.1,0.9], # This puts struts at top and bottom
# bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
# bshapez = shapeZ,
# angularOffset = -pi/2) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

# PyPlot.figure()
# for icon = 1:length(mymesh.conn[:,1])
#     idx1 = mymesh.conn[icon,1]
#     idx2 = mymesh.conn[icon,2]
#     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
#     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
#     # sleep(0.1)
# end

# for ijoint = 1:length(myjoint[:,1])
#     idx2 = Int(myjoint[ijoint,2])
#     idx1 = Int(myjoint[ijoint,3])
#     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
#     sleep(0.1)
#     PyPlot.scatter3D(1,1,1,".")
# end

# #########################################
# ### Set up Sectional Properties
# #########################################
# println("Calculate/Set up sectional properties")
# #Tower
# NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
# numadIn_twr = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

# #Add the full path
# for (i,airfoil) in enumerate(numadIn_twr.airfoil)
#     numadIn_twr.airfoil[i] = "$path/airfoils/$airfoil"
# end

# NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
# plyprops_twr = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = OWENS.getOWENSPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
# nTwrElem = Int(mymesh.meshSeg[1])+1
# sectionPropsArray_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
# stiff_twr, mass_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

# #Blades
# NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
# numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

# for (i,airfoil) in enumerate(numadIn_bld.airfoil)
#     numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
# end

# NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
# plyprops_bld = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# # Get blade spanwise position
# bld1start = Int(mymesh.structuralNodeNumbers[1,1])
# bld1end = Int(mymesh.structuralNodeNumbers[1,end])
# spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

# bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getOWENSPreCompOutput(numadIn_bld;plyprops = plyprops_bld)
# sectionPropsArray_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
# stiff_bld, mass_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

# #Struts
# # They are the same as the end properties of the blades

# # Combined Section Props
# bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
# Nremain = mymesh.numEl-length(sectionPropsArray_twr)-length(bldssecprops) #strut elements remain
# strutssecprops = fill(sectionPropsArray_bld[end],Nremain)

# sectionPropsArray = [sectionPropsArray_twr; bldssecprops; strutssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
        
# # GXBeam sectional properties
# stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
# stiff_struts = fill(stiff_bld[end],Nremain)#collect(Iterators.flatten(fill(stiff_strut, Nstrutperbld*Nbld)))
# stiff_array = [stiff_twr; stiff_blds; stiff_struts]

# mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
# mass_struts = fill(mass_bld[end],Nremain)#collect(Iterators.flatten(fill(mass_strut, Nstrutperbld*Nbld)))
# mass_array = [mass_twr; mass_blds; mass_struts]

# rotationalEffects = ones(mymesh.numEl)

# #store data in element object
# myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

# println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
# system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_array,mass_array)

# #########################################
# ### Create Aero Functions
# #########################################

# # # Example aero forces file input
# # d1 = "$path/data/Simplified5MW2.geom"
# # d2 = "$path/data/Simplified5MW2_ElementData.csv"
# # d3 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.bld"
# # d4 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.el"
# # d5 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.ort"
# # d6 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.mesh"
# #
# # aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = OWENS.mapCactusLoadsFile(d1,d2,d3,d4,d5,d6)
# # aeroForcesCACTUS(t,azi) = OWENS.externalForcing(t,aerotimeArray,aeroForceValHist,aeroForceDof)

# # aeroForcesDMSFile(t) = OWENSAero.mapCACTUSFILE_minimalio(t,mymesh,myel,"$path/data/DMS_Simplified5MW_ElementData.csv")
# # Example realtime aero calculation setup
# aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)


mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    RPM,
    Vinf,
    eta,
    B = Nbld,
    H = Blade_Height,
    R = Blade_Radius,
    shapeZ,
    shapeX,
    ifw,
    delta_t,
    numTS,
    adi_lib,
    adi_rootname,
    AD15hubR = 0.0,
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
    # ncelem,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.1,0.9],
    strut_bld_mountpoint = [0.1,0.9],
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = pi/2,
    meshtype = turbineType)


######################################
#### Perform Aerostructural One Way Test
#######################################

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

model = OWENS.Inputs(;analysisType = "GX",
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS,
delta_t,
aeroLoadsOn = 2)

feamodel = OWENSFEA.FEAModel(;analysisType = "GX",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
numNodes = mymesh.numNodes,
RayleighAlpha = 0.0,
RayleighBeta = 0.0,
iterationType = "DI")

# deformaero3(x;newOmega=-1,newVinf=-1) = 0

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
topFexternal_hist,rbDataHist = OWENS.Unsteady(model;
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,system,assembly)

saveName = "$path/vtk/SNL5MW_timedomain"
OWENS.OWENSVTK(saveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FReactionHist,topFexternal_hist)

##########################################
#### Get strain values at the blades #####
##########################################

# massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
# SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL = OWENS.extractSF(bld_precompinput,
# bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
# twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
# mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
# kappa_x_hist,epsilon_y_hist;verbosity=1, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor
# # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
# LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
# LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
# Twr_LE_U_idx=1,Twr_LE_L_idx=1)

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
