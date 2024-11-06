using Statistics:mean
using Test
import HDF5
import PyPlot
import DelimitedFiles
import QuadGK
import FLOWMath
import OWENSFEA
import Composites
import OWENS
import OWENSAero
import RollingFunctions
import GXBeam

steady = false

path = splitdir(@__FILE__)[1]

t_offset = 5.0 #sec

verbosity = 1

turbineType = "Darrieus"
eta = 0.5
Nbld = 2
towerHeight = 0.5
rho = 0.94
Nslices = 35
ntheta = 30
ntelem = 20
nbelem = 60
ncelem = 10
nselem = 5
ifw = false
AeroModel = "DMS"
windINPfilename = nothing
ifw_libfile = nothing
Blade_Height = 41.9
Blade_Radius = 17.1
numTS = 2000
delta_t = 0.05
NuMad_geom_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mGeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mMaterials.csv"
NuMad_geom_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mGeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mMaterials.csv"
adi_lib = nothing
adi_rootname = "$(path)/SNL34m"

println("Running Gravity Loading")


controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = safeakima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
SNL34_unit_xz = [x_shape;;z_shape]
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])
SNL34Z = SNL34z.*Blade_Height
SNL34X = SNL34x.*Blade_Radius

shapeZ = SNL34Z#collect(LinRange(0,H,Nslices+1))
shapeX = SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
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
    ncelem,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.03,0.97],
    strut_bld_mountpoint = [0.03,0.97],
    AeroModel, #AD, DMS, AC
    DynamicStallModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = pi/2,
    meshtype = turbineType)


# thickness_flap is distance from shear center x to top
# thickness_lag is distance from shear center y to trailing edge
# shear center is relative to the blade reference axis
# blade reference axis is from leading edge to le_loc*chord and along the chord line
# reference axes Y is along the chord, and X is perpendicular to the chord

thickness_precomp_lag = zeros(length(bld_precompinput))
thickness_precomp_flap = zeros(length(bld_precompinput))
for ipc = 1:length(bld_precompinput)
    refY = bld_precompinput[ipc].le_loc*bld_precompinput[ipc].chord
                                # Negative distance for lag, to align with SAND-88-1144
    thickness_precomp_lag[ipc] = -(bld_precompinput[ipc].chord-(refY+bld_precompoutput[ipc].y_sc))
    thickness_precomp_flap[ipc] = maximum(bld_precompinput[ipc].ynode)*bld_precompinput[ipc].chord - bld_precompoutput[ipc].x_sc
end
bld1start = Int(mymesh.structuralNodeNumbers[1,1])
bld1end = Int(mymesh.structuralNodeNumbers[1,end])
spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]
spanposmid = cumsum(diff(spanpos))
thickness = safeakima(numadIn_bld.span,thickness_precomp_flap,spanposmid)
thickness_lag = safeakima(numadIn_bld.span,thickness_precomp_lag,spanposmid)
# thickness = thicknessGX[1:end-1]

PyPlot.figure()
PyPlot.plot(mymesh.x,mymesh.z,"b-")
    for myi = 1:length(mymesh.x)
        PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
        PyPlot.draw()
        #sleep(0.1)
    end
PyPlot.xlabel("x")
PyPlot.ylabel("y")

top_idx = ntelem+1+Nbld#Int(myjoint[7,2])
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
top_idx 1 0
top_idx 2 0
top_idx 3 0
top_idx 4 0
top_idx 5 0]

model = OWENS.Inputs(;verbosity,analysisType = "TNB",
outFilename = "none",
tocp = [0.0, 1e6],#SNL34m_5_3_RPM[:,1],#[0.0,10.0,100000.1],
Omegaocp = [0.0, 0.0]./60,#SNL34m_5_3_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
tocp_Vinf = [0.0, 1e6],
Vinfocp = [0.0, 0.0],
numTS = 600,
delta_t = 0.03,
turbineStartup = 0,
aeroLoadsOn = 2,
generatorOn = false,
ratedTorque = 100000.0,
OmegaInit = 0.0/60,
zeroTorqueGenSpeed = 0.010,
pulloutRatio = 0.90,
ratedGenSlipPerc = 9000.0)

feamodel = OWENSFEA.FEAModel(;analysisType = "TNB",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn=true,
numNodes = mymesh.numNodes,
numModes = 200,
RayleighAlpha = 0.1,
RayleighBeta = 0.1,
gravityOn = true,
# nodalTerms = OWENSFEA.readNodalTerms(;data=[94 "F6" 2 2 1e6;145 "F6" 2 2 1e6]),
iterationType = "DI")

deformTurb(azi_j;newOmega=0,newVinf=0,bld_x=0,bld_z=0,bld_twist=0) = 0

eps_xG,eps_zG,eps_yG,kappa_xG,kappa_yG,kappa_zG,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist = OWENS.run34m(model,feamodel,mymesh,myel,
aeroForces,deformTurb;steady)

##################################################################
########### FIG 4.1 Gravity Only #############
##################################################################

# Load data
experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_bothbladesreversed.csv",',',skipstart = 0)
predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1 = (kappa_yG[1,end-1,2:end] .* thickness .+ 0*eps_xG[1,end-1,2:end]) .* Ealuminum
flatwise_stress2 = (kappa_yG[2,end-1,1:end-1] .* thickness .+ 0*eps_xG[2,end-1,1:end-1]) .* Ealuminum

edgewise_stress1 = (kappa_zG[1,end-1,2:end] .* thickness .+ 0*eps_xG[1,end-1,2:end]) .* Ealuminum

using PyPlot
PyPlot.pygui(true)
PyPlot.close("all")
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

PyPlot.figure("Grav")
PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "VAWT Test Data")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1[2:end-5]./1e6,spanposmid[2:end-5].+toweroffset,".-",color=plot_cycle[1],label = "OWENS Timoshenko Beam")
PyPlot.plot(flatwise_stress2[2:end-5]./1e6,spanposmid[2:end-5].+toweroffset,".-",color=plot_cycle[3],label = "OWENS Timoshenko Beam blade2")
# PyPlot.plot(spanposmid),-flatwise_stress2./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
PyPlot.legend(loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig4_1_GravityOnly_flapwise_BladeGXNL.pdf",transparent = true)


PyPlot.figure(11)
# PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(edgewise_stress1[1:end]./1e6,spanposmid[1:end].+toweroffset,".-",color=plot_cycle[1],label = "OWENS Timoshenko Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Edgewise Stress (MPa)")
PyPlot.legend(loc = (0.06,1.0),ncol=2)

##################################################################
########### FIG 4.3 28RPM with Gravity Loads Tared out  ##########
##################################################################

model.Omegaocp = [28.0, 28.0]./60
model.OmegaInit = 28.0/60

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist  = OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformTurb;steady)


# Load data
experimentalCent28 = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.3_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1_centrifugal28 = (kappa_y[1,end-1,2:end].* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum
flatwise_stress2_centrifugal28 = (kappa_y[2,end-1,1:end-1].* thickness .+ 0*eps_x[2,end-1,1:end-1]) .* Ealuminum

# Zero out gravity loads
flatwise_stress1_centrifugal28 .-= flatwise_stress1
flatwise_stress2_centrifugal28 .-= flatwise_stress2

PyPlot.figure("Cent28")
PyPlot.plot(experimentalCent28[:,2],experimentalCent28[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1_centrifugal28[5:end-5]./1e6,spanposmid[5:end-5].+toweroffset,".-",color=plot_cycle[1],label = "OWENS Timoshenko Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2_centrifugal28./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
PyPlot.legend(loc = (0.06,0.95),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig4_3_Centrifugal_28RPM_flapwise_BladeGXNL.pdf",transparent = true)



println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","Fx","Fy","Fz","Mx","My","Mz"]
# userPointData[iname,it,ipt] = Float64

# map el props to points using con
userPointData = zeros(length(userPointNames),length(t),mymesh.numNodes)
EA_points = zeros(mymesh.numNodes)
EIyy_points = zeros(mymesh.numNodes)
EIzz_points = zeros(mymesh.numNodes)

# Time-invariant data
for iel = 1:length(myel.props)
    # iel = 1
    nodes = mymesh.conn[iel,:]
    EA_points[Int.(nodes)] = myel.props[iel].EA
    EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
    EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
end

# fill in the big matrix
for it = 1:length(t)

    userPointData[1,it,:] = EA_points
    userPointData[2,it,:] = EIyy_points
    userPointData[3,it,:] = EIzz_points
    # userPointData[4,it,:] = FReactionHist[it,1:6:end]
    # userPointData[5,it,:] = FReactionHist[it,2:6:end]
    # userPointData[6,it,:] = FReactionHist[it,3:6:end]
    # userPointData[7,it,:] = FReactionHist[it,4:6:end]
    # userPointData[8,it,:] = FReactionHist[it,5:6:end]
    # userPointData[9,it,:] = FReactionHist[it,6:6:end]
end

azi=aziHist#./aziHist*1e-6
VTKsaveName = "$path/vtk/two_blade"
OWENS.OWENSFEA_VTK(VTKsaveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)

##################################################################
########### FIG 4.4 40RPM with Gravity Loads Tared out  ##########
##################################################################

model.Omegaocp = [40.0, 40.0]./60
model.OmegaInit = 40.0/60

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist= OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformTurb;steady,system,assembly)

# println("writing")
# OWENS.OWENSFEA_VTK("SNL34m_test",t,uHist,system,assembly,sections;scaling=10)#,azi=aziHist)
# println("done")

# Load data
experimentalCent = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.4_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1_centrifugal40 = (kappa_y[1,end-1,2:end].* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum
flatwise_stress2_centrifugal40 = (kappa_y[2,end-1,1:end-1].* thickness .+ 0*eps_x[2,end-1,1:end-1]) .* Ealuminum
edgewise_stress1_centrifugal40 = (kappa_z[1,end-1,2:end] .* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum



# Zero out gravity loads
flatwise_stress1_centrifugal40 .-= flatwise_stress1
flatwise_stress2_centrifugal40 .-= flatwise_stress2
edgewise_stress1_centrifugal40 .-= edgewise_stress1

PyPlot.figure("Cent40")
PyPlot.plot(experimentalCent[:,2],experimentalCent[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1_centrifugal40[5:end-5]./1e6,spanposmid[5:end-5].+toweroffset,".-",color=plot_cycle[1],label = "OWENS Timoshenko Beam")
PyPlot.plot(flatwise_stress2_centrifugal40[5:end-5]./1e6,spanposmid[5:end-5].+toweroffset,".-",color=plot_cycle[3],label = "OWENS Timoshenko Beam2")
# PyPlot.plot(spanposmid),-flatwise_stress2_centrifugal./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
PyPlot.legend(loc = (0.06,0.95),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig4_4_Centrifugal_40RPM_flapwise_BladeGXNL.pdf",transparent = true)

PyPlot.figure("Cent40lag")
# PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(edgewise_stress1_centrifugal40[1:end-4]./1e6,spanposmid[1:end-4].+toweroffset,".-",color=plot_cycle[1],label = "OWENS Timoshenko Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Edgewise Stress (MPa)")
PyPlot.legend(loc = (0.06,1.0),ncol=2)

PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,6])

############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################
############### GX ################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################

# Vinf = 1e-2#mean(SNL34m_5_3_Vinf[:,2])
# TSR = 0.001
# ntheta = 200#176

# OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
#     eta = 0.25,
#     rho,
#     mu = 1.7894e-5,
#     ntheta,
#     Nslices,
#     ifw = false,
#     RPI = true,
#     DynamicStallModel = "BV",
#     AeroModel = "DMS",
#     tau = [1e-5,1e-5],
#     afname = airfoils)

# dt = 1/(RPM/60*ntheta)

# aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.advanceTurb)
# deformTurb(azi_j;newOmega=0,newVinf=0,bld_x=0,bld_z=0,bld_twist=0) = 0

model = OWENS.Inputs(;analysisType = "GX",
outFilename = "none",
tocp = [0.0, 1e6],#SNL34m_5_3_RPM[:,1],#[0.0,10.0,100000.1],
Omegaocp = [0.0, 0.0]./60,#SNL34m_5_3_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
tocp_Vinf = [0.0, 1e6],
Vinfocp = [0.0, 0.0],
numTS = 600,
delta_t = 0.03,
turbineStartup = 0,
aeroLoadsOn = 2,
generatorOn = false,
ratedTorque = 100000.0,
OmegaInit = 0.0/60,
zeroTorqueGenSpeed = 0.010,
pulloutRatio = 0.90,
ratedGenSlipPerc = 9000.0)

feamodel = OWENSFEA.FEAModel(;analysisType = "GX",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn=false,
numNodes = mymesh.numNodes,
numModes = 200,
RayleighAlpha = 0.1,
RayleighBeta = 0.1,
gravityOn = true,
# nodalTerms = OWENSFEA.readNodalTerms(;data=[94 "F6" 2 2 1e6;145 "F6" 2 2 1e6]),
iterationType = "DI")


eps_xG,eps_zG,eps_yG,kappa_xG,kappa_yG,kappa_zG,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist = OWENS.run34m(model,feamodel,mymesh,myel,
aeroForces,deformTurb;steady,system,assembly)

##################################################################
########### FIG 4.1 Gravity Only #############
##################################################################

# Load data
experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_bothbladesreversed.csv",',',skipstart = 0)
predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1GX = (kappa_yG[1,end-1,2:end] .* thickness .+ 0*eps_xG[1,end-1,2:end]) .* Ealuminum
flatwise_stress2GX = (kappa_yG[2,end-1,1:end-1] .* thickness .+ 0*eps_xG[2,end-1,1:end-1]) .* Ealuminum

edgewise_stress1GX = (kappa_zG[1,end-1,2:end] .* thickness .+ 0*eps_xG[1,end-1,2:end]) .* Ealuminum

# # Test Gravity
# for ipt = 1:length(flatwise_stress1GX)
#     atol = abs(flatwise_stress1GX[ipt]*0.05)
#     println( isapprox(flatwise_stress1GX[ipt],flatwise_stress1[ipt];atol))
# end

PyPlot.figure("Grav")
# PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1GX[2:end-5]./1e6,spanposmid[2:end-5].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2GX./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
# PyPlot.legend(loc = (0.06,1.0),ncol=2)
PyPlot.legend(loc = (.0,0.83))
# PyPlot.savefig("$(path)/../figs/34m_fig4_1_GravityOnly_flapwise_BladeGXNL.pdf",transparent = true)


PyPlot.figure()
# PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
# PyPlot.plot(edgewise_stress1[0.8[8[30--]]][1:end]./1e6,spanposmid[1:end].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2GX./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Edgewise Stress (MPa)")
PyPlot.legend(loc = (0.06,1.0),ncol=2)

##################################################################
########### FIG 4.3 28RPM with Gravity Loads Tared out  ##########
##################################################################

model.Omegaocp = [28.0, 28.0]./60
model.OmegaInit = 28.0/60

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist  = OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformTurb;steady,system,assembly)


# Load data
experimentalCent28 = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.3_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1GX_centrifugal28 = (kappa_y[1,end-1,2:end].* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum
flatwise_stress2GX_centrifugal28 = (kappa_y[2,end-1,1:end-1].* thickness .+ 0*eps_x[2,end-1,1:end-1]) .* Ealuminum

# Zero out gravity loads
flatwise_stress1GX_centrifugal28 .-= flatwise_stress1GX
flatwise_stress2GX_centrifugal28 .-= flatwise_stress2GX

# # Test 28 RPM
# for ipt = 1:length(flatwise_stress1GX_centrifugal28)
#     atol = abs(flatwise_stress1GX_centrifugal28[ipt]*0.05)
#     println( isapprox(flatwise_stress1GX_centrifugal28[ipt],flatwise_stress1_centrifugal28[ipt];atol))
# end

PyPlot.figure("Cent28")
# PyPlot.plot(experimentalCent28[:,2],experimentalCent28[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1GX_centrifugal28[5:end-5]./1e6,spanposmid[5:end-5].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2GX_centrifugal28./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
PyPlot.legend(loc = (.41,0.81))
# PyPlot.legend(loc=4)
# PyPlot.savefig("$(path)/../figs/34m_fig4_3_Centrifugal_28RPM_flapwise_BladeGXNL.pdf",transparent = true)


##################################################################
########### FIG 4.4 40RPM with Gravity Loads Tared out  ##########
##################################################################

model.Omegaocp = [40.0, 40.0]./60
model.OmegaInit = 40.0/60

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist= OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformTurb;steady,system,assembly)

# println("writing")
# OWENS.OWENSFEA_VTK("SNL34m_test",t,uHist,system,assembly,sections;scaling=10)#,azi=aziHist)
# println("done")

# Load data
experimentalCent = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.4_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1GX_centrifugal = (kappa_y[1,end-1,2:end].* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum
flatwise_stress2GX_centrifugal = (kappa_y[2,end-1,1:end-1].* thickness .+ 0*eps_x[2,end-1,1:end-1]) .* Ealuminum
edgewise_stress1GX_centrifugal = (kappa_z[1,end-1,2:end] .* thickness .+ 0*eps_x[1,end-1,2:end]) .* Ealuminum



# Zero out gravity loads
flatwise_stress1GX_centrifugal .-= flatwise_stress1GX
flatwise_stress2GX_centrifugal .-= flatwise_stress2GX
edgewise_stress1GX_centrifugal .-= edgewise_stress1GX

# # Test 40 RPM
# for ipt = 1:length(flatwise_stress1GX_centrifugal)
#     atol = abs(flatwise_stress1GX_centrifugal[ipt]*0.05)
#     println( isapprox(flatwise_stress1GX_centrifugal[ipt],flatwise_stress1_centrifugal40[ipt];atol))
# end

PyPlot.figure("Cent40")
# PyPlot.plot(experimentalCent[:,2],experimentalCent[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(flatwise_stress1GX_centrifugal[5:end-5]./1e6,spanposmid[5:end-5].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2GX_centrifugal./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Flapwise Stress (MPa)")
# PyPlot.legend(loc = (0.06,0.95),ncol=2)
PyPlot.legend(loc = (.41,0.81))
# PyPlot.savefig("$(path)/../figs/34m_fig4_4_Centrifugal_40RPM_flapwise_BladeGXNL.pdf",transparent = true)

PyPlot.figure("Cent40lag")
# PyPlot.plot(experimental_grav[:,2],experimental_grav[:,1],"ko",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(edgewise_stress1GX_centrifugal[1:end-4]./1e6,spanposmid[1:end-4].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
# PyPlot.plot(spanposmid),-flatwise_stress2GX./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.ylabel("Blade Span (m)")
PyPlot.xlabel("Edgewise Stress (MPa)")
# PyPlot.legend(loc = (0.06,1.0),ncol=2)
PyPlot.legend()

