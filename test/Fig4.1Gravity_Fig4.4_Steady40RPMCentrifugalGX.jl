using Statistics:mean
using Test
import HDF5
import PyPlot
import DelimitedFiles
import QuadGK
import FLOWMath
import ModelGen
import GyricFEA
import Composites
import OWENS
import VAWTAero
import RollingFunctions
import GXBeam

steady = false

path = splitdir(@__FILE__)[1]

t_offset = 5.0 #sec

println("Running Gravity Loading")

#Put in one place so its not repeated for all of the analyses
include("$(path)/34mSetup.jl")
# include("$(path)/speedupdebugger.jl")

Vinf = 1e-2#mean(SNL34m_5_3_Vinf[:,2])
TSR = omega*radius./Vinf
ntheta = 200#176

VAWTAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
    eta = 0.25,
    rho,
    mu = 1.7894e-5,
    ntheta,
    Nslices,
    ifw = false,
    RPI = true,
    DSModel = "BV",
    AModel = "DMS",
    tau = [1e-5,1e-5],
    afname = airfoils)

dt = 1/(RPM/60*ntheta)

aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.advanceTurb)

model = OWENS.Inputs(;analysisType = "TNB",
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

feamodel = GyricFEA.FEAModel(;analysisType = "TNB",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn=true,
numNodes = mymesh.numNodes,
numModes = 200,
RayleighAlpha = 0.1,
RayleighBeta = 0.1,
gravityOn = true,
# nodalTerms = GyricFEA.readNodalTerms(;data=[94 "F6" 2 2 1e6;145 "F6" 2 2 1e6]),
iterationType = "DI")

deformTurb(azi_j;newOmega=0,newVinf=0,bld_x=0,bld_z=0,bld_twist=0) = 0

eps_xG,eps_zG,eps_yG,kappa_xG,kappa_yG,kappa_zG,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist = runowens(model,feamodel,mymesh,myel,
aeroForcesDMS,deformTurb;steady,system,assembly)

##################################################################
########### FIG 4.1 Gravity Only #############
##################################################################

# Load data
experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_bothbladesreversed.csv",',',skipstart = 0)
predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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
torqueDriveShaft,aziHist,uHist  = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb;steady,system,assembly)


# Load data
experimentalCent28 = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.3_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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
saveName = "$path/vtk/two_blade"
ModelGen.gyricFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)

##################################################################
########### FIG 4.4 40RPM with Gravity Loads Tared out  ##########
##################################################################

model.Omegaocp = [40.0, 40.0]./60
model.OmegaInit = 40.0/60

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist= runowens(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb;steady,system,assembly)

# println("writing")
# ModelGen.gyricFEA_VTK("SNL34m_test",t,uHist,system,assembly,sections;scaling=10)#,azi=aziHist)
# println("done")

# Load data
experimentalCent = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.4_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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
# include("$(path)/../../../OWENSKevin.jl/src/OWENS.jl")
#Put in one place so its not repeated for all of the analyses
include("$(path)/34mSetup.jl")
# include("$(path)/speedupdebugger.jl")

Vinf = 1e-2#mean(SNL34m_5_3_Vinf[:,2])
TSR = omega*radius./Vinf
ntheta = 200#176

VAWTAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
    eta = 0.25,
    rho,
    mu = 1.7894e-5,
    ntheta,
    Nslices,
    ifw = false,
    RPI = true,
    DSModel = "BV",
    AModel = "DMS",
    tau = [1e-5,1e-5],
    afname = airfoils)

dt = 1/(RPM/60*ntheta)

aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.advanceTurb)
deformTurb(azi_j;newOmega=0,newVinf=0,bld_x=0,bld_z=0,bld_twist=0) = 0

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

feamodel = GyricFEA.FEAModel(;analysisType = "GX",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn=false,
numNodes = mymesh.numNodes,
numModes = 200,
RayleighAlpha = 0.1,
RayleighBeta = 0.1,
gravityOn = true,
# nodalTerms = GyricFEA.readNodalTerms(;data=[94 "F6" 2 2 1e6;145 "F6" 2 2 1e6]),
iterationType = "DI")


eps_xG,eps_zG,eps_yG,kappa_xG,kappa_yG,kappa_zG,t,FReactionHist,OmegaHist,genTorque,
torqueDriveShaft,aziHist,uHist = runowens(model,feamodel,mymesh,myel,
aeroForcesDMS,deformTurb;steady,system,assembly)

##################################################################
########### FIG 4.1 Gravity Only #############
##################################################################

# Load data
experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_bothbladesreversed.csv",',',skipstart = 0)
predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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
PyPlot.plot(edgewise_stress1GX[1:end]./1e6,spanposmid[1:end].+toweroffset,".-",color=plot_cycle[2],label = "OWENS GX Beam")
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
torqueDriveShaft,aziHist,uHist  = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb;steady,system,assembly)


# Load data
experimentalCent28 = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.3_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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
torqueDriveShaft,aziHist,uHist= runowens(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb;steady,system,assembly)

# println("writing")
# ModelGen.gyricFEA_VTK("SNL34m_test",t,uHist,system,assembly,sections;scaling=10)#,azi=aziHist)
# println("done")

# Load data
experimentalCent = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.4_bothbladesreversed.csv",',',skipstart = 0)
# predicted = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1_predicted.csv",',',skipstart = 0)


Ealuminum = plyprops.plies[end].e1
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

