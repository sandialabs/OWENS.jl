
using Statistics:mean
using Test
import HDF5
import DelimitedFiles
import QuadGK
import FLOWMath
import GyricFEA
import OWENS
import VAWTAero
import RollingFunctions
# import FFTW
import GXBeam

path = splitdir(@__FILE__)[1]

import PyPlot
PyPlot.pygui(true)
# close("all")
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

##############################################
# Setup
#############################################

SNL34m_5_11_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.11EmergStopVinf.csv",',',skipstart = 0)
SNL34m_5_11_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.11EmergStopRPM.csv",',',skipstart = 0)
SNL34m_5_11_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.11EmergStopTorque.csv",',',skipstart = 0)

# Plot
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("axes.spines", right=true, top=false)
PyPlot.rc("figure.subplot", left=.15, bottom=.17, top=0.85, right=.9)
PyPlot.figure()
PyPlot.plot(SNL34m_5_11_RPM[:,1],SNL34m_5_11_RPM[:,2],color=plot_cycle[1],".-",label="RPM")
# PyPlot.xlim([0,0.35])
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("RPM")
PyPlot.ylim([-2,50])
PyPlot.legend(loc=(.08,.25))
ax = PyPlot.gca()
ax2 = ax.twinx() # Create another axis on top of the current axis
# new_position = [0.06;0.06;0.77;0.91] # Position Method 2
# ax2.set_position(new_position) # Position Method 2
PyPlot.xlabel("Time (s)")
PyPlot.plot(SNL34m_5_11_Vinf[:,1],SNL34m_5_11_Vinf[:,2],"k.-",label="Windspeed (m/s)")
PyPlot.xlabel("X (m)")
PyPlot.ylabel("Windspeed (m/s)")
PyPlot.ylim([0,13.5])
PyPlot.legend(loc=(.08,.15))
# PyPlot.savefig("$(path)/../figs/Emerg_stop_Vinf_RPM.pdf",transparent = true)


# smooth out the input data
SNL34m_5_11_RPM = SNL34m_5_11_RPM[1:6:end,:]
SNL34m_5_11_Vinf = SNL34m_5_11_Vinf[1:2:end,:]
SNL34m_5_11_RPM[SNL34m_5_11_RPM[:,2].<0.0,2] .*= 0

new_t = LinRange(SNL34m_5_11_RPM[1,1],SNL34m_5_11_RPM[end,1],100).-SNL34m_5_11_RPM[1,1].+1e-6
new_RPM = FLOWMath.akima(SNL34m_5_11_RPM[:,1].-SNL34m_5_11_RPM[1,1],SNL34m_5_11_RPM[:,2],new_t)

new_Torque = FLOWMath.akima(SNL34m_5_11_Torque[:,1],SNL34m_5_11_Torque[:,2],new_t)

t_Vinf = LinRange(SNL34m_5_11_Vinf[1,1],SNL34m_5_11_Vinf[end,1],60).-SNL34m_5_11_RPM[1,1].+1e-6
Vinf_spec = FLOWMath.akima(SNL34m_5_11_Vinf[:,1].-SNL34m_5_11_RPM[1,1],SNL34m_5_11_Vinf[:,2],t_Vinf)

offsetTime = 20.0
t_Vinf = [0;t_Vinf.+offsetTime;t_Vinf[end]*10]
Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
PyPlot.figure()
PyPlot.plot(t_Vinf,Vinf_spec,".-",label="Orig")
PyPlot.xlabel("t")
PyPlot.ylabel("V")

PyPlot.figure()
PyPlot.plot(new_t,new_RPM,".-",label="Orig")
PyPlot.plot(new_t[1:end-1],diff(new_RPM),".-",label="Orig")
PyPlot.xlabel("t")
PyPlot.ylabel("RPM")
#
# PyPlot.figure()
# PyPlot.plot(new_t,new_Torque,".-",label="Orig")
# PyPlot.xlabel("t")
# PyPlot.ylabel("Torque")
#
# PyPlot.figure()
# PyPlot.plot(new_t,new_RPM,".-",label="Orig")
# PyPlot.xlabel("t")
# PyPlot.ylabel("RPM")
#
# PyPlot.figure()
# PyPlot.plot(new_t,Vinf_spec,".-",label="Orig")
# PyPlot.xlabel("t")
# PyPlot.ylabel("Vinf")
# PyPlot.plot(t_Vinf,Vinf_spec,label="New")
# PyPlot.legend()

#Put in one place so its not repeated for all of the analyses
include("$(path)/34mSetup.jl")

Vinf = mean(SNL34m_5_11_Vinf[:,2])
TSR = omega*radius./Vinf
windpower = 0.5*rho*Vinf.^3*RefArea

ntheta = 30

VAWTAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
    eta = 0.25,
    rho,
    mu = 1.7894e-5,
    ntheta,
    Nslices,
    ifw = false,
    turbsim_filename = "$path/data/40mx40mVinf10_41ms10percturb.bts",
    RPI = true,
    DSModel = "BV",
    AModel = "DMS",
    tau = [1e-5,1e-5],
    afname = airfoils)

# # UnSteady
# CP,
# Rp,
# Tp,
# Zp,
# alpha,
# cl_af,
# cd_af,
# Vloc,
# Re,
# thetavec,
# nstep,
# Fx_base,
# Fy_base,
# Fz_base,
# Mx_base,
# My_base,
# Mz_base,
# power,
# power2 = VAWTAero.AdvanceTurbineInterpolate(maximum(SNL34m_5_11_Vinf[:,1]))#VAWTAero.steadyTurb(omega,Vinf)

aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
# tnew = 50.0
#
# CP,Rp,Tp,Zp,alpha,cl_af,cd_af,Vloc,Re,thetavec,nstep,Fx_base,Fy_base,Fz_base,
# Mx_base,My_base,Mz_base = VAWTAero.AdvanceTurbineInterpolate(tnew)
#
# filterwindow = 20
# PyPlot.rc("figure.subplot", left=.22, bottom=.17, top=0.9, right=.9)
# PyPlot.figure()
# PyPlot.plot(LinRange(0,tnew,nstep),Mz_base/1000,color=plot_cycle[1],label="DMS Aero Only")
# rollingave = RollingFunctions.runmean(Mz_base/1000,filterwindow)
# PyPlot.plot(LinRange(0,tnew,nstep),rollingave,"r.-",label="Rolling Average")
# PyPlot.xlabel("Time (s)")
# PyPlot.ylabel("Base Torque (kN-m)")
# PyPlot.legend()
# PyPlot.savefig("$(path)/../figs/34m_baseTorqueSteadyFig5.3_test.pdf",transparent = true)

Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]]./60

model = OWENS.Inputs(;analysisType = "ROM",
    outFilename = "none",
    tocp = [0.0; new_t.+offsetTime; 1e6],#SNL34m_5_11_RPM[:,1],#[0.0,10.0,100000.1],
    Omegaocp,#SNL34m_5_11_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
    tocp_Vinf = t_Vinf,
    Vinfocp = Vinf_spec,
    numTS = 2000,
    delta_t = 0.05,
    aeroLoadsOn = 2,
    turbineStartup = 0,
    generatorOn = false,
    useGeneratorFunction = true,
    OmegaInit = new_RPM[1]/60)

feamodel = GyricFEA.FEAModel(;analysisType = "ROM",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn=false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.01,
    RayleighBeta = 0.01,
    iterationType = "DI")

# Get Gravity Loads
model.Omegaocp = model.Omegaocp.*0.0
model.OmegaInit = model.OmegaInit.*0.0
model.Vinfocp = model.Vinfocp.*0.0
feamodel.nlOn = true

# Returns data filled with e.g. eps[Nbld,N_ts,Nel_bld]
eps_x_grav,eps_z_grav,eps_y_grav,kappa_x_grav,kappa_y_grav,kappa_z_grav,t,FReactionHist_grav = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,VAWTAero.deformTurb;steady=true)

Ealuminum = plyprops.plies[end].e1
flatwise_stress1grav = (kappa_y_grav[1,end,1:end-1].* thickness .+ 0*eps_x_grav[1,end,1:end-1]) .* Ealuminum
flatwise_stress2grav = (kappa_y_grav[2,end,1:end-1].* thickness .+ 0*eps_x_grav[2,end,1:end-1]) .* Ealuminum
lag_stress1grav = (kappa_z_grav[2,end,1:end-1].* thickness_lag .+ 0*eps_x_grav[2,end,1:end-1]) .* Ealuminum
# experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_b1.csv",',',skipstart = 0)
# experimental_grav = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.1_bothbladesreversed.csv",',',skipstart = 0)
# PyPlot.figure()
# PyPlot.plot(experimental_grav[:,1],experimental_grav[:,2],"ko",label = "Experimental")
# # PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
# PyPlot.plot((spanposmid),flatwise_stress1grav./1e6,".-",color=plot_cycle[1],label = "OWENS Blade 1")
# PyPlot.plot((spanposmid),flatwise_stress2grav./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
# PyPlot.xlabel("Blade Position (m)")
# PyPlot.ylabel("Flapwise Stress (MPa)")
# PyPlot.legend()
# # PyPlot.savefig("$(path)/../figs/34m_fig4_1_GravityOnly_flapwise_Blade.pdf",transparent = true)

println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld)


model.Omegaocp = Omegaocp
model.OmegaInit = Omegaocp[1]
model.Vinfocp = Vinf_spec
feamodel.nlOn = false
eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,VAWTAero.deformTurb;steady=false,system,assembly,VTKFilename="$path/vtk/EmergencyStop")

# Get stress and "zero" out the loads from the initial 0-RPM
flatwise_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
flatwise_stress2 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
lag_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
for its = 1:length(eps_x[1,:,1])
    flatwise_stress1[its,:] = (kappa_y[1,its,1:end-1].* thickness .+ 0*eps_x[1,its,1:end-1]) .* Ealuminum .- flatwise_stress1grav
    flatwise_stress2[its,:] = (kappa_y[2,its,1:end-1].* thickness .+ 0*eps_x[2,its,1:end-1]) .* Ealuminum .- flatwise_stress2grav

    lag_stress1[its,:] = (kappa_z[1,its,1:end-1].* thickness_lag .+ 0*eps_x[1,its,1:end-1]) .* Ealuminum .- lag_stress1grav
end

# PyPlot.figure()
# for its = length(eps_x[1,:,1])-100:length(eps_x[1,:,1])
#     PyPlot.cla()
#     PyPlot.plot((spanposmid),flatwise_stress1[its,:]./1e6,".-",color=plot_cycle[1],label = "OWENS Blade 1")
#     # PyPlot.plot((spanposmid),flatwise_stress2grav./1e6,".-",color=plot_cycle[2],label = "OWENS Blade 2")
#     PyPlot.xlabel("Blade Position (m)")
#     PyPlot.ylabel("Edgewise Stress (MPa)")
#     PyPlot.title("Time: $(t[its])")
#     PyPlot.ylim([-40,40])
#     # PyPlot.legend()
#     sleep(0.0001)
#     # PyPlot.savefig("$(path)/../figs/34m_fig4_1_GravityOnly_flapwise_Blade.pdf",transparent = true)
# end

# Load in experimentat data
SNL34m_5_11_FlatwiseStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.12_EmergStopFlap.csv",',',skipstart = 1)

# Plots
PyPlot.figure()
PyPlot.plot(SNL34m_5_11_FlatwiseStress[:,1].-SNL34m_5_11_FlatwiseStress[1,1],SNL34m_5_11_FlatwiseStress[:,2],"k-",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(t[1:end].-offsetTime,flatwise_stress1[:,6]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
# PyPlot.plot(t[1:end],-flatwise_stress2[:,3]./1e6,"-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Flapwise Stress (MPa)")
PyPlot.xlim([0,80])
PyPlot.legend(loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig5_11_EmergStop_flapwise_Blade.pdf",transparent = true)

SNL34m_5_11_LeadLagStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.12_EmergStopLeadLag.csv",',',skipstart = 1)

PyPlot.figure()
PyPlot.plot(SNL34m_5_11_LeadLagStress[:,1].-SNL34m_5_11_LeadLagStress[1,1],SNL34m_5_11_LeadLagStress[:,2],"k-",label = "Experimental")
# PyPlot.plot(predicted[:,1],predicted[:,2],"k-",label = "Predicted")
PyPlot.plot(t[1:end].-offsetTime,lag_stress1[:,6]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
# PyPlot.plot(t[1:end],-flatwise_stress2[:,3]./1e6,"-",color=plot_cycle[2],label = "OWENS Blade 2")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Lead-Lag Stress (MPa)")
PyPlot.xlim([0,80])
PyPlot.legend(loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig5_11_EmergStop_LeadLag_Blade.pdf",transparent = true)


##########################################
#### Torque Plot
##########################################

SNL34m_5_11_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.11EmergStopTorque.csv",',',skipstart = 0)

filterwindow = 100
PyPlot.ion()
PyPlot.figure()
# PyPlot.plot(thetavec[1,:]/omega,Mz_base/1000,color=plot_cycle[1],label="DMS Aero Only")
# PyPlot.plot(t,FReactionHist[:,1]/1000,color=plot_cycle[1],label="One-Way Aer1")
# PyPlot.plot(t,-FReactionHist[:,2]/1000,color=plot_cycle[2],label="One-Way Aer2")
# PyPlot.plot(t,-FReactionHist[:,3]/1000,color=plot_cycle[3],label="One-Way Aer3")
# PyPlot.plot(t,-FReactionHist[:,4]/1000,color=plot_cycle[4],label="One-Way Aer4")
# PyPlot.plot(t,FReactionHist[:,5]/1000,color=plot_cycle[5],label="One-Way Aero5")
PyPlot.plot(t.-offsetTime,-FReactionHist[:,6]/1000 ,color=plot_cycle[1],label="OWENS")
# PyPlot.plot(t,OmegaHist.*60,color=plot_cycle[1],label="RPM")
# rollingave = RollingFunctions.runmean(FReactionHist[:,6]/1000,filterwindow)
# PyPlot.plot(t,rollingave,"r.-",label="Rolling Average")
# usedLogic = SNL34m_5_11_Torque[:,1].<100
PyPlot.plot(SNL34m_5_11_Torque[:,1].-SNL34m_5_11_Torque[1,1],SNL34m_5_11_Torque[:,2],"k-",label="Experimental")
PyPlot.xlabel("Time (s)")
PyPlot.xlim([0,80])
PyPlot.ylabel("Torque (kN-m)")
PyPlot.legend()
# PyPlot.savefig("$(path)/../figs/34m_fig5_11EmergStop.pdf",transparent = true)

# # Extract and Plot Frequency - Amplitude
# signal = -flatwise_stress1[:,4]./1e6
# signal .-= mean(signal)
# Y = FFTW.fft(signal)

# L = length(signal)-1
# Fs = 1/(t[2]-t[1])
# P2 = abs.(Y./L)
# P1 = P2[1:Int(L/2)+1]
# P1[2:end] = 2*P1[2:end]

# f = Fs.*(0:Int(L/2))./L
# PyPlot.figure()
# PyPlot.plot(f,P1)
# PyPlot.title("Single-Sided Amplitude Spectrum of X(t)")
# PyPlot.xlabel("f (Hz)")
# PyPlot.ylabel("|P1(f)|")
# PyPlot.xlim([0.0,10.0])
