using PyPlot
PyPlot.pygui(true)
using Statistics:mean
using Statistics
# close("all")
using Test
import PyPlot
import DelimitedFiles
import FLOWMath
import OWENSFEA
import OWENS
import OWENSAero
# import FFTW

path = splitdir(@__FILE__)[1]

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

SNL34m_5_3_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Vinf.csv",',',skipstart = 0)
SNL34m_5_3_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_RPM.csv",',',skipstart = 0)
SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque.csv",',',skipstart = 0)


# Plot
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.85, right=.9)
PyPlot.figure()
PyPlot.plot(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],"k.-",label="Windspeed (m/s)")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Windspeed (m/s)")
PyPlot.ylim([0,12.0])
PyPlot.xlim([0,100.0])
# PyPlot.legend(loc=(.4,.3))
# PyPlot.savefig("$(path)/../figs/NormalOperation_Vinf.pdf",transparent = true)

new_t = LinRange(SNL34m_5_3_RPM[1,1],SNL34m_5_3_RPM[end,1],100)
new_RPM = FLOWMath.akima(SNL34m_5_3_RPM[:,1],SNL34m_5_3_RPM[:,2],new_t)

new_Torque = FLOWMath.akima(SNL34m_5_3_Torque[:,1],SNL34m_5_3_Torque[:,2],new_t)

Vinf_spec = FLOWMath.akima(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],new_t)

t_Vinf = [0;new_t;1e6]
Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
# PyPlot.figure()
# PyPlot.plot(new_RPM,new_Torque,".",label="Orig")
# PyPlot.xlabel("RPM")
# PyPlot.ylabel("Torque")
#
# PyPlot.figure()
# PyPlot.plot(Vinf_spec,new_Torque,".",label="Orig")
# PyPlot.xlabel("Vinf")
# PyPlot.ylabel("Torque")
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

Vinf = mean(SNL34m_5_3_Vinf[:,2])
TSR = omega*radius./Vinf
windpower = 0.5*rho*Vinf.^3*RefArea
ntheta = 30#176

OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
    eta = 0.5,
    rho,
    mu = 1.7894e-5,
    ntheta,
    Nslices,
    ifw = false,
    wind_filename = "$path/data/40mx40mVinf10_41ms10percturb.bts",
    RPI = true,
    DSModel = "BV",
    AModel = "DMS",
    tau = [1e-5,1e-5],
    afname = airfoils)

dt = 1/(RPM/60*ntheta)

aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)

offsetTime = 20.0
Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]]./60 .*0 .+33.92871/60
tocp_Vinf = [0.0;t_Vinf.+offsetTime; 1e6]
Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]].*1e-6

model = OWENS.Inputs(;analysisType = "ROM",
    outFilename = "none",
    tocp = [0.0;new_t.+offsetTime; 1e6],#SNL34m_5_3_RPM[:,1],#[0.0,10.0,100000.1],
    Omegaocp,#SNL34m_5_3_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
    tocp_Vinf,
    Vinfocp,
    numTS = 2000,
    delta_t = 0.05,#dt,
    aeroLoadsOn = 2,
    turbineStartup = 1,
    generatorOn = true,
    useGeneratorFunction = true,
    driveTrainOn = true,
    JgearBox = 250.0,#(2.15e3+25.7)/12*1.35582*100,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    driveShaftProps = OWENS.DriveShaftProps(10000,1.5e2), #8.636e5*1.35582*0.6
    OmegaInit = Omegaocp[1]/60)

println(sqrt(model.driveShaftProps.k/model.JgearBox)*60/2/pi/2)

feamodel = OWENSFEA.FEAModel(;analysisType = "ROM",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn=false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")

# Get Gravity Loads
model.Omegaocp = model.Omegaocp.*0.0
model.OmegaInit = model.OmegaInit.*0.0
model.Vinfocp = model.Vinfocp.*0.0
feamodel.nlOn = true

# Returns data filled with e.g. eps[Nbld,N_ts,Nel_bld]
eps_x_grav,eps_z_grav,eps_y_grav,kappa_x_grav,kappa_y_grav,kappa_z_grav,t,FReactionHist_grav = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,OWENSAero.deformTurb;steady=true)

#####
###****** SAND-88-1144 Specifies Bending Strains and Axial Strains Separate ****
#####

Ealuminum = plyprops.plies[end].e1
flatwise_stress1grav = (kappa_y_grav[1,end-1,2:end].* thickness .+ 0*eps_x_grav[1,end-1,2:end]) .* Ealuminum
flatwise_stress2grav = (kappa_y_grav[2,end-1,1:end-1].* thickness .+ 0*eps_x_grav[2,end-1,1:end-1]) .* Ealuminum
lag_stress1grav = (kappa_z_grav[1,end-1,2:end].* thickness_lag .+ 0*eps_x_grav[1,end-1,2:end]) .* Ealuminum
lag_stress2grav = (kappa_z_grav[2,end-1,1:end-1].* thickness_lag .+ 0*eps_x_grav[2,end-1,1:end-1]) .* Ealuminum

# println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld)#;damp_coef=0.05)


model.Omegaocp = Omegaocp
model.OmegaInit = Omegaocp[1]
model.Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
feamodel.nlOn = false
feamodel.analysisType = "ROM"
model.analysisType = "ROM"

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,omegaHist,genTorque,torqueDriveShaft,aziHist,uHist = runowens(model,feamodel,mymesh,myel,
aeroForcesDMS,OWENSAero.deformTurb;steady=false,system,assembly,VTKFilename="$path/vtk/NormalOperation")

# Get stress and "zero" out the loads from the initial 0-RPM
flatwise_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
flatwise_stress2 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
lag_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
lag_stress2 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
for its = 1:length(eps_x[1,:,1])
    flatwise_stress1[its,:] = (kappa_y[1,its,2:end].* thickness .+ 0*eps_x[1,its,2:end]) .* Ealuminum .- flatwise_stress1grav
    flatwise_stress2[its,:] = (kappa_y[2,its,1:end-1].* thickness .+ 0*eps_x[2,its,1:end-1]) .* Ealuminum .- flatwise_stress2grav

    lag_stress1[its,:] = (kappa_z[1,its,2:end].* thickness_lag .+ 0*eps_x[1,its,2:end]) .* Ealuminum .- lag_stress1grav
    lag_stress2[its,:] = (kappa_z[2,its,1:end-1].* thickness_lag .+ 0*eps_x[2,its,1:end-1]) .* Ealuminum .- lag_stress2grav
end

# Load in experimental data
SNL34m_5_4_FlatwiseStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.4AMF.csv",',',skipstart = 1)

# Plots
PyPlot.figure()
PyPlot.plot(t.-offsetTime,flatwise_stress1[:,end-5]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
PyPlot.plot(SNL34m_5_4_FlatwiseStress[:,1].-0.8,SNL34m_5_4_FlatwiseStress[:,2],"k-",label = "Experimental")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Flapwise Stress (MPa)")
PyPlot.xlim([0,SNL34m_5_4_FlatwiseStress[end,1]])
# PyPlot.ylim([-60,0])
PyPlot.legend(loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig5_4_NormalOperation_flapwise_Blade2Way.pdf",transparent = true)

println("here")
exp_std_flap = Statistics.std(SNL34m_5_4_FlatwiseStress[:,2])
println("exp_std_flap $exp_std_flap")
exp_mean_flap = Statistics.mean(SNL34m_5_4_FlatwiseStress[:,2])
println("exp_mean_flap $exp_mean_flap")
sim_std_flap = Statistics.std(flatwise_stress1[:,end-4]./1e6)
println("sim_std_flap $sim_std_flap")
sim_mean_flap = Statistics.mean(flatwise_stress1[:,end-4]./1e6)
println("sim_mean_flap $sim_mean_flap")

SNL34m_5_4_LeadLagStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.4AML.csv",',',skipstart = 1)

PyPlot.figure()
PyPlot.plot(t.-offsetTime,lag_stress1[:,end-5]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
PyPlot.plot(SNL34m_5_4_LeadLagStress[:,1],SNL34m_5_4_LeadLagStress[:,2],"k-",label = "Experimental")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Lead-Lag Stress (MPa)")
PyPlot.xlim([0,SNL34m_5_4_LeadLagStress[end,1]])
PyPlot.legend(loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig5_4_NormalOperation_LeadLag_Blade2Way.pdf",transparent = true)

exp_std_lag = Statistics.std(SNL34m_5_4_LeadLagStress[:,2])
println("exp_std_lag $exp_std_lag")
exp_mean_lag = Statistics.mean(SNL34m_5_4_LeadLagStress[:,2])
println("exp_mean_lag $exp_mean_lag")
sim_std_lag = Statistics.std(lag_stress1[:,end-4]./1e6)
println("sim_std_lag $sim_std_lag")
sim_mean_lag = Statistics.mean(lag_stress1[:,end-4]./1e6)
println("sim_mean_lag $sim_mean_lag")

##########################################
#### Torque Plot
##########################################

SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque2.csv",',',skipstart = 0)

filterwindow = 100
PyPlot.ion()
PyPlot.figure()
PyPlot.plot(t.-offsetTime,torqueDriveShaft/1000 ,color=plot_cycle[1],label="Simulated Drive Shaft")
usedLogic = SNL34m_5_3_Torque[:,1].<100
PyPlot.plot(SNL34m_5_3_Torque[usedLogic,1],SNL34m_5_3_Torque[usedLogic,2],"k-",label="Experimental")
PyPlot.xlabel("Time (s)")
PyPlot.xlim([0,100])
PyPlot.ylabel("Torque (kN-m)")
PyPlot.legend()#loc = (0.06,1.0),ncol=2)
# PyPlot.savefig("$(path)/../figs/34m_fig5_32Way.pdf",transparent = true)

# #Extract and Plot Frequency - Amplitude

# lensig = length(flatwise_stress1[:,1])
# ipt = length(flatwise_stress1[1,:])-5
# istart = 1#500
# iend = lensig
# signal = -flatwise_stress1[istart:iend,ipt]./1e6
# L = length(signal)
# if L%2 != 0
#     signal = signal[1:end-1]
#     L = length(signal)
# end
# # signal .-= mean(signal)
# Y = FFTW.fft(signal)


# Fs = 1/(t[2]-t[1])
# P2 = abs.(Y./L)
# P1 = P2[1:Int(L/2)+1]
# P1[2:end-1] = 2*P1[2:end-1]

# f = Fs.*(0:Int(L/2))./L
# PyPlot.figure()
# PyPlot.plot(f,P1)
# PyPlot.title("Single-Sided Amplitude Spectrum of X(t)")
# PyPlot.xlabel("f (Hz)")
# PyPlot.ylabel("|P1(f)|")
# PyPlot.xlim([0.0,10.0])
# PyPlot.ylim([0.0,7.0])
# # PyPlot.savefig("$(path)/../figs/34m_fig5_6_upperRootFlatwiseStressSpectrum2Way.pdf",transparent = true)

################################################################
################ SAVE VTK TIME DOMAIN OUTPUT ###################
################################################################

# println("Saving VTK time domain files")
# OWENS.OWENSFEA_VTK("$path/vtk/SNL34m_timedomain",t,uHist,system,assembly,sections;scaling=10)#,azi=aziHist)

# Open Paraview, open animation pane, adjust as desired, export animation (which exports frames)
# ffmpeg -i Ux.%04d.png -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 -y -an -pix_fmt yuv420p video34m34RPM_Ux.mp4
