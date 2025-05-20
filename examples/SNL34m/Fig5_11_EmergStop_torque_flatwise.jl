
using Statistics:mean
using Test
import HDF5
import DelimitedFiles
import QuadGK
import FLOWMath
import OWENSFEA
import OWENS
import OWENSAero
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

nothing

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
new_RPM = safeakima(SNL34m_5_11_RPM[:,1].-SNL34m_5_11_RPM[1,1],SNL34m_5_11_RPM[:,2],new_t)

new_Torque = safeakima(SNL34m_5_11_Torque[:,1],SNL34m_5_11_Torque[:,2],new_t)

t_Vinf = LinRange(SNL34m_5_11_Vinf[1,1],SNL34m_5_11_Vinf[end,1],60).-SNL34m_5_11_RPM[1,1].+1e-6
Vinf_spec = safeakima(SNL34m_5_11_Vinf[:,1].-SNL34m_5_11_RPM[1,1],SNL34m_5_11_Vinf[:,2],t_Vinf)

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

Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]]./60

model = OWENS.Inputs(;verbosity,analysisType = "ROM",
    dataOutputFilename = "none",
    tocp = [0.0; new_t.+offsetTime; 1e6],#SNL34m_5_11_RPM[:,1],#[0.0,10.0,100000.1],
    Omegaocp,#SNL34m_5_11_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
    tocp_Vinf = t_Vinf,
    Vinfocp = Vinf_spec,
    numTS,
    delta_t,
    aeroLoadsOn = 2,
    turbineStartup = 0,
    generatorOn = false,
    useGeneratorFunction = true,
    OmegaInit = new_RPM[1]/60)

feamodel = OWENSFEA.FEAModel(;analysisType = "ROM",
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
eps_x_grav,eps_z_grav,eps_y_grav,kappa_x_grav,kappa_y_grav,kappa_z_grav,t,FReactionHist_grav = OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformAero;steady=true)

Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1grav = (kappa_y_grav[1,end,1:end-1].* thickness .+ 0*eps_x_grav[1,end,1:end-1]) .* Ealuminum
flatwise_stress2grav = (kappa_y_grav[2,end,1:end-1].* thickness .+ 0*eps_x_grav[2,end,1:end-1]) .* Ealuminum
lag_stress1grav = (kappa_z_grav[2,end,1:end-1].* thickness_lag .+ 0*eps_x_grav[2,end,1:end-1]) .* Ealuminum

model.Omegaocp = Omegaocp
model.OmegaInit = Omegaocp[1]
model.Vinfocp = Vinf_spec
feamodel.nlOn = false
model.analysisType = "ROM"
feamodel.analysisType = "ROM"
eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist = OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformAero;steady=false,system,assembly,VTKFilename="$path/vtk/EmergencyStop")

# Get stress and "zero" out the loads from the initial 0-RPM
flatwise_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
flatwise_stress2 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
lag_stress1 = zeros(length(eps_x[1,:,1]),length(eps_x[1,1,1:end-1]))
for its = 1:length(eps_x[1,:,1])
    flatwise_stress1[its,:] = (kappa_y[1,its,1:end-1].* thickness .+ 0*eps_x[1,its,1:end-1]) .* Ealuminum .- flatwise_stress1grav
    flatwise_stress2[its,:] = (kappa_y[2,its,1:end-1].* thickness .+ 0*eps_x[2,its,1:end-1]) .* Ealuminum .- flatwise_stress2grav

    lag_stress1[its,:] = (kappa_z[1,its,1:end-1].* thickness_lag .+ 0*eps_x[1,its,1:end-1]) .* Ealuminum .- lag_stress1grav
end



# UnitFilename = "$path/data/EmergencyStop34m_UNIT.h5"

# HDF5.h5open(UnitFilename, "w") do file
#     HDF5.write(file,"t",t)
#     HDF5.write(file,"flatwise_stress1",flatwise_stress1)
#     HDF5.write(file,"lag_stress1",lag_stress1)
# end

# t_UNIT = HDF5.h5read(UnitFilename,"t")
# flatwise_stress1_UNIT = HDF5.h5read(UnitFilename,"flatwise_stress1")
# lag_stress1_UNIT = HDF5.h5read(UnitFilename,"lag_stress1")

# for (i_t,t) in enumerate(t)
#     @test isapprox(t[i_t],t_UNIT[i_t];atol = abs(t_UNIT[i_t])*0.001)
#     for iel = 1:length(flatwise_stress1[i_t,:])
#         @test isapprox(flatwise_stress1[i_t,iel],flatwise_stress1_UNIT[i_t,iel];atol = abs(flatwise_stress1_UNIT[i_t,iel])*0.001)
#         @test isapprox(lag_stress1[i_t,iel],lag_stress1_UNIT[i_t,iel];atol = abs(lag_stress1_UNIT[i_t,iel])*0.001)
#     end
# end

# Load in experimentat data
SNL34m_5_11_FlatwiseStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.12_EmergStopFlap.csv",',',skipstart = 1)

# Stress Plots
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
