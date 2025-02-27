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
AModel = "DMS"
windINPfilename = nothing
ifw_libfile = nothing
Blade_Height = 41.9
Blade_Radius = 17.1
numTS = 200
delta_t = 0.05
NuMad_geom_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mGeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mMaterials.csv"
NuMad_geom_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mGeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mMaterials.csv"
adi_lib = nothing
adi_rootname = "$(path)/SNL34m"

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
new_RPM = OWENS.safeakima(SNL34m_5_3_RPM[:,1],SNL34m_5_3_RPM[:,2],new_t)

new_Torque = OWENS.safeakima(SNL34m_5_3_Torque[:,1],SNL34m_5_3_Torque[:,2],new_t)

Vinf_spec = OWENS.safeakima(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],new_t)

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



controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = OWENS.safeakima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
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
    shapeY = zero(shapeX),
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
    AModel, #AD, DMS, AC
    DSModel="BV",
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
thickness = OWENS.safeakima(numadIn_bld.span,thickness_precomp_flap,spanposmid)
thickness_lag = OWENS.safeakima(numadIn_bld.span,thickness_precomp_lag,spanposmid)
# thickness = thicknessGX[1:end-1]


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
    numTS,
    delta_t,#dt,
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
eps_x_grav,eps_z_grav,eps_y_grav,kappa_x_grav,kappa_y_grav,kappa_z_grav,t,FReactionHist_grav = OWENS.run34m(model,feamodel,mymesh,myel,aeroForces,deformAero;steady=true)

#####
###****** SAND-88-1144 Specifies Bending Strains and Axial Strains Separate ****
#####

Ealuminum = plyprops_bld.plies[end].e1
flatwise_stress1grav = (kappa_y_grav[1,end-1,2:end].* thickness .+ 0*eps_x_grav[1,end-1,2:end]) .* Ealuminum
flatwise_stress2grav = (kappa_y_grav[2,end-1,1:end-1].* thickness .+ 0*eps_x_grav[2,end-1,1:end-1]) .* Ealuminum
lag_stress1grav = (kappa_z_grav[1,end-1,2:end].* thickness_lag .+ 0*eps_x_grav[1,end-1,2:end]) .* Ealuminum
lag_stress2grav = (kappa_z_grav[2,end-1,1:end-1].* thickness_lag .+ 0*eps_x_grav[2,end-1,1:end-1]) .* Ealuminum

model.Omegaocp = Omegaocp
model.OmegaInit = Omegaocp[1]
model.Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
feamodel.nlOn = false
feamodel.analysisType = "ROM"
model.analysisType = "ROM"

eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,omegaHist,genTorque,torqueDriveShaft,aziHist,uHist = OWENS.run34m(model,feamodel,mymesh,myel,
aeroForces,deformAero;steady=false,system,assembly,VTKFilename="$path/vtk/NormalOperation")

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
