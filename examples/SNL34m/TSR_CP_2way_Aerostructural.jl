
using Statistics:mean
import HDF5
import DelimitedFiles
import FLOWMath
import OWENSFEA
import GXBeam
import OWENS
import OWENSAero
import QuadGK


path = splitdir(@__FILE__)[1]

using PyPlot
PyPlot.pygui(true)
close("all")
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

# ##############################################
# # Setup
# #############################################
# #Put in one place so its not repeated for all of the analyses
# include("$(path)/34mSetup.jl")
# numTS = 125
# RPMsetpoint = 34.0
# omega = RPMsetpoint/60*2*pi
# nTSR = 5
# # TSRvec = LinRange(2.43,15.22,nTSR) #[5.843]
# offsetTime = 1e-6
# t_Vinf = LinRange(0,1e6,10)
# # Vinf_array = omega*radius./TSRvec
# Vinf_array = [5.0]#collect(LinRange(5.0,25.0,nTSR))
# TSRvec = omega*radius./Vinf_array
# mytorque = zeros(nTSR,numTS-1)

# # for (iTSR,TSRspec) in enumerate(TSRvec)
# iTSR = 1
# TSRspec = TSRvec[iTSR]
# Vinf_spec = ones(10).*Vinf_array[iTSR]

# Vinf = Vinf_spec[1]#mean(SNL34m_5_3_Vinf[:,2])
# TSR = omega*radius./Vinf
# # windpower = 0.5*rho*Vinf.^3*RefArea
# ntheta = 30#176
# global rho = 1.225 #override since it was corrected for sea level
# OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
#     eta = 0.5,
#     rho,
#     mu = 1.7894e-5,
#     ntheta,
#     Nslices,
#     ifw = false,
#     turbsim_filename = "$path/data/40mx40mVinf10_41ms10percturb.bts",
#     RPI = true,
#     DSModel = "BV",
#     AModel = "DMS",
#     tau = [1e-5,1e-5],
#     afname = airfoils)

# dt = 1/(RPM/60*ntheta)

# aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=false)

#######################
# New Setup
#######################
verbosity = 1
RPMsetpoint = 34.0
omega = RPMsetpoint/60*2*pi
turbineType = "Darrieus"
eta = 0.5
Nbld = 2
towerHeight = 0.5
rho = 1.225
Vinf = 10.1
controlStrategy = "constantRPM"
RPM = 34.0
Nslices = 35
ntheta = 30
structuralModel = "GX"
ntelem = 10
nbelem = 60
ncelem = 10
nselem = 5
ifw = false
AModel = "DMS"
windINPfilename = "$(path)/data/turbsim/115mx115m_30x30_20.0msETM.bts"
ifw_libfile = nothing
Blade_Height = 41.9
Blade_Radius = 17.1

Vinf_array = [5.0]#collect(LinRange(5.0,25.0,nTSR))
TSRvec = omega*Blade_Radius./Vinf_array

numTS = 50
delta_t = 0.05
NuMad_geom_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mGeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mMaterials.csv"
NuMad_geom_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mGeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mMaterials.csv"
adi_lib = nothing

adi_rootname = "/SNL34m"

##############################################
# Setup
#############################################

SNL34m_5_3_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Vinf.csv",',',skipstart = 0)
SNL34m_5_3_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_RPM.csv",',',skipstart = 0)
SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque.csv",',',skipstart = 0)


new_t = LinRange(SNL34m_5_3_RPM[1,1],SNL34m_5_3_RPM[end,1],100)
new_RPM = FLOWMath.akima(SNL34m_5_3_RPM[:,1],SNL34m_5_3_RPM[:,2],new_t)

Vinf_spec = FLOWMath.akima(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],new_t)

offsetTime = 20.0 # seconds
t_Vinf = [0;new_t;1e6]
tocp_Vinf = [0.0;t_Vinf.+offsetTime; 1e6]
Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
Omegaocp = zero(tocp_Vinf) .+RPMsetpoint/60

controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = FLOWMath.akima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
SNL34_unit_xz = [x_shape;;z_shape]
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])
SNL34Z = SNL34z.*Blade_Height
SNL34X = SNL34x.*Blade_Radius

shapeZ = SNL34Z#collect(LinRange(0,H,Nslices+1))
shapeX = SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

shapeX_spline = FLOWMath.Akima(SNL34Z, SNL34X)
RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, Blade_Height, atol=1e-10)
RefArea = RefArea_half*2

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
    # ncelem=3,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.04,0.96],
    strut_bld_mountpoint = [0.04,0.96],
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = pi/2,
    meshtype = turbineType)

top_idx = 23#Int(myjoint[7,2])
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

model = OWENS.Inputs(;analysisType = structuralModel,
    outFilename = "none",
    tocp = tocp_Vinf,#SNL34m_5_3_RPM[:,1],#[0.0,10.0,100000.1],
    Omegaocp,#SNL34m_5_3_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
    tocp_Vinf,
    Vinfocp,
    numTS,
    delta_t,#dt,
    aeroLoadsOn = 2,
    turbineStartup = 0,
    generatorOn = false,
    useGeneratorFunction = false,
    driveTrainOn = false,
    JgearBox = 250.0,#(2.15e3+25.7)/12*1.35582*100,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    driveShaftProps = OWENS.DriveShaftProps(10000,1.5e2), #8.636e5*1.35582*0.6
    OmegaInit = Omegaocp[1])

feamodel = OWENSFEA.FEAModel(;analysisType = structuralModel,
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.1,
    RayleighBeta = 0.1,
    iterationType = "DI")

# eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,torqueDriveShaft,aziHist,uHist = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,OWENSAero.deformTurb;steady=false,system=system,assembly=assembly)

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(model;
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,system,assembly)

mytorque = zeros(1,numTS-1)
mytorque[1,:] = FReactionHist[:,6]#torqueDriveShaft
# end

# import HDF5
# filename = "$(path)/data/SNL34m34RPM_2way_TSR_Sweep4.h5"
# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"t",t)
#     HDF5.write(file,"mytorque",mytorque)
#     HDF5.write(file,"RPMsetpoint",RPMsetpoint)
#     HDF5.write(file,"Vinf_array",Vinf_array)
#     HDF5.write(file,"RefArea",RefArea)
#     HDF5.write(file,"rho",rho)
#     HDF5.write(file,"radius",radius)
# end

# # Reset and run aero only
# OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSRvec[1],Vinf_array[1];
#     eta = 0.5,
#     rho,
#     mu = 1.7894e-5,
#     ntheta = 30,
#     Nslices,
#     RPI=true,
#     ifw = false,
#     DSModel = "BV",
#     AModel = "DMS",
#     tau = [1e-5,1e-5],
#     afname = airfoils)

# t = 0:0.05:30
Mz_base = zero(t)
Xpbase = zero(t)
Ypbase = zero(t)
Xpbase2 = zero(t)
Ypbase2 = zero(t)
myazi = zero(t)
for (i,myt) in enumerate(t)
    azi = omega*myt + 270/360*2*pi +0.1780235837034216
    myazi[i] = azi
    CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,_,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base[i],power,power2,_,z3Dnorm,delta,Xp,Yp = OWENSAero.AdvanceTurbineInterpolate(myt;azi,alwaysrecalc=false)
    for ibld = 1:length(Xp[:,1,1])
        Xpbase[i] += OWENSAero.trapz(z3Dnorm.*Blade_Height,Xp[ibld,:,end])
        Ypbase[i] += OWENSAero.trapz(z3Dnorm.*Blade_Height,Yp[ibld,:,end])
        Xpbase2[i] = Xpbase[i]*cos(-azi) + Ypbase[i]*sin(-azi)
        Ypbase2[i] = -Xpbase[i]*sin(-azi) + Ypbase[i]*cos(-azi)
    end
end


CPsteady,Rpsteady,Tpsteady,Zpsteady,alphasteady,cl_afsteady,cd_afsteady,Vlocsteady,Resteady,thetavecsteady,nstepsteady,Fx_basesteady,Fy_basesteady,Fz_basesteady,
Mx_basesteady,My_basesteady,Mz_basesteady,powersteady,power2steady,torquesteady = OWENSAero.steadyTurb()
tsteady = (thetavecsteady.+(270/360*2*pi))./omega

# PyPlot.figure()
# PyPlot.plot(myazi,Xpbase,"k-",label="Xp")
# PyPlot.plot(myazi,Xpbase2,"k--",label="Xp2")
# PyPlot.plot(myazi,Ypbase,"b-",label="Yp")
# PyPlot.plot(myazi,Ypbase2,"b--",label="Yp2")
# PyPlot.legend()

freac = mytorque[1,:]
istart = 14
istart2 = 14
iend2 = numTS-1
PyPlot.ion()
PyPlot.figure()
PyPlot.plot(t[istart:iend2],-freac[istart:iend2]/1000 ,color=plot_cycle[1],label="Freac")
PyPlot.plot([t[1],t[end]],mean(-freac[istart:iend2]/1000).*ones(2),"--" ,color=plot_cycle[1],label="Freacmean")
# PyPlot.plot(t,torqueDriveShaft/1000 ,color=plot_cycle[4],label="torqueDriveShaft")
# PyPlot.plot([t[1],t[end]],mean(torqueDriveShaft/1000).*ones(2),"--" ,color=plot_cycle[4],label="torqueDriveShaftmean")
PyPlot.plot(t[istart2:end],Mz_base[istart2:end]./1000,"-" ,color=plot_cycle[2],label="aeroOnly")
PyPlot.plot([t[1],t[end]],mean(Mz_base[istart2:end]./1000).*ones(2),"--" ,color=plot_cycle[2],label="aeroOnlymean")
PyPlot.plot(tsteady,Mz_basesteady./1000 ,color=plot_cycle[3],label="aeroOnlysteady")
PyPlot.plot([tsteady[1],tsteady[end]],mean(Mz_basesteady[:]./1000).*ones(2),"--" ,color=plot_cycle[3],label="aeroOnlymeansteady")
PyPlot.xlabel("Time (s)")
# PyPlot.xlim([0,5])
PyPlot.ylabel("Torque (kN-m)")
PyPlot.legend()


# #######################################################################################
# #######################################################################################
# ############################# Now the other way #######################################
# #######################################################################################
# #######################################################################################
# import QuadGK
# import OWENS
# import OWENSFEA
# import OWENSAero
# import FLOWMath
# import DelimitedFiles
# using Statistics:mean
# using Statistics

# import PyPlot
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(4.5, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# # function runprofilefunction()
# path = runpath = splitdir(@__FILE__)[1]

# # Inp = OWENS.MasterInput("$path/SNL34m_Inputs.yml")
# # Inp = OWENS.MasterInput("$path/SNL34m_InputsAeroDyn.yml")

# nothing

# # Unpack inputs, or you could directly input them here and bypass the file 

# verbosity = 1

# turbineType = "Darrieus"
# eta = 0.5
# Nbld = 2
# towerHeight = 0.5
# rho = 1.225
# Vinf = 10.1
# controlStrategy = "constantRPM"
# RPM = 34.0
# Nslices = 35
# ntheta = 30
# structuralModel = "ROM"
# ntelem = 10
# nbelem = 60
# ncelem = 10
# nselem = 5
# ifw = false
# AModel = "DMS"
# windINPfilename = "$(path)/data/turbsim/115mx115m_30x30_20.0msETM.bts"
# ifw_libfile = nothing
# Blade_Height = 41.9
# Blade_Radius = 17.1
# numTS = 400
# delta_t = 0.05
# NuMad_geom_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerGeom.csv"
# NuMad_mat_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerMaterials.csv"
# NuMad_geom_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mGeomBlades.csv"
# NuMad_mat_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mMaterials.csv"
# NuMad_geom_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mGeomStruts.csv"
# NuMad_mat_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mMaterials.csv"
# adi_lib = nothing

# adi_rootname = "/SNL34m"

# ##############################################
# # Setup
# #############################################

# SNL34m_5_3_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Vinf.csv",',',skipstart = 0)
# SNL34m_5_3_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_RPM.csv",',',skipstart = 0)
# SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque.csv",',',skipstart = 0)


# new_t = LinRange(SNL34m_5_3_RPM[1,1],SNL34m_5_3_RPM[end,1],100)
# new_RPM = FLOWMath.akima(SNL34m_5_3_RPM[:,1],SNL34m_5_3_RPM[:,2],new_t)

# Vinf_spec = FLOWMath.akima(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],new_t)

# offsetTime = 20.0 # seconds
# tocp = [0.0;new_t.+offsetTime; 1e6]
# Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]]./60 .*0 .+33.92871/60
# t_Vinf = [0;new_t;1e6]
# Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
# tocp_Vinf = [0.0;t_Vinf.+offsetTime; 1e6]
# Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]

# controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# # z_shape = collect(LinRange(0,41.9,length(x_shape)))
# z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
# x_shape1 = [0.0;controlpts;0.0]
# z_shape = collect(LinRange(0,41.9,60))
# x_shape = FLOWMath.akima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
# toweroffset = 4.3953443986241725
# SNL34_unit_xz = [x_shape;;z_shape]
# SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
# SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])
# SNL34Z = SNL34z.*Blade_Height
# SNL34X = SNL34x.*Blade_Radius

# shapeZ = SNL34Z#collect(LinRange(0,H,Nslices+1))
# shapeX = SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

# shapeX_spline = FLOWMath.Akima(SNL34Z, SNL34X)
# RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, Blade_Height, atol=1e-10)
# RefArea = RefArea_half*2

# mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
# stiff_twr, stiff_bld,bld_precompinput,
# bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
# twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
# mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
#     rho,
#     Nslices,
#     ntheta,
#     RPM,
#     Vinf,
#     eta,
#     B = Nbld,
#     H = Blade_Height,
#     R = Blade_Radius,
#     shapeZ,
#     shapeX,
#     ifw,
#     delta_t,
#     numTS,
#     adi_lib,
#     adi_rootname,
#     AD15hubR = 0.0,
#     windINPfilename,
#     ifw_libfile,
#     NuMad_geom_xlscsv_file_twr,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS_Cables.csv",
#     NuMad_mat_xlscsv_file_twr,# = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv",
#     NuMad_geom_xlscsv_file_bld,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS.csv",
#     NuMad_mat_xlscsv_file_bld,# = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
#     NuMad_geom_xlscsv_file_strut,
#     NuMad_mat_xlscsv_file_strut,
#     Htwr_base=towerHeight,
#     ntelem, 
#     nbelem, 
#     ncelem,
#     nselem,
#     joint_type = 0,
#     strut_twr_mountpoint = [0.03,0.97],
#     strut_bld_mountpoint = [0.03,0.97],
#     AModel, #AD, DMS, AC
#     DSModel="BV",
#     RPI=true,
#     cables_connected_to_blade_base = true,
#     angularOffset = pi/2,
#     meshtype = turbineType)

# ##########################################
# ############ AeroElastic #################
# ##########################################

# top_idx = 23#Int(myjoint[7,2])
# pBC = [1 1 0
# 1 2 0
# 1 3 0
# 1 4 0
# 1 5 0
# 1 6 0
# top_idx 1 0
# top_idx 2 0
# top_idx 3 0
# top_idx 4 0
# top_idx 5 0]

# if AModel=="AD"
#     AD15On = true
# else
#     AD15On = false
# end

# inputs = OWENS.Inputs(;analysisType = structuralModel,
#     tocp,
#     Omegaocp,
#     tocp_Vinf,
#     Vinfocp,
#     numTS,
#     delta_t,
#     AD15On,
#     aeroLoadsOn = 2,
#     turbineStartup = 0,
#     generatorOn = false,
#     useGeneratorFunction = false,
#     driveTrainOn = false,
#     JgearBox = 250.0,#(2.15e3+25.7)/12*1.35582*100,
#     gearRatio = 1.0,
#     gearBoxEfficiency = 1.0,
#     driveShaftProps = OWENS.DriveShaftProps(10000,1.5e2), #8.636e5*1.35582*0.6
#     OmegaInit = Omegaocp[1])

# nothing

# # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

# feamodel = OWENS.FEAModel(;analysisType = structuralModel,
#     joint = myjoint,
#     platformTurbineConnectionNodeNumber = 1,
#     pBC,
#     nlOn=false,
#     numNodes = mymesh.numNodes,
#     numModes = 200,
#     RayleighAlpha = 0.05,
#     RayleighBeta = 0.05,
#     iterationType = "DI")



# # inputs.Omegaocp = Omegaocp
# # inputs.OmegaInit = Omegaocp[1]
# # inputs.Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
# # feamodel.nlOn = false
# # feamodel.analysisType = "ROM"
# # inputs.analysisType = "ROM"

# # eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,omegaHist,genTorque,torqueDriveShaft,aziHist,uHist,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.run34m(inputs,feamodel,mymesh,myel,
# # aeroForces,deformAero;steady=false,system,assembly,VTKFilename="$path/vtk/NormalOperation")

# t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
# FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
# epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
# topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;
# topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,system,assembly)


# # Reset and run aero only
# # OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSRvec[1],Vinf_array[1];
# #     eta = 0.5,
# #     rho,
# #     mu = 1.7894e-5,
# #     ntheta = 30,
# #     Nslices,
# #     RPI=true,
# #     ifw = false,
# #     DSModel = "BV",
# #     AModel = "DMS",
# #     tau = [1e-5,1e-5],
# #     afname = airfoils)

# # t = 0:0.05:30

# RPMsetpoint = 34.0
# omega = RPMsetpoint/60*2*pi
# Mz_base = zero(t)
# Xpbase = zero(t)
# Ypbase = zero(t)
# Xpbase2 = zero(t)
# Ypbase2 = zero(t)
# myazi = zero(t)
# for (i,myt) in enumerate(t)
#     azi = omega*myt + 270/360*2*pi +0.1780235837034216
#     myazi[i] = azi
#     CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,_,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base[i],power,power2,_,z3Dnorm,delta,Xp,Yp = OWENSAero.AdvanceTurbineInterpolate(myt;azi,alwaysrecalc=false)
#     for ibld = 1:length(Xp[:,1,1])
#         Xpbase[i] += OWENSAero.trapz(z3Dnorm.*Blade_Height,Xp[ibld,:,end])
#         Ypbase[i] += OWENSAero.trapz(z3Dnorm.*Blade_Height,Yp[ibld,:,end])
#         Xpbase2[i] = Xpbase[i]*cos(-azi) + Ypbase[i]*sin(-azi)
#         Ypbase2[i] = -Xpbase[i]*sin(-azi) + Ypbase[i]*cos(-azi)
#     end
# end


# CPsteady,Rpsteady,Tpsteady,Zpsteady,alphasteady,cl_afsteady,cd_afsteady,Vlocsteady,Resteady,thetavecsteady,nstepsteady,Fx_basesteady,Fy_basesteady,Fz_basesteady,
# Mx_basesteady,My_basesteady,Mz_basesteady,powersteady,power2steady,torquesteady = OWENSAero.steadyTurb()
# tsteady = (thetavecsteady.+(270/360*2*pi))./omega

# # PyPlot.figure()
# # PyPlot.plot(myazi,Xpbase,"k-",label="Xp")
# # PyPlot.plot(myazi,Xpbase2,"k--",label="Xp2")
# # PyPlot.plot(myazi,Ypbase,"b-",label="Yp")
# # PyPlot.plot(myazi,Ypbase2,"b--",label="Yp2")
# # PyPlot.legend()


# ##########################################
# #### Torque Plot
# ##########################################

# SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque2.csv",',',skipstart = 0)

# PyPlot.ion()
# PyPlot.figure()
# # PyPlot.plot(t.-offsetTime,torqueDriveShaft/1000 ,color=plot_cycle[1],label="Simulated Drive Shaft")
# # PyPlot.plot([-20,80],ones(2).*mean(torqueDriveShaft/1000) ,color=plot_cycle[2],label="Simulated Drive Shaft Mean")
# PyPlot.plot(t.-offsetTime,-FReactionHist[:,6]/1000 ,color=plot_cycle[3],label="Reaction Force")
# PyPlot.plot([-20,80],ones(2).*mean(-FReactionHist[:,6]/1000) ,color=plot_cycle[3],label="Reaction Force Mean")
# usedLogic = SNL34m_5_3_Torque[:,1].<100
# # PyPlot.plot(SNL34m_5_3_Torque[usedLogic,1],SNL34m_5_3_Torque[usedLogic,2],"k-",label="Experimental")
# PyPlot.plot(t.-offsetTime,topFexternal_hist[:,6]/1000,"k--",label="Aero Torque1")
# PyPlot.plot([-20,80],ones(2).*mean(topFexternal_hist[:,6]/1000),"k--",label="Aero Mean1")
# # PyPlot.plot(t.-offsetTime,OmegaHist.*60,"k--",label="OmegaHist RPM")
# PyPlot.plot(t.-offsetTime,Mz_base./1000,"+-" ,color=plot_cycle[2],label="aeroOnlydirect")
# PyPlot.plot([t[1],t[end]].-offsetTime,mean(Mz_base[15:end]./1000).*ones(2),"+-" ,color=plot_cycle[2],label="aeroOnlymeandirect")
# PyPlot.plot(tsteady.-offsetTime,Mz_basesteady./1000 ,color=plot_cycle[1],label="aeroOnlysteadydirect")
# PyPlot.plot([tsteady[1],tsteady[end]].-offsetTime,mean(Mz_basesteady[:]./1000).*ones(2),"--" ,color=plot_cycle[1],label="aeroOnlymeansteadydirect")
# PyPlot.xlabel("Time (s)")
# PyPlot.xlim([-20,0])
# PyPlot.ylabel("Torque (kN-m)")
# PyPlot.legend()#loc = (0.06,1.0),ncol=2)
# # PyPlot.savefig("$(path)/../figs/34m_fig5_32Way.pdf",transparent = true)
