
using Statistics:mean
import HDF5
import DelimitedFiles
import FLOWMath
import OWENSFEA
import GXBeam
import OWENS
import OWENSAero


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

##############################################
# Setup
#############################################
#Put in one place so its not repeated for all of the analyses
include("$(path)/34mSetup.jl")
numTS = 100
RPMsetpoint = 34.0
omega = RPMsetpoint/60*2*pi
nTSR = 5
# TSRvec = LinRange(2.43,15.22,nTSR) #[5.843]
offsetTime = 1e-6
t_Vinf = LinRange(0,1e6,10)
# Vinf_array = omega*radius./TSRvec
Vinf_array = [5.0]#collect(LinRange(5.0,25.0,nTSR))
TSRvec = omega*radius./Vinf_array
mytorque = zeros(nTSR,numTS-1)
uHist = []
torqueDriveShaft = []
t = 0.0
for (iTSR,TSRspec) in enumerate(TSRvec)
    TSRspec = TSRvec[iTSR]
    Vinf_spec = ones(10).*Vinf_array[iTSR]
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



    Vinf = Vinf_spec[1]#mean(SNL34m_5_3_Vinf[:,2])
    TSR = omega*radius./Vinf
    # windpower = 0.5*rho*Vinf.^3*RefArea
    ntheta = 30#176
    global rho = 1.225 #override since it was corrected for sea level
    OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSR,Vinf;
        eta = 0.5,
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

    dt = 1/(RPM/60*ntheta)

    aeroForcesDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=false)

    tocp_Vinf = [0.0;t_Vinf.+offsetTime; 1e6]
    Omegaocp = zero(tocp_Vinf) .+RPMsetpoint/60
    Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]

    model = OWENS.Inputs(;analysisType = "ROM",
        outFilename = "none",
        tocp = tocp_Vinf,#SNL34m_5_3_RPM[:,1],#[0.0,10.0,100000.1],
        Omegaocp,#SNL34m_5_3_RPM[:,2]./ 60,#[RPM,RPM,RPM] ./ 60,
        tocp_Vinf,
        Vinfocp,
        numTS,
        delta_t = 0.02,#dt,
        aeroLoadsOn = 2,
        turbineStartup = 0,
        generatorOn = false,
        useGeneratorFunction = false,
        driveTrainOn = false,
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
        nlOn = false,
        numNodes = mymesh.numNodes,
        numModes = 200,
        RayleighAlpha = 0.1,
        RayleighBeta = 0.1,
        iterationType = "DI")

    # Get Gravity Loads
    model.Omegaocp = model.Omegaocp.*0.0
    model.OmegaInit = model.OmegaInit.*0.0
    model.Vinfocp = model.Vinfocp.*0.0
    feamodel.nlOn = true

    # Returns data filled with e.g. eps[Nbld,N_ts,Nel_bld]
    eps_x_grav,eps_z_grav,eps_y_grav,kappa_x_grav,kappa_y_grav,kappa_z_grav,tin,FReactionHist_grav = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,OWENSAero.deformTurb;steady=true)
    global t = tin
    #####
    ###****** SAND-88-1144 Specifies Bending Strains and Axial Strains Separate ****
    #####

    Ealuminum = plyprops.plies[end].e1
    flatwise_stress1grav = (kappa_y_grav[1,end,2:end].* thickness .+ eps_x_grav[1,end,2:end]) .* Ealuminum
    flatwise_stress2grav = (kappa_y_grav[2,end,1:end-1].* thickness .+ eps_x_grav[2,end,1:end-1]) .* Ealuminum
    lag_stress1grav = (kappa_z_grav[1,end,2:end].* thickness_lag .+ eps_x_grav[1,end,2:end]) .* Ealuminum
    lag_stress2grav = (kappa_z_grav[2,end,1:end-1].* thickness_lag .+ eps_x_grav[2,end,1:end-1]) .* Ealuminum

    model.Omegaocp = Omegaocp
    model.OmegaInit = Omegaocp[1]
    model.Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
    feamodel.nlOn = false
    feamodel.analysisType = "ROM"
    model.analysisType = "ROM"

    global uHist
    global torqueDriveShaft
    eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,torqueDriveShaft,aziHist,uHist = runowens(model,feamodel,mymesh,myel,aeroForcesDMS,OWENSAero.deformTurb;steady=false,system=system,assembly=assembly)
    # Juno.@enter runowens(model,feamodel,mymesh,myel,aeroForcesDMS,OWENSAero.deformTurb;steady=false,system=system,assembly=assembly)

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

    # # Plots
    # PyPlot.figure()
    # PyPlot.plot(t[1:end-1].-offsetTime,flatwise_stress1[:,end-5]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
    # PyPlot.plot(SNL34m_5_4_FlatwiseStress[:,1],SNL34m_5_4_FlatwiseStress[:,2],"k-",label = "Experimental")
    # PyPlot.xlabel("Time (s)")
    # PyPlot.ylabel("Flatwise Stress (MPa)")
    # # PyPlot.xlim([0,SNL34m_5_4_FlatwiseStress[end,1]])
    # PyPlot.legend()
    # PyPlot.savefig("$(path)/../figs/34m_fig5_4_NormalOperation_flapwise_Blade2Way_TSR_$TSRspec.pdf",transparent = true)
    #
    # SNL34m_5_4_LeadLagStress = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.4AML.csv",',',skipstart = 1)
    #
    # PyPlot.figure()
    # PyPlot.plot(t[1:end-1].-offsetTime,lag_stress1[:,end-5]./1e6,"-",color=plot_cycle[1],label = "OWENS Blade 1")
    # PyPlot.plot(SNL34m_5_4_LeadLagStress[:,1],SNL34m_5_4_LeadLagStress[:,2],"k-",label = "Experimental")
    # PyPlot.xlabel("Time (s)")
    # PyPlot.ylabel("Lead-Lag Stress (MPa)")
    # # PyPlot.xlim([0,SNL34m_5_4_LeadLagStress[end,1]])
    # PyPlot.legend()
    # PyPlot.savefig("$(path)/../figs/34m_fig5_4_NormalOperation_LeadLag_Blade2Way_TSR_$TSRspec.pdf",transparent = true)
    #
    #
    # ##########################################
    # #### Torque Plot
    # ##########################################
    #
    # PyPlot.ion()
    # PyPlot.figure()
    # PyPlot.plot(t.-offsetTime,-FReactionHist[:,6]/1000 ,color=plot_cycle[1],label="OWENSFreact")
    # PyPlot.plot(t.-offsetTime,mean(-FReactionHist[1000:end,6]/1000).+zero(t),"--" ,color=plot_cycle[1],label="OWENSFreactMean")
    # # PyPlot.plot(t.-offsetTime,FhatHist/1000 ,color=plot_cycle[2],label="FhatHist")
    # # PyPlot.plot(t.-offsetTime,torqueDriveShaft/1000 ,color=plot_cycle[2],label="torqueDriveShaft")
    # # PyPlot.plot(t.-offsetTime,mean(torqueDriveShaft[1000:end]/1000).+zero(t),"--" ,color=plot_cycle[2],label="torqueDriveShaftMean")
    # # PyPlot.plot(t.-offsetTime,genTorque/1000 ,color=plot_cycle[4],label="genTorque")
    # PyPlot.xlabel("Time (s)")
    # PyPlot.xlim([0,100])
    # PyPlot.ylabel("Torque (kN-m)")
    # PyPlot.legend()
    # PyPlot.savefig("$(path)/../figs/34m_fig5_32Way_TSR_$TSRspec.pdf",transparent = true)

    mytorque[iTSR,:] = FReactionHist[:,6]#torqueDriveShaft
end

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

# Reset and run aero only
OWENSAero.setupTurb(SNL34X,SNL34Z,B,chord,TSRvec[1],Vinf_array[1];
    eta = 0.5,
    rho,
    mu = 1.7894e-5,
    ntheta = 30,
    Nslices,
    RPI=true,
    ifw = false,
    DSModel = "BV",
    AModel = "DMS",
    tau = [1e-5,1e-5],
    afname = airfoils)

# t = 0:0.05:30
Mz_base = zero(t)
Xpbase = zero(t)
Ypbase = zero(t)
Xpbase2 = zero(t)
Ypbase2 = zero(t)
Zp = zero(t)
myazi = zero(t)
for (i,myt) in enumerate(t)
    azi = omega*myt + 270/360*2*pi +0.1780235837034216
    myazi[i] = azi
    CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,ntheta,Fx_base,Fy_base,Fz_base,Mx_base,My_base,Mz_base[i],power,power2,_,z3Dnorm,delta,Xp,Yp = OWENSAero.AdvanceTurbineInterpolate(myt;azi,alwaysrecalc=false)
    for ibld = 1:length(Xp[:,1,1])
        Xpbase[i] += OWENSAero.trapz(z3Dnorm.*height,Xp[ibld,:,end])
        Ypbase[i] += OWENSAero.trapz(z3Dnorm.*height,Yp[ibld,:,end])
        Xpbase2[i] = Xpbase[i]*cos(-azi) + Ypbase[i]*sin(-azi)
        Ypbase2[i] = -Xpbase[i]*sin(-azi) + Ypbase[i]*cos(-azi)
    end
end



CPsteady,Rpsteady,Tpsteady,Zpsteady,alphasteady,cl_afsteady,cd_afsteady,Vlocsteady,Resteady,thetavecsteady,nstepsteady,Fx_basesteady,Fy_basesteady,Fz_basesteady,
Mx_basesteady,My_basesteady,Mz_basesteady,powersteady,power2steady,torquesteady = OWENSAero.steadyTurb()
tsteady = (thetavecsteady.+(270/360*2*pi))./omega

PyPlot.figure()
PyPlot.plot(myazi,Xpbase,"k-",label="Xp")
PyPlot.plot(myazi,Xpbase2,"k--",label="Xp2")
PyPlot.plot(myazi,Ypbase,"b-",label="Yp")
PyPlot.plot(myazi,Ypbase2,"b--",label="Yp2")
PyPlot.legend()

freac = mytorque[1,:]
istart = 55
PyPlot.ion()
PyPlot.figure()
PyPlot.plot(t[istart:end],-freac[istart:end]/1000 ,color=plot_cycle[1],label="Freac")
PyPlot.plot([t[1],t[end]],mean(-freac[istart:end]/1000).*ones(2),"--" ,color=plot_cycle[1],label="Freacmean")
PyPlot.plot(t,torqueDriveShaft/1000 ,color=plot_cycle[4],label="torqueDriveShaft")
PyPlot.plot([t[1],t[end]],mean(torqueDriveShaft/1000).*ones(2),"--" ,color=plot_cycle[4],label="torqueDriveShaftmean")
PyPlot.plot(t,Mz_base./1000,"-" ,color=plot_cycle[2],label="aeroOnly")
PyPlot.plot([t[1],t[end]],mean(Mz_base[15:end]./1000).*ones(2),"--" ,color=plot_cycle[2],label="aeroOnlymean")
PyPlot.plot(tsteady,Mz_basesteady./1000 ,color=plot_cycle[3],label="aeroOnlysteady")
PyPlot.plot([tsteady[1],tsteady[end]],mean(Mz_basesteady[:]./1000).*ones(2),"--" ,color=plot_cycle[3],label="aeroOnlymeansteady")
PyPlot.xlabel("Time (s)")
# PyPlot.xlim([0,5])
PyPlot.ylabel("Torque (kN-m)")
PyPlot.legend()
# PyPlot.savefig("$(path)/../figs/34m_fig5_32Way_TSR_$TSRspec.pdf",transparent = true)
#
# PyPlot.figure()
# for its = 1:length(uHist[1,:])
# # its = length(uHist[1,:])
#     PyPlot.cla()

#     displ = uHist[:,its]
#     deformFact = 1.0
#     disp_x_OW = [displ[i] for i = 1:6:length(displ)]
#     disp_y_OW = [displ[i] for i = 2:6:length(displ)]
#     disp_z_OW = [displ[i] for i = 3:6:length(displ)]
#     disp_curv1_OW = [displ[i] for i = 4:6:length(displ)]
#     disp_curv2_OW = [displ[i] for i = 5:6:length(displ)]
#     disp_curv3_OW = [displ[i] for i = 6:6:length(displ)]
#     # disp_r = sqrt.(disp_x_OW.^2 .+ disp_y_OW.^2)
#     PyPlot.plot(mymesh.x,mymesh.z,"k-",label="Undeformed xz")
#     PyPlot.plot(mymesh.x,mymesh.y,"k-",label="Undeformed zy")
#     PyPlot.plot(mymesh.y,mymesh.z,"k-",label="Undeformed yz")
#     # PyPlot.plot((mymesh.x+disp_r*deformFact),(mymesh.z+disp_z_OW*deformFact),"--",color=plot_cycle[1],label="OWENSxz")
#     # PyPlot.plot(((mymesh.x+disp_r)./mymesh.x),(mymesh.z+disp_z_OW*deformFact),"--",color=plot_cycle[2],label="Diff")
#     PyPlot.plot((mymesh.x+disp_x_OW*deformFact),(mymesh.z+disp_z_OW*deformFact),"--",color=plot_cycle[1],label="OWENSxz")
#     PyPlot.plot((mymesh.x+disp_x_OW*deformFact),(mymesh.y+disp_y_OW*deformFact),"--",color=plot_cycle[2],label="OWENSxy")
#     PyPlot.plot((mymesh.y+disp_y_OW*deformFact),(mymesh.z+disp_z_OW*deformFact),"--",color=plot_cycle[3],label="OWENSyz")
#     PyPlot.legend()
#     PyPlot.xlabel("x-position (m)")
#     PyPlot.ylabel("y-position (m)")
#     sleep(0.01)

#     # println(maximum(sqrt.(disp_curv3_OW.^2 .+disp_curv2_OW.^2 .+disp_curv1_OW.^2))*180/pi)
# end
# PyPlot.axis("equal")
# PyPlot.xlim([0,0.5])
# PyPlot.savefig("$(path)/../figs/GravCentDisp.pdf",transparent = true)
