import OWENS
import OWENSAero
import DelimitedFiles
using Statistics:mean
using Test
import FLOWMath
import HDF5

# path = runpath = splitdir(@__FILE__)[1]
runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]

nothing

turbineType = "H-VAWT" # turbine type, for the automatic meshing
Vinf = 1.2 # inflow velocity
TSRrange = [3.5]#LinRange(1.0,5.0,2) range of tip speed ratios
Nslices = 20 # vertical discretizations if DMS or AC aero model
ntheta = 30 # azimuthal discretizations if DMS or AC aero model
structuralModel = "TNB"
ntelem = 100 # tower elements
nbelem = 30 # blade elements
nselem = 10 # strut elements
ifw = false # use inflow wind, if DMS or AC aero model
numTS = 21#321 # number of simulation time steps
delta_t = 0.01 # simulation time step spacing
adi_lib = nothing#"$path/../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding"
adi_rootname = "$path/RM2" # path and name that all the aerodyn files are saved with
VTKsaveName = "$path/vtk/RM2_medium" # path and name that all OWENS VTK files are saved with
tsave_idx = 1:1:numTS-1 #you don't have to save every timestep in VTK
ifw_libfile = nothing#"$path/../../../openfast/build/modules/inflowwind/libifw_c_binding"
fluid_density = 1000.0
fluid_dyn_viscosity = 1.792E-3
AddedMass_Coeff_Ca = 1.0 #For structural side added mass
Aero_Buoyancy_Active = true # For buoyancy forcing, handled by the OWENSAero module
AeroModel = "DMS"
verbosity = 1 # verbosity level where higher is more
if AeroModel=="AD"
    AD15On = true
else
    AD15On = false
end

eta = 0.5 # blade-strut mount point ratio, fraction from leading edge
number_of_blades = 3 # number of blades
towerHeight = 0.2165 # height of the tower past the blades
Blade_Height = 0.807 # height of the blades
Blade_Radius = 0.5375 # radius of the turbine

NuMad_geom_xlscsv_file_twr = "$path/data_RM2/TowerGeom.csv" #
NuMad_mat_xlscsv_file_twr = "$path/data_RM2/TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$path/data_RM2/GeomBlades$AeroModel.csv"
NuMad_mat_xlscsv_file_bld = "$path/data_RM2/materials_NuMAD.csv"
NuMad_geom_xlscsv_file_strut = "$path/data_RM2/GeomStruts$AeroModel.csv"
NuMad_mat_xlscsv_file_strut = "$path/data_RM2/materials_NuMAD.csv"

nothing

WindType = 1
windINPfilename = "$path/data_RM2/3mx3m1pt2msNTM.bts"

nothing

CP = zeros(length(TSRrange))
mymesh = []
myel = []
system = []
assembly = []
sections = []
myjoint = []
pBC = []
# for (iTSR,TSR) in enumerate(collect(TSRrange))
TSR = TSRrange[1]
iTSR = 1
    global Vinf
    global mymesh
    global myel
    global system
    global assembly
    global sections
    global myjoint
    global pBC

    omega = Vinf/Blade_Radius*TSR
    RPM = omega * 60 / (2*pi)

    println(RPM)

    nothing

    shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
    helix_angle = 0.0#-pi/4
    shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius
    shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius

    mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
    mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
        rho=fluid_density,
        mu=fluid_dyn_viscosity,
        Nslices,
        ntheta,
        RPM,
        Vinf,
        eta,
        B=number_of_blades,
        H = Blade_Height,
        R = Blade_Radius,
        shapeZ,
        shapeX,
        shapeY,
        ifw,
        WindType,
        delta_t,
        numTS,
        adi_lib,
        adi_rootname,
        windINPfilename,
        ifw_libfile,
        NuMad_geom_xlscsv_file_twr,
        NuMad_mat_xlscsv_file_twr,
        NuMad_geom_xlscsv_file_bld,
        NuMad_mat_xlscsv_file_bld,
        NuMad_geom_xlscsv_file_strut,
        NuMad_mat_xlscsv_file_strut,
        Htwr_base=towerHeight,
        Htwr_blds = Blade_Height+towerHeight,
        ntelem,
        nbelem,
        nselem,
        joint_type = 0,
        strut_twr_mountpoint = [0.5],
        strut_bld_mountpoint = [0.5],
        AeroModel, #AD, DMS, AC
        DynamicStallModel="BV",
        Aero_AddedMass_Active = false,
        Aero_RotAccel_Active = false,
        AddedMass_Coeff_Ca,
        Aero_Buoyancy_Active,
        verbosity,
        meshtype = turbineType)

    nothing

    # # PyPlot.figure()
    # # for idot = 1:length(sectionPropsArray[170].xaf)
    # #     PyPlot.scatter(sectionPropsArray[170].xaf[idot],sectionPropsArray[170].yaf[idot])
    # #     sleep(0.001)
    # # end
    #
    # ## This plots the mesh and node numbering of the resulting mesh and overlays the joint connections
    #
    # PyPlot.figure()
    # for icon = 1:length(mymesh.conn[:,1])
        # idx1 = mymesh.conn[icon,1]
        # idx2 = mymesh.conn[icon,2]
        # PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
        # PyPlot.plot3D([1,1],[1,1],[1,1],"k.-")
        # PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    # end
    # for ijoint = 1:length(myjoint[:,1])
        # idx2 = Int(myjoint[ijoint,2])
        # idx1 = Int(myjoint[ijoint,3])
        # PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
        # PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",color="r",ha="center",va="center")
        # PyPlot.text3D(mymesh.x[idx2].+rand()/30,mymesh.y[idx2].+rand()/30,mymesh.z[idx2].+rand()/30,"$idx2",color="r",ha="center",va="center")
    # end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.zlabel("z")
    # PyPlot.axis("equal")

    nothing

    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0
    AD15bldElIdxRng[1,2]-1 1 0
    AD15bldElIdxRng[1,2]-1 2 0
    AD15bldElIdxRng[1,2]-1 3 0
    AD15bldElIdxRng[1,2]-1 4 0
    AD15bldElIdxRng[1,2]-1 5 0]

    nothing

    tocp = [0.0;10.0; 1e6]
    Omegaocp = [RPM; RPM; RPM]./60 #control inputs

    tocp_Vinf = [0.0;10.0; 1e6]
    Vinfocp = [Vinf;Vinf;Vinf]

    inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
    tocp,
    dataOutputFilename = "./InitialDataOutputs_scripting.out",
    Omegaocp,
    tocp_Vinf,
    Vinfocp,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = 2,
    MAXITER = 20)

    nothing

    FEAinputs = OWENS.FEAModel(;analysisType = structuralModel,
    dataOutputFilename = "./InitialDataOutputs_scripting.out",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    numModes = 70,
    gravityOn = [0,0,9.81], #positive since the turbine is upside down
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.005,
    RayleighBeta = 0.005,
    AddedMass_Coeff_Ca,
    iterationType = "DI")

    nothing

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
    topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=FEAinputs,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,turbsimfile = windINPfilename)

    area = Blade_Height*2*Blade_Radius
    full_rev_N_timesteps = round(Int,RPM/60/delta_t)*2
    if full_rev_N_timesteps>numTS
        idx_start = 1
    else
        idx_start = numTS-full_rev_N_timesteps
    end
    CP[iTSR] = mean(FReactionHist[idx_start:end,6].*OmegaHist[idx_start:end]*2*pi)/(0.5*fluid_density*mean(Vinfocp)^3*area)
    TSR = mean(OmegaHist*2*pi*Blade_Radius/mean(Vinfocp))
    ReD = fluid_density*mean(Vinfocp)*Blade_Radius*2/fluid_dyn_viscosity

    nothing

    ##########################################
    #### Ultimate Failure #####
    ##########################################

    massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
    SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
    topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
    topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,number_of_blades,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,delta_t,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing)

    OWENS.outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)

    # OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
    #     epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    #     FReactionHist,topFexternal_hist;tsave_idx)


# end

# PyPlot.figure("CP2")
# PyPlot.plot(TSRrange,CP,".-",color=plot_cycle[2],label="OWENS Aero, 1.3 RE_d (No Added Mass)") #,color=color_cycle[2]
# RM2_0_538D_RE_D_1_3E6 = DelimitedFiles.readdlm("$(path)/data_RM2/RM2_0.538D_RE_D_1.3E6.csv", ',',Float64)
# PyPlot.plot(RM2_0_538D_RE_D_1_3E6[:,1],RM2_0_538D_RE_D_1_3E6[:,2],"k-",label="Exp. 1.3e6 RE_d")
# PyPlot.legend()
# PyPlot.xlabel("TSR")
# PyPlot.ylabel("Cp")

Qinst = -FReactionHist[idx_start:end,6]
Qinst2 = topFexternal_hist[idx_start:end,6]

drag = FReactionHist[idx_start:end,1] #./ (0.5*fluid_density*mean(Vinfocp)^2*area)
drag2 = topFexternal_hist[idx_start:end,1] #./ (0.5*fluid_density*mean(Vinfocp)^2*area)

dat_strt = round(Int,160800/80*11.0)
dat_end = round(Int,160800/80*14.0)

# # Open the HDF5 file in read mode
# carriage_pos = nothing
# drag_left = nothing
# drag_right = nothing
# time = nothing
# torque_arm = nothing
# torque_trans = nothing
# turbine_angle = nothing
# # This data is located at https://figshare.com/articles/dataset/UNH_RM2_tow_tank_experiment_Raw_data/1302029, with the file matching that here
# c = HDF5.h5open("$path/data_RM2/Perf1.2b_16_nidata.h5", "r") do file
#     global carriage_pos = read(file,"data/carriage_pos")
#     global drag_left = read(file,"data/drag_left")
#     global drag_right = read(file,"data/drag_right")
#     global time = read(file,"data/time")
#     global torque_arm = read(file,"data/torque_arm")
#     global torque_trans = read(file,"data/torque_trans")
#     global turbine_angle = read(file,"data/turbine_angle")
# end



# PyPlot.figure()
# PyPlot.plot((t[idx_start:end].-t[idx_start]),Qinst,".-",color=plot_cycle[1],label="Reaction") #,color=color_cycle[2]
# # PyPlot.plot((t[idx_start:end].-t[idx_start]),-Qinst2,"x-",color=plot_cycle[2],label="Applied") #,color=color_cycle[2]
# PyPlot.plot(time[dat_strt:dat_end].-time[dat_strt],torque_trans[dat_strt:dat_end],"k-",label="Exp. ")
# PyPlot.legend()
# PyPlot.xlabel("Time (s)")
# PyPlot.ylabel("Q (instantaneous)")
#
# PyPlot.figure()
# # PyPlot.plot((t[idx_start:end].-t[idx_start]),drag,".-",color=plot_cycle[1],label="Reaction") #,color=color_cycle[2]
# PyPlot.plot((t[idx_start:end].-t[idx_start]),drag2,"x-",color=plot_cycle[2],label="Applied") #,color=color_cycle[2]
# PyPlot.plot(time[dat_strt:dat_end].-time[dat_strt],drag_left[dat_strt:dat_end],"k-",label="Exp. L")
# PyPlot.plot(time[dat_strt:dat_end].-time[dat_strt],drag_right[dat_strt:dat_end],"k--",label="Exp. R")
# PyPlot.legend()
# PyPlot.xlabel("Time (s)")
# PyPlot.ylabel("Drag (instantaneous)")

# nothing

# # Calculate mean Q and compare
# mean_Q_owens = mean(Qinst)
# mean_Q_exp = mean(torque_trans[dat_strt:dat_end])
# println("Percent Difference in Torque: $((mean_Q_owens-mean_Q_exp)/mean_Q_exp*100)")

# dat_strt = round(Int,160800/80*10)
# dat_end = round(Int,160800/80*10.5)
# amp_Q_exp = maximum(torque_trans[dat_strt:dat_end])-minimum(torque_trans[dat_strt:dat_end])
# amp_Q_owens = maximum(Qinst[end-round(Int,0.5/delta_t):end])-minimum(Qinst[end-round(Int,0.5/delta_t):end])
# println("Percent Difference in Amplitude: $((amp_Q_owens-amp_Q_exp)/amp_Q_exp*100)")

rotSpdArrayRPM = [0.0, 42.64]

FEAinputs = OWENS.FEAModel(;analysisType = "GX",
dataOutputFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
gravityOn = [0,0,9.81], #positive since the turbine is upside down
numNodes = mymesh.numNodes,
RayleighAlpha = 0.05,
RayleighBeta = 0.05,
AddedMass_Coeff_Ca,
iterationType = "DI")

freq2 = OWENS.AutoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    rotSpdArrayRPM,
    VTKsavename=VTKsaveName,
    saveModes = [1,3,5], #must be int
    saveRPM = [2], #must be int
    mode_scaling = 500.0,
    )
freqGX = [freq2[:,i] for i=1:2:FEAinputs.numModes-6-2]

nothing

# Now the Campbell diagram can be generated
# NperRevLines = 8
# PyPlot.figure()
# for i=1:NperRevLines
#     linex=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5]
#     liney=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5].*i./60.0
#     PyPlot.plot(linex,liney,"--k",linewidth=0.5)
#     PyPlot.annotate("$i P",xy=(0.95*linex[2],liney[2]+.05+(i-1)*.01))
# end
# PyPlot.grid()
# PyPlot.xlabel("Rotor Speed (RPM)")
# PyPlot.ylabel("Frequency (Hz)")
# PyPlot.plot(0,0,"k-",label="Experimental")
# PyPlot.plot(0,0,color=plot_cycle[1],"-",label="OWENS")
# PyPlot.legend()

# for i=1:1:FEAinputs.numModes
#        PyPlot.plot(rotSpdArrayRPM,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
# end
# PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
# PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
# PyPlot.ylim([0,40.0])

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
