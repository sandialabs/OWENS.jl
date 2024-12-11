
import OWENS
import OWENSAero
import DelimitedFiles
using Statistics:mean
using Test
import FLOWMath

import PyPlot
PyPlot.pygui(true)
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

# function runprofilefunction()
path = runpath = splitdir(@__FILE__)[1]

# Unpack inputs, or you could directly input them here and bypass the file 

verbosity = 1


naca0021_coords = DelimitedFiles.readdlm("$path/airfoils/circular.csv",',',skipstart=0)

PyPlot.figure()
PyPlot.plot(naca0021_coords[:,1],naca0021_coords[:,2],".-")
PyPlot.axis("equal")

analysisType = "unsteady"
turbineType = "H-VAWT"
eta = 0.5
Nbld = 3
towerHeight = 0.2165
Vinf = 1.2
controlStrategy = "constantRPM"
Nslices = 20
ntheta = 30
structuralModel = "GX"
ntelem = 100
nbelem = 30
ncelem = 10
nselem = 10
ifw = false
AeroModel = "DMS"
windINPfilename = "$path/3mx3m1pt2msNTM.bts"
# run(`$(OWENS.OWENSOpenFASTWrappers.turbsim()) $(windINPfilename[1:end-3])inp`)
ifw_libfile = nothing#"$path/../../openfast/build/modules/inflowwind/libifw_c_binding"
Blade_Height = 0.807
Blade_Radius = 0.5375
area = Blade_Height*2*Blade_Radius
numTS = 10
delta_t = 0.01
NuMad_geom_xlscsv_file_twr = "$path/TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$path/TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$path/GeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$path/Materials34m.csv"
NuMad_geom_xlscsv_file_strut = "$path/GeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$path/Materials34m.csv"
adi_lib = nothing#"$path/../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" 
adi_rootname = "$path/helical"

fluid_density = 1000.0
fluid_dyn_viscosity = 1.792E-3
number_of_blades = Nbld
WindType = 3

AddedMass_Coeff_Ca=0.0 #For structural side

TSRrange = LinRange(1.0,5.0,2)
CP = zeros(length(TSRrange))
mymesh = []
myel = []
system = []
assembly = []
sections = []
myjoint = []
pBC = []
for (iTSR,TSR) in enumerate(collect(TSRrange))
    global Vinf
    global mymesh
    global myel
    global system
    global assembly
    global sections
    global myjoint
    global pBC
    # global TSR 
    omega = Vinf/Blade_Radius*TSR  
    RPM = omega * 60 / (2*pi)

    println(RPM)
    ##############################################
    # Setup
    #############################################

    shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
    helix_angle = 0.0#-pi/4
    shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius#SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ) ones(length(shapeZ)).*Blade_Radius#
    shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius # zeros(length(shapeX))#

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
        H = Blade_Height, #windio
        R = Blade_Radius, #windio
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
        NuMad_geom_xlscsv_file_twr,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS_Cables.csv",
        NuMad_mat_xlscsv_file_twr,# = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv",
        NuMad_geom_xlscsv_file_bld,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS.csv",
        NuMad_mat_xlscsv_file_bld,# = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
        NuMad_geom_xlscsv_file_strut,
        NuMad_mat_xlscsv_file_strut,
        Htwr_base=towerHeight,
        Htwr_blds = Blade_Height+towerHeight,
        ntelem, 
        nbelem, 
        ncelem,
        nselem,
        joint_type = 0,
        c_mount_ratio = 0.05,
        strut_twr_mountpoint = [0.5],
        strut_bld_mountpoint = [0.5],
        AeroModel, #AD, DMS, AC
        DynamicStallModel="BV",
        Aero_AddedMass_Active = false,
        Aero_RotAccel_Active = false,
        AddedMass_Coeff_Ca,
        Aero_Buoyancy_Active = true,
        verbosity,
        RPI=true,
        cables_connected_to_blade_base = true,
        meshtype = turbineType)

    nothing

    # Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
    # this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
    # includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
    # to work with OWENS inputs as well.

    tower_base_props = myel.props[1]
    blade_tip_props = myel.props[AD15bldElIdxRng[1,1]]


    PyPlot.figure()
    for icon = 1:length(mymesh.conn[:,1])
        idx1 = mymesh.conn[icon,1]
        idx2 = mymesh.conn[icon,2]
        PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
        PyPlot.plot3D([1,1],[1,1],[1,1],"k.-")
        PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
        # sleep(0.1)
    end

    for ijoint = 1:length(myjoint[:,1])
        idx2 = Int(myjoint[ijoint,2])
        idx1 = Int(myjoint[ijoint,3])
        PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
        PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",color="r",ha="center",va="center")
        PyPlot.text3D(mymesh.x[idx2].+rand()/30,mymesh.y[idx2].+rand()/30,mymesh.z[idx2].+rand()/30,"$idx2",color="r",ha="center",va="center")
        # sleep(0.1)
    end
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.zlabel("z")
    PyPlot.axis("equal")

    nothing

    # Here we apply the boundary conditions.  For this case, with a regular cantelever tower, the tower base node which is 
    # 1 is constrained in all 6 degrees of freedom to have a displacement of 0.  You can change this displacement to allow for things
    # like pretension, and you can apply boundary conditions to any node.

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

    # There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    # PyPlot.figure()
    # xtest = LinRange(0,1,10)
    # PyPlot.plot(xtest,FLOWMath.akima([0,0.5,1],[10,12,12],xtest))



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
    MAXITER = 10)

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

    FEAinputs = OWENS.FEAModel(;analysisType = structuralModel,
    dataOutputFilename = "./InitialDataOutputs_scripting.out",
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

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time.

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
    topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=FEAinputs,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,turbsimfile = windINPfilename)

    CP[iTSR] = mean(FReactionHist[:,6].*OmegaHist*2*pi)/(0.5*fluid_density*mean(Vinfocp)^3*area)
    TSR = mean(OmegaHist*2*pi*Blade_Radius/mean(Vinfocp))
    ReD = fluid_density*mean(Vinfocp)*Blade_Radius*2/fluid_dyn_viscosity



    nothing

    # Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

    # This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
    # it also calculates analytical simply supported buckling safety factors

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
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

    OWENS.outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)

    if iTSR == 2

        azi=aziHist#./aziHist*1e-6
        VTKsaveName = "$path/vtk/RM2_medium"
        tsave_idx=1:1:numTS-1
        OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
            epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
            FReactionHist,topFexternal_hist;tsave_idx)

        nothing
    end

end

# tsr_aeroelastic = [1.0, 2.0, 3.0, 4.0, 5.0]
# cp_aeroelastic = [0.005925193025218385, 0.061526662120420406, 0.3977425517719758, 0.40540280139050694, 0.2658708292801572]

PyPlot.figure("CP2")
PyPlot.plot(TSRrange,-CP,"-",color=plot_cycle[2],label="OWENS Aero, 1.3 RE_d (No Added Mass)") #,color=color_cycle[2]
# # end
# RM2_0_538D_RE_D_0_9E6 = DelimitedFiles.readdlm("$(path)/RM2_0.538D_RE_D_0.9E6.csv", ',',Float64)
# RM2_0_538D_RE_D_1_3E6 = DelimitedFiles.readdlm("$(path)/RM2_0.538D_RE_D_1.3E6.csv", ',',Float64)
# # PyPlot.plot(RM2_0_538D_RE_D_0_9E6[:,1],RM2_0_538D_RE_D_0_9E6[:,2],"k--",label="Exp. 0.9e6 RE_d")
# PyPlot.plot(RM2_0_538D_RE_D_1_3E6[:,1],RM2_0_538D_RE_D_1_3E6[:,2],"k-",label="Exp. 1.3e6 RE_d")
# # PyPlot.plot(right_TSR,right_CP,"ko",label="Right Only Exp.")
PyPlot.legend()
PyPlot.xlabel("TSR")
PyPlot.ylabel("Cp")


# ffmpeg -i smallerhelical.%04d.png -vf palettegen=reserve_transparent=1 palette.png
# ffmpeg -framerate 30 -i smallerhelical.%04d.png -i palette.png -lavfi paletteuse=alpha_threshold=30 -gifflags -offsetting smallerhelical.gif

nothing



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
    VTKsavename="$path/campbellVTK/RM2",
    saveModes = [1,3,5], #must be int
    saveRPM = [2], #must be int
    mode_scaling = 500.0,
    )
freqGX = [freq2[:,i] for i=1:2:FEAinputs.numModes-6-2]


import PyPlot
PyPlot.close("all")
PyPlot.pygui(true)
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

#plot per rev lines
NperRevLines = 8
PyPlot.figure()
for i=1:NperRevLines
    linex=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5]
    liney=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5].*i./60.0
    PyPlot.plot(linex,liney,"--k",linewidth=0.5)
    PyPlot.annotate("$i P",xy=(0.95*linex[2],liney[2]+.05+(i-1)*.01))
end
PyPlot.grid()
PyPlot.xlabel("Rotor Speed (RPM)")
PyPlot.ylabel("Frequency (Hz)")
PyPlot.plot(0,0,"k-",label="Experimental")
PyPlot.plot(0,0,color=plot_cycle[1],"-",label="OWENS")
PyPlot.legend()
# PyPlot.ylim([0.0,0.8])
# PyPlot.savefig("$(path)/../figs/34mCampbell.pdf",transparent = true)

# Add to figure
for i=1:1:FEAinputs.numModes
       PyPlot.plot(rotSpdArrayRPM,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
end
PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
PyPlot.ylim([0,40.0])
# PyPlot.savefig("$(path)/../figs/34mCampbellWGX.pdf",transparent = true)


# Frequency differences
FreqSolidWorks = [12.548,12.549,12.896,21.848,21.854]
for ifreq = 1:length(FreqSolidWorks)
    println("Frequency $ifreq Solidworks $(FreqSolidWorks[ifreq]) OWENS $(freq2[2,ifreq]) Percent Diff $((FreqSolidWorks[ifreq]-freq2[2,ifreq])/FreqSolidWorks[ifreq]*100)")
end
