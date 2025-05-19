
import OWENS
import OWENSFEA
import OWENSAero
import FLOWMath
import DelimitedFiles
using Statistics:mean
using Test
import HDF5
import YAML

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

analysisType = "unsteady"
turbineType = "H-VAWT"
eta = 0.5
Nbld = 3
towerHeight = 0.5
rho = 0.94
Vinf = 10.1
controlStrategy = "constantRPM"
RPM = 8.0
Nslices = 35
ntheta = 30
structuralModel = "GX"
ntelem = 10
nbelem = 60
ncelem = 10
nselem = 5
ifw = false
AeroModel = "AD"
windINPfilename = "$path/300mx300m12msETM_Coarse.inp"
run(`$(OWENS.OWENSOpenFASTWrappers.turbsim()) $windINPfilename`)

ifw_libfile = nothing#"$path/../../openfast/build/modules/inflowwind/libifw_c_binding"
Blade_Height = 41.9
Blade_Radius = 6.0
numTS = 10
delta_t = 0.05
NuMad_geom_xlscsv_file_twr = "$path/TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$path/TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$path/GeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$path/Materials34m.csv"
NuMad_geom_xlscsv_file_strut = "$path/GeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$path/Materials34m.csv"
adi_lib = nothing#"$path/../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" 
adi_rootname = "$path/helical"

##############################################
# Setup
#############################################

tocp = [0.0;10.0; 1e6]
Omegaocp = [34.0; 34.0; 34.0]./60 #control inputs

tocp_Vinf = [0.0;10.0; 1e6]
Vinfocp = [10.0;10.0;10.0]

shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
helix_angle = -pi/4
shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius#SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ) ones(length(shapeZ)).*Blade_Radius#
shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius # zeros(length(shapeX))#

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng  = OWENS.setupOWENS(OWENSAero,path;
    rho, #windio
    Nslices, #modeling options
    ntheta, #modeling options
    eta, #windio
    B = Nbld, #windio
    H = Blade_Height, #windio
    R = Blade_Radius, #windio
    shapeZ, #windio
    shapeX, #windio
    shapeY,
    ifw, #modeling options
    delta_t, #modeling options
    numTS, #modeling options
    adi_lib, #remove - make smart enough to find
    adi_rootname, #modeling options
    AD15hubR = 0.0, #modeling options
    windINPfilename, #modeling options
    ifw_libfile, #remove - make smart enough to find
    NuMad_geom_xlscsv_file_twr,#windio
    NuMad_mat_xlscsv_file_twr,#windio
    NuMad_geom_xlscsv_file_bld,#windio
    NuMad_mat_xlscsv_file_bld,#windio
    NuMad_geom_xlscsv_file_strut,#windio
    NuMad_mat_xlscsv_file_strut,#windio
    Htwr_base=towerHeight,
    ntelem, 
    nbelem, 
    ncelem,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.1,0.5,0.9],
    strut_bld_mountpoint = [0.05,0.5,0.95],
    AeroModel, #AD, DMS, AC
    DynamicStallModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0,#pi/2,
    meshtype = turbineType)


# PyPlot.figure()
# PyPlot.plot(mymesh.x,mymesh.z,"b-")
#     for myi = 1:length(mymesh.x)
#         PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
#         PyPlot.draw()
#         #sleep(0.1)
#     end
# PyPlot.xlabel("x")
# PyPlot.ylabel("y")

PyPlot.figure()
PyPlot.scatter3D(mymesh.x,mymesh.y,mymesh.z)
PyPlot.scatter3D(shapeY,shapeX,shapeZ,color="red")
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")

top_idx = Int(myjoint[6,1])
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

if AeroModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
    tocp,
    Omegaocp,
    tocp_Vinf,
    Vinfocp,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = 2,
    turbineStartup = 1,
    generatorOn = true,
    useGeneratorFunction = true,
    driveTrainOn = true,
    JgearBox = 250.0,#(2.15e3+25.7)/12*1.35582*100,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    driveShaftProps = OWENS.DriveShaftProps(10000,1.5e2), #8.636e5*1.35582*0.6
    OmegaInit = Omegaocp[1])

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn=false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")



t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist = OWENS.Unsteady_Land(inputs;
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero,system,assembly)

################################################################
################ SAVE VTK TIME DOMAIN OUTPUT ###################
################################################################


println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","e_x","e_y","e_z","k_x","k_y","k_z","Fx_Reaction","Fy_Reaction","Fz_Reaction","Mx_Reaction","My_Reaction","Mz_Reaction"]#,"Fx","Fy","Fz","Mx","My","Mz"]
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


epsilon_x_histused = mean(epsilon_x_hist;dims=1)
epsilon_y_histused = mean(epsilon_y_hist;dims=1)
epsilon_z_histused = mean(epsilon_z_hist;dims=1)
kappa_x_histused = mean(kappa_x_hist;dims=1)
kappa_y_histused = mean(kappa_y_hist;dims=1)
kappa_z_histused = mean(kappa_z_hist;dims=1)

# fill in the big matrix
for it = 1:length(t)

    userPointData[1,it,:] = EA_points
    userPointData[2,it,:] = EIyy_points
    userPointData[3,it,:] = EIzz_points
    for iel = 1:length(myel.props)
        nodes = mymesh.conn[iel,:]
        userPointData[4,it,Int.(nodes)] .= epsilon_x_histused[1,iel,it] 
        userPointData[5,it,Int.(nodes)] .= epsilon_y_histused[1,iel,it] 
        userPointData[6,it,Int.(nodes)] .= epsilon_z_histused[1,iel,it] 
        userPointData[7,it,Int.(nodes)] .= kappa_x_histused[1,iel,it] 
        userPointData[8,it,Int.(nodes)] .= kappa_y_histused[1,iel,it] 
        userPointData[9,it,Int.(nodes)] .= kappa_z_histused[1,iel,it] 
    end
    userPointData[10,it,:] .= FReactionHist[it,1:6:end]
    userPointData[11,it,:] .= FReactionHist[it,2:6:end]
    userPointData[12,it,:] .= FReactionHist[it,3:6:end]
    userPointData[13,it,:] .= FReactionHist[it,4:6:end]
    userPointData[14,it,:] .= FReactionHist[it,5:6:end]
    userPointData[15,it,:] .= FReactionHist[it,6:6:end]
    
    # userPointData[4,it,:] = FReactionHist[it,1:6:end]
    # userPointData[5,it,:] = FReactionHist[it,2:6:end]
    # userPointData[6,it,:] = FReactionHist[it,3:6:end]
    # userPointData[7,it,:] = FReactionHist[it,4:6:end]
    # userPointData[8,it,:] = FReactionHist[it,5:6:end]
    # userPointData[9,it,:] = FReactionHist[it,6:6:end]
end

azi=aziHist#./aziHist*1e-6
VTKsaveName = "$path/vtk/helical"
OWENS.OWENSFEA_VTK(VTKsaveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)

# Open Paraview, open animation pane, adjust as desired, export animation (which exports frames)
# ffmpeg -i Ux.%04d.png -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 -y -an -pix_fmt yuv420p video34m34RPM_Ux.mp4

#########################################
### Ultimate Failure #####
#########################################

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
topstrainout_tower_U,topstrainout_tower_L = OWENS.extractSF(bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
composite_station_idx_U_strut = [1,6,3,2,5],
composite_station_name_U_strut = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_L_strut = [1,6,3,2,5],
composite_station_name_L_strut = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_U_bld = [1,6,3,2,5],
composite_station_name_U_bld = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_L_bld = [1,6,3,2,5],
composite_station_name_L_bld = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
Twr_LE_U_idx=1,Twr_LE_L_idx=1,AD15bldNdIdxRng,AD15bldElIdxRng) #TODO: add in ability to have material safety factors and load safety factors

nothing


# ffmpeg -i smallerhelical.%04d.png -vf palettegen=reserve_transparent=1 palette.png
# ffmpeg -framerate 30 -i smallerhelical.%04d.png -i palette.png -lavfi paletteuse=alpha_threshold=30 -gifflags -offsetting smallerhelical.gif
