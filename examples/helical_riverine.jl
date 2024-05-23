
import OWENS
import OWENSFEA
import OWENSAero
import FLOWMath
import DelimitedFiles
using Statistics:mean
using Test
import HDF5
import YAML
using StructTypes
import OrderedCollections

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

Inp = OWENS.MasterInput("$path/RM2_OWENS_modeling_options.yaml")

windio = YAML.load_file("$path/RM1.yaml"; dicttype=OrderedCollections.OrderedDict{Symbol,Any})


# mutable struct myTest
#     can_log
#     sys_log
#     version
#     info
#     myfloat
#     myTest() = new()
# end

# yamlInputtest = YAML.load_file("$path/testdata.yaml"; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
# test_dict = yamlInputtest[:tests][1]
# println(test_dict)

# StructTypes.StructType(::Type{myTest}) = StructTypes.Mutable()
# test = StructTypes.constructfrom(myTest, test_dict)



# Unpack inputs, or you could directly input them here and bypass the file 

verbosity = 1

analysisType = Inp.analysisType
turbineType = "H-VAWT"#Inp.turbineType
eta = Inp.eta
Nbld = 3#Inp.Nbld
towerHeight = Inp.towerHeight
rho = Inp.rho
Vinf = Inp.Vinf
controlStrategy = Inp.controlStrategy
RPM = 8.0#Inp.RPM
Nslices = Inp.Nslices
ntheta = Inp.ntheta
structuralModel = Inp.structuralModel
ntelem = Inp.ntelem
nbelem = Inp.nbelem
ncelem = Inp.ncelem
nselem = Inp.nselem
ifw = Inp.ifw
AModel = "AD"#Inp.AModel
windINPfilename = Inp.windINPfilename
ifw_libfile = Inp.ifw_libfile
Blade_Height = Inp.Blade_Height
Blade_Radius = 6.0#Inp.Blade_Radius
numTS = 1000#Inp.numTS
delta_t = Inp.delta_t
NuMad_geom_xlscsv_file_twr = Inp.NuMad_geom_xlscsv_file_twr
NuMad_mat_xlscsv_file_twr = Inp.NuMad_mat_xlscsv_file_twr
NuMad_geom_xlscsv_file_bld = Inp.NuMad_geom_xlscsv_file_bld
NuMad_mat_xlscsv_file_bld = Inp.NuMad_mat_xlscsv_file_bld
NuMad_geom_xlscsv_file_strut = Inp.NuMad_geom_xlscsv_file_strut
NuMad_mat_xlscsv_file_strut = Inp.NuMad_mat_xlscsv_file_strut
adi_lib = Inp.adi_lib
adi_rootname = "$path/helical"#Inp.adi_rootname

##############################################
# Setup
#############################################

tocp = [0.0;10.0; 1e6]
Omegaocp = [34.0; 34.0; 34.0]./60 #control inputs

tocp_Vinf = [0.0;10.0; 1e6]
Vinfocp = [10.0;10.0;10.0]

# windio
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
SNL34Z = SNL34z.*Blade_Height #windio
SNL34X = SNL34x.*Blade_Radius #windio

shapeY = SNL34Z#collect(LinRange(0,H,Nslices+1))
helix_angle = pi/2
shapeX = sin.(shapeY/maximum(shapeY)*helix_angle).*Blade_Radius#SNL34X#R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY) ones(length(shapeY)).*Blade_Radius#
bshapey =cos.(shapeY/maximum(shapeY)*helix_angle).*Blade_Radius # zeros(length(shapeX))#

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng  = OWENS.setupOWENS(OWENSAero,path;
    rho, #windio
    Nslices, #modeling options
    ntheta, #modeling options
    RPM, #remove
    Vinf, #remove
    eta, #windio
    B = Nbld, #windio
    H = Blade_Height, #windio
    R = Blade_Radius, #windio
    shapeY, #windio
    shapeX, #windio
    bshapey,
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
    Ht=towerHeight,
    ntelem, 
    nbelem, 
    ncelem,
    nselem,
    joint_type = 0,
    strut_mountpointbot = 0.03,
    strut_mountpointtop = 0.03,
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = pi/2,
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
PyPlot.scatter3D(shapeX,bshapey,shapeY,color="red")

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

if AModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;analysisType = structuralModel,
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
saveName = "$path/vtk/RM2"
OWENS.OWENSFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)

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
LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
Twr_LE_U_idx=1,Twr_LE_L_idx=1,AD15bldNdIdxRng,AD15bldElIdxRng) #TODO: add in ability to have material safety factors and load safety factors

nothing


# ffmpeg -i smallerhelical.%04d.png -vf palettegen=reserve_transparent=1 palette.png
# ffmpeg -framerate 30 -i smallerhelical.%04d.png -i palette.png -lavfi paletteuse=alpha_threshold=30 -gifflags -offsetting smallerhelical.gif
