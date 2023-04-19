import PyPlot
PyPlot.close("all")
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

import Statistics:mean
import DelimitedFiles
import Dierckx
import QuadGK
import FLOWMath
import HDF5

import ModelGen
import GyricFEA
import OWENS
import VAWTAero
import Composites

path,_ = splitdir(@__FILE__)

# include("$(path)/../../../../OWENS.jl/src/OWENS.jl")
# include("$(path)/../../../../ModelGen.jl/src/ModelGen.jl")
println("Set up Macro Geometry/Inputs")
rho = 1.225
Nslices = 30
ntheta = 30
RPM = 7.2
Vinf = 22.26786 * 0.3048
slc1 = 3
slc2 = 5
eta = 0.5
B = Nbld = 2
R = 177.2022*0.3048 #m
H = 1.02*R*2 #m
chord = 5.0*ones(Nslices)
omega = RPM / 60 * 2 * pi
tsr = omega*R/Vinf

shapeY = collect(LinRange(0,H,Nslices+1))
shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)
shapeX_spline = FLOWMath.Akima(shapeY, shapeX)
h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
h_elem = (shapeY[2:end] - shapeY[1:end-1])
h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
RefArea = RefArea_half*2
delta_xs = shapeX[2:end] - shapeX[1:end-1]
delta_zs = shapeY[2:end] - shapeY[1:end-1]

delta3D = atan.(delta_xs./delta_zs)

#########################################
### Set up aero forces
#########################################
println("Initialize Aerodynamics")
VAWTAero.setupTurb(shapeX,shapeY,B,chord,tsr,Vinf;AModel="DMS",DSModel="BV",
afname = "$(path)/../airfoils/NACA_0021.dat",
ifw=true,
ifw_libfile = joinpath(dirname(@__FILE__), "../openfast/build/modules/inflowwind/libifw_c_binding"),
turbsim_filename="$(path)/data/300mx300m12msETM_Coarse.bts",
ntheta,Nslices,rho,eta,RPI=true)

#########################################
### Set up mesh
#########################################
println("Create Mesh")
mymesh,myort,myjoint = ModelGen.create_mesh_struts(;Ht=15.0,
Hb = H, #blade height
R, # m bade radius
nblade = 2,
ntelem = 30, #tower elements
nbelem = 60, #blade elements
nselem = 10,
strut_mountpointbot = 0.1,
strut_mountpointtop = 0.1,
bshapex = bshapex=shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
bshapez = shapeY,
angularOffset = -pi/2) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

#########################################
### Set up Sectional Properties
#########################################
println("Calculate/Set up sectional properties")
#Tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn_twr = ModelGen.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

#Add the full path
for (i,airfoil) in enumerate(numadIn_twr.airfoil)
    numadIn_twr.airfoil[i] = "$path/../airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops_twr = ModelGen.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = ModelGen.getPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
nTwrElem = Int(mymesh.meshSeg[1])+1
sectionPropsArray_twr = ModelGen.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
stiff_twr, mass_twr = ModelGen.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

#Blades
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = ModelGen.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/../airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops_bld = ModelGen.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# Get blade spanwise position
bld1start = Int(mymesh.structuralNodeNumbers[1,1])
bld1end = Int(mymesh.structuralNodeNumbers[1,end])
spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = ModelGen.getPreCompOutput(numadIn_bld;plyprops = plyprops_bld)
sectionPropsArray_bld = ModelGen.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
stiff_bld, mass_bld = ModelGen.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = sum(Int,mymesh.meshSeg[Nbld+1+1:end]) #strut elements remain
# Nremain = 8 #strut elements remain
bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
sectionPropsArray = [sectionPropsArray_twr;bldssecprops; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = ModelGen.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld;VTKmeshfilename="$path/vtk/SNL5MW")

#########################################
### Create Aero Functions
#########################################

# # Example aero forces file input
# d1 = "$path/data/Simplified5MW2.geom"
# d2 = "$path/data/Simplified5MW2_ElementData.csv"
# d3 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.bld"
# d4 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.el"
# d5 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.ort"
# d6 = "$path/data/oldfiles/_15mTower_transient_dvawt_c_2_lcdt.mesh"
#
# aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = OWENS.mapCactusLoadsFile(d1,d2,d3,d4,d5,d6)
# aeroForcesCACTUS(t,azi) = OWENS.externalForcing(t,aerotimeArray,aeroForceValHist,aeroForceDof)

# aeroForcesDMSFile(t) = VAWTAero.mapCACTUSFILE_minimalio(t,mymesh,myel,"$path/data/DMS_Simplified5MW_ElementData.csv")
# Example realtime aero calculation setup
aeroForcesDMS(t,azi) = ModelGen.mapACDMS(t,azi,mymesh,myel,VAWTAero.AdvanceTurbineInterpolate;alwaysrecalc=true)

######################################
#### Perform Aerostructural One Way Test
#######################################

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

model = OWENS.Model(;analysisType = "ROM",
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS = 40.0,
delta_t = 0.05,
aeroLoadsOn = 2)

feamodel = GyricFEA.FEAModel(;analysisType = "ROM",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
numNodes = mymesh.numNodes,
RayleighAlpha = 0.0,
RayleighBeta = 0.0,
iterationType = "DI")

# Choose which aeroforces
# aeroForces = aeroForcesCACTUS
aeroForces = aeroForcesDMS

# deformaero3(x;newOmega=-1,newVinf=-1) = 0

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,epsilon_x_hist,kappa_y_hist,
kappa_z_hist,epsilon_z_hist,kappa_x_hist,epsilon_y_hist = OWENS.Unsteady(model,feamodel,mymesh,myel,aeroForces,VAWTAero.deformTurb)

println("Saving VTK time domain files")
ModelGen.gyricFEA_VTK("$path/vtk/SNL5MW_timedomain",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

##########################################
#### Get strain values at the blades #####
##########################################

meanepsilon_z_hist = mean(epsilon_z_hist,dims=1)
meanepsilon_y_hist = mean(epsilon_y_hist,dims=1)

N_ts = length(epsilon_x_hist[1,1,:])
eps_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
eps_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
eps_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
kappa_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
kappa_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
kappa_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
mesh_span_bld = zeros(length(mymesh.structuralNodeNumbers[1,:]))
composites_span_bld = zeros(length(bld_precompinput))
for ibld = 1:Nbld
    start = Int(mymesh.structuralNodeNumbers[ibld,1])
    stop = Int(mymesh.structuralNodeNumbers[ibld,end])
    x = mymesh.z[start:stop]
    x = x.-x[1] #zero
    x = x./x[end] #normalize
    mesh_span_bld[:] = mymesh.z[start:stop].-mymesh.z[start]
    global composites_span_bld = FLOWMath.akima(LinRange(0,1,length(mesh_span_bld)),mesh_span_bld,LinRange(0,1,length(bld_precompinput)))
    for its = 1:N_ts
        # Interpolate to the composite inputs
        eps_x_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,epsilon_x_hist[1,start:stop,its],composites_span_bld)
        eps_z_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,meanepsilon_z_hist[1,start:stop,its],composites_span_bld)
        eps_y_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,meanepsilon_y_hist[1,start:stop,its],composites_span_bld)
        kappa_x_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_x_hist[1,start:stop,its],composites_span_bld)
        kappa_y_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_y_hist[1,start:stop,its],composites_span_bld)
        kappa_z_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_z_hist[1,start:stop,its],composites_span_bld)
    end
end


##########################################
#### Get strain values at the tower #####
##########################################

N_ts = length(epsilon_x_hist[1,1,:])
eps_x_twr = zeros(1,N_ts,length(twr_precompinput))
eps_z_twr = zeros(1,N_ts,length(twr_precompinput))
eps_y_twr = zeros(1,N_ts,length(twr_precompinput))
kappa_x_twr = zeros(1,N_ts,length(twr_precompinput))
kappa_y_twr = zeros(1,N_ts,length(twr_precompinput))
kappa_z_twr = zeros(1,N_ts,length(twr_precompinput))

start = 1
stop = length(mymesh.type[mymesh.type.==1])
x = mymesh.z[start:stop]
x = x.-x[1] #zero
x = x./x[end] #normalize
mesh_span_twr = mymesh.z[start:stop].-mymesh.z[start]
composites_span_twr = FLOWMath.akima(LinRange(0,1,length(mesh_span_twr)),mesh_span_twr,LinRange(0,1,length(twr_precompinput)))
for its = 1:N_ts
    # Interpolate to the composite inputs
    eps_x_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,epsilon_x_hist[1,start:stop,its],composites_span_twr)
    eps_z_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,meanepsilon_z_hist[1,start:stop,its],composites_span_twr)
    eps_y_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,meanepsilon_y_hist[1,start:stop,its],composites_span_twr)
    kappa_x_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_x_hist[1,start:stop,its],composites_span_twr)
    kappa_y_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_y_hist[1,start:stop,its],composites_span_twr)
    kappa_z_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_z_hist[1,start:stop,its],composites_span_twr)

end


##########################################
#### Composite Failure & Buckling ########
##########################################

stress_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]),3)
SF_ult_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]))
SF_buck_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]))

ModelGen.calcSF(stress_U,SF_ult_U,SF_buck_U,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_U_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "tsaiwu",upper=true)

println("\nUPPER SURFACE")
if !isempty(SF_buck_U[SF_buck_U.>0.0])
    SF_buck_U[SF_buck_U.<0.0] .= 1e6
    SF_buck_U[:,:,1] .= 1e6 #ignore leading edge
    SF_buck_U[:,:,6] .= 1e6 #ignore trailing edge
    minbuck_sf,minbuck_sfidx = findmin(SF_buck_U)
    println("Worst buckling safety factor $(minbuck_sf)")
    println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_bld)) at lam $(minbuck_sfidx[3]) of $(length(lam_U_bld[minbuck_sfidx[2],:]))")
else
    println("Buckling not a factor, no sections in compression")
end
println("Minimum Safety Factor on Blade Surface: $(minimum(SF_ult_U))")
mymin,idx = findmin(SF_ult_U)
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld)) at lam $(idx[3]) of $(length(lam_U_bld[idx[2],:]))")
mymin,idx = findmin(SF_ult_U[:,:,3])
println("Spar Cap SF min: $mymin")
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
mymin,idx = findmin(SF_ult_U[:,:,1])
println("Leading Edge SF min: $mymin")
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")

println("Upper Spar")
for SF in SF_ult_U[idx[1],:,3]
    println(SF)
end

println("Upper Leading Edge")
for SF in SF_ult_U[idx[1],:,1]
    println(SF)
end

println("Upper Trailing Edge")
for SF in SF_ult_U[idx[1],:,6]
    println(SF)
end

println("Upper Fore Panel")
for SF in SF_ult_U[idx[1],:,2]
    println(SF)
end

println("Upper Aft Panel")
for SF in SF_ult_U[idx[1],:,5]
    println(SF)
end

println("Upper Buckling")
for istation = 1:length(composites_span_bld)
    println(minimum(SF_buck_U[idx[1],istation,:]))
end


stress_L = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]),3)
SF_ult_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))
SF_buck_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))

ModelGen.calcSF(stress_L,SF_ult_L,SF_buck_L,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_L_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "tsaiwu",upper=false)

println("\nLOWER SURFACE")
if !isempty(SF_buck_L[SF_buck_L.>0.0])
    SF_buck_L[SF_buck_L.<0.0] .= 1e6
    SF_buck_L[:,:,1] .= 1e6 #ignore leading edge
    SF_buck_L[:,:,6] .= 1e6 #ignore trailing edge
    minbuck_sf,minbuck_sfidx = findmin(SF_buck_L)
    println("Worst buckling safety factor $(minbuck_sf)")
    println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_bld)) at lam $(minbuck_sfidx[3]) of $(length(lam_L_bld[minbuck_sfidx[2],:]))")
else
    println("Buckling not a factor, no sections in compression")
end
println("Minimum Safety Factor on Blade Surface: $(minimum(SF_ult_L))")
mymin,idx = findmin(SF_ult_L)
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld)) at lam $(idx[3]) of $(length(lam_L_bld[idx[2],:]))")
mymin,idx = findmin(SF_ult_L[:,:,3])
println("Spar Cap SF min: $mymin")
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
mymin,idx = findmin(SF_ult_L[:,:,1])
println("Leading Edge SF min: $mymin")
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")

println("Lower Spar")
for SF in SF_ult_L[idx[1],:,3]
    println(SF)
end

println("Lower Leading Edge")
for SF in SF_ult_L[idx[1],:,1]
    println(SF)
end

println("Lower Trailing Edge")
for SF in SF_ult_L[idx[1],:,6]
    println(SF)
end

println("Lower Fore Panel")
for SF in SF_ult_L[idx[1],:,2]
    println(SF)
end

println("Lower Aft Panel")
for SF in SF_ult_L[idx[1],:,5]
    println(SF)
end

println("Lower Buckling")
for istation = 1:length(composites_span_bld)
    println(minimum(SF_buck_L[idx[1],istation,:]))
end

stress_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
SF_ult_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
SF_buck_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

ModelGen.calcSF(stress_TU,SF_ult_TU,SF_buck_TU,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "tsaiwu",upper=true)


println("\nUPPER TOWER")
if !isempty(SF_buck_TU[SF_buck_TU.>0.0])
    SF_buck_TU[SF_buck_TU.<0.0] .= 1e6
    # SF_buck_TU[:,:,1] .= 1e6 #ignore leading edge
    # SF_buck_TU[:,:,6] .= 1e6 #ignore trailing edge
    minbuck_sf,minbuck_sfidx = findmin(SF_buck_TU)
    println("Worst buckling safety factor $(minbuck_sf)")
    println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_U_twr[minbuck_sfidx[2],:]))")
else
    println("Buckling not a factor, no sections in compression")
end
println("Minimum Safety Factor on tower Surface: $(minimum(SF_ult_TU))")
mymin,idx = findmin(SF_ult_TU)
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_U_twr[idx[2],:]))")

println("Leading Edge")
for SF in SF_ult_TU[idx[1],:,1]
    println(SF)
end

println("Buckling")
for istation = 1:length(composites_span_twr)
    println(minimum(SF_buck_TU[minbuck_sfidx[1],istation,:]))
end


stress_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
SF_ult_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
SF_buck_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

ModelGen.calcSF(stress_TL,SF_ult_TL,SF_buck_TL,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "tsaiwu",upper=false)


println("\nLOWER TOWER")
if !isempty(SF_buck_TL[SF_buck_TL.>0.0])
    SF_buck_TL[SF_buck_TL.<0.0] .= 1e6
    # SF_buck_TL[:,:,1] .= 1e6 #ignore leading edge
    # SF_buck_TL[:,:,6] .= 1e6 #ignore trailing edge
    minbuck_sf,minbuck_sfidx = findmin(SF_buck_TL)
    println("Worst buckling safety factor $(minbuck_sf)")
    println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_L_twr[minbuck_sfidx[2],:]))")
else
    println("Buckling not a factor, no sections in compression")
end
println("Minimum Safety Factor on tower Surface: $(minimum(SF_ult_TL))")
mymin,idx = findmin(SF_ult_TL)
println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_L_twr[idx[2],:]))")

println("Leading Edge")
for SF in SF_ult_TL[idx[1],:,1]
    println(SF)
end

println("Buckling")
for istation = 1:length(composites_span_twr)
    println(minimum(SF_buck_TL[minbuck_sfidx[1],istation,:]))
end

function massOWENS(sectionPropsArray,myort)
    mass = 0.0
    for (i,sectionProp) in enumerate(sectionPropsArray)
        lenEl = 0
        try
            lenEl = myort.Length[i]
        catch
            lenEl = myort.Length[i-1]
        end
        rhoA = sectionProp.rhoA[1]
        mass += lenEl*rhoA
    end
return mass
end

massOwens = massOWENS(sectionPropsArray,myort)
println("Mass Turbine: $massOwens")

######################################
#### Plot
#######################################

# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,1])
# PyPlot.ylabel("FReaction Hist 1")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction1.pdf",transparent = true)
#
# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,2])
# PyPlot.ylabel("FReaction Hist 2")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction2.pdf",transparent = true)
#
# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,3])
# PyPlot.ylabel("FReaction Hist 3")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction3.pdf",transparent = true)
#
# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,4])
# PyPlot.ylabel("FReaction Hist 4")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction4.pdf",transparent = true)
#
# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,5])
# PyPlot.ylabel("FReaction Hist 5")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction5.pdf",transparent = true)
#
# PyPlot.figure()
# PyPlot.plot(t,FReactionHist[:,6])
# PyPlot.ylabel("FReaction Hist 6")
# # PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_FReaction6.pdf",transparent = true)
#
# PyPlot.figure()
# for ii = 1:length(uHist[1,:])
# PyPlot.plot(1:length(old_uHist[ii,:]),old_uHist[ii,:],"k-")
#     PyPlot.plot(1:length(uHist[ii,:]),uHist[ii,:],"k--")
#     if ii%10 == 0.0
#         PyPlot.savefig("$(path)/figs/AeroOnly/SNL5MW_Uhist$ii.pdf",transparent = true)
#         PyPlot.figure()
#     end
# end
