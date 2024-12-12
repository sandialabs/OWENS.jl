
import OWENS
import OWENSAero
import DelimitedFiles
using Statistics:mean
using Test
import FFTW

# import PyPlot
# PyPlot.close("all")
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

# function runprofilefunction()
path = runpath = splitdir(@__FILE__)[1]

# Unpack inputs, or you could directly input them here and bypass the file 

verbosity = 1


analysisType = "unsteady"
turbineType = "H-VAWT"
eta = 0.5
Nbld = 2
towerHeight = 0.5
Vinf = 1.2
controlStrategy = "constantRPM"
Nslices = 20
ntheta = 30
structuralModel = "TNB"
ntelem = 20
nbelem = 60
ncelem = 10
nselem = 10
ifw = false
AeroModel = "DMS"
windINPfilename = "$path/300mx300m12msETM_Coarse.bts"
ifw_libfile = nothing#"$path/../../openfast/build/modules/inflowwind/libifw_c_binding"
Blade_Height = 20.0
Blade_Radius = 4.5
area = Blade_Height*2*Blade_Radius
numTS = 200
delta_t = 0.005
NuMad_geom_xlscsv_file_twr = "$path/TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$path/TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$path/GeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$path/Materials34m.csv"
NuMad_geom_xlscsv_file_strut = "$path/GeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$path/Materials34m.csv"
adi_lib = nothing#"$path/../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" 
adi_rootname = "$path/helical"

# TSR = 3.0

# omega = Vinf/Blade_Radius*TSR
RPM = 1e-6#omega / (2*pi) * 60

fluid_density = 998.0
fluid_dyn_viscosity = 1.792E-3
number_of_blades = Nbld
WindType = 3

Aero_AddedMass_Active = false # if on aero side
AddedMass_Coeff_Ca=1.0 #For either side, set to 0 or don't include for it to be off on structures
Aero_RotAccel_Active = false
Aero_Buoyancy_Active = false

##############################################
# Setup
#############################################

tocp = [0.0;10.0; 1e6]
Omegaocp = [RPM; RPM; RPM]./60 #control inputs

tocp_Vinf = [0.0;10.0; 1e6]
Vinfocp = [Vinf;Vinf;Vinf]

shapeZ = collect(LinRange(0,Blade_Height,Nslices+1))
helix_angle = 0.0#-pi/4
shapeX = cos.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius#SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ) ones(length(shapeZ)).*Blade_Radius#
shapeY = sin.(shapeZ/maximum(shapeZ)*helix_angle).*Blade_Radius # zeros(length(shapeX))#

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAeroreal,
mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng,AD15bldElIdxRng,custom_mesh_outputs,stiff_array,mass_array = OWENS.setupOWENS(OWENSAero,path;
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
    Aero_AddedMass_Active,
    AddedMass_Coeff_Ca,
    Aero_RotAccel_Active,
    Aero_Buoyancy_Active,
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)



# PyPlot.figure()
# for icon = 1:length(mymesh.conn[:,1])
#     idx1 = mymesh.conn[icon,1]
#     idx2 = mymesh.conn[icon,2]
#     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
#     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
#     # sleep(0.1)
# end

# for ijoint = 1:length(myjoint[:,1])
#     idx2 = Int(myjoint[ijoint,2])
#     idx1 = Int(myjoint[ijoint,3])
#     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
#     sleep(0.1)
# end

# PyPlot.scatter3D(1,1,1,"b.")

nothing

# Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
# this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
# includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
# to work with OWENS inputs as well.

nothing

# If the sectional properties material files includes cost information, that is combined with the density 
# to estimate the overall material cost of of materials in the blades

if verbosity>0
    
    println("\nBlades' Mass Breakout")
    for (i,name) in enumerate(plyprops_bld.names)
        println("$name $(mass_breakout_blds[i]) kg, $(plyprops_bld.costs[i]) \$/kg: \$$(mass_breakout_blds[i]*plyprops_bld.costs[i])")
    end
    
    println("\nTower Mass Breakout")
    for (i,name) in enumerate(plyprops_twr.names)
        println("$name $(mass_breakout_twr[i]) kg, $(plyprops_twr.costs[i]) \$/kg: \$$(mass_breakout_twr[i]*plyprops_twr.costs[i])")
    end
    
    println("Total Material Cost Blades: \$$(sum(mass_breakout_blds.*plyprops_bld.costs))")
    println("Total Material Cost Tower: \$$(sum(mass_breakout_twr.*plyprops_twr.costs))")
    println("Total Material Cost: \$$(sum(mass_breakout_blds.*plyprops_bld.costs)+ sum(mass_breakout_twr.*plyprops_twr.costs))")
    
end

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
24 1 0
24 2 0
24 3 0
24 4 0
24 5 0
24 6 0
87 1 0
87 2 0
87 3 0
87 4 0
87 5 0
87 6 0
146 1 0
146 2 0
146 3 0
146 4 0
146 5 0
146 6 0
]

nothing

# There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

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
aeroLoadsOn = 1)


inputs.iteration_parameters.MAXITER = 5

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.
forced_node = 85
function flappingForces(t,azi)
    if t<0.4
        Fexternal = [-5000000.0*t^2]
        Fdof = [(forced_node-1)*6+1]
    else#if t<0.5
        Fexternal = [0.0]
        Fdof = [(forced_node-1)*6+1]
    # else
    #     Fexternal,Fdof = aeroForces(t,azi)
    end
    return Fexternal, Fdof
end



function deformTurb(azi;newOmega=-1,newVinf=-1,bld_x=-1,
bld_z=-1,
bld_twist=-1,
accel_flap_in=-1,
accel_edge_in=-1,
gravity = [0.0,0.0,-9.81],
steady=false) 

end

deformAero = deformTurb
# deformAero = deformAeroreal

Keff = 6e5 #N/m
L = 9.5 #m 

EI = 2e8 #this is what was in the paper, which is not the same as the equation #Keff*L^3/3
# K = 3EI/L^3
# Set the element properties
for iprop = 1:length(myel.props)
    # do just the non-super stiff elements
    if myel.props[iprop].EIyy[2] < 1e9
        myel.props[iprop].EIyy .= EI
        myel.props[iprop].EIzz .= EI
        myel.props[iprop].GJ .= EI*5
        myel.props[iprop].EA .= EI*100
        stiff_array[iprop][1,1] = myel.props[iprop].EA[1]
        stiff_array[iprop][2,2] = myel.props[iprop].EA[1]/2.6*5/6
        stiff_array[iprop][3,3] = myel.props[iprop].EA[1]/2.6*5/6
        stiff_array[iprop][4,4] = myel.props[iprop].GJ[1]
        stiff_array[iprop][5,5] = myel.props[iprop].EIzz[1]
        stiff_array[iprop][6,6] = myel.props[iprop].EIyy[1]
    end
end

system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,myel.props,stiff_array,mass_array;damp_coef=0.0001)

"initCond`: array containing initial conditions.
    initCond(i,1) node number for init cond i.
    initCond(i,2) local DOF number for init cond i.
    initCond(i,3) value for init cond i.`"


ihalf = round(Int,size(mymesh.structuralNodeNumbers)[2]/2)
halfbladenodes = mymesh.structuralNodeNumbers[1,ihalf:end]
displace_x = LinRange(0,Blade_Height/2,length(halfbladenodes))
analyticalforce = 600000*0
displace_y = analyticalforce*displace_x.^2.0./(6*2e8).*(3*Blade_Height/2.0.-displace_x)
initTopConditions = zeros(length(displace_y)*2,3)
initTopConditions[1:length(displace_y),1] = halfbladenodes
initTopConditions[1:length(displace_y),2] .= 1
initTopConditions[1:length(displace_y),3] = displace_y

curve_y = atan.(displace_y,displace_x)
initTopConditions[length(displace_y)+1:end,1] = halfbladenodes
initTopConditions[length(displace_y)+1:end,2] .= 5
initTopConditions[length(displace_y)+1:end,3] = curve_y

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
dataOutputFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
initCond = initTopConditions,
gravityOn = [0,0,0.0],
numNodes = mymesh.numNodes,
numModes=200,
AddedMass_Coeff_Ca,
RayleighAlpha = 0.00,
RayleighBeta = 0.00,

iterationType = "DI")

nothing



println("Running Unsteady")
topdata = OWENS.Unsteady_Land(inputs;system,assembly,returnold=false,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=flappingForces,deformAero)
t = topdata.t
aziHist = topdata.aziHist
OmegaHist = topdata.OmegaHist
OmegaDotHist = topdata.OmegaDotHist
gbHist = topdata.gbHist
gbDotHist = topdata.gbDotHist
gbDotDotHist = topdata.gbDotDotHist
FReactionHist = topdata.FReactionHist
FTwrBsHist = topdata.FTwrBsHist
genTorque = topdata.genTorque
genPower = topdata.genPower
torqueDriveShaft = topdata.torqueDriveShaft
uHist = topdata.uHist
uHist_prp = topdata.uHist_prp
epsilon_x_hist = topdata.epsilon_x_hist
epsilon_y_hist = topdata.epsilon_y_hist
epsilon_z_hist = topdata.epsilon_z_hist
kappa_x_hist = topdata.kappa_x_hist
kappa_y_hist = topdata.kappa_y_hist
kappa_z_hist = topdata.kappa_z_hist
FPtfmHist = topdata.FPtfmHist
FHydroHist = topdata.FHydroHist
FMooringHist = topdata.FMooringHist
topFexternal_hist = topdata.topFexternal_hist
rbDataHist = topdata.rbDataHist
uddotHist = topdata.uddotHist
udotHist = topdata.udotHist

nothing

OWENS_tip_displ = uHist[:,(forced_node-1)*6+1]

if AddedMass_Coeff_Ca == 1.0
    omega_OF = 4.5105 * 2*pi
else
    omega_OF = 13.53* 2*pi
end

ofast_tdispl = cos.(t.*omega_OF)

ts_start = 111
# PyPlot.figure()
# PyPlot.plot(t,ofast_tdispl,".",label="OpenFAST")
# PyPlot.plot(t[ts_start:end].-0.32,OWENS_tip_displ[ts_start:end],".",label="OWENS")
# PyPlot.legend()
# PyPlot.savefig("$(path)AddedMassOff.pdf",transparent = true)

signal = OWENS_tip_displ[ts_start:end]
L = length(signal)
if L%2 != 0
    signal = signal[ts_start:end-1]
    L = length(signal)
end
# signal .-= mean(signal)
Y = FFTW.fft(signal)


Fs = 1/(t[2]-t[1])
P2 = abs.(Y./L)
P1 = P2[1:round(Int,L/2)+1]
P1[2:end-1] = 2*P1[2:end-1]

f = Fs.*(0:round(Int,L/2))./L

# PyPlot.figure()
# PyPlot.plot(f,P1,color=plot_cycle[1],label="OWENS")
# PyPlot.plot([omega_OF,omega_OF]./(2*pi),[0.0,1.0],color=plot_cycle[2],label="OpenFAST")
# # PyPlot.title("Single-Sided Amplitude Spectrum of X(t)")
# PyPlot.xlabel("Frequency (Hz)")
# PyPlot.ylabel("Frequency Amplitude Norm")
# PyPlot.legend()
# PyPlot.xlim([0.0,30.0])
# PyPlot.ylim([0.0,1.0])
# PyPlot.savefig("$(path)/added_Mass_off_bode.pdf",transparent = true)

val,idxmax = findmax(P1)
f_natural = f[idxmax]
percentdiff = (omega_OF/(2*pi)-f_natural)/(omega_OF/(2*pi))*100
println("percent_diff: $percentdiff%")
@test abs(percentdiff)<5.0

# import FLOWMath
# displace_spl = FLOWMath.Akima(t,OWENS_tip_displ)
# velocity = [FLOWMath.derivative(displace_spl,t[i]) for i = 1:length(t)]
# velocity_spl = FLOWMath.Akima(t,velocity)
# accel = [FLOWMath.derivative(velocity_spl,t[i]) for i = 1:length(t)]

# PyPlot.figure("Vel")
# PyPlot.plot(t,velocity)
# PyPlot.plot(t,udotHist[:,(forced_node-1)*6+1])

# PyPlot.figure("Accel")
# PyPlot.plot(t,accel)
# PyPlot.plot(t,uddotHist[:,(forced_node-1)*6+1])

# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

# azi=aziHist#./aziHist*1e-6
# VTKsaveName = "$path/vtk/flapping_added_mass"
# tsave_idx=1:1:numTS-1
# OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
#     epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
#     FReactionHist,topFexternal_hist;tsave_idx)

nothing