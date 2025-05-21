import OWENS
import OWENSAero
#### import PyPlot
runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
# runpath = path = splitdir(@__FILE__)[1]
Inp = OWENS.MasterInput("$runpath/sampleOWENS.yml")

nothing

verbosity = 1

analysisType = Inp.analysisType
turbineType = Inp.turbineType
eta = Inp.eta
Nbld = Inp.Nbld
towerHeight = Inp.towerHeight
rho = Inp.rho
Vinf = Inp.Vinf
controlStrategy = Inp.controlStrategy
RPM = Inp.RPM
Nslices = Inp.Nslices
ntheta = Inp.ntheta
structuralModel = Inp.structuralModel
ntelem = Inp.ntelem
nbelem = Inp.nbelem
ncelem = Inp.ncelem
nselem = Inp.nselem
ifw = Inp.ifw
WindType = Inp.WindType
AeroModel = Inp.AeroModel
windINPfilename = "$(path)$(Inp.windINPfilename)"
ifw_libfile = Inp.ifw_libfile
if ifw_libfile == "nothing"
    ifw_libfile = nothing
end
Blade_Height = Inp.Blade_Height
Blade_Radius = Inp.Blade_Radius
numTS = Inp.numTS
delta_t = Inp.delta_t
NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld)"
NuMad_mat_xlscsv_file_bld = "$(path)$(Inp.NuMad_mat_xlscsv_file_bld)"
NuMad_geom_xlscsv_file_strut = "$(path)$(Inp.NuMad_geom_xlscsv_file_strut)"
NuMad_mat_xlscsv_file_strut = "$(path)$(Inp.NuMad_mat_xlscsv_file_strut)"
adi_lib = Inp.adi_lib
if adi_lib == "nothing"
    adi_lib = nothing
end
adi_rootname = "$(path)$(Inp.adi_rootname)"

B = Nbld
R = Blade_Radius#177.2022*0.3048 #m
H = Blade_Height#1.02*R*2 #m

shapeZ = collect(LinRange(0,H,Nslices+1))
shapeX = R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

nothing

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    RPM,
    Vinf,
    eta,
    B,
    H,
    R,
    shapeZ,
    shapeX,
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
    strut_twr_mountpoint = [0.2,0.8],
    strut_bld_mountpoint = [0.2,0.8],
    ntelem,
    nbelem,
    ncelem,
    nselem,
    joint_type = 0,
    c_mount_ratio = 0.05,
    AeroModel, #AD, DMS, AC
    DynamicStallModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)

nothing

nothing

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

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

nothing

if AeroModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2)

nothing

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
dataOutputFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = true,
numNodes = mymesh.numNodes,
RayleighAlpha = 0.05,
RayleighBeta = 0.05,
iterationType = "DI")

nothing

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

if AD15On #TODO: move this into the run functions
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

nothing

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

nothing

##########################################
#### Ultimate Failure #####
##########################################

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
topstrainout_tower_U,topstrainout_tower_LtopDamage_blade_U,
topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
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
Twr_LE_U_idx=1,Twr_LE_L_idx=1,
AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
