import OWENS
import OWENSFEA
import OWENSAero
import FLOWMath
import DelimitedFiles
using Statistics:mean
using Statistics
using GXBeam
using Test

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

# function runprofilefunction()
path = runpath = splitdir(@__FILE__)[1]

nothing

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
AeroModel = "DMS"
windINPfilename = nothing
ifw_libfile = nothing
Blade_Height = 41.9
Blade_Radius = 17.1
numTS = 2000
delta_t = 0.05
NuMad_geom_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerGeom.csv"
NuMad_mat_xlscsv_file_twr = "$(path)/data/NuMAD_34m_TowerMaterials.csv"
NuMad_geom_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mGeomBlades.csv"
NuMad_mat_xlscsv_file_bld = "$(path)/data/NuMAD_SNL34mMaterials.csv"
NuMad_geom_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mGeomStruts.csv"
NuMad_mat_xlscsv_file_strut = "$(path)/data/NuMAD_SNL34mMaterials.csv"
adi_lib = nothing
adi_rootname = "$(path)/SNL34m"

##############################################
# Setup
#############################################

SNL34m_5_3_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Vinf.csv",',',skipstart = 0)
SNL34m_5_3_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_RPM.csv",',',skipstart = 0)
SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque.csv",',',skipstart = 0)


new_t = LinRange(SNL34m_5_3_RPM[1,1],SNL34m_5_3_RPM[end,1],100)
new_RPM = OWENS.safeakima(SNL34m_5_3_RPM[:,1],SNL34m_5_3_RPM[:,2],new_t)

Vinf_spec = OWENS.safeakima(SNL34m_5_3_Vinf[:,1],SNL34m_5_3_Vinf[:,2],new_t)

offsetTime = 20.0 # seconds
tocp = [0.0;new_t.+offsetTime; 1e6]
Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]]./60 .*0 .+33.92871/60
t_Vinf = [0;new_t;1e6]
Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
tocp_Vinf = [0.0;t_Vinf.+offsetTime; 1e6]
Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]].*1e-6

controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = OWENS.safeakima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
SNL34_unit_xz = [x_shape z_shape]
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
mass_breakout_blds,mass_breakout_twr,system, assembly, sections = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    eta,
    B = Nbld,
    H = Blade_Height,
    R = Blade_Radius,
    shapeZ,
    shapeX,
    shapeY=zero(shapeX),
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
    AeroModel, #AD, DMS, AC
    DynamicStallModel="BV",
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

# node, dof, bc
top_idx = Int(myjoint[4,2])
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
top_idx 1 0
top_idx 2 0]

# filename = "$(path)/data/newmesh_34m"
# OWENS.saveOWENSfiles(filename,mymesh,myort,myjoint,myel,pBC,numadIn_bld)


# if testModal
##############################################
# Modal Test
#############################################

FEAinputs = OWENSFEA.FEAModel(;analysisType = "M",
        # dataOutputFilename = "$path/data/outplat.out",
        joint = myjoint,
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        gravityOn=true,
        nlOn = true,
        tolerance = 1e-6,
        spinUpOn = true,
        iterationType = "NR",
        numNodes = mymesh.numNodes,
        numModes = 16)  # number of modes to calculate)


starttime = time()
freq = OWENS.AutoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    rotSpdArrayRPM = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
    )
freqOWENS = [freq[:,i] for i=1:4:FEAinputs.numModes-2]
elapsedtime = time() - starttime

starttime2 = time()
FEAinputs.analysisType = "GX"
freq2 = OWENS.AutoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    rotSpdArrayRPM = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
    VTKsavename="$path/campbellVTK/SNL34m",
    saveModes = [1,3,5], #must be int
    mode_scaling = 500.0,
    )
freqGX = [freq2[:,i] for i=1:2:FEAinputs.numModes-6-2]
elapsedtime2 = time() - starttime2

println("OWENS: $elapsedtime")
println("GX: $elapsedtime2")

for imode = 1:length(freqOWENS)
    for ifreq = 1:length(freqOWENS[imode])
        # println("imode $imode, ifreq $ifreq")
        atol = freqGX[imode][ifreq]*0.065
        @test isapprox(freqOWENS[imode][ifreq],freqGX[imode][ifreq];atol)
    end
end

# import PyPlot
# PyPlot.close("all")
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(4.5, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=4.0)
# PyPlot.rc("legend", frameon=true)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# NperRevLines = 8
# rotSpdArrayRPM = LinRange(0.0, 40.0, 9) # int
# PyPlot.figure()
# for i=1:1:FEAinputs.numModes-2
#        PyPlot.plot(rotSpdArrayRPM,freq[:,i],color=plot_cycle[1],"-") #plot mode i at various rotor speeds
# end

# # Plot 34m experimental data
# SNL34_flap = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_Flatwise.csv",',',skipstart = 0)
# SNL34_lead = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_LeadLag.csv",',',skipstart = 0)

# PyPlot.plot(SNL34_flap[:,1],SNL34_flap[:,2],"k.",label="Flapwise Gauges")
# PyPlot.plot(SNL34_lead[:,1],SNL34_lead[:,2],"kx",label="Lead-Lag Gauges")
# PyPlot.plot(0,0,"k--",label="Tower-Mode")

# SNL34_1F = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1F.csv",',',skipstart = 0)
# SNL34_1BE = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1BE.csv",',',skipstart = 0)
# SNL34_2FA = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_2FA.csv",',',skipstart = 0)
# SNL34_2FS = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_2FS.csv",',',skipstart = 0)
# SNL34_1TO = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1TO.csv",',',skipstart = 0)
# SNL34_3F = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_3F.csv",',',skipstart = 0)

# PyPlot.plot(SNL34_1F[:,1],SNL34_1F[:,2],"k-")
# PyPlot.plot(SNL34_1BE[:,1],SNL34_1BE[:,2],"k-")
# PyPlot.plot(SNL34_2FA[:,1],SNL34_2FA[:,2],"k-")
# PyPlot.plot(SNL34_2FS[:,1],SNL34_2FS[:,2],"k-")
# PyPlot.plot(SNL34_1TO[:,1],SNL34_1TO[:,2],"k--")
# PyPlot.plot(SNL34_3F[:,1],SNL34_3F[:,2],"k-")

# #plot per rev lines
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
# # PyPlot.ylim([0.0,0.8])
# # PyPlot.savefig("$(path)/../figs/34mCampbell.pdf",transparent = true)

# # Add to figure
# for i=1:2:FEAinputs.numModes-6-2
#        PyPlot.plot(rotSpdArrayRPM,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
# end
# PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
# PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
# PyPlot.ylim([0,6.01])
# # PyPlot.savefig("$(path)/../figs/34mCampbellWGX.pdf",transparent = true)
