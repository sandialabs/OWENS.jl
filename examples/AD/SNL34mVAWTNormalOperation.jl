using Revise
import QuadGK
import OWENS
import OWENSFEA
import OWENSAero
import FLOWMath
import DelimitedFiles
using Statistics: mean
using Statistics

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize = (4.5, 3))
PyPlot.rc("font", size = 10.0)
PyPlot.rc("lines", linewidth = 1.5)
PyPlot.rc("lines", markersize = 3.0)
PyPlot.rc("legend", frameon = false)
PyPlot.rc("axes.spines", right = false, top = false)
PyPlot.rc("figure.subplot", left = 0.18, bottom = 0.17, top = 0.9, right = 0.9)
PyPlot.rc("figure", max_open_warning = 500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle = ["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# function runprofilefunction()
# path = runpath = splitdir(@__FILE__)[1]
path = runpath = "/home/fredrik/dev/sandia-OWENS/OWENS/examples/AD"

Inp = OWENS.MasterInput("$path/SNL34m_Inputs.yml")
# Inp = OWENS.MasterInput("$path/SNL34m_InputsAeroDyn.yml")

nothing

# Unpack inputs, or you could directly input them here and bypass the file

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
AModel = Inp.AModel
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

##############################################
# Setup
#############################################

SNL34m_5_3_Vinf = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Vinf.csv", ',', skipstart = 0)
SNL34m_5_3_RPM = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_RPM.csv", ',', skipstart = 0)
SNL34m_5_3_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/5.3_Torque.csv", ',', skipstart = 0)


new_t = LinRange(SNL34m_5_3_RPM[1, 1], SNL34m_5_3_RPM[end, 1], 100)
new_RPM = OWENS.safeakima(SNL34m_5_3_RPM[:, 1], SNL34m_5_3_RPM[:, 2], new_t)

Vinf_spec = OWENS.safeakima(SNL34m_5_3_Vinf[:, 1], SNL34m_5_3_Vinf[:, 2], new_t)

offsetTime = 20.0 # seconds
tocp = [0.0;new_t .+ offsetTime; 1.0e6]
Omegaocp = [new_RPM[1]; new_RPM; new_RPM[end]] ./ 60 .* 0 .+ 33.92871 / 60
t_Vinf = [0;new_t;1.0e6]
Vinf_spec = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]]
tocp_Vinf = [0.0;t_Vinf .+ offsetTime; 1.0e6]
Vinfocp = [Vinf_spec[1];Vinf_spec;Vinf_spec[end]] .* 1.0e-6

controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0, 41.9, length(controlpts) + 2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0, 41.9, 60))
x_shape = OWENS.safeakima(z_shape1, x_shape1, z_shape) #[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
SNL34_unit_xz = [x_shape;;z_shape]
SNL34x = SNL34_unit_xz[:, 1] ./ maximum(SNL34_unit_xz[:, 1])
SNL34z = SNL34_unit_xz[:, 2] ./ maximum(SNL34_unit_xz[:, 2])
SNL34Z = SNL34z .* Blade_Height
SNL34X = SNL34x .* Blade_Radius

shapeZ = SNL34Z #collect(LinRange(0,H,Nslices+1))
shapeX = SNL34X #R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

shapeX_spline = FLOWMath.Akima(SNL34Z, SNL34X)
RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, Blade_Height, atol = 1.0e-10)
RefArea = RefArea_half * 2

mymesh, myel, myort, myjoint, sectionPropsArray, mass_twr, mass_bld,
    stiff_twr, stiff_bld, bld_precompinput,
    bld_precompoutput, plyprops_bld, numadIn_bld, lam_U_bld, lam_L_bld,
    twr_precompinput, twr_precompoutput, plyprops_twr, numadIn_twr, lam_U_twr, lam_L_twr, aeroForces, deformAero,
    mass_breakout_blds, mass_breakout_twr, system, assembly, sections, AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(
    OWENSAero, path;
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
    shapeY = zero(shapeX),
    ifw,
    delta_t,
    numTS,
    adi_lib,
    adi_rootname,
    AD15hubR = 0.0,
    windINPfilename,
    ifw_libfile,
    NuMad_geom_xlscsv_file_twr, # = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS_Cables.csv",
    NuMad_mat_xlscsv_file_twr, # = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv",
    NuMad_geom_xlscsv_file_bld, # = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS.csv",
    NuMad_mat_xlscsv_file_bld, # = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
    NuMad_geom_xlscsv_file_strut,
    NuMad_mat_xlscsv_file_strut,
    Htwr_base = towerHeight,
    ntelem,
    nbelem,
    ncelem,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.03, 0.97],
    strut_bld_mountpoint = [0.03, 0.97],
    AModel, #AD, DMS, AC
    DSModel = "BV",
    RPI = true,
    cables_connected_to_blade_base = true,
    angularOffset = pi / 2,
    meshtype = turbineType
)

# Insert Dual's

using ForwardDiff

seed = ForwardDiff.Dual(1.0, 1.0)

sectionPropsArray = [
    OWENSFEA.SectionPropsArray(
            x.ac, x.twist, x.rhoA, seed * x.EIyy, seed * x.EIzz, x.GJ, seed * x.EA, x.rhoIyy,
            x.rhoIzz, x.rhoJ, x.zcm, x.ycm, x.a, seed * x.EIyz, x.alpha1, x.alpha2, x.alpha3,
            x.alpha4, x.alpha5, x.alpha6, x.rhoIyz, x.b, x.a0, x.aeroCenterOffset, x.xaf, x.yaf,
            x.added_M22, x.added_M33
        ) for x in sectionPropsArray
]
myel.props = sectionPropsArray


error()

# PyPlot.figure()
# for icon in 1:length(mymesh.conn[:, 1])
#     idx1 = mymesh.conn[icon, 1]
#     idx2 = mymesh.conn[icon, 2]
#     PyPlot.plot3D([mymesh.x[idx1], mymesh.x[idx2]], [mymesh.y[idx1], mymesh.y[idx2]], [mymesh.z[idx1], mymesh.z[idx2]], "k.-")
#     PyPlot.plot3D([1, 1], [1, 1], [1, 1], "k.-")
#     PyPlot.text3D(mymesh.x[idx1] .+ rand() / 30, mymesh.y[idx1] .+ rand() / 30, mymesh.z[idx1] .+ rand() / 30, "$idx1", ha = "center", va = "center")
#     # sleep(0.1)
# end

# for ijoint in 1:length(myjoint[:, 1])
#     idx2 = Int(myjoint[ijoint, 2])
#     idx1 = Int(myjoint[ijoint, 3])
#     PyPlot.plot3D([mymesh.x[idx1], mymesh.x[idx2]], [mymesh.y[idx1], mymesh.y[idx2]], [mymesh.z[idx1], mymesh.z[idx2]], "r.-")
#     PyPlot.text3D(mymesh.x[idx1] .+ rand() / 30, mymesh.y[idx1] .+ rand() / 30, mymesh.z[idx1] .+ rand() / 30, "$idx1", color = "r", ha = "center", va = "center")
#     PyPlot.text3D(mymesh.x[idx2] .+ rand() / 30, mymesh.y[idx2] .+ rand() / 30, mymesh.z[idx2] .+ rand() / 30, "$idx2", color = "r", ha = "center", va = "center")
#     # sleep(0.1)
# end
# PyPlot.xlabel("x")
# PyPlot.ylabel("y")
# PyPlot.zlabel("z")

RPM = 34.0
omega = RPM * 2 * pi / 60
nVinf = 20
Vinfvec = collect(LinRange(4.0, 25.0, nVinf))
torque = zeros(nVinf)
thrust = zeros(nVinf)
torque2 = zeros(nVinf)
power = zeros(nVinf)
power2 = zeros(nVinf)
TSRvec = omega * Blade_Radius ./ Vinfvec
windpower = 0.5 * rho * Vinfvec .^ 3 * RefArea
for (i, Vinf) in enumerate(Vinfvec)
    # Vinf = 10.0
    println("$i of $(length(Vinfvec))")
    CP, Rp, Tp, Zp, alpha, cl_af, cd_af, Vloc, Re, thetavec, nstep, Fx_base, Fy_base, Fz_base,
        Mx_base, My_base, Mz_base, power[i], power2[i], torque[i] = OWENSAero.steadyTurb(; omega, Vinf)
    thrust[i] = mean(Fx_base)
    torque2[i] = mean(Mz_base)
end

##########################################
#### Torque Plot
##########################################

# SNL34m_3_6_Torque = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/3.6.csv", ',', skipstart = 0)

# exptorque = SNL34m_3_6_Torque[:, 2]
# greyidx = exptorque .< 0.25 * maximum(exptorque)
# regidx = exptorque .>= 0.25 * maximum(exptorque)
# torquegrey = exptorque[greyidx]
# torquereg = exptorque[regidx]

# Vinfarray = SNL34m_3_6_Torque[:, 1]
# Vinfgrey = Vinfarray[greyidx]
# Vinfreg = Vinfarray[regidx]

# PyPlot.figure()
# PyPlot.plot(Vinfvec, torque ./ 1000, color = plot_cycle[1], label = "OWENS Predictions")
# # PyPlot.plot(Vinfvec,torque2./1000,"+",label="DMS Aero Only2")
# # PyPlot.plot(Vinfgrey,torquegrey,".",color="0.5",label="Experimental < 25% Peak")
# # PyPlot.plot(Vinfreg,torquereg,"k.",label="Experimental > 25% Peak")
# PyPlot.plot(Vinfarray, exptorque, "k.", label = "VAWT Test Data")
# PyPlot.xlabel("Wind Speed (m/s)")
# PyPlot.ylabel("Torque (kN-m)")
# # PyPlot.grid(linestyle="--")
# PyPlot.xlim([0.0, 25.0])
# PyPlot.ylim([-25.0, 150.0])
# PyPlot.legend()
# PyPlot.savefig("$(path)/figs/34m_fig3_6.pdf", transparent = true)

##########################################
#### Power Plot
##########################################

# SNL34m_3_7_Power = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/3.7.csv", ',', skipstart = 0)

# PyPlot.figure()
# PyPlot.plot(Vinfvec, power ./ 1000, color = plot_cycle[1], label = "DMS Aero Only 1")
# # PyPlot.plot(Vinfvec,power2./1000,label="DMS Aero Only 2")
# PyPlot.plot(SNL34m_3_7_Power[:, 1], SNL34m_3_7_Power[:, 2], "k.", label = "Experimental")
# PyPlot.xlabel("Wind Speed (m/s)")
# PyPlot.ylabel("Power (kW)")
# PyPlot.xlim([0.0, 25.0])
# PyPlot.ylim([-100.0, 600.0])
# PyPlot.grid(linestyle = "--")
# PyPlot.legend()
# PyPlot.savefig("$(path)/figs/34m_fig3_7.pdf", transparent = true)

##########################################
#### CP Plot
##########################################

# TSRvecExpGrey = omega * Blade_Radius ./ SNL34m_3_6_Torque[greyidx, 1]
# windpowerexp = 0.5 * rho * SNL34m_3_6_Torque[greyidx, 1] .^ 3 * RefArea
# Cp2expGrey = SNL34m_3_6_Torque[greyidx, 2] * omega * 1000.0 ./ windpowerexp

# TSRvecExpReg = omega * Blade_Radius ./ SNL34m_3_6_Torque[regidx, 1]
# windpowerexp = 0.5 * rho * SNL34m_3_6_Torque[regidx, 1] .^ 3 * RefArea
# Cp2expReg = SNL34m_3_6_Torque[regidx, 2] * omega * 1000.0 ./ windpowerexp

# SNL34m_3_7_Cp = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/3.8.csv", ',', skipstart = 0)

# PyPlot.rc("figure.subplot", left = 0.15, bottom = 0.17, top = 0.9, right = 0.9)
# PyPlot.figure()
# PyPlot.plot(TSRvec, power ./ windpower, color = plot_cycle[2], label = "Aero Only (DMS)")
# PyPlot.plot(SNL34m_3_7_Cp[:, 1], SNL34m_3_7_Cp[:, 2] / 1000, "ko", label = "Exp. From CP Plot")
# # PyPlot.plot(TSRvecExpGrey,Cp2expGrey,".",color="0.5",label="Exp. From Torque < 0.25% Peak Torque")
# # PyPlot.plot(TSRvecExpReg,Cp2expReg,"k.",label="Exp. From Torque > 0.25% Peak Torque")
# PyPlot.xlabel("Tip Speed Ratio")
# PyPlot.ylabel("Power Coefficient")
# PyPlot.xlim([0, 13])
# PyPlot.ylim([0, 0.8])
# PyPlot.legend()
# PyPlot.savefig("$(path)/figs/34m_fig3_8.pdf", transparent = true)

##########################################
#### CT Plot
##########################################


# PyPlot.figure()
# PyPlot.plot(TSRvec, -thrust ./ (windpower ./ Vinfvec), color = plot_cycle[1], label = "DMS Aero Only")
# PyPlot.xlabel("Tip Speed Ratio")
# PyPlot.ylabel("Thrust Coefficient")
# PyPlot.xlim([0, 15])
# PyPlot.ylim([0, 0.9])
# PyPlot.legend()
# PyPlot.savefig("$(path)/figs/34m_figCT_34RPM.pdf", transparent = true)

##########################################
############ AeroElastic #################
##########################################

top_idx = 23 #Int(myjoint[7,2])
pBC = [
    1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0
    top_idx 1 0
    top_idx 2 0
    top_idx 3 0
    top_idx 4 0
    top_idx 5 0
]

if AModel == "AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;
    analysisType = structuralModel,
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
    JgearBox = 250.0, #(2.15e3+25.7)/12*1.35582*100,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    driveShaftProps = OWENS.DriveShaftProps(10000, 1.5e2), #8.636e5*1.35582*0.6
    OmegaInit = Omegaocp[1] / 60
)

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

feamodel = OWENS.FEAModel(;
    analysisType = structuralModel,
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI"
)


# Get Gravity Loads
inputs.Omegaocp = inputs.Omegaocp .* 0.0
inputs.OmegaInit = inputs.OmegaInit .* 0.0
inputs.Vinfocp = inputs.Vinfocp .* 0.0
feamodel.nlOn = true

## Returns data filled with e.g. eps[Nbld,N_ts,Nel_bld]
eps_x_grav, eps_z_grav, eps_y_grav, kappa_x_grav, kappa_y_grav, kappa_z_grav, t, FReactionHist_grav = OWENS.run34m_ad(inputs, feamodel, mymesh, myel, aeroForces, deformAero; steady = true)
