# using PyPlot
# close("all")
using Test
import GyricFEA
import OWENS
import MAT

path = splitdir(@__FILE__)[1]

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

##############################################
# Setup
#############################################


file = MAT.matopen("$path/data/ROMInputs.mat")
mymesh2 = MAT.read(file,"mymesh")
myel2 = MAT.read(file,"myel")
myjoint = MAT.read(file,"myjoint")
close(file)

mymesh = GyricFEA.Mesh(mymesh2["nodeNum"],mymesh2["numEl"],mymesh2["numNodes"],
mymesh2["x"],mymesh2["y"],mymesh2["z"],mymesh2["elNum"],mymesh2["conn"],
mymesh2["type"],mymesh2["meshSeg"],mymesh2["structuralSpanLocNorm"],
mymesh2["structuralNodeNumbers"],mymesh2["structuralElNumbers"])

sectionPropsArray = Array{GyricFEA.SectionPropsArray, 1}(undef, length(myel2["props"]))

pr2 = myel2["props"]

for i=1:length(myel2["props"])
    sectionPropsArray[i] = GyricFEA.SectionPropsArray(pr2[i]["ac"],pr2[i]["twist"],
    pr2[i]["rhoA"],pr2[i]["EIyy"],pr2[i]["EIzz"],pr2[i]["GJ"],pr2[i]["EA"],pr2[i]["rhoIyy"],
    pr2[i]["rhoIzz"],pr2[i]["rhoJ"],pr2[i]["zcm"],pr2[i]["ycm"],pr2[i]["a"],pr2[i]["EIyz"],
    pr2[i]["alpha1"],pr2[i]["alpha2"],pr2[i]["alpha3"],pr2[i]["alpha4"],pr2[i]["alpha5"],
    pr2[i]["alpha6"],pr2[i]["rhoIyz"],pr2[i]["b"],pr2[i]["a0"],pr2[i]["aeroCenterOffset"])
end
myel = GyricFEA.El(sectionPropsArray,myel2["elLen"],myel2["psi"],myel2["theta"],myel2["roll"],myel2["rotationalEffects"])

deformTurb2(a;newOmega=-1,newVinf=-1) = 0

delta_t = 0.1
function aeroForcesDMS(t,azi)
    if t+delta_t < 0.2
        Fexternal = 1e6;
        Fdof = 20*6+1;
    else
        Fexternal = [];
        Fdof = [];
    end
    return Fexternal, Fdof
end

top_idx = Int(myjoint[7,2])
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
top_idx 5 0
top_idx 6 0]

model = OWENS.Inputs(;analysisType = "ROM",
    outFilename = "none",
    tocp = [0.0,100000.1],
    numTS = 100,
    delta_t = 0.1,
    turbineStartup = 0,
    generatorOn = false,
    ratedTorque = 100000.0,
    OmegaInit = 34.0/60,
    zeroTorqueGenSpeed = 0.010,
    pulloutRatio = 0.90,
    ratedGenSlipPerc = 9000.0,
    Omegaocp = [34.0,34.0] ./ 60)

feamodel = GyricFEA.FEAModel(;analysisType = "ROM",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    numModes = 400,
    numNodes = mymesh.numNodes)

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,eps_xx_0_hist,eps_xx_z_hist,
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb2)


##########################################
#### Plot
##########################################

file = "$(path)/data/newmesh_34mout34m_ROMtransient.mat"
# mfile = MAT.matopen(file, "w")
# MAT.write(mfile, "t", t)
# MAT.write(mfile, "FReactionHist", FReactionHist)
# MAT.close(mfile)
vars = MAT.matread(file)
told = vars["t"]
FReactionHistold = vars["FReactionHist"]


for iDof = 1:6
    for i_t = round(Int,length(told)/20):length(told)
        atol = max(abs(FReactionHistold[i_t,iDof] * 0.01),1e-2)
        @test isapprox(FReactionHistold[i_t,iDof],FReactionHist[i_t,iDof];atol)
    end
end

# using PyPlot
# PyPlot.close("all")
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
#
# PyPlot.ion()
# PyPlot.figure()
# # PyPlot.plot(thetavec[1,:]/omega,Mz_base/1000,color=plot_cycle[1],label="DMS Aero Only")
# PyPlot.plot(told',FReactionHistold[:,3]/1000,color=plot_cycle[2],label="Old")
# PyPlot.plot(t,FReactionHist[:,3]/1000,color=plot_cycle[1],label="New")
# PyPlot.xlabel("Time (s)")
# # PyPlot.xlim([0,50])
# PyPlot.ylabel("Torque (kN-m)")
# PyPlot.legend()
# # PyPlot.savefig("$(path)/figs/34mRomTest.pdf",transparent = true)
