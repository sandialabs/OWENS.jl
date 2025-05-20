using Test
import OWENSFEA
import OWENS
import MAT

path = splitdir(@__FILE__)[1]

##############################################
# Setup
#############################################

file = MAT.matopen("$path/data/ROMInputs.mat")
mymesh2 = MAT.read(file,"mymesh")
myel2 = MAT.read(file,"myel")
myjoint = MAT.read(file,"myjoint")
close(file)

mymesh = OWENSFEA.Mesh(round.(Int,mymesh2["nodeNum"]),round.(Int,mymesh2["numEl"]),Int.(mymesh2["numNodes"]),
mymesh2["x"],mymesh2["y"],mymesh2["z"],round.(Int,mymesh2["elNum"]),Int.(mymesh2["conn"]),
Int.(mymesh2["type"]),Int.(mymesh2["meshSeg"]),mymesh2["structuralSpanLocNorm"],
Int.(mymesh2["structuralNodeNumbers"]),Int.(mymesh2["structuralElNumbers"]))

sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, length(myel2["props"]))

pr2 = myel2["props"]

for i=1:length(myel2["props"])
    sectionPropsArray[i] = OWENSFEA.SectionPropsArray(pr2[i]["ac"],pr2[i]["twist"],
    pr2[i]["rhoA"],pr2[i]["EIyy"],pr2[i]["EIzz"],pr2[i]["GJ"],pr2[i]["EA"],pr2[i]["rhoIyy"],
    pr2[i]["rhoIzz"],pr2[i]["rhoJ"],pr2[i]["zcm"],pr2[i]["ycm"],pr2[i]["a"],pr2[i]["EIyz"],
    pr2[i]["alpha1"],pr2[i]["alpha2"],pr2[i]["alpha3"],pr2[i]["alpha4"],pr2[i]["alpha5"],
    pr2[i]["alpha6"],pr2[i]["rhoIyz"],pr2[i]["b"],pr2[i]["a0"],pr2[i]["aeroCenterOffset"])
end
myel = OWENSFEA.El(sectionPropsArray,myel2["elLen"],myel2["psi"],myel2["theta"],myel2["roll"],myel2["rotationalEffects"])

deformTurb2(a;newOmega=-1,newVinf=-1, bld_x=-1, bld_z=-1, bld_twist=-1,accel_flap_in=-1,accel_edge_in=-1,gravity=-1) = 0

delta_t = 0.1
function aeroForcesDMS(t,azi)
    if t+delta_t < 0.2
        Fexternal = [1e6]
        Fdof = [20*6+1]
    else
        Fexternal = [0]
        Fdof = [1]
    end
    return Fexternal, Fdof, nothing, nothing, nothing, nothing
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
    dataOutputFilename = "none",
    tocp = [0.0,100000.1],
    numTS = 100,
    delta_t = 0.1,
    turbineStartup = 0,
    aeroLoadsOn = 1,
    generatorOn = false,
    ratedTorque = 100000.0,
    OmegaInit = 34.0/60,
    zeroTorqueGenSpeed = 0.010,
    pulloutRatio = 0.90,
    ratedGenSlipPerc = 9000.0,
    Omegaocp = [34.0,34.0] ./ 60)

feamodel = OWENSFEA.FEAModel(;analysisType = "ROM",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    numModes = 400,
    numNodes = mymesh.numNodes)

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist = OWENS.Unsteady(model;
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesDMS,deformAero=deformTurb2)
##########################################
#### Plot
##########################################

file = "$(path)/data/newmesh_34mout34m_ROMtransient.mat"
# mfile = MAT.matopen(file, "w")
# MAT.write(mfile, "t", collect(t))
# MAT.write(mfile, "FReactionHist", FReactionHist)
# MAT.close(mfile)
vars = MAT.matread(file)
told = vars["t"]
FReactionHistold = vars["FReactionHist"]



for i_t = round(Int,length(t)/20):length(t)
    # println(i_t)
    digits = 1e-4 
    @test isapprox(FReactionHistold[i_t,1],FReactionHist[i_t,1];atol= max(abs(FReactionHistold[i_t,1] * digits),digits))
    @test isapprox(FReactionHistold[i_t,2],FReactionHist[i_t,2];atol= max(abs(FReactionHistold[i_t,2] * digits),digits))
    @test isapprox(FReactionHistold[i_t,3],FReactionHist[i_t,3];atol= max(abs(FReactionHistold[i_t,3] * digits),digits))
    @test isapprox(FReactionHistold[i_t,4],FReactionHist[i_t,4];atol= max(abs(FReactionHistold[i_t,4] * digits),digits))
    @test isapprox(FReactionHistold[i_t,5],FReactionHist[i_t,5];atol= max(abs(FReactionHistold[i_t,5] * digits),digits))
    @test isapprox(FReactionHistold[i_t,6],FReactionHist[i_t,6];atol= max(abs(FReactionHistold[i_t,6] * digits),digits))
end


# using PyPlot
# PyPlot.pygui(true)
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

# for idof = 1:6
#     PyPlot.figure()
#     # PyPlot.plot(thetavec[1,:]/omega,Mz_base/1000,color=plot_cycle[1],label="DMS Aero Only")
#     PyPlot.plot(told',FReactionHistold[:,idof]/1000,".-",color=plot_cycle[2],label="Old")
#     PyPlot.plot(t,FReactionHist[:,idof]/1000,"+-",color=plot_cycle[1],label="New")
#     PyPlot.xlabel("Time (s)")
#     # PyPlot.xlim([0,50])
#     PyPlot.ylabel("dof $idof (kN-m)")
#     PyPlot.legend()
#     # PyPlot.savefig("$(path)/figs/34mRomTest.pdf",transparent = true)
# end