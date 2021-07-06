# using PyPlot
# close("all")
using Test
import HDF5
import PyPlot
import DelimitedFiles
import PrePostOWENS
PyPlot.close("all")
# import OWENS
path = splitdir(@__FILE__)[1]
include("$(path)/../src/OWENS.jl")

##############################################
# Setup
#############################################
# start = time()
# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL34_unit_xz = DelimitedFiles.readdlm("$(path)/data/SNL34m/SNL34m_unit_blade_shape.txt",'\t',skipstart = 0)

#Ensure the data is fully unitized since the hand picking process is only good to one or two significant digits.
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])

#Scale the turbine to the full dimensions
height = 41.9 #m
radius = 17.1 #m
SNL34Z = SNL34z*height
SNL34X = SNL34x*radius
SNL5MW_bld_z = SNL34Z
SNL5MW_bld_x = SNL34X

# Juno.@enter PrePostOWENS.create_mesh(;Ht = 15.0, #tower height before blades attach
# Hb = 147.148-15.0, #blade height
# R = 54.014, # m bade radius
# nstrut = 2,
# strut_mout_ratio = 0.1, #distance from top/bottom
# ntelem = 20, #tower elements
# nbelem = 20, #blade elements
# nselem = 2,  #strut elements
# bshapex=SNL5MW_bld_x,
# bshapez=SNL5MW_bld_z) #use defaults

mymesh,myort,myjoint = PrePostOWENS.create_mesh(;Ht = 15.0, #tower height before blades attach
Hb = height, #blade height
R = radius, # m bade radius
nstrut = 2,
strut_mout_ratio = 0.05, #distance from top/bottom
ntelem = 20, #tower elements
nbelem = 20, #blade elements
nselem = 2,  #strut elements
bshapex=SNL5MW_bld_x,
bshapez=SNL5MW_bld_z) #use defaults


nTwrElem = Int(mymesh.meshSeg[1])+1
nBldElem = Int(mymesh.meshSeg[2])+1
#Blades
NuMad_geom_xlscsv_file = "$path/data/SNL34m/SNL34mGeom.csv"
numadIn_bld = PrePostOWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/SNL34m/SNL34mMaterials.csv"
plyprops = PrePostOWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld_precompoutput,bld_precompinput = PrePostOWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = PrePostOWENS.getSectPropsFromPreComp(LinRange(0,1,nBldElem),numadIn_bld,bld_precompoutput)

#Tower #TODO: use actual tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn = PrePostOWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops = PrePostOWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = PrePostOWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = PrePostOWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn,precompoutput)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = 8 #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;sectionPropsArray_bld;sectionPropsArray_bld; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

# node, dof, bc.  Constrain the top and bottom of the turbine TODO: figure out if there is a less trial and error approach to hit the end boundary condition
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
70 1 0
70 2 0
70 3 0
70 4 0
70 5 0
70 6 0]
# owensfile = "$path/SNL34m/SNL34.owens", #controls output filenames: TODO: clean up
model = OWENS.Model(;analysisType = "TNB",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC = pBC,
bladeData = numadIn_bld,
numNodes = mymesh.numNodes)

##############################################
# Unsteady Test
#############################################

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,eps_xx_0_hist,eps_xx_z_hist,
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(model,mymesh,myel;getLinearizedMatrices=false)

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,1]),FReactionHist[:,1])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 1")

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,2]),FReactionHist[:,2])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 2")

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,3]),FReactionHist[:,3])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 3")

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,4]),FReactionHist[:,4])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 4")

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,5]),FReactionHist[:,5])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 5")

PyPlot.figure()
PyPlot.plot(1:length(FReactionHist[:,6]),FReactionHist[:,6])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 6")


PyPlot.figure()
PyPlot.ylabel("uhist")
for ii = 1:length(uHist[1,:])
    PyPlot.plot(1:length(uHist[ii,:]),uHist[ii,:],"k--")
    if ii%10 == 0.0
        PyPlot.figure()
    end
end

##############################################
# Modal Test
#############################################
maxRPM = 10
Omega=0.5*maxRPM*2*pi/60
OmegaStart = 0.0
displInitGuess = zeros(mymesh.numNodes*6)

mymodel = OWENS.Model(;analysisType = "M",
        outFilename = "none",
        joint = myjoint,
        owensfile = "$path/SNL34m/SNL34.owens", #controls output filenames: TODO: clean up
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        numNodes = mymesh.numNodes)

freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENS.Modal(mymodel,mymesh,myel,displInitGuess,Omega,OmegaStart)



numNodes = 82#mesh.numNodes
new_filename = "$path/data/input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out" #TODO: fix this filename
if true
    # PyPlot.close("all")
    println("Plotting Modes")
    Ndof = 10
    savePlot = true


    for df = 1:Ndof
        PrePostOWENS.viz("$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh",new_filename,df,10)
        if savePlot # save the plot
            PyPlot.savefig(string(new_filename[1:end-4],"_MODE$(df)newplot.pdf"),transparent = true)
        else # flip through the plots visually
            sleep(0.1)
        end
        # PyPlot.close("all")
    end

println("MODAL PLOTTING COMPLETE")

end

# ac
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].ac[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("ac")
PyPlot.xlabel("Element")

# twist
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].twist[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("twist")
PyPlot.xlabel("Element")

# rhoA
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoA[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoA")
PyPlot.xlabel("Element")

# EIyy
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIyy[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIyy")
PyPlot.xlabel("Element")

# EIzz
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIzz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIzz")
PyPlot.xlabel("Element")

# GJ
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].GJ[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("GJ")
PyPlot.xlabel("Element")

# EA
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EA[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EA")
PyPlot.xlabel("Element")

# rhoIyy
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIyy[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIyy")
PyPlot.xlabel("Element")

# rhoIzz
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIzz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIzz")
PyPlot.xlabel("Element")

# rhoJ
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoJ[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoJ")
PyPlot.xlabel("Element")

# zcm
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].zcm[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("zcm")
PyPlot.xlabel("Element")

# ycm
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].ycm[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("ycm")
PyPlot.xlabel("Element")

# a
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].a[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("a")
PyPlot.xlabel("Element")

# EIyz
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIyz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIyz")
PyPlot.xlabel("Element")

# rhoIyz
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIyz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIyz")
PyPlot.xlabel("Element")

# b
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].b[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("b")
PyPlot.xlabel("Element")

# a0
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].a0[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("a0")
PyPlot.xlabel("Element")

# aeroCenterOffset
PyPlot.figure()
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].aeroCenterOffset[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("aeroCenterOffset")
PyPlot.xlabel("Element")

# PyPlot.figure()
# PyPlot.plot(mymesh.structuralSpanLocNorm[1,:],zero(mymesh.structuralSpanLocNorm[1,:]),"r.")
