# using PyPlot
# close("all")
using Test
import HDF5
import PyPlot
import DelimitedFiles
import FLOWMath
PyPlot.close("all")
# import OWENS
path = splitdir(@__FILE__)[1]
include("$(path)/../src/OWENS.jl")

##############################################
# Setup
#############################################
# start = time()
# Use the SNL5MW as the baseline check

SNL34_unit_xz = DelimitedFiles.readdlm("$(path)/data/SNL34m/SNL34m_unit_blade_shape.txt",'\t',skipstart = 0)

#Ensure the data is fully unitized since the hand picking process is only good to one or two significant digits.
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])

#Scale the turbine to the full dimensions
height = 41.9 #m
radius = 17.1 #m
SNL34Z = SNL34z*height
SNL34X = SNL34x*radius

# Transitions
idx = [1,6,11,24,31,40]
z_transitions = [SNL34Z[idx[1]],SNL34Z[idx[2]],SNL34Z[idx[3]],SNL34Z[idx[4]],SNL34Z[idx[5]],SNL34Z[idx[6]]]
x_transitions = [SNL34X[idx[1]],SNL34X[idx[2]],SNL34X[idx[3]],SNL34X[idx[4]],SNL34X[idx[5]],SNL34X[idx[6]]]

# Create the spline, and for simplicity, we will not use the exact transition points in the mesh
npt = 30
spl_z = LinRange(0,height,npt)
spl_x = FLOWMath.akima(SNL34Z,SNL34X,spl_z)

mymesh, myort, myjoint = OWENS.create_mesh(;Ht = 0.1, #tower height before blades attach
    Hb = height, #blade height
    R = radius, # m bade radius
    nstrut = 2, #strut elements
    strut_mout_ratio = 0.025, #distance from top/bottom
    ntelem = 20, #tower elements
    nbelem = 20, #blade elements
    nselem = 2,  #strut elements
    bshapex = SNL34X, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = SNL34Z) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

nTwrElem = Int(mymesh.meshSeg[1])
nBldElem = Int(mymesh.meshSeg[2])
#Blades
NuMad_geom_xlscsv_file = "$path/data/SNL34m/SNL34mGeom.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/SNL34m/SNL34mMaterials.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld_precompoutput,bld_precompinput = OWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(LinRange(0,1,nBldElem),numadIn_bld,bld_precompoutput)

#Tower #TODO: use actual tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = OWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn,precompoutput)



#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = 8 #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;sectionPropsArray_bld;sectionPropsArray_bld; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = OWENS.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)


nodalinputdata = [1 "M6" 1 1 9.8088e6
1 "M6" 2 2 9.7811e6
1 "M6" 3 3 1.8914e7
1 "M6" 4 4 3.6351e9
1 "M6" 5 5 3.6509e9
1 "M6" 6 6 2.4362e9
1 "K6" 1 1 132900.0
1 "K6" 2 2 132900.0
1 "K6" 3 3 1.985e6
1 "K6" 4 4 2.2878204759573773e8
1 "K6" 5 5 2.2889663915476388e8
1 "K6" 6 6 6.165025875607658e7]

mynodalTerms = OWENS.readNodalTerms(data = nodalinputdata)

# node, dof, bc
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

model = OWENS.Model(;analysisType = "TNB",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC = pBC,
nodalTerms = mynodalTerms,
numNodes = mymesh.numNodes)

##############################################
# Unsteady Test
#############################################

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,eps_xx_0_hist,eps_xx_z_hist,
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(mymodel,mymesh,myel;getLinearizedMatrices=false)


PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,1]),old_FReactionHist[:,1])
PyPlot.plot(1:length(FReactionHist[:,1]),FReactionHist[:,1])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 1")

PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,2]),old_FReactionHist[:,2])
PyPlot.plot(1:length(FReactionHist[:,2]),FReactionHist[:,2])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 2")

PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,3]),old_FReactionHist[:,3])
PyPlot.plot(1:length(FReactionHist[:,3]),FReactionHist[:,3])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 3")

PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,4]),old_FReactionHist[:,4])
PyPlot.plot(1:length(FReactionHist[:,4]),FReactionHist[:,4])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 4")

PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,5]),old_FReactionHist[:,5])
PyPlot.plot(1:length(FReactionHist[:,5]),FReactionHist[:,5])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 5")

PyPlot.figure()
# PyPlot.plot(1:length(old_FReactionHist[:,6]),old_FReactionHist[:,6])
PyPlot.plot(1:length(FReactionHist[:,6]),FReactionHist[:,6])
#PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 6")

# for ii = 1:length(old_uHist)
#     if isapprox(old_uHist[ii],uHist[ii],atol=tol)
#         @test isapprox(old_uHist[ii],uHist[ii],atol=tol)
#     else
#         @warn "$ii tolerance is 1000%, error is $(old_uHist[ii]-uHist[ii]), old: $(old_uHist[ii]), new:$(uHist[ii]), percent error: $((old_uHist[ii]-uHist[ii])/old_uHist[ii]*100)%"
#         @test isapprox(old_uHist[ii],uHist[ii],atol=abs(old_uHist[ii])*10000.0)
#     end
# end
PyPlot.figure()
for ii = 1:length(uHist[1,:])
    # PyPlot.plot(1:length(old_uHist[ii,:]),old_uHist[ii,:],"k-")
    PyPlot.plot(1:length(uHist[ii,:]),uHist[ii,:],"k--")
    if ii%10 == 0.0
        PyPlot.figure()
    end
end

# if testModal
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
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        nodalTerms = mynodalTerms,
        numNodes = mymesh.numNodes)

freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENS.Modal(mymodel,mesh,myel,displInitGuess,Omega,OmegaStart)

numNodes = 82#mesh.numNodes

# freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)

el=myel #TODO...

# ac
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.ac[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].ac[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("ac")
PyPlot.xlabel("Element")

# twist
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.twist[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].twist[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("twist")
PyPlot.xlabel("Element")

# rhoA
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.rhoA[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoA[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoA")
PyPlot.xlabel("Element")

# EIyy
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.EIyy[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIyy[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIyy")
PyPlot.xlabel("Element")

# EIzz
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.EIzz[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIzz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIzz")
PyPlot.xlabel("Element")

# GJ
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.GJ[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].GJ[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("GJ")
PyPlot.xlabel("Element")

# EA
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.EA[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EA[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EA")
PyPlot.xlabel("Element")

# rhoIyy
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.rhoIyy[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIyy[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIyy")
PyPlot.xlabel("Element")

# rhoIzz
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.rhoIzz[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIzz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIzz")
PyPlot.xlabel("Element")

# rhoJ
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.rhoJ[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoJ[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoJ")
PyPlot.xlabel("Element")

# zcm
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.zcm[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].zcm[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("zcm")
PyPlot.xlabel("Element")

# ycm
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.ycm[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].ycm[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("ycm")
PyPlot.xlabel("Element")

# a
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.a[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].a[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("a")
PyPlot.xlabel("Element")

# EIyz
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.EIyz[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].EIyz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("EIyz")
PyPlot.xlabel("Element")

# rhoIyz
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.rhoIyz[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].rhoIyz[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("rhoIyz")
PyPlot.xlabel("Element")

# b
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.b[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].b[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("b")
PyPlot.xlabel("Element")

# a0
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.a0[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].a0[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("a0")
PyPlot.xlabel("Element")

# aeroCenterOffset
PyPlot.figure()
elplot = zeros(length(el.props))
for ii = 1:length(el.props)
    secprops = el.props[ii]
    elplot[ii] = secprops.aeroCenterOffset[1]
end
PyPlot.plot(1:length(elplot),elplot,"k.-")
myelplot = zeros(length(sectionPropsArray))
for ii = 1:length(sectionPropsArray)
    myelplot[ii] = sectionPropsArray[ii].aeroCenterOffset[1]
end
PyPlot.plot(1:length(myelplot),myelplot,"r.-")
PyPlot.ylabel("aeroCenterOffset")
PyPlot.xlabel("Element")
