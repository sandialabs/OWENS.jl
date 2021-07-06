using Test
import HDF5
import PyPlot
import DelimitedFiles
import PrePostOWENS
import GyricFEA
PyPlot.close("all")
# import OWENS
path = splitdir(@__FILE__)[1]
include("$(path)/../src/OWENS.jl")
#TODO:  initial conditions

##############################################
# Setup
#############################################

# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL5MW_bld_z = [15.0, 21.61004296, 28.20951408, 28.2148, 34.81955704, 41.4296, 48.03964296, 54.63911408, 61.24915704, 67.8592, 74.46924296, 81.06871408, 87.67875704, 94.2888, 100.89884296, 107.49831408, 114.10835704, 120.7184, 127.32844296, 133.92791408, 133.9332, 140.53795704, 147.148].-15.0
SNL5MW_bld_x = -[0.0, -10.201, -20.361, -20.368290684, -29.478, -36.575, -42.579, -47.177, -50.555, -52.809, -53.953, -54.014, -53.031, -51.024, -47.979, -43.942, -38.768, -32.91, -25.587, -17.587, -17.580079568, -8.933, 8.0917312607e-15]

H_plat = 10.0

mymesh,myort,myjoint = PrePostOWENS.create_mesh(;Ht = 15.0+H_plat, #tower height before blades attach
Hb = 147.148-15.0, #blade height
R = 54.014, # m bade radius
nstrut = 2,
strut_mout_ratio = 0.1, #distance from top/bottom
ntelem = 20, #tower elements
nbelem = 20, #blade elements
nselem = 2,  #strut elements
bshapex=SNL5MW_bld_x,
bshapez=SNL5MW_bld_z) #use defaults

PyPlot.figure()
PyPlot.plot(mymesh.x,mymesh.z,"k.-",markersize=10.0)
PyPlot.axis("equal")

PyPlot.figure()
PyPlot.plot(mymesh.y,mymesh.z,"k.-",markersize=10.0)
PyPlot.legend(["mymesh"])
# PyPlot.axis("equal")

#Tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn = PrePostOWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops = PrePostOWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = PrePostOWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = PrePostOWENS.getSectPropsFromPreComp(mymesh.z[1:25],numadIn,precompoutput)

#Blades
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = PrePostOWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops = PrePostOWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld_precompoutput,bld_precompinput = PrePostOWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = PrePostOWENS.getSectPropsFromPreComp(mymesh.z[26:48].-15.0,numadIn_bld,bld_precompoutput)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = 8 #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;sectionPropsArray_bld;sectionPropsArray_bld; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

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
1 "K6" 6 6 6.165025875607658e7
2 "K6" 1 1 Inf
2 "K6" 2 2 Inf
2 "K6" 3 3 Inf
2 "K6" 4 4 Inf
2 "K6" 5 5 Inf
2 "K6" 6 6 Inf]

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
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(model,mymesh,myel;getLinearizedMatrices=false)

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

freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENS.Modal(mymodel,mymesh,myel,displInitGuess,Omega,OmegaStart)

old_filename = "$path/data/input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
new_filename = "$path/data/input_files_test/NEW_1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
#Reading function

numNodes = mymesh.numNodes

freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)

if true
    PyPlot.close("all")
    println("Plotting Modes")
    Ndof = 10
    savePlot = false


    for df = 1:Ndof
        PrePostOWENS.viz("$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh",new_filename,df,10)
        if savePlot # save the plot
            PyPlot.savefig(string(new_filename[1:end-4],"_MODE$(df)newplot.pdf"),transparent = true)
        else # flip through the plots visually
            sleep(0.1)
        end
        PyPlot.close("all")
    end

    for df = 1:Ndof
        PrePostOWENS.viz("$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh",old_filename,df,10)
        if savePlot # save the plot
            PyPlot.savefig(string(old_filename[1:end-4],"_MODE$(df)newplot.pdf"),transparent = true)
        else # flip through the plots visually
            sleep(0.1)
        end
        PyPlot.close("all")
    end
println("MODAL PLOTTING COMPLETE")

end
