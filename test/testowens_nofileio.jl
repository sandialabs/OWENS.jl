# using PyPlot
# close("all")
using Test
import HDF5
import PyPlot
# import OWENS
path,_ = splitdir(@__FILE__)
include("$(path)/../src/OWENS.jl")
#TODO: ort file, nodal file, element file, initial conditions, and blade file

##############################################
# Setup
#############################################

# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL5MW_bld_z = [15.0, 21.61004296, 28.20951408, 28.2148, 34.81955704, 41.4296, 48.03964296, 54.63911408, 61.24915704, 67.8592, 74.46924296, 81.06871408, 87.67875704, 94.2888, 100.89884296, 107.49831408, 114.10835704, 120.7184, 127.32844296, 133.92791408, 133.9332, 140.53795704, 147.148].-15.0
SNL5MW_bld_x = -[0.0, -10.201, -20.361, -20.368290684, -29.478, -36.575, -42.579, -47.177, -50.555, -52.809, -53.953, -54.014, -53.031, -51.024, -47.979, -43.942, -38.768, -32.91, -25.587, -17.587, -17.580079568, -8.933, 8.0917312607e-15]


mymesh,myort,myjoint = OWENS.create_mesh(;Ht = 15.0, #tower height before blades attach
Hb = 147.148-15.0, #blade height
R = 54.014, # m bade radius
nstrut = 2,
strut_mout_ratio = 0.1, #distance from top/bottom
ntelem = 20, #tower elements
nbelem = 20, #blade elements
nselem = 2,  #strut elements
bshapex=SNL5MW_bld_x,
bshapez=SNL5MW_bld_z) #use defaults


bladeData,_,_,_ = OWENS.readBladeData("$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.bld") #reads overall blade data file
el = OWENS.readElementData(mymesh.numEl,"$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.el","$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.ort",bladeData) #read element data file (also reads orientation and blade data file associated with elements)

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
bladeData,
pBC = pBC,
nodalTerms = mynodalTerms,
numNodes = mymesh.numNodes)


##############################################
# Unsteady Test
#############################################

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,eps_xx_0_hist,eps_xx_z_hist,
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(model,mymesh,el)

#PROCEED WITH UNIT TEST
old_filename = "$path/data/input_files_test/UNIT_TEST_15mTower_transient_dvawt_c_2_lcdt.h5"

old_t = HDF5.h5read(old_filename,"t")
old_aziHist = HDF5.h5read(old_filename,"aziHist")
old_OmegaHist = HDF5.h5read(old_filename,"OmegaHist")
old_OmegaDotHist = HDF5.h5read(old_filename,"OmegaDotHist")
old_gbHist = HDF5.h5read(old_filename,"gbHist")
old_gbDotHist = HDF5.h5read(old_filename,"gbDotHist")
old_gbDotDotHist = HDF5.h5read(old_filename,"gbDotDotHist")
old_FReactionHist = HDF5.h5read(old_filename,"FReactionHist")
old_rigidDof = HDF5.h5read(old_filename,"rigidDof")
old_genTorque = HDF5.h5read(old_filename,"genTorque")
old_genPower = HDF5.h5read(old_filename,"genPower")
old_torqueDriveShaft = HDF5.h5read(old_filename,"torqueDriveShaft")
old_uHist = HDF5.h5read(old_filename,"uHist")
old_eps_xx_0_hist = HDF5.h5read(old_filename,"eps_xx_0_hist")
old_eps_xx_z_hist = HDF5.h5read(old_filename,"eps_xx_z_hist")
old_eps_xx_y_hist = HDF5.h5read(old_filename,"eps_xx_y_hist")
old_gam_xz_0_hist = HDF5.h5read(old_filename,"gam_xz_0_hist")
old_gam_xz_y_hist = HDF5.h5read(old_filename,"gam_xz_y_hist")
old_gam_xy_0_hist = HDF5.h5read(old_filename,"gam_xy_0_hist")
old_gam_xy_z_hist = HDF5.h5read(old_filename,"gam_xy_z_hist")

tol = 1e-4

@test isapprox(old_t,t,atol = tol)
@test isapprox(old_aziHist,aziHist,atol = tol)
@test isapprox(old_OmegaHist,OmegaHist,atol = tol)
@test isapprox(old_OmegaDotHist,OmegaDotHist,atol = tol)
@test isapprox(old_gbHist,gbHist,atol = tol)
@test isapprox(old_gbDotHist,gbDotHist,atol = tol)
@test isapprox(old_gbDotDotHist,gbDotDotHist,atol = tol)
for ii = 1:length(FReactionHist)
    if isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.01)
        @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.01)
    elseif isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.02)

        @warn "$ii tolerance is 2%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
        @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.02)
    elseif isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.1)
        @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.1)
        @warn "$ii tolerance is 10%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
    else
        @warn "$ii tolerance is 30%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
        @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.3)

    end
end
@test isapprox(old_rigidDof,rigidDof,atol = tol)
@test isapprox(old_genTorque,genTorque,atol = tol)
@test isapprox(old_genPower,genPower,atol = tol)
@test isapprox(old_torqueDriveShaft,torqueDriveShaft,atol = tol)
for ii = 1:length(old_uHist)
    @test isapprox(old_uHist[ii],uHist[ii],atol=tol)
end
@test isapprox(old_eps_xx_0_hist,eps_xx_0_hist,atol = tol)
@test isapprox(old_eps_xx_z_hist,eps_xx_z_hist,atol = tol)
@test isapprox(old_eps_xx_y_hist,eps_xx_y_hist,atol = tol)
@test isapprox(old_gam_xz_0_hist,gam_xz_0_hist,atol = tol)
@test isapprox(old_gam_xz_y_hist,gam_xz_y_hist,atol = tol)
@test isapprox(old_gam_xy_0_hist,gam_xy_0_hist,atol = tol)
@test isapprox(old_gam_xy_z_hist,gam_xy_z_hist,atol = tol)

# if testModal
##############################################
# Modal Test
#############################################
maxRPM = 10
Omega=0.5*maxRPM*2*pi/60
OmegaStart = 0.0
displInitGuess = zeros(mymesh.numNodes*6)
meshFile = "$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh"
mesh = OWENS.readMesh(meshFile)

mymodel = OWENS.Model(;analysisType = "M",
        outFilename = "none",
        joint = myjoint,
        platformTurbineConnectionNodeNumber = 1,
        bladeData,
        pBC = pBC,
        nodalTerms = mynodalTerms,
        numNodes = mymesh.numNodes)

freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENS.Modal(mymodel,mesh,el,displInitGuess,Omega,OmegaStart)

old_filename = "$path/data/input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
#Reading function

numNodes = 82#mesh.numNodes

freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)

# if true
#     PyPlot.close("all")
#     println("Plotting Modes")
#     Ndof = 10
#     savePlot = true
#
#
#     for df = 1:Ndof
#         OWENS.viz("$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh",new_filename,df,10)
#         if savePlot # save the plot
#             PyPlot.savefig(string(new_filename[1:end-4],"_MODE$(df)newplot.pdf"),transparent = true)
#         else # flip through the plots visually
#             sleep(0.1)
#         end
#         PyPlot.close("all")
#     end
#
#     for df = 1:Ndof
#         OWENS.viz("$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh",old_filename,df,10)
#         if savePlot # save the plot
#             PyPlot.savefig(string(old_filename[1:end-4],"_MODE$(df)newplot.pdf"),transparent = true)
#         else # flip through the plots visually
#             sleep(0.1)
#         end
#         PyPlot.close("all")
#     end
# println("MODAL PLOTTING COMPLETE")
#
# end

tol = 1e-6
for imode = 1:length(freq)
    used_tol = max(tol*freq[imode],tol) #don't enforce 1e-6 precision on a 1e6 number when we want 6 digits and not 12 digits of precision, also limit it for small number errors
    @test isapprox(freqOLD[imode],freq[imode],atol = used_tol)
    used_tol = max(tol*damp[imode],tol)
    @test isapprox(dampOLD[imode],damp[imode],atol = used_tol)
end

tol = 1e-2
U_x_0pass = 0
U_y_0pass = 0
U_z_0pass = 0
theta_x_0pass = 0
theta_y_0pass = 0
theta_z_0pass = 0
U_x_90pass = 0
U_y_90pass = 0
U_z_90pass = 0
theta_x_90pass = 0
theta_y_90pass = 0
theta_z_90pass = 0

for imode = 1:length(U_x_0OLD[1,:])
    # println("mode: $imode")
    # for inode = 1:numNodes
    # println("node $inode")
    if isapprox(abs.(U_x_0OLD[:,imode]),abs.(U_x_0[:,imode]),atol = tol)
        global U_x_0pass += 1
    end
    if isapprox(abs.(U_y_0OLD[:,imode]),abs.(U_y_0[:,imode]),atol = tol)
        global U_y_0pass += 1
    end
    if isapprox(abs.(U_z_0OLD[:,imode]),abs.(U_z_0[:,imode]),atol = tol)
        global U_z_0pass += 1
    end
    if isapprox(abs.(theta_x_0OLD[:,imode]),abs.(theta_x_0[:,imode]),atol = tol)
        global theta_x_0pass += 1
    end
    if isapprox(abs.(theta_y_0OLD[:,imode]),abs.(theta_y_0[:,imode]),atol = tol)
        global theta_y_0pass += 1
    end
    if isapprox(abs.(theta_z_0OLD[:,imode]),abs.(theta_z_0[:,imode]),atol = tol)
        global theta_z_0pass += 1
    end
    if isapprox(abs.(U_x_90OLD[:,imode]),abs.(U_x_90[:,imode]),atol = tol)
        global U_x_90pass += 1
    end
    if isapprox(abs.(U_y_90OLD[:,imode]),abs.(U_y_90[:,imode]),atol = tol)
        global U_y_90pass += 1
    end
    if isapprox(abs.(U_z_90OLD[:,imode]),abs.(U_z_90[:,imode]),atol = tol)
        global U_z_90pass += 1
    end
    if isapprox(abs.(theta_x_90OLD[:,imode]),abs.(theta_x_90[:,imode]),atol = tol)
        global theta_x_90pass += 1
    end
    if isapprox(abs.(theta_y_90OLD[:,imode]),abs.(theta_y_90[:,imode]),atol = tol)
        global theta_y_90pass += 1
    end
    if isapprox(abs.(theta_z_90OLD[:,imode]),abs.(theta_z_90[:,imode]),atol = tol)
        global theta_z_90pass += 1
    end
    # end
end

# at least 70 percent of the modeshapes are identical indicates (despite the recripocity of the solutions) that the analysis is adequate

@test U_x_0pass/length(U_x_0OLD[1,:])>0.70
@test U_y_0pass/length(U_x_0OLD[1,:])>0.89
@test U_z_0pass/length(U_x_0OLD[1,:])>0.88
@test theta_x_0pass/length(U_x_0OLD[1,:])>0.90
@test theta_y_0pass/length(U_x_0OLD[1,:])>0.90
@test theta_z_0pass/length(U_x_0OLD[1,:])>0.90
@test U_x_90pass/length(U_x_0OLD[1,:])>0.70
@test U_y_90pass/length(U_x_0OLD[1,:])>0.70
@test U_z_90pass/length(U_x_0OLD[1,:])>0.70
@test theta_x_90pass/length(U_x_0OLD[1,:])>0.80
@test theta_y_90pass/length(U_x_0OLD[1,:])>0.80
@test theta_z_90pass/length(U_x_0OLD[1,:])>0.80
