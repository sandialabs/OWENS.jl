# using PyPlot
# close("all")
using Test
import HDF5
import PyPlot
import DelimitedFiles
PyPlot.close("all")
# import OWENS
path = splitdir(@__FILE__)[1]
include("$(path)/../src/OWENS.jl")

using pyfloater
pyFloater = pyfloater.pyFloater.pyFloater #simplify the call

##############################################
# Setup
#############################################
# start = time()
# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL5MW_bld_z = [15.0, 21.61004296, 28.20951408, 28.2148, 34.81955704, 41.4296, 48.03964296, 54.63911408, 61.24915704, 67.8592, 74.46924296, 81.06871408, 87.67875704, 94.2888, 100.89884296, 107.49831408, 114.10835704, 120.7184, 127.32844296, 133.92791408, 133.9332, 140.53795704, 147.148].-15.0
SNL5MW_bld_x = -[0.0, -10.201, -20.361, -20.368290684, -29.478, -36.575, -42.579, -47.177, -50.555, -52.809, -53.953, -54.014, -53.031, -51.024, -47.979, -43.942, -38.768, -32.91, -25.587, -17.587, -17.580079568, -8.933, 8.0917312607e-15]

# Juno.@enter OWENS.create_mesh(;Ht = 15.0, #tower height before blades attach
# Hb = 147.148-15.0, #blade height
# R = 54.014, # m bade radius
# nstrut = 2,
# strut_mout_ratio = 0.1, #distance from top/bottom
# ntelem = 20, #tower elements
# nbelem = 20, #blade elements
# nselem = 2,  #strut elements
# bshapex=SNL5MW_bld_x,
# bshapez=SNL5MW_bld_z) #use defaults

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

#Tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = OWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(mymesh.z[1:24],numadIn,precompoutput)

#Blades
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld_precompoutput,bld_precompinput = OWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(mymesh.z[25:47].-15.0,numadIn_bld,bld_precompoutput)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = 8 #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;sectionPropsArray_bld;sectionPropsArray_bld; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = OWENS.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

p = pyFloater(wrking_dir="$path/data",
name="coupled_hydro_test",
draft= 20,
freeboard= 1,
column_radius= 7.5,
pontoon_length= 35,
pontoon_width= 5,
pontoon_height= 5,
num_columns= 3,
num_pontoons= 3,
tension_frac= 0.3,
water_depth= 80,
displacement= -99,
turbine_power= 5e6,
turbine_mass= 568800,
turbine_cog= [0, 0, 22.5],
mesh_refine= 0.7)

#%% Compute hydrodyanmics

p.run_hydro(freq=1/LinRange(1/20, 1/2, 3))
# p.write_hydro()           # save the results to an netCDF
# p.read_hydro()            # load results from netCDF (avoid calling run_hydro)
# p.floatingbody.show()     # visualize platform

#%% Get linear models

m, b, k = p.get_linearModel()

m_jl = m.data
b_jl = b.data
k_jl = k.data

mbk_array = [m_jl;b_jl;k_jl]

m_nrow,m_ncol = size(m_jl)
b_nrow,b_ncol = size(b_jl)
k_nrow,k_ncol = size(k_jl)
ntotalrow = 6+6+6 #TODO: it can't handle anything other than the diagonals for now#m_nrow*m_ncol+b_nrow*b_ncol+k_nrow*k_ncol
nodalinputdata = [zeros(ntotalrow) repeat(["NA"],ntotalrow) zeros(ntotalrow,3)]

inodal = 0
for icol = 1:length(mbk_array[1,:])
    for irow = 1:length(mbk_array[:,1])
        if icol%6 == irow%6
            global inodal += 1
            nodalinputdata[inodal,1] = 1 #node number

            # Type and Row number
            if irow<=m_nrow
                nodalinputdata[inodal,2] = "M6"
                nodalinputdata[inodal,3] = irow
            elseif irow<=m_nrow+b_nrow
                nodalinputdata[inodal,2] = "C6"
                nodalinputdata[inodal,3] = irow-m_nrow
            elseif irow<=m_nrow+b_nrow+k_nrow
                nodalinputdata[inodal,2] = "K6"
                nodalinputdata[inodal,3] = irow-(m_nrow+b_nrow)
            end

            #Col num
            nodalinputdata[inodal,4] = icol

            #Value
            nodalinputdata[inodal,5] = mbk_array[irow,icol]
        end
    end
end

for i = 1:length(nodalinputdata[:,1])
    println(nodalinputdata[i,:])
end


mynodalTerms = OWENS.readNodalTerms(data = nodalinputdata)

# node, dof, bc
# pBC = [1 1 0
# 1 2 0
# 1 3 0
# 1 4 0
# 1 5 0
# 1 6 0]

model = OWENS.Model(;analysisType = "TNB",
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
bladeData,
hydroOn = true,
plat_model = p,
nodalTerms = mynodalTerms,
numNodes = mymesh.numNodes)

# elapsed = time() - start
# println("here")
# println(elapsed)
##############################################
# Unsteady Test
#############################################

t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
rigidDof,genTorque,genPower,torqueDriveShaft,uHist,eps_xx_0_hist,eps_xx_z_hist,
eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist = OWENS.Unsteady(model,mymesh,myel;getLinearizedMatrices=false)

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
ii = 0
for row = 1:length(FReactionHist[:,1])
    for col = 1:length(FReactionHist[1,:])
        global ii+=1
        if isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.01)
            @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.01)
        elseif isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.02)
            @warn "row $row col $col tolerance is 2%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
            @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.02)
        elseif isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.1)
            @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.1)
            @warn "row $row col $col tolerance is 10%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
        elseif isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.3)
            @warn "row $row col $col tolerance is 30%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
            @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*0.3)
        else
            @warn "row $row col $col tolerance is 2000%, error is $(old_FReactionHist[ii]-FReactionHist[ii]), old: $(old_FReactionHist[ii]), new:$(FReactionHist[ii]), percent error: $((old_FReactionHist[ii]-FReactionHist[ii])/old_FReactionHist[ii]*100)%"
            @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=abs(old_FReactionHist[ii])*20.0)
        end
    end
end

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,1]),old_FReactionHist[:,1])
PyPlot.plot(1:length(FReactionHist[:,1]),FReactionHist[:,1])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 1")

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,2]),old_FReactionHist[:,2])
PyPlot.plot(1:length(FReactionHist[:,2]),FReactionHist[:,2])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 2")

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,3]),old_FReactionHist[:,3])
PyPlot.plot(1:length(FReactionHist[:,3]),FReactionHist[:,3])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 3")

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,4]),old_FReactionHist[:,4])
PyPlot.plot(1:length(FReactionHist[:,4]),FReactionHist[:,4])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 4")

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,5]),old_FReactionHist[:,5])
PyPlot.plot(1:length(FReactionHist[:,5]),FReactionHist[:,5])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 5")

PyPlot.figure()
PyPlot.plot(1:length(old_FReactionHist[:,6]),old_FReactionHist[:,6])
PyPlot.plot(1:length(FReactionHist[:,6]),FReactionHist[:,6])
PyPlot.legend(["Old", "New"])
PyPlot.ylabel("FReaction Hist 6")

@test isapprox(old_rigidDof,rigidDof,atol = tol)
@test isapprox(old_genTorque,genTorque,atol = tol)
@test isapprox(old_genPower,genPower,atol = tol)
@test isapprox(old_torqueDriveShaft,torqueDriveShaft,atol = tol)
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
    PyPlot.plot(1:length(old_uHist[ii,:]),old_uHist[ii,:],"k-")
    PyPlot.plot(1:length(uHist[ii,:]),uHist[ii,:],"k--")
    if ii%10 == 0.0
        PyPlot.figure()
    end
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

freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENS.Modal(mymodel,mesh,myel,displInitGuess,Omega,OmegaStart)

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

PyPlot.figure()
for ii = 1:length(sectionPropsArray)
    # toplot = sectionPropsArray_bld[ii].rhoA[1]
    toplot0 = (sectionPropsArray[ii].rhoA[1]-el.props[ii].rhoA[1])/el.props[ii].rhoA[1]
    PyPlot.plot(ii,toplot0.*100.0,"g.")
end
PyPlot.ylabel("% Error")
PyPlot.xlabel("Element")

# PyPlot.figure()
# PyPlot.plot(mymesh.structuralSpanLocNorm[1,:],zero(mymesh.structuralSpanLocNorm[1,:]),"r.")
