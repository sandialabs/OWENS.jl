path = splitdir(@__FILE__)[1]

import OWENS

using Test
import HDF5

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


test_transient = true
test_modal = true
test_flutter = false
println("Starting")

## ****************** COORDINATE SYSTEM DEFINITION *********************** #
# 1 - x -  surge (pointed downstream)
# 2 - y -  sway  (right hand rule)
# 3 - z -  heave (pointed upwards)
# 4 - Ox - roll
# 5 - Oy - pitch
# 6 - Oz - yaw
# *********************************************************************** #

# use this benchmark file
bmOwens = "$path/data/_15mTower_transient_dvawt_c_2_lcdt"
# append this name to the end of the saved files
outFileExt = "_15mTowerExt_NOcentStiff"

# convert rotational stiffness to N-m/rad
StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]
convRotStiff = [1 1 1 180/pi 180/pi 180/pi]
platStiffDiag = StiffDiag_Nm_deg .* convRotStiff

# filename root to save the created nodal file
platfileRoot = "$path/data/1_FourColumnSemi_2ndPass"
platMassDiag = [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9]
platStiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]

# external dependencies
hydrodynLib = []#"$path/../bin/HydroDyn_c_lib_x64"
moordynLib = []#"$path/../bin/MoorDyn_c_lib_x64"

# define the filename saving convention
fname = string(platfileRoot,outFileExt)

# *********************************************************************
# perform operations for the nodal file generation
# *********************************************************************
MassVal = platMassDiag
StiffVal = platStiffDiag

# nodes = [1 1]
# cmkType = ["M6" "K6"]
# cmkValues = zeros(6,6,2)
# for dd = 1:6
#     # set up mass matrix
#     cmkValues[dd,dd,1] = MassVal[dd]
#     # set up stiffness matrix
#     cmkValues[dd,dd,2] = StiffVal[dd]
# end
# OWENS.writeOwensNDL(fname, nodes, cmkType, cmkValues)

# *********************************************************************
# perform operations for the aerodynamic forces file generation
# *********************************************************************
CACTUSfileRoot = "$path/data/DVAWT_2B_LCDT"
OWENSfileRoot = bmOwens

operatingRPM = 7.2 # rpm
Nrpm = 10    # number of rpm stations
maxRPM = 10
Nmodes = 4  # number of modes to calculate/extract:
timeStep = 2e-3
timeSim = 0.1       # [sec]
n_t = timeSim/timeStep # length of time vector
timeArray = [0, timeSim+1]
rpmArray  = [operatingRPM, operatingRPM]
omegaArrayHz = rpmArray ./ 60
omegaArrayHz2 = [7.1]/60

## ************************************************************************
# perform the transient simulations using OWENS
# *************************************************************************
if test_transient
    println("Running Transient")

    OWENS.owens(string(fname, ".owens"),"TNB";iterationType="DI", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz, hydrodynLib=hydrodynLib, moordynLib=moordynLib)

    # Perform Tests
    tol = 1e-5
    n_t = 50
    old_filename = "$path/data/UNIT_TEST_15mTower_transient_dvawt_c_2_lcdt.h5"
    new_filename = "$path/data/_15mTower_transient_dvawt_c_2_lcdt.h5"

    if (time() - mtime(new_filename)) >0.70
        println("Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.")
    end

    old_t = HDF5.h5read(old_filename,"t")
    old_aziHist = HDF5.h5read(old_filename,"aziHist")
    old_OmegaHist = HDF5.h5read(old_filename,"OmegaHist")
    old_OmegaDotHist = HDF5.h5read(old_filename,"OmegaDotHist")
    old_gbHist = HDF5.h5read(old_filename,"gbHist")
    old_gbDotHist = HDF5.h5read(old_filename,"gbDotHist")
    old_gbDotDotHist = HDF5.h5read(old_filename,"gbDotDotHist")
    old_FReactionHist = HDF5.h5read(old_filename,"FReactionHist")
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

    t = HDF5.h5read(new_filename,"t")
    aziHist = HDF5.h5read(new_filename,"aziHist")
    OmegaHist = HDF5.h5read(new_filename,"OmegaHist")
    OmegaDotHist = HDF5.h5read(new_filename,"OmegaDotHist")
    gbHist = HDF5.h5read(new_filename,"gbHist")
    gbDotHist = HDF5.h5read(new_filename,"gbDotHist")
    gbDotDotHist = HDF5.h5read(new_filename,"gbDotDotHist")
    FReactionHist = HDF5.h5read(new_filename,"FReactionHist")
    genTorque = HDF5.h5read(new_filename,"genTorque")
    genPower = HDF5.h5read(new_filename,"genPower")
    torqueDriveShaft = HDF5.h5read(new_filename,"torqueDriveShaft")
    uHist = HDF5.h5read(new_filename,"uHist")
    eps_xx_0_hist = HDF5.h5read(new_filename,"epsilon_x_hist")
    eps_xx_z_hist = HDF5.h5read(new_filename,"kappa_y_hist")
    eps_xx_y_hist = HDF5.h5read(new_filename,"kappa_z_hist")
    gam_xz_0_hist = HDF5.h5read(new_filename,"epsilon_z_hist")
    gam_xz_y_hist = HDF5.h5read(new_filename,"kappa_x_hist")
    gam_xy_0_hist = HDF5.h5read(new_filename,"epsilon_y_hist")
    maxT_idx = length(t)
    @test isapprox(old_t[1:maxT_idx],t,atol = tol)
    @test isapprox(old_aziHist[1:maxT_idx],aziHist,atol = tol)
    @test isapprox(old_OmegaHist[1:maxT_idx],OmegaHist,atol = tol)
    @test isapprox(old_OmegaDotHist[1:maxT_idx],OmegaDotHist,atol = tol)
    @test isapprox(old_gbHist[1:maxT_idx],gbHist,atol = tol)
    @test isapprox(old_gbDotHist[1:maxT_idx],gbDotHist,atol = tol)
    @test isapprox(old_gbDotDotHist[1:maxT_idx],gbDotDotHist,atol = tol)
    # for ii = 1:6
    #     PyPlot.figure()
    #     PyPlot.plot(LinRange(0,1,length(old_FReactionHist[:,1])),old_FReactionHist[:,ii],label="old")
    #     PyPlot.plot(LinRange(0,1,length(FReactionHist[:,1])),FReactionHist[:,ii],label="new")
    #     PyPlot.ylabel("Freaction $ii")
    #     PyPlot.legend()
    # end
    for ii = 1:length(FReactionHist[:,1])
        for jj = 1:6#length(FReactionHist[1,:])
            local digits = floor(log10(abs(old_FReactionHist[ii,jj]))) #this way if the tol is 1e-5, then we are actually looking at significant digits, much better than comparing 1e-5 on a 1e6 large number, that's 11 significant digits!
            @test isapprox(old_FReactionHist[ii,jj],FReactionHist[ii,jj],atol=tol*10^(digits+3))
        end
    end
    @test isapprox(old_genTorque[1:maxT_idx],genTorque,atol = tol)
    @test isapprox(old_genPower[1:maxT_idx],genPower,atol = tol)
    @test isapprox(old_torqueDriveShaft[1:maxT_idx],torqueDriveShaft,atol = tol)
    @test isapprox(old_uHist[:,1:maxT_idx],collect(uHist'),atol = 1e-4)
    @test isapprox(old_eps_xx_0_hist[:,:,1:maxT_idx],eps_xx_0_hist,atol = tol)
    @test isapprox(old_eps_xx_z_hist[:,:,1:maxT_idx],eps_xx_z_hist,atol = tol)
    @test isapprox(old_eps_xx_y_hist[:,:,1:maxT_idx],eps_xx_y_hist,atol = tol)
    @test isapprox(old_gam_xz_0_hist[:,:,1:maxT_idx],gam_xz_0_hist,atol = tol)
    @test isapprox(old_gam_xz_y_hist[:,:,1:maxT_idx],gam_xz_y_hist,atol = tol)
    @test isapprox(old_gam_xy_0_hist[:,:,1:maxT_idx],gam_xy_0_hist,atol = tol)
end

# *********************************************************************
# run a modal analysis of the platform design
# *********************************************************************
if test_modal

    OWENS.owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)

    # if verify_modal
    old_filename = "$path/data/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
    new_filename = "$path/data/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out"

    #Reading function

    numNodes = 82#mesh.numNodes

    freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)
    freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90 = OWENS.readResultsModalOut(new_filename,numNodes)

    # tol = 1e-6 # This is covered by the campbell diagram test
    # for imode = 1:length(freq)
    #     used_tol = max(tol*freq[imode],tol) #don't enforce 1e-6 precision on a 1e6 number when we want 6 digits and not 12 digits of precision, also limit it for small number errors
    #     @test isapprox(freqOLD[imode],freq[imode],atol = used_tol)
    #     used_tol = max(tol*damp[imode],tol)
    #     @test isapprox(dampOLD[imode],damp[imode],atol = used_tol)
    # end

    tol = 1e-1
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
        for inode = 1:length(U_x_0OLD[:,1])
            if isapprox(abs.(U_x_0OLD[inode,imode]),abs.(U_x_0[inode,imode]),atol = tol)
                global U_x_0pass += 1
            end
            if isapprox(abs.(U_y_0OLD[inode,imode]),abs.(U_y_0[inode,imode]),atol = tol)
                global U_y_0pass += 1
            end
            if isapprox(abs.(U_z_0OLD[inode,imode]),abs.(U_z_0[inode,imode]),atol = tol)
                global U_z_0pass += 1
            end
            if isapprox(abs.(theta_x_0OLD[inode,imode]),abs.(theta_x_0[inode,imode]),atol = tol)
                global theta_x_0pass += 1
            end
            if isapprox(abs.(theta_y_0OLD[inode,imode]),abs.(theta_y_0[inode,imode]),atol = tol)
                global theta_y_0pass += 1
            end
            if isapprox(abs.(theta_z_0OLD[inode,imode]),abs.(theta_z_0[inode,imode]),atol = tol)
                global theta_z_0pass += 1
            end
            if isapprox(abs.(U_x_90OLD[inode,imode]),abs.(U_x_90[inode,imode]),atol = tol)
                global U_x_90pass += 1
            end
            if isapprox(abs.(U_y_90OLD[inode,imode]),abs.(U_y_90[inode,imode]),atol = tol)
                global U_y_90pass += 1
            end
            if isapprox(abs.(U_z_90OLD[inode,imode]),abs.(U_z_90[inode,imode]),atol = tol)
                global U_z_90pass += 1
            end
            if isapprox(abs.(theta_x_90OLD[inode,imode]),abs.(theta_x_90[inode,imode]),atol = tol)
                global theta_x_90pass += 1
            end
            if isapprox(abs.(theta_y_90OLD[inode,imode]),abs.(theta_y_90[inode,imode]),atol = tol)
                global theta_y_90pass += 1
            end
            if isapprox(abs.(theta_z_90OLD[inode,imode]),abs.(theta_z_90[inode,imode]),atol = tol)
                global theta_z_90pass += 1
            end
        end
    end

    # at least 80 percent of the modeshapes are identical indicates (despite the recripocity of the solutions) that the analysis is adequate
    tol2 = 0.8
    @test U_x_0pass/length(U_x_0OLD)>tol2
    @test U_y_0pass/length(U_x_0OLD)>tol2
    @test U_z_0pass/length(U_x_0OLD)>tol2
    @test theta_x_0pass/length(U_x_0OLD)>tol2
    @test theta_y_0pass/length(U_x_0OLD)>tol2
    @test theta_z_0pass/length(U_x_0OLD)>tol2
    @test U_x_90pass/length(U_x_0OLD)>tol2
    @test U_y_90pass/length(U_x_0OLD)>tol2
    @test U_z_90pass/length(U_x_0OLD)>tol2
    @test theta_x_90pass/length(U_x_0OLD)>tol2
    @test theta_y_90pass/length(U_x_0OLD)>tol2
    @test theta_z_90pass/length(U_x_0OLD)>tol2
    # end

end

if test_flutter
    # freq,damp=owens(string(fname ".owens"),"FA", omegaArrayHz2, true, 1.2041, Nmodes)
end

println("Function Finished")
