path = splitdir(@__FILE__)[1]

# import OWENS
include("$path/../src/OWENS.jl")
# mat_path = string("$path/../src/","/Matlab_Cxx")
# using MATLAB #if this is used, must run from src location
# mat"addpath($mat_path)"
using Test
import HDF5
import PyPlot

mutable struct PlatformProp
    fileRoot
    MassDiag
    StiffDiag_Nm_deg
    StiffDiag
end

function test_owens(test_transient,test_modal,test_flutter)

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
    bmOwens = "$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt"
    # append this name to the end of the saved files
    outFileExt = "_15mTowerExt_NOcentStiff"

    # convert rotational stiffness to N-m/rad
    StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]
    convRotStiff = [1 1 1 180/pi 180/pi 180/pi]
    StiffDiag = StiffDiag_Nm_deg .* convRotStiff

    # filename root to save the created nodal file
    platformProp = PlatformProp("$path/data/input_files_test/1_FourColumnSemi_2ndPass",
    [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9],
    [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6],
    StiffDiag)

    ## ************************************************************************
    # perform the transient simulations using OWENS
    # *************************************************************************


    # define the filename saving convention
    fname = string(platformProp.fileRoot,outFileExt)

    # *********************************************************************
    # perform operations for the nodal file generation
    # *********************************************************************
    MassVal = platformProp.MassDiag
    StiffVal = platformProp.StiffDiag

    nodes = [1 1]
    cmkType = ["M6" "K6"]
    cmkValues = zeros(6,6,2)
    for dd = 1:6
        # set up mass matrix
        cmkValues[dd,dd,1] = MassVal[dd]
        # set up stiffness matrix
        cmkValues[dd,dd,2] = StiffVal[dd]
    end
    OWENS.writeOwensNDL(fname, nodes, cmkType, cmkValues) #TODO: make this a direct pass

    # *********************************************************************
    # perform operations for the aerodynamic forces file generation
    # *********************************************************************
    CACTUSfileRoot = "$path/data/input_files_test/DVAWT_2B_LCDT"
    OWENSfileRoot = bmOwens

    # *********************************************************************
    # run a modal analysis of the platform design
    # *********************************************************************
    operatingRPM = 7.2 # rpm
    Nrpm = 10    # number of rpm stations
    maxRPM = 10
    Nmodes = 4  # number of modes to calculate/extract: TODO: Since we can only use eig, push this change through
    timeStep = 2e-3
    timeSim = 0.1       # [sec]
    n_t = timeSim/timeStep # length of time vector
    timeArray = [0, timeSim+1]
    rpmArray  = [operatingRPM, operatingRPM]
    omegaArrayHz = rpmArray ./ 60
    omegaArrayHz2 = [7.1]/60

    if test_transient
        println("Running Transient")
        # Juno.@enter OWENS.owens(string(fname, ".owens"),"TNB", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz)
        OWENS.owens(string(fname, ".owens"),"TNB", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz)
    end

    if test_modal
        # Juno.@enter owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)
        OWENS.owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)
    end

    if test_flutter
        # freq,damp=owens(string(fname ".owens"),"FA", omegaArrayHz2, true, 1.2041, Nmodes)
    end

    println("Function Finished")
end

verify_transient = true
verify_modal = true
plotmodal = false
# Juno.@profiler test_owens(verify_transient,verify_modal,false)
test_owens(verify_transient,verify_modal,false)
# mat"verifyOWENSPost"

if verify_transient

    tol = 1e-5
    n_t = 50
    old_filename = "$path/data/input_files_test/UNIT_TEST_15mTower_transient_dvawt_c_2_lcdt.h5"
    new_filename = "$path/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.h5"

    if (time() - mtime(new_filename)) > 10000 #milliseconds
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

    t = HDF5.h5read(new_filename,"t")
    aziHist = HDF5.h5read(new_filename,"aziHist")
    OmegaHist = HDF5.h5read(new_filename,"OmegaHist")
    OmegaDotHist = HDF5.h5read(new_filename,"OmegaDotHist")
    gbHist = HDF5.h5read(new_filename,"gbHist")
    gbDotHist = HDF5.h5read(new_filename,"gbDotHist")
    gbDotDotHist = HDF5.h5read(new_filename,"gbDotDotHist")
    FReactionHist = HDF5.h5read(new_filename,"FReactionHist")
    rigidDof = HDF5.h5read(new_filename,"rigidDof")
    genTorque = HDF5.h5read(new_filename,"genTorque")
    genPower = HDF5.h5read(new_filename,"genPower")
    torqueDriveShaft = HDF5.h5read(new_filename,"torqueDriveShaft")
    uHist = HDF5.h5read(new_filename,"uHist")
    eps_xx_0_hist = HDF5.h5read(new_filename,"eps_xx_0_hist")
    eps_xx_z_hist = HDF5.h5read(new_filename,"eps_xx_z_hist")
    eps_xx_y_hist = HDF5.h5read(new_filename,"eps_xx_y_hist")
    gam_xz_0_hist = HDF5.h5read(new_filename,"gam_xz_0_hist")
    gam_xz_y_hist = HDF5.h5read(new_filename,"gam_xz_y_hist")
    gam_xy_0_hist = HDF5.h5read(new_filename,"gam_xy_0_hist")
    gam_xy_z_hist = HDF5.h5read(new_filename,"gam_xy_z_hist")

    @test isapprox(old_t,t,atol = tol)
    @test isapprox(old_aziHist,aziHist,atol = tol)
    @test isapprox(old_OmegaHist,OmegaHist,atol = tol)
    @test isapprox(old_OmegaDotHist,OmegaDotHist,atol = tol)
    @test isapprox(old_gbHist,gbHist,atol = tol)
    @test isapprox(old_gbDotHist,gbDotHist,atol = tol)
    @test isapprox(old_gbDotDotHist,gbDotDotHist,atol = tol)
    for ii = 1:length(FReactionHist)
        # println(ii)
        local digits = floor(log10(abs(old_FReactionHist[ii]))) #this way if the tol is 1e-5, then we are actually looking at significant digits, much better than comparing 1e-5 on a 1e6 large number, that's 11 significant digits!
        @test isapprox(old_FReactionHist[ii],FReactionHist[ii],atol=tol*10^digits)
    end
    @test isapprox(old_rigidDof,rigidDof,atol = tol)
    @test isapprox(old_genTorque,genTorque,atol = tol)
    @test isapprox(old_genPower,genPower,atol = tol)
    @test isapprox(old_torqueDriveShaft,torqueDriveShaft,atol = tol)
    @test isapprox(old_uHist,uHist,atol = tol)
    @test isapprox(old_eps_xx_0_hist,eps_xx_0_hist,atol = tol)
    @test isapprox(old_eps_xx_z_hist,eps_xx_z_hist,atol = tol)
    @test isapprox(old_eps_xx_y_hist,eps_xx_y_hist,atol = tol)
    @test isapprox(old_gam_xz_0_hist,gam_xz_0_hist,atol = tol)
    @test isapprox(old_gam_xz_y_hist,gam_xz_y_hist,atol = tol)
    @test isapprox(old_gam_xy_0_hist,gam_xy_0_hist,atol = tol)
    @test isapprox(old_gam_xy_z_hist,gam_xy_z_hist,atol = tol)

end

#MODAL
# @testset "modal" begin
# if verify_modal
old_filename = "$path/data/input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
new_filename = "$path/data/input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out"

if (time() - mtime(new_filename)) > 10 #seconds
    println("It is this many seconds old")
    println((time() - mtime(new_filename)))
    println("Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.")
end

#Reading function

numNodes = 82#mesh.numNodes

freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)
freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90 = OWENS.readResultsModalOut(new_filename,numNodes)

if plotmodal
    PyPlot.close("all")
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
# end
