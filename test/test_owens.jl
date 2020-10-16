path,_ = splitdir(@__FILE__)


mat_path = string(path,"/../src/Matlab_Cxx")
using MATLAB #if this is used, must run from src location
mat"addpath($mat_path)"
include("$path/../old_src/owens.jl") #TODO: organize correctly

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
    bmOwens = "$path/../old_src/input_files_test/_15mTower_transient_dvawt_c_2_lcdt"
    # append this name to the end of the saved files
    outFileExt = "_15mTowerExt_NOcentStiff"

    # convert rotational stiffness to N-m/rad
    StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]
    convRotStiff = [1 1 1 180/pi 180/pi 180/pi]
    StiffDiag = StiffDiag_Nm_deg .* convRotStiff

    # filename root to save the created nodal file
    platformProp = PlatformProp("$path/../old_src/input_files_test/1_FourColumnSemi_2ndPass",
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
    # writeOwensNDL(fname, nodes, cmkType, cmkValues) #TODO: translate this

    # *********************************************************************
    # perform operations for the aerodynamic forces file generation
    # *********************************************************************
    CACTUSfileRoot = "$path/input_files_test/DVAWT_2B_LCDT"
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
        # Juno.@enter owens(string(fname, ".owens"),"TNB", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz)
        freq,damp=owens(string(fname, ".owens"),"TNB", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz)
    end

    if test_modal
        # Juno.@enter owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)
        freq,damp=owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)
    end

    if test_flutter
        # freq,damp=owens(string(fname ".owens"),"FA", omegaArrayHz2, true, 1.2041, Nmodes)
    end

    println("Function Finished")
end

test_owens(true,false,false)
mat"verifyOWENSPost"
