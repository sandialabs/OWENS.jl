function runOWENSWINDIO(modelopt, windio, path)
    if isa(windio, String)
        windio = Design_Data(windio)
    end

    if isa(modelopt, String)
        modelopt = ModelingOptions(modelopt)
    end

    # Create SetupOptions instance using the preprocessing function
    setup_opts = preprocess_windio_setup(modelopt, windio, path)

    # Extract environment parameters needed for gravity
    gravity = Float64.(windio[:environment][:gravity])
    if isa(gravity, Float64)
        gravityOn = [0, 0, gravity]
    else
        gravityOn = gravity
    end

    # Extract control strategy parameters
    controlStrategy = modelopt.OWENS_Options.controlStrategy
    DLCs = modelopt.DLC_Options.DLCs
    normalRPM = modelopt.OWENS_Options.Prescribed_RPM_RPM_controlpoints[1] ./ 60
    slowRPM = 1.0 #RPM

    # Handle control strategy
    if DLCs != "none"
        tocp, Omegaocp, turbineStartup, generatorOn, useGeneratorFunction =
            handle_control_strategy(
                controlStrategy,
                normalRPM,
                slowRPM,
                modelopt.OWENS_Options.Prescribed_RPM_time_controlpoints,
                modelopt.OWENS_Options.Prescribed_RPM_RPM_controlpoints,
            )
    else
        tocp = modelopt.OWENS_Options.Prescribed_RPM_time_controlpoints
        Omegaocp = modelopt.OWENS_Options.Prescribed_RPM_RPM_controlpoints ./ 60
        turbineStartup = modelopt.Drivetrain_Options.turbineStartup
        generatorOn = modelopt.Drivetrain_Options.generatorOn
        useGeneratorFunction = modelopt.Drivetrain_Options.useGeneratorFunction
    end

    println("controlStrategy: $controlStrategy")

    # Create drive shaft properties
    driveShaftProps = DriveShaftProps(
        modelopt.Drivetrain_Options.driveShaft_K,
        modelopt.Drivetrain_Options.driveShaft_C,
    )

    # Create inputs for the simulation
    inputs = OWENS.Inputs(;
        analysisType = modelopt.OWENS_Options.structuralModel,
        tocp,
        Omegaocp,
        tocp_Vinf = modelopt.OWENS_Options.Prescribed_Vinf_time_controlpoints,
        Vinfocp = modelopt.OWENS_Options.Prescribed_Vinf_Vinf_controlpoints,
        numTS = modelopt.OWENS_Options.numTS,
        delta_t = modelopt.OWENS_Options.delta_t,
        AD15On = modelopt.OWENS_Options.AeroModel == "AD",
        aeroLoadsOn = modelopt.OWENS_Options.aeroLoadsOn,
        turbineStartup,
        usingRotorSpeedFunction = modelopt.Drivetrain_Options.usingRotorSpeedFunction,
        driveTrainOn = modelopt.Drivetrain_Options.driveTrainOn,
        generatorOn,
        platformActive = modelopt.OWENS_Options.platformActive,
        topsideOn = modelopt.OWENS_Options.topsideOn,
        interpOrder = modelopt.OWENS_Options.interpOrder,
        hd_input_file = modelopt.OWENSOpenFASTWrappers_Options.hd_input_file,
        ss_input_file = modelopt.OWENSOpenFASTWrappers_Options.ss_input_file,
        md_input_file = modelopt.OWENSOpenFASTWrappers_Options.md_input_file,
        JgearBox = modelopt.Drivetrain_Options.JgearBox,
        gearRatio = modelopt.Drivetrain_Options.gearRatio,
        gearBoxEfficiency = modelopt.Drivetrain_Options.gearBoxEfficiency,
        useGeneratorFunction,
        generatorProps = modelopt.Drivetrain_Options.generatorProps,
        ratedTorque = modelopt.Drivetrain_Options.ratedTorque,
        zeroTorqueGenSpeed = modelopt.Drivetrain_Options.zeroTorqueGenSpeed,
        pulloutRatio = modelopt.Drivetrain_Options.pulloutRatio,
        ratedGenSlipPerc = modelopt.Drivetrain_Options.ratedGenSlipPerc,
        OmegaGenStart = modelopt.Drivetrain_Options.OmegaGenStart,
        omegaControl = false,
        OmegaInit = 7.2/60,
        aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv",
        owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens",
        potflowfile = "$(path)$(modelopt.OWENSOpenFASTWrappers_Options.potflowfile)",
        dataOutputFilename = modelopt.OWENS_Options.dataOutputFilename,
        numDOFPerNode = 6,
        bladeData = [],
        rigid = modelopt.OWENS_Options.rigid,
        driveShaftProps,
        TOL = modelopt.OWENS_Options.TOL,
        MAXITER = modelopt.OWENS_Options.MAXITER,
    )

    # Call setupOWENS with the setup options
    mymesh,
    myel,
    myort,
    myjoint,
    sectionPropsArray,
    mass_twr,
    mass_bld,
    stiff_twr,
    stiff_bld,
    bld_precompinput,
    bld_precompoutput,
    plyprops_bld,
    numadIn_bld,
    lam_U_bld,
    lam_L_bld,
    twr_precompinput,
    twr_precompoutput,
    plyprops_twr,
    numadIn_twr,
    lam_U_twr,
    lam_L_twr,
    aeroForces,
    deformAero,
    mass_breakout_blds,
    mass_breakout_twr,
    system,
    assembly,
    sections,
    AD15bldNdIdxRng,
    AD15bldElIdxRng = OWENS.setupOWENS(path; setup_options = setup_opts)

    # Apply boundary conditions
    pBC = [
        1 1 0
        1 2 0
        1 3 0
        1 4 0
        1 5 0
        1 6 0
    ]

    # Create FEA model
    feamodel = OWENS.FEAModel(;
        analysisType = modelopt.OWENS_Options.structuralModel,
        airDensity = setup_opts.aero.rho,
        pBC,
        nlOn = modelopt.OWENSFEA_Options.nlOn,
        gravityOn,
        numNodes = mymesh.numNodes,
        RayleighAlpha = modelopt.OWENSFEA_Options.RayleighAlpha,
        RayleighBeta = modelopt.OWENSFEA_Options.RayleighBeta,
        joint = myjoint,
        iterationType = modelopt.OWENSFEA_Options.iterationType,
        initCond = [],
        aeroElasticOn = modelopt.OWENSFEA_Options.aeroElasticOn,
        guessFreq = modelopt.OWENSFEA_Options.guessFreq,
        dataOutputFilename = modelopt.OWENS_Options.dataOutputFilename,
        jointTransform = 0.0,
        reducedDOFList = zeros(Int, 2),
        numDOFPerNode = 6,
        platformTurbineConnectionNodeNumber = modelopt.OWENSFEA_Options.platformTurbineConnectionNodeNumber,
        spinUpOn = modelopt.OWENSFEA_Options.spinUpOn,
        numModes = modelopt.OWENSFEA_Options.numModes,
        nlParams = 0,
        adaptiveLoadSteppingFlag = modelopt.OWENSFEA_Options.adaptiveLoadSteppingFlag,
        tolerance = modelopt.OWENSFEA_Options.tolerance,
        maxIterations = modelopt.OWENSFEA_Options.maxIterations,
        maxNumLoadSteps = modelopt.OWENSFEA_Options.maxNumLoadSteps,
        minLoadStepDelta = modelopt.OWENSFEA_Options.minLoadStepDelta,
        minLoadStep = modelopt.OWENSFEA_Options.minLoadStep,
        prescribedLoadStep = modelopt.OWENSFEA_Options.prescribedLoadStep,
        predef = modelopt.OWENSFEA_Options.predef,
        elementOrder = modelopt.OWENSFEA_Options.elementOrder,
        alpha = modelopt.OWENSFEA_Options.alpha,
        gamma = modelopt.OWENSFEA_Options.gamma,
        AddedMass_Coeff_Ca = modelopt.OWENSFEA_Options.AddedMass_Coeff_Ca,
        nodalTerms = 0.0,
    )

    # Run unsteady simulation
    println("Running Unsteady")
    t,
    aziHist,
    OmegaHist,
    OmegaDotHist,
    gbHist,
    gbDotHist,
    gbDotDotHist,
    FReactionHist,
    FTwrBsHist,
    genTorque,
    genPower,
    torqueDriveShaft,
    uHist,
    uHist_prp,
    epsilon_x_hist,
    epsilon_y_hist,
    epsilon_z_hist,
    kappa_x_hist,
    kappa_y_hist,
    kappa_z_hist,
    FPtfmHist,
    FHydroHist,
    FMooringHist,
    topFexternal_hist,
    rbDataHist = OWENS.Unsteady_Land(
        inputs;
        system,
        assembly,
        topModel = feamodel,
        topMesh = mymesh,
        topEl = myel,
        aero = aeroForces,
        deformAero,
    )

    # Output VTK files
    tsave_idx = 1:3:(modelopt.OWENS_Options.numTS-1)
    OWENS.OWENSVTK(
        modelopt.OWENS_Options.VTKsaveName,
        t,
        uHist,
        system,
        assembly,
        sections,
        aziHist,
        mymesh,
        myel,
        epsilon_x_hist,
        epsilon_y_hist,
        epsilon_z_hist,
        kappa_x_hist,
        kappa_y_hist,
        kappa_z_hist,
        FReactionHist,
        topFexternal_hist;
        tsave_idx,
    )

    # Calculate safety factors
    massOwens,
    stress_U,
    SF_ult_U,
    SF_buck_U,
    stress_L,
    SF_ult_L,
    SF_buck_L,
    stress_TU,
    SF_ult_TU,
    SF_buck_TU,
    stress_TL,
    SF_ult_TL,
    SF_buck_TL,
    topstrainout_blade_U,
    topstrainout_blade_L,
    topstrainout_tower_U,
    topstrainout_tower_L,
    topDamage_blade_U,
    topDamage_blade_L,
    topDamage_tower_U,
    topDamage_tower_L = OWENS.extractSF(
        bld_precompinput,
        bld_precompoutput,
        plyprops_bld,
        numadIn_bld,
        lam_U_bld,
        lam_L_bld,
        twr_precompinput,
        twr_precompoutput,
        plyprops_twr,
        numadIn_twr,
        lam_U_twr,
        lam_L_twr,
        mymesh,
        myel,
        myort,
        setup_opts.blade.B,
        epsilon_x_hist,
        kappa_y_hist,
        kappa_z_hist,
        epsilon_z_hist,
        kappa_x_hist,
        epsilon_y_hist;
        verbosity = modelopt.OWENS_Options.verbosity,
        Twr_LE_U_idx = 1,
        Twr_LE_L_idx = 1,
        delta_t = modelopt.OWENS_Options.delta_t,
        AD15bldNdIdxRng,
        AD15bldElIdxRng,
        strut_precompoutput = nothing,
    )

    # Output data
    outputData(;
        mymesh,
        inputs,
        t,
        aziHist,
        OmegaHist,
        OmegaDotHist,
        gbHist,
        gbDotHist,
        gbDotDotHist,
        FReactionHist,
        genTorque,
        genPower,
        torqueDriveShaft,
        uHist,
        uHist_prp,
        epsilon_x_hist,
        epsilon_y_hist,
        epsilon_z_hist,
        kappa_x_hist,
        kappa_y_hist,
        kappa_z_hist,
        FTwrBsHist,
        massOwens,
        stress_U,
        SF_ult_U,
        SF_buck_U,
        stress_L,
        SF_ult_L,
        SF_buck_L,
        stress_TU,
        SF_ult_TU,
        SF_buck_TU,
        stress_TL,
        SF_ult_TL,
        SF_buck_TL,
        topstrainout_blade_U,
        topstrainout_blade_L,
        topstrainout_tower_U,
        topstrainout_tower_L,
        topDamage_blade_U,
        topDamage_blade_L,
        topDamage_tower_U,
        topDamage_tower_L,
    )
end

"""
    handle_control_strategy(controlStrategy, normalRPM, slowRPM, prescribed_time, prescribed_RPM)

Handles different control strategies and returns the appropriate control parameters.

# Arguments
- `controlStrategy`: The control strategy to use
- `normalRPM`: The normal operating RPM
- `slowRPM`: The slow RPM for startup/shutdown
- `prescribed_time`: Time control points for prescribed RPM
- `prescribed_RPM`: RPM control points for prescribed RPM

# Returns
- `tocp`: Time control points
- `Omegaocp`: Angular velocity control points
- `turbineStartup`: Turbine startup flag
- `generatorOn`: Generator on flag
- `useGeneratorFunction`: Use generator function flag
"""
function handle_control_strategy(
    controlStrategy,
    normalRPM,
    slowRPM,
    prescribed_time,
    prescribed_RPM,
)
    if controlStrategy == "normal"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [normalRPM, normalRPM, normalRPM, normalRPM] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "freewheelatNormalOperatingRPM"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [normalRPM, normalRPM, normalRPM, normalRPM] ./ 60 #hz
        turbineStartup = 2
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "startup"
        tocp = [0.0, 30.0, 1000.0]
        Omegaocp = [slowRPM, normalRPM, normalRPM] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "shutdown"
        tocp = [0.0, 10.0, 15.0, 30.0, 1000.0]
        Omegaocp = [normalRPM, normalRPM, normalRPM/2, slowRPM, 0.0] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "emergencyshutdown"
        tocp = [0.0, 10.0, 15.0, 30.0, 1000.0]
        Omegaocp = [normalRPM, normalRPM, slowRPM, 0.0, 0.0] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "parked"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [slowRPM, slowRPM, slowRPM, slowRPM] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "parked_idle"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [slowRPM, slowRPM, slowRPM, slowRPM] ./ 60 #hz
        turbineStartup = 2
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "parked_yaw"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [slowRPM, slowRPM, slowRPM, slowRPM] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "transport"
        tocp = [0.0, 10.0, 30.0, 1000.0]
        Omegaocp = [slowRPM, slowRPM, slowRPM, slowRPM] ./ 60 #hz
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    elseif controlStrategy == "prescribedRPM"
        tocp = prescribed_time
        Omegaocp = prescribed_RPM ./ 60
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    else
        @warn "ControlStrategy $controlStrategy not recognized, using prescribed RPM spline control points"
        tocp = prescribed_time
        Omegaocp = prescribed_RPM ./ 60
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
    end

    return tocp, Omegaocp, turbineStartup, generatorOn, useGeneratorFunction
end
