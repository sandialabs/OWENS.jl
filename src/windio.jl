
function runOWENSWINDIO(modelopt,windio,path)
    if isa(windio, String)
        windio = Design_Data(windio)
    end

    if isa(modelopt, String)
        modelopt = ModelingOptions(modelopt)
    end

    # Assembly
    # turbine_class = windio[:assembly][:turbine_class]
    # turbulence_class = windio[:assembly][:turbulence_class]
    # drivetrain = windio[:assembly][:drivetrain]
    # rotor_orientation = windio[:assembly][:rotor_orientation]
    number_of_blades = windio[:assembly][:number_of_blades] #Used
    # number_of_struts_per_blade = windio[:assembly][:number_of_struts_per_blade] #NEW?
    hub_height = windio[:assembly][:hub_height]
    # rotor_diameter = windio[:assembly][:rotor_diameter]
    # rated_power = windio[:assembly][:rated_power]
    # lifetime = windio[:assembly][:lifetime]
    # marine_hydro = windio[:assembly][:marine_hydro]
    
    # Components: blade, strut, tower, cable all reuse the blade format
    blade_mountpoint = windio[:components][:blade][:outer_shape_bem][:blade_mountpoint]
    # blade_x_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:grid]
    blade_x = windio[:components][:blade][:outer_shape_bem][:reference_axis][:x][:values] #Used

    # blade_y_grid = windio[:components][:blade][:outer_shape_bem][:reference_axis][:y][:grid]
    blade_y = windio[:components][:blade][:outer_shape_bem][:reference_axis][:y][:values] #Used

    # blade_z_grid = windio[:components][:blade][:outer_shape_bem][:reference_axis][:z][:grid]
    blade_z = windio[:components][:blade][:outer_shape_bem][:reference_axis][:z][:values] #Used

    tower_z = windio[:components][:tower][:outer_shape_bem][:reference_axis][:z][:values] #Used

    
    Blade_Height = maximum(blade_z) #TODO: resolve DLC dependence
    Blade_Radius = maximum(sqrt.(blade_x.^2 .+ blade_y.^2))
    
    Htwr_base = hub_height-Blade_Height/2
    Htwr_blds = maximum(tower_z)-Htwr_base

    # Struts
    #TODO: multiple struts
    tower_strut_connection = windio[:components][:struts][1][:mountfraction_tower]
    blade_strut_connection = windio[:components][:struts][1][:mountfraction_blade]

    # Airfoils
    #TODO

    # materials

    # Control
    # Vin = windio[:control][:supervisory][:Vin] 
    # Vout = windio[:control][:supervisory][:Vout] 
    # maxTS = windio[:control][:supervisory][:maxTS] 
 
    # PC_zeta = windio[:control][:pitch][:PC_zeta]
    # PC_omega = windio[:control][:pitch][:PC_omega]
    # ps_percent = windio[:control][:pitch][:ps_percent]
    # max_pitch = windio[:control][:pitch][:max_pitch]
    # max_pitch_rate = windio[:control][:pitch][:max_pitch_rate]
    # min_pitch = windio[:control][:pitch][:min_pitch]

    # control_type = windio[:control][:torque][:control_type]
    # tsr = windio[:control][:torque][:tsr]
    # VS_zeta = windio[:control][:torque][:VS_zeta]
    # VS_omega = windio[:control][:torque][:VS_omega]
    # max_torque_rate = windio[:control][:torque][:max_torque_rate]
    # VS_minspd = windio[:control][:torque][:VS_minspd]
    # VS_maxspd = windio[:control][:torque][:VS_maxspd]
    
    # ss_vsgain = windio[:control][:setpoint_smooth][:ss_vsgain]
    # ss_pcgain = windio[:control][:setpoint_smooth][:ss_pcgain]
    
    # limit_type = windio[:control][:shutdown][:limit_type]
    # limit_value = windio[:control][:shutdown][:limit_value]
    
    # Environment
    air_density = windio[:environment][:air_density] #used
    air_dyn_viscosity = windio[:environment][:air_dyn_viscosity] #used
    # air_speed_sound = windio[:environment][:air_speed_sound]
    # shear_exp = windio[:environment][:shear_exp]
    gravity = Float64.(windio[:environment][:gravity]) #used
    # weib_shape_parameter = windio[:environment][:weib_shape_parameter]
    # water_density = windio[:environment][:water_density]
    # water_dyn_viscosity = windio[:environment][:water_dyn_viscosity]
    # soil_shear_modulus = windio[:environment][:soil_shear_modulus]
    # soil_poisson = windio[:environment][:soil_poisson]
    # water_depth = windio[:environment][:water_depth]
    # air_pressure = windio[:environment][:air_pressure]
    # air_vapor_pressure = windio[:environment][:air_vapor_pressure]
    # significant_wave_height = windio[:environment][:significant_wave_height]
    # significant_wave_period = windio[:environment][:significant_wave_period]
    
    # BOS
    # plant_turbine_spacing = windio[:bos][:plant_turbine_spacing]
    # plant_row_spacing = windio[:bos][:plant_row_spacing]
    # commissioning_pct = windio[:bos][:commissioning_pct]
    # decommissioning_pct = windio[:bos][:decommissioning_pct]
    # distance_to_substation = windio[:bos][:distance_to_substation]
    # distance_to_interconnection = windio[:bos][:distance_to_interconnection]
    # interconnect_voltage = windio[:bos][:interconnect_voltage]
    # distance_to_site = windio[:bos][:distance_to_site]
    # distance_to_landfall = windio[:bos][:distance_to_landfall]
    # port_cost_per_month = windio[:bos][:port_cost_per_month]
    # site_auction_price = windio[:bos][:site_auction_price]
    # site_assessment_plan_cost = windio[:bos][:site_assessment_plan_cost]
    # site_assessment_cost = windio[:bos][:site_assessment_cost]
    # construction_operations_plan_cost = windio[:bos][:construction_operations_plan_cost]
    # boem_review_cost = windio[:bos][:boem_review_cost]
    # design_install_plan_cost = windio[:bos][:design_install_plan_cost]

    # Costs
    # wake_loss_factor = windio[:costs][:wake_loss_factor]
    # fixed_charge_rate = windio[:costs][:fixed_charge_rate]
    # bos_per_kW = windio[:costs][:bos_per_kW]
    # opex_per_kW = windio[:costs][:opex_per_kW]
    # turbine_number = windio[:costs][:turbine_number]
    # labor_rate = windio[:costs][:labor_rate]
    # painting_rate = windio[:costs][:painting_rate]
    # blade_mass_cost_coeff = windio[:costs][:blade_mass_cost_coeff]
    # hub_mass_cost_coeff = windio[:costs][:hub_mass_cost_coeff]
    # pitch_system_mass_cost_coeff = windio[:costs][:pitch_system_mass_cost_coeff]
    # spinner_mass_cost_coeff = windio[:costs][:spinner_mass_cost_coeff]
    # lss_mass_cost_coeff = windio[:costs][:lss_mass_cost_coeff]
    # bearing_mass_cost_coeff = windio[:costs][:bearing_mass_cost_coeff]
    # gearbox_mass_cost_coeff = windio[:costs][:gearbox_mass_cost_coeff]
    # hss_mass_cost_coeff = windio[:costs][:hss_mass_cost_coeff]
    # generator_mass_cost_coeff = windio[:costs][:generator_mass_cost_coeff]
    # bedplate_mass_cost_coeff = windio[:costs][:bedplate_mass_cost_coeff]
    # yaw_mass_cost_coeff = windio[:costs][:yaw_mass_cost_coeff]
    # converter_mass_cost_coeff = windio[:costs][:converter_mass_cost_coeff]
    # transformer_mass_cost_coeff = windio[:costs][:transformer_mass_cost_coeff]
    # hvac_mass_cost_coeff = windio[:costs][:hvac_mass_cost_coeff]
    # cover_mass_cost_coeff = windio[:costs][:cover_mass_cost_coeff]
    # elec_connec_machine_rating_cost_coeff = windio[:costs][:elec_connec_machine_rating_cost_coeff]
    # platforms_mass_cost_coeff = windio[:costs][:platforms_mass_cost_coeff]
    # tower_mass_cost_coeff = windio[:costs][:tower_mass_cost_coeff]
    # controls_machine_rating_cost_coeff = windio[:costs][:controls_machine_rating_cost_coeff]
    # crane_cost = windio[:costs][:crane_cost]
    # electricity_price = windio[:costs][:electricity_price]
    # reserve_margin_price = windio[:costs][:reserve_margin_price]
    # capacity_credit = windio[:costs][:capacity_credit]
    # benchmark_price = windio[:costs][:benchmark_price]

    if isa(gravity, Float64)
        gravityOn = [0,0,gravity]
    else
        gravityOn = gravity
    end

    # Top Level OWENS Options
    analysisType = modelopt.OWENS_Options.analysisType
    AeroModel = modelopt.OWENS_Options.AeroModel
    controlStrategy = modelopt.OWENS_Options.controlStrategy
    numTS = modelopt.OWENS_Options.numTS
    delta_t = modelopt.OWENS_Options.delta_t
    platformActive = modelopt.OWENS_Options.platformActive
    topsideOn = modelopt.OWENS_Options.topsideOn
    interpOrder = modelopt.OWENS_Options.interpOrder
    dataOutputFilename = modelopt.OWENS_Options.dataOutputFilename
    rigid = modelopt.OWENS_Options.rigid
    TOL = modelopt.OWENS_Options.TOL
    MAXITER = modelopt.OWENS_Options.MAXITER
    verbosity = modelopt.OWENS_Options.verbosity
    VTKsaveName = modelopt.OWENS_Options.VTKsaveName
    aeroLoadsOn = modelopt.OWENS_Options.aeroLoadsOn
    structuralModel = modelopt.OWENS_Options.structuralModel
    Prescribed_RPM_time_controlpoints = modelopt.OWENS_Options.Prescribed_RPM_time_controlpoints
    Prescribed_RPM_RPM_controlpoints = modelopt.OWENS_Options.Prescribed_RPM_RPM_controlpoints
    Prescribed_Vinf_time_controlpoints = modelopt.OWENS_Options.Prescribed_Vinf_time_controlpoints
    Prescribed_Vinf_Vinf_controlpoints = modelopt.OWENS_Options.Prescribed_Vinf_Vinf_controlpoints

    # OWENSAero Options
    Nslices = modelopt.OWENSAero_Options.Nslices
    ntheta = modelopt.OWENSAero_Options.ntheta
    ifw = modelopt.OWENSAero_Options.ifw
    DynamicStallModel = modelopt.OWENSAero_Options.DynamicStallModel
    RPI = modelopt.OWENSAero_Options.RPI
    Aero_Buoyancy_Active = modelopt.OWENSAero_Options.Aero_Buoyancy_Active
    Aero_AddedMass_Active = modelopt.OWENSAero_Options.Aero_AddedMass_Active
    Aero_RotAccel_Active = modelopt.OWENSAero_Options.Aero_RotAccel_Active

    # DLC Options
    DLCs = modelopt.DLC_Options.DLCs
    Vinf_range = modelopt.DLC_Options.Vinf_range
    IEC_std = modelopt.DLC_Options.IEC_std
    WindChar = modelopt.DLC_Options.WindChar
    WindClass = modelopt.DLC_Options.WindClass
    turbsimsavepath = modelopt.DLC_Options.turbsimsavepath
    pathtoturbsim = modelopt.DLC_Options.pathtoturbsim
    NumGrid_Z = modelopt.DLC_Options.NumGrid_Z
    NumGrid_Y = modelopt.DLC_Options.NumGrid_Y
    Vref = modelopt.DLC_Options.Vref
    Vdesign = modelopt.DLC_Options.Vdesign
    grid_oversize = modelopt.DLC_Options.grid_oversize
    regenWindFiles = modelopt.DLC_Options.regenWindFiles
    delta_t_turbsim = modelopt.DLC_Options.delta_t_turbsim
    simtime_turbsim = modelopt.DLC_Options.simtime_turbsim

    # OWENSFEA Options
    nlOn = modelopt.OWENSFEA_Options.nlOn
    RayleighAlpha = modelopt.OWENSFEA_Options.RayleighAlpha
    RayleighBeta = modelopt.OWENSFEA_Options.RayleighBeta
    iterationType = modelopt.OWENSFEA_Options.iterationType
    guessFreq = modelopt.OWENSFEA_Options.guessFreq
    numModes = modelopt.OWENSFEA_Options.numModes
    adaptiveLoadSteppingFlag = modelopt.OWENSFEA_Options.adaptiveLoadSteppingFlag
    minLoadStepDelta = modelopt.OWENSFEA_Options.minLoadStepDelta
    minLoadStep = modelopt.OWENSFEA_Options.minLoadStep
    prescribedLoadStep = modelopt.OWENSFEA_Options.prescribedLoadStep
    maxNumLoadSteps = modelopt.OWENSFEA_Options.maxNumLoadSteps
    tolerance = modelopt.OWENSFEA_Options.tolerance
    maxIterations = modelopt.OWENSFEA_Options.maxIterations
    elementOrder = modelopt.OWENSFEA_Options.elementOrder
    alpha = modelopt.OWENSFEA_Options.alpha
    gamma = modelopt.OWENSFEA_Options.gamma
    AddedMass_Coeff_Ca = modelopt.OWENSFEA_Options.AddedMass_Coeff_Ca
    platformTurbineConnectionNodeNumber = modelopt.OWENSFEA_Options.platformTurbineConnectionNodeNumber
    aeroElasticOn = modelopt.OWENSFEA_Options.aeroElasticOn
    spinUpOn = modelopt.OWENSFEA_Options.spinUpOn
    predef = modelopt.OWENSFEA_Options.predef

    # OWENSOpenFASTWrappers Options
    windINPfilename = modelopt.OWENSOpenFASTWrappers_Options.windINPfilename
    ifw_libfile = modelopt.OWENSOpenFASTWrappers_Options.ifw_libfile
    hd_lib = modelopt.OWENSOpenFASTWrappers_Options.hd_lib
    md_lib = modelopt.OWENSOpenFASTWrappers_Options.md_lib
    adi_lib = modelopt.OWENSOpenFASTWrappers_Options.adi_lib
    adi_rootname = modelopt.OWENSOpenFASTWrappers_Options.adi_rootname
    hd_input_file = modelopt.OWENSOpenFASTWrappers_Options.hd_input_file
    ss_input_file = modelopt.OWENSOpenFASTWrappers_Options.ss_input_file
    md_input_file = modelopt.OWENSOpenFASTWrappers_Options.md_input_file
    potflowfile = modelopt.OWENSOpenFASTWrappers_Options.potflowfile
    WindType = modelopt.OWENSOpenFASTWrappers_Options.WindType
    
    if adi_lib == "nothing"
        adi_lib = nothing
    end

    windINPfilename = "$(path)$windINPfilename"
    potflowfile = "$(path)$potflowfile"
    if ifw_libfile == "nothing"
        ifw_libfile = nothing
    end

    # Mesh Options
    ntelem = modelopt.Mesh_Options.ntelem
    nbelem = modelopt.Mesh_Options.nbelem
    ncelem = modelopt.Mesh_Options.ncelem
    nselem = modelopt.Mesh_Options.nselem
    angularOffset = modelopt.Mesh_Options.angularOffset
    joint_type = modelopt.Mesh_Options.joint_type
    c_mount_ratio = modelopt.Mesh_Options.c_mount_ratio
    AD15hubR = modelopt.Mesh_Options.AD15hubR
    cables_connected_to_blade_base = modelopt.Mesh_Options.cables_connected_to_blade_base
    turbineType = modelopt.Mesh_Options.turbineType

    # Drivetrain Options
    turbineStartup = modelopt.Drivetrain_Options.turbineStartup
    usingRotorSpeedFunction = modelopt.Drivetrain_Options.usingRotorSpeedFunction
    driveTrainOn = modelopt.Drivetrain_Options.driveTrainOn
    JgearBox = modelopt.Drivetrain_Options.JgearBox
    gearRatio = modelopt.Drivetrain_Options.gearRatio
    gearBoxEfficiency = modelopt.Drivetrain_Options.gearBoxEfficiency
    generatorOn = modelopt.Drivetrain_Options.generatorOn
    useGeneratorFunction = modelopt.Drivetrain_Options.useGeneratorFunction
    generatorProps = modelopt.Drivetrain_Options.generatorProps
    ratedTorque = modelopt.Drivetrain_Options.ratedTorque
    zeroTorqueGenSpeed = modelopt.Drivetrain_Options.zeroTorqueGenSpeed
    pulloutRatio = modelopt.Drivetrain_Options.pulloutRatio
    ratedGenSlipPerc = modelopt.Drivetrain_Options.ratedGenSlipPerc
    OmegaGenStart = modelopt.Drivetrain_Options.OmegaGenStart
    driveShaft_K = modelopt.Drivetrain_Options.driveShaft_K
    driveShaft_C = modelopt.Drivetrain_Options.driveShaft_C

    strut_twr_mountpoint = tower_strut_connection
    strut_bld_mountpoint = blade_strut_connection

    # Inputs that are part of the overall options, but which are not yet available at the top level yaml input method
    nlParams = 0 # derived struct, we aren't going to pass in the nlParams struct, but rather use the detailed inputs above, so hard code here.
    bladeData = [] # same as above
    numDOFPerNode = 6 #while much of the model can operate with fewer dofs, too much is hard coded on the full 6 dof.
    omegaControl = false # this is a derived parameter
    meshtype = turbineType #derived, should probably be cleaned up TODO
    initCond = [] #OWENSFEA initial conditions array, will be derived at this level, if using the OWENS scripting method, this can be used to initalize the structure
    jointTransform = 0.0 #OWENSFEA, derived matrix to transform from total matrix to the reduced matrix and vice-versa
    reducedDOFList = zeros(Int,2) #OWENSFEA, derived array that maps the joint-reduced dofs
    nodalTerms = 0.0 #OWENSFEA the ability to apply concentrated nodal masses and forces etc., currently only available at the scripting level of analysis
    stack_layers_bld = nothing #, enables direct specification of the numbers of stack layers in the numad format
    stack_layers_scale = [1.0,1.0] #, simple scaling across the blade span with a linear interpolation between
    chord_scale = [1.0,1.0] #, simple scaling across the blade span with a linear interpolation between
    thickness_scale = [1.0,1.0] #, simple scaling across the blade span with a linear interpolation between

    OmegaInit = 7.2/60 #TODO: simplify this in the code since it is redundant
    aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv"
    owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens"
    
    
    driveShaftProps = DriveShaftProps(driveShaft_K,driveShaft_C)

    custommesh = nothing

    tsave_idx=1:3:numTS-1

    NuMad_geom_xlscsv_file_twr = windio
    NuMad_mat_xlscsv_file_twr = windio
    NuMad_geom_xlscsv_file_bld = windio
    NuMad_mat_xlscsv_file_bld = windio
    NuMad_geom_xlscsv_file_strut = windio
    NuMad_mat_xlscsv_file_strut = windio

    adi_rootname = "$(path)$(adi_rootname)"

    shapeZ = blade_z#collect(LinRange(0,H,Nslices+1))
    shapeX = blade_x#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)
    shapeY = blade_y

    R = maximum(blade_x) #m 
    H = maximum(blade_z) #m
    
    nothing

    # Call the helper function that builds the mesh, calculates the sectional properties,
    # and aligns the sectional properties to the mesh elements, 

    mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
    mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
        rho=air_density,
        mu=air_dyn_viscosity,
        Nslices,
        ntheta,
        RPM=Prescribed_RPM_RPM_controlpoints[1],
        Vinf=Prescribed_Vinf_Vinf_controlpoints[1],
        eta=blade_mountpoint,
        B=number_of_blades,
        H,
        R,
        AD15hubR,
        shapeZ, 
        shapeX,
        shapeY, 
        ifw,
        WindType,
        delta_t,
        numTS,
        adi_lib,
        adi_rootname,
        windINPfilename,
        ifw_libfile,
        NuMad_geom_xlscsv_file_twr,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS_Cables.csv",
        NuMad_mat_xlscsv_file_twr,# = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv",
        NuMad_geom_xlscsv_file_bld,# = "$path/data/NuMAD_Geom_SNL_5MW_ARCUS.csv",
        NuMad_mat_xlscsv_file_bld,# = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
        NuMad_geom_xlscsv_file_strut,
        NuMad_mat_xlscsv_file_strut,
        Htwr_base,
        ntelem, 
        nbelem, 
        ncelem,
        nselem,
        joint_type,
        c_mount_ratio,
        strut_twr_mountpoint,
        strut_bld_mountpoint,
        AeroModel, #AD, DMS, AC
        DynamicStallModel,
        RPI,
        cables_connected_to_blade_base,
        meshtype,
        Htwr_blds,
        stack_layers_bld,
        stack_layers_scale,
        chord_scale,
        thickness_scale,
        angularOffset,
        Aero_AddedMass_Active,
        Aero_RotAccel_Active,
        Aero_Buoyancy_Active,
        custommesh)

    nothing

    # Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
    # this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
    # includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
    # to work with OWENS inputs as well.

    nothing

    # Here we apply the boundary conditions.  For this case, with a regular cantelever tower, the tower base node which is 
    # 1 is constrained in all 6 degrees of freedom to have a displacement of 0.  You can change this displacement to allow for things
    # like pretension, and you can apply boundary conditions to any node.

    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0]

    nothing

    # There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

     # Handle the control strategy #TODO: function and hook up discon controller
    if DLCs != "none"
        normalRPM = Prescribed_RPM_RPM_controlpoints[1] ./ 60
        slowRPM = 1.0 #RPM
        tocp = [0.0,10.0,30.0,1000.0]
        turbineStartup = 0
        generatorOn = false
        useGeneratorFunction = false
        throwawayTimeSteps=1#round(Int,10.0/delta_t) # 10 seconds
        if controlStrategy == "normal"
            Omegaocp = [normalRPM,normalRPM,normalRPM,normalRPM]./60 #hz
            turbineStartup = 0
            generatorOn = false
            useGeneratorFunction = false
        elseif controlStrategy == "freewheelatNormalOperatingRPM"
            Omegaocp = [normalRPM,normalRPM,normalRPM,normalRPM]./60 #hz
            turbineStartup = 2
            generatorOn = false
            useGeneratorFunction = false
        elseif controlStrategy == "startup"
            throwawayTimeSteps = 1
            tocp = [0.0,30.0,1000.0]
            Omegaocp = [slowRPM,normalRPM,normalRPM]./60 #hz
            turbineStartup = 0
            generatorOn = false
            useGeneratorFunction = false
        elseif controlStrategy == "shutdown"
            tocp = [0.0,10.0,15.0,30.0,1000.0]
            Omegaocp = [normalRPM,normalRPM,normalRPM/2,slowRPM,0.0]./60 #hz
        elseif controlStrategy == "emergencyshutdown"
            tocp = [0.0,10.0,15.0,30.0,1000.0]
            Omegaocp = [normalRPM,normalRPM,slowRPM,0.0,0.0]./60 #hz
        elseif controlStrategy == "parked"
            Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
        elseif controlStrategy == "parked_idle"
            Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
            turbineStartup = 2
            generatorOn = false
            useGeneratorFunction = false
        elseif controlStrategy == "parked_yaw"
            Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
        elseif controlStrategy == "parked"
            Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
        elseif controlStrategy == "transport"
            Omegaocp = [slowRPM,slowRPM,slowRPM,slowRPM]./60 #hz
        elseif controlStrategy == "prescribedRPM"
            tocp = Prescribed_RPM_time_controlpoints
            Omegaocp = Prescribed_RPM_RPM_controlpoints ./ 60
        else
            @warn "ControlStrategy $controlStrategy not recognized, using prescribed RPM spline control points"
            tocp = Prescribed_RPM_time_controlpoints
            Omegaocp = Prescribed_RPM_RPM_controlpoints ./ 60
        end
    else
        tocp = Prescribed_RPM_time_controlpoints
        Omegaocp = Prescribed_RPM_RPM_controlpoints ./ 60

    end


    println("controlStrategy: $controlStrategy")

    inputs = OWENS.Inputs(;analysisType = structuralModel,
    tocp,
    Omegaocp,
    tocp_Vinf = Prescribed_Vinf_time_controlpoints,
    Vinfocp = Prescribed_Vinf_Vinf_controlpoints,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn,
    turbineStartup,
    usingRotorSpeedFunction,
    driveTrainOn,
    generatorOn,
    platformActive,
    topsideOn,
    interpOrder,
    hd_input_file,
    ss_input_file,
    md_input_file,
    JgearBox,
    gearRatio,
    gearBoxEfficiency,
    useGeneratorFunction,
    generatorProps,
    ratedTorque,
    zeroTorqueGenSpeed,
    pulloutRatio,
    ratedGenSlipPerc,
    OmegaGenStart,
    omegaControl,
    OmegaInit,
    aeroloadfile,
    owensfile,
    potflowfile,
    dataOutputFilename,
    numDOFPerNode,
    bladeData,
    rigid,
    driveShaftProps,
    TOL,
    MAXITER)

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

    feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    airDensity=air_density,
    pBC,
    nlOn,
    gravityOn,
    numNodes = mymesh.numNodes,
    RayleighAlpha,
    RayleighBeta,
    joint = myjoint,
    iterationType,
    initCond,
    aeroElasticOn, #for the automated flutter model
    guessFreq,
    dataOutputFilename,
    jointTransform,
    reducedDOFList,
    numDOFPerNode,
    platformTurbineConnectionNodeNumber,
    spinUpOn,
    numModes,
    nlParams, # can pass in strut, or leave at 0 to use other inputs
    adaptiveLoadSteppingFlag,
    tolerance,
    maxIterations,
    maxNumLoadSteps,
    minLoadStepDelta,
    minLoadStep,
    prescribedLoadStep,
    predef,
    elementOrder,
    alpha,
    gamma,
    AddedMass_Coeff_Ca,
    nodalTerms)

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time.

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
    topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

    nothing

    # Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

    OWENS.OWENSVTK(VTKsaveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
        epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
        FReactionHist,topFexternal_hist;tsave_idx)

    nothing

    # This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
    # it also calculates analytical simply supported buckling safety factors

    ##########################################
    #### Ultimate Failure #####
    ##########################################

    massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
    SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
    topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
    topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,number_of_blades,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,delta_t,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

    outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)
end