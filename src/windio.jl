
function runOWENSWINDIO(windio,modelopt,path)
    if typeof(windio) == String
        windio = YAML.load_file(windio; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
        println("Running: $(windio[:name])")
    end

    if typeof(modelopt) == String
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

    
    Blade_Height = maximum(blade_z) #TODO: resolve DLC dependence
    Blade_Radius = maximum(sqrt.(blade_x.^2 .+ blade_y.^2))
    
    Htwr_base = hub_height-Blade_Height/2
    if turbineType == "Darrieus"
        Htwr_blds = Blade_Height
    else
        Htwr_blds = Blade_Height*0.6 #TODO: finer grained inputs
    end

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

    if typeof(gravity) == Float64
        gravityOn = [0,0,gravity]
    else
        gravityOn = gravity
    end


    analysisType = modelopt.analysisType
    turbineType = modelopt.turbineType
    Vinf = modelopt.Vinf
    controlStrategy = modelopt.controlStrategy
    RPM = modelopt.RPM
    Nslices = modelopt.Nslices
    ntheta = modelopt.ntheta
    structuralModel = modelopt.structuralModel
    ntelem = modelopt.ntelem
    nbelem = modelopt.nbelem
    ncelem = modelopt.ncelem
    nselem = modelopt.nselem
    ifw = modelopt.ifw
    WindType = modelopt.WindType
    AeroModel = modelopt.AeroModel
    windINPfilename = "$(path)$(modelopt.windINPfilename)"
    ifw_libfile = modelopt.ifw_libfile
    if ifw_libfile == "nothing"
        ifw_libfile = nothing
    end
    numTS = modelopt.numTS
    delta_t = modelopt.delta_t

    turbineStartup = 0
    usingRotorSpeedFunction = false
    driveTrainOn = false
    generatorOn = false
    platformActive = false
    topsideOn = true
    interpOrder = 2
    hd_input_file = "none"
    ss_input_file = "none"
    md_input_file = "none"
    JgearBox = 0.0
    gearRatio = 1.0
    gearBoxEfficiency = 1.0
    useGeneratorFunction = false
    generatorProps = 0.0
    ratedTorque = 0.0
    zeroTorqueGenSpeed = 0.0
    pulloutRatio = 0.0
    ratedGenSlipPerc = 0.0
    OmegaGenStart = 0.0
    omegaControl = false
    OmegaInit = 7.2/60 #TODO: simplify this in the code since it is redundant
    aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv"
    owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens"
    potflowfile = "$module_path/../test/data/potential_flow_data"
    dataOutputFilename = "none"
    numDofPerNode = 6
    bladeData = []
    rigid = false
    driveShaftProps = DriveShaftProps(0.0,0.0)
    TOl = 1e-4
    MAXITER = 300
    verbosity = 2
    joint_type = 0
    c_mount_ratio = 0.05
    strut_twr_mountpoint = tower_strut_connection #TODO: multiple struts
    strut_bld_mountpoint = blade_strut_connection
    DynamicStallModel="BV"
    RPI=true
    cables_connected_to_blade_base = true
    meshtype = turbineType
    VTKsaveName = "$path/vtk/windio"
    aeroLoadsOn = 2
    nlOn = true
    RayleighAlpha = 0.05
    RayleighBeta = 0.05
    iterationType = "DI"
    initCond = []
    aeroElasticOn = false #for the automated flutter model
    guessFreq = 0.0
    dataOutputFilename = "none"
    jointTransform = 0.0
    reducedDOFList = zeros(Int,2)
    numDOFPerNode = 6
    platformTurbineConnectionNodeNumber = 1
    spinUpOn = false
    numModes = 20
    nlParams = 0 # can pass in strut, or leave at 0 to use other inputs
    adaptiveLoadSteppingFlag = true
    tolerance = 1.0000e-06
    maxIterations = 50
    maxNumLoadSteps = 20
    minLoadStepDelta = 0.0500
    minLoadStep = 0.0500
    prescribedLoadStep = 0.0
    predef = false
    elementOrder = 1
    alpha = 0.5
    gamma = 0.5
    nodalTerms = 0.0

    AD15hubR = 0.1
    # Htwr_base = 3.0#maximum(windio[:components][:tower][:outer_shape_bem][:reference_axis][:z][:values])
    stack_layers_bld = nothing
    stack_layers_scale = [1.0,1.0]
    chord_scale = [1.0,1.0]
    thickness_scale = [1.0,1.0]
    angularOffset = -pi/2
    Aero_AddedMass_Active = false
    Aero_RotAccel_Active = false
    Aero_Buoyancy_Active = false

    custommesh = nothing

    tsave_idx=1:3:numTS

    NuMad_geom_xlscsv_file_twr = windio #"$(path)$(modelopt.NuMad_geom_xlscsv_file_twr)"
    NuMad_mat_xlscsv_file_twr = windio #"$(path)$(modelopt.NuMad_mat_xlscsv_file_twr)"
    NuMad_geom_xlscsv_file_bld = windio #"$(path)$(modelopt.NuMad_geom_xlscsv_file_bld)"
    NuMad_mat_xlscsv_file_bld = windio #"$(path)$(modelopt.NuMad_mat_xlscsv_file_bld)"
    NuMad_geom_xlscsv_file_strut = windio #"$(path)$(modelopt.NuMad_geom_xlscsv_file_strut)"
    NuMad_mat_xlscsv_file_strut = windio #"$(path)$(modelopt.NuMad_mat_xlscsv_file_strut)"
    adi_lib = modelopt.adi_lib
    if adi_lib == "nothing"
        adi_lib = nothing
    end
    adi_rootname = "$(path)$(modelopt.adi_rootname)"

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
        RPM,
        Vinf,
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

    # If the sectional properties material files includes cost information, that is combined with the density 
    # to estimate the overall material cost of of materials in the blades

    if verbosity>0
        
        println("\nBlades' Mass Breakout")
        for (i,name) in enumerate(plyprops_bld.names)
            println("$name $(mass_breakout_blds[i]) kg, $(plyprops_bld.costs[i]) \$/kg: \$$(mass_breakout_blds[i]*plyprops_bld.costs[i])")
        end
        
        println("\nTower Mass Breakout")
        for (i,name) in enumerate(plyprops_twr.names)
            println("$name $(mass_breakout_twr[i]) kg, $(plyprops_twr.costs[i]) \$/kg: \$$(mass_breakout_twr[i]*plyprops_twr.costs[i])")
        end
        
        println("Total Material Cost Blades: \$$(sum(mass_breakout_blds.*plyprops_bld.costs))")
        println("Total Material Cost Tower: \$$(sum(mass_breakout_twr.*plyprops_twr.costs))")
        println("Total Material Cost: \$$(sum(mass_breakout_blds.*plyprops_bld.costs)+ sum(mass_breakout_twr.*plyprops_twr.costs))")
        
    end

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

    inputs = OWENS.Inputs(;analysisType = structuralModel,
    tocp = [0.0,100000.1],
    Omegaocp = [RPM,RPM] ./ 60,
    tocp_Vinf = [0.0,100000.1],
    Vinfocp = [Vinf,Vinf],
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
    numDofPerNode,
    bladeData,
    rigid,
    driveShaftProps,
    TOl,
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

    azi=aziHist#./aziHist*1e-6
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
    LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
    LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

    outputData(;mymesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FTwrBsHist,massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,topDamage_blade_L,topDamage_tower_U,topDamage_tower_L)
end