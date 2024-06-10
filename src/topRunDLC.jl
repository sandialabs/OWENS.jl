
mutable struct MasterInput
    analysisType
    turbineType
    eta
    Nbld #WindIO
    towerHeight
    rho
    Vinf
    controlStrategy
    RPM
    Nslices
    ntheta
    structuralModel
    ntelem
    nbelem
    ncelem
    nselem
    AModel
    ifw
    WindType
    windINPfilename
    ifw_libfile
    adi_lib
    adi_rootname
    Blade_Height
    Blade_Radius
    numTS
    delta_t
    NuMad_geom_xlscsv_file_twr
    NuMad_mat_xlscsv_file_twr
    NuMad_geom_xlscsv_file_bld
    NuMad_mat_xlscsv_file_bld
    NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut
end

function MasterInput(;
    analysisType =  "unsteady", # unsteady, steady, modal
    turbineType =  "Darrieus", #Darrieus, H-VAWT, ARCUS
    eta =  0.5, # blade mount point ratio, 0.5 is the blade half chord is perpendicular with the axis of rotation, 0.25 is the quarter chord, etc
    Nbld =  3, # number of blade
    Blade_Height = 54.01123056,
    Blade_Radius = 110.1829092,
    towerHeight =  3.0, # m tower extension height below blades
    rho =  1.225, # air density
    Vinf =  17.2, # m/s
    controlStrategy = "constantRPM", # TODO: incorporate the others
    RPM =  17.2, #RPM
    Nslices =  30, # number of VAWTAero discritizations 
    ntheta =  30, # number of VAWTAero azimuthal discretizations
    structuralModel = "GX", #GX, TNB, ROM
    ntelem =  10, #tower elements in each 
    nbelem =  60, #blade elements in each 
    ncelem =  10, #central cable elements in each if turbineType is ARCUS
    nselem =  5, #strut elements in each if turbineType has struts
    AModel = "AD",
    ifw = false,
    WindType = 1,
    ifw_libfile = "./../openfast/build/modules/inflowwind/libifw_c_binding",
    adi_lib = "./../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    numTS = 100,
    delta_t = 0.01,
    windINPfilename ="$module_path/../test/data/turbsim/115mx115m_30x30_20.0msETM.bts",
    NuMad_geom_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv",
    NuMad_mat_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
    NuMad_mat_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_Struts.csv",
    NuMad_mat_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv"
    )

    return MasterInput(analysisType,turbineType,eta,Nbld,towerHeight,rho,Vinf,controlStrategy,
    RPM,Nslices,ntheta,structuralModel,ntelem,nbelem,ncelem,nselem,AModel,ifw,WindType,windINPfilename,ifw_libfile,adi_lib,adi_rootname,
    Blade_Height,Blade_Radius,numTS,delta_t,NuMad_geom_xlscsv_file_twr,NuMad_mat_xlscsv_file_twr,
    NuMad_geom_xlscsv_file_bld,NuMad_mat_xlscsv_file_bld,NuMad_geom_xlscsv_file_strut,NuMad_mat_xlscsv_file_strut)
end

function MasterInput(yamlInputfile)
    yamlInput = YAML.load_file(yamlInputfile)
    # Unpack YAML
    general = yamlInput["general"]
        analysisType = general["analysisType"]
        turbineType = general["turbineType"]

    designParameters = yamlInput["designParameters"]
        eta = designParameters["eta"]
        Nbld = designParameters["Nbld"]
        Blade_Height = designParameters["Blade_Height"]
        Blade_Radius = designParameters["Blade_Radius"]
        towerHeight = designParameters["towerHeight"]

    operationParameters = yamlInput["operationParameters"]
        rho = operationParameters["rho"]
        Vinf = operationParameters["Vinf"]

    controlParameters = yamlInput["controlParameters"]
        controlStrategy = controlParameters["controlStrategy"]
        RPM = controlParameters["RPM"]
        numTS = controlParameters["numTS"]
        delta_t = controlParameters["delta_t"]

    AeroParameters = yamlInput["AeroParameters"]
        Nslices = AeroParameters["Nslices"]
        ntheta = AeroParameters["ntheta"]
        AModel = AeroParameters["AModel"]
        adi_lib = AeroParameters["adi_lib"]
        adi_rootname = AeroParameters["adi_rootname"]

    turbulentInflow = yamlInput["turbulentInflow"]
        ifw = turbulentInflow["ifw"]
        WindType = turbulentInflow["WindType"]
        windINPfilename = turbulentInflow["windINPfilename"]
        ifw_libfile = turbulentInflow["ifw_libfile"]

    structuralParameters = yamlInput["structuralParameters"]
        structuralModel = structuralParameters["structuralModel"]
        ntelem = structuralParameters["ntelem"]
        nbelem = structuralParameters["nbelem"]
        ncelem = structuralParameters["ncelem"]
        nselem = structuralParameters["nselem"]
        NuMad_geom_xlscsv_file_twr = structuralParameters["NuMad_geom_xlscsv_file_twr"]
        NuMad_mat_xlscsv_file_twr = structuralParameters["NuMad_mat_xlscsv_file_twr"]
        NuMad_geom_xlscsv_file_bld = structuralParameters["NuMad_geom_xlscsv_file_bld"]
        NuMad_mat_xlscsv_file_bld = structuralParameters["NuMad_mat_xlscsv_file_bld"]
        NuMad_geom_xlscsv_file_strut = structuralParameters["NuMad_geom_xlscsv_file_strut"]
        NuMad_mat_xlscsv_file_strut = structuralParameters["NuMad_mat_xlscsv_file_strut"]

    return MasterInput(analysisType,turbineType,eta,Nbld,towerHeight,rho,Vinf,
    controlStrategy,RPM,Nslices,ntheta,structuralModel,ntelem,nbelem,ncelem,
    nselem,AModel,ifw,WindType,windINPfilename,ifw_libfile,adi_lib,adi_rootname,Blade_Height,Blade_Radius,numTS,
    delta_t,NuMad_geom_xlscsv_file_twr,NuMad_mat_xlscsv_file_twr,
    NuMad_geom_xlscsv_file_bld,NuMad_mat_xlscsv_file_bld,NuMad_geom_xlscsv_file_strut,NuMad_mat_xlscsv_file_strut)
end

function runOWENS(Inp,path;verbosity=2)
    analysisType = Inp.analysisType
    turbineType = Inp.turbineType
    eta = Inp.eta
    Nbld = Inp.Nbld
    towerHeight = Inp.towerHeight
    rho = Inp.rho
    Vinf = Inp.Vinf
    controlStrategy = Inp.controlStrategy
    RPM = Inp.RPM
    Nslices = Inp.Nslices
    ntheta = Inp.ntheta
    structuralModel = Inp.structuralModel
    ntelem = Inp.ntelem
    nbelem = Inp.nbelem
    ncelem = Inp.ncelem
    nselem = Inp.nselem
    ifw = Inp.ifw
    WindType = Inp.WindType
    AModel = Inp.AModel
    windINPfilename = "$(path)$(Inp.windINPfilename)"
    ifw_libfile = Inp.ifw_libfile
    if ifw_libfile == "nothing"
        ifw_libfile = nothing
    end
    Blade_Height = Inp.Blade_Height
    Blade_Radius = Inp.Blade_Radius
    numTS = Inp.numTS
    delta_t = Inp.delta_t
    NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
    NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
    NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld)"
    NuMad_mat_xlscsv_file_bld = "$(path)$(Inp.NuMad_mat_xlscsv_file_bld)"
    NuMad_geom_xlscsv_file_strut = "$(path)$(Inp.NuMad_geom_xlscsv_file_strut)"
    NuMad_mat_xlscsv_file_strut = "$(path)$(Inp.NuMad_mat_xlscsv_file_strut)"
    adi_lib = Inp.adi_lib
    if adi_lib == "nothing"
        adi_lib = nothing
    end
    adi_rootname = "$(path)$(Inp.adi_rootname)"

    B = Nbld
    R = Blade_Radius#177.2022*0.3048 #m
    H = Blade_Height#1.02*R*2 #m

    shapeY = collect(LinRange(0,H,Nslices+1))
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)

    nothing

    # Call the helper function that builds the mesh, calculates the sectional properties,
    # and aligns the sectional properties to the mesh elements, 


mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    RPM,
    Vinf,
    eta,
    B,
    H,
    R,
    shapeY,
    shapeX,
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
    Ht=towerHeight,
    ntelem, 
    nbelem, 
    ncelem,
    nselem,
    joint_type = 0,
    c_mount_ratio = 0.05,
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = turbineType)

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

    if AModel=="AD"
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
    aeroLoadsOn = 2)

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

    feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    outFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = true,
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time.

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

    nothing

    # Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

    println("Saving VTK time domain files")
    OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

    nothing

    # This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
    # it also calculates analytical simply supported buckling safety factors

    ##########################################
    #### Ultimate Failure #####
    ##########################################

    massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
    SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
    topstrainout_tower_U,topstrainout_tower_LtopDamage_blade_U,
    topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
    LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
    LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,
    AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

    return [1.0,2.0,3.0]

end


function runOWENSWINDIO(WINDIO_filename,Inp,path;verbosity=2)

    windio = YAML.load_file(WINDIO_filename; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
    println("Running: $(windio[:name])")

    # Assembly
    # turbine_class = windio[:assembly][:turbine_class]
    # turbulence_class = windio[:assembly][:turbulence_class]
    # drivetrain = windio[:assembly][:drivetrain]
    # rotor_orientation = windio[:assembly][:rotor_orientation]
    number_of_blades = windio[:assembly][:number_of_blades] #Used
    # hub_height = windio[:assembly][:hub_height]
    # rotor_diameter = windio[:assembly][:rotor_diameter]
    # rated_power = windio[:assembly][:rated_power]
    # lifetime = windio[:assembly][:lifetime]
    # marine_hydro = windio[:assembly][:marine_hydro]
    
    # Components
    # blade_x_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:grid]
    blade_x = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:x][:values] #Used

    # blade_y_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:y][:grid]
    blade_y = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:y][:values] #Used

    # blade_z_grid = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:z][:grid]
    blade_z = windio[:components][:blade][:internal_structure_2d_fem][:reference_axis][:z][:values] #Used

    # Airfoils

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
    gravity = windio[:environment][:gravity] #used
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



    analysisType = Inp.analysisType
    turbineType = Inp.turbineType
    eta = Inp.eta
    towerHeight = Inp.towerHeight
    Vinf = Inp.Vinf
    controlStrategy = Inp.controlStrategy
    RPM = Inp.RPM
    Nslices = Inp.Nslices
    ntheta = Inp.ntheta
    structuralModel = Inp.structuralModel
    ntelem = Inp.ntelem
    nbelem = Inp.nbelem
    ncelem = Inp.ncelem
    nselem = Inp.nselem
    ifw = Inp.ifw
    WindType = Inp.WindType
    AModel = Inp.AModel
    windINPfilename = "$(path)$(Inp.windINPfilename)"
    ifw_libfile = Inp.ifw_libfile
    if ifw_libfile == "nothing"
        ifw_libfile = nothing
    end
    Blade_Height = Inp.Blade_Height # WindIO TODO: resolve DLC dependence
    Blade_Radius = Inp.Blade_Radius # WindIO TODO: resolve DLC dependence
    numTS = Inp.numTS
    delta_t = Inp.delta_t
    NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
    NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
    NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld)"
    NuMad_mat_xlscsv_file_bld = "$(path)$(Inp.NuMad_mat_xlscsv_file_bld)"
    NuMad_geom_xlscsv_file_strut = "$(path)$(Inp.NuMad_geom_xlscsv_file_strut)"
    NuMad_mat_xlscsv_file_strut = "$(path)$(Inp.NuMad_mat_xlscsv_file_strut)"
    adi_lib = Inp.adi_lib
    if adi_lib == "nothing"
        adi_lib = nothing
    end
    adi_rootname = "$(path)$(Inp.adi_rootname)"

    shapeY = blade_z#collect(LinRange(0,H,Nslices+1))
    shapeX = blade_x#R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)
    bshapey = blade_y

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
        eta,
        B=number_of_blades,
        H,
        R,
        shapeY, #TODO: rename to shape Z
        shapeX,
        bshapey, #TODO: rename to shapeY
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
        Ht=towerHeight,
        ntelem, 
        nbelem, 
        ncelem,
        nselem,
        joint_type = 0,
        c_mount_ratio = 0.05,
        AModel, #AD, DMS, AC
        DSModel="BV",
        RPI=true,
        cables_connected_to_blade_base = true,
        meshtype = turbineType)

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

    if AModel=="AD"
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
    aeroLoadsOn = 2)

    nothing

    # Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

    feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    outFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = true,
    gravityOn = [0,0,gravity],
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")

    nothing

    # Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
    # and propogates things in time.

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

    nothing

    # Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
    # deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
    # for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

    println("Saving VTK time domain files")
    OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

    # println("Saving VTK time domain files")
    # OWENS.OWENSFEA_VTK("$path/vtk/XFlowDeflections",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)


    # println("Saving VTK time domain files")
    # userPointNames=["EA","EIyy","EIzz","e_x","e_y","e_z","k_x","k_y","k_z","Fx_Reaction","Fy_Reaction","Fz_Reaction","Mx_Reaction","My_Reaction","Mz_Reaction"]#,"Fx","Fy","Fz","Mx","My","Mz"]
    # # userPointData[iname,it,ipt] = Float64

    # # map el props to points using con
    # userPointData = zeros(length(userPointNames),length(t),mymesh.numNodes)
    # EA_points = zeros(mymesh.numNodes)
    # EIyy_points = zeros(mymesh.numNodes)
    # EIzz_points = zeros(mymesh.numNodes)

    # # Time-invariant data
    # for iel = 1:length(myel.props)
    #     # iel = 1
    #     nodes = mymesh.conn[iel,:]
    #     EA_points[Int.(nodes)] = myel.props[iel].EA
    #     EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
    #     EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
    # end
    

    # epsilon_x_histused = mean(epsilon_x_hist;dims=1)
    # epsilon_y_histused = mean(epsilon_y_hist;dims=1)
    # epsilon_z_histused = mean(epsilon_z_hist;dims=1)
    # kappa_x_histused = mean(kappa_x_hist;dims=1)
    # kappa_y_histused = mean(kappa_y_hist;dims=1)
    # kappa_z_histused = mean(kappa_z_hist;dims=1)

    # # fill in the big matrix
    # for it = 1:length(t)

    #     userPointData[1,it,:] = EA_points
    #     userPointData[2,it,:] = EIyy_points
    #     userPointData[3,it,:] = EIzz_points
    #     for iel = 1:length(myel.props)
    #         nodes = mymesh.conn[iel,:]
    #         userPointData[4,it,Int.(nodes)] .= epsilon_x_histused[1,iel,it] 
    #         userPointData[5,it,Int.(nodes)] .= epsilon_y_histused[1,iel,it] 
    #         userPointData[6,it,Int.(nodes)] .= epsilon_z_histused[1,iel,it] 
    #         userPointData[7,it,Int.(nodes)] .= kappa_x_histused[1,iel,it] 
    #         userPointData[8,it,Int.(nodes)] .= kappa_y_histused[1,iel,it] 
    #         userPointData[9,it,Int.(nodes)] .= kappa_z_histused[1,iel,it] 
    #     end
    #     userPointData[10,it,:] .= FReactionHist[it,1:6:end]
    #     userPointData[11,it,:] .= FReactionHist[it,2:6:end]
    #     userPointData[12,it,:] .= FReactionHist[it,3:6:end]
    #     userPointData[13,it,:] .= FReactionHist[it,4:6:end]
    #     userPointData[14,it,:] .= FReactionHist[it,5:6:end]
    #     userPointData[15,it,:] .= FReactionHist[it,6:6:end]
        
    #     # userPointData[4,it,:] = FReactionHist[it,1:6:end]
    #     # userPointData[5,it,:] = FReactionHist[it,2:6:end]
    #     # userPointData[6,it,:] = FReactionHist[it,3:6:end]
    #     # userPointData[7,it,:] = FReactionHist[it,4:6:end]
    #     # userPointData[8,it,:] = FReactionHist[it,5:6:end]
    #     # userPointData[9,it,:] = FReactionHist[it,6:6:end]
    # end

    # azi=aziHist#./aziHist*1e-6
    # saveName = "$path/vtk/$(windINPfilename[1:end-4])"
    # OWENS.OWENSFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)


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

    dataDumpFilename = "$path/InitialDataOutputs.h5"

    HDF5.h5open(dataDumpFilename, "w") do file
        HDF5.write(file,"t",collect(t))
        HDF5.write(file,"aziHist",aziHist)
        HDF5.write(file,"OmegaHist",OmegaHist)
        HDF5.write(file,"OmegaDotHist",OmegaDotHist)
        HDF5.write(file,"gbHist",gbHist)
        HDF5.write(file,"gbDotHist",gbDotHist)
        HDF5.write(file,"gbDotDotHist",gbDotDotHist)
        HDF5.write(file,"FReactionHist",FReactionHist)
        HDF5.write(file,"FTwrBsHist",FTwrBsHist)
        HDF5.write(file,"genTorque",genTorque)
        HDF5.write(file,"genPower",genPower)
        HDF5.write(file,"torqueDriveShaft",torqueDriveShaft)
        HDF5.write(file,"uHist",uHist)
        HDF5.write(file,"uHist_prp",uHist_prp)
        HDF5.write(file,"epsilon_x_hist",epsilon_x_hist)
        HDF5.write(file,"epsilon_y_hist",epsilon_y_hist)  
        HDF5.write(file,"epsilon_z_hist",epsilon_z_hist)
        HDF5.write(file,"kappa_x_hist",kappa_x_hist)
        HDF5.write(file,"kappa_y_hist",kappa_y_hist)
        HDF5.write(file,"kappa_z_hist",kappa_z_hist) 
        HDF5.write(file,"massOwens",massOwens)
        HDF5.write(file,"stress_U",stress_U)
        HDF5.write(file,"SF_ult_U",SF_ult_U)
        HDF5.write(file,"SF_buck_U",SF_buck_U)
        HDF5.write(file,"stress_L",stress_L)
        HDF5.write(file,"SF_ult_L",SF_ult_L)
        HDF5.write(file,"SF_buck_L",SF_buck_L)
        HDF5.write(file,"stress_TU",stress_TU)
        HDF5.write(file,"SF_ult_TU",SF_ult_TU)
        HDF5.write(file,"SF_buck_TU",SF_buck_TU)
        HDF5.write(file,"stress_TL",stress_TL)
        HDF5.write(file,"SF_ult_TL",SF_ult_TL)
        HDF5.write(file,"SF_buck_TL",SF_buck_TL)
        HDF5.write(file,"topstrainout_blade_U",topstrainout_blade_U)
        HDF5.write(file,"topstrainout_blade_L",topstrainout_blade_L)
        HDF5.write(file,"topstrainout_tower_U",topstrainout_tower_U)
        HDF5.write(file,"topstrainout_tower_L",topstrainout_tower_L)
        HDF5.write(file,"topDamage_blade_U",topDamage_blade_U)
        HDF5.write(file,"topDamage_blade_L",topDamage_blade_L)
        HDF5.write(file,"topDamage_tower_U",topDamage_tower_U)
        HDF5.write(file,"topDamage_tower_L",topDamage_tower_L)
    end
end

# # Test
# Inp = OWENS.MasterInput(;numTS=3,ifw_libfile="$localpath/../../openfastandy/build/modules/inflowwind/libifw_c_binding")
# OWENS.runDLC(["1_1","1_2"],Inp,localpath;Vinf_range=LinRange(5,20,2),regenWindFiles=true,pathtoturbsim="$localpath/../../openfastandy/build/modules/turbsim/turbsim")
"""

runDLC(DLCs,Inp,path;
    Vinf_range=LinRange(5,20,16),
    IEC_std="\"2\"",
    WindChar="\"A\"",
    WindClass=1,
    turbsimpath="./turbsimfiles",
    templatefile="./templateTurbSim.inp",
    pathtoturbsim="../../openfast/build/modules/turbsim/turbsim",
    NumGrid_Z=100,
    NumGrid_Y=100,
    Vref=10.0,
    Vdesign=11.0,
    grid_oversize=1.1,
    regenWindFiles=false)

    # Input
    * `DLCs`: ["1_1","1_2"]
    * `Inp::MasterInput`: see ?OWENS.MasterInput
    * `path`: desired path to run everything
    * `Vinf_range`: =LinRange(5,20,16),
    * `IEC_std`: ="\"2\"",
    * `WindChar`: ="\"A\"",
    * `WindClass`: =1,
    * `turbsimpath`: ="./turbsimfiles", path where it dumps the turbsim files
    * `templatefile`: ="./template_files/templateTurbSim.inp",
    * `pathtoturbsim`: ="../../openfast/build/modules/turbsim/turbsim",
    * `NumGrid_Z`: =100,
    * `NumGrid_Y`: =100,
    * `Vref`: =10.0,
    * `Vdesign`: =11.0, # Design or rated speed
    * `grid_oversize`: =1.1,
    * `regenWindFiles`: =false

    # Output
    * `nothing`: 
    """
function runDLC(DLCs,Inp,path;
    Vinf_range=LinRange(5,20,16),
    IEC_std="\"1-ED3\"",
    WindChar="\"A\"",
    WindClass=1,
    turbsimpath="./turbsimfiles",
    templatefile="$module_path/template_files/templateTurbSim.inp",
    pathtoturbsim="../../openfast/build/modules/turbsim/turbsim",
    NumGrid_Z=nothing,
    NumGrid_Y=nothing,
    Vref=10.0,
    Vdesign=11.0, # Design or rated speed
    grid_oversize=1.1,
    regenWindFiles=false,
    delta_t_turbsim=nothing,
    simtime_turbsim=nothing,
    runScript = OWENS.runOWENS)

    if !isdir(turbsimpath)
        mkdir(turbsimpath)
    end

    # Fill in DLC parameters based on model inputs
    DLCParams = Array{DLCParameters, 1}(undef, length(DLCs))

    for (iDLC, DLC) in enumerate(DLCs) #TODO parallelize this

        DLCParams[iDLC] = getDLCparams(DLC, Inp, Vinf_range, Vdesign, Vref, WindChar,WindClass, IEC_std;grid_oversize,simtime_turbsim,delta_t_turbsim,NumGrid_Z,NumGrid_Y)


        # Run Simulation at each Wind Speed
        for windspeed in DLCParams[iDLC].Vinf_range_used #TODO: parallelize this

            

            DLCParams[iDLC].URef = windspeed
            # Check if turbulent inflow file exists, if not create it
            windspeedStr = round(windspeed;digits=2)
            windspeedStr = lpad(windspeedStr,4,"0")
            println("Running DLC $DLC at Vinf $windspeedStr m/s")
            windINPfilename = "$turbsimpath/DLC$(DLC)Vinf$(windspeedStr).inp"
            
            if contains(DLCParams[iDLC].IEC_WindType, "NTM") || contains(DLCParams[iDLC].IEC_WindType, "ETM") || contains(DLCParams[iDLC].IEC_WindType, "EWM")
                if !isfile(windINPfilename) || regenWindFiles
                    generateTurbsimBTS(DLCParams[iDLC],windINPfilename,pathtoturbsim;templatefile)
                end
                Inp.WindType = 3
                Inp.windINPfilename = "$(windINPfilename[1:end-4]).bts"
            else
                if !isfile(windINPfilename) || regenWindFiles
                    generateUniformwind(DLCParams[iDLC],windINPfilename)
                end
                Inp.windINPfilename = windINPfilename
                Inp.WindType = 2
            end

            Inp.ifw = true
            Inp.controlStrategy = DLCParams[iDLC].controlStrategy
            # run owens simulation
            runScript(Inp,path)
        end
    end
end

mutable struct DLCParameters
    Vinf_range_used
    analysis_type # "U", "F", "UF"
    controlStrategy # "constRPM", function handle
    RandSeed1 # Turbulent Random Seed Number
    NumGrid_Z # Vertical grid-point matrix dimension
    NumGrid_Y # Horizontal grid-point matrix dimension
    TimeStepSim # Time step [s]
    TimeStep # Time step [s]
    HubHt # Hub height [m] (should be > 0.5*GridHeight)
    AnalysisTime # Length of analysis time series [s] (program will add time if necessary)
    GridHeight # Grid height [m]
    GridWidth # Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))
    VFlowAng # Vertical mean flow (uptilt) angle [degrees]
    HFlowAng # Horizontal mean flow (skew) angle [degrees]
    TurbModel # Turbulence model (see Table 4 for valid codes)
    IECstandard # Number of the IEC standard (61400-x, x=1,2,3) with optional 61400-1 ed. number
    IECturbc # IEC turbulence characteristic ("A", "B", "C" or TI in %) or KHTEST
    IEC_WindType # IEC turbulence type ("NTM", "xETM", "xEWM1", or "xEWM50" for x=class 1, 2, or 3)
    RefHt # Height of the reference wind speed [m]
    URef # Mean wind speed at the reference height [m/s] [must be 1-hr mean for API model]
    time
    windvel
    winddir
    windvertvel
    horizshear
    pwrLawVertShear
    LinVertShear
    gustvel
    UpflowAngle
end


function getDLCparams(DLC, Inp, Vinf_range, Vdesign, Vref, WindChar, WindClass, IEC_std;grid_oversize=1.2,simtime_turbsim=nothing,delta_t_turbsim=nothing,NumGrid_Z=nothing,NumGrid_Y=nothing)

    Ve50 = 50.0 #TODO change by class etc
    Ve1 = 30.0 #TODO

    numTS = Inp.numTS
    delta_t = Inp.delta_t
    simtime = numTS*delta_t

    GridHeight = (Inp.towerHeight-Inp.Blade_Height/2+Inp.Blade_Height)*grid_oversize
    GridWidth = Inp.Blade_Radius * 2.0 * grid_oversize
    HubHt = GridHeight*2/3

    if !isnothing(NumGrid_Z)
        NumGrid_Z = NumGrid_Z #Inp.ntelem+Inp.nbelem
        NumGrid_Y = NumGrid_Y #Inp.ntelem+Inp.nbelem
    else
        NumGrid_Z = Inp.ntelem+Inp.nbelem
        NumGrid_Y = Inp.nbelem
    end

    RandSeed1 = 40071 #TODO
    
    if !isnothing(simtime_turbsim)
        AnalysisTime = simtime_turbsim
    else
        AnalysisTime = simtime
    end

    VFlowAng = 0.0
    HFlowAng = 0.0
    
    IECstandard = IEC_std
    IECturbc = WindChar
    TurbModel = "\"IECKAI\""
    
    RefHt = round(Inp.towerHeight) #TODO: what if tower doesn't extend into blade z level
    URef = 0.0 #gets filled in later from the Vinf_range when the .bst is generated

    TimeStepSim = delta_t
    if !isnothing(delta_t_turbsim)
        TimeStep = delta_t_turbsim
    else
        TimeStep = delta_t
    end

    time = LinRange(0,10,10)
    windvel = nothing # gets supersceded   
    winddir = nothing  
    windvertvel = nothing  
    horizshear = nothing  
    pwrLawVertShear = nothing  
    LinVertShear = nothing  
    gustvel = nothing  
    UpflowAngle = nothing  

    if contains(IEC_std,"1-")
        if DLC == "1_1" || DLC == "1_2"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""                        

        elseif DLC == "1_3"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ETM\""
            
        elseif DLC == "1_4"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign-2.0,Vdesign+2.0]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = [0,0,45,90,45,0,0.0] 
            windvertvel = zeros(length(time))   
            horizshear = zeros(length(time))  
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = zeros(length(time))  
            gustvel = [0,0,7.0,15,7.5,0.0,0.0]  
            UpflowAngle = zeros(length(time))    

        elseif DLC == "1_5"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWS\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = zeros(length(time))  
            windvertvel = zeros(length(time))   
            horizshear = [0,0,5.0,0,0,0,0]#ones(length(time)).*10.0   
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = [0,0,0,0,5.0,0,0]#ones(length(time)).*10.0 
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))  
            

        elseif DLC == "2_1" || DLC == "2_2" || DLC == "2_4"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "2_3"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = zeros(length(time))      
            LinVertShear = zeros(length(time))   
            UpflowAngle = zeros(length(time))     

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
        elseif DLC == "3_1"
            ControlStrategy = "startup"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0,30,10)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = ones(length(time)).*0.2  
            LinVertShear = zeros(length(time))   
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))     
            
        elseif DLC == "3_2"
            ControlStrategy = "startup"
            Vinf_range_used = [Vinf_range[1];collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = zeros(length(time))      
            LinVertShear = zeros(length(time))   
            UpflowAngle = zeros(length(time))     

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
        elseif DLC == "3_3"
            ControlStrategy = "startup"
            Vinf_range_used = [Vinf_range[1];collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD\""

            time = [0,10,15,20,25,30,10000.0]#LinRange(0,10,10)
            winddir = [0,0,45,90,45,0,0.0] 
            windvertvel = zeros(length(time))   
            horizshear = zeros(length(time))  
            pwrLawVertShear = zeros(length(time))   
            LinVertShear = zeros(length(time))  
            gustvel = [0,0,7.0,15,7.5,0.0,0.0]  
            UpflowAngle = zeros(length(time)) 
            
        elseif DLC == "4_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NWP\""

            time = LinRange(0,30,10)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = ones(length(time)).*0.2  
            LinVertShear = zeros(length(time))   
            gustvel = zeros(length(time))   
            UpflowAngle = zeros(length(time))     
            
        elseif DLC == "4_2"
            ControlStrategy = "shutdown"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG\""

            time = LinRange(0,30,100)
            winddir = zeros(length(time))     
            windvertvel = zeros(length(time))      
            horizshear = zeros(length(time))   
            pwrLawVertShear = zeros(length(time))      
            LinVertShear = zeros(length(time))   
            UpflowAngle = zeros(length(time))     

            time_delay = 10.0 #sec
            time_delay2 = 20.0
            G_amp = 15.0 #m/s
            gustT = 10.0
            gustvel = simpleGustVel.(time, time_delay, G_amp,gustT) .+ simpleGustVel.(time, time_delay2, G_amp,gustT)
            
            
        elseif DLC == "5_1"
            ControlStrategy = "emergencyshutdown"
            Vinf_range_used = [collect(LinRange(Vdesign-2.0,Vdesign+2.0,2));Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""
            
        elseif DLC == "6_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM50\""
            
        elseif DLC == "6_2"
            ControlStrategy = "parked_idle"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM50\""

        elseif DLC == "6_3"
            ControlStrategy = "parked_yaw"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""

        elseif DLC == "6_4"
            ControlStrategy = "parked"
            Vinf_range_used = [0.7*Ve50]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""

        elseif DLC == "7_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""
            
        elseif DLC == "8_1" #Startup
            ControlStrategy = "transport"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM1\""
            
        else
            error("IEC61400_1 DLCs such as 1_1, 1_2 defined, you requested $DLC")
        end

    elseif contains(IEC_std,"2")
        error("IEC61400_2 DLCs are not fully defined")
        if DLC == "1_1"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "1_2"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)ECD\""
            

        elseif DLC == "1_3"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG50\""
            

        elseif DLC == "1_4"
            ControlStrategy = "normal"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECD50\""
            

        elseif DLC == "1_5"
            ControlStrategy = "normal"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)ECG\""
            

        elseif DLC == "2_1"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NWP\""
            

        elseif DLC == "2_2"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "UF"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "2_3"
            ControlStrategy = "freewheelatNormalOperatingRPM"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG1\""
            

        elseif DLC == "3_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = Vinf_range
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "3_2"
            ControlStrategy = "shutdown"
            Vinf_range_used = [Vinf_range[end]]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EOG1\""
            

        elseif DLC == "4_1"
            ControlStrategy = "shutdown"
            Vinf_range_used = [Vdesign]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "5_1"
            ControlStrategy = "freewheelatIdle"
            Vinf_range_used = [Ve50]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""
            

        elseif DLC == "5_2"
            ControlStrategy = "idle"
            Vinf_range_used = [Vdesign]
            analysis_type = "F"
            IEC_WindType = "\"$(WindClass)NTM\""
            

        elseif DLC == "6_1"
            ControlStrategy = "parked"
            Vinf_range_used = [Ve1]
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""
            
        elseif DLC == "8_1" #Startup
            ControlStrategy = "startup"
            Vinf_range_used = Vinf_range
            analysis_type = "U"
            IEC_WindType = "\"$(WindClass)EWM\""
            
        else
            error("IEC61400_2 DLCs [1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,3.1,3.2,4.1,5.1,5.2,6.1] defined, you requested $DLC")
        end
    else
        error("IEC_std 61400 1-ED3 and 2 defined, you requested $IEC_std")
    end

    return DLCParameters(
        Vinf_range_used,
        analysis_type, # array of windspeeds m/s
        ControlStrategy, # "constRPM", function handle
        RandSeed1, # Turbulent Random Seed Number
        NumGrid_Z, # Vertical grid-point matrix dimension
        NumGrid_Y, # Horizontal grid-point matrix dimension
        TimeStepSim, # Time step [s]
        TimeStep, # Turbsim time step [s]
        HubHt, # Hub height [m] (should be > 0.5*GridHeight)
        AnalysisTime, # Length of analysis time series [s] (program will add time if necessary)
        GridHeight, # Grid height [m]
        GridWidth, # Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))
        VFlowAng, # Vertical mean flow (uptilt) angle [degrees]
        HFlowAng, # Horizontal mean flow (skew) angle [degrees]
        TurbModel, # Turbulence model (see Table 4 for valid codes)
        IECstandard, # Number of the IEC standard (61400-x, x=1,2,3) with optional 61400-1 ed. number
        IECturbc, # IEC turbulence characteristic ("A", "B", "C" or TI in %) or KHTEST
        IEC_WindType, # IEC turbulence type ("NTM", "xETM", "xEWM1", or "xEWM50" for x=class 1, 2, or 3)
        RefHt, # Height of the reference wind speed [m]
        URef, # Mean wind speed at the reference height [m/s] [must be 1-hr mean for API model]
        time,
        windvel,
        winddir,
        windvertvel,
        horizshear,
        pwrLawVertShear,
        LinVertShear,
        gustvel,
        UpflowAngle,
    )
end

function generateUniformwind(DLCParams,windINPfilename)

    time = DLCParams.time
    windvel = ones(length(DLCParams.time)) .* DLCParams.URef
    winddir = DLCParams.winddir
    windvertvel = DLCParams.windvertvel
    horizShear = DLCParams.horizshear
    pwrLawVertShear = DLCParams.pwrLawVertShear
    LinVertShear = DLCParams.LinVertShear
    gustvel = DLCParams.gustvel
    UpflowAngle = DLCParams.UpflowAngle

    lines = ["! OpenFAST Deterministic Wind File",
    "#",
    "# Comment lines begin with \"!\" or \"#\" or \"%\", then the data lines must contain the following columns:",
    "#",
    "# If there are only 8 columns, upflow is assumed to be 0.",
    "#",
    "# Parameters are interpolated linearly between time steps; using nearest neighbor before the first time ",
    "# listed in this file and after the last time listed in the file. ",
    "#",
    "! Time     Wind    Wind    Vertical    Horiz.      Pwr.Law     Lin.Vert.   Gust     Upflow",
    "!          Speed   Dir     Speed       Shear       Vert.Shr    Shear       Speed    Angle ",
    "! (sec)    (m/s)   (Deg)   (m/s)                                            (m/s)   (deg)"]
    # "0.000000   10   0.000000   0          0.000000   0.300000   0.000000   0.000000      8",
    # "10.000000   12   0.000000   0          0.000000   0.300000   0.000000   0.000000      8",

    for itime = 1:length(time)
        lines = [lines; "$(time[itime]) $(windvel[itime]) $(winddir[itime]) $(windvertvel[itime]) $(horizShear[itime]) $(pwrLawVertShear[itime]) $(LinVertShear[itime]) $(gustvel[itime]) $(UpflowAngle[itime])"]
    end

    
    # Write the new file
    open(windINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end
end

function generateTurbsimBTS(DLCParams,windINPfilename,pathtoturbsim;templatefile="$localpath/templateTurbSim.inp") 

    lines = readlines(templatefile)

    for fieldname in fieldnames(typeof(DLCParams))
        turbsimKeyName = String(fieldname)
        myvalue = getfield(DLCParams,fieldname)
        if turbsimKeyName != "Vinf_range_used" || turbsimKeyName != "analysis_type" || turbsimKeyName != "ControlStrategy"# || other strings
            for (iline,line) in enumerate(lines)
                if contains(line," - ") #TODO: this assumes that the keys aren't in the comments
                    linenocomments,comments = split(line," - ")
                    if contains(linenocomments,turbsimKeyName)
                        value,descriptor = split(linenocomments)
                        newline = "$myvalue $turbsimKeyName - $comments"
                        lines[iline] = newline
                        break
                    end
                end
            end
        end
    end
    
    # Write the new file
    open(windINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end

    run(`$pathtoturbsim $windINPfilename`)
end

"""
* `time::TF`: in seconds
* `nominalVinf::TF`: Nominal velocity used to calculate the IEC gust size (m/s)
* `R::TF`: Turbine Radius (m)
* `G_amp::TF`: IEC gust amplitude (m/s)
* `gustT::TF`: IEC gust duration (s)
* `gustDelayT::TF`: IEC gust delay time
"""
function getGustVel(time,nominalVinf,R,G_amp,gustT,gustDelayT)
    ele_x = 0.0 #TODO: I don't think inflowwind takes in account the 3D nature of a vawt

    gustT = gustT * nominalVinf / R
    tr = time .- ele_x .- gustDelayT / R
    if (tr >= 0) && (tr<=gustT)
        IECGustFactor = 1.0 - 0.37 * G_amp/nominalVinf * sin(3*pi*tr/gustT)  * (1.0 - cos(2*pi*tr/gustT))
        return nominalVinf*IECGustFactor
    else
        return nominalVinf
    end

end

function simpleGustVel(time, time_delay, G_amp,gustT)
    timeused = time - time_delay
    if (timeused >= 0) && (timeused<=gustT)
        gustV = -0.37 * G_amp * sin(3*pi*timeused/gustT)  * (1.0 - cos(2*pi*timeused/gustT))
    else
        gustV = 0.0
    end
    return gustV
end

# # Test
# Inp = OWENS.MasterInput(;numTS=3,ifw_libfile="$localpath/../../openfastandy/build/modules/inflowwind/libifw_c_binding")
# OWENS.runDLC(["1_1","1_2"],Inp,localpath;Vinf_range=LinRange(5,20,2),regenWindFiles=true,pathtoturbsim="$localpath/../../openfastandy/build/modules/turbsim/turbsim")
