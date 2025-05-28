using OrderedCollections
import ..OWENS: MeshConfig, TowerConfig, BladeConfig, MaterialConfig, AeroConfig

function setupOWENS(
    OWENSAero,
    path;
    rho = 1.225,
    mu = 1.7894e-5,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    shapeZ = collect(LinRange(0, H, Nslices+1)),
    shapeX = R .* (1.0 .- 4.0 .* (shapeZ/H .- 0.5) .^ 2),#shapeX_spline(shapeZ)
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    nselem = 5,
    shapeY = zeros(Nslices+1),
    ifw = false,
    AD15hubR = 0.1,
    WindType = 1,
    delta_t = 0.01,
    numTS = 100,
    adi_lib = "$(path)../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    windINPfilename = "$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_mat_xlscsv_file_twr::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_geom_xlscsv_file_bld::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_mat_xlscsv_file_bld::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_geom_xlscsv_file_strut::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_mat_xlscsv_file_strut::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_geom_xlscsv_file_intra_blade_cable::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_mat_xlscsv_file_intra_blade_cable::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_geom_xlscsv_file_guys::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    NuMad_mat_xlscsv_file_guys::Union{Nothing,String,OrderedDict{Symbol,Any}} = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0, 1.0],
    chord_scale = [1.0, 1.0],
    thickness_scale = [1.0, 1.0],
    Htwr_base = 2.0,
    Htwr_blds = H,
    strut_twr_mountpoint = [0.25, 0.75],
    strut_bld_mountpoint = [0.25, 0.75],
    joint_type = 2,
    c_mount_ratio = 0.05,
    angularOffset = -pi/2,
    AeroModel = "DMS",
    DynamicStallModel = "BV",
    RPI = true,
    Aero_AddedMass_Active = false,
    Aero_RotAccel_Active = false,
    Aero_Buoyancy_Active = false,
    cables_connected_to_blade_base = true,
    meshtype = "Darrieus",
    custommesh = nothing,
    AddedMass_Coeff_Ca = 0.0,
    verbosity = 1,
    VTKmeshfilename = nothing,
    return_componentized = false,
)

    # Bundle the inputs into config structs
    mesh_config = MeshConfig(
        Nslices = Nslices,
        ntheta = ntheta,
        ntelem = ntelem,
        nbelem = nbelem,
        ncelem = ncelem,
        nselem = nselem,
        meshtype = meshtype,
        custommesh = custommesh,
        connectBldTips2Twr = meshtype == "Darrieus", # Connect blade tips to tower for Darrieus type
        AD15_ccw = true,
    )

    tower_config = TowerConfig(
        Htwr_base = Htwr_base,
        Htwr_blds = Htwr_blds,
        strut_twr_mountpoint = strut_twr_mountpoint,
        strut_bld_mountpoint = strut_bld_mountpoint,
        joint_type = joint_type,
        c_mount_ratio = c_mount_ratio,
        angularOffset = angularOffset,
        NuMad_geom_xlscsv_file_twr = NuMad_geom_xlscsv_file_twr,
        NuMad_mat_xlscsv_file_twr = NuMad_mat_xlscsv_file_twr,
    )

    blade_config = BladeConfig(
        B = B,
        H = H,
        R = R,
        shapeZ = shapeZ,
        shapeX = shapeX,
        shapeY = shapeY,
        NuMad_geom_xlscsv_file_bld = NuMad_geom_xlscsv_file_bld,
        NuMad_mat_xlscsv_file_bld = NuMad_mat_xlscsv_file_bld,
    )

    material_config = MaterialConfig(
        stack_layers_bld = stack_layers_bld,
        stack_layers_scale = stack_layers_scale,
        chord_scale = chord_scale,
        thickness_scale = thickness_scale,
        AddedMass_Coeff_Ca = AddedMass_Coeff_Ca,
    )

    aero_config = AeroConfig(
        rho = rho,
        mu = mu,
        RPM = RPM,
        Vinf = Vinf,
        eta = eta,
        delta_t = delta_t,
        AD15hubR = AD15hubR,
        WindType = WindType,
        AeroModel = AeroModel,
        DynamicStallModel = DynamicStallModel,
        numTS = numTS,
        adi_lib = adi_lib,
        adi_rootname = adi_rootname,
        windINPfilename = windINPfilename,
        ifw_libfile = ifw_libfile,
        ifw = ifw,
        RPI = RPI,
        Aero_AddedMass_Active = Aero_AddedMass_Active,
        Aero_RotAccel_Active = Aero_RotAccel_Active,
        Aero_Buoyancy_Active = Aero_Buoyancy_Active,
        centrifugal_force_flag = AddedMass_Coeff_Ca > 0.0,
    )

    # Run the setup function with the config structs
    return setupOWENS_config(
        OWENSAero,
        path;
        verbosity = verbosity,
        VTKmeshfilename = VTKmeshfilename,
        return_componentized = return_componentized,
        mesh_config = mesh_config,
        tower_config = tower_config,
        blade_config = blade_config,
        material_config = material_config,
        aero_config = aero_config,
    )
end