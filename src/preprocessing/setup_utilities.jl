using OrderedCollections: OrderedDict

export MeshSetupOptions,
    TowerSetupOptions,
    BladeSetupOptions,
    MaterialSetupOptions,
    AeroSetupOptions,
    SetupOptions

"""
    MeshSetupOptions

Contains the configuration for the mesh.

# Fields
- `Nslices::Int`: The number of slices.
- `ntheta::Int`: The number of theta points.
- `ntelem::Int`: The number of elements.
- `nbelem::Int`: The number of blade elements.
- `ncelem::Int`: The number of chord elements.
- `nselem::Int`: The number of span elements.
- `meshtype::String`: The type of mesh.
- `custommesh::Union{Nothing,Function}`: The custom mesh function.
- `connectBldTips2Twr::Bool`: Whether to connect the blade tips to the tower.
- `AD15_ccw::Bool`: Whether to use the AD15 convention of VAWT counter-clockwise with blade root at top (blade points down).
"""
mutable struct MeshSetupOptions
    Nslices::Int
    ntheta::Int
    ntelem::Int
    nbelem::Int
    ncelem::Int
    nselem::Int
    meshtype::String
    custommesh::Union{Nothing,Function} # this doesn't appear to be used anywhere
    connectBldTips2Twr::Bool
    AD15_ccw::Bool

    function MeshSetupOptions(
        Nslices,
        ntheta,
        ntelem,
        nbelem,
        ncelem,
        nselem,
        meshtype,
        custommesh,
        connectBldTips2Twr,
        AD15_ccw,
    )
        new(
            Nslices,
            ntheta,
            ntelem,
            nbelem,
            ncelem,
            nselem,
            meshtype,
            custommesh,
            connectBldTips2Twr,
            AD15_ccw,
        )
    end
    function MeshSetupOptions(;
        Nslices::Int,
        ntheta::Int,
        ntelem::Int,
        nbelem::Int,
        ncelem::Int,
        nselem::Int,
        meshtype::String,
        custommesh::Union{Nothing,Function} = nothing,
        connectBldTips2Twr::Bool = false,
        AD15_ccw::Bool = true,
    )
        new(
            Nslices,
            ntheta,
            ntelem,
            nbelem,
            ncelem,
            nselem,
            meshtype,
            custommesh,
            connectBldTips2Twr,
            AD15_ccw,
        )
    end
    # Default constructor
    function MeshSetupOptions()
        new(
            30,  # Nslices
            30,  # ntheta
            10,  # ntelem
            60,  # nbelem
            10,  # ncelem
            5,   # nselem
            "Darrieus",  # meshtype
            nothing,  # custommesh
            false,  # connectBldTips2Twr
            true,   # AD15_ccw
        )
    end
    function MeshSetupOptions(dict::OrderedDict{Symbol,Any})
        new(
            get(dict, :Nslices, 30),
            get(dict, :ntheta, 30),
            get(dict, :ntelem, 10),
            get(dict, :nbelem, 60),
            get(dict, :ncelem, 10),
            get(dict, :nselem, 5),
            get(dict, :meshtype, "Darrieus"),
            get(dict, :custommesh, nothing),
            get(dict, :connectBldTips2Twr, false),
            get(dict, :AD15_ccw, true),
        )
    end
end

"""
    TowerSetupOptions

Contains the configuration for the tower.

# Fields
- `Htwr_base::Float64`: The height of the tower base.
- `Htwr_blds::Float64`: The height of the tower blades.
- `strut_twr_mountpoint::Vector{Float64}`: The mount points of the struts on the tower.
- `strut_bld_mountpoint::Vector{Float64}`: The mount points of the struts on the blades.
- `joint_type::Int`: The type of joint between the struts and tower.
- `c_mount_ratio::Float64`: The mount point of the struts to the tower.
- `angularOffset::Float64`: The angular offset of the tower.
- `NuMad_geom_xlscsv_file_twr::Any`: The path to the tower geometry file. Can be nothing, String, or Vector{String}.
- `NuMad_mat_xlscsv_file_twr::Any`: The path to the tower material file. Can be nothing, String, or Vector{String}.
- `NuMad_geom_xlscsv_file_strut::Any`: The path to the strut geometry file. Can be nothing, String, or Vector{String}.
- `NuMad_mat_xlscsv_file_strut::Any`: The path to the strut material file. Can be nothing, String, or Vector{String}.
- `strut_tower_joint_type::Int`: The type of joint between the struts and tower.
"""
mutable struct TowerSetupOptions
    Htwr_base::Float64
    Htwr_blds::Float64
    strut_twr_mountpoint::Vector{Float64}
    strut_bld_mountpoint::Vector{Float64}
    joint_type::Int
    c_mount_ratio::Float64
    angularOffset::Float64
    NuMad_geom_xlscsv_file_twr::Any
    NuMad_mat_xlscsv_file_twr::Any
    NuMad_geom_xlscsv_file_strut::Any
    NuMad_mat_xlscsv_file_strut::Any
    strut_tower_joint_type::Int

    function TowerSetupOptions(
        Htwr_base,
        Htwr_blds,
        strut_twr_mountpoint,
        strut_bld_mountpoint,
        joint_type,
        c_mount_ratio,
        angularOffset,
        NuMad_geom_xlscsv_file_twr,
        NuMad_mat_xlscsv_file_twr,
        NuMad_geom_xlscsv_file_strut,
        NuMad_mat_xlscsv_file_strut,
        strut_tower_joint_type,
    )
        new(
            Htwr_base,
            Htwr_blds,
            strut_twr_mountpoint,
            strut_bld_mountpoint,
            joint_type,
            c_mount_ratio,
            angularOffset,
            NuMad_geom_xlscsv_file_twr,
            NuMad_mat_xlscsv_file_twr,
            NuMad_geom_xlscsv_file_strut,
            NuMad_mat_xlscsv_file_strut,
            strut_tower_joint_type,
        )
    end
    function TowerSetupOptions(;
        Htwr_base::Float64,
        Htwr_blds::Float64,
        strut_twr_mountpoint::Vector{Float64},
        strut_bld_mountpoint::Vector{Float64},
        joint_type::Int,
        c_mount_ratio::Float64,
        angularOffset::Float64,
        NuMad_geom_xlscsv_file_twr::Any = nothing,
        NuMad_mat_xlscsv_file_twr::Any = nothing,
        NuMad_geom_xlscsv_file_strut::Any = nothing,
        NuMad_mat_xlscsv_file_strut::Any = nothing,
        strut_tower_joint_type::Int = 2,
    )
        new(
            Htwr_base,
            Htwr_blds,
            strut_twr_mountpoint,
            strut_bld_mountpoint,
            joint_type,
            c_mount_ratio,
            angularOffset,
            NuMad_geom_xlscsv_file_twr,
            NuMad_mat_xlscsv_file_twr,
            NuMad_geom_xlscsv_file_strut,
            NuMad_mat_xlscsv_file_strut,
            strut_tower_joint_type,
        )
    end
    # Default constructor
    function TowerSetupOptions()
        new(
            2.0,           # Htwr_base
            5.0,           # Htwr_blds
            [0.25, 0.75],  # strut_twr_mountpoint
            [0.25, 0.75],  # strut_bld_mountpoint
            2,             # joint_type
            0.05,          # c_mount_ratio
            -pi/2,         # angularOffset
            nothing,       # NuMad_geom_xlscsv_file_twr
            nothing,       # NuMad_mat_xlscsv_file_twr
            nothing,       # NuMad_geom_xlscsv_file_strut
            nothing,       # NuMad_mat_xlscsv_file_strut
            2,              # strut_tower_joint_type
        )
    end
    function TowerSetupOptions(dict::OrderedDict{Symbol,Any})
        new(
            get(dict, :Htwr_base, 2.0),
            get(dict, :Htwr_blds, 5.0),
            get(dict, :strut_twr_mountpoint, [0.25, 0.75]),
            get(dict, :strut_bld_mountpoint, [0.25, 0.75]),
            get(dict, :joint_type, 2),
            get(dict, :c_mount_ratio, 0.05),
            get(dict, :angularOffset, -pi/2),
            get(dict, :NuMad_geom_xlscsv_file_twr, nothing),
            get(dict, :NuMad_mat_xlscsv_file_twr, nothing),
            get(dict, :NuMad_geom_xlscsv_file_strut, nothing),
            get(dict, :NuMad_mat_xlscsv_file_strut, nothing),
            get(dict, :strut_tower_joint_type, 2),
        )
    end
end

"""
    BladeSetupOptions

Contains the configuration for the blades.

# Fields
- `B::Int`: The number of blades.
- `H::Float64`: The height of the blades.
- `R::Float64`: The radius of the blades.
- `shapeZ::Vector{Float64}`: The z-coordinates of the blade.
- `shapeX::Vector{Float64}`: The x-coordinates of the blade.
- `shapeY::Vector{Float64}`: The y-coordinates of the blade.
- `NuMad_geom_xlscsv_file_bld::Any`: The path to the blade geometry file. Can be nothing, String, or Vector{String}.
- `NuMad_mat_xlscsv_file_bld::Any`: The path to the blade material file. Can be nothing, String, or Vector{String}.
- `strut_blade_joint_type::Int`: The type of joint between the struts and blades.
- `blade_joint_angle_Degrees::Float64`: The angle of the blade joint in degrees.
"""
mutable struct BladeSetupOptions
    B::Int
    H::Float64
    R::Float64
    shapeZ::Vector{Float64}
    shapeX::Vector{Float64}
    shapeY::Vector{Float64}
    NuMad_geom_xlscsv_file_bld::Any
    NuMad_mat_xlscsv_file_bld::Any
    strut_blade_joint_type::Int
    blade_joint_angle_Degrees::Float64

    function BladeSetupOptions(
        B,
        H,
        R,
        shapeZ,
        shapeX,
        shapeY,
        NuMad_geom_xlscsv_file_bld,
        NuMad_mat_xlscsv_file_bld,
        strut_blade_joint_type,
        blade_joint_angle_Degrees,
    )
        new(
            B,
            H,
            R,
            shapeZ,
            shapeX,
            shapeY,
            NuMad_geom_xlscsv_file_bld,
            NuMad_mat_xlscsv_file_bld,
            strut_blade_joint_type,
            blade_joint_angle_Degrees,
        )
    end
    function BladeSetupOptions(;
        B::Int,
        H::Float64,
        R::Float64,
        shapeZ::Vector{Float64},
        shapeX::Vector{Float64},
        shapeY::Vector{Float64},
        NuMad_geom_xlscsv_file_bld::Any = nothing,
        NuMad_mat_xlscsv_file_bld::Any = nothing,
        strut_blade_joint_type::Int = 0,
        blade_joint_angle_Degrees::Float64 = 0.0,
    )
        new(
            B,
            H,
            R,
            shapeZ,
            shapeX,
            shapeY,
            NuMad_geom_xlscsv_file_bld,
            NuMad_mat_xlscsv_file_bld,
            strut_blade_joint_type,
            blade_joint_angle_Degrees,
        )
    end
    # Default constructor
    function BladeSetupOptions()
        new(
            3,  # B
            5.0,  # H
            2.5,  # R
            collect(LinRange(0, 5.0, 31)),  # shapeZ
            2.5 .* (1.0 .- 4.0 .* (collect(LinRange(0, 5.0, 31))/5.0 .- 0.5) .^ 2),  # shapeX
            zeros(31),  # shapeY
            nothing,  # NuMad_geom_xlscsv_file_bld
            nothing,   # NuMad_mat_xlscsv_file_bld
            0,  # strut_blade_joint_type
            0.0,  # blade_joint_angle_Degrees
        )
    end
    function BladeSetupOptions(dict::OrderedDict{Symbol,Any})
        new(
            get(dict, :B, 3),
            get(dict, :H, 5.0),
            get(dict, :R, 2.5),
            get(dict, :shapeZ, collect(LinRange(0, 5.0, 31))),
            get(
                dict,
                :shapeX,
                2.5 .* (1.0 .- 4.0 .* (collect(LinRange(0, 5.0, 31))/5.0 .- 0.5) .^ 2),
            ),
            get(dict, :shapeY, zeros(31)),
            get(dict, :NuMad_geom_xlscsv_file_bld, nothing),
            get(dict, :NuMad_mat_xlscsv_file_bld, nothing),
            get(dict, :strut_blade_joint_type, 0),
            get(dict, :blade_joint_angle_Degrees, 0.0),
        )
    end
end

"""
    MaterialSetupOptions

Contains the material properties for the blades, struts, and tower.

# Fields
- `stack_layers_bld::Union{Nothing,Matrix{Float64}}`: The stack layers for the blades.
- `stack_layers_scale::Vector{Float64}`: The scale factors for the stack layers.
- `chord_scale::Vector{Float64}`: The scale factors for the chord.
- `thickness_scale::Vector{Float64}`: The scale factors for the thickness.
- `AddedMass_Coeff_Ca::Float64`: The added mass coefficient.
"""
mutable struct MaterialSetupOptions
    stack_layers_bld::Union{Nothing,Matrix{Float64}}
    stack_layers_scale::Vector{Float64}
    chord_scale::Vector{Float64}
    thickness_scale::Vector{Float64}
    AddedMass_Coeff_Ca::Float64

    function MaterialSetupOptions(
        stack_layers_bld,
        stack_layers_scale,
        chord_scale,
        thickness_scale,
        AddedMass_Coeff_Ca,
    )
        new(
            stack_layers_bld,
            stack_layers_scale,
            chord_scale,
            thickness_scale,
            AddedMass_Coeff_Ca,
        )
    end
    function MaterialSetupOptions(;
        stack_layers_bld::Union{Nothing,Matrix{Float64}} = nothing,
        stack_layers_scale::Vector{Float64} = [1.0, 1.0],
        chord_scale::Vector{Float64} = [1.0, 1.0],
        thickness_scale::Vector{Float64} = [1.0, 1.0],
        AddedMass_Coeff_Ca::Float64 = 0.0,
    )
        new(
            stack_layers_bld,
            stack_layers_scale,
            chord_scale,
            thickness_scale,
            AddedMass_Coeff_Ca,
        )
    end
    # Default constructor
    function MaterialSetupOptions()
        new(
            nothing,      # stack_layers_bld
            [1.0, 1.0],   # stack_layers_scale
            [1.0, 1.0],   # chord_scale
            [1.0, 1.0],   # thickness_scale
            0.0,           # AddedMass_Coeff_Ca
        )
    end
    function MaterialSetupOptions(dict::OrderedDict{Symbol,Any})
        new(
            get(dict, :stack_layers_bld, nothing),
            get(dict, :stack_layers_scale, [1.0, 1.0]),
            get(dict, :chord_scale, [1.0, 1.0]),
            get(dict, :thickness_scale, [1.0, 1.0]),
            get(dict, :AddedMass_Coeff_Ca, 0.0),
        )
    end
end

"""
    AeroSetupOptions

Contains the configuration for the aerodynamic model.

# Fields
- `rho::Float64`: The density of the air.
- `mu::Float64`: The dynamic viscosity of the air.
- `RPM::Float64`: The rotational speed of the blades.
- `Vinf::Float64`: The free stream velocity.
- `eta::Float64`: The tip speed ratio.
- `delta_t::Float64`: The time step for the aerodynamic model.
- `AD15hubR::Float64`: The hub radius for the AD15 model.
- `WindType::Int`: The type of wind model.
- `AeroModel::String`: The type of aerodynamic model.
- `DynamicStallModel::String`: The type of dynamic stall model.
- `numTS::Int`: The number of time steps for the aerodynamic model.
- `adi_lib::Union{Nothing,String}`: The path to the AeroDyn library.
- `adi_rootname::Union{Nothing,String}`: The root name for the AeroDyn files.
- `windINPfilename::Union{Nothing,String}`: The path to the wind input file.
- `ifw_libfile::Union{Nothing,String}`: The path to the inflow wind library.
- `ifw::Bool`: Whether to use the inflow wind model.
- `RPI::Bool`: Whether to use the RPI model.
- `Aero_AddedMass_Active::Bool`: Whether to use the added mass model.
- `Aero_RotAccel_Active::Bool`: Whether to use the rotational acceleration model.
- `Aero_Buoyancy_Active::Bool`: Whether to use the buoyancy model.
- `centrifugal_force_flag::Bool`: Whether to use the centrifugal force model.
"""
mutable struct AeroSetupOptions
    rho::Float64
    mu::Float64
    RPM::Float64
    Vinf::Float64
    eta::Float64
    delta_t::Float64
    AD15hubR::Float64
    WindType::Int
    AeroModel::String
    DynamicStallModel::String
    numTS::Int
    adi_lib::Union{Nothing,String}
    adi_rootname::Union{Nothing,String}
    windINPfilename::Union{Nothing,String}
    ifw_libfile::Union{Nothing,String}
    ifw::Bool
    RPI::Bool
    Aero_AddedMass_Active::Bool
    Aero_RotAccel_Active::Bool
    Aero_Buoyancy_Active::Bool
    centrifugal_force_flag::Bool
    AD15On::Bool

    function AeroSetupOptions(
        rho,
        mu,
        RPM,
        Vinf,
        eta,
        delta_t,
        AD15hubR,
        WindType,
        AeroModel,
        DynamicStallModel,
        numTS,
        adi_lib,
        adi_rootname,
        windINPfilename,
        ifw_libfile,
        ifw,
        RPI,
        Aero_AddedMass_Active,
        Aero_RotAccel_Active,
        Aero_Buoyancy_Active,
        centrifugal_force_flag,
    )
        new(
            rho,
            mu,
            RPM,
            Vinf,
            eta,
            delta_t,
            AD15hubR,
            WindType,
            AeroModel,
            DynamicStallModel,
            numTS,
            adi_lib,
            adi_rootname,
            windINPfilename,
            ifw_libfile,
            ifw,
            RPI,
            Aero_AddedMass_Active,
            Aero_RotAccel_Active,
            Aero_Buoyancy_Active,
            centrifugal_force_flag,
            false,
        )
    end
    function AeroSetupOptions(;
        rho::Float64 = 1.225,
        mu::Float64 = 1.7894e-5,
        RPM::Float64 = 1e-6,
        Vinf::Float64 = 25.0,
        eta::Float64 = 0.5,
        delta_t::Float64 = 0.01,
        AD15hubR::Float64 = 0.1,
        WindType::Int = 1,
        AeroModel::String = "DMS",
        DynamicStallModel::String = "BV",
        numTS::Int = 100,
        adi_lib::Union{Nothing,String} = nothing,
        adi_rootname::Union{Nothing,String} = nothing,
        windINPfilename::Union{Nothing,String} = nothing,
        ifw_libfile::Union{Nothing,String} = nothing,
        ifw::Bool = false,
        RPI::Bool = true,
        Aero_AddedMass_Active::Bool = false,
        Aero_RotAccel_Active::Bool = false,
        Aero_Buoyancy_Active::Bool = false,
        centrifugal_force_flag::Bool = false,
        AD15On::Bool = false,
    )
        new(
            rho,
            mu,
            RPM,
            Vinf,
            eta,
            delta_t,
            AD15hubR,
            WindType,
            AeroModel,
            DynamicStallModel,
            numTS,
            adi_lib,
            adi_rootname,
            windINPfilename,
            ifw_libfile,
            ifw,
            RPI,
            Aero_AddedMass_Active,
            Aero_RotAccel_Active,
            Aero_Buoyancy_Active,
            centrifugal_force_flag,
            AD15On,
        )
    end
    # Default constructor
    function AeroSetupOptions()
        new(
            1.225,        # rho
            1.7894e-5,    # mu
            1e-6,         # RPM
            25.0,         # Vinf
            0.5,          # eta
            0.01,         # delta_t
            0.1,          # AD15hubR
            1,            # WindType
            "DMS",        # AeroModel
            "BV",         # DynamicStallModel
            100,          # numTS
            "",           # adi_lib
            "",           # adi_rootname
            "",           # windINPfilename
            "",           # ifw_libfile
            false,        # ifw
            true,         # RPI
            false,        # Aero_AddedMass_Active
            false,        # Aero_RotAccel_Active
            false,        # Aero_Buoyancy_Active
            false,        # centrifugal_force_flag
            false,         # AD15On
        )
    end
    function AeroSetupOptions(dict::OrderedDict{Symbol,Any})
        new(
            get(dict, :rho, 1.225),
            get(dict, :mu, 1.7894e-5),
            get(dict, :RPM, 1e-6),
            get(dict, :Vinf, 25.0),
            get(dict, :eta, 0.5),
            get(dict, :delta_t, 0.01),
            get(dict, :AD15hubR, 0.1),
            get(dict, :WindType, 1),
            get(dict, :AeroModel, "DMS"),
            get(dict, :DynamicStallModel, "BV"),
            get(dict, :numTS, 100),
            get(dict, :adi_lib, ""),
            get(dict, :adi_rootname, ""),
            get(dict, :windINPfilename, ""),
            get(dict, :ifw_libfile, ""),
            get(dict, :ifw, false),
            get(dict, :RPI, true),
            get(dict, :Aero_AddedMass_Active, false),
            get(dict, :Aero_RotAccel_Active, false),
            get(dict, :Aero_Buoyancy_Active, false),
            get(dict, :centrifugal_force_flag, false),
            get(dict, :AD15On, false),
        )
    end
end

"""
    SetupOptions

Contains all setup options for the OWENS turbine model.

# Fields
- `mesh::MeshSetupOptions`: Mesh configuration options
- `tower::TowerSetupOptions`: Tower configuration options
- `blade::BladeSetupOptions`: Blade configuration options
- `material::MaterialSetupOptions`: Material configuration options
- `aero::AeroSetupOptions`: Aerodynamic configuration options
"""
mutable struct SetupOptions
    mesh::MeshSetupOptions
    tower::TowerSetupOptions
    blade::BladeSetupOptions
    material::MaterialSetupOptions
    aero::AeroSetupOptions

    function SetupOptions(;
        mesh::MeshSetupOptions = MeshSetupOptions(),
        tower::TowerSetupOptions = TowerSetupOptions(),
        blade::BladeSetupOptions = BladeSetupOptions(),
        material::MaterialSetupOptions = MaterialSetupOptions(),
        aero::AeroSetupOptions = AeroSetupOptions(),
    )
        new(mesh, tower, blade, material, aero)
    end

    # Default constructor
    function SetupOptions()
        new(
            MeshSetupOptions(),
            TowerSetupOptions(),
            BladeSetupOptions(),
            MaterialSetupOptions(),
            AeroSetupOptions(),
        )
    end

    # Constructor that takes a YAML dictionary
    function SetupOptions(dict_in::OrderedCollections.OrderedDict{Symbol,Any})
        # Extract each section from the dictionary, using default constructors if not present
        mesh =
            haskey(dict_in, :mesh) ? MeshSetupOptions(dict_in[:mesh]) : MeshSetupOptions()
        tower =
            haskey(dict_in, :tower) ? TowerSetupOptions(dict_in[:tower]) :
            TowerSetupOptions()
        blade =
            haskey(dict_in, :blade) ? BladeSetupOptions(dict_in[:blade]) :
            BladeSetupOptions()
        material =
            haskey(dict_in, :material) ? MaterialSetupOptions(dict_in[:material]) :
            MaterialSetupOptions()
        aero =
            haskey(dict_in, :aero) ? AeroSetupOptions(dict_in[:aero]) : AeroSetupOptions()

        new(mesh, tower, blade, material, aero)
    end
end

###
# Setup functions
###

struct MeshProperties
    mymesh::Any  # Type depends on mesh implementation
    myort::Any   # Type depends on orientation implementation
    myjoint::Any # Type depends on joint implementation
    AD15bldNdIdxRng::Matrix{Int}
    AD15bldElIdxRng::Matrix{Int}
    custom_mesh_outputs::Any  # Type depends on custom mesh implementation

    function MeshProperties(;
        mymesh,
        myort,
        myjoint,
        AD15bldNdIdxRng,
        AD15bldElIdxRng,
        custom_mesh_outputs,
    )
        new(mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs)
    end
end

# Component setup functions
function setup_mesh(
    mesh_config::MeshSetupOptions,
    blade_config::BladeSetupOptions,
    tower_config::TowerSetupOptions,
    aero_config::AeroSetupOptions,
    verbosity::Int64 = 1,
)
    mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs =
        OWENS.create_mesh_struts(;
            Htwr_base = tower_config.Htwr_base,
            Htwr_blds = tower_config.Htwr_blds,
            Hbld = blade_config.H,
            R = blade_config.R,
            AD15hubR = aero_config.AD15hubR,
            nblade = blade_config.B,
            ntelem = mesh_config.ntelem,
            nbelem = mesh_config.nbelem,
            nselem = mesh_config.nselem,
            strut_twr_mountpoint = tower_config.strut_twr_mountpoint,
            strut_bld_mountpoint = tower_config.strut_bld_mountpoint,
            bshapex = blade_config.shapeX,
            bshapez = blade_config.shapeZ,
            bshapey = blade_config.shapeY,
            angularOffset = tower_config.angularOffset,
            AD15_ccw = mesh_config.AD15_ccw,
            verbosity = verbosity,
            blade_joint_angle_Degrees = blade_config.blade_joint_angle_Degrees,
            connectBldTips2Twr = mesh_config.connectBldTips2Twr,
            strut_tower_joint_type = tower_config.strut_tower_joint_type,
            strut_blade_joint_type = blade_config.strut_blade_joint_type,
        )
    return MeshProperties(
        mymesh = mymesh,
        myort = myort,
        myjoint = myjoint,
        AD15bldNdIdxRng = AD15bldNdIdxRng,
        AD15bldElIdxRng = AD15bldElIdxRng,
        custom_mesh_outputs = custom_mesh_outputs,
    )
end

"""
Contains the sectional properties of the components.

# Fields
- `sectionPropsArray::Vector{Any}`: The sectional properties of the components.
- `stiff_array::Vector{Any}`: The stiffness of the components.
- `mass_array::Vector{Any}`: The mass of the components.
- `rotationalEffects::Vector{Float64}`: The rotational effects of the components.
"""
struct SectionalProperties
    sectionPropsArray::Vector{Any}
    stiff_array::Vector{Any}
    mass_array::Vector{Any}
    rotationalEffects::Vector{Float64}

    function SectionalProperties(
        sectionPropsArray,
        stiff_array,
        mass_array,
        rotationalEffects,
    )
        new(sectionPropsArray, stiff_array, mass_array, rotationalEffects)
    end
end

function setup_sectional_props(
    mesh_props,
    material_config,
    aero_config,
    tower_config,
    blade_config,
    path,
    NuMad_geom_xlscsv_file_intra_blade_cable,
    NuMad_mat_xlscsv_file_intra_blade_cable,
    NuMad_geom_xlscsv_file_guys,
    NuMad_mat_xlscsv_file_guys,
    verbosity,
)
    # Unpack the mesh properties
    mymesh = mesh_props.mymesh
    myort = mesh_props.myort
    myjoint = mesh_props.myjoint
    AD15bldNdIdxRng = mesh_props.AD15bldNdIdxRng
    AD15bldElIdxRng = mesh_props.AD15bldElIdxRng

    components,
    s2t_idx,
    intra_blade_cable_startidx_el,
    intra_blade_cable_endidx_el,
    topel_idx,
    s2b_idx = mesh_props.custom_mesh_outputs

    # Unpack the aero config
    rho = aero_config.rho

    # Unpack the material config
    AddedMass_Coeff_Ca = material_config.AddedMass_Coeff_Ca

    # Unpack the tower config
    NuMad_geom_xlscsv_file_twr = tower_config.NuMad_geom_xlscsv_file_twr
    NuMad_mat_xlscsv_file_twr = tower_config.NuMad_mat_xlscsv_file_twr

    # Unpack the blade config
    NuMad_geom_xlscsv_file_bld = blade_config.NuMad_geom_xlscsv_file_bld
    NuMad_mat_xlscsv_file_bld = blade_config.NuMad_mat_xlscsv_file_bld

    # Unpack the strut config
    NuMad_geom_xlscsv_file_strut = tower_config.NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut = tower_config.NuMad_mat_xlscsv_file_strut

    sectionPropsArray = []
    stiff_array = []
    mass_array = []
    rotationalEffects = ones(mymesh.numEl)

    for icomponent = 1:size(components)[1]
        if contains(components[icomponent].name, "tower")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_twr
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_twr
        elseif contains(components[icomponent].name, "blade")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_bld
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_bld
        elseif contains(components[icomponent].name, "strut")
            istrut = parse(Int, components[icomponent].name[end]) #This assumes that you have a numad file for each strut, and that you have 9 or fewer struts
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_strut[istrut]
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_strut
        elseif contains(components[icomponent].name, "intra_cable")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_intra_blade_cable
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_intra_blade_cable
        elseif contains(components[icomponent].name, "guy")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_guys
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_guys
            rotationalEffects[components[icomponent].elNumbers] .= 0.0 #turn rotational effects off for guy wires
        end

        # Skip components with missing layup or materials
        if isnothing(components[icomponent].input_layup) ||
           isnothing(components[icomponent].input_materials)
            continue
        end

        if verbosity>1
            println(
                "Calculating Sectional Properties for Component $icomponent $(components[icomponent].name)",
            )
        elseif verbosity>2
            println(
                "Using for the layup, $(components[icomponent].input_layup) and for the materials, $(components[icomponent].input_materials)",
            )
        end

        sectionPropsArray, stiff_array, mass_array, components[icomponent] =
            addSectionalPropertiesComponent!(
                sectionPropsArray,
                stiff_array,
                mass_array,
                components[icomponent],
                path,
                rho,
                AddedMass_Coeff_Ca;
                name = nothing,
            )

        components[icomponent].mass = OWENS.get_material_mass(
            components[icomponent].plyProps,
            components[icomponent].nuMadIn,
        )
        components[icomponent].cost = [
            "$name $(components[icomponent].mass[i]) kg, $(components[icomponent].plyProps.costs[i]) /kg: \$( $(components[icomponent].mass[i]*components[icomponent].plyProps.costs[i]) )"
            for (i, name) in enumerate(components[icomponent].plyProps.names)
        ]
    end
    return SectionalProperties(
        sectionPropsArray,
        stiff_array,
        mass_array,
        rotationalEffects,
    )
end

mutable struct AeroProperties
    aeroForcesAD::Any
    deformAeroAD::Any
    aeroForcesACDMS::Any
    deformAeroACDMS::Any

    function AeroProperties(aeroForcesAD, deformAeroAD, aeroForcesACDMS, deformAeroACDMS)
        new(aeroForcesAD, deformAeroAD, aeroForcesACDMS, deformAeroACDMS)
    end
end

function setup_aerodynamic_model(
    blade_config::BladeSetupOptions,
    aero_config::AeroSetupOptions,
    tower_config::TowerSetupOptions,
    mesh_config::MeshSetupOptions,
    mesh_props::MeshProperties,
    components::Vector{OWENS.Component},
    path::String,
    myel::OWENSFEA.El,
    verbosity::Int64 = 1,
)

    blade_component_index = findfirst(s -> s.name[1:(end-1)] == "blade", components) #This assumes 9 or fewer blades, not including struts
    numadIn_bld = components[blade_component_index].nuMadIn

    strut_component_index = findall(s -> contains(s.name, "strut"), components) #This assumes 9 or fewer blades, not including struts
    numadIn_strut = [
        components[strut_component_index1].nuMadIn for
        strut_component_index1 in strut_component_index[1:blade_config.B:end]
    ]
    # Initialize aerodynamic variables to nothing
    # TODO: Add types to the variables for proper allocation
    aeroForcesAD = nothing
    deformAeroAD = nothing
    aeroForcesACDMS = nothing
    deformAeroACDMS = nothing

    # Unpack the config
    AD15On = aero_config.AeroModel == "AD"
    WindType = aero_config.WindType
    windINPfilename = aero_config.windINPfilename
    ifw_libfile = aero_config.ifw_libfile
    AeroModel = aero_config.AeroModel
    DynamicStallModel = aero_config.DynamicStallModel
    Aero_AddedMass_Active = aero_config.Aero_AddedMass_Active
    Aero_RotAccel_Active = aero_config.Aero_RotAccel_Active
    Aero_Buoyancy_Active = aero_config.Aero_Buoyancy_Active
    centrifugal_force_flag = aero_config.centrifugal_force_flag
    ntheta = mesh_config.ntheta
    Nslices = mesh_config.Nslices
    meshtype = mesh_config.meshtype
    RPI = aero_config.RPI
    rho = aero_config.rho
    Vinf = aero_config.Vinf
    mu = aero_config.mu
    eta = aero_config.eta
    ifw = aero_config.ifw
    Nbld = blade_config.B
    mymesh = mesh_props.mymesh
    myort = mesh_props.myort
    myjoint = mesh_props.myjoint
    shapeX = blade_config.shapeX
    shapeZ = blade_config.shapeZ
    shapeY = blade_config.shapeY
    AD15bldNdIdxRng = mesh_props.AD15bldNdIdxRng
    AD15bldElIdxRng = mesh_props.AD15bldElIdxRng
    H = blade_config.H
    Htwr_base = tower_config.Htwr_base

    Nstrutperbld = length(tower_config.strut_twr_mountpoint)

    # Calculate tsr locally
    omega = aero_config.RPM / 60 * 2 * pi
    tsr = omega * blade_config.R / aero_config.Vinf

    # Set up AeroDyn if used
    # Here we create AeroDyn the files, first by specifying the names, then by creating the files, TODO: hook up the direct sectionPropsArray_str
    # Then by initializing AeroDyn and grabbing the backend functionality with a function handle
    if AD15On
        ad_input_file = "$path/ADInputFile_SingleTurbine2.dat"
        ifw_input_file = "$path/IW2.dat"
        OLAF_filename = "$path/OLAF2.dat"

        NumADBldNds = NumADStrutNds = 10

        bldchord_spl = OWENS.safeakima(
            numadIn_bld.span ./ maximum(numadIn_bld.span),
            numadIn_bld.chord,
            LinRange(0, 1, NumADBldNds),
        )

        # Discretely assign the blade airfoils based on the next closest neighbor
        bld_airfoil_filenames = fill("nothing", NumADBldNds) #TODO: cable drag?
        for (ispan_numad, span_numad) in
            enumerate(numadIn_bld.span ./ maximum(numadIn_bld.span))
            for (ispan, span_slices) in enumerate(collect(LinRange(0, 1, NumADBldNds)))
                if bld_airfoil_filenames[ispan]=="nothing" && span_slices<=span_numad
                    bld_airfoil_filenames[ispan] = "$(numadIn_bld.airfoil[ispan_numad]).dat"
                end
            end
        end

        if meshtype == "ARCUS"
            blade_filenames = ["$path/blade$i.dat" for i = 1:Nbld]
            blade_chords = [bldchord_spl for i = 1:Nbld]
            blade_Nnodes = [NumADBldNds for i = 1:Nbld]
            airfoil_filenames = [bld_airfoil_filenames for i = 1:Nbld]

        else
            blade_filenames = ["$path/blade$i.dat" for i = 1:Nbld]
            blade_chords = [bldchord_spl for i = 1:Nbld]
            blade_Nnodes = [NumADBldNds for i = 1:Nbld]
            airfoil_filenames =
                collect(Iterators.flatten([bld_airfoil_filenames for i = 1:Nbld]))

            for istrut = 1:Nstrutperbld
                strutchord_spl = OWENS.safeakima(
                    numadIn_strut[istrut].span ./ maximum(numadIn_strut[istrut].span),
                    numadIn_strut[istrut].chord,
                    LinRange(0, 1, NumADStrutNds),
                )
                for ibld = 1:Nbld
                    blade_filenames = [blade_filenames; "$path/strut$(istrut)_bld$ibld.dat"]
                    blade_chords = [blade_chords; [strutchord_spl]]
                    blade_Nnodes = [blade_Nnodes; NumADStrutNds]

                    # Discretely assign the strut airfoils based on the next closest neighbor
                    strut_airfoil_filenames = fill("nothing", NumADStrutNds)
                    for (ispan_numad, span_numad) in enumerate(
                        numadIn_strut[istrut].span ./ maximum(numadIn_strut[istrut].span),
                    )
                        for (ispan, span_slices) in
                            enumerate(collect(LinRange(0, 1, NumADBldNds)))
                            if strut_airfoil_filenames[ispan]=="nothing" &&
                               span_slices<=span_numad
                                strut_airfoil_filenames[ispan] = "$(numadIn_strut[istrut].airfoil[ispan_numad]).dat"
                            end
                        end
                    end

                    airfoil_filenames = [airfoil_filenames; strut_airfoil_filenames]

                end
            end
        end

        OWENSOpenFASTWrappers.writeADinputFile(
            ad_input_file,
            blade_filenames,
            airfoil_filenames,
            OLAF_filename;
            rho,
        )

        NumADBody = length(AD15bldNdIdxRng[:, 1])
        bld_len = zeros(NumADBody)
        # bladefileissaved = false
        for (iADBody, filename) in enumerate(blade_filenames)
            strt_idx = AD15bldNdIdxRng[iADBody, 1]
            end_idx = AD15bldNdIdxRng[iADBody, 2]
            if end_idx<strt_idx
                tmp_end = end_idx
                end_idx = strt_idx
                strt_idx = tmp_end
            end

            #Get the blade length
            x1 = mymesh.x[strt_idx]
            x2 = mymesh.x[end_idx]
            y1 = mymesh.y[strt_idx]
            y2 = mymesh.y[end_idx]
            z1 = mymesh.z[strt_idx]
            z2 = mymesh.z[end_idx]
            bld_len[iADBody] = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)

            #Get the blade shape
            ADshapeZ = collect(LinRange(0, H, NumADBldNds))
            xmesh = mymesh.x[strt_idx:end_idx]
            ymesh = mymesh.y[strt_idx:end_idx]
            ADshapeX = sqrt.(xmesh .^ 2 .+ ymesh .^ 2)
            ADshapeX .-= ADshapeX[1] #get it starting at zero #TODO: make robust for blades that don't start at 0
            ADshapeXspl =
                OWENS.safeakima(LinRange(0, H, length(ADshapeX)), ADshapeX, ADshapeZ)

            if iADBody<=Nbld #&& !bladefileissaved#Note that the blades can be curved and are assumed to be oriented vertically
                # bladefileissaved = true
                BlSpn0=ADshapeZ
                BlCrvAC0=ADshapeXspl

                bladeangle = (iADBody-1)*2.0*pi/Nbld + tower_config.angularOffset #TODO: pitch offset and twist offset that isn't from the helical

                BlSpn = ADshapeZ
                blade_twist = atan.(xmesh, ymesh) .- bladeangle

                #TODO: reevalueate these equations and make sure they are robust against varying designs
                # BlCrvACinput = xmesh.*cos(bladeangle) .+ ymesh.*sin(bladeangle)
                # BlCrvACinput = BlCrvACinput .- BlCrvACinput[1]
                # BlCrvAC = OWENS.safeakima(LinRange(0,H,length(BlCrvACinput)),BlCrvACinput,ADshapeZ)

                # BlSwpACinput = -ymesh.*cos(bladeangle) .+ xmesh.*sin(bladeangle)
                # BlSwpACinput = BlSwpACinput .- BlSwpACinput[1]
                # BlSwpAC = OWENS.safeakima(LinRange(0,H,length(BlSwpACinput)),BlSwpACinput,ADshapeZ)

                # BlCrvAng = zeros(blade_Nnodes[iADBody])

                # BlTwistinput =(blade_twist.-blade_twist[1])*180/pi
                # BlTwist = OWENS.safeakima(LinRange(0,H,length(BlTwistinput)),BlTwistinput,ADshapeZ)


                BlCrvACinput = -ymesh .* sin(bladeangle) .+ xmesh .* cos(bladeangle)
                BlCrvACinput = BlCrvACinput .- BlCrvACinput[1]
                BlSwpAC =
                    -safeakima(LinRange(0, H, length(BlCrvACinput)), BlCrvACinput, ADshapeZ)

                BlSwpACinput = xmesh .* sin(bladeangle) .+ ymesh .* cos(bladeangle)
                BlSwpACinput = BlSwpACinput .- BlSwpACinput[1]
                BlCrvAC =
                    safeakima(LinRange(0, H, length(BlSwpACinput)), BlSwpACinput, ADshapeZ)

                BlCrvAng = zeros(blade_Nnodes[iADBody])

                BlTwistinput = (blade_twist .- blade_twist[2])*180/pi
                BlTwist =
                    safeakima(LinRange(0, H, length(BlTwistinput)), BlTwistinput, ADshapeZ)

                BlChord=blade_chords[iADBody]

                BlAFID=collect(((iADBody-1)*NumADBldNds+1):(iADBody*NumADBldNds))

            elseif iADBody>Nbld # while the arms/struts are assumed to be straight and are oriented by the mesh angle
                BlSpn=collect(LinRange(0, bld_len[iADBody], blade_Nnodes[iADBody]))
                BlCrvAC=zeros(blade_Nnodes[iADBody])
                BlSwpAC=zeros(blade_Nnodes[iADBody])
                BlCrvAng=zeros(blade_Nnodes[iADBody])
                BlTwist=zeros(blade_Nnodes[iADBody])
                BlChord=blade_chords[iADBody]
                BlAFID=collect(((iADBody-1)*NumADStrutNds+1):(iADBody*NumADStrutNds))
            end
            OWENSOpenFASTWrappers.writeADbladeFile(
                filename;
                NumBlNds = blade_Nnodes[iADBody],
                BlSpn,
                BlCrvAC,
                BlSwpAC,
                BlCrvAng,
                BlTwist,
                BlChord,
                BlAFID,
            )
        end

        OWENSOpenFASTWrappers.writeOLAFfile(OLAF_filename; nNWPanel = 200, nFWPanels = 10)

        OWENSOpenFASTWrappers.writeIWfile(Vinf, ifw_input_file; WindType, windINPfilename)

        OWENSOpenFASTWrappers.setupTurb(
            aero_config.adi_lib,
            ad_input_file,
            ifw_input_file,
            aero_config.adi_rootname,
            [shapeX],
            [shapeZ],
            [Nbld],
            [Htwr_base],
            [mymesh],
            [myort],
            [AD15bldNdIdxRng],
            [AD15bldElIdxRng];
            rho = rho,
            adi_dt = aero_config.delta_t,
            adi_tmax = aero_config.numTS*aero_config.delta_t,
            omega = [omega],
            adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
            adi_DT_Outs = aero_config.delta_t,   # output frequency
            numTurbines = 1,
            refPos = [[0, 0, 0]],
            hubPos = [[0, 0, 0.0]],
            hubAngle = [[0, 0, 0]],
            nacPos = [[0, 0, 0]],
            adi_nstrut = [Nstrutperbld],
            adi_debug = 0,
            isHAWT = false,     # true for HAWT, false for crossflow or VAWT
        )

        aeroForcesAD(t, azi) = OWENS.mapAD15(
            t,
            azi,
            [mymesh],
            OWENSOpenFASTWrappers.advanceAD15;
            alwaysrecalc = true,
            verbosity,
        )
        deformAeroAD=OWENSOpenFASTWrappers.deformAD15

    else
        #########################################
        ### Set up aero forces
        #########################################
        #### translate from blade span to blade height between the numad definition and the vertical slice positions
        #### First get the angles from the overall geometry npoints and go to the numad npoints
        delta_xs = shapeX[2:end] - shapeX[1:(end-1)]
        delta_zs = shapeZ[2:end] - shapeZ[1:(end-1)]
        delta3D = atan.(delta_xs ./ delta_zs)
        delta3D_spl = OWENS.safeakima(
            shapeZ[1:(end-1)] ./ maximum(shapeZ[1:(end-1)]),
            delta3D,
            LinRange(0, 1, length(numadIn_bld.span)-1),
        )
        #### now convert the numad span to a height
        bld_height_numad =
            cumsum(diff(numadIn_bld.span) .* (1.0 .- abs.(sin.(delta3D_spl))))
        bld_height_numad_unit = bld_height_numad ./ maximum(bld_height_numad)
        #### now we can use it to access the numad data 
        chord = OWENS.safeakima(
            bld_height_numad_unit,
            numadIn_bld.chord,
            LinRange(bld_height_numad_unit[1], 1, Nslices),
        )
        airfoils = fill("nothing", Nslices)

        twist = OWENS.safeakima(
            bld_height_numad_unit,
            numadIn_bld.twist_d .* pi/180,
            LinRange(bld_height_numad_unit[1], 1, Nslices),
        )

        # Discretely assign the airfoils
        for (iheight_numad, height_numad) in enumerate(bld_height_numad_unit)
            for (iheight, height_slices) in enumerate(collect(LinRange(0, 1, Nslices)))
                if airfoils[iheight]=="nothing" && height_slices<=height_numad
                    airfoils[iheight] = "$(numadIn_bld.airfoil[iheight_numad]).dat"
                end
            end
        end

        # Map the element wise mass to the input aero shape
        mass_bld = []
        for component in components
            if contains(component.name, "blade1")
                mass_bld = component.mass_matrix
                break
            end
        end

        rhoA_el = [mass_bld[i][1, 1] for i = 1:length(mass_bld)]
        if AD15bldNdIdxRng[1, 2]<AD15bldNdIdxRng[1, 1]
            bld_z_node = mymesh.z[AD15bldNdIdxRng[1, 2]:AD15bldNdIdxRng[1, 1]]
        else
            bld_z_node = mymesh.z[AD15bldNdIdxRng[1, 1]:AD15bldNdIdxRng[1, 2]]
        end
        bld_z_el =
            bld_z_node[1:(end-1)] .+ (bld_z_node[2:end] .- bld_z_node[1:(end-1)]) ./ 2
        rhoA_in = FLOWMath.akima(bld_z_el, rhoA_el, shapeZ)

        OWENSAero.setupTurb(
            shapeX,
            shapeZ,
            Nbld,
            chord,
            tsr,
            Vinf;
            AeroModel,
            DynamicStallModel,
            afname = airfoils,
            bld_y = shapeY,
            rho,
            twist, #TODO: verify twist is in same direction
            mu,
            eta,
            ifw, #TODO: propogate WindType
            turbsim_filename = windINPfilename,
            ifw_libfile,
            tau = [1e-5, 1e-5],
            Aero_AddedMass_Active,
            Aero_RotAccel_Active,
            Aero_Buoyancy_Active,
            centrifugal_force_flag,
            ntheta,
            Nslices,
            RPI,
            rhoA_in,
        )

        aeroForcesACDMS(t, azi) = OWENS.mapACDMS(
            t,
            azi,
            mymesh,
            myel,
            OWENSAero.AdvanceTurbineInterpolate;
            alwaysrecalc = true,
        )
        deformAeroACDMS = OWENSAero.deformTurb
    end

    # Always return all four values, even if some are nothing
    return AeroProperties(aeroForcesAD, deformAeroAD, aeroForcesACDMS, deformAeroACDMS)
end

function addSectionalPropertiesComponent!(
    sectionPropsArray,
    stiff_array,
    mass_array,
    component,
    path,
    rho,
    AddedMass_Coeff_Ca;
    name = nothing,
)

    input_layup = component.input_layup
    input_materials = component.input_materials
    start_el_idx = component.elNumbers[1]
    end_el_idx = component.elNumbers[end]

    nElem = end_el_idx-start_el_idx+1+1
    subsection = nothing
    if !isnothing(input_layup)
        if contains(component.name, "tower")
            section=:tower
        elseif contains(component.name, "blade")
            section=:blade
        elseif contains(component.name, "strut")
            section=:struts
            if typeof(input_layup) == OrderedCollections.OrderedDict{Symbol,Any}
                if length(input_layup[:components][:struts])==1
                    subsection = 1
                else
                    subsection = parse(Int, component.name[end]) #This assumes that you have a numad file for each strut, and that you have 9 or fewer struts
                end
            else
                subsection = nothing
            end
        elseif contains(component.name, "intra_cable")
            section=:intra_cable
        elseif contains(component.name, "guy")
            section=:guy
        end
        numadIn = OWENS.readNuMadGeomCSV(input_layup; section, subsection) #note that this function will run either the numad input or the windio input depending on the input type
    else
        @error "Input data or numad file must be defined for this method"
    end

    #### Add the full path
    for (i, airfoil) in enumerate(numadIn.airfoil)
        numadIn.airfoil[i] = "$path/airfoils/$airfoil"
    end

    # Here is where the material properties for the tower are either read in from the file, or directly input
    if !isnothing(input_materials)
        plyprops = OWENS.readNuMadMaterialsCSV(input_materials)
    else
        @error "Input data or numad file must be defined for this method"
    end

    # Then this is where precomp.jl is called to get first the precomp outputs, then formatting those into the OWENS format, and then in the GXBeam.jl format for if GXBeam is used as the structural solver.
    precompoutput, precompinput, lam_U, lam_L, lam_W =
        OWENS.getOWENSPreCompOutput(numadIn; plyprops = plyprops)
    sectionPropsArray_component = OWENS.getSectPropsFromOWENSPreComp(
        LinRange(0, 1, nElem),
        numadIn,
        precompoutput;
        precompinputs = precompinput,
        fluid_density = rho,
        AddedMass_Coeff_Ca,
    )
    stiff, mass = OWENS.getSectPropsFromOWENSPreComp(
        LinRange(0, 1, nElem),
        numadIn,
        precompoutput;
        GX = true,
        precompinputs = precompinput,
        fluid_density = rho,
        AddedMass_Coeff_Ca,
    )

    sectionPropsArray = [sectionPropsArray; sectionPropsArray_component]
    stiff_array = [stiff_array; stiff]
    mass_array = [mass_array; mass]


    component.lam_U = lam_U
    component.lam_L = lam_L
    component.lam_W = lam_W
    component.preCompInput = precompinput
    component.preCompOutput = precompoutput
    component.plyProps = plyprops
    component.nuMadIn = numadIn
    component.sectionProps = sectionPropsArray_component
    component.stiff_matrix = stiff
    component.mass_matrix = mass

    return sectionPropsArray, stiff_array, mass_array, component
end

function get_material_mass(
    plyprops_in,
    numadIn;
    int_start = numadIn.span[1],
    int_stop = numadIn.span[end],
)
    # Get Relative contribution to mass by setting all but one of the materials to zero.
    mass_component_material = zeros(length(plyprops_in.names))
    for imat = 1:length(plyprops_in.names)
        # Initialize array since it is immutable
        plies = Array{Composites.Material}(undef, length(plyprops_in.names))

        # Fill it in, setting all the rho's to zero except the one matching imat
        for imat2 = 1:length(plyprops_in.names)
            if imat != imat2
                plies[imat2] = Composites.Material(
                    plyprops_in.plies[imat2].e1,
                    plyprops_in.plies[imat2].e2,
                    plyprops_in.plies[imat2].g12,
                    plyprops_in.plies[imat2].nu12,
                    0.0, #rho
                    plyprops_in.plies[imat2].xt,
                    plyprops_in.plies[imat2].xc,
                    plyprops_in.plies[imat2].yt,
                    plyprops_in.plies[imat2].yc,
                    plyprops_in.plies[imat2].s,
                    plyprops_in.plies[imat2].t,
                )
            else
                plies[imat2] = plyprops_in.plies[imat2]
            end
        end

        plyprops = plyproperties(plyprops_in.names, plies)
        # Get the precomp output
        precompoutput, _, _, _, _ = getOWENSPreCompOutput(numadIn; plyprops)
        mass_array = [precompoutput[iter].mass for iter = 1:length(precompoutput)]
        # Spline and integrate that output across the span
        mass_spl = FLOWMath.Akima(numadIn.span, mass_array)
        mass_component_material[imat], error =
            QuadGK.quadgk(mass_spl, int_start, int_stop, atol = 1e-10)

    end
    return mass_component_material
end

"""
    preprocess_windio_setup(modelopt, windio, path)

Creates and populates a SetupOptions instance from the windio and modelopt inputs.

# Arguments
- `modelopt`: The modeling options containing configuration parameters
- `windio`: The windio data containing turbine geometry and properties
- `path`: The base path for file operations

# Returns
- `SetupOptions`: A populated SetupOptions instance containing all configuration parameters
"""
function preprocess_windio_setup(modelopt, windio, path)
    # Extract basic geometry parameters
    number_of_blades = windio[:assembly][:number_of_blades]
    hub_height = windio[:assembly][:hub_height]

    # Extract blade geometry
    blade_x = windio[:components][:blade][:outer_shape_bem][:reference_axis][:x][:values]
    blade_y = windio[:components][:blade][:outer_shape_bem][:reference_axis][:y][:values]
    blade_z = windio[:components][:blade][:outer_shape_bem][:reference_axis][:z][:values]
    tower_z = windio[:components][:tower][:outer_shape_bem][:reference_axis][:z][:values]

    # Calculate derived geometry
    Blade_Height = maximum(blade_z)
    Blade_Radius = maximum(sqrt.(blade_x .^ 2 .+ blade_y .^ 2))
    Htwr_base = hub_height - Blade_Height/2
    Htwr_blds = maximum(tower_z) - Htwr_base

    # Extract strut parameters
    tower_strut_connection = windio[:components][:struts][1][:mountfraction_tower]
    blade_strut_connection = windio[:components][:struts][1][:mountfraction_blade]

    # Extract environment parameters
    air_density = windio[:environment][:air_density]
    air_dyn_viscosity = windio[:environment][:air_dyn_viscosity]

    # Create mesh options
    mesh_opts = MeshSetupOptions(
        Nslices = modelopt.OWENSAero_Options.Nslices,
        ntheta = modelopt.OWENSAero_Options.ntheta,
        ntelem = modelopt.Mesh_Options.ntelem,
        nbelem = modelopt.Mesh_Options.nbelem,
        ncelem = modelopt.Mesh_Options.ncelem,
        nselem = modelopt.Mesh_Options.nselem,
        meshtype = modelopt.Mesh_Options.turbineType,
        custommesh = nothing,
        connectBldTips2Twr = modelopt.Mesh_Options.cables_connected_to_blade_base,
        AD15_ccw = true,
    )

    # Create tower options
    tower_opts = TowerSetupOptions(
        Htwr_base = Htwr_base,
        Htwr_blds = Htwr_blds,
        strut_twr_mountpoint = tower_strut_connection,
        strut_bld_mountpoint = blade_strut_connection,
        joint_type = modelopt.Mesh_Options.joint_type,
        c_mount_ratio = modelopt.Mesh_Options.c_mount_ratio,
        angularOffset = modelopt.Mesh_Options.angularOffset,
        NuMad_geom_xlscsv_file_twr = windio,
        NuMad_mat_xlscsv_file_twr = windio,
        NuMad_geom_xlscsv_file_strut = windio,
        NuMad_mat_xlscsv_file_strut = windio,
        strut_tower_joint_type = 2,
    )

    # Create blade options
    blade_opts = BladeSetupOptions(
        B = number_of_blades,
        H = Blade_Height,
        R = Blade_Radius,
        shapeZ = blade_z,
        shapeX = blade_x,
        shapeY = blade_y,
        NuMad_geom_xlscsv_file_bld = windio,
        NuMad_mat_xlscsv_file_bld = windio,
        strut_blade_joint_type = 0,
        blade_joint_angle_Degrees = 0.0,
    )

    # Create material options
    material_opts = MaterialSetupOptions(
        stack_layers_bld = nothing,
        stack_layers_scale = [1.0, 1.0],
        chord_scale = [1.0, 1.0],
        thickness_scale = [1.0, 1.0],
        AddedMass_Coeff_Ca = modelopt.OWENSFEA_Options.AddedMass_Coeff_Ca,
    )

    # Create aero options
    aero_opts = AeroSetupOptions(
        rho = air_density,
        mu = air_dyn_viscosity,
        RPM = modelopt.OWENS_Options.Prescribed_RPM_RPM_controlpoints[1],
        Vinf = modelopt.OWENS_Options.Prescribed_Vinf_Vinf_controlpoints[1],
        eta = windio[:components][:blade][:outer_shape_bem][:blade_mountpoint],
        delta_t = modelopt.OWENS_Options.delta_t,
        AD15hubR = modelopt.Mesh_Options.AD15hubR,
        WindType = modelopt.OWENSOpenFASTWrappers_Options.WindType,
        AeroModel = modelopt.OWENS_Options.AeroModel,
        DynamicStallModel = modelopt.OWENSAero_Options.DynamicStallModel,
        numTS = modelopt.OWENS_Options.numTS,
        adi_lib = modelopt.OWENSOpenFASTWrappers_Options.adi_lib == "nothing" ?
                  nothing : modelopt.OWENSOpenFASTWrappers_Options.adi_lib,
        adi_rootname = "$(path)$(modelopt.OWENSOpenFASTWrappers_Options.adi_rootname)",
        windINPfilename = "$(path)$(modelopt.OWENSOpenFASTWrappers_Options.windINPfilename)",
        ifw_libfile = modelopt.OWENSOpenFASTWrappers_Options.ifw_libfile == "nothing" ?
                      nothing : modelopt.OWENSOpenFASTWrappers_Options.ifw_libfile,
        ifw = modelopt.OWENSAero_Options.ifw,
        RPI = modelopt.OWENSAero_Options.RPI,
        Aero_AddedMass_Active = modelopt.OWENSAero_Options.Aero_AddedMass_Active,
        Aero_RotAccel_Active = modelopt.OWENSAero_Options.Aero_RotAccel_Active,
        Aero_Buoyancy_Active = modelopt.OWENSAero_Options.Aero_Buoyancy_Active,
        centrifugal_force_flag = false,
        AD15On = modelopt.OWENS_Options.AeroModel == "AD",
    )

    # Create and return the complete SetupOptions instance
    return SetupOptions(
        mesh = mesh_opts,
        tower = tower_opts,
        blade = blade_opts,
        material = material_opts,
        aero = aero_opts,
    )
end
