using OrderedCollections
using Test
import OWENS

@testset "Setup option constructors and dictionaries" begin
    mesh = OWENS.MeshSetupOptions(;
        Nslices = 8,
        ntheta = 12,
        ntelem = 4,
        nbelem = 5,
        ncelem = 2,
        nselem = 3,
        meshtype = "H-VAWT",
        connectBldTips2Twr = true,
        AD15_ccw = false,
    )
    @test mesh.Nslices == 8
    @test mesh.ntheta == 12
    @test mesh.ntelem == 4
    @test mesh.nbelem == 5
    @test mesh.ncelem == 2
    @test mesh.nselem == 3
    @test mesh.meshtype == "H-VAWT"
    @test mesh.custommesh === nothing
    @test mesh.connectBldTips2Twr === true
    @test mesh.AD15_ccw === false

    mesh_dict = OWENS.MeshSetupOptions(
        OrderedDict{Symbol,Any}(
            :Nslices => 9,
            :ntheta => 18,
            :meshtype => "ARCUS",
            :connectBldTips2Twr => true,
        ),
    )
    @test mesh_dict.Nslices == 9
    @test mesh_dict.ntheta == 18
    @test mesh_dict.ntelem == 10
    @test mesh_dict.nbelem == 60
    @test mesh_dict.meshtype == "ARCUS"
    @test mesh_dict.connectBldTips2Twr === true
    @test mesh_dict.AD15_ccw === true

    tower = OWENS.TowerSetupOptions(;
        Htwr_base = 3.5,
        Htwr_blds = 12.0,
        strut_twr_mountpoint = [0.2, 0.8],
        strut_bld_mountpoint = [0.25, 0.75],
        joint_type = 4,
        c_mount_ratio = 0.08,
        angularOffset = 0.25,
        NuMad_geom_xlscsv_file_twr = "tower-geom.csv",
        NuMad_mat_xlscsv_file_twr = "tower-mat.csv",
        NuMad_geom_xlscsv_file_strut = ["strut1.csv", "strut2.csv"],
        NuMad_mat_xlscsv_file_strut = "strut-mat.csv",
        strut_tower_joint_type = 3,
    )
    @test tower.Htwr_base == 3.5
    @test tower.Htwr_blds == 12.0
    @test tower.strut_twr_mountpoint == [0.2, 0.8]
    @test tower.strut_bld_mountpoint == [0.25, 0.75]
    @test tower.joint_type == 4
    @test tower.c_mount_ratio == 0.08
    @test tower.angularOffset == 0.25
    @test tower.NuMad_geom_xlscsv_file_twr == "tower-geom.csv"
    @test tower.NuMad_mat_xlscsv_file_strut == "strut-mat.csv"
    @test tower.strut_tower_joint_type == 3

    tower_dict = OWENS.TowerSetupOptions(
        OrderedDict{Symbol,Any}(:Htwr_base => 4.0, :strut_tower_joint_type => 5),
    )
    @test tower_dict.Htwr_base == 4.0
    @test tower_dict.Htwr_blds == 5.0
    @test tower_dict.strut_twr_mountpoint == [0.25, 0.75]
    @test tower_dict.strut_tower_joint_type == 5

    blade = OWENS.BladeSetupOptions(;
        B = 2,
        H = 14.0,
        R = 3.2,
        shapeZ = [0.0, 7.0, 14.0],
        shapeX = [0.0, 3.2, 0.0],
        shapeY = [0.0, 0.1, 0.0],
        NuMad_geom_xlscsv_file_bld = "blade-geom.csv",
        NuMad_mat_xlscsv_file_bld = "blade-mat.csv",
        strut_blade_joint_type = 2,
        blade_joint_angle_Degrees = 12.5,
    )
    @test blade.B == 2
    @test blade.H == 14.0
    @test blade.R == 3.2
    @test blade.shapeZ == [0.0, 7.0, 14.0]
    @test blade.shapeX == [0.0, 3.2, 0.0]
    @test blade.shapeY == [0.0, 0.1, 0.0]
    @test blade.NuMad_geom_xlscsv_file_bld == "blade-geom.csv"
    @test blade.NuMad_mat_xlscsv_file_bld == "blade-mat.csv"
    @test blade.strut_blade_joint_type == 2
    @test blade.blade_joint_angle_Degrees == 12.5

    blade_dict = OWENS.BladeSetupOptions(
        OrderedDict{Symbol,Any}(:B => 4, :shapeY => [0.0, 0.0, 0.0]),
    )
    @test blade_dict.B == 4
    @test blade_dict.H == 5.0
    @test length(blade_dict.shapeZ) == 31
    @test blade_dict.shapeY == [0.0, 0.0, 0.0]

    material = OWENS.MaterialSetupOptions(;
        stack_layers_bld = [0.0 45.0; -45.0 90.0],
        stack_layers_scale = [1.0, 1.2],
        chord_scale = [0.9, 1.1],
        thickness_scale = [0.8, 1.3],
        AddedMass_Coeff_Ca = 0.75,
        sectional_property_source = "gxbeam",
    )
    @test material.stack_layers_bld == [0.0 45.0; -45.0 90.0]
    @test material.stack_layers_scale == [1.0, 1.2]
    @test material.chord_scale == [0.9, 1.1]
    @test material.thickness_scale == [0.8, 1.3]
    @test material.AddedMass_Coeff_Ca == 0.75
    @test material.sectional_property_source == :gxbeam

    material_dict = OWENS.MaterialSetupOptions(
        OrderedDict{Symbol,Any}(:AddedMass_Coeff_Ca => 1.1, :sectional_property_source => "precomp"),
    )
    @test material_dict.stack_layers_bld === nothing
    @test material_dict.AddedMass_Coeff_Ca == 1.1
    @test material_dict.sectional_property_source == :precomp

    aero = OWENS.AeroSetupOptions(;
        rho = 1.2041,
        mu = 1.8e-5,
        RPM = 42.0,
        Vinf = 11.5,
        eta = 0.42,
        delta_t = 0.025,
        AD15hubR = 0.3,
        WindType = 3,
        AeroModel = "AC",
        DynamicStallModel = "LB",
        numTS = 123,
        adi_lib = "adi-lib",
        adi_rootname = "adi-root",
        windINPfilename = "wind.bts",
        ifw_libfile = "ifw-lib",
        ifw = true,
        RPI = false,
        Aero_AddedMass_Active = true,
        Aero_RotAccel_Active = true,
        Aero_Buoyancy_Active = true,
        centrifugal_force_flag = true,
        AD15On = true,
    )
    @test aero.rho == 1.2041
    @test aero.mu == 1.8e-5
    @test aero.RPM == 42.0
    @test aero.Vinf == 11.5
    @test aero.eta == 0.42
    @test aero.delta_t == 0.025
    @test aero.AD15hubR == 0.3
    @test aero.WindType == 3
    @test aero.AeroModel == "AC"
    @test aero.DynamicStallModel == "LB"
    @test aero.numTS == 123
    @test aero.adi_lib == "adi-lib"
    @test aero.adi_rootname == "adi-root"
    @test aero.windINPfilename == "wind.bts"
    @test aero.ifw_libfile == "ifw-lib"
    @test aero.ifw === true
    @test aero.RPI === false
    @test aero.Aero_AddedMass_Active === true
    @test aero.Aero_RotAccel_Active === true
    @test aero.Aero_Buoyancy_Active === true
    @test aero.centrifugal_force_flag === true
    @test aero.AD15On === true

    aero_dict = OWENS.AeroSetupOptions(
        OrderedDict{Symbol,Any}(:RPM => 36.0, :AeroModel => "AD", :AD15On => true),
    )
    @test aero_dict.RPM == 36.0
    @test aero_dict.AeroModel == "AD"
    @test aero_dict.adi_rootname == "/aerodyn"
    @test aero_dict.AD15On === true

    setup = OWENS.SetupOptions(; mesh, tower, blade, material, aero)
    @test setup.mesh === mesh
    @test setup.tower === tower
    @test setup.blade === blade
    @test setup.material === material
    @test setup.aero === aero

    setup_dict = OWENS.SetupOptions(
        OrderedDict{Symbol,Any}(
            :mesh => OrderedDict{Symbol,Any}(:Nslices => 5),
            :material => OrderedDict{Symbol,Any}(:sectional_property_source => "gxbeam"),
            :aero => OrderedDict{Symbol,Any}(:Vinf => 9.5),
        ),
    )
    @test setup_dict.mesh.Nslices == 5
    @test setup_dict.tower.Htwr_base == 2.0
    @test setup_dict.blade.B == 3
    @test setup_dict.material.sectional_property_source == :gxbeam
    @test setup_dict.aero.Vinf == 9.5
end

@testset "Setup aggregate containers" begin
    mesh_props = OWENS.MeshProperties(
        mymesh = :mesh,
        myort = :orientations,
        myjoint = :joint,
        AD15bldNdIdxRng = [1 2; 3 4],
        AD15bldElIdxRng = [5 6; 7 8],
        custom_mesh_outputs = (:components, :s2t),
    )
    @test mesh_props.mymesh === :mesh
    @test mesh_props.myort === :orientations
    @test mesh_props.myjoint === :joint
    @test mesh_props.AD15bldNdIdxRng == [1 2; 3 4]
    @test mesh_props.AD15bldElIdxRng == [5 6; 7 8]
    @test mesh_props.custom_mesh_outputs == (:components, :s2t)

    sectional = OWENS.SectionalProperties([:section], [:stiff], [:mass], [1.0, 0.0])
    @test sectional.sectionPropsArray == [:section]
    @test sectional.stiff_array == [:stiff]
    @test sectional.mass_array == [:mass]
    @test sectional.rotationalEffects == [1.0, 0.0]

    aero = OWENS.AeroProperties(:ad_forces, :ad_deform, :acdms_forces, :acdms_deform)
    @test aero.aeroForcesAD === :ad_forces
    @test aero.deformAeroAD === :ad_deform
    @test aero.aeroForcesACDMS === :acdms_forces
    @test aero.deformAeroACDMS === :acdms_deform
end

@testset "Modeling option and design defaults" begin
    modeling = OWENS.ModelingOptions()
    @test modeling.OWENS_Options.analysisType == "Unsteady"
    @test modeling.OWENS_Options.numTS == 10
    @test modeling.DLC_Options.Vdesign == 11.0
    @test modeling.OWENSAero_Options.Nslices == 20
    @test modeling.OWENSFEA_Options.sectional_property_source == :precomp
    @test modeling.OWENSOpenFASTWrappers_Options.WindType == 3
    @test modeling.Mesh_Options.Nbld == 3
    @test modeling.Drivetrain_Options.gearRatio == 1.0

    defaults = OWENS.Design_Data()
    @test defaults isa OrderedDict{Symbol,Any}
    @test defaults[:name] == "WINDIO Example"
    @test defaults[:assembly][:number_of_blades] == 3
end
