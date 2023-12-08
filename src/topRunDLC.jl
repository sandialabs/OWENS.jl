
mutable struct MasterInput
    analysisType
    turbineType
    eta
    Nbld
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
    turbsim_filename
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
    ifw_libfile = "./../openfast/build/modules/inflowwind/libifw_c_binding",
    numTS = 100,
    delta_t = 0.01,
    turbsim_filename ="$module_path/../test/data/turbsim/115mx115m_30x30_20.0msETM.bts",
    NuMad_geom_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv",
    NuMad_mat_xlscsv_file_twr = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv",
    NuMad_mat_xlscsv_file_bld = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv",
    NuMad_geom_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Geom_SNL_5MW_Struts.csv",
    NuMad_mat_xlscsv_file_strut = "$module_path/../test/examples/data/NuMAD_Materials_SNL_5MW.csv"
    )

    return MasterInput(analysisType,turbineType,eta,Nbld,towerHeight,rho,Vinf,controlStrategy,
    RPM,Nslices,ntheta,structuralModel,ntelem,nbelem,ncelem,nselem,AModel,ifw,turbsim_filename,ifw_libfile,adi_lib,adi_rootname,
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
        turbsim_filename = turbulentInflow["turbsim_filename"]
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
    nselem,AModel,ifw,turbsim_filename,ifw_libfile,adi_lib,adi_rootname,Blade_Height,Blade_Radius,numTS,
    delta_t,NuMad_geom_xlscsv_file_twr,NuMad_mat_xlscsv_file_twr,
    NuMad_geom_xlscsv_file_bld,NuMad_mat_xlscsv_file_bld,NuMad_geom_xlscsv_file_strut,NuMad_mat_xlscsv_file_strut)
end

function runOWENS(Inp,path;verbosity=2)
    # Unpack inputs
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
    AModel = Inp.AModel
    ifw = Inp.ifw
    turbsim_filename = Inp.turbsim_filename
    ifw_libfile = Inp.ifw_libfile
    Blade_Height = Inp.Blade_Height
    Blade_Radius = Inp.Blade_Radius
    numTS = Inp.numTS
    delta_t = Inp.delta_t
    NuMad_geom_xlscsv_file_twr = Inp.NuMad_geom_xlscsv_file_twr
    NuMad_mat_xlscsv_file_twr = Inp.NuMad_mat_xlscsv_file_twr
    NuMad_geom_xlscsv_file_bld = Inp.NuMad_geom_xlscsv_file_bld
    NuMad_mat_xlscsv_file_bld = Inp.NuMad_mat_xlscsv_file_bld
    NuMad_geom_xlscsv_file_strut = Inp.NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut = Inp.NuMad_mat_xlscsv_file_strut
    adi_lib = Inp.adi_lib
    adi_rootname = Inp.adi_rootname

    println("Set up Turbine")

    B = Nbld
    R = Blade_Radius#177.2022*0.3048 #m
    H = Blade_Height#1.02*R*2 #m

    shapeY = collect(LinRange(0,H,Nslices+1))
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2)#shapeX_spline(shapeY)

    mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
    mass_breakout_blds,mass_breakout_twr = OWENS.setupOWENS(VAWTAero,path;
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
        delta_t,
        numTS,
        adi_lib,
        adi_rootname,
        turbsim_filename,
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

    println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
    system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld)#;VTKmeshfilename="ARCUS5MW")

    if verbosity>0
        # Print Blades and Tower Materials Breakdown and Costs
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
        
        # println("\nBlades' Material Max Strain")
        # for (i,name) in enumerate(plyprops_bld.names)
        #     println("$name $(plyprops_bld.plies[i].xt/plyprops_bld.plies[i].e1) xt $(plyprops_bld.plies[i].xc/plyprops_bld.plies[i].e1) xc $(plyprops_bld.plies[i].yt/plyprops_bld.plies[i].e2) yt $(plyprops_bld.plies[i].yc/plyprops_bld.plies[i].e2) yc")
        # end
    end

    ######################################
    #### Perform Aerostructural One Way Test
    #######################################

    pBC = [1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0]


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

    feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    outFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = true,
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI",
    predef = "update")

    println("Running Unsteady")
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)
    
    println("Saving VTK time domain files")
    OWENS.gyricFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

    ##########################################
    #### Ultimate Failure #####
    ##########################################

    massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
    SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
    topstrainout_tower_U,topstrainout_tower_L = OWENS.extractSF(bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
    kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor
    # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
    LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
    LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1) #TODO: add in ability to have material safety factors and load safety factors

    ##########################################
    #### Fatigue #####
    ##########################################

    # DEL

    ##########################################
    #### Data Dump in OpenFAST Format #####
    ##########################################

end

# # Test
# Inp = OWENS.MasterInput(;numTS=3,ifw_libfile="$localpath/../../openfastandy/build/modules/inflowwind/libifw_c_binding")
# OWENS.runDLC(["1_1","1_2"],Inp,localpath;Vinf_range=LinRange(5,20,2),regenTurbsim=true,pathtoturbsim="$localpath/../../openfastandy/build/modules/turbsim/turbsim")
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
    regenTurbsim=false)

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
    * `Vdesign`: =11.0,
    * `grid_oversize`: =1.1,
    * `regenTurbsim`: =false

    # Output
    * `nothing`: 
    """
function runDLC(DLCs,Inp,path;
    Vinf_range=LinRange(5,20,16),
    IEC_std="\"2\"",
    WindChar="\"A\"",
    WindClass=1,
    turbsimpath="./turbsimfiles",
    templatefile="$module_path/template_files/templateTurbSim.inp",
    pathtoturbsim="../../openfast/build/modules/turbsim/turbsim",
    NumGrid_Z=100,
    NumGrid_Y=100,
    Vref=10.0,
    Vdesign=11.0,
    grid_oversize=1.1,
    regenTurbsim=false,
    runScript = OWENS.runOWENS)

    if !isdir(turbsimpath)
        mkdir(turbsimpath)
    end

    # Fill in DLC parameters based on model inputs
    DLCParams = Array{DLCParameters, 1}(undef, length(DLCs))

    for (iDLC, DLC) in enumerate(DLCs) #TODO parallelize this
        if IEC_std=="\"2\""
            DLCParams[iDLC] = getDLCparams(DLC, Inp, Vinf_range, Vdesign, Vref, WindChar,WindClass, IEC_std;grid_oversize)
        else
            error("Only IEC61400_2 is currently implemented")
        end

        # Run Simulation at each Wind Speed
        for windspeed in DLCParams[iDLC].Vinf_range_used #TODO: parallelize this

            

            DLCParams[iDLC].URef = windspeed
            # Check if turbulent inflow file exists, if not create it
            windspeedStr = round(windspeed;digits=2)
            windspeedStr = lpad(windspeedStr,4,"0")
            println("Running DLC $DLC at Vinf $windspeedStr m/s")
            turbsimINPfilename = "$turbsimpath/DLC$(DLC)Vinf$(windspeedStr).inp"
            if !isfile(turbsimINPfilename) || regenTurbsim
                generateTurbsimBTS(DLCParams[iDLC],turbsimINPfilename,pathtoturbsim;templatefile)
            end

            Inp.ifw = true #TODO: steady non-ifw?
            Inp.turbsim_filename = "$(turbsimINPfilename[1:end-4]).bts"

            # run owens simulation
            runScript(Inp,path)
        end
    end
end

mutable struct DLCParameters
    Vinf_range_used
    analysis_type # "U", "F", "UF"
    TorqController # "constRPM", function handle
    RandSeed1 # Turbulent Random Seed Number
    NumGrid_Z # Vertical grid-point matrix dimension
    NumGrid_Y # Horizontal grid-point matrix dimension
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
end


function getDLCparams(DLC, Inp, Vinf_range, Vdesign, Vref, WindChar, WindClass, IEC_std;grid_oversize=1.2)

    Ve50 = 50.0 #TODO
    Ve1 = 30.0 #TODO

    numTS = Inp.numTS
    delta_t = Inp.delta_t
    simtime = numTS*delta_t

    Blade_Radius = Inp.Blade_Radius
    Blade_Height = Inp.Blade_Height

    NumGrid_Z = Inp.ntelem+Inp.nbelem
    NumGrid_Y = Inp.ntelem+Inp.nbelem

    RandSeed1 = 40071 #TODO
    HubHt = (Inp.towerHeight+Inp.Blade_Height)*grid_oversize/2 #TODO
    AnalysisTime = simtime
    GridHeight = (Inp.towerHeight+Inp.Blade_Height)*grid_oversize
    GridWidth = ceil((Blade_Radius) * 2.0 * grid_oversize)
    VFlowAng = 0.0
    HFlowAng = 0.0
    
    IECstandard = IEC_std
    IECturbc = WindChar
    TurbModel = "\"IECKAI\""
    
    RefHt = round(Blade_Height) #TODO
    URef = 0.0 #gets filled in later from the Vinf_range when the .bst is generated

    TimeStep = delta_t

    if DLC == "1_1"
        TorqController = "normal"
        Vinf_range_used = Vinf_range
        analysis_type = "UF"
        IEC_WindType = "\"$(WindClass)NTM\""
        

    elseif DLC == "1_2"
        TorqController = "normal"
        Vinf_range_used = [Vdesign]
        analysis_type = "F"
        IEC_WindType = "\"$(WindClass)ECD\""
        

    elseif DLC == "1_3"
        TorqController = "normal"
        Vinf_range_used = Vinf_range
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EOG50\""
        

    elseif DLC == "1_4"
        TorqController = "normal"
        Vinf_range_used = Vinf_range
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)ECD50\""
        

    elseif DLC == "1_5"
        TorqController = "normal"
        Vinf_range_used = [Vdesign]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)ECG\""
        

    elseif DLC == "2_1"
        TorqController = "freewheelatNormalOperatingRPM"
        Vinf_range_used = [Vdesign]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)NWP\""
        

    elseif DLC == "2_2"
        TorqController = "freewheelatNormalOperatingRPM"
        Vinf_range_used = Vinf_range
        analysis_type = "UF"
        IEC_WindType = "\"$(WindClass)NTM\""
        

    elseif DLC == "2_3"
        TorqController = "freewheelatNormalOperatingRPM"
        Vinf_range_used = Vinf_range
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EOG1\""
        

    elseif DLC == "3_1"
        TorqController = "shutdown"
        Vinf_range_used = Vinf_range
        analysis_type = "F"
        IEC_WindType = "\"$(WindClass)NTM\""
        

    elseif DLC == "3_2"
        TorqController = "shutdown"
        Vinf_range_used = [Vinf_range[end]]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EOG1\""
        

    elseif DLC == "4_1"
        TorqController = "shutdown"
        Vinf_range_used = [Vdesign]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)NTM\""
        

    elseif DLC == "5_1"
        TorqController = "freewheelatIdle"
        Vinf_range_used = [Ve50]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EWM\""
        

    elseif DLC == "5_2"
        TorqController = "idle"
        Vinf_range_used = [Vdesign]
        analysis_type = "F"
        IEC_WindType = "\"$(WindClass)NTM\""
        

    elseif DLC == "6_1"
        TorqController = "parked"
        Vinf_range_used = [Ve1]
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EWM\""
        
    elseif DLC == "8_1" #Startup
        TorqController = "startup"
        Vinf_range_used = Vinf_range
        analysis_type = "U"
        IEC_WindType = "\"$(WindClass)EWM\""
        
    else
        error("IEC61400_2 DLCs [1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,3.1,3.2,4.1,5.1,5.2,6.1] defined, you requested $DLC")
    end

    return DLCParameters(
        Vinf_range_used,
        analysis_type, # array of windspeeds m/s
        TorqController, # "constRPM", function handle
        RandSeed1, # Turbulent Random Seed Number
        NumGrid_Z, # Vertical grid-point matrix dimension
        NumGrid_Y, # Horizontal grid-point matrix dimension
        TimeStep, # Time step [s]
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
    )
end

function generateTurbsimBTS(DLCParams,turbsimINPfilename,pathtoturbsim;templatefile="$localpath/templateTurbSim.inp") 

    lines = readlines(templatefile)

    for fieldname in fieldnames(typeof(DLCParams))
        turbsimKeyName = String(fieldname)
        myvalue = getfield(DLCParams,fieldname)
        if turbsimKeyName != "Vinf_range_used" || turbsimKeyName != "analysis_type" || turbsimKeyName != "TorqController"# || other strings
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
    open(turbsimINPfilename, "w") do file
        # Write new data to file
        for line in lines
            write(file, "$(line)\n")
        end
    end

    run(`$pathtoturbsim $turbsimINPfilename`)
end

# # Test
# Inp = OWENS.MasterInput(;numTS=3,ifw_libfile="$localpath/../../openfastandy/build/modules/inflowwind/libifw_c_binding")
# OWENS.runDLC(["1_1","1_2"],Inp,localpath;Vinf_range=LinRange(5,20,2),regenTurbsim=true,pathtoturbsim="$localpath/../../openfastandy/build/modules/turbsim/turbsim")
