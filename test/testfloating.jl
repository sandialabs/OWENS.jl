using Test
import DelimitedFiles
path = splitdir(@__FILE__)[1]
import OWENS

function runSim(;
    potflowfile = "$path/data/potential_flow_data/marin_semi",
    hd_input_file = "$path/data/HydroDyn_CCT2_test.dat",
    ss_input_file = "$path/data/HydroDyn_CCT2_SeaState_test.dat",
    md_input_file = "$path/data/MoorDyn_CCT2_test.dat",
    topForcingFile = "$(path)/data/PrescribedForcesMoments.csv",
    hd_lib = nothing,
    md_lib = nothing,
    moordyn_on = true,
    bcDOFs = [])

    ##############################################
    # Setup
    #############################################
    dt = .00625 # seconds
    t_max = 1 # seconds
    num_ts = Int(round(t_max/dt) + 1) # +1 since time 0 is considered a time step
    numDOFPerNode = 6

    # create points
    n_ptfm_elem = 1
    n_twr_elem = 10

    # Platform Properties

    bottomMesh = OWENS.OWENSFEA.Mesh(
    [1.0, 2.0], #nodeNum
    1, #numElx
    2, #numNodes
    [0.0, 0.0], #x
    [0.0, 0.0], #y
    [0.0, 10.0], #z
    [1], #elNum
    [1 2], #conn
    [0], #type
    [0, 1], #meshSeg
    Any[], #structuralSpanLocNorm
    Any[], #structuralNodeNumbers
    Any[] #structuralElNumbers
    )

    topMesh = OWENS.OWENSFEA.Mesh(
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0], #nodeNum
    10, #numEl
    11, #numNodes
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #x
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #y
    [0.0, 7.76, 15.52, 23.28, 31.04, 38.8, 46.56, 54.32, 62.08, 69.84, 77.6], #z
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], #elNum
    [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11], #conn
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #type
    [0, 10], #meshSeg
    Any[], #structuralSpanLocNorm
    Any[], #structuralNodeNumbers
    Any[] #structuralElNumbers
    )

    bottom_ort = OWENS.OWENSFEA.Ort(
    [180.0], #Psi_d
    [-90.0],  #Theta_d
    [90.0], #Twist_d
    [10.0], #Length
    [1.0], #elNum
    [0.0; #Offset
     0.0;
     0.0]
    )

    top_ort = OWENS.OWENSFEA.Ort(
    [180.0, 0.0, 180.0, 0.0, 180.0, 0.0, 180.0, 0.0, 180.0, 0.0], #Psi_d
    [-90.0, -90.0, -90.0, -90.0, -90.0, -90.0, -90.0, -90.0, -90.0, -90.0],  #Theta_d
    [90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0], #Twist_d
    [7.76, 7.76, 7.76, 7.76, 7.76, 7.76, 7.76, 7.76, 7.76, 7.76], #Length
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], #elNum
    [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; #Offset
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 7.76 15.52 23.28 31.04 38.8 46.56 54.32 62.08 69.84] #TODO fix this too
    )

    joint = [] #Joint Number, Master Node, Slave Node, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D

    numBottomNodes = length(bottomMesh.z)
    numTopNodes = length(topMesh.z)

    # Set initial displacements
    initDisps = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    initBottomConditions = OWENS.createInitCondArray(initDisps, numBottomNodes, numDOFPerNode) #TODO: generalize this so it's not dependent on numNodes
    initTopConditions = OWENS.createInitCondArray(initDisps, numTopNodes, numDOFPerNode)

    ## Create Sectional Properties

    # Bottom Side
    bottomSectionProps = Array{OWENS.OWENSFEA.SectionPropsArray, 1}(undef, n_ptfm_elem)

    for ii = 1:n_ptfm_elem
        bottomSectionProps[ii] = OWENS.OWENSFEA.SectionPropsArray(
            [0.0, 0.0], #ac
            [0.0, 0.0], #twist_d
            [0.0, 0.0], #rhoA
            [1e18, 1e18], #EIyy
            [1e18, 1e18], #EIzz
            [1e18, 1e18], #GJ
            [1e18, 1e18], #EA
            [0.0, 0.0], #rhoIyy
            [0.0, 0.0], #rhoIzz
            [0.0, 0.0], #rhoJ
            [-8.6588, 0.0], #zcm
            [0.0, 0.0], #ycm,
            [0.0, 0.0], #a
            [0.0, 0.0], #EIyz
            [0.0, 0.0], #alpha1
            [0.0, 0.0], #alpha2
            [0.0, 0.0], #alpha3
            [0.0, 0.0], #alpha4
            [0.0, 0.0], #alpha5
            [0.0, 0.0], #alpha6
            [0.0, 0.0], #rhoIyz
            [0.0, 0.0], #b
            [0.0, 0.0], #a0
            [0.0, 0.0], #aeroCenterOffset
        )
    end

    ptfm_mass = 3.85218e6
    ptfm_roll_iner = 2.56193e9
    ptfm_pitch_iner = 2.56193e9
    ptfm_yaw_iner = 4.24265e9
    ptfm_Ixx = ptfm_roll_iner
    ptfm_Iyy = ptfm_pitch_iner
    ptfm_Izz = ptfm_yaw_iner
    ptfm_Ixy = 0.0
    ptfm_Iyz = 0.0
    ptfm_Ixz = 0.0
    ptfm_mass = [
        ptfm_mass                0.0                      0.0                      0.0             0.0        0.0
        0.0                       ptfm_mass                0.0                      0.0             0.0        0.0
        0.0                       0.0                      ptfm_mass                0.0             0.0        0.0
        0.0                       0.0                      0.0                      ptfm_Ixx        ptfm_Ixy   ptfm_Ixz
        0.0                       0.0                      0.0                      ptfm_Ixy        ptfm_Iyy  -ptfm_Iyz
        0.0                       0.0                      0.0                      ptfm_Ixz       -ptfm_Iyz   ptfm_Izz
        ]

    # Top Side
    topSectionProps = Array{OWENS.OWENSFEA.SectionPropsArray, 1}(undef, n_twr_elem)

    twr_rhoA = [4.667e3, 4.34528e3, 4.03476e3, 3.73544e3, 3.44732e3, 3.1704e3, 2.90469e3, 2.65018e3, 2.40688e3, 2.17477e3, 1.95387e3]
    twr_EA = [115.302e9, 107.354e9, 99.682e9, 92.287e9,  85.169e9, 78.328e9, 71.763e9, 65.475e9, 59.464e9, 53.730e9, 48.272e9]
    twr_GJ = [464.718e9, 398.339e9, 339.303e9, 287.049e9, 241.043e9, 200.767e9, 165.729e9, 135.458e9, 109.504e9, 87.441e9, 68.863e9]
    twr_EIyy = [603.903e9, 517.644e9, 440.925e9, 373.022e9, 313.236e9, 260.897e9, 215.365e9, 176.028e9, 142.301e9, 113.630e9, 89.488e9]
    twr_rhoIyy = [24443.7, 20952.2, 17847.0, 15098.5, 12678.6, 10560.1, 8717.2, 7124.9, 5759.8, 4599.3, 3622.1] #kg-m
    twr_rhoJ = copy(twr_rhoIyy)
    twr_EIzz = copy(twr_EIyy)
    twr_rhoIzz = copy(twr_rhoIyy)
    twr_zcm = [0.0, 0.0]
    twr_ycm = [0.0, 0.0]

    for ii = 1:n_twr_elem
        rhoA = [twr_rhoA[ii], twr_rhoA[ii+1]]
        EIyy = [twr_EIyy[ii], twr_EIyy[ii+1]]
        EIzz = [twr_EIzz[ii], twr_EIzz[ii+1]]
        GJ = [twr_GJ[ii], twr_GJ[ii+1]]
        EA = [twr_EA[ii], twr_EA[ii+1]]
        rhoIyy = [twr_rhoIyy[ii], twr_rhoIyy[ii+1]]
        rhoIzz = [twr_rhoIzz[ii], twr_rhoIzz[ii+1]]
        rhoJ = [twr_rhoJ[ii], twr_rhoJ[ii+1]]

        topSectionProps[ii] = OWENS.OWENSFEA.SectionPropsArray(
            [0.0, 0.0], #ac
            [0.0, 0.0], #twist_d
            rhoA,
            EIyy,
            EIzz,
            GJ,
            EA,
            rhoIyy,
            rhoIzz,
            rhoJ,
            [0.0, 0.0], #zcm
            [0.0, 0.0], #ycm,
            [0.0, 0.0], #a
            [0.0, 0.0], #EIyz
            [0.0, 0.0], #alpha1
            [0.0, 0.0], #alpha2
            [0.0, 0.0], #alpha3
            [0.0, 0.0], #alpha4
            [0.0, 0.0], #alpha5
            [0.0, 0.0], #alpha6
            [0.0, 0.0], #rhoIyz
            [0.0, 0.0], #b
            [0.0, 0.0], #a0
            [0.0, 0.0], #aeroCenterOffset
        )
    end

    bottomRotationalEffects = ones(bottomMesh.numEl)
    topRotationalEffects = ones(topMesh.numEl)

    #store data in element object
    bottomEl = OWENS.OWENSFEA.El(bottomSectionProps,bottom_ort.Length,bottom_ort.Psi_d,bottom_ort.Theta_d,bottom_ort.Twist_d,bottomRotationalEffects)
    topEl = OWENS.OWENSFEA.El(topSectionProps,top_ort.Length,top_ort.Psi_d,top_ort.Theta_d,top_ort.Twist_d,topRotationalEffects)

    # Add the platform properties
    bottomConcInputs = [
    1 "M6" 1 1 ptfm_mass[1,1]
    1 "M6" 1 2 ptfm_mass[1,2]
    1 "M6" 1 3 ptfm_mass[1,3]
    1 "M6" 1 4 ptfm_mass[1,4]
    1 "M6" 1 5 ptfm_mass[1,5]
    1 "M6" 1 6 ptfm_mass[1,6]
    1 "M6" 2 1 ptfm_mass[2,1]
    1 "M6" 2 2 ptfm_mass[2,2]
    1 "M6" 2 3 ptfm_mass[2,3]
    1 "M6" 2 4 ptfm_mass[2,4]
    1 "M6" 2 5 ptfm_mass[2,5]
    1 "M6" 2 6 ptfm_mass[2,6]
    1 "M6" 3 1 ptfm_mass[3,1]
    1 "M6" 3 2 ptfm_mass[3,2]
    1 "M6" 3 3 ptfm_mass[3,3]
    1 "M6" 3 4 ptfm_mass[3,4]
    1 "M6" 3 5 ptfm_mass[3,5]
    1 "M6" 3 6 ptfm_mass[3,6]
    1 "M6" 4 1 ptfm_mass[4,1]
    1 "M6" 4 2 ptfm_mass[4,2]
    1 "M6" 4 3 ptfm_mass[4,3]
    1 "M6" 4 4 ptfm_mass[4,4]
    1 "M6" 4 5 ptfm_mass[4,5]
    1 "M6" 4 6 ptfm_mass[4,6]
    1 "M6" 5 1 ptfm_mass[5,1]
    1 "M6" 5 2 ptfm_mass[5,2]
    1 "M6" 5 3 ptfm_mass[5,3]
    1 "M6" 5 4 ptfm_mass[5,4]
    1 "M6" 5 5 ptfm_mass[5,5]
    1 "M6" 5 6 ptfm_mass[5,6]
    1 "M6" 6 1 ptfm_mass[6,1]
    1 "M6" 6 2 ptfm_mass[6,2]
    1 "M6" 6 3 ptfm_mass[6,3]
    1 "M6" 6 4 ptfm_mass[6,4]
    1 "M6" 6 5 ptfm_mass[6,5]
    1 "M6" 6 6 ptfm_mass[6,6]
    ]

    topConcInputs = [
    length(topMesh.z) "M6" 1 1 350000
    length(topMesh.z) "M6" 2 2 350000
    length(topMesh.z) "M6" 3 3 350000
    ]

    bottomConcTerms = OWENS.OWENSFEA.applyConcentratedTerms(bottomMesh.numNodes, numDOFPerNode, data=bottomConcInputs, jointData=joint)
    topConcTerms = OWENS.OWENSFEA.applyConcentratedTerms(topMesh.numNodes, numDOFPerNode, data=topConcInputs, jointData=joint)

    # node, dof, bc
    fixedTopNodes = [1]
    fixedBottomNodes = []
    topBCs = OWENS.setBCs(bcDOFs, fixedTopNodes, numTopNodes, numDOFPerNode)
    bottomBCs = OWENS.setBCs(bcDOFs, fixedBottomNodes, numBottomNodes, numDOFPerNode)

    # Setting library paths
    if !moordyn_on
        md_lib = nothing
    end
    bin = OWENS.Bin(hd_lib, md_lib)

    inputs = OWENS.Inputs(analysisType = "TNB",
    outFilename = "none",
    tocp = [0.0,num_ts],
    Omegaocp = [0.0, 0.0],
    OmegaInit = 0.0,
    hydroOn = true,
    aeroLoadsOn = 2,
    interpOrder = 2,
    numTS = num_ts,
    delta_t = dt,
    hd_input_file = hd_input_file,
    ss_input_file = ss_input_file,
    md_input_file = md_input_file,
    potflowfile = potflowfile)

    if topBCs == []
        topModel = OWENS.OWENSFEA.FEAModel(analysisType = "TNB",
        joint = joint,
        nlOn = true,
        platformTurbineConnectionNodeNumber = 1,
        initCond = initTopConditions,
        iterationType = "DI",
        gamma = 1.0,
        alpha = 0.5,
        gravityOn = [0.0,0.0,-9.80665],
        numNodes = topMesh.numNodes,
        nodalTerms = topConcTerms,
        RayleighAlpha = 0.01/pi)
    else
        topModel = OWENS.OWENSFEA.FEAModel(analysisType = "TNB",
        joint = joint,
        nlOn = true,
        platformTurbineConnectionNodeNumber = 1,
        initCond = initTopConditions,
        iterationType = "DI",
        numNodes = topMesh.numNodes,
        nodalTerms = topConcTerms,
        RayleighAlpha = 0.01/pi,
        gamma = 1.0,
        alpha = 0.5,
        gravityOn = [0.0,0.0,-9.80665],
        pBC = topBCs)
    end

    if bottomBCs == []
        bottomModel = OWENS.OWENSFEA.FEAModel(analysisType = "TNB",
        joint = joint,
        nlOn = true,
        platformTurbineConnectionNodeNumber = 1,
        initCond = initBottomConditions,
        iterationType = "DI",
        gamma = 1.0,
        alpha = 0.5,
        gravityOn = [0.0,0.0,-9.80665],
        numNodes = bottomMesh.numNodes,
        nodalTerms = bottomConcTerms
        )
    else
        bottomModel = OWENS.OWENSFEA.FEAModel(analysisType = "TNB",
        joint = joint,
        nlOn = true,
        platformTurbineConnectionNodeNumber = 1,
        initCond = initBottomConditions,
        iterationType = "DI",
        gamma = 1.0,
        alpha = 0.5,
        gravityOn = [0.0,0.0,-9.80665],
        numNodes = bottomMesh.numNodes,
        nodalTerms = bottomConcTerms,
        pBC = bottomBCs
        )
    end

    ##############################################
    # Unsteady Test
    #############################################

    topForcing = DelimitedFiles.readdlm(topForcingFile, ',', Float32, '\n', header=false)
    topFrcValArray = topForcing[:, 2:end]
    topFrcDOFs = (numTopNodes-1)*numDOFPerNode+1:numTopNodes*numDOFPerNode

    deformAero(omega) = 0.0 #placeholder function
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,FTwrBsHist,
    genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,
    epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FPtfmHist,FHydroHist,FMooringHist = OWENS.Unsteady(inputs,
                                                topModel=topModel,
                                                topMesh=topMesh,
                                                topEl=topEl,
                                                aeroVals=topFrcValArray,
                                                aeroDOFs=topFrcDOFs,
                                                deformAero=deformAero,
                                                bottomModel=bottomModel,
                                                bottomMesh=bottomMesh,
                                                bottomEl=bottomEl,
                                                bin=bin)

    return uHist_prp, FPtfmHist, FHydroHist, FMooringHist, uHist[:,topFrcDOFs], FReactionHist
end

ptfm_disps, ptfm_forces, hydro_forces, mooring_forces, tt_disps, FReaction = runSim(
    potflowfile = "$path/data/potential_flow_data/marin_semi",
    hd_input_file = "$path/data/HydroDyn_CCT2_test.dat",
    ss_input_file = "$path/data/HydroDyn_CCT2_SeaState_test.dat",
    md_input_file = "$path/data/MoorDyn_CCT2_test.dat",
    topForcingFile = "$(path)/data/PrescribedForcesMoments.csv",
    hd_lib = nothing,
    md_lib = nothing)


ptfm_disps_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_xPtfm.csv", ',')
ptfm_forces_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_FPtfm.csv", ',')
hydro_forces_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_FHydro.csv", ',')
mooring_forces_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_FMooring.csv", ',')
tt_disps_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_xTowerTop.csv", ',')
FReaction_UNIT = DelimitedFiles.readdlm("$path/data/UNIT_TEST_cct2_FReaction.csv", ',')

mytol = 0.00001

for iel = 1:length(ptfm_disps_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(ptfm_disps_UNIT[iel],ptfm_disps[iel];atol=abs(ptfm_disps_UNIT[iel])*mytol)
end

for iel = 1:length(ptfm_forces_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(ptfm_forces_UNIT[iel],ptfm_forces[iel];atol=abs(ptfm_forces_UNIT[iel])*mytol)
end

for iel = 1:length(hydro_forces_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(hydro_forces_UNIT[iel],hydro_forces[iel];atol=abs(hydro_forces_UNIT[iel])*mytol)
end

for iel = 1:length(mooring_forces_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(mooring_forces_UNIT[iel],mooring_forces[iel];atol=abs(mooring_forces_UNIT[iel])*mytol)
end

for iel = 1:length(tt_disps_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(tt_disps_UNIT[iel],tt_disps[iel];atol=abs(tt_disps_UNIT[iel])*mytol)
end

for iel = 1:length(FReaction_UNIT) #note, purposely iterating over multidimensional array with a single for loop
    @test isapprox(FReaction_UNIT[iel],FReaction[iel];atol=abs(FReaction_UNIT[iel])*mytol)
end


# FReactionHist = FReaction
# dt = .00625 # seconds
# t_max = 1 # seconds
# FReaction_tvec = collect(dt:dt:t_max)

# Fz_orig = DelimitedFiles.readdlm("$path/data/Fz_orig.csv", ',')
# Fz_oldowens = DelimitedFiles.readdlm("$path/data/Fz_MD_newsim.csv", ',')

# Fy_oldowens = DelimitedFiles.readdlm("$path/data/Fy_MD_newsim.csv", ',')

# Mz_oldowens = DelimitedFiles.readdlm("$path/data/Mz_MD_newsim.csv", ',')

# Fx = FReactionHist[:,1]
# Fy = FReactionHist[:,2]
# Fz = FReactionHist[:,3]
# Mx = FReactionHist[:,4]
# My = FReactionHist[:,5]
# Mz = FReactionHist[:,6]


# Fx_UNIT = FReaction_UNIT[:,1]
# Fy_UNIT = FReaction_UNIT[:,2]
# Fz_UNIT = FReaction_UNIT[:,3]
# Mx_UNIT = FReaction_UNIT[:,4]
# My_UNIT = FReaction_UNIT[:,5]
# Mz_UNIT = FReaction_UNIT[:,6]

# using PyPlot
# PyPlot.pygui(true)

# PyPlot.figure("Fx")
# PyPlot.plot(FReaction_tvec,Fx)
# PyPlot.plot(FReaction_tvec,Fx_UNIT)

# PyPlot.figure("Fy")
# # PyPlot.plot(Fy_oldowens[:,1],Fy_oldowens[:,2],label="Old OWENS")
# PyPlot.plot(FReaction_tvec,Fy,label="Latest OWENS")
# PyPlot.plot(FReaction_tvec,Fy_UNIT,label="MD OWENS")
# # PyPlot.xlim([300,600])
# # PyPlot.ylim([-150000,100000])
# PyPlot.legend()

# PyPlot.figure("Fz")
# # PyPlot.plot(Fz_orig[:,1],Fz_orig[:,2].*1e6,"k",label="OpenFAST")
# # PyPlot.plot(Fz_oldowens[:,1],Fz_oldowens[:,2].*1e6,label="Old OWENS")
# PyPlot.plot(FReaction_tvec,Fz,label="Latest OWENS")
# PyPlot.plot(FReaction_tvec,Fz_UNIT,label="MD OWENS")
# # PyPlot.xlim([300,600])
# # PyPlot.ylim([-5.5e6,-6.0e6])
# PyPlot.legend()

# PyPlot.figure("Mx")
# PyPlot.plot(FReaction_tvec,Mx)

# PyPlot.figure("My")
# PyPlot.plot(FReaction_tvec,My)

# PyPlot.figure("Mz")
# # PyPlot.plot(Mz_oldowens[:,1],Mz_oldowens[:,2],label="Old OWENS")
# PyPlot.plot(FReaction_tvec,Mz_UNIT,label="MD OWENS")
# PyPlot.plot(FReaction_tvec,Mz,label="Latest OWENS")
# # PyPlot.xlim([300,600])
# # PyPlot.ylim([-5.5e6,-6.0e6])
# PyPlot.legend()
