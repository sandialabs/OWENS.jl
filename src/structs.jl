
mutable struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
    type
    meshSeg
end

mutable struct Ort
    Psi_d
    Theta_d
    Twist_d
    Length
    elNum
    Offset
end

mutable struct TimeInt
    delta_t
    a1
    a2
    a3
    a4
    a5
    a6
    a7
    a8
end

mutable struct BC_struct
    numpBC
    pBC
    numsBC
    nummBC
    isConstrained
    map
    redVectorMap
end

mutable struct ElInput
    elementOrder
    modalFlag
    timeInt
    xloc
    sectionProps
    sweepAngle
    coneAngle
    rollAngle
    aeroSweepAngle
    iterationType
    useDisp
    preStress
    aeroElasticOn
    aeroForceOn
    loadStepPrev
    loadStep
    maxNumLoadSteps
    MAXIT
    tolerance
    analysisType
    disp
    dispdot
    dispddot
    displ_iter
    concMass
    concStiff
    concLoad
    dispm1
    x
    y
    z
    gravityOn
    RayleighAlpha
    RayleighBeta
    accelVec
    omegaVec
    omegaDotVec
    Omega
    OmegaDot
    CN2H
    airDensity
    freq
    firstIteration
end

mutable struct ElOutput
    FhatLessConc
    Ke
    Fe
    Me
    Ce
end

mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
end

mutable struct CactusGeom
    NBlade
    NStrut
    RotN
    RotP
    RefAR
    RefR
    blade
    strut
end

mutable struct Blade
    NElem
    FlipN
    QCx
    QCy
    QCz
    tx
    ty
    tz
    CtoR
    PEx
    PEy
    PEz
    tEx
    tEy
    tEz
    nEx
    nEy
    nEz
    sEx
    sEy
    sEz
    ECtoR
    EAreaR
    iSect
end

mutable struct Strut
    NElem
    TtoC
    MCx
    MCy
    MCz
    CtoR
    PEx
    PEy
    PEz
    sEx
    sEy
    sEz
    ECtoR
    EAreaR
    BIndS
    EIndS
    BIndE
    EIndE
end

mutable struct SectionPropsArray
    ac
    twist
    rhoA
    EIyy
    EIzz
    GJ
    EA
    rhoIyy
    rhoIzz
    rhoJ
    zcm
    ycm
    a
    EIyz
    alpha1
    alpha2
    alpha3
    alpha4
    alpha5
    alpha6
    rhoIyz
    b
    a0
    aeroCenterOffset
end

mutable struct El
    props
    elLen
    psi
    theta
    roll
    rotationalEffects
end

mutable struct ElStorage
      K11
      K12
      K13
      K14
      K15
      K16
      K22
      K23
      K24
      K25
      K26
      K33
      K34
      K35
      K36
      K44
      K45
      K46
      K55
      K56
      K66
      M11
      M15
      M16
      M22
      M24
      M33
      M34
      M44
      M55
      M56
      M66
      S11
      S12
      S13
      S15
      S16
      S22
      S23
      S25
      S26
      S33
      S35
      S36
      S55
      S56
      S66
      S14_1
      S14_2
      S24_1
      S24_2
      S34_1
      S34_2
      S45_1
      S45_2
      S46_1
      S46_2
      S44_1
      S44_2
      S44_3
      C12
      C13
      C23
      C24
      C25
      C26
      C34
      C35
      C36
      C14_1
      C14_2
      C45_1
      C45_2
      C46_1
      C46_2
      mel
      moiel
      xmel
end

mutable struct NodalTerms
    concLoad
    concStiff
    concMass
    concStiffGen
    concMassGen
    concDampGen
end

mutable struct ConcNDL
    nodeNum
    dof
    val
end

mutable struct ConcNDLGen
    nodeNum
    dof1
    dof2
    val
end

mutable struct ElStrain
    eps_xx_0
    eps_xx_z
    eps_xx_y
    gam_xz_0
    gam_xz_y
    gam_xy_0
    gam_xy_z
end

mutable struct DispOut
    elStrain
    displ_sp1
    displddot_sp1
    displdot_sp1
end

mutable struct DispData
    displ_s
    displdot_s
    displddot_s
end

mutable struct Model
    analysisType
    turbineStartup
    usingRotorSpeedFunction
    tocp
    initCond
    numTS
    delta_t
    Omegaocp
    aeroElasticOn
    aeroForceOn
    aeroLoadsOn
    driveTrainOn
    airDensity
    guessFreq
    gravityOn
    generatorOn
    hydroOn
    JgearBox
    gearRatio
    gearBoxEfficiency
    useGeneratorFunction
    generatorProps
    OmegaGenStart
    omegaControl
    OmegaInit
    totalNumDof
    spinUpOn
    nlOn
    numModesToExtract
    aeroloadfile
    owensfile
    outFilename
    RayleighAlpha
    RayleighBeta
    elementOrder
    joint
    platformTurbineConnectionNodeNumber
    jointTransform
    reducedDOFList
    bladeData
    nlParams
    BC
    nodalTerms
    driveShaftProps
end

# this way you can use defaults and pass in what is different, and it's mapped
# by keyword so it doesn't have to be in order.
function Model(;analysisType = "TNB",
    turbineStartup = 0,
    usingRotorSpeedFunction = false,
    tocp = [0.0,1.1],
    initCond = [],
    numTS = 50.0,
    delta_t = 2e-3,
    Omegaocp = [7.2,7.2] ./ 60,
    aeroElasticOn = false,
    aeroForceOn = true, #this need to get cleaned up in the code
    aeroLoadsOn = false, #this need to get cleaned up in the code
    driveTrainOn = false,
    airDensity=1.2041,
    guessFreq = 0,
    gravityOn = true,
    generatorOn = false,
    hydroOn = false,
    JgearBox = 0.0,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    OmegaGenStart = 0.0,
    omegaControl = false,
    OmegaInit = 7.2/60, #TODO: simplify this in the code since it is redundant
    totalNumDof = 0.0,
    spinUpOn = false,
    nlOn = false,
    numModesToExtract = 20,
    aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv",
    owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens",
    outFilename = "none",
    RayleighAlpha = 0.0,
    RayleighBeta = 0.0,
    elementOrder = 1,
    joint = [0,0],
    platformTurbineConnectionNodeNumber = 1,
    jointTransform = zeros(2,2),
    reducedDOFList = zeros(Int,2),
    numDofPerNode = 6,
    numNodes = 0,
    bladeData = [],
    nlParams = [],
    pBC = 0,
    BC = [],
    nodalTerms = [],
    driveShaftProps = DriveShaftProps(0.0,0.0),
    iterationType = "NR", # nlParams
    adaptiveLoadSteppingFlag = true, # nlParams
    tolerance = 1.0000e-06,# nlParams
    maxIterations = 50,# nlParams
    maxNumLoadSteps = 20,# nlParams
    minLoadStepDelta = 0.0500,# nlParams
    minLoadStep = 0.0500,# nlParams
    prescribedLoadStep = 0.0)# nlParams

    if pBC!=0
        BC = makeBCdata(pBC,numNodes,numDofPerNode)
    end

    nlParams = NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
    maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

    return Model(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,initCond,numTS,delta_t,Omegaocp,
    aeroElasticOn,aeroForceOn,aeroLoadsOn,driveTrainOn,airDensity,
    guessFreq,gravityOn,generatorOn,hydroOn,JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,totalNumDof,
    spinUpOn,nlOn,numModesToExtract,aeroloadfile,owensfile,outFilename,RayleighAlpha,
    RayleighBeta,elementOrder,joint,platformTurbineConnectionNodeNumber,jointTransform,
    reducedDOFList,bladeData,nlParams,BC,nodalTerms,driveShaftProps)
end

mutable struct NlParams
    iterationType
    adaptiveLoadSteppingFlag
    tolerance
    maxIterations
    maxNumLoadSteps
    minLoadStepDelta
    minLoadStep
    prescribedLoadStep
end

mutable struct DriveShaftProps
    k
    c
end
