
struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
    type
end

mutable struct Ort
    Psi_d
    Theta_d
    Twist_d
    Length
    elNum
    Offset
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

mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
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
