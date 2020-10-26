
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
