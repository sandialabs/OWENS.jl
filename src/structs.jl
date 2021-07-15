
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
    aeroLoadsOn
    driveTrainOn
    airDensity
    guessFreq
    gravityOn
    generatorOn
    hydroOn
    plat_model
    JgearBox
    gearRatio
    gearBoxEfficiency
    useGeneratorFunction
    generatorProps
    OmegaGenStart
    omegaControl
    OmegaInit
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
    aeroLoadsOn = false, #this need to get cleaned up in the code
    driveTrainOn = false,
    airDensity=1.2041,
    guessFreq = 0,
    gravityOn = true,
    generatorOn = false,
    hydroOn = false,
    plat_model = [],
    JgearBox = 0.0,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    OmegaGenStart = 0.0,
    omegaControl = false,
    OmegaInit = 7.2/60, #TODO: simplify this in the code since it is redundant
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
    jointTransform = 0.0,
    reducedDOFList = zeros(Int,2),
    numDofPerNode = 6,
    numNodes = 0,
    bladeData = [],
    nlParams = 0,
    pBC = 0,
    BC = [],
    nodalTerms = 0.0,
    driveShaftProps = DriveShaftProps(0.0,0.0),
    iterationType = "NR",
    adaptiveLoadSteppingFlag = true,
    tolerance = 1.0000e-06,
    maxIterations = 50,
    maxNumLoadSteps = 20,
    minLoadStepDelta = 0.0500,
    minLoadStep = 0.0500,
    prescribedLoadStep = 0.0)

    if jointTransform==0.0
        jointTransform, reducedDOFList = GyricFEA.createJointTransform(joint,numNodes,numDofPerNode) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    end
    if nodalTerms == 0.0
        nodalTerms = OWENS.readNodalTerms() #Fill in the data structure with nothing
    end

    if pBC!=0
        BC = makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)
    else
        BC = GyricFEA.BC_struct(0,
        0,
        0,
        0,
        0,
        0,
        0)
    end

    if nlParams==0
        nlParams = NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
        maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)
    end

    return Model(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,initCond,numTS,delta_t,Omegaocp,
    aeroElasticOn,aeroLoadsOn,driveTrainOn,airDensity,
    guessFreq,gravityOn,generatorOn,hydroOn,plat_model,JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,
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

# Cactus Related Structs

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
