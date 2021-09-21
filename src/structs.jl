mutable struct Model
    analysisType
    turbineStartup
    usingRotorSpeedFunction
    tocp
    numTS
    delta_t
    Omegaocp
    aeroLoadsOn
    driveTrainOn
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
    numModesToExtract
    aeroloadfile
    owensfile
    outFilename
    bladeData
    driveShaftProps
end

# this way you can use defaults and pass in what is different, and it's mapped
# by keyword so it doesn't have to be in order.
function Model(;analysisType = "TNB",
    turbineStartup = 0,
    usingRotorSpeedFunction = false,
    tocp = [0.0,1.1],
    numTS = 50.0,
    delta_t = 2e-3,
    Omegaocp = [7.2,7.2] ./ 60,
    aeroLoadsOn = false, #this need to get cleaned up in the code
    driveTrainOn = false,
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
    numModesToExtract = 20,
    aeroloadfile = "$module_path/../test/data/input_files_test/DVAWT_2B_LCDT_ElementData.csv",
    owensfile = "$module_path/../test/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.owens",
    outFilename = "none",
    numDofPerNode = 6,
    bladeData = [],
    driveShaftProps = DriveShaftProps(0.0,0.0)
    )

    return Model(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,numTS,delta_t,Omegaocp,
    aeroLoadsOn,driveTrainOn,generatorOn,hydroOn,JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,
    numModesToExtract,aeroloadfile,owensfile,outFilename,bladeData,driveShaftProps)
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

mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
end
