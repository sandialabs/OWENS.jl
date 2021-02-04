
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
    structuralSpanLocNorm
    structuralNodeNumbers
    structuralElNumbers
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
    jointTransform = 0.0,
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

    if jointTransform==0.0
        jointTransform, reducedDOFList = createJointTransform(joint,numNodes,numDofPerNode) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    end

    if pBC!=0
        BC = makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)
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


"""
NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin,web_stack,web_dp)

Parameters defining the blade composite layup. See NuMad user guide SAND2012_7028 appendix B for more details

**Arguments**
- `n_web::Int64`: number of shear webs
- `n_stack::Int64`: number of predefined composite stacks
- `n_segments::Int64`: number of segments around the airfoil
- `span::Vector{Float64}`: span-wise position
- `airfoil::Vector{String}`: airfoil name
- `te_type::Vector{String}`: trailing edge type
- `twist_d::Vector{Float64}`: twist_d in degrees
- `chord::Vector{Float64}`: chord length
- `xoffset::Vector{Float64}`: The distance from the “nose” of a station to the blade reference axis.
- `aerocenter::Vector{Float64}`: This is an aerodynamic parameter that is an output from aerodynamic performance analysis of a two-dimensional airfoil section. The aerodynamic center is the point along the chord where the aerodynamic pitching moment does not vary with changes in angle of attack.
- `stack_mat_types::Vector{Int64}`: Material numbers used that correspond to each stack number
- `stack_layers::Array{Int64,2}`: number of layers at each span used corresponding to each material type (first index corresponds to spanwise position, second index corresponds to the stack number)
- `segments::Array{Float64,2}`: normalized starting and stopping points of each section (i.e. leading edge, sparcap, etc).  First index corresponds to  spanwise position, second index corresponds to the section, except there is an extra first column starting at -1 for the trailing edge. There must be a leading edge position at 0, and the last column must be 1 corresponding to the trailing edge again.  Positions are fractions of the chord, lower (HP) is negative, upper (LP) is positive
- `DPtypes::Array{Int64,2}`: division point types (NOTE THAT THIS ISN'T IMPLEMENTED AND DOES NOTHING CURRENTLY, i.e. only SINGLE is being used). First index corresponds to  spanwise positoin, second corresponds to section number
- `skin_seq::Array{Seq,2}`: stack sequence, is an array of structures, each containing a Vector{Int64} of the sequence (i.e. skin[2,5].seq). First index corresponds to spanwise positoin, second index the section
- `web_seq::Array{Seq,2}`: same format and meaning as skin sequence, but for the webs with the second index corresponding to the web number
- `web_dp::Array{Seq,2}`: same format as skin sequence, but this corresponds to the section numbers the web connects to at the top and bottom at both edges. There are always four entries in the CSV list and the order goes as follows: inboard LP, inboard HP, outboard HP, outboard LP.
"""

mutable struct NuMad
    n_web
    n_stack
    n_segments
    span
    airfoil
    te_type
    twist_d
    chord
    xoffset
    aerocenter
    stack_mat_types
    stack_layers
    segments
    DPtypes
    skin_seq
    web_seq
    web_dp
end

mutable struct Seq
    seq
end

struct plyproperties
    names#::Array{String,1}
    plies#::Array{Composites.Material,1}
end



##------- MATERIAL PROPERTIES -------##
function plyproperties()
    #TODO: read in numad xls materials file
    names = ["highmodulus_uni",
    "standard_uni",
    "MR60H",
    "T3900_uni",
    "T700_uni",
    "ELT5500",
    "UDCarbon",
    "highmodulus_weave",
    "standard_weave",
    "T3900_weave",
    "T700_weave",
    "Gelcoat",
    "Triax",
    "Saertex",
    "taylor_foam",
    "SNL_foam",
    "aluminum6061",
    "aluminum6063",
    "CLA_5500",
    "CBX_2400",
    "ETLX_2400",
    "Airex_C70_55",
    "EBX_2400_x10",
    "ETLX_2400_x10",
    "Airex_C70_55_x10"]

    n_materials = 18
    e1  =zeros(n_materials) #pa
    e2  =zeros(n_materials) #pa
    g12 =zeros(n_materials) #pa
    anu =zeros(n_materials) #ratio
    rho =zeros(n_materials) #kg/m3
    xt  =zeros(n_materials) #pa
    xc  =zeros(n_materials) #pa
    yt  =zeros(n_materials) #pa
    yc  =zeros(n_materials) #pa
    s   =zeros(n_materials) #pa
    plythickness =zeros(n_materials) #meters

    # "highmodulus_uni"
    e1[1]  = 175.0e9
    e2[1]  = 8.0e9
    g12[1] = 5.0e9
    anu[1] = 0.30
    rho[1] = 1600.0
    xt[1]  = min(1.0,1355.656/1682.011)*1000e6 #mean -> A-basis
    xc[1]  = min(1.0,1103.943/1396.504)*850e6 #mean -> A-basis
    yt[1]  = min(1.0,39.226/52.975)*40e6 #mean -> A-basis
    yc[1]  = min(1.0,235.434/282.439)*200e6 #mean -> A-basis
    s[1]   = min(1.0,142.411/159.516)*60e6 #mean -> A-basis
    plythickness[1] = 0.152e-3
    # "standard_uni"
    e1[2]  = 135.0e9
    e2[2]  = 10.0e9
    g12[2] = 5.0e9
    anu[2] = 0.30
    rho[2] = 1600.0
    xt[2]  = min(1.0,1355.656/1682.011)*1500e6 #mean -> A-basis
    xc[2]  = min(1.0,1103.943/1396.504)*1200e6 #mean -> A-basis
    yt[2]  = min(1.0,39.226/52.975)*50e6 #mean -> A-basis
    yc[2]  = min(1.0,235.434/282.439)*250e6 #mean -> A-basis
    s[2]   = min(1.0,142.411/159.516)*70e6 #mean -> A-basis
    plythickness[2] = 0.152e-3
    # "MR60H"
    e1[3]  = (165e9+150e9)/2.0
    e2[3]  = 8.56e9
    g12[3] = 4.39e9
    anu[3] = 0.326
    rho[3] = 1810.0
    xt[3]  = min(1.0,1355.656/1682.011)*3190e6 #mean -> A-basis
    xc[3]  = min(1.0,1103.943/1396.504)*1440e6 #mean -> A-basis
    yt[3]  = min(1.0,39.226/52.975)*82.0e6 #mean -> A-basis
    yc[3]  = min(1.0,235.434/282.439)*200.0e6 #mean -> A-basis
    s[3]   = min(1.0,142.411/159.516)*141e6 #mean -> A-basis
    plythickness[3] = 0.152e-3
    # "T3900_uni"
    e1[4]  = (148e9+131e9)/2.0
    e2[4]  = (9.7e9+9.7e9)/2.0
    g12[4] = 4.83e9
    anu[4] = 0.33
    rho[4] = 1573.0
    xt[4]  = min(1.0,1355.656/1682.011)*2830e6 #CTD mean -> A-basis
    xc[4]  = min(1.0,1103.943/1396.504)*1772e6 #CTD mean -> A-basis
    yt[4]  = min(1.0,39.226/52.975)*56.9e6 #CTD mean -> A-basis
    yc[4]  = min(1.0,235.434/282.439)*303e6 #CTD mean -> A-basis
    s[4]   = min(1.0,142.411/159.516)*89.6e6 #CTD mean -> A-basis
    plythickness[4] = 0.191e-3
    # "T700_uni"
    e1[5]=120.8e9
    e2[5]=11.57e9
    g12[5]=5.219e9
    anu[5]=0.350
    rho[5]=1525.0
    xt[5]=1356e6
    xc[5]=1104e6
    yt[5]=39.23e6
    yc[5]=235.4e6
    s[5]=142.4e6
    plythickness[5] = 0.152e-3
    # "ELT5500"
    e1[6]=41.8e9
    e2[6]=14.00e9
    g12[6]=2.630e9
    anu[6]=0.280
    rho[6]=1920.0 #g/cc * 1000
    xt[6]=972.0e6
    xc[6]=702.0e6
    yt[6]=100.0e6 #made up
    yc[6]=100.0e6 #made up
    s[6]=100.0e6 #made up
    plythickness[6] = 0.47e-3
    # "UDCarbon"
    e1[7]=114.5e9
    e2[7]=8.39e9
    g12[7]=5.990e9
    anu[7]=0.270
    rho[7]=1220.0
    xt[7]=1546.0e6
    xc[7]=1047.0e6
    yt[7]=100.0e6 #made up
    yc[7]=100.0e6 #made up
    s[7]=100.0e6 #made up
    plythickness[7] = 0.47e-3

    # FABRICS
    # "highmodulus_weave"
    e1[8]  = 85.0e9
    e2[8]  = 85.0e9
    g12[8] = 5.0e9
    anu[8] = 0.10
    rho[8] = 1600.0
    xt[8]  = min(1.0,1355.656/1682.011)*350e6 #mean -> A-basis
    xc[8]  = min(1.0,1103.943/1396.504)*150e6 #mean -> A-basis
    yt[8]  = min(1.0,39.226/52.975)*350e6 #mean -> A-basis
    yc[8]  = min(1.0,235.434/282.439)*150e6 #mean -> A-basis
    s[8]   = min(1.0,142.411/159.516)*35e6 #mean -> A-basis
    plythickness[8] = 0.218e-3
    # "standard_weave"
    e1[9]  = 70.0e9
    e2[9]  = 70.0e9
    g12[9] = 5.0e9
    anu[9] = 0.10
    rho[9] = 1600.0
    xt[9]  = min(1.0,1355.656/1682.011)*600e6 #mean -> A-basis
    xc[9]  = min(1.0,1103.943/1396.504)*570e6 #mean -> A-basis
    yt[9]  = min(1.0,39.226/52.975)*600e6 #mean -> A-basis
    yc[9]  = min(1.0,235.434/282.439)*570e6 #mean -> A-basis
    s[9]   = min(1.0,142.411/159.516)*90e6 #mean -> A-basis
    plythickness[9] = 0.218e-3
    # "T3900_weave"
    e1[10]  = (70.3e9+71.0e9)/2.0
    e2[10]  = (68.9e9+67.6e9)/2.0
    g12[10] = 4.6e9
    anu[10] = 0.032
    rho[10] = 1551.0
    xt[10]  = min(1.0,701.302/803.236)*1055e6 #CTD mean -> A-basis
    xc[10]  = min(1.0,549.748/749.955)*676e6 #CTD mean -> A-basis
    yt[10]  = min(1.0,557.575/722.602)*945e6 #CTD mean -> A-basis
    yc[10]  = min(1.0,604.067/741.866)*614e6 #CTD mean -> A-basis
    s[10]   = min(1.0,138.440/154.888)*79.3e6 #CTD mean -> A-basis
    plythickness[10] = 0.218e-3
    # "T700_weave"
    e1[11]=55.82e9
    e2[11]=52.10e9
    g12[11]=4.295e9
    anu[11]=0.085
    rho[11]=1501.0
    xt[11]=701.32e6
    yt[11]=557.59e6
    xc[11]=549.77e6
    yc[11]=604.08e6
    s[11]=138.44e6
    plythickness[11] = 0.218e-3
    # "Gelcoat"
    e1[12]=3.44e9
    e2[12]=3.44e9
    g12[12]=1.38e9
    anu[12]=0.3
    rho[12]=1235.0
    xt[12]=100.0e6 #made up
    xc[12]=100.0e6 #made up
    yt[12]=100.0e6 #made up
    yc[12]=100.0e6 #made up
    s[12]=100.0e6 #made up
    plythickness[12] = 0.05e-3
    # "Triax"
    e1[13]=27.7e9
    e2[13]=13.65e9
    g12[13]=7.2e9
    anu[13]=0.39
    rho[13]=1850.0
    xt[13]=700.0e6
    xc[13]=100.0e6 #made up
    yt[13]=700.0e6
    yc[13]=100.0e6 #made up
    s[13]=100.0e6 #made up
    plythickness[13] = 0.94e-3
    # "Saertex"
    e1[14]=13.6e9
    e2[14]=13.3e9
    g12[14]=11.8e9
    anu[14]=0.49
    rho[14]=1780.0
    xt[14]=144.0e6
    xc[14]=213.0e6
    yt[14]=144.0e6
    yc[14]=213.0e6
    s[14]=100.0e6 #made up
    plythickness[14] = 1.0e-3

    #FOAMS
    # "taylor_foam"
    e1[15]=48.0e6
    e2[15]=48.0e6
    g12[15]=28.0e6
    anu[15]=0.3
    rho[15]=75.0
    xt[15]=100.0e6 #made up
    yt[15]=100.0e6 #made up
    xc[15]=100.0e6 #made up
    yc[15]=100.0e6 #made up
    s[15]=100.0e6 #made up
    plythickness[15]=1.0E-3
    # "SNL_foam"
    e1[16]=256.0e6
    e2[16]=256.0e6
    g12[16]=22.0e6
    anu[16]=0.3
    rho[16]=200.0
    xt[16]=100.0e6 #made up
    yt[16]=100.0e6 #made up
    xc[16]=100.0e6 #made up
    yc[16]=100.0e6 #made up
    s[16]=100.0e6 #made up
    plythickness[16]=1.0E-3

    # Metals
    # "aluminum6061" http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061T6
    e1[17]=68.9e9
    e2[17]=68.9e9
    g12[17]=26.0e9
    anu[17]=0.33
    rho[17]=2700.0 #g/cc * 1000
    xt[17]=276.0e6 #use yield for metal?
    yt[17]=276.0e6 #use yield for metal?
    xc[17]=386.0e6 #use yield for metal?
    yc[17]=386.0e6 #use yield for metal?
    s[17]=207.0e6
    plythickness[17]=1.0E-3

    # "aluminum6063" http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T6
    e1[18]=68.9e9
    e2[18]=68.9e9
    g12[18]=25.8e9
    anu[18]=0.33
    rho[18]=2700.0 #g/cc * 1000
    xt[18]=214.0e6 #use yield for metal?
    yt[18]=214.0e6 #use yield for metal?
    xc[18]=276.0e6 #use yield for metal?
    yc[18]=276.0e6 #use yield for metal?
    s[18]=152.0e6
    plythickness[18]=1.0E-3

    # "CLA_5500"
    e1[19]=98.240e9 #pa
    e2[19]=5.102e9 #pa
    g12[19]=4.274e9 #pa
    anu[19]=0.3 #ratio
    rho[19]=1540.0 #g/cc * 1000 #kg/m3
    xt[19]=875.634139e6 #pa
    xc[19]=592.949102e6 #pa
    yt[19]=100.0e6 #made up
    yc[19]=100.0e6 #made up
    s[19]=100.0e6 #made up
    plythickness[19]=0.66E-3 #meters

    # "CBX_2400"
    e1[20]=14.931e9 #pa
    e2[20]=14.931e9 #pa
    g12[20]=23.890e9 #pa
    anu[20]=0.3 #ratio
    rho[20]=1530.0 #g/cc * 1000 #kg/m3
    xt[20]=455.053962e6 #pa
    xc[20]=455.053962e6 #pa
    yt[20]=100.0e6 #made up
    yc[20]=100.0e6 #made up
    s[20]=100.0e6 #made up
    plythickness[20]=0.81E-3 #meters

    # "ETLX_2400"
    e1[21]=20.333e9 #pa
    e2[21]=9.305e9 #pa
    g12[21]=4.756e9 #pa
    anu[21]=0.3 #ratio
    rho[21]=1900.0 #g/cc * 1000 #kg/m3
    xt[21]=530.896289e6 #pa
    xc[21]=530.896289e6 #pa
    yt[21]=100.0e6 #made up
    yc[21]=100.0e6 #made up
    s[21]=100.0e6 #made up
    plythickness[21]=0.66E-3 #meters

    # "Airex_C70_55"
    e1[22]=0.045e9 #pa
    e2[22]=0.045e9 #pa
    g12[22]=0.022e9 #pa
    anu[22]=0.2 #ratio
    rho[22]=59.0 #g/cc * 1000 #kg/m3
    xt[22]=100.0e6 #pa #made up
    xc[22]=100.0e6 #pa
    yt[22]=100.0e6 #made up
    yc[22]=100.0e6 #made up
    s[22]=100.0e6 #made up
    plythickness[22]=1.0E-3 #meters

    # "EBX_2400_x10"
    e1[23]=982.400e9 #pa
    e2[23]=51.020e9 #pa
    g12[23]=42.740e9 #pa
    anu[23]=0.3 #ratio
    rho[23]=15300.0 #g/cc * 1000 #kg/m3
    xt[23]=4550.53962e6 #pa
    xc[23]=4550.53962e6 #pa
    yt[23]=100.0e6 #made up
    yc[23]=100.0e6 #made up
    s[23]=100.0e6 #made up
    plythickness[23]=0.07E-3 #meters

    # "ETLX_2400_x10"
    e1[24]=149.310e9 #pa
    e2[24]=149.310e9 #pa
    g12[24]=238.900e9 #pa
    anu[24]=0.3 #ratio
    rho[24]=19000.0 #g/cc * 1000 #kg/m3
    xt[24]=5308.96289e6 #pa
    xc[24]=5308.96289e6 #pa
    yt[24]=100.0e6 #made up
    yc[24]=100.0e6 #made up
    s[24]=100.0e6 #made up
    plythickness[24]=0.08E-3 #meters

    # "Airex_C70_55_x10"
    e1[25]=203.335e9 #pa
    e2[25]=93.051e9 #pa
    g12[25]=47.560e9 #pa
    anu[25]=0.2 #ratio
    rho[25]=590.0 #g/cc * 1000 #kg/m3
    xt[25]=100.0e6 #pa #made up
    xc[25]=100.0e6 #pa
    yt[25]=100.0e6 #made up
    yc[25]=100.0e6 #made up
    s[25]=100.0e6 #made up
    plythickness[25]=0.07E-3 #meters

    return plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness))
end
