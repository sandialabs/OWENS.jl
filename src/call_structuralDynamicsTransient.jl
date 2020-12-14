
function call_structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t_in,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)
    #Break Out Model
    # start = time()
    analysisType = model.analysisType
    turbineStartup = model.turbineStartup
    usingRotorSpeedFunction = model.usingRotorSpeedFunction
    tocp = model.tocp
    initCond = model.initCond
    initCond_size = zeros(Int,2)
    initCond_size[:] .= size(initCond)
    numTS = model.numTS
    delta_t = model.delta_t
    Omegaocp = model.Omegaocp
    aeroElasticOn = model.aeroElasticOn
    aeroForceOn = model.aeroForceOn
    aeroLoadsOn = model.aeroLoadsOn
    driveTrainOn = model.driveTrainOn
    airDensity = model.airDensity
    guessFreq = model.guessFreq
    gravityOn = model.gravityOn
    generatorOn = model.generatorOn
    hydroOn = model.hydroOn
    JgearBox = model.JgearBox
    gearRatio = model.gearRatio
    gearBoxEfficiency = model.gearBoxEfficiency
    useGeneratorFunction = model.useGeneratorFunction
    generatorProps = model.generatorProps
    OmegaGenStart = model.OmegaGenStart
    omegaControl = model.omegaControl
    OmegaInit = model.OmegaInit
    totalNumDof = model.totalNumDof
    spinUpOn = model.spinUpOn
    nlOn = model.nlOn
    numModesToExtract = model.numModesToExtract
    aeroloadfile = model.aeroloadfile
    owensfile = model.owensfile
    outFilename = model.outFilename
    RayleighAlpha = model.RayleighAlpha
    RayleighBeta = model.RayleighBeta
    elementOrder = model.elementOrder
    joint = model.joint
    platformTurbineConnectionNodeNumber = model.platformTurbineConnectionNodeNumber
    jointTransform = model.jointTransform
    reducedDOFList = model.reducedDOFList

    bladeData = model.bladeData #Struct
    numBlades = bladeData.numBlades
    bladeNum = bladeData.bladeNum
    h = bladeData.h
    nodeNum = bladeData.nodeNum
    elementNum = bladeData.elementNum
    remaining = bladeData.remaining

    nlParams = model.nlParams #Struct
    iterationType = nlParams.iterationType
    adaptiveLoadSteppingFlag = nlParams.adaptiveLoadSteppingFlag
    tolerance = nlParams.tolerance
    maxIterations = nlParams.maxIterations
    maxNumLoadSteps = nlParams.maxNumLoadSteps
    minLoadStepDelta = nlParams.minLoadStepDelta
    minLoadStep = nlParams.minLoadStep
    prescribedLoadStep = nlParams.prescribedLoadStep

    BC = model.BC #Struct
    numpBC = BC.numpBC
    pBC = BC.pBC
    numsBC = BC.numsBC
    nummBC = BC.nummBC
    isConstrained = BC.isConstrained
    map = BC.map
    redVectorMap = BC.redVectorMap

    nodalTerms = model.nodalTerms #Struct
    concLoad = nodalTerms.concLoad
    mysize = length(concLoad)
    concLoadnodeNum = zeros(mysize)
    concLoaddof = zeros(mysize)
    concLoadval = zeros(mysize)
    for ii = 1:mysize
        concLoadnodeNum[ii] = concLoad[ii].nodeNum
        concLoaddof[ii] = concLoad[ii].dof
        concLoadval[ii] = concLoad[ii].val
    end

    concStiff = nodalTerms.concStiff
    mysize = length(concStiff)
    concStiffnodeNum = zeros(mysize)
    concStiffdof = zeros(mysize)
    concStiffval = zeros(mysize)
    for ii = 1:mysize
        concStiffnodeNum[ii] = concStiff[ii].nodeNum
        concStiffdof[ii] = concStiff[ii].dof
        concStiffval[ii] = concStiff[ii].val
    end

    concMass = nodalTerms.concMass
    mysize = length(concMass)
    concMassnodeNum = zeros(mysize)
    concMassdof = zeros(mysize)
    concMassval = zeros(mysize)
    for ii = 1:mysize
        concMassnodeNum[ii] = concMass[ii].nodeNum
        concMassdof[ii] = concMass[ii].dof
        concMassval[ii] = concMass[ii].val
    end

    concStiffGen = nodalTerms.concStiffGen
    mysize = length(concStiffGen)
    StiffGennodeNum = zeros(mysize)
    StiffGendof1 = zeros(mysize)
    StiffGendof2 = zeros(mysize)
    StiffGenval = zeros(mysize)
    for ii = 1:mysize
        StiffGennodeNum[ii] = concStiffGen[ii].nodeNum
        StiffGendof1[ii] = concStiffGen[ii].dof1
        StiffGendof2[ii] = concStiffGen[ii].dof2
        StiffGenval[ii] = concStiffGen[ii].val
    end

    concMassGen = nodalTerms.concMassGen
    mysize = length(concMassGen)
    MassGennodeNum = zeros(mysize)
    MassGendof1 = zeros(mysize)
    MassGendof2 = zeros(mysize)
    MassGenval = zeros(mysize)
    for ii = 1:mysize
        MassGennodeNum[ii] = concMassGen[ii].nodeNum
        MassGendof1[ii] = concMassGen[ii].dof1
        MassGendof2[ii] = concMassGen[ii].dof2
        MassGenval[ii] = concMassGen[ii].val
    end


    concDampGen = nodalTerms.concDampGen
    mysize = length(concDampGen)
    DampGennodeNum = zeros(mysize)
    DampGendof1 = zeros(mysize)
    DampGendof2 = zeros(mysize)
    DampGenval = zeros(mysize)
    for ii = 1:mysize
        DampGennodeNum[ii] = concDampGen[ii].nodeNum
        DampGendof1[ii] = concDampGen[ii].dof1
        DampGendof2[ii] = concDampGen[ii].dof2
        DampGenval[ii] = concDampGen[ii].val
    end


    driveShaftProps = model.driveShaftProps #Struct
    k = driveShaftProps.k
    c = driveShaftProps.c

    #Break out Mesh
    nodeNummesh = mesh.nodeNum
    numEl = mesh.numEl
    numNodes = mesh.numNodes
    x = mesh.x
    y = mesh.y
    z = mesh.z
    elNum = mesh.elNum
    conn = mesh.conn

    #Break out el
    props = el.props #struct
    ac = zeros(length(props),2)
    twist = zeros(length(props),2)
    rhoA = zeros(length(props),2)
    EIyy = zeros(length(props),2)
    EIzz = zeros(length(props),2)
    GJ = zeros(length(props),2)
    EA = zeros(length(props),2)
    rhoIyy = zeros(length(props),2)
    rhoIzz = zeros(length(props),2)
    rhoJ = zeros(length(props),2)
    zcm = zeros(length(props),2)
    ycm = zeros(length(props),2)
    a = zeros(length(props),2)
    EIyz = zeros(length(props),2)
    alpha1 = zeros(length(props),2)
    alpha2 = zeros(length(props),2)
    alpha3 = zeros(length(props),2)
    alpha4 = zeros(length(props),2)
    alpha5 = zeros(length(props),2)
    alpha6 = zeros(length(props),2)
    rhoIyz = zeros(length(props),2)
    b = zeros(length(props),2)
    a0 = zeros(length(props),2)
    aeroCenterOffset = zeros(length(props),2)

    for jj = 1:length(props)
        ac[jj,:] = props[jj].ac
        twist[jj,:] = props[jj].twist
        rhoA[jj,:] = props[jj].rhoA
        EIyy[jj,:] = props[jj].EIyy
        EIzz[jj,:] = props[jj].EIzz
        GJ[jj,:] = props[jj].GJ
        EA[jj,:] = props[jj].EA
        rhoIyy[jj,:] = props[jj].rhoIyy
        rhoIzz[jj,:] = props[jj].rhoIzz
        rhoJ[jj,:] = props[jj].rhoJ
        zcm[jj,:] = props[jj].zcm
        ycm[jj,:] = props[jj].ycm
        a[jj,:] = props[jj].a
        EIyz[jj,:] = props[jj].EIyz
        alpha1[jj,:] = props[jj].alpha1
        alpha2[jj,:] = props[jj].alpha2
        alpha3[jj,:] = props[jj].alpha3
        alpha4[jj,:] = props[jj].alpha4
        alpha5[jj,:] = props[jj].alpha5
        alpha6[jj,:] = props[jj].alpha6
        rhoIyz[jj,:] = props[jj].rhoIyz
        b[jj,:] = props[jj].b
        a0[jj,:] = props[jj].a0
        aeroCenterOffset[jj,:] = props[jj].aeroCenterOffset
    end

    elLen = el.elLen
    psi = el.psi
    theta = el.theta
    roll = el.roll
    rotationalEffects = el.rotationalEffects

    # Break out dispData
    displ_s = dispData.displ_s
    displdot_s = dispData.displdot_s
    displddot_s = dispData.displddot_s

    # Break out elStorage
    K11 = zeros(length(elStorage),4)
    K12 = zeros(length(elStorage),4)
    K13 = zeros(length(elStorage),4)
    K14 = zeros(length(elStorage),4)
    K15 = zeros(length(elStorage),4)
    K16 = zeros(length(elStorage),4)
    K22 = zeros(length(elStorage),4)
    K23 = zeros(length(elStorage),4)
    K24 = zeros(length(elStorage),4)
    K25 = zeros(length(elStorage),4)
    K26 = zeros(length(elStorage),4)
    K33 = zeros(length(elStorage),4)
    K34 = zeros(length(elStorage),4)
    K35 = zeros(length(elStorage),4)
    K36 = zeros(length(elStorage),4)
    K44 = zeros(length(elStorage),4)
    K45 = zeros(length(elStorage),4)
    K46 = zeros(length(elStorage),4)
    K55 = zeros(length(elStorage),4)
    K56 = zeros(length(elStorage),4)
    K66 = zeros(length(elStorage),4)
    M11 = zeros(length(elStorage),4)
    M15 = zeros(length(elStorage),4)
    M16 = zeros(length(elStorage),4)
    M22 = zeros(length(elStorage),4)
    M24 = zeros(length(elStorage),4)
    M33 = zeros(length(elStorage),4)
    M34 = zeros(length(elStorage),4)
    M44 = zeros(length(elStorage),4)
    M55 = zeros(length(elStorage),4)
    M56 = zeros(length(elStorage),4)
    M66 = zeros(length(elStorage),4)
    S11 = zeros(length(elStorage),4)
    S12 = zeros(length(elStorage),4)
    S13 = zeros(length(elStorage),4)
    S15 = zeros(length(elStorage),4)
    S16 = zeros(length(elStorage),4)
    S22 = zeros(length(elStorage),4)
    S23 = zeros(length(elStorage),4)
    S25 = zeros(length(elStorage),4)
    S26 = zeros(length(elStorage),4)
    S33 = zeros(length(elStorage),4)
    S35 = zeros(length(elStorage),4)
    S36 = zeros(length(elStorage),4)
    S55 = zeros(length(elStorage),4)
    S56 = zeros(length(elStorage),4)
    S66 = zeros(length(elStorage),4)
    S14_1 = zeros(length(elStorage),4)
    S14_2 = zeros(length(elStorage),4)
    S24_1 = zeros(length(elStorage),4)
    S24_2 = zeros(length(elStorage),4)
    S34_1 = zeros(length(elStorage),4)
    S34_2 = zeros(length(elStorage),4)
    S45_1 = zeros(length(elStorage),4)
    S45_2 = zeros(length(elStorage),4)
    S46_1 = zeros(length(elStorage),4)
    S46_2 = zeros(length(elStorage),4)
    S44_1 = zeros(length(elStorage),4)
    S44_2 = zeros(length(elStorage),4)
    S44_3 = zeros(length(elStorage),4)
    C12 = zeros(length(elStorage),4)
    C13 = zeros(length(elStorage),4)
    C23 = zeros(length(elStorage),4)
    C24 = zeros(length(elStorage),4)
    C25 = zeros(length(elStorage),4)
    C26 = zeros(length(elStorage),4)
    C34 = zeros(length(elStorage),4)
    C35 = zeros(length(elStorage),4)
    C36 = zeros(length(elStorage),4)
    C14_1 = zeros(length(elStorage),4)
    C14_2 = zeros(length(elStorage),4)
    C45_1 = zeros(length(elStorage),4)
    C45_2 = zeros(length(elStorage),4)
    C46_1 = zeros(length(elStorage),4)
    C46_2 = zeros(length(elStorage),4)
    mel = zeros(length(elStorage))
    moiel = zeros(length(elStorage),9)
    xmel = zeros(length(elStorage),3)

    for jj = 1:length(elStorage)
        K11[jj,:] = [elStorage[jj].K11[1,:];elStorage[jj].K11[2,:]]
        K12[jj,:] = [elStorage[jj].K12[1,:];elStorage[jj].K12[2,:]]
        K13[jj,:] = [elStorage[jj].K13[1,:];elStorage[jj].K13[2,:]]
        K14[jj,:] = [elStorage[jj].K14[1,:];elStorage[jj].K14[2,:]]
        K15[jj,:] = [elStorage[jj].K15[1,:];elStorage[jj].K15[2,:]]
        K16[jj,:] = [elStorage[jj].K16[1,:];elStorage[jj].K16[2,:]]
        K22[jj,:] = [elStorage[jj].K22[1,:];elStorage[jj].K22[2,:]]
        K23[jj,:] = [elStorage[jj].K23[1,:];elStorage[jj].K23[2,:]]
        K24[jj,:] = [elStorage[jj].K24[1,:];elStorage[jj].K24[2,:]]
        K25[jj,:] = [elStorage[jj].K25[1,:];elStorage[jj].K25[2,:]]
        K26[jj,:] = [elStorage[jj].K26[1,:];elStorage[jj].K26[2,:]]
        K33[jj,:] = [elStorage[jj].K33[1,:];elStorage[jj].K33[2,:]]
        K34[jj,:] = [elStorage[jj].K34[1,:];elStorage[jj].K34[2,:]]
        K35[jj,:] = [elStorage[jj].K35[1,:];elStorage[jj].K35[2,:]]
        K36[jj,:] = [elStorage[jj].K36[1,:];elStorage[jj].K36[2,:]]
        K44[jj,:] = [elStorage[jj].K44[1,:];elStorage[jj].K44[2,:]]
        K45[jj,:] = [elStorage[jj].K45[1,:];elStorage[jj].K45[2,:]]
        K46[jj,:] = [elStorage[jj].K46[1,:];elStorage[jj].K46[2,:]]
        K55[jj,:] = [elStorage[jj].K55[1,:];elStorage[jj].K55[2,:]]
        K56[jj,:] = [elStorage[jj].K56[1,:];elStorage[jj].K56[2,:]]
        K66[jj,:] = [elStorage[jj].K66[1,:];elStorage[jj].K66[2,:]]
        M11[jj,:] = [elStorage[jj].M11[1,:];elStorage[jj].M11[2,:]]
        M15[jj,:] = [elStorage[jj].M15[1,:];elStorage[jj].M15[2,:]]
        M16[jj,:] = [elStorage[jj].M16[1,:];elStorage[jj].M16[2,:]]
        M22[jj,:] = [elStorage[jj].M22[1,:];elStorage[jj].M22[2,:]]
        M24[jj,:] = [elStorage[jj].M24[1,:];elStorage[jj].M24[2,:]]
        M33[jj,:] = [elStorage[jj].M33[1,:];elStorage[jj].M33[2,:]]
        M34[jj,:] = [elStorage[jj].M34[1,:];elStorage[jj].M34[2,:]]
        M44[jj,:] = [elStorage[jj].M44[1,:];elStorage[jj].M44[2,:]]
        M55[jj,:] = [elStorage[jj].M55[1,:];elStorage[jj].M55[2,:]]
        M56[jj,:] = [elStorage[jj].M56[1,:];elStorage[jj].M56[2,:]]
        M66[jj,:] = [elStorage[jj].M66[1,:];elStorage[jj].M66[2,:]]
        S11[jj,:] = [elStorage[jj].S11[1,:];elStorage[jj].S11[2,:]]
        S12[jj,:] = [elStorage[jj].S12[1,:];elStorage[jj].S12[2,:]]
        S13[jj,:] = [elStorage[jj].S13[1,:];elStorage[jj].S13[2,:]]
        S15[jj,:] = [elStorage[jj].S15[1,:];elStorage[jj].S15[2,:]]
        S16[jj,:] = [elStorage[jj].S16[1,:];elStorage[jj].S16[2,:]]
        S22[jj,:] = [elStorage[jj].S22[1,:];elStorage[jj].S22[2,:]]
        S23[jj,:] = [elStorage[jj].S23[1,:];elStorage[jj].S23[2,:]]
        S25[jj,:] = [elStorage[jj].S25[1,:];elStorage[jj].S25[2,:]]
        S26[jj,:] = [elStorage[jj].S26[1,:];elStorage[jj].S26[2,:]]
        S33[jj,:] = [elStorage[jj].S33[1,:];elStorage[jj].S33[2,:]]
        S35[jj,:] = [elStorage[jj].S35[1,:];elStorage[jj].S35[2,:]]
        S36[jj,:] = [elStorage[jj].S36[1,:];elStorage[jj].S36[2,:]]
        S55[jj,:] = [elStorage[jj].S55[1,:];elStorage[jj].S55[2,:]]
        S56[jj,:] = [elStorage[jj].S56[1,:];elStorage[jj].S56[2,:]]
        S66[jj,:] = [elStorage[jj].S66[1,:];elStorage[jj].S66[2,:]]
        S14_1[jj,:] = [elStorage[jj].S14_1[1,:];elStorage[jj].S14_1[2,:]]
        S14_2[jj,:] = [elStorage[jj].S14_2[1,:];elStorage[jj].S14_2[2,:]]
        S24_1[jj,:] = [elStorage[jj].S24_1[1,:];elStorage[jj].S24_1[2,:]]
        S24_2[jj,:] = [elStorage[jj].S24_2[1,:];elStorage[jj].S24_2[2,:]]
        S34_1[jj,:] = [elStorage[jj].S34_1[1,:];elStorage[jj].S34_1[2,:]]
        S34_2[jj,:] = [elStorage[jj].S34_2[1,:];elStorage[jj].S34_2[2,:]]
        S45_1[jj,:] = [elStorage[jj].S45_1[1,:];elStorage[jj].S45_1[2,:]]
        S45_2[jj,:] = [elStorage[jj].S45_2[1,:];elStorage[jj].S45_2[2,:]]
        S46_1[jj,:] = [elStorage[jj].S46_1[1,:];elStorage[jj].S46_1[2,:]]
        S46_2[jj,:] = [elStorage[jj].S46_2[1,:];elStorage[jj].S46_2[2,:]]
        S44_1[jj,:] = [elStorage[jj].S44_1[1,:];elStorage[jj].S44_1[2,:]]
        S44_2[jj,:] = [elStorage[jj].S44_2[1,:];elStorage[jj].S44_2[2,:]]
        S44_3[jj,:] = [elStorage[jj].S44_3[1,:];elStorage[jj].S44_3[2,:]]
        C12[jj,:] = [elStorage[jj].C12[1,:];elStorage[jj].C12[2,:]]
        C13[jj,:] = [elStorage[jj].C13[1,:];elStorage[jj].C13[2,:]]
        C23[jj,:] = [elStorage[jj].C23[1,:];elStorage[jj].C23[2,:]]
        C24[jj,:] = [elStorage[jj].C24[1,:];elStorage[jj].C24[2,:]]
        C25[jj,:] = [elStorage[jj].C25[1,:];elStorage[jj].C25[2,:]]
        C26[jj,:] = [elStorage[jj].C26[1,:];elStorage[jj].C26[2,:]]
        C34[jj,:] = [elStorage[jj].C34[1,:];elStorage[jj].C34[2,:]]
        C35[jj,:] = [elStorage[jj].C35[1,:];elStorage[jj].C35[2,:]]
        C36[jj,:] = [elStorage[jj].C36[1,:];elStorage[jj].C36[2,:]]
        C14_1[jj,:] = [elStorage[jj].C14_1[1,:];elStorage[jj].C14_1[2,:]]
        C14_2[jj,:] = [elStorage[jj].C14_2[1,:];elStorage[jj].C14_2[2,:]]
        C45_1[jj,:] = [elStorage[jj].C45_1[1,:];elStorage[jj].C45_1[2,:]]
        C45_2[jj,:] = [elStorage[jj].C45_2[1,:];elStorage[jj].C45_2[2,:]]
        C46_1[jj,:] = [elStorage[jj].C46_1[1,:];elStorage[jj].C46_1[2,:]]
        C46_2[jj,:] = [elStorage[jj].C46_2[1,:];elStorage[jj].C46_2[2,:]]
        mel[jj] = elStorage[jj].mel
        moiel[jj,:] = [elStorage[jj].moiel[1,:];elStorage[jj].moiel[2,:];elStorage[jj].moiel[3,:]]
        xmel[jj,:] = elStorage[jj].xmel[:]
    end

    # Initialize outputs, which are being changed by reference

    displ_sp1 = zeros(model.totalNumDof)
    displddot_sp1 = zeros(model.totalNumDof)
    displdot_sp1 = zeros(model.totalNumDof)
    eps_xx_0 = zeros(mesh.numEl*4)
    eps_xx_z = zeros(mesh.numEl*4)
    eps_xx_y = zeros(mesh.numEl*4)
    gam_xz_0 = zeros(mesh.numEl*4)
    gam_xz_y = zeros(mesh.numEl*4)
    gam_xy_0 = zeros(mesh.numEl*4)
    gam_xy_z = zeros(mesh.numEl*4)
    FReaction_j = zeros(6)

    # ccall((:structuralDynamicsTransient,"/Users/kevmoor/Documents/coderepos/OWENS.jl/src/Matlab_Cxx/codegen/dll/structuralDynamicsTransient/structuralDynamicsTransient"),Cvoid,
    # (Cint, #const boolean_T rotationalEffects
    # Cint, #boolean_T gravityOn
    # Cint, #boolean_T nlOn
    # Cdouble, #double maxNumLoadSteps
    # Cdouble, #double airDensity
    # Cdouble, #double RayleighAlpha
    # Cdouble, #double RayleighBeta
    # Cdouble, #double elementOrder
    # Cdouble, #double delta_t
    # Cdouble, #double c_platformTurbineConnectionNode
    # Cdouble, #double tolerance
    # Cdouble, #double maxIterations
    # Cdouble, #double numpBC
    # Cdouble, #double numEl
    # Cdouble, #double Omega
    # Cdouble, #double OmegaDot
    # Ptr{Cdouble}, #const double mel[75]
    # Ptr{Cdouble}, #const double moiel[675]
    # Ptr{Cdouble}, #const double xmel[225]
    # Ptr{Cdouble}, #const double Fexternal[273]
    # Ptr{Cdouble}, #const double Fdof[273]
    # Ptr{Cdouble}, #const double CN2H[9]
    # Ptr{Cdouble}, #const double rbData[9]
    # Ptr{Cdouble}, #const double jointTransform[206640]
    # Ptr{Cdouble}, #const double conn[150]
    # Ptr{Cdouble}, #const double joint[96]
    # Ptr{Cdouble}, #const double pBC[18]
    # Ptr{Cdouble}, #const double x[82]
    # Ptr{Cdouble}, #const double y[82]
    # Ptr{Cdouble}, #const double z[82]
    # Ptr{Cdouble}, #const double ac[150]
    # Ptr{Cdouble}, #const double twist[150]
    # Ptr{Cdouble}, #const double rhoA[150]
    # Ptr{Cdouble}, #const double EIyy[150]
    # Ptr{Cdouble}, #const double EIzz[150]
    # Ptr{Cdouble}, #const double GJ[150]
    # Ptr{Cdouble}, #const double EA[150]
    # Ptr{Cdouble}, #const double rhoIyy[150]
    # Ptr{Cdouble}, #const double rhoIzz[150]
    # Ptr{Cdouble}, #const double rhoJ[150]
    # Ptr{Cdouble}, #const double zcm[150]
    # Ptr{Cdouble}, #const double ycm[150]
    # Ptr{Cdouble}, #const double a[150]
    # Ptr{Cdouble}, #const double EIyz[150]
    # Ptr{Cdouble}, #const double alpha1[150]
    # Ptr{Cdouble}, #const double alpha2[150]
    # Ptr{Cdouble}, #const double alpha3[150]
    # Ptr{Cdouble}, #const double alpha4[150]
    # Ptr{Cdouble}, #const double alpha5[150]
    # Ptr{Cdouble}, #const double alpha6[150]
    # Ptr{Cdouble}, #const double rhoIyz[150]
    # Ptr{Cdouble}, #const double b[150]
    # Ptr{Cdouble}, #const double a0[150]
    # Ptr{Cdouble}, #const double aeroCenterOffset[150]
    # Ptr{Cdouble}, #const double elLen[75]
    # Ptr{Cdouble}, #const double psi[75]
    # Ptr{Cdouble}, #const double theta[75]
    # Ptr{Cdouble}, #const double roll[75]
    # Ptr{Cdouble}, #const double displ_s[492]
    # Ptr{Cdouble}, #const double displdot_s[492]
    # Ptr{Cdouble}, #const double displddot_s[492]
    # Ptr{Cdouble}, # concLoadnodeNum
    # Ptr{Cdouble}, # concLoaddof
    # Ptr{Cdouble}, # concLoadval
    # Ptr{Cdouble}, # concStiffnodeNum
    # Ptr{Cdouble}, # concStiffdof
    # Ptr{Cdouble}, # concStiffval
    # Ptr{Cdouble}, # concMassnodeNum
    # Ptr{Cdouble}, # concMassdof
    # Ptr{Cdouble}, # concMassval
    # Ptr{Cdouble}, # StiffGennodeNum
    # Ptr{Cdouble}, # StiffGendof1
    # Ptr{Cdouble}, # StiffGendof2
    # Ptr{Cdouble}, # StiffGenval
    # Ptr{Cdouble}, # MassGennodeNum
    # Ptr{Cdouble}, # MassGendof1
    # Ptr{Cdouble}, # MassGendof2
    # Ptr{Cdouble}, # MassGenval
    # Ptr{Cdouble}, # DampGennodeNum
    # Ptr{Cdouble}, # DampGendof1
    # Ptr{Cdouble}, # DampGendof2
    # Ptr{Cdouble}, # DampGenval
    # Ptr{Cdouble}, #const double K11[300]
    # Ptr{Cdouble}, #const double K12[300]
    # Ptr{Cdouble}, #const double K13[300]
    # Ptr{Cdouble}, #const double K14[300]
    # Ptr{Cdouble}, #const double K15[300]
    # Ptr{Cdouble}, #const double K16[300]
    # Ptr{Cdouble}, #const double K22[300]
    # Ptr{Cdouble}, #const double K23[300]
    # Ptr{Cdouble}, #const double K24[300]
    # Ptr{Cdouble}, #const double K25[300]
    # Ptr{Cdouble}, #const double K26[300]
    # Ptr{Cdouble}, #const double K33[300]
    # Ptr{Cdouble}, #const double K34[300]
    # Ptr{Cdouble}, #const double K35[300]
    # Ptr{Cdouble}, #const double K36[300]
    # Ptr{Cdouble}, #const double K44[300]
    # Ptr{Cdouble}, #const double K45[300]
    # Ptr{Cdouble}, #const double K46[300]
    # Ptr{Cdouble}, #const double K55[300]
    # Ptr{Cdouble}, #const double K56[300]
    # Ptr{Cdouble}, #const double K66[300]
    # Ptr{Cdouble}, #const double M11[300]
    # Ptr{Cdouble}, #const double M15[300]
    # Ptr{Cdouble}, #const double M16[300]
    # Ptr{Cdouble}, #const double M22[300]
    # Ptr{Cdouble}, #const double M24[300]
    # Ptr{Cdouble}, #const double M33[300]
    # Ptr{Cdouble}, #const double M34[300]
    # Ptr{Cdouble}, #const double M44[300]
    # Ptr{Cdouble}, #const double M55[300]
    # Ptr{Cdouble}, #const double M56[300]
    # Ptr{Cdouble}, #const double M66[300]
    # Ptr{Cdouble}, #const double S11[300]
    # Ptr{Cdouble}, #const double S12[300]
    # Ptr{Cdouble}, #const double S13[300]
    # Ptr{Cdouble}, #const double S15[300]
    # Ptr{Cdouble}, #const double S16[300]
    # Ptr{Cdouble}, #const double S22[300]
    # Ptr{Cdouble}, #const double S23[300]
    # Ptr{Cdouble}, #const double S25[300]
    # Ptr{Cdouble}, #const double S26[300]
    # Ptr{Cdouble}, #const double S33[300]
    # Ptr{Cdouble}, #const double S35[300]
    # Ptr{Cdouble}, #const double S36[300]
    # Ptr{Cdouble}, #const double S55[300]
    # Ptr{Cdouble}, #const double S56[300]
    # Ptr{Cdouble}, #const double S66[300]
    # Ptr{Cdouble}, #const double S14_1[300]
    # Ptr{Cdouble}, #const double S14_2[300]
    # Ptr{Cdouble}, #const double S24_1[300]
    # Ptr{Cdouble}, #const double S24_2[300]
    # Ptr{Cdouble}, #const double S34_1[300]
    # Ptr{Cdouble}, #const double S34_2[300]
    # Ptr{Cdouble}, #const double S45_1[300]
    # Ptr{Cdouble}, #const double S45_2[300]
    # Ptr{Cdouble}, #const double S46_1[300]
    # Ptr{Cdouble}, #const double S46_2[300]
    # Ptr{Cdouble}, #const double S44_1[300]
    # Ptr{Cdouble}, #const double S44_2[300]
    # Ptr{Cdouble}, #const double S44_3[300]
    # Ptr{Cdouble}, #const double C12[300]
    # Ptr{Cdouble}, #const double C13[300]
    # Ptr{Cdouble}, #const double C23[300]
    # Ptr{Cdouble}, #const double C24[300]
    # Ptr{Cdouble}, #const double C25[300]
    # Ptr{Cdouble}, #const double C26[300]
    # Ptr{Cdouble}, #const double C34[300]
    # Ptr{Cdouble}, #const double C35[300]
    # Ptr{Cdouble}, #const double C36[300]
    # Ptr{Cdouble}, #const double C14_1[300]
    # Ptr{Cdouble}, #const double C14_2[300]
    # Ptr{Cdouble}, #const double C45_1[300]
    # Ptr{Cdouble}, #const double C45_2[300]
    # Ptr{Cdouble}, #const double C46_1[300]
    # Ptr{Cdouble}, #const double C46_2[300]
    # Ptr{Cdouble}, #double displ_sp1[492]
    # Ptr{Cdouble}, #double displddot_sp1[492]
    # Ptr{Cdouble}, #double displdot_sp1[492]
    # Ptr{Cdouble}, #double eps_xx_0[300]
    # Ptr{Cdouble}, #double eps_xx_z[300]
    # Ptr{Cdouble}, #double eps_xx_y[300]
    # Ptr{Cdouble}, #double gam_xz_0[300]
    # Ptr{Cdouble}, #double gam_xz_y[300]
    # Ptr{Cdouble}, #double gam_xy_0[300]
    # Ptr{Cdouble}, #double gam_xy_z[300]
    # Ptr{Cdouble}), #double FReaction_sp1[6]
    # rotationalEffects,
    # gravityOn,
    # nlOn,
    # maxNumLoadSteps,
    # airDensity,
    # RayleighAlpha,
    # RayleighBeta,
    # elementOrder,
    # delta_t,
    # platformTurbineConnectionNodeNumber,
    # tolerance,
    # maxIterations,
    # numpBC,
    # numEl,
    # Omega_j,
    # OmegaDot_j,
    # mel,
    # moiel,
    # xmel,
    # Fexternal,
    # Fdof,
    # CN2H,
    # rbData,
    # jointTransform,
    # conn,
    # joint,
    # float(pBC),
    # x,
    # y,
    # z,
    # ac,
    # twist,
    # rhoA,
    # EIyy,
    # EIzz,
    # GJ,
    # EA,
    # rhoIyy,
    # rhoIzz,
    # rhoJ,
    # zcm,
    # ycm,
    # a,
    # EIyz,
    # alpha1,
    # alpha2,
    # alpha3,
    # alpha4,
    # alpha5,
    # alpha6,
    # rhoIyz,
    # b,
    # a0,
    # aeroCenterOffset,
    # elLen,
    # psi,
    # theta,
    # roll,
    # displ_s,
    # displdot_s,
    # displddot_s,
    # concLoadnodeNum,
    # concLoaddof,
    # concLoadval,
    # concStiffnodeNum,
    # concStiffdof,
    # concStiffval,
    # concMassnodeNum,
    # concMassdof,
    # concMassval,
    # StiffGennodeNum,
    # StiffGendof1,
    # StiffGendof2,
    # StiffGenval,
    # MassGennodeNum,
    # MassGendof1,
    # MassGendof2,
    # MassGenval,
    # DampGennodeNum,
    # DampGendof1,
    # DampGendof2,
    # DampGenval,
    # K11,
    # K12,
    # K13,
    # K14,
    # K15,
    # K16,
    # K22,
    # K23,
    # K24,
    # K25,
    # K26,
    # K33,
    # K34,
    # K35,
    # K36,
    # K44,
    # K45,
    # K46,
    # K55,
    # K56,
    # K66,
    # M11,
    # M15,
    # M16,
    # M22,
    # M24,
    # M33,
    # M34,
    # M44,
    # M55,
    # M56,
    # M66,
    # S11,
    # S12,
    # S13,
    # S15,
    # S16,
    # S22,
    # S23,
    # S25,
    # S26,
    # S33,
    # S35,
    # S36,
    # S55,
    # S56,
    # S66,
    # S14_1,
    # S14_2,
    # S24_1,
    # S24_2,
    # S34_1,
    # S34_2,
    # S45_1,
    # S45_2,
    # S46_1,
    # S46_2,
    # S44_1,
    # S44_2,
    # S44_3,
    # C12,
    # C13,
    # C23,
    # C24,
    # C25,
    # C26,
    # C34,
    # C35,
    # C36,
    # C14_1,
    # C14_2,
    # C45_1,
    # C45_2,
    # C46_1,
    # C46_2,
    # displ_sp1,
    # displddot_sp1,
    # displdot_sp1,
    # eps_xx_0,
    # eps_xx_z,
    # eps_xx_y,
    # gam_xz_0,
    # gam_xz_y,
    # gam_xy_0,
    # gam_xy_z,
    # FReaction_j)

    mat"""
    [$displ_sp1,$displddot_sp1,$displdot_sp1,$eps_xx_0,$eps_xx_z,$eps_xx_y,$gam_xz_0,$gam_xz_y,$gam_xy_0,$gam_xy_z,$FReaction_j] = structuralDynamicsTransient($rotationalEffects,...
    $gravityOn,...
    $nlOn,...
    $maxNumLoadSteps,...
    $airDensity,...
    $RayleighAlpha,...
    $RayleighBeta,...
    $elementOrder,...
    $delta_t,...
    $platformTurbineConnectionNodeNumber,...
    $tolerance,...
    $maxIterations,...
    $numpBC,...
    $numEl,...
    $Omega_j,...
    $OmegaDot_j,...
    $mel,...
    $moiel,...
    $xmel,...
    $Fexternal,...
    $Fdof,...
    $CN2H,...
    $rbData,...
    $jointTransform,...
    $conn,...
    $joint,...
    $pBC,...
    $x,...
    $y,...
    $z,...
    $ac,...
    $twist,...
    $rhoA,...
    $EIyy,...
    $EIzz,...
    $GJ,...
    $EA,...
    $rhoIyy,...
    $rhoIzz,...
    $rhoJ,...
    $zcm,...
    $ycm,...
    $a,...
    $EIyz,...
    $alpha1,...
    $alpha2,...
    $alpha3,...
    $alpha4,...
    $alpha5,...
    $alpha6,...
    $rhoIyz,...
    $b,...
    $a0,...
    $aeroCenterOffset,...
    $elLen,...
    $psi,...
    $theta,...
    $roll,...
    $displ_s,...
    $displdot_s,...
    $displddot_s,...
    $concLoadnodeNum,...
    $concLoaddof,...
    $concLoadval,...
    $concStiffnodeNum,...
    $concStiffdof,...
    $concStiffval,...
    $concMassnodeNum,...
    $concMassdof,...
    $concMassval,...
    $StiffGennodeNum,...
    $StiffGendof1,...
    $StiffGendof2,...
    $StiffGenval,...
    $MassGennodeNum,...
    $MassGendof1,...
    $MassGendof2,...
    $MassGenval,...
    $DampGennodeNum,...
    $DampGendof1,...
    $DampGendof2,...
    $DampGenval,...
    $K11,...
    $K12,...
    $K13,...
    $K14,...
    $K15,...
    $K16,...
    $K22,...
    $K23,...
    $K24,...
    $K25,...
    $K26,...
    $K33,...
    $K34,...
    $K35,...
    $K36,...
    $K44,...
    $K45,...
    $K46,...
    $K55,...
    $K56,...
    $K66,...
    $M11,...
    $M15,...
    $M16,...
    $M22,...
    $M24,...
    $M33,...
    $M34,...
    $M44,...
    $M55,...
    $M56,...
    $M66,...
    $S11,...
    $S12,...
    $S13,...
    $S15,...
    $S16,...
    $S22,...
    $S23,...
    $S25,...
    $S26,...
    $S33,...
    $S35,...
    $S36,...
    $S55,...
    $S56,...
    $S66,...
    $S14_1,...
    $S14_2,...
    $S24_1,...
    $S24_2,...
    $S34_1,...
    $S34_2,...
    $S45_1,...
    $S45_2,...
    $S46_1,...
    $S46_2,...
    $S44_1,...
    $S44_2,...
    $S44_3,...
    $C12,...
    $C13,...
    $C23,...
    $C24,...
    $C25,...
    $C26,...
    $C34,...
    $C35,...
    $C36,...
    $C14_1,...
    $C14_2,...
    $C45_1,...
    $C45_2,...
    $C46_1,...
    $C46_2);
    """

    # Reconstruct dispOut.elStrain

    elStrain = Array{ElStrain, 1}(undef, mesh.numEl)

    for jj = 1:mesh.numEl
        elStrain[jj] =  ElStrain(eps_xx_0[jj*4-3:jj*4], eps_xx_z[jj*4-3:jj*4], eps_xx_y[jj*4-3:jj*4], gam_xz_0[jj*4-3:jj*4], gam_xz_y[jj*4-3:jj*4], gam_xy_0[jj*4-3:jj*4], gam_xy_z[jj*4-3:jj*4])
    end

    dispOut = DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1)
    # dispOut.elStrain = elStrain
    # println("Breakout Struct Dyn Time Elapsed: $(time() - start)")
    return elStrain,dispOut,FReaction_j

end
