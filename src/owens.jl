include("readMesh.jl")
include("readBladeData.jl")
include("readBCdata.jl")
include("readElementData.jl")
include("readJointData.jl")
include("readNodalTerms.jl")
mutable struct Model
    analysisType
    turbineStartup
    aeroElasticOn
    aeroForceOn
    airDensity
    guessFreq
    gravityOn
    generatorOn
    hydroOn
    OmegaGenStart
    omegaControl
    totalNumDof
    spinUpOn
    nlOn
    numModesToExtract
    aeroloadfile
    owensfile
    RayleighAlpha
    RayleighBeta
    elementOrder
    joint
    platformTurbineConnectionNodeNumber
    bladeData
    nlParams
    BC
    nodalTerms
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

function owens(owensfile,analysisType;Omega=0.0,spinUpOn=false,numModesToExtract=20,displInitGuess=0.0,airDensity=1.2041)

    #spinUpOn flag for pre-stressed modal analysis
    #lowest freq modes extracted first
    # displInitGuess sets initial guess for nonlinear calculations

    #owens Startup function for the OWENS toolkit
    # **********************************************************************
    # *                   Part of the SNL OWENS toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #   [freq,damp]=owens(varargin)
    #
    #   This function is a start up function for launching various analysis
    #   modes of the OWENS toolkit.
    #
    #      input:
    #      varargin      = input parameter list
    #         varargin{1} is the .owens file associated with analysis
    #         varargin{2} is a string describing analysis type
    #                     "S" = static analysis
    #                     "M" = modal analysis
    #                     "TNB" = transient analysis with Newmark-Beta time
    #                     integration
    #                     "ROM" = reduced order model transient analysis
    #      output:
    #      freq         = array of modal frequencies (when applicable to
    #                     analysis type)
    #      damp         = array of modal damping (when applicable to analysis
    #                     type)
    #      displ        = array containing converged solution for static
    #                     displacement

    turbineStartup = 0           #initialization of turbine startup,,
    aeroElasticOn = false        # aeroElastic flags, and air density,
    aeroForceOn = true
    airDensity = 0
    guessFreq = 0          #``guess"" modal frequency
    gravityOn = true             #flag to activate gravity loading in structural dynamics/static simulations,
    generatorOn = false #Initialize only, gets changed later on,
    OmegaGenStart = 0.0 #Initialize only, gets changed later on,
    omegaControl = false #Initialize only, gets changed later on,
    totalNumDof = 0.0 #Initialize only, gets changed later on,
    nlOn = true

    # nlParams
    iterationType = "NR"
    adaptiveLoadSteppingFlag = true
    tolerance = 1.0000e-06
    maxIterations = 50
    maxNumLoadSteps = 20
    minLoadStepDelta = 0.0500
    minLoadStep = 0.0500
    prescribedLoadStep = 0.0

    # if(strcmp(analysisType,"S")) #STATIC ANALYSIS
    #     Omega = varargin{3}            #initialization of rotor speed (Hz)
    #     model.nlOn= varargin{4}        #flag for nonlinear elastic calculation
    #     if(length(varargin)>4)                #sets initial guess for nonlinear calculations
    #         displInitGuess = varargin{5}
    #     end
    #     #    if(length(varargin)>5)                #sets air density if simple thin
    #     #        model.airDensity = varargin{6}   # airfoil theory loading desired
    #     #    else
    #     model.airDensity = 1.2041
    #     #    end



    # elseif(strcmp(analysisType,"TNB")||strcmp(analysisType,"TD")) #TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration)
    #     model.delta_t = varargin{3}  # time step size
    #     model.numTS = varargin{4}    # number of time steps
    #     model.nlOn = varargin{5}     # flag for nonlinear elastic calculation
    #     turbineOpFlag = varargin{6}           #turbine operation flag
    #     if(turbineOpFlag == 1) #generator start up operation mode
    #         model.turbineStartup = turbineOpFlag
    #         model.OmegaInit = varargin{7}   #initial rotor speed (Hz)
    #     elseif(turbineOpFlag == 2) #self starting operation mode
    #         model.turbineStartup = turbineOpFlag
    #         model.OmegaInit = varargin{7}   #initial rotor speed (Hz)
    #         model.OmegaGenStart = varargin{8} #rotor speed at which generator activates (Hz)
    #     else                       #specified rotor speed profile
    #         model.turbineStartup = 0
    #         if(length(varargin) == 6)
    #             model.usingRotorSpeedFunction = true #set flag to use user specified rotor speed function
    #             [~,model.OmegaInit,~] = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
    #         else
    #             #this option uses a discretely specified rotor speed profile
    #             model.usingRotorSpeedFunction = false #set flag to not use user specified rotor speed function
    #             model.tocp = varargin{7} #time points for rotor speed provfile
    #             Omegaocp = varargin{8} #rotor speed value at time points (Hz)
    #             model.Omegaocp = Omegaocp
    #             model.OmegaInit = Omegaocp(1)
    #         end
    #     end

    # elseif(strcmp(analysisType,"ROM")) #REDUCED ORDER MODEL FOR TRANSIENT ANALYSIS
    #     model.delta_t = varargin{3} #time step size
    #     model.numTS = varargin{4}   #number of time steps
    #     model.numModesForROM = varargin{5} #number of lower system modes to include in ROM
    #     model.nlOn = varargin{6}    #flag for nonlinear elastic calculation
    #     turbineOpFlag = varargin{7}
    #     if(turbineOpFlag == 1) #generator start up operation mode
    #         model.turbineStartup = turbineOpFlag
    #         model.OmegaInit = varargin{8} #initial rotor speed
    #     elseif(turbineOpFlag == 2) # self starting operation mode
    #         model.turbineStartup = turbineOpFlag
    #         model.OmegaInit = varargin{8} #initial rotor speed (Hz)
    #         model.OmegaGenStart = varargin{9} #rotor speed at which generator activates (Hz)
    #     else                         #specified rotor speed profile
    #         model.turbineStartup = 0
    #         if(length(varargin) == 7)
    #             model.usingRotorSpeedFunction = true #set flag to use user specified rotor speed function
    #             [~,model.OmegaInit,~] = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
    #         else
    #             #this option uses a discretely specified rotor speed profile
    #             model.usingRotorSpeedFunction = false #set flag to not use user specified rotor speed function
    #             model.tocp = varargin{8} #time points for rotor speed provfile
    #             Omegaocp = varargin{9} #rotor speed value at time points (Hz)
    #             model.Omegaocp = Omegaocp
    #             model.OmegaInit = Omegaocp(1)
    #         end
    #     end
    # elseif(strcmp(analysisType,"F"))  #MANUAL FLUTTER ANALYSIS
    #     Omega = varargin{3}   #rotor speed (Hz)
    #     model.spinUpOn = varargin{4} #flag for pre-stressed modal analysis
    #     model.guessFreq = varargin{5} #``guess"" modal frequency
    #     model.aeroElasticOn = true
    #     model.nlOn = true
    #
    #     if(length(varargin)>5)   #air density initialization
    #         model.airDensity = varargin{6}
    #     else
    #         model.airDensity = 1.2041
    #     end
    #     if(length(varargin)>6)   #number of lower system modes to extract
    #         model.numModesToExtract = varargin{7}
    #     else
    #         model.numModesToExtract = 20
    #     end
    #
    # elseif(strcmp(analysisType,"FA")) #AUTOMATED FLUTTER ANALYSIS
    #     omegaArray = varargin{3}    #array of rotor speed values(Hz)
    #     model.spinUpOn = varargin{4} #flag for pre-stressed modal analysis
    #     model.aeroElasticOn = true
    #     model.nlOn = true
    #
    #     if(length(varargin)>4)    #air density initializatio
    #         model.airDensity = varargin{5}
    #     else
    #         model.airDensity = 1.2041
    #     end
    #     if(length(varargin)>5)    #number of lower system modes to extract
    #         model.numModesToExtract = varargin{6}
    #     else
    #         model.numModesToExtract = 20
    #     end
    #
    # else
    #     error("Analysis type not recognized.")
    # end

    fid = open(owensfile,"r") #reads in model file names from .owens file

    last_delimiter   = findall(r"\W", owensfile) #Find the file directory
    fdirectory       = owensfile[1:last_delimiter[end-1][1]]

    meshfilename     = string(fdirectory, readline(fid)) #mesh file name
    eldatafilename   = string(fdirectory, readline(fid)) #element data file name
    ortdatafilename  = string(fdirectory, readline(fid)) #element orientation file name
    jntdatafilename  = string(fdirectory, readline(fid)) #joint data file name
    ndldatafilename  = string(fdirectory, readline(fid)) #concentrated nodal data file name
    bcdatafilename   = string(fdirectory, readline(fid)) #boundary condition file name
    line             = readline(fid)
    platformFlag     = real(parse(Int,line[1:2]))
    platfilename     = string(fdirectory, line[3:end])

    initcondfilename = string(fdirectory, readline(fid)) #initial condition filename

    line             = readline(fid)
    delimiter_idx    = findall(" ",line)

    aeroFlag         = real(parse(Int,line[1])) #flag for activating aerodynamic analysis
    blddatafilename  = string(fdirectory, line[delimiter_idx[1][1]+1:delimiter_idx[2][1]-1]) #blade data file name
    aeroloadfile = string(fdirectory, line[delimiter_idx[2][1]+1:end]) #.csv file containing CACTUS aerodynamic loads

    line             = readline(fid) #flag to include drive shaft effects
    driveShaftFlag   = real(parse(Int,line[1:2]))
    driveshaftfilename = string(fdirectory, line[3:end]) #drive shaft file name

    generatorfilename = string(fdirectory, readline(fid)) #generator file name
    # rayleighDamping   = getSplitLine((fid),"	") #read in alpha/beta for rayleigh damping
    line = readline(fid)
    rayleighDamping = split(line)

    if (isempty(rayleighDamping))
        RayleighAlpha = 0.0
        RayleighBeta = 0.0
    else
        RayleighAlpha = parse(Float64,rayleighDamping[1])
        RayleighBeta = parse(Float64,rayleighDamping[2])
    end

    close(fid) # close .owens file

    #model definitions
    elementOrder = 1 #linear element order
    #--------------------------------------------
    mesh = readMesh(meshfilename) #read mesh file
    numDofPerNode = 6
    bladeData = readBladeData(blddatafilename) #reads overall blade data file
    BC = readBCdata(bcdatafilename,mesh.numNodes,numDofPerNode) #read boundary condition file
    el = readElementData(mesh.numEl,eldatafilename,ortdatafilename,bladeData) #read element data file (also reads orientation and blade data file associated with elements)
    joint = readJointData(jntdatafilename) #read joint data file
    # rbarFileName = [owensfile(1:end-6),".rbar"] #setrbarfile
    # [model.joint] = readRBarFile(rbarFileName,model.joint,mesh) #read rbar file name
    nodalTerms = readNodalTerms(ndldatafilename) #read concentrated nodal terms file
    # [model] = readPlatformFile(model,platformFlag,platfilename)
    hydroOn = false
    platformTurbineConnectionNodeNumber = 1

    # if(strcmp(analysisType,"TNB")||strcmp(analysisType,"TD")||strcmp(analysisType,"ROM")) #for transient analysis...
    #
    #     [model.initCond] = readInitCond(initcondfilename) #read initial conditions
    #
    #     if(aeroFlag)
    #         model.aeroLoadsOn = true
    #     else
    #         model.aeroLoadsOn = false
    #     end
    #
    #     #     [model] = readDriveShaftProps(model,driveShaftFlag,driveshaftfilename) #reads drive shaft properties
    #     model.driveTrainOn = false          #set drive shaft unactive
    #
    #     model.driveShaftProps.k = 0.0       #set drive shat properties to 0
    #     model.driveShaftProps.c = 0.0
    #     model.JgearBox =0.0
    #
    #     model.gearRatio = 1.0             #set gear ratio and efficiency to 1
    #     model.gearBoxEfficiency = 1.0
    #
    #     if(real(str2double(generatorfilename))==1.0)
    #         model.useGeneratorFunction = true
    #         model.generatorProps = 0.0
    #     else
    #         model.useGeneratorFunction = false
    #         [model.generatorProps] = readGeneratorProps(generatorfilename) #reads generator properties
    #     end
    #
    # end
    # outFilename = generateOutputFilename(owensfile,analysisType) #generates an output filename for analysis results #TODO: map to the output location instead of input
    #
    # jnt_struct = createJointTransform(joint,mesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    # model.jointTransform = jnt_struct.jointTransform
    # model.reducedDOFList = jnt_struct.reducedDOF
    # [model.BC.map] = calculateBCMap(BC.numpBC,BC.pBC,numDofPerNode,jnt_struct.reducedDOF) #create boundary condition map from original DOF numbering to reduced/constrained DOF numbering
    # numReducedDof = length(jnt_struct.jointTransform(1,:))
    # model.BC.redVectorMap = zeros(numReducedDof,1)
    # [model.BC.redVectorMap] = constructReducedDispVectorMap(mesh.numNodes,numDofPerNode,numReducedDof,model.BC) #create a map between reduced and full DOF lists

    nlParams = NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
        maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

    model = Model(analysisType,turbineStartup,aeroElasticOn,aeroForceOn,airDensity,
        guessFreq,gravityOn,generatorOn,hydroOn,OmegaGenStart,omegaControl,totalNumDof,
        spinUpOn,nlOn,numModesToExtract,aeroloadfile,owensfile,RayleighAlpha,
        RayleighBeta,elementOrder,joint,platformTurbineConnectionNodeNumber,
        bladeData,nlParams,BC,nodalTerms)

#     if(strcmp(analysisType,"S")) #EXECUTE STATIC ANALYSIS
#         [model.nlParams] = readNLParamsFile(owensfile)
#         if(length(varargin)<=4 || ~model.nlOn)                #sets initial guess for nonlinear calculations
#             displInitGuess = zeros(mesh.numNodes*6,1)
#         end
#
#         OmegaStart = 0.0
#         staticExec(model,mesh,el,displInitGuess,Omega,OmegaStart)
#     end
#
#     if(strcmp(analysisType,"M") || strcmp(analysisType,"F")) #EXECUTE MODAL OR MANUAL FLUTTER ANALYSIS
#         [model.nlParams] = readNLParamsFile(owensfile)
#         if(length(varargin)<=5 || ~model.nlOn)
#             displInitGuess = zeros(mesh.numNodes*6,1)
#         end
#         OmegaStart = 0.0
#         [freq,damp]=modalExec(model,mesh,el,displInitGuess,Omega,OmegaStart)
#     end
#
#     if(strcmp(analysisType,"FA")) #EXECUTE AUTOMATED FLUTTER ANALYSIS
#         displ = zeros(mesh.numNodes*6,1)
#         OmegaStart = 0.0
#         [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
#     end
#
#     if(strcmp(analysisType,"TNB")||strcmp(analysisType,"TD")||strcmp(analysisType,"ROM")) #EXECUTE TRANSIENT ANALYSIS
#         [model.nlParams] = readNLParamsFile(owensfile)
#         model.analysisType = analysisType
#         transientExec(model,mesh,el)
#         freq = 0
#         damp = 0
#     end
#
# end
#
# function [outputfilename] = generateOutputFilename(owensfilename,analysisType)
#     #This function generates an output file name depending on the analysis type
#
#     #find the last "." in owensfilename - helps to extract the prefix in the .owens
#     index_all = find(owensfilename == ".")
#     index = index_all(end)
#
#     if(strcmp(analysisType,"M")||strcmp(analysisType,"F")||strcmp(analysisType,"FA")) #output filename (*.out) for modal/flutter analysis
#         outputfilename = [owensfilename(1:index-1),".out"]
#     elseif(strcmp(analysisType,"S")) #output file name (*_static.mat) for static analysis
#         outputfilename = [owensfilename(1:index-1),"_static.mat"]
#     elseif(strcmp(analysisType,"TNB")||strcmp(analysisType,"TD")||strcmp(analysisType,"ROM")) #output filename (*.mat) for transient analysis
#         outputfilename = [owensfilename(1:index-1),".mat"]
#     end
#
# end
#
# function [redVectorMap] = constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,BC)
#     #This function creates a map of unconstrained DOFs between a full
#     #listing and reduced listing (aftger constraints have been applied)
#
#     bcdoflist=zeros(BC.numpBC,1)
#
#     #create a listing of constrained DOFs from boundary condition file
#     for i=1:BC.numpBC
#         bcnodenum = BC.pBC(i,1)
#         bcdofnum = BC.pBC(i,2)
#         bcdoflist(i) = (bcnodenum-1)*numDofPerNode + bcdofnum
#     end
#
#     dofList = calculateReducedDOFVector(numNodes,numDofPerNode,BC.isConstrained) #calculate a reduced (unconstrained) DOF vector
#
#     redVectorMap = zeros(numReducedDof,1)
#
#     for i=1:numReducedDof
#
#         if(ismember(i,bcdoflist))              #creates a map of unconstrained reduced DOFs
#             redVectorMap(i) = -1.0
#         else
#             index = find(ismember(dofList,i))
#             redVectorMap(i) = index
#         end
#
#     end
#
# end
#
# function [dofVector] = calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained)
#     #This function searches over all DOFs in a structural model and
#     #determines and returns "dofVector" containing only unconstrained DOFs
#
#     #loop over all DOFs in the model checking if constrained by BC or not
#     index = 1
#     for i=1:numNodes
#         for j=1:numDofPerNode
#             if~(isConstrained((i-1)*numDofPerNode + j))
#                 #             dofVector(index) = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
#                 index = index + 1
#             end
#         end
#     end
#
#     dofVector = zeros(1,index)
#     index = 1
#     for i=1:numNodes
#         for j=1:numDofPerNode
#             if~(isConstrained((i-1)*numDofPerNode + j))
#                 dofVector(index) = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
#                 index = index + 1
#             end
#         end
#     end
#
#     return freq,damp
return 1,2
end
