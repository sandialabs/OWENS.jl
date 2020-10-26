import DelimitedFiles
import LinearAlgebra
import FLOWMath
import VAWTAero
using MATLAB #if this is used, must run from src location
include("../src/meshing_utilities.jl")
include("../src/structs.jl")
include("../src/file_reading.jl")
include("../src/structural_utilities.jl")

include("transientExec.jl")

function owens(owensfile,analysisType;
    delta_t=0.01,
    numTS=100,
    tocp=[0.0,1.1],
    Omegaocp=[0.0,1.0],
    OmegaInit=0.0,
    OmegaGenStart=0.0,
    usingRotorSpeedFunction=false,
    nlOn=true,
    Omega=0.0,
    turbineStartup=0,
    spinUpOn=false,
    numModesToExtract=20,
    displInitGuess=0.0, #TODO: clean this up, is overwritten below
    airDensity=1.2041,
    aeroElasticOn = false,        # aeroElastic flags, and air density,
    aeroForceOn = true,
    guessFreq = 0,          #``guess"" modal frequency
    gravityOn = true,             #flag to activate gravity loading in structural dynamics/static simulations,
    generatorOn = false, #Initialize only, gets changed later on,
    omegaControl = false, #Initialize only, gets changed later on,
    totalNumDof = 0.0, #Initialize only, gets changed later on,
    iterationType = "NR", # nlParams
    adaptiveLoadSteppingFlag = true,
    tolerance = 1.0000e-06,
    maxIterations = 50,
    maxNumLoadSteps = 20,
    minLoadStepDelta = 0.0500,
    minLoadStep = 0.0500,
    prescribedLoadStep = 0.0)

    # if(occursin("S",analysisType)) #STATIC ANALYSIS
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



    if (occursin("TNB",analysisType) || occursin("TD",analysisType)) #TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration)
        if usingRotorSpeedFunction
            _,OmegaInit,_ = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
        else
            #this option uses a discretely specified rotor speed profile
            OmegaInit = Omegaocp[1] #TODO: simplify this downstream since the data doesn't need to be repeated
        end
    end

    # elseif(occursin("ROM",analysisType)) #REDUCED ORDER MODEL FOR TRANSIENT ANALYSIS
    #     model.delta_t = varargin{3} #time step size
    #     model.numTS = varargin{4}   #number of time steps
    #     model.numModesForROM = varargin{5} #number of lower system modes to include in ROM
    #     model.nlOn = varargin{6}    #flag for nonlinear elastic calculation
    #     turbineOpFlag = varargin{7}
    #     if(turbineOpFlag == 1) #generator start up operation mode
    #         model.OmegaInit = varargin{8} #initial rotor speed
    #     elseif(turbineOpFlag == 2) # self starting operation mode
    #         model.OmegaInit = varargin{8} #initial rotor speed (Hz)
    #         model.OmegaGenStart = varargin{9} #rotor speed at which generator activates (Hz)
    #     else                         #specified rotor speed profile
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
    # elseif(occursin("F",analysisType))  #MANUAL FLUTTER ANALYSIS
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
    # elseif(occursin("FA",analysisType)) #AUTOMATED FLUTTER ANALYSIS
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

    aeroLoadsOn      = Bool(real(parse(Int,line[1]))) #flag for activating aerodynamic analysis

    blddatafilename  = string(fdirectory, line[delimiter_idx[1][1]+1:delimiter_idx[2][1]-1]) #blade data file name
    aeroloadfile = string(fdirectory, line[delimiter_idx[2][1]+1:end]) #.csv file containing CACTUS aerodynamic loads
    owensfile = string(blddatafilename[1:end-4], ".owens") #TODO: this is redundant and confusing since it is specified at the beginning, clean up
    line             = readline(fid) #flag to include drive shaft effects
    driveShaftFlag   = real(parse(Int,line[1:2]))
    driveshaftfilename = string(fdirectory, line[3:end]) #drive shaft file name

    generatorfilename = string(fdirectory, readline(fid)) #generator file name
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
    bladeData,_,_,_ = readBladeData(blddatafilename) #reads overall blade data file
    BC = readBCdata(bcdatafilename,mesh.numNodes,numDofPerNode) #read boundary condition file
    el = readElementData(mesh.numEl,eldatafilename,ortdatafilename,bladeData) #read element data file (also reads orientation and blade data file associated with elements)
    joint = DelimitedFiles.readdlm(jntdatafilename,'\t',skipstart = 0) #readJointData(jntdatafilename) #read joint data file
    # rbarFileName = [owensfile(1:end-6),".rbar"] #setrbarfile
    # [model.joint] = readRBarFile(rbarFileName,model.joint,mesh) #read rbar file name
    nodalTerms = readNodalTerms(ndldatafilename) #read concentrated nodal terms file
    # [model] = readPlatformFile(model,platformFlag,platfilename)
    hydroOn = false
    platformTurbineConnectionNodeNumber = 1
    initCond = []
    JgearBox =0.0
    gearRatio = 1.0             #set gear ratio and efficiency to 1
    gearBoxEfficiency = 1.0
    useGeneratorFunction = false
    generatorProps = 0.0

    #     [model] = readDriveShaftProps(model,driveShaftFlag,driveshaftfilename) #reads drive shaft properties
    driveTrainOn = false          #set drive shaft unactive
    driveShaftProps = DriveShaftProps(0.0,0.0)       #set drive shat properties to 0

    if (occursin("TNB",analysisType)||occursin("TD",analysisType)||occursin("ROM",analysisType)) #for transient analysis...

        initCond = readInitCond(initcondfilename) #read initial conditions

        if !(occursin("[",generatorfilename)) #If there isn't a file
            useGeneratorFunction = true
            generatorProps = 0.0
        else
            useGeneratorFunction = false
            generatorProps = readGeneratorProps(generatorfilename) #reads generator properties
        end

    end

    outFilename = generateOutputFilename(owensfile,analysisType) #generates an output filename for analysis results #TODO: map to the output location instead of input

    jointTransform, reducedDOFList = createJointTransform(joint,mesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    BC.map = calculateBCMap(BC.numpBC,BC.pBC,numDofPerNode,reducedDOFList) #TODO: put this before BC so the structs can be made immutable in the future. #create boundary condition map from original DOF numbering to reduced/constrained DOF numbering

    numReducedDof = length(jointTransform[1,:])
    # BC.redVectorMap = zeros(numReducedDof,1)
    BC.redVectorMap = constructReducedDispVectorMap(mesh.numNodes,numDofPerNode,numReducedDof,BC) #TODO: put this before BC so the structs can be made immutable in the future. #create a map between reduced and full DOF lists

    nlParams = NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
    maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

    model = Model(analysisType,turbineStartup,usingRotorSpeedFunction,tocp,initCond,numTS,delta_t,Omegaocp,
    aeroElasticOn,aeroForceOn,aeroLoadsOn,driveTrainOn,airDensity,
    guessFreq,gravityOn,generatorOn,hydroOn,JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,totalNumDof,
    spinUpOn,nlOn,numModesToExtract,aeroloadfile,owensfile,outFilename,RayleighAlpha,
    RayleighBeta,elementOrder,joint,platformTurbineConnectionNodeNumber,jointTransform,
    reducedDOFList,bladeData,nlParams,BC,nodalTerms,driveShaftProps)

    #     if(occursin("S",analysisType)) #EXECUTE STATIC ANALYSIS
    #         [model.nlParams] = readNLParamsFile(owensfile)
    #         if(length(varargin)<=4 || ~model.nlOn)                #sets initial guess for nonlinear calculations
    #             displInitGuess = zeros(mesh.numNodes*6,1)
    #         end
    #
    #         OmegaStart = 0.0
    #         staticExec(model,mesh,el,displInitGuess,Omega,OmegaStart)
    #     end
    #
    if (occursin(analysisType,"M") || occursin(analysisType,"F")) #EXECUTE MODAL OR MANUAL FLUTTER ANALYSIS
        # [model.nlParams] = readNLParamsFile(owensfile) #TODO: clean this up, is redundant
        if (displInitGuess!=0 || !nlOn)
            displInitGuess = zeros(mesh.numNodes*6)
        end
        OmegaStart = 0.0
        # freq,damp=modalExec(model,mesh,el,displInitGuess,Omega,OmegaStart)
    end
    #
    #     if(occursin("FA",analysisType)) #EXECUTE AUTOMATED FLUTTER ANALYSIS
    #         displ = zeros(mesh.numNodes*6,1)
    #         OmegaStart = 0.0
    #         [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
    #     end
    #
    if (occursin("TNB",analysisType)||occursin("TD",analysisType)||occursin("ROM",analysisType)) #EXECUTE TRANSIENT ANALYSIS
        # [model.nlParams] = readNLParamsFile(owensfile) #TODO: this isn't really used, clean up
        # Juno.@enter transientExec(model,mesh,el)
        transientExec(model,mesh,el)
        freq = 0.0
        damp = 0.0
    end
    #
    # end
    #



    return freq,damp
end



function calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
    #calculateBCMap   calculates a boundary condition map
    #   [bcMap] = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
    #
    #   This function creates a boundary condition map between full and reduced
    #   dof listing as a result of constraints.
    #
    #      input:
    #      numpBC            = number of boundary conditions
    #      pBC               = array of boundary  condition data
    #      numDofPerNode     = number of degrees of freedom per node
    #      reducedDofList    = array of reduced DOF numbering
    #
    #      output:
    #      elStorage         = map for boundary conditions between full and
    #                          reduced dof list


    constrainedDof = zeros(numpBC)
    for i=1:numpBC
        constrainedDof[i] = (pBC[i,1]-1)*numDofPerNode + pBC[i,2]  #creates an array of constrained DOFs
    end
    constrainedDof = sort(constrainedDof)

    reducedDOFCount = length(reducedDofList)

    bcMap = zeros(reducedDOFCount)
    index = 1
    for i=1:reducedDOFCount
        if reducedDofList[i] in constrainedDof  #searches reduced DOF for constrained DOFs
            bcMap[i] = -1
        else
            bcMap[i] = index
            index = index + 1
        end
    end

    return bcMap

end


function generateOutputFilename(owensfilename,analysisType)
    #This function generates an output file name depending on the analysis type

    #find the last "." in owensfilename - helps to extract the prefix in the .owens
    index = findlast(".",owensfilename)[1]

    if (occursin("M",analysisType)||occursin("F",analysisType)||occursin("FA",analysisType)) #output filename (*.out) for modal/flutter analysis
        outputfilename = string(owensfilename[1:index-1],".out")
    elseif (occursin("S",analysisType)) #output file name (*_static.mat) for static analysis
        outputfilename = string(owensfilename[1:index-1],"_static.mat")
    elseif (occursin("TNB",analysisType)||occursin("TD",analysisType)||occursin("ROM",analysisType)) #output filename (*.mat) for transient analysis
        outputfilename = string(owensfilename[1:index-1],".mat")
    end

    return outputfilename

end

function calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained)
    #This function searches over all DOFs in a structural model and
    #determines and returns "dofVector" containing only unconstrained DOFs

    #loop over all DOFs in the model checking if constrained by BC or not
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                #             dofVector(index) = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    dofVector = zeros(Int,index)
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                dofVector[index] = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    return dofVector
end

function constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,BC)
    #This function creates a map of unconstrained DOFs between a full
    #listing and reduced listing (aftger constraints have been applied)

    bcdoflist=zeros(Int, BC.numpBC)

    #create a listing of constrained DOFs from boundary condition file
    for i=1:BC.numpBC
        bcnodenum = BC.pBC[i,1]
        bcdofnum = BC.pBC[i,2]
        bcdoflist[i] = (bcnodenum-1)*numDofPerNode + bcdofnum
    end

    dofList = calculateReducedDOFVector(numNodes,numDofPerNode,BC.isConstrained) #calculate a reduced (unconstrained) DOF vector

    redVectorMap = zeros(numReducedDof)

    for i=1:numReducedDof

        if (i in bcdoflist)              #creates a map of unconstrained reduced DOFs
            redVectorMap[i] = -1.0
        else
            index = findall(x->x==i,dofList)[1]
            redVectorMap[i] = index
        end

    end
    return redVectorMap
end
