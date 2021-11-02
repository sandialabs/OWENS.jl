"""
    owens(owensfile,analysisType;
        delta_t=2e-3,
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
        displInitGuess=0.0,
        airDensity=1.2041,
        aeroElasticOn = false,
        guessFreq = 0,
        gravityOn = true,
        generatorOn = false,
        omegaControl = false,
        iterationType = "NR", # nlParams
        adaptiveLoadSteppingFlag = true,
        tolerance = 1.0000e-06,
        maxIterations = 50,
        maxNumLoadSteps = 20,
        minLoadStepDelta = 0.0500,
        minLoadStep = 0.0500,
        prescribedLoadStep = 0.0,
        elementOrder = 1,
        numDofPerNode = 6,
        hydroOn = false,
        platformTurbineConnectionNodeNumber = 1,
        JgearBox =0.0,
        gearRatio = 1.0,
        gearBoxEfficiency = 1.0,
        useGeneratorFunction = false,
        generatorProps = 0.0,
        driveTrainOn = false)

Original serial and file reading method of running an analysis.

#Inputs
See ?OWENS.Model and ?GyricFEA.FEAModel

#Outputs
See ?OWENS.Unsteady, ?GyricFEA.Modal

"""
function owens(owensfile,analysisType;
    delta_t=2e-3,
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
    displInitGuess=0.0,
    airDensity=1.2041,
    aeroElasticOn = false,        # aeroElastic flags, and air density,
    guessFreq = 0,          #``guess"" modal frequency
    gravityOn = true,             #flag to activate gravity loading in structural dynamics/static simulations,
    generatorOn = false,
    omegaControl = false,
    iterationType = "NR", # nlParams
    adaptiveLoadSteppingFlag = true,
    tolerance = 1.0000e-06,
    maxIterations = 50,
    maxNumLoadSteps = 20,
    minLoadStepDelta = 0.0500,
    minLoadStep = 0.0500,
    prescribedLoadStep = 0.0,
    elementOrder = 1, #linear element order
    numDofPerNode = 6,
    hydroOn = false,
    interpOrder = 2,
    platformTurbineConnectionNodeNumber = 1,
    JgearBox =0.0,
    gearRatio = 1.0,             #set gear ratio and efficiency to 1
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    driveTrainOn = false,          #set drive shaft unactive
    hydrodynLib = "none",
    moordynLib = "none")

    # if(analysisType=="S") #STATIC ANALYSIS
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

    if (analysisType=="TNB" || analysisType=="TD") #TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration)
        if usingRotorSpeedFunction
            _,OmegaInit,_ = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
        else
            #this option uses a discretely specified rotor speed profile
            OmegaInit = Omegaocp[1] #TODO: simplify this downstream since the data doesn't need to be repeated
        end
    end

    # elseif(analysisType=="ROM") #REDUCED ORDER MODEL FOR TRANSIENT ANALYSIS
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
    # elseif(analysisType=="F")  #MANUAL FLUTTER ANALYSIS
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
    # elseif(analysisType=="FA") #AUTOMATED FLUTTER ANALYSIS
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
    if analysisType!="M"
        owensfile = string(blddatafilename[1:end-4], ".owens") #TODO: this is redundant and confusing since it is specified at the beginning, clean up
    end
    line             = readline(fid) #flag to include drive shaft effects
    driveShaftFlag   = real(parse(Int,line[1:2]))
    driveshaftfilename = string(fdirectory, line[3:end]) #drive shaft file name

    generatorfilename = string(fdirectory, readline(fid)) #generator file name
    
    line = readline(fid)
    delimiter_idx = findall(" ",line)
    hydroOn = Bool(real(parse(Int,line[1]))) #flag for activating hydrodynamic analysis
    potflowfile  = string(fdirectory, line[delimiter_idx[1][1]+1:delimiter_idx[2][1]-1]) # potential flow file prefix
    interpOrder = real(parse(Int,line[delimiter_idx[2][1]+1])) # interpolation order for HD/MD libraries
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
    #--------------------------------------------
    mesh = readMesh(meshfilename) #read mesh file
    bladeData,_,_,_ = readBladeData(blddatafilename) #reads overall blade data file
    BC = readBCdata(bcdatafilename,mesh.numNodes,numDofPerNode) #read boundary condition file
    el = readElementData(mesh.numEl,eldatafilename,ortdatafilename,bladeData) #read element data file (also reads orientation and blade data file associated with elements)
    joint = DelimitedFiles.readdlm(jntdatafilename,'\t',skipstart = 0) #readJointData(jntdatafilename) #read joint data file
    nodalTerms = GyricFEA.readNodalTerms(filename=ndldatafilename) #read concentrated nodal terms file
    initCond = []

    #     [model] = readDriveShaftProps(model,driveShaftFlag,driveshaftfilename) #reads drive shaft properties

    driveShaftProps = DriveShaftProps(0.0,0.0)       #set drive shat properties to 0

    if (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #for transient analysis...

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
    jointTransform, reducedDOFList = GyricFEA.createJointTransform(joint,mesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    numReducedDof = length(jointTransform[1,:])

    nlParams = GyricFEA.NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
    maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

    model = Model(;analysisType,turbineStartup,usingRotorSpeedFunction,tocp,numTS,delta_t,Omegaocp,
    aeroLoadsOn,driveTrainOn,generatorOn,hydroOn,interpOrder,JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,
    numModesToExtract,aeroloadfile,owensfile,potflowfile,outFilename,bladeData,driveShaftProps)

    feamodel = GyricFEA.FEAModel(;analysisType,
    initCond,
    aeroElasticOn,
    guessFreq,
    airDensity,
    gravityOn,
    nlOn,
    spinUpOn,
    outFilename,
    RayleighAlpha,
    RayleighBeta,
    elementOrder,
    joint,
    platformTurbineConnectionNodeNumber,
    jointTransform,
    reducedDOFList,
    nlParams,
    numNodes = mesh.numNodes,
    numModes = numModesToExtract,
    pBC=BC.pBC,
    nodalTerms)

    #     if(analysisType=="S") #EXECUTE STATIC ANALYSIS
    #         if(length(varargin)<=4 || ~model.nlOn)                #sets initial guess for nonlinear calculations
    #             displInitGuess = zeros(mesh.numNodes*6,1)
    #         end
    #
    #         OmegaStart = 0.0
    #         staticExec(model,mesh,el,displInitGuess,Omega,OmegaStart)
    #     end
    #
    if (analysisType == "M" || analysisType == "F") #EXECUTE MODAL OR MANUAL FLUTTER ANALYSIS
        if (displInitGuess==0 || !nlOn)
            displInitGuess = zeros(mesh.numNodes*6)
        end
        OmegaStart = 0.0
        freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=Modal(feamodel,mesh,el;displ=displInitGuess,Omega,OmegaStart)
        return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
    end
    #
    #     if(analysisType=="FA") #EXECUTE AUTOMATED FLUTTER ANALYSIS
    #         displ = zeros(mesh.numNodes*6,1)
    #         OmegaStart = 0.0
    #         [freq,damp]=ModalAuto(model,mesh,el,displ,omegaArray,OmegaStart)
    #     end
    #

    if ((hydrodynLib!="none")||(moordynLib!="none")) # Map shared library paths for external dependencies
        bin = Bin(hydrodynLib, moordynLib)
    end

    if (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #EXECUTE TRANSIENT ANALYSIS
        aeroLoadsFile_root = model.aeroloadfile[1:end-16] #cut off the _ElementData.csv
        OWENSfile_root = model.owensfile[1:end-6] #cut off the .owens

        geomFn = string(aeroLoadsFile_root, ".geom")
        loadsFn = string(aeroLoadsFile_root, "_ElementData.csv")
        bldFn = string(OWENSfile_root, ".bld")
        elFn = string(OWENSfile_root, ".el")
        ortFn = string(OWENSfile_root, ".ort")
        meshFn = string(OWENSfile_root, ".mesh")

        aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)

        aeroForces(t) = externalForcing(t+delta_t,aerotimeArray,aeroForceValHist,aeroForceDof)

        deformAero(omega) = 0.0 #placeholder function
        Unsteady(model,feamodel,mesh,el,bin,aeroForces,deformAero)

        return model
    end

    #
    # end
    #

end

"""

    readMesh(filename)

Reads the mesh file and stores data in the mesh object.

input:
* `filename::string` string containing mesh path/to/filename.mesh

output:
* `mesh::GyricFEA.Mesh` see GyricFEA.Mesh
"""
function readMesh(filename)

    fid = open(filename,"r")   #open mesh file

    # temp = fscanf(fid,'#i',2)   #read in number of nodes and number of elements
    line = readline(fid)
    temp = split(line)

    numNodes = parse(Int,temp[1])
    numEl = parse(Int,temp[2])

    nodeNum = zeros(numNodes,1)
    x = zeros(numNodes,1)
    y = zeros(numNodes,1)
    z = zeros(numNodes,1)

    conn = zeros(numEl,2)
    elNum = zeros(numEl,1)

    for i=1:numNodes            # read in node number and node coordinates
        line = readline(fid)
        temp = split(line)
        nodeNum[i] = parse(Float64,temp[1])
        x[i] = parse(Float64,temp[2])
        y[i] = parse(Float64,temp[3])
        z[i] = parse(Float64,temp[4])
    end

    for i=1:numEl               # read in element number and connectivity list
        line = readline(fid)
        temp = split(line)
        elNum[i] = parse(Float64,temp[1])

        conn[i,1] = parse(Float64,temp[3])
        conn[i,2] = parse(Float64,temp[4])
    end

    line = readline(fid) #get blank line
    line = readline(fid)
    temp = split(line)
    numComponents = parse(Int,temp[1])
    meshSeg = zeros(Int,numComponents)
    for i=1:numComponents
        meshSeg[i] = parse(Int,temp[i+1])
    end

    close(fid)  #close mesh file

    mesh = GyricFEA.Mesh(nodeNum,
    numEl,
    numNodes,
    x,
    y,
    z,
    elNum,
    conn,
    zeros(Int,numEl),
    meshSeg,
    0,
    0,
    0)

    return mesh

end

"""

    readBCdata(bcfilename,numNodes,numDofPerNode)

This function reads the boundray condition file and stores data in the
boundary condition object.

#Input
* `bcfilename::string`:    string containing boundary condition filename
* `numNodes::int`:      number of nodes in structural model
* `numDofPerNode::int`: number of degrees of freedom per node

#Output
* `BC::GyricFEA.BC_struct`:   see GyricFEA.BC_struct, object containing boundary condition data
"""
function readBCdata(bcfilename,numNodes,numDofPerNode)

    fid = open(bcfilename)       #open boundary condition file
    numpBC = parse(Int,readline(fid)) #read in number of boundary conditions (displacement boundary conditions)
    pBC = zeros(Int,numpBC,3)         #initialize boundary conditions
    for i=1:numpBC

        line = readline(fid)

        # Find where all of the delimiters are
        #first two are boundary condition node number and local DOF number
        #third is boundary condition value (typically zero)
        delimiter_idx = [0;collect.(Int,findall(" ",line));length(line)+1]
        # Extract the data from the beginning to the last delimiter
        for k = 2:length(delimiter_idx)
            pBC[i,k-1] = Int(parse(Float64,line[delimiter_idx[k-1][1]+1:delimiter_idx[k][1]-1]))
        end

    end

    totalNumDof = numNodes*numDofPerNode

    numsBC = 0
    nummBC = 0

    close(fid)

    #create a vector denoting constrained DOFs in the model (0 unconstrained, 1
    #constrained)


    #calculate constrained dof vector
    isConstrained = zeros(totalNumDof,1)
    constDof = (pBC[:,1].-1)*numDofPerNode + pBC[:,2]
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if ((i-1)*numDofPerNode + j in constDof)
                isConstrained[index] = 1
            end
            index = index + 1
        end
    end

    BC = GyricFEA.BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    [],
    [])

    return BC

end

"""

    readBladeData(filename)

This function reads blade data from file

#Input
* `filename::string`:   string containing /path/to/bladedata.bld

#Output
* `bladeData::BladeData`:  see ?BladeData object containing blade data
"""
function readBladeData(filename)

    a = DelimitedFiles.readdlm(filename,'\t',skipstart = 0)

    bladeNum = a[:,1]

    numBlades = maximum(bladeNum)
    #     numStruts = min(bladeNum)
    #     if (numStruts>0)
    #         numStruts = 0
    #     else
    #         numStruts = abs(numStruts)
    #     end

    strutStartIndex = 0
    for i=1:length(bladeNum)
        if (isempty(a[i,end]))
            strutStartIndex = i
            break
        end
    end



    if (strutStartIndex!=0)
        #         strutDataBlock = a(strutStartIndex:end,:)
        #         [strutEntries, _] = size(strutDataBlock)
        #         numNodesPerStrut = strutEntries/numStruts
        #         numElPerStrut = numNodesPerStrut - 1
    else
        temp=size(a)

        strutStartIndex = temp[1] + 1
    end

    bladeDataBlock = a[1:strutStartIndex-1,:]
    bladeEntries, _ = size(bladeDataBlock)
    numNodesPerBlade = round(Int,bladeEntries/numBlades)

    structuralSpanLocNorm = zeros(numBlades,numNodesPerBlade)
    structuralNodeNumbers = zeros(numBlades,numNodesPerBlade)
    structuralElNumbers = zeros(numBlades,numNodesPerBlade)
    for i=1:numBlades
        structuralSpanLocNorm[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,2]./bladeDataBlock[i*numNodesPerBlade,2]
        structuralNodeNumbers[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,3]
        structuralElNumbers[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,4]
    end

    bladeData = BladeData(numBlades,  #assign data to bladeData object
    bladeDataBlock[:,1],
    bladeDataBlock[:,2],
    bladeDataBlock[:,3],
    bladeDataBlock[:,4],
    bladeDataBlock[:,5:end])

    return bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers

end

"""

    readElementData(numElements,elfile,ortfile,bldfile

Reads element data and stores data in the element data
object.

#Input
* `numElements::int`:  number of elements in structural mesh
* `elfile::string`:       element data path/to/filename
* `ortfile::string`:      element orientation path/to/filename
* `bldfile::string`:      blade data path/to/filename

#Output
* `el::GyricFEA.El`:       see GyricFEA.El    element data object
"""
function readElementData(numElements,elfile,ortfile,bladeData_struct)

    fid = open(elfile,"r") #open element data file

    ac = zeros(2)
    twist = zeros(2)
    rhoA = zeros(2)
    EIyy = zeros(2)
    EIzz = zeros(2)
    GJ = zeros(2)
    EA = zeros(2)
    rhoIyy = zeros(2)
    rhoIzz = zeros(2)
    rhoJ = zeros(2)
    zcm = zeros(2)
    ycm = zeros(2)
    a = zeros(2)
    EIyz = zeros(2)
    alpha1 = zeros(2)
    alpha2 = zeros(2)
    alpha3 = zeros(2)
    alpha4 = zeros(2)
    alpha5 = zeros(2)
    alpha6 = zeros(2)
    rhoIyz = zeros(2)
    b = zeros(2)
    a0 = zeros(2)
    aeroCenterOffset = zeros(2)

    sectionPropsArray = Array{GyricFEA.SectionPropsArray, 1}(undef, numElements)

    data1 = zeros(1,17)
    data2 = zeros(1,17)
    for i=1:numElements
        data1=parse.(Float64,split(readline(fid))) #read element data
        data2=parse.(Float64,split(readline(fid)))

        #structural properties
        ac = -([data1[2], data2[2]].-0.5) #TODO: why are we doing it this way???
        twist=[data1[3], data2[3]]
        rhoA = [data1[4], data2[4]]
        EIyy = [data1[5], data2[5]]
        EIzz = [data1[6], data2[6]]
        if (minimum(abs.(EIyy - EIzz)) < 1.0e-3)
            EIzz = EIzz.*1.0001
        end
        GJ = [data1[7], data2[7]]
        EA = [data1[8], data2[8]]
        alpha1 = [data1[9], data2[9]]

        rhoIyy = [data1[10], data2[10]]
        rhoIzz = [data1[11], data2[11]]
        rhoJ = [(data1[10]+data1[11]), (data2[10]+data2[11])]
        zcm = [data1[14], data2[14]]
        ycm = [data1[15], data2[15]]
        a = [data1[17], data2[17]]

        #coupling factors
        EIyz = [0.0, 0.0]
        alpha1 = [0.0, 0.0]
        alpha2 = [0.0, 0.0]
        alpha3 = [0.0, 0.0]
        alpha4 = [0.0, 0.0]
        alpha5 = [0.0, 0.0]
        alpha6 = [0.0, 0.0]
        rhoIyz = [0.0, 0.0]
        b = [0.0, 0.0]
        a0 = [2*pi, 2*pi]

        sectionPropsArray[i] = GyricFEA.SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

    end
    close(fid) #close element file

    nodeNum = bladeData_struct.nodeNum  #node number associated with blade section
    elNum = bladeData_struct.elementNum    #element number associated with blade section
    bladeData = bladeData_struct.remaining  #blade data

    chord = zeros(maximum(nodeNum),1)
    for i=1:length(elNum)
        chord[nodeNum[i]] = bladeData[i,10]  #store chord of blade sections
    end

    for i=1:length(elNum)
        if (elNum[i]!=-1)

            sectionPropsArray[elNum[i]].b = 0.5.*[chord[nodeNum[i]], chord[nodeNum[i+1]]] #element semi chord
            sectionPropsArray[elNum[i]].a0 = [bladeData[i,12], bladeData[i+1,12]]         #element lift curve slope (needed for flutter analysis)

            #convert "a" to semichord fraction aft of halfchord
            sectionPropsArray[elNum[i]].a = (sectionPropsArray[elNum[i]].a + 0.25*(2*sectionPropsArray[elNum[i]].b) - sectionPropsArray[elNum[i]].b)./sectionPropsArray[elNum[i]].b

            #convert "ac" to semichord fraction foreward of halfchord TODO: why are we doing it this way???
            sectionPropsArray[elNum[i]].ac = (sectionPropsArray[elNum[i]].ac).*2

            #physical aero center offset from elastic axis
            sectionPropsArray[elNum[i]].aeroCenterOffset = (sectionPropsArray[elNum[i]].ac).*sectionPropsArray[elNum[i]].b - sectionPropsArray[elNum[i]].a
        end
    end


    println("EIyz, rhoIyz deactivated")

    #read element orientation data
    elLen = zeros(numElements)
    psi = zeros(numElements)
    theta = zeros(numElements)
    roll = zeros(numElements)
    fid = open(ortfile,"r")
    for i=1:numElements
        temp = parse.(Float64,split(readline(fid)))
        elLen[i]=temp[5]
        psi[i]=temp[2]
        theta[i]=temp[3]
        roll[i]=temp[4]
    end
    close(fid) #close ort file

    rotationalEffects = ones(numElements)

    #store data in element object
    el = GyricFEA.El(sectionPropsArray,elLen,psi,theta,roll,rotationalEffects)

    return el

end

"""

    readGeneratorProps(generatorfilename)

This function reads generator properties from file.

#Input
* `generatorfilenanme::string`:  = string containing path/to/generatorfile

#Output
* `genprops`:    = model object containing generator properties
"""
function readGeneratorProps(generatorfilename)

    # fid = fopen(generatorfilename) #open generator property file
    # if (fid!=-1) #if file can be opened
    #         genprops.ratedTorque = fscanf(fid,'#f',1) #store rated torque
    #         dum = fgetl(fid)
    #         genprops.zeroTorqueGenSpeed = fscanf(fid,'#f',1) #store zero torque generator zpeed
    #         dum = fgetl(fid)
    #         genprops.pulloutRatio = fscanf(fid,'#f',1) #store pullout ratio
    #         dum = fgetl(fid)
    #         genprops.ratedGenSlipPerc= fscanf(fid,'#f',1) #store rated generator slip percentage
    #         dum = fgetl(fid)
    #
    #         fclose(fid) #close generator propery file
    genprops = 0.0
    println("GENERATOR NOT FULLY ENABLED")

    genprops = 0.0 #if generator property file does not exist, set object to null

    return genprops
end

"""

    readInitCond(filename)

Reads initial conditions from file

#Input
* `filename`      string containing file name for initial conditions file

#Output
* `initCond`      array containing initial conditions
"""
function readInitCond(filename)

    initCond =[] #initialize intial condition to null

    # fid = open(filename) #open initial  conditions file

    #         index = 1
    #         while(~feof(fid))
    #             temp1 = fscanf(fid,'#i',2) #read node number and local DOF number for initial cond.
    #             temp2 = fscanf(fid,'#f',1) #read value for initial cond.
    #
    #             #place node number, dof number and value into array
    #             initCond(index,1:3) = [temp1(1), temp1(2), temp2(1)]
    #
    #             index = index + 1
    #         end

    println("INITIAL CONDITIONS NOT FULLY ENABLED")
    return initCond
end

"""

    writeOwensNDL(fileRoot, nodes, cmkType, cmkValues)

writes a nodal input file

#Intput
* `fileRoot::string`: string path to desired location with name but no extension
* `nodes::int`: node numbers for C/M/K
* `cmkType::string`: "C" "M" or "K"
* `cmkValues::float`: C/M/K value

#Output
* `none`:

"""
function writeOwensNDL(fileRoot, nodes, cmkType, cmkValues)

    # open the BC file to save the boundary conditions to
    BCfile = string(fileRoot, ".ndl")    #construct file name

    open(BCfile, "w") do file
        # write out the boundary conditions into the file
        for nn = 1:length(nodes)
            # [row, col, val] = find(cmkValues[nn])
            indices = findall(x->x!=0,cmkValues[:,:,nn])
            # println(indices)
            for ii = 1:length(indices)
                row = indices[ii][1]
                col = indices[ii][2]
                write(file, "$(nodes[nn]) $(cmkType[nn]) $(row) $(col) $(cmkValues[row,col,nn])\n")
            end
        end
    end
end

"""
Internal, reads cactus .geom file and stores each column in an array within the CactusGeom struct
"""
function readCactusGeom(geom_fn)

    data = DelimitedFiles.readdlm(geom_fn,skipstart = 0)

    NBlade = Int(data[1,2])
    NStrut = Int(data[2,2])
    RotN = Float64.(data[3,2:4])
    RotP = Float64.(data[4,2:4])
    RefAR = Float64(data[5,2])
    RefR = Float64(data[6,2])

    blade = Array{Blade, 1}(undef, NBlade)

    idx = 9
    for bld_idx = 1:NBlade

        NElem = Int(data[idx+0,2])
        FlipN = Int(data[idx+1,2])
        QCx = Float64.(data[idx+2,2:NElem+2])
        QCy = Float64.(data[idx+3,2:NElem+2])
        QCz = Float64.(data[idx+4,2:NElem+2])
        tx = Float64.(data[idx+5,2:NElem+2])
        ty = Float64.(data[idx+6,2:NElem+2])
        tz = Float64.(data[idx+7,2:NElem+2])
        CtoR = Float64.(data[idx+8,2:NElem+2])
        PEx = Float64.(data[idx+9,2:NElem+1])
        PEy = Float64.(data[idx+10,2:NElem+1])
        PEz = Float64.(data[idx+11,2:NElem+1])
        tEx = Float64.(data[idx+12,2:NElem+1])
        tEy = Float64.(data[idx+13,2:NElem+1])
        tEz = Float64.(data[idx+14,2:NElem+1])
        nEx = Float64.(data[idx+15,2:NElem+1])
        nEy = Float64.(data[idx+16,2:NElem+1])
        nEz = Float64.(data[idx+17,2:NElem+1])
        sEx = Float64.(data[idx+18,2:NElem+1])
        sEy = Float64.(data[idx+19,2:NElem+1])
        sEz = Float64.(data[idx+20,2:NElem+1])
        ECtoR = Float64.(data[idx+21,2:NElem+1])
        EAreaR = Float64.(data[idx+22,2:NElem+1])
        iSect = Float64.(data[idx+23,2:NElem+1])
        idx += 25
        blade[bld_idx] = Blade(NElem,FlipN,QCx,QCy,QCz,tx,ty,tz,CtoR,PEx,PEy,PEz,tEx,tEy,tEz,nEx,nEy,nEz,sEx,sEy,sEz,ECtoR,EAreaR,iSect)
    end

    strut = Array{Strut, 1}(undef, NStrut)

    for strut_idx = 1:NStrut

        NElem = Int(data[idx+0,2])
        TtoC = Float64(data[idx+1,2])
        MCx = Float64.(data[idx+2,2:NElem+2])
        MCy = Float64.(data[idx+3,2:NElem+2])
        MCz = Float64.(data[idx+4,2:NElem+2])
        CtoR = Float64.(data[idx+5,2:NElem+2])
        PEx = Float64.(data[idx+6,2:NElem+1])
        PEy = Float64.(data[idx+7,2:NElem+1])
        PEz = Float64.(data[idx+8,2:NElem+1])
        sEx = Float64.(data[idx+9,2:NElem+1])
        sEy = Float64.(data[idx+10,2:NElem+1])
        sEz = Float64.(data[idx+11,2:NElem+1])
        ECtoR = Float64.(data[idx+12,2:NElem+1])
        EAreaR = Float64.(data[idx+13,2:NElem+1])
        BIndS = Int(data[idx+14,2])
        EIndS = Int(data[idx+15,2])
        BIndE = Int(data[idx+16,2])
        EIndE = Int(data[idx+17,2])

        idx += 19
        strut[strut_idx] = Strut(NElem,TtoC,MCx,MCy,MCz,CtoR,PEx,PEy,PEz,sEx,sEy,sEz,ECtoR,EAreaR,BIndS,EIndS,BIndE,EIndE)
    end

    return CactusGeom(NBlade,NStrut,RotN,RotP,RefAR,RefR,blade,strut)
end

"""
Internal, takes cactus loads and geometry and OWENS mesh and maps the loads to the blades' FEA dofs
"""
function mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)

    cactusGeom = readCactusGeom(geomFn)
    blade = cactusGeom.blade

    # aero_data = importCactusFile(loadsFn,1,2002,22,',')
    aero_data = DelimitedFiles.readdlm(loadsFn,',',skipstart = 1)
    #define these from params file
    ft2m = 1 / 3.281
    rho = 1.225
    #     RefAR = cactusGeom.RefAR*ft2m*ft2m
    RefR = cactusGeom.RefR*ft2m
    V = 25.0#6.787243728 #m/s #27.148993200000003
    @warn "Velocity is hardcoded here at $V"

    normTime = aero_data[:,1]

    numAeroEl = 0
    for i=1:cactusGeom.NBlade
        numAeroEl = numAeroEl + cactusGeom.blade[i].NElem
    end

    len,_ = size(aero_data)

    numAeroTS = Int(len/numAeroEl)
    # time = zeros(len/numAeroEl,1)
    time = normTime[1:Int(numAeroEl):end,1].*RefR[1]./V[1]

    urel = aero_data[:,15]
    uloc = urel.*V

    cn = aero_data[:,24]
    ct = aero_data[:,25]
    cm25 = aero_data[:,22]

    NperSpan = zeros(len)
    TperSpan = zeros(len)
    M25perSpan = zeros(len)
    Mecc = zeros(len)

    for i=1:len

        NperSpan[i] =  cn[i]  * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
        TperSpan[i] =  ct[i]  * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
        M25perSpan[i] = cm25[i] * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)*blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR
        momentArm = 0.0
        Mecc[i] = NperSpan[i] * momentArm
    end

    # Initialize bladeForces
    N = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)
    T = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)
    M25 = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)

    index = 1
    for i=1:numAeroTS
        for j=1:cactusGeom.NBlade
            for k=1:blade[j].NElem
                N[j,i,k] = NperSpan[index]
                T[j,i,k] = TperSpan[index]
                M25[j,i,k] = M25perSpan[index]
                index = index + 1
            end
        end
    end

    spanLocNorm = zeros(cactusGeom.NBlade,blade[1].NElem)
    for i=1:cactusGeom.NBlade
        spanLocNorm[i,:] = blade[i].PEy[1:blade[1].NElem[1,1],1].*RefR[1,1]/(blade[i].QCy[blade[1].NElem[1,1]+1,1]*RefR[1,1])
    end

    bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers = readBladeData(bldFn)

    #Initialize structuralLoad
    struct_N = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_T = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_M25 = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))

    # TODO: the interpolation input here is blade height vs blade span length, should be updated to be consistent
    for i=1:cactusGeom.NBlade
        for j=1:numAeroTS
            struct_N[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],N[i,j,:],structuralSpanLocNorm[i,:])
            struct_T[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],T[i,j,:],structuralSpanLocNorm[i,:])
            struct_M25[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],M25[i,j,:],structuralSpanLocNorm[i,:])
        end
    end

    _,numNodesPerBlade = size(structuralNodeNumbers)

    #integrate over elements

    #read element aero_data in
    mesh = readMesh(meshFn)
    el = readElementData(mesh.numEl,elFn,ortFn,bladeData)
    numDofPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6),numAeroTS)
    for i=1:numAeroTS
        for j = 1:cactusGeom.NBlade
            for k = 1:numNodesPerBlade-1
                #get element aero_data
                # orientation angle,xloc,sectionProps,element order]
                elNum = Int(structuralElNumbers[j,k])
                #get dof map
                node1 = Int(structuralNodeNumbers[j,k])
                node2 = Int(structuralNodeNumbers[j,k+1])
                dofList = [(node1-1)*numDofPerNode.+(1:6) (node2-1)*numDofPerNode.+(1:6)]

                elementOrder = 1
                x = [mesh.x[node1], mesh.x[node2]]
                elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
                xloc = [0 elLength]
                twist = el.props[elNum].twist
                sweepAngle = el.psi[elNum]
                coneAngle = el.theta[elNum]
                rollAngle = el.roll[elNum]

                extDistF2Node =  [struct_T[j,i,k]    struct_T[j,i,k+1]]
                extDistF3Node = -[struct_N[j,i,k]    struct_N[j,i,k+1]]
                extDistF4Node = -[struct_M25[j,i,k]  struct_M25[j,i,k+1]]

                Fe = GyricFEA.calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)

                #asssembly
                for m = 1:length(dofList)
                    Fg[dofList[m],i] =  Fg[dofList[m],i]+Fe[m]
                end

            end
        end
    end

    #reduce Fg to nonzero components
    #assumes any loaded DOF will never be identically zero throughout time
    #history
    # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
    # ForceDof = zeros(sum(Fg[:,1].!=0),1)
    ForceValHist = zeros(length(Fg[:,1]),length(Fg[1,:]))
    ForceDof = zeros(length(Fg[:,1]),1)
    index = 1
    for i=1:Int(maximum(maximum(structuralNodeNumbers))*6)
        # if !isempty(findall(x->x!=0,Fg[i,:]))

            ForceValHist[index,:] = Fg[i,:]
            ForceDof[index] = i
            index = index + 1
        # end
    end
    return time,ForceValHist,ForceDof,cactusGeom
end

"""
Internal, generates an output file name depending on the analysis type
"""
function generateOutputFilename(owensfilename,analysisType)

    #find the last "." in owensfilename - helps to extract the prefix in the .owens
    index = findlast(".",owensfilename)[1]

    if (analysisType=="M"||analysisType=="F"||analysisType=="FA") #output filename (*.out) for modal/flutter analysis
        outputfilename = string(owensfilename[1:index-1],".out")
    elseif (analysisType=="S") #output file name (*_static.mat) for static analysis
        outputfilename = string(owensfilename[1:index-1],"_static.mat")
    elseif (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #output filename (*.mat) for transient analysis
        outputfilename = string(owensfilename[1:index-1],".mat")
    end

    return outputfilename

end

"""
Internal, reads modal file and returns freq, damp, and modeshapes
"""
function readResultsModalOut(resultsFile,numNodes)
    data = DelimitedFiles.readdlm(resultsFile,'\t',skipstart = 0)

    nmodes = round(Int,min(length(data[:,1])/(numNodes*2+6),30))

    freq = zeros(nmodes)
    damp = zeros(nmodes)
    U_x_0 = zeros(numNodes,nmodes)
    U_y_0 = zeros(numNodes,nmodes)
    U_z_0 = zeros(numNodes,nmodes)
    theta_x_0 = zeros(numNodes,nmodes)
    theta_y_0 = zeros(numNodes,nmodes)
    theta_z_0 = zeros(numNodes,nmodes)
    U_x_90 = zeros(numNodes,nmodes)
    U_y_90 = zeros(numNodes,nmodes)
    U_z_90 = zeros(numNodes,nmodes)
    theta_x_90 = zeros(numNodes,nmodes)
    theta_y_90 = zeros(numNodes,nmodes)
    theta_z_90 = zeros(numNodes,nmodes)

    for i_mode = 1:nmodes
        i_line = (i_mode-1)*(numNodes*2+7)+1

        freq[i_mode] = parse(Float64,(split(data[i_line+1,1])[2])[1:end-1])
        damp[i_mode] = parse(Float64,(split(data[i_line+2,1])[2])[1:end-1])

        # 0 degree shapes, with the max value scaled to 1

        temp = Float64.(data[i_line+5:i_line+4+numNodes,1])
        U_x_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,2])
        U_y_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,3])
        U_z_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,4])
        theta_x_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,5])
        theta_y_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,6])
        theta_z_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())

        i_line = i_line+numNodes+2 #90 degree shapes, with the max value scaled to 1

        temp = Float64.(data[i_line+5:i_line+4+numNodes,1])
        U_x_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,2])
        U_y_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,3])
        U_z_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,4])
        theta_x_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,5])
        theta_y_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,6])
        theta_z_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
    end
    return freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
end

"""
    simpleGenerator(generatorProps,genSpeed)

Caclulates generator torque for simple induction generator

#Input
* `generatorProps` object containing generator properties, see ?model
* `genSpeed::float`       generator speed (Hz)

#Output
* `genTorque::float`      generator torque
"""
function simpleGenerator(model,genSpeed)

    #assign generator properties form model object
    ratedTorque        = model.ratedTorque;
    ratedGenSlipPerc   = model.ratedGenSlipPerc;
    zeroTorqueGenSpeed = model.zeroTorqueGenSpeed;
    pulloutRatio       = model.pulloutRatio;

    #calculate rated generator speed
    ratedGenSpeed = zeroTorqueGenSpeed*(1.0 + 0.01*ratedGenSlipPerc);

    #calculate slope between lower and upper torque limits for simple induction generator
    midSlope = (ratedTorque/(ratedGenSpeed-zeroTorqueGenSpeed));

    #calculate lower and upper torque limits of generator
    upperTorqueLimit =ratedTorque*pulloutRatio;
    lowerTorqueLimit = -upperTorqueLimit;

    #calculate upper and lower generator speeds at which linear torque vs. speed region begins/ends
    upperGenSpeed = zeroTorqueGenSpeed + upperTorqueLimit/midSlope;
    lowerGenSpeed = zeroTorqueGenSpeed - upperTorqueLimit/midSlope;

    #calculate generator torque
    if genSpeed<0.0#lowerGenSpeed
        genTorque = 0.0#lowerTorqueLimit;
    elseif genSpeed>upperGenSpeed
        genTorque = upperTorqueLimit;
    else
        genTorque = midSlope*(genSpeed-zeroTorqueGenSpeed);
    end

    return genTorque

end

#
# function stress_wrapper(usedmaterials,plyprops,spanlocy,orientation,strain,shear)
#
#     c_stress = []
#     stress = [] #used in vtk output
#     stresslayer = zeros(length(usedmaterials),length(spanlocy),length(strain[1,:]))
#     j=1
#     for i = 1:length(usedmaterials)
#
#         if !contains(usedmaterials[i],"foam") #don't check stress on foam
#             matnames = plyprops.names
#             idx = find(matnames -> matnames == usedmaterials[i],matnames)
#
#             material = plyprops.plies[idx[1]]
#             orien = orientation[j]
#             j+=1 #so we don't have to have extra design variables in orientation
#
#             stressi = stresscalc(strain,shear,material,orien)
#
#             if i==1 #only calc once since it's the same for the whole structure (assuming thin layups)
#                 stress = sqrt.(stressi[:,:,1].^2+stressi[:,:,2].^2+stressi[:,:,3].^2) #TODO break up or calculate von-mises etc
#             end
#
#             # determine contraint values
#             c_stressi = stresscon(stressi,material)
#             stresslayer[i,:,:] = c_stressi
#
#             c_stress = [c_stress;c_stressi]
#         end
#     end
#
#     return c_stress, stress, stresslayer
# end #stress_wrapper
#
# """
#     stresscalc(strain::Array{Float64,2},shear::Array{Float64,1},mat::Composites.material,theta_d::Float64)
#
# #Inputs
# * `strain::Array{Float64,2}`:: strain at every span location and every chord location around the airfoil
# * `shear::Array{Float64,1}`:: shear stress, needs to be included, is zeros for now
# * `mat::Composites.material`:: material properties
# * `theta_d::Float64`:: orientation of the ply in degrees
#
# #Outputs
# * `plystress::Array{Float64,2}`: stress at every span location and every chord location around the airfoil
# """
# function stresscalc(strain,shear,mat,theta_d)
#     # Get rotated material stiffness matrix
#     qbar = Composites.getQ(mat,theta_d)
#     # Determine Stresses from Strains in individual plys
#     nnodes = size(strain,1)
#     nloc = size(strain,2)
#     plystress = zeros(nnodes,nloc,3)
#     for i = 1:nnodes #Run entire length of wing #
#         for iloc = 1:nloc # Run at every specified strain location
#             # Convert shear stress to shear strain
#             gam12 = shear[i]/qbar[3,3]-qbar[1,3]/qbar[3,3]*strain[i,iloc]
#             # Calculate laminate stresses
#             stress = qbar*[strain[i,iloc],0,gam12]
#             # Transform stresses to ply axes
#             plystress[i,iloc,:] = Composites.rotstress(stress,theta_d)
#         end
#     end
#     return plystress
# end #stress_calc
