function Modal(model,mesh,el,displ,Omega,OmegaStart)
    #Modal  Executive function for modal analysis
    #   [freq,damp]=Modal(model,mesh,el,displ,Omega,OmegaStart,fid)
    #
    #   This function executes modal analysis.
    #
    #   input:
    #   model          = object containing model information
    #   mesh           = object containing mesh information
    #   el             = object containing element information
    #   displ          = displacement vector for use in pre-stressed analysis
    #   Omega          = rotor speed (Hz)
    #   OmegaStart     = rotor speed (Hz) from previous analysis if stepping
    #                    through various rotor speeds, may be useful in load
    #                    stepping
    #
    #   output:
    #   freq           = modal frequencies of system
    #   damp           = modal damping ratios of system

    elStorage = initialElementCalculations(model,el,mesh) #performs initial element calculations

    # [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage) #calculate mass properties of structure

    #Do nonlinear iteration if needed
    if model.spinUpOn
        @warn "static analysis not yet coded"
        # model.aeroElasticOn = false
        # model.aeroForceOn = true
        # displ,_,staticAnalysisSuccessful=staticAnalysis(model,mesh,el,displ,Omega,OmegaStart,elStorage) #performs static analysis about specified operating condition
        # #save displnl displ #saves displacement vector from static analysis TODO: This doesn't appear to be used
    else
        staticAnalysisSuccessful = true
    end

    if staticAnalysisSuccessful
        freq,damp,phase1,phase2,imagCompSign= linearAnalysisModal(model,mesh,el,displ,Omega,elStorage) #performs modal analysis
    else
        error("Static analysis unsuccessful. Exiting")
    end
    return freq,damp,phase1,phase2,imagCompSign
end

function  linearAnalysisModal(model,mesh,el,displ,Omega,elStorage)
    #linearAnalysisModal performs modal analysis
    #   [freq,damp,phase1,phase2,imagCompSign] = linearAnalysisModal(model,mesh,
    #                                                 el,displ,Omega,elStorage)
    #
    #   This function performs a modal analysis of a structural dynamics
    #   system and returns freq, damping, and mode shapes
    #
    #   input:
    #   model          = object containing model information
    #   mesh           = object containing mesh information
    #   el             = object containing element information
    #   displ          = displacement vector for use in pre-stressed analysis
    #   Omega          = rotor speed (Hz)
    #   elStorage      = previously calculated element system matrices
    #
    #
    #   output:
    #   freq         = modal frequencies
    #   damp         = modal damping ratios
    #   phase1       = in phase mode shapes (real part of mode shape)
    #   phase2       = out of phase mode shapes (imaginary part of mode shape)
    #   imagCompSign = sign of imaginary component of eigenvalues
    # tic

    numEl = mesh.numEl  #extract mesh information from mesh object
    x = mesh.x
    y = mesh.y
    z = mesh.z
    conn = Int.(mesh.conn) #TODO fix this at the source
    numNodes = length(x)

    elementOrder = model.elementOrder  #extract element order from model

    BC = model.BC #extract boundary  condition object from model

    numNodesPerEl = elementOrder + 1  #do initialization
    numDOFPerNode = 6
    totalNumDOF = numNodes * numDOFPerNode

    Kg = zeros(totalNumDOF,totalNumDOF)
    Mg = zeros(totalNumDOF,totalNumDOF)
    Cg = zeros(totalNumDOF,totalNumDOF)

    nodalTerms = model.nodalTerms #extract concentrated nodal terms from model

    elInput = ElInput(elementOrder,
    true, #modalFlag,
    TimeInt(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
    zeros(2), #xloc,
    el.props[1], #sectionProps,
    0.0, #sweepAngle,
    0.0, #coneAngle,
    0.0, #rollAngle,
    0.0, #aeroSweepAngle,
    "none",#iterationType,
    model.nlOn, #useDisp,
    false, #preStress,
    model.aeroElasticOn, #aeroElasticOn,
    false, #aeroForceOn,
    0.0, #loadStepPrev,
    1.0, #loadStep,
    model.nlParams.maxNumLoadSteps,
    model.nlParams.maxIterations, #MAXIT
    model.nlParams.tolerance, #tolerance
    model.analysisType,
    zeros(numNodesPerEl*numDOFPerNode), #disp,
    zeros(numNodesPerEl*numDOFPerNode), #dispdot,
    zeros(numNodesPerEl*numDOFPerNode), #dispddot,
    zeros(numNodesPerEl*numDOFPerNode), #displ_iter,
    zeros(4,2), #concMass,
    zeros(6,2), #concStiff,
    zeros(6,2), #concLoad,
    zeros(numNodesPerEl*numDOFPerNode), #dispm1,
    zeros(numNodesPerEl), #x,
    zeros(numNodesPerEl), #y,
    zeros(numNodesPerEl), #z,
    model.gravityOn,
    model.RayleighAlpha,
    model.RayleighBeta,
    zeros(3), #accelVec,
    zeros(3), #omegaVec,
    zeros(3), #omegaDotVec,
    Omega,
    0.0,
    1.0*LinearAlgebra.I(3), #CN2H
    model.airDensity,
    0.0, #freq,
    true) #firstIteration

    for i=1:numEl   #element loop
        eldisp = zeros(numNodesPerEl*numDOFPerNode)
        #Calculate Ke and Fe for element i
        index = 1
        elInput.analysisType = "M"      #define element input object flags and element properties from el object
        elInput.iterationType = "NONE"
        elInput.displ_iter = eldisp
        elInput.useDisp = false
        elInput.preStress = true
        elInput.elementOrder = elementOrder
        elInput.modalFlag = true
        elInput.xloc = [0.0 el.elLen[i]]
        elInput.sectionProps = el.props[i]
        elInput.sweepAngle = el.psi[i]
        elInput.coneAngle = el.theta[i]
        elInput.rollAngle = el.roll[i]
        elInput.aeroSweepAngle = 0.0

        for j=1:numNodesPerEl       #define element coordinates and displacements associated with element
            elInput.x[j] = x[conn[i,j]]
            elInput.y[j] = y[conn[i,j]]
            elInput.z[j] = z[conn[i,j]]
            for k=1:numDOFPerNode
                eldisp[index] = displ[(conn[i,j]-1)*numDOFPerNode + k]
                index = index + 1
            end
        end

        #retrieve concentrated nodal terms associated with element
        massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff = ConcMassAssociatedWithElement(conn[i,:],model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad)

        #assign concentrated nodal terms and coordinates to element input
        #object
        elInput.concMass = massConc
        elInput.concStiff = stiffConc
        elInput.concLoad = loadConc
        elInput.disp = eldisp

        if el.rotationalEffects #TODO: vectorize this by i to make it element specific #activate or deactivate rotational effects for element
            elInput.Omega = Omega
            elInput.OmegaDot = 0.0
        else
            elInput.Omega = 0.0
            elInput.OmegaDot = 0.0
        end

        elInput.aeroElasticOn = model.aeroElasticOn   #set aeroelastic flag
        elInput.aeroForceOn = false
        elInput.airDensity = model.airDensity
        elInput.gravityOn = model.gravityOn
        elInput.RayleighAlpha = model.RayleighAlpha
        elInput.RayleighBeta = model.RayleighBeta

        if model.aeroElasticOn
            elInput.freq = model.guessFreq*2.0*pi     #set guess frequency if aeroelastic analysis
        else
            elInput.freq = 0.0 #Declare variable on all execution paths
        end

        elOutput = calculateTimoshenkoElementNL(elInput,elStorage[i]) #do element calculation

        Kg = assemblyMatrixOnly(elOutput.Ke,conn[i,:],numNodesPerEl,numDOFPerNode,Kg) #assemble element into global stiffness matrix
        Mg = assemblyMatrixOnly(elOutput.Me,conn[i,:],numNodesPerEl,numDOFPerNode,Mg) #assemble element into global mass matrix
        Cg = assemblyMatrixOnly(elOutput.Ce,conn[i,:],numNodesPerEl,numDOFPerNode,Cg) #assemble element into global damping matrix

    end

    #apply general 6x6  mass, damping, and stiffness matrices to nodes
    Kg,Mg,Cg = applyGeneralConcentratedTerms(Kg,Mg,Cg,model.nodalTerms.concStiffGen,model.nodalTerms.concMassGen,model.nodalTerms.concDampGen)

    #----------------------------------------------------------------------
    #APPLY CONSTRAINT
    Kg = applyConstraints(Kg,model.jointTransform)  #modify global matrices for joint constraints using joint transform
    Mg = applyConstraints(Mg,model.jointTransform)
    Cg = applyConstraints(Cg,model.jointTransform)

    #APPLY BOUNDARY CONDITIONS
    KgTotal = applyBCModal(Kg,BC.numpBC,BC.map)     #apply boundary conditions to global matrices
    MgTotal = applyBCModal(Mg,BC.numpBC,BC.map)
    CgTotal = applyBCModal(Cg,BC.numpBC,BC.map)

    if Omega==0.0 #set eigensolver flag
        solveFlag = 2
    else
        solveFlag = 1
    end
    # eigVec,eigVal = eigSolve(MgTotal,CgTotal,KgTotal)#,... #eigensolve of global system

    #constructs state space form (with mass matrix inverted)
    matwidth = length(MgTotal[:,1])
    sysMat = [zeros(matwidth,matwidth) 1.0*LinearAlgebra.I(matwidth)
    -MgTotal\KgTotal -MgTotal\CgTotal]
    # eigVal, eigVec = ArnoldiMethod.partialeigen(ArnoldiMethod.partialschur(sysMat))# nev=min(nx,nev), which=ArnoldiMethod.LM())[1])
    eigVal, eigVec = ArnoldiMethod.partialeigen(ArnoldiMethod.partialschur(sysMat; nev=size(sysMat)[1], which=ArnoldiMethod.LM())[1])
    # println(maximum(abs.(sysMat*eigVec .- eigVec*(eigVal.*LinearAlgebra.I(length(eigVal))))))
    #model.numModesToExtract,solveFlag)
    #save eigVectors eigVec #save eigenvector for later use (if needed) TODO: Doesn't appear to be used


    #extract frequency, damping, mode shapes from eigenvalues and vectors
    len3 = length(eigVal)
    freq = zeros(len3)
    damp = zeros(len3)
    len1 = Int(length(displ)/numDOFPerNode)
    phase1 = zeros(len1,numDOFPerNode,len3)
    phase2 = zeros(len1,numDOFPerNode,len3)
    sortedModes0 = zeros(len1,numDOFPerNode,len3)
    sortedModes = complex.(sortedModes0,zero(sortedModes0))
    imagCompSign = zeros(len3)
    for i=1:len3
        freq[i],damp[i],phase1[:,:,i],phase2[:,:,i],sortedModes[:,:,i] = extractFreqDamp(eigVal[i],eigVec[:,i],numDOFPerNode,model.jointTransform,model.reducedDOFList,model.BC,model.analysisType)
        imagCompSign[i] = sign(imag(eigVal[i]))
    end

    # #write output
    # t_modal = toc
    # disp('Elapsed time for modal analysis(s):')
    # disp(t_modal)

    if model.analysisType !="FA"
        freq,damp,imagCompSign = ModalOutput(freq,damp,phase1,phase2,imagCompSign,model.outFilename)
    end

    return freq,damp,phase1,phase2,imagCompSign

end

function applyBCModal(K,numpBC,bcMap)
    #applyBCModal Applies boundary conditions to system for modal analysis
    #   [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
    #
    #   This function applies boundary conditions to a system matrix for modal
    #   analysis
    #
    #      input:
    #      K             = assembled global system matrix
    #      BC            = struct of boundary condition information
    #      numDofPerNode = number of degrees of freedom per node

    #      output:
    #      K             = global system matrix with boundary conditions
    #      dofVector     = reduced DOF vector after imposing BCs

    numEq=length(bcMap)
    # indVec = zeros(numEq-numpBC,1) #can't pre-initialize...puts zeros in map
    # that causes problems when creating Knew

    index = 1
    indVec = zeros(Int,sum(bcMap .!= -1))
    for i=1:numEq
        if bcMap[i] != -1
            indVec[index] = i
            index = index +1
        end
    end

    if numpBC > 0
        Knew = K[indVec,indVec]
    else
        Knew = K
    end
    return Knew
end

function assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)
    #assemblyMatrixOnly Assembles element matrices into global sys of equations
    #   [Kg] = assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)
    #
    #   This function assembles the element matrix into the
    #   global system of equations
    #
    #      input:
    #      Ke             = element matrix
    #      conn           = element connectivity
    #      numNodesPerEl  = number of nodes per element
    #      numDofPerNode  = number of degrees of freedom per node
    #      Kg             = global system matrix

    #      output:
    #      Kg             = global system matrix with assembled element

    count = 1
    dofList = zeros(Int,numNodesPerEl*numDOFPerNode)
    for i=1:numNodesPerEl
        for j=1:numDOFPerNode
            dofList[count] = (conn[i]-1)*numDOFPerNode + j
            count = count +1
        end
    end

    Kg[dofList,dofList] = Kg[dofList,dofList] + Ke
    # numDOFPerEl = length(dofList)
    # #Assemble element i into global system
    #         for j=1:numDOFPerEl
    #             J = dofList(j)
    #             for m=1:numDOFPerEl
    #                 M = dofList(m)
    #                 Kg(J,M) = Kg(J,M) + Ke(j,m)
    #             end
    #         end

    return Kg

end

function extractFreqDamp(val,vec,numDOFPerNode,jointTransform,reducedDOFList,BC,analysisType)
    #extractFreqDamp   extract frequency, damping, mode shapes from eigsolution
    #   [freq,damp,phase1,phase2,sortedModes] = extractFreqDamp(val,vec,...
    #                                         numDOFPerNode,jointTransform,...
    #                                         reducedDOFList,BC,analysisType)
    #
    #   This function calculates the eigenvalues and vectors of a structural
    #   dynamic system.
    #
    #   input:
    #   val            = eigenvalue
    #   vec            = eigenvector
    #   numDOFPerNode  = number of degrees of freedom per node
    #   jointTransform = joint transformation matrix from reduced to full DOF
    #                    list
    #   reducedDOFList = listing of reduced DOFs
    #   BC             = boundary condition object containing boundary
    #                    condition info
    #   analysisType   = analysis type
    #
    #   output:
    #   freq        = modal frequency
    #   damp        = modal damping
    #   phase1      = in phase mode shape (real part of mode shape)
    #   phase2      = out of phase mode shape (imaginary part of mode shape)
    #   sortedModes = total, complex mode shape

    freq = abs(imag(val))/(2*pi)   #damped frequency for state space system
    damp = -real(val)/abs(imag(val)) #modal damping

    if (abs(imag(val)) < 1.0e-4)     #if imaginary part of eigenvalue is numeric zero treat as spring-mass system
        freq = sqrt(abs(real(val)))/(2*pi)
        damp = 0.0
    end

    # if (~strcmp(analysisType,'FA'))  #for all but automated flutter analysis
    #          [len,numModeShapes] = size(vec)
    dispReduc = constructReducedDispVecFromEigVec(vec,reducedDOFList,BC) #construct mode shape vector with boundary conditinos
    dispOrig = jointTransform*dispReduc #transform from reduced DOF list to full DOF list
    lenOrig=length(dispOrig)

    sortedModes0=zeros(Int(lenOrig/numDOFPerNode),numDOFPerNode)
    sortedModes = complex.(sortedModes0,0)
    for i=1:Int(lenOrig/numDOFPerNode)     #construct matrix of nodal DOF values from full DOF eigenvector
        for j=1:numDOFPerNode
            sortedModes[i,j] = dispOrig[(i-1)*6+j]
        end
    end

    phase1 = real(sortedModes)  #phase 1 is real part of modeshape (0 deg in phase)
    phase2 = imag(sortedModes)  #phase 2 is imag part of modeshape (90 deg out of phase)

    max1=maximum(abs.(phase1)) #find maximum values for modeshape normalization
    max2=maximum(abs.(phase2))


    if (max1>max2)
        maxOverall = max1
    else
        maxOverall = max2
    end

    if (maxOverall == 0)
        maxOverall = 1
    end

    phase1 = phase1./maxOverall  #normalize modeshapes
    phase2 = phase2./maxOverall

    if (abs(minimum(phase1))+1)<1.0e-4 || (abs(minimum(phase2))+1)<1.0e-4
        phase1 = -1*phase1
        phase2 = -1*phase2

    end

    # else  #return null mode shapes if mode shapes not requested
    #     phase1 = []
    #     phase2 = []
    #     sortedModes = []
    # end
    return freq,damp,phase1,phase2,sortedModes
end


function constructReducedDispVecFromEigVec(vec1,reducedDOFList,BC)
    #This function takes the original mode shape and modifies it to
    #account for boundary conditions
    bcdoflist=zeros(1,BC.numpBC)
    #form pBC DOF list
    for i=1:BC.numpBC
        bcnodenum = BC.pBC[i,1]
        bcdofnum = BC.pBC[i,2]
        bcdoflist[i] = (bcnodenum-1)*6 + bcdofnum
    end

    index = 1
    len = length(vec1)
    vec1Red0 = zeros(length(reducedDOFList))
    vec1Red = complex.(vec1Red0,0)
    for i=1:length(reducedDOFList)
        if isnothing(findfirst(x->x==reducedDOFList[i],bcdoflist))#(!ismember(reducedDOFList[i],bcdoflist))
            vec1Red[i] = vec1[Int(len/2+index),1]
            index = index + 1
        else
            vec1Red[i] = 0
        end
    end
    return vec1Red
end
