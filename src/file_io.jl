
function readMesh(filename)
    #readMesh  reads mesh file and stores data in mesh object
    #   [mesh] = readMesh(filename)
    #
    #   This function reads the mesh file and stores data in the mesh object.
    #
    #      input:
    #      filename      = string containing mesh filename

    #      output:
    #      mesh          = object containing mesh data

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

    mesh = Mesh(nodeNum,
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


function readBCdata(bcfilename,numNodes,numDofPerNode)
    #readBDdata  reads boundary condition file
    #   [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
    #
    #   This function reads the boundray condition file and stores data in the
    #   boundary condition object.
    #
    #      input:
    #      bcfilename    = string containing boundary condition filename
    #      numNodes      = number of nodes in structural model
    #      numDofPerNode = number of degrees of freedom per node

    #      output:
    #      BC            = object containing boundary condition data

    fid = open(bcfilename)       #open boundary condition file
    numpBC = real(parse(Int,readline(fid))) #read in number of boundary conditions (displacement boundary conditions)
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

    BC = BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    [],
    [])

    return BC

end


function readBladeData(filename)
    #readBladeDAta reads blade data
    #   [bladeData] = readBladeData(filename)
    #
    #   This function reads blade data from file
    #
    #   input:
    #   filename      = string containing file name for blade data file
    #
    #   output:
    #   bladeData     = object containing blade data
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


function readElementData(numElements,elfile,ortfile,bladeData_struct)
    #readElementData  reads element data
    #   [el] = readElementData(numElements,elfile,ortfile,bldfile
    #
    #   This function reads element data and stores data in the element data
    #   object.
    #
    #      input:
    #      numElements   = number of elements in structural mesh
    #      elfile        = element data filename
    #      ortfile       = element orientation filename
    #      bldfile       = blade data filename

    #      output:
    #      el            = element data object

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

    sectionPropsArray = Array{SectionPropsArray, 1}(undef, numElements)

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

        sectionPropsArray[i] = SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

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
    el = El(sectionPropsArray,elLen,psi,theta,roll,rotationalEffects)

    return el

end

function readGeneratorProps(generatorfilename)
    #readGeneratorProps reads generator properties from file
    #   [genprops] = readGeneratorProps(generatorfilename)
    #
    #   This function reads generator properties from file.
    #
    #   input:
    #   generatorfilenanme  = string containing generator property file name
    #
    #   output:
    #   genprops          = model object containing generator properties


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

function readInitCond(filename)
    #readInitCond reads initial conditions
    #   [initCond] = readInitCond(filename)
    #
    #   This function reads initial conditions from file
    #
    #   input:
    #   filename      = string containing file name for initial conditions file
    #
    #   output:
    #   initCond      = array containing initial conditions

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

function writeOwensNDL(fileRoot, nodes, cmkType, cmkValues)
    # writeOwensNDL writes a nodal input file for the OWENS Toolkit
    #   This function writes a boundary condition file for the OWENS Toolkit
    #      input:
    #      fileRoot     = string containing input prefix of file name
    #      output:     (NONE)
    # **********************************************************************

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

function readNodalTerms(;filename="none",data=zeros(2))
    #readNodalTerms reads concentrated nodal terms file
    #   [nodalTerms] = readNodalTerms(filename)
    #
    #   This function reads the nodal terms file and stores data in the nodal
    #   terms object.
    #
    #      input:
    #      filename      = string containing nodal terms filename
    #
    #      output:
    #      nodalTerms    = object containing concentrated nodal data
    if filename!="none"
        data = DelimitedFiles.readdlm(filename,' ',skipstart = 0)
    end

    n_M = sum(count.(x->(x=='M'), data[:,2]))
    n_K = sum(count.(x->(x=='K'), data[:,2]))
    n_C = sum(count.(x->(x=='C'), data[:,2]))
    n_F = sum(count.(x->(x=='F'), data[:,2]))

    concMnodeNum = zeros(n_M)
    concMdof1 = zeros(n_M)
    concMdof2 = zeros(n_M)
    concMval = zeros(n_M)

    concKnodeNum = zeros(n_K)
    concKdof1 = zeros(n_K)
    concKdof2 = zeros(n_K)
    concKval = zeros(n_K)

    concCnodeNum = zeros(n_C)
    concCdof1 = zeros(n_C)
    concCdof2 = zeros(n_C)
    concCval = zeros(n_C)

    concFnodeNum = zeros(n_F)
    concFdof1 = zeros(n_F)
    concFdof2 = zeros(n_F)
    concFval = zeros(n_F)

    i_M = 1
    i_K = 1
    i_C = 1
    i_F = 1
    for i_data = 1:length(data[:,1])
        if data[i_data,2][1] == 'M'
            concMnodeNum[i_M] = data[i_data,1]
            concMdof1[i_M] = data[i_data,3]
            if length(data[1,:])==5 #If 6x6 general method
                concMdof2[i_M] = data[i_data,4]
            end
            concMval[i_M] = data[i_data,end]
            i_M += 1
        elseif data[i_data,2][1] == 'K'
            concKnodeNum[i_K] = data[i_data,1]
            concKdof1[i_K] = data[i_data,3]
            if length(data[1,:])==5 #If 6x6 general method
                concKdof2[i_K] = data[i_data,4]
            end
            concKval[i_K] = data[i_data,end]
            i_K += 1
        elseif data[i_data,2][1] == 'C'
            concCnodeNum[i_C] = data[i_data,1]
            concCdof1[i_C] = data[i_data,3]
            if length(data[1,:])==5 #If 6x6 general method
                concCdof2[i_C] = data[i_data,4]
            end
            concCval[i_C] = data[i_data,end]
            i_C += 1
        elseif data[i_data,2][1] == 'F'
            concFnodeNum[i_F] = data[i_data,1]
            concFdof1[i_F] = data[i_data,3]
            if length(data[1,:])==5 #If 6x6 general method
                concFdof2[i_F] = data[i_data,4]
            end
            concFval[i_F] = data[i_data,end]
            i_F += 1
        else
            error("Unknown Nodal Data Type")
        end
    end

    if length(data[1,:])==5
        @warn "General 6x6 concentrated diagonal terms are being applied to the old diagonal method with coulping (e.g. mass and force through acceleration), and no coupling is happening for the non-diagonal terms"
        #TODO: implement the 6x6 terms since they will be necessary for the linearized platform since there is strong cross coupling
        concLoad = Array{ConcNDL, 1}(undef, n_F)
        for i_F = 1:n_F
            if concFnodeNum[i_F] == concFdof2[i_F]
                concLoad[i_F]= ConcNDL(concFnodeNum[i_F], concFdof1[i_F], concFval[i_F])
            end
        end

        concStiff = Array{ConcNDL, 1}(undef, n_K)
        for i_K = 1:n_K
            if concKdof1[i_K] == concKdof2[i_K]
                concStiff[i_K] = ConcNDL(concKnodeNum[i_K], concKdof1[i_K], concKval[i_K])
            end
        end

        concMass = Array{ConcNDL, 1}(undef, n_M)
        for i_M = 1:n_M
            if concMdof1[i_M] == concMdof2[i_M]
                concMass[i_M] = ConcNDL(concMnodeNum[i_M], concMdof1[i_M], concMval[i_M])
            end
        end

        concStiffGen = Array{ConcNDLGen, 1}(undef, n_K)
        for i_K = 1:n_K
            concStiffGen[i_K] = ConcNDLGen(concKnodeNum[i_K],concKdof1[i_K],concKdof2[i_K],concKval[i_K])
        end

        concMassGen = Array{ConcNDLGen, 1}(undef, n_M)
        for i_M = 1:n_M
            concMassGen[i_M] = ConcNDLGen(concMnodeNum[i_M], concMdof1[i_M], concMdof2[i_M], concMval[i_M])
        end

        concDampGen = Array{ConcNDLGen, 1}(undef, n_C)
        for i_C = 1:n_C
            concDampGen[i_C] = ConcNDLGen(concCnodeNum[i_C], concCdof1[i_C], concCdof2[i_C], concCval[i_C])
        end

    elseif length(data[1,:])==4
        @warn "Only diagonal terms being used, there are no cross terms"
        # This portion is different in that it uses the nongeneral terms and applies them to the general just at the diagonal, TODO: once the general terms are implemented, this needs to be updated
        concLoad = Array{ConcNDL, 1}(undef, n_F)
        for i_F = 1:n_F
            concLoad[i_F]= ConcNDL(concFnodeNum[i_F], concFdof1[i_F], concFval[i_F])
        end

        concStiff = Array{ConcNDL, 1}(undef, n_K)
        for i_K = 1:n_K
            concStiff[i_K] = ConcNDL(concKnodeNum[i_K], concKdof1[i_K], concKval[i_K])
        end

        concMass = Array{ConcNDL, 1}(undef, n_M)
        for i_M = 1:n_M
            concMass[i_M] = ConcNDL(concMnodeNum[i_M], concMdof1[i_M], concMval[i_M])
        end

        concStiffGen = Array{ConcNDLGen, 1}(undef, n_K)
        for i_K = 1:n_K #NOTE dof1 is being used twice since in this case we didn't read in any cross terms!
            concStiffGen[i_K] = ConcNDLGen(concKnodeNum[i_K],concKdof1[i_K],concKdof1[i_K],concKval[i_K])
        end

        concMassGen = Array{ConcNDLGen, 1}(undef, n_M)
        for i_M = 1:n_M #NOTE dof1 is being used twice since in this case we didn't read in any cross terms!
            concMassGen[i_M] = ConcNDLGen(concMnodeNum[i_M], concMdof1[i_M], concMdof1[i_M], concMval[i_M])
        end

        concDampGen = Array{ConcNDLGen, 1}(undef, n_C)
        for i_C = 1:n_C #NOTE dof1 is being used twice since in this case we didn't read in any cross terms!
            concDampGen[i_C] = ConcNDLGen(concCnodeNum[i_C], concCdof1[i_C], concCdof1[i_C], concCval[i_C])
        end
    else
        error("Wrong number of terms in the .ndl file")
    end

    #store concentrated nodal term data in nodalTerms object
    return NodalTerms(concLoad,concStiff,concMass,concStiffGen,concMassGen,concDampGen)

end

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

function generateOutputFilename(owensfilename,analysisType)
    #This function generates an output file name depending on the analysis type

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

function ModalOutput(freq,damp,phase1,phase2,imagComponentSign,filename)
    #writeOutput writes output from modal analysis
    #   [freqSorted,dampSorted,imagCompSignSorted] = writeOutput(freq,damp,...
    #                                                phase1,phase2,...
    #                                                imagComponentSign,fid)
    #
    #   This function writes an output file for modal analysis.
    #
    #      input:
    #      freq               = array of modal frequencies
    #      damp               = array of modal damping ratios
    #      phase1             = array of in phase mode shapes
    #      phase2             = array of out of phase mode shapes
    #      imagComponentSign  = array of sign of imaginary components
    #      fid                = file identifier for output
    #
    #      output:
    #      freqSorted         = array of sorted(by frequency) modal frequencies
    #      dampSorted         = array of sorted(by frequency) modal damping ratios
    #      imagCompSignSorted = array of sorted(by frequency) of imaginarycomponentSign array
    if filename!="none"
        println(filename)
        fid=open(filename,"w")
    end
    Nnodes = size(phase1[:,:,1])[1] #gets number of nodes for mode shape printing

    freq,map,posIndex = bubbleSort(freq) #sorts frequency #TODO: get rid of this

    Nmodes = length(freq)

    dampSorted = zeros(Nmodes)
    freqSorted = zeros(Nmodes)
    imagCompSignSorted = zeros(Nmodes)
    U_x_0 = zeros(Nnodes,Nmodes)
    U_y_0 = zeros(Nnodes,Nmodes)
    U_z_0 = zeros(Nnodes,Nmodes)
    theta_x_0 = zeros(Nnodes,Nmodes)
    theta_y_0 = zeros(Nnodes,Nmodes)
    theta_z_0 = zeros(Nnodes,Nmodes)
    U_x_90 = zeros(Nnodes,Nmodes)
    U_y_90 = zeros(Nnodes,Nmodes)
    U_z_90 = zeros(Nnodes,Nmodes)
    theta_x_90 = zeros(Nnodes,Nmodes)
    theta_y_90 = zeros(Nnodes,Nmodes)
    theta_z_90 = zeros(Nnodes,Nmodes)

    index = 1
    for i=posIndex:1:posIndex+(Nmodes-1) #prints mode frequency, damping and in/out of phase mode shapes
        if filename!="none"
            Printf.@printf(fid,"MODE # %0.0f \n\n",index)
            Printf.@printf(fid,"Frequency: %e: \n",freq[i])
            Printf.@printf(fid,"Damping %e: \n",damp[map[i]])
            Printf.@printf(fid,"0 deg Mode Shape:\n")
            Printf.@printf(fid,"U_x          U_y          U_z          theta_x     theta_y     theta_z \n")
        end
        for j=1:Nnodes
            U_x_0[j,index] = phase1[j,1,map[i]]
            U_y_0[j,index] = phase1[j,2,map[i]]
            U_z_0[j,index] = phase1[j,3,map[i]]
            theta_x_0[j,index] = phase1[j,4,map[i]]
            theta_y_0[j,index] = phase1[j,5,map[i]]
            theta_z_0[j,index] = phase1[j,6,map[i]]
            if filename!="none"
                Printf.@printf(fid,"%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",U_x_0[j,index],U_y_0[j,index],U_z_0[j,index],theta_x_0[j,index],theta_y_0[j,index],theta_z_0[j,index])
                Printf.@printf(fid,"\n")
            end
        end

        if filename!="none"
            Printf.@printf(fid,"\n")

            Printf.@printf(fid,"90 deg Mode Shape:\n")
            Printf.@printf(fid,"U_x          U_y          U_z          theta_x     theta_y     theta_z \n")
        end

        for j=1:Nnodes
            U_x_90[j,index] = phase2[j,1,map[i]]
            U_y_90[j,index] = phase2[j,2,map[i]]
            U_z_90[j,index] = phase2[j,3,map[i]]
            theta_x_90[j,index] = phase2[j,4,map[i]]
            theta_y_90[j,index] = phase2[j,5,map[i]]
            theta_z_90[j,index] = phase2[j,6,map[i]]
            if filename!="none"
                Printf.@printf(fid,"%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",U_x_90[j,index],U_y_90[j,index],U_z_90[j,index],theta_x_90[j,index],theta_y_90[j,index],theta_z_90[j,index])
                Printf.@printf(fid,"\n")
            end
        end

        if (i<posIndex+(Nmodes-1)) && filename!="none"
            Printf.@printf(fid,"\n\n")
        end



        dampSorted[i] = damp[map[i]]
        freqSorted[i] = freq[i]
        imagCompSignSorted[i] = imagComponentSign[map[i]]

        index = index + 1
    end

    if filename!="none"
        close(fid)
    end

    return freqSorted,dampSorted,imagCompSignSorted,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
end

function bubbleSort(A)
    #bubbleSort Sorts the vector A in ascending order
    #   [A,origMap,posIndex] = bubbleSort(A)
    #
    #   This function accepts a vector A, sorts the vector in ascending order,
    #   outputting the sorted vector, a map to the original ordering, and the
    #   index at which a positive value first occurs.
    #
    #      input:
    #      A        = vector to be sorted
    #
    #      output:
    #      A        = sorted vector
    #      origMap  = map of sorted vector to original ordering
    #      posIndex = index at which positive value first occurs

    Aorig=deepcopy(A)
    len = length(A)
    swapped = true
    origMap = collect(1:len)
    while swapped
        swapped = false

        for i=1:len-1
            if (A[i+1] < A[i])
                temp = A[i]
                A[i] = A[i+1]
                A[i+1] = temp
                swapped = true

                temp2 = origMap[i]
                origMap[i] = origMap[i+1]
                origMap[i+1] = temp2
            end
        end
    end

    #     posIndex = length(A)/2+1
    #     [posIndex] = findPositiveCrossOver(A)
    posIndex = 1

    return A,origMap,posIndex

end

function readResultsModalOut(resultsFile,numNodes)
    data = DelimitedFiles.readdlm(resultsFile,'\t',skipstart = 0)

    nmodes = Int(length(data[:,1])/171) #This requires the specific formatting of the file to remain the same
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

    i_line = 1
    for i_mode = 1:nmodes
        i_line = (i_mode-1)*171+1

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

        i_line = i_line+84 #90 degree shapes, with the max value scaled to 1

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
readNuMadGeomCSV(Rhub, Rtip, B; precone=0.0, turbine=false,
mach=nothing, re=nothing, rotation=nothing, tip=PrandtlTipHub())

Parameters defining the rotor (apply to all sections).

**Arguments**
- `NuMad_geom_xlscsv_file::String`: name of the numad excel CSV file being read (!!! THE NUMAD TAB MUST BE SAVED AS A CSV FOR THIS TO WORK !!!)


**Returns**
- `Output::NuMad`: numad structure as defined in the NuMad structure docstrings.
"""

function readNuMadGeomCSV(NuMad_geom_xlscsv_file)
    #TODO: add composite orientation
    csvdata = DelimitedFiles.readdlm(NuMad_geom_xlscsv_file,',',skipstart = 0)

    n_station = length(csvdata[4:end,1])- sum(isempty.(csvdata[4:end,1]))
    n_web = Int(csvdata[1,6])
    n_stack = Int(csvdata[1,8])
    n_segments = Int(csvdata[2,8])
    span = Float64.(csvdata[4:n_station+3,2])
    #TODO: interpolations
    airfoil = csvdata[4:n_station+3,3]
    te_type = csvdata[4:n_station+3,4]
    twist_d = Float64.(csvdata[4:n_station+3,5])
    chord = Float64.(csvdata[4:n_station+3,6])
    xoffset = Float64.(csvdata[4:n_station+3,7])
    aerocenter = Float64.(csvdata[4:n_station+3,8])

    # Read stack info
    stack_idx_end = 10+n_stack-1
    stack_mat_types = Int.(csvdata[2,10:stack_idx_end])
    stack_layers_tmp = csvdata[4:n_station+3,10:stack_idx_end]
    if mod(stack_layers_tmp[1],1) != 0.0
        @warn "Stack layers not integer value, rounding"
    end
    stack_layers = round.(Int,stack_layers_tmp)

    seg_idx_end = stack_idx_end+n_segments+1
    segments = Float64.(csvdata[4:n_station+3,stack_idx_end+1:seg_idx_end])

    DP_idx_end = seg_idx_end+n_segments+1
    DPtypes = csvdata[4:n_station+3,seg_idx_end+1:DP_idx_end]

    skin_seq = Array{Seq, 2}(undef, n_station,n_segments) #can be any number of stack nums, so we have to make non-square containers

    skin_idx_end = DP_idx_end+n_segments

    for sta_idx = 1:n_station
        # sta_idx = 1
        for seg_idx = 1:n_segments
            # seg_idx = 1
            str = split(csvdata[3+sta_idx,DP_idx_end+seg_idx],",")
            if length(str)>1 && str[2]=="" #allow for single number
                str = str[1]
            end
            skin_seq[sta_idx,seg_idx] = Seq([parse(Int,x) for x in str])
        end
    end

    web_seq = Array{Seq, 2}(undef, n_station,n_web) #can be any number of stack nums, so we have to make non-square containers
    web_dp = Array{Seq, 2}(undef, n_station,n_web) #this is fixed size square, but it's easier to do it this way

    for web_idx = 1:n_web
        # web_idx = 1
        for sta_idx = 1:n_station
            # sta_idx = 1
            str = split(csvdata[3+sta_idx,skin_idx_end+web_idx*2-1],",")
            if !isempty(str[1])
                if str[2]=="" #allow for single number
                    str = str[1]
                end
                web_seq[sta_idx,web_idx] = Seq([parse(Int,x) for x in str])

                str = split(csvdata[3+sta_idx,skin_idx_end+web_idx*2],",")
                web_dp[sta_idx,web_idx] = Seq([parse(Int,x) for x in str])
            end
        end
    end

    return NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
end

function readNuMadMaterialsCSV(NuMad_geom_xlscsv_file)


    csvdata = DelimitedFiles.readdlm(NuMad_geom_xlscsv_file,',',skipstart = 0)

    data_start = 0
    data_end = length(csvdata[:,1])
    for ii = 1:length(csvdata[:,1])
        if csvdata[ii,1]=="Material ID"
            data_start = ii+1
        end

        if data_start!=0 && ii>data_start && !isa(csvdata[ii,1],Number) # in case they have a file with excess rows
            data_end = ii-1
            break
        end
    end

    names = csvdata[data_start:data_end,2]
    e1 = Float64.(csvdata[data_start:data_end,5]) .* 1e6
    e2 = Float64.(csvdata[data_start:data_end,6]) .* 1e6
    g12 = Float64.(csvdata[data_start:data_end,8]) .* 1e6
    anu = Float64.(csvdata[data_start:data_end,11])  #ratio
    rho = Float64.(csvdata[data_start:data_end,14])  #g/cc * 1000 #kg/m3
    xt = Float64.(csvdata[data_start:data_end,15]) .* 1e6  #pa
    xc = Float64.(csvdata[data_start:data_end,16]) .* 1e6  #pa
    yt = ones(length(e1)) .* 100.0e6 #made up
    yc = ones(length(e1)) .* 100.0e6  #made up
    s = ones(length(e1)) .* 100.0e6  #made up
    plythickness = Float64.(csvdata[data_start:data_end,4]) .* 1e-3 #meters

    return OWENS.plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness))

end
