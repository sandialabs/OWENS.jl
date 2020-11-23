
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

    close(fid)  #close mesh file

    mesh = Mesh(nodeNum,
    numEl,
    numNodes,
    x,
    y,
    z,
    elNum,
    conn,
    zeros(Int,numEl))

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
    #     if(numStruts>0)
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

    bladeData = BladeData(numBlades,  #assign data to bladeData object #TODO: Should not be loading this file in multiple times
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
        ac = -([data1[2], data2[2]].-0.5)
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


    nodeNum = bladeData_struct.nodeNum  #node number associated with blade section
    elNum = bladeData_struct.elementNum    #element number associated with blade sectino
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

            #convert "ac" to semichord fraction foreward of halfchord
            sectionPropsArray[elNum[i]].ac = (sectionPropsArray[elNum[i]].ac).*2

            #physical aero center offset from elastic axis
            sectionPropsArray[elNum[i]].aeroCenterOffset = (sectionPropsArray[elNum[i]].ac).*sectionPropsArray[elNum[i]].b - sectionPropsArray[elNum[i]].a
        end
    end


    println("EIyz, rhoIyz deactivated")
    close(fid) #close element file

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

    #store data in element object
    el = El(sectionPropsArray,elLen,psi,theta,roll,true)

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
    # if(fid!=-1) #if file can be opened
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


function readNodalTerms(filename)
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

    concStiff=[] #initialize stiffness, load, and mass arrays to null
    concLoad=[]
    concMass=[]
    concStiffGen = []
    concMassGen = []
    concDampGen = []
    fid=open(filename,"r")  #open nodal terms file
    #     index = 1

    #         while(~feof(fid))
    #
    #             temp = strsplit(fgetl(fid))
    #
    #             nodeNum(index,1) = real(str2double(temp{1}))   #read in node number
    #             termType{index} = temp{2}  #read in concentrated term type (M=mass, K = stiffness, F = force)
    #             if(length(temp) == 4)
    #                 localDof{index} = real(str2double(temp{3}))
    #                 val(index) = real(str2double(temp{4})
    #             else
    #                 localDof{index} = [real(str2double(temp{3})) real(str2double(temp{4}))]
    #                 val(index) = real(str2double(temp{5}))
    #             end
    #             #         val(index,1) = fscanf(fid,'#f',1)       #read in value for concentrated nodal term
    #             index = index + 1
    #         end
    #
    #         fIndex = 1    #counters for concentrated nodal force, stiffness, and mass
    #         kIndex = 1
    #         mIndex = 1
    #         mIndex2 = 1
    #         kIndex2 = 1
    #         cIndex2 = 1
    #         for i=1:length(nodeNum)
    #             if(strcmpi(termType{i},'F'))    #read in concentrated load data
    #                 concLoad(fIndex).nodeNum = nodeNum(i)
    #                 concLoad(fIndex).dof = localDof{i}
    #                 concLoad(fIndex).val = val(i)
    #                 fIndex = fIndex + 1
    #             end
    #
    #             if(strcmpi(termType{i},'K'))    #read in concentrated stiffness data
    #                 concStiff(kIndex).nodeNum = nodeNum(i)
    #                 concStiff(kIndex).dof = localDof{i}
    #                 concStiff(kIndex).val = val(i)
    #                 kIndex = kIndex + 1
    #             end
    #
    #             if(strcmpi(termType{i},'M'))   #read in concentrated mass data
    #                 concMass(mIndex).nodeNum = nodeNum(i)
    #                 concMass(mIndex).dof = localDof{i}
    #                 concMass(mIndex).val = val(i)
    #                 mIndex = mIndex + 1
    #             end
    #
    #             if(strcmpi(termType{i},'M6'))   #read in concentrated mass data
    #                 concMassGen(mIndex2).nodeNum = nodeNum(i)
    #                 temp = localDof{i}
    #                 concMassGen(mIndex2).dof1 = temp(1)
    #                 concMassGen(mIndex2).dof2 = temp(2)
    #                 concMassGen(mIndex2).val = val(i)
    #                 mIndex2 = mIndex2 + 1
    #             end
    #
    #             if(strcmpi(termType{i},'K6'))   #read in concentrated stiffness data
    #                 concStiffGen(kIndex2).nodeNum = nodeNum(i)
    #                 temp = localDof{i}
    #                 concStiffGen(kIndex2).dof1 = temp(1)
    #                 concStiffGen(kIndex2).dof2 = temp(2)
    #                 concStiffGen(kIndex2).val = val(i)
    #                 kIndex2 = kIndex2 + 1
    #             end
    #
    #             if(strcmpi(termType{i},'C6'))   #read in concentrated damp data
    #                 concDampGen(cIndex2).nodeNum = nodeNum(i)
    #                 temp = localDof{i}
    #                 concDampGen(cIndex2).dof1 = temp(1)
    #                 concDampGen(cIndex2).dof2 = temp(2)
    #                 concDampGen(cIndex2).val = val(i)
    #                 cIndex2 = cIndex2 + 1
    #             end
    #         end
    println("NODAL FILE NOT LOADED, CODE UNFINISHED IN readNodalTerms") #TODO


    close(fid)

    #store concentrated nodal term data in nodalTerms object
    nodalTerms = NodalTerms(concLoad,concStiff,concMass,concStiffGen,concMassGen,concDampGen)

    return nodalTerms

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
