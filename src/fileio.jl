
"""
readNuMadGeomCSV(NuMad_geom_xlscsv_file)

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
    stack_layers_tmp = Float64.(csvdata[4:n_station+3,10:stack_idx_end])
    # if mod(stack_layers_tmp[1],1) != 0.0
    #     @warn "Stack layers not integer value, rounding"
    # end

    stack_layers = stack_layers_tmp#round.(Float64.(stack_layers_tmp))

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

"""
readNuMadMaterialsCSV(NuMad_materials_xlscsv_file)

Parameters defining the rotor materials.

**Arguments**
- `NuMad_materials_xlscsv_file::String`: name of the numad excel CSV file being read (!!! THE NUMAD TAB MUST BE SAVED AS A CSV FOR THIS TO WORK !!!)


**Returns**
- `Output::plyproperties`: plyproperties structure as defined in the plyproperties structure docstrings.
"""
function readNuMadMaterialsCSV(NuMad_materials_xlscsv_file)


    csvdata = DelimitedFiles.readdlm(NuMad_materials_xlscsv_file,',',skipstart = 0)

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
    plythickness = Float64.(csvdata[data_start:data_end,4]) .* 1e-3 #meters
    e1 = Float64.(csvdata[data_start:data_end,5]) .* 1e6
    e2 = Float64.(csvdata[data_start:data_end,6]) .* 1e6
    g12 = Float64.(csvdata[data_start:data_end,8]) .* 1e6
    anu = Float64.(csvdata[data_start:data_end,11])  #ratio
    rho = Float64.(csvdata[data_start:data_end,14])  #g/cc * 1000 #kg/m3
    xt = abs.(Float64.(csvdata[data_start:data_end,15]) .* 1e6)  #pa
    xc = abs.(Float64.(csvdata[data_start:data_end,16]) .* 1e6)  #pa, abs since composites is looking for positive failure values and handles the negative.
    if length(csvdata[4,:])>17
        yt = abs.(Float64.(csvdata[data_start:data_end,17]) .* 1e6)  #pa
        yc = abs.(Float64.(csvdata[data_start:data_end,18]) .* 1e6)  #pa, abs since composites is looking for positive failure values and handles the negative.
        costs = Float64.(csvdata[data_start:data_end,21]) #$/kg
        try
            SN_stressMpa = Float64.(csvdata[data_start:data_end,23:28])
            Log_SN_cycles2Fail = Float64.(csvdata[data_start:data_end,29:34])
        catch
            SN_stressMpa = collect(LinRange(1e12,0,6))
            Log_SN_cycles2Fail = collect(LinRange(0,7,6))
            @warn "Data for SN curve control points not found in material file columns 23:28 for stress in Mpa, 29:33 for cycles in log10"
        end
    else
        yt = ones(length(e1)) .* 100.0e6 #made up
        yc = ones(length(e1)) .* 100.0e6  #made up, abs since composites is looking for positive failure values and handles the negative.
        costs = zeros(length(e1)) #$/kg
        SN_stressMpa = collect(LinRange(1e12,0,6))
        Log_SN_cycles2Fail = collect(LinRange(0,7,6))
    end

    s = abs.(ones(length(e1)) .* 100.0e6)  #made up, abs since composites.jl is expecting positive failure values

    return plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness),costs,SN_stressMpa,Log_SN_cycles2Fail)

end


function saveOWENSfiles(filename,mymesh,myort,myjoint,myEl,pBC,numadIn_bld)
    # Create Filenames

    # Save mesh
    DelimitedFiles.open("$filename.mesh", "w") do io
           DelimitedFiles.writedlm(io, [mymesh.numNodes mymesh.numEl], ' ')
           DelimitedFiles.writedlm(io, [round.(Int,mymesh.nodeNum) mymesh.x mymesh.y mymesh.z], ' ')
           DelimitedFiles.writedlm(io, [round.(Int,mymesh.elNum) zeros(Int,length(mymesh.elNum)).+2 round.(Int,mymesh.conn)], ' ')
           DelimitedFiles.writedlm(io, " ", ' ')
           DelimitedFiles.writedlm(io, [length(mymesh.meshSeg) [meshSeg for meshSeg in mymesh.meshSeg]'], ' ')
       end

    # Save El

    DelimitedFiles.open("$filename.el", "w") do io
        for ii = 1:mymesh.numEl
            DelimitedFiles.writedlm(io, [mymesh.z[ii]/maximum(mymesh.z) -(myEl.props[ii].ac[1].+0.5)./2 myEl.props[ii].twist[1] myEl.props[ii].rhoA[1] myEl.props[ii].EIyy[1] myEl.props[ii].EIzz[1] myEl.props[ii].GJ[1] myEl.props[ii].EA[1] myEl.props[ii].alpha1[1] myEl.props[ii].rhoIyy[1] myEl.props[ii].rhoIzz[1] 0.0 0.0 myEl.props[ii].zcm[1] myEl.props[ii].ycm[1] 0.0 myEl.props[ii].a[1]], ' ')
            DelimitedFiles.writedlm(io, [mymesh.z[ii]/maximum(mymesh.z) -(myEl.props[ii].ac[2].+0.5)./2 myEl.props[ii].twist[2] myEl.props[ii].rhoA[2] myEl.props[ii].EIyy[2] myEl.props[ii].EIzz[2] myEl.props[ii].GJ[2] myEl.props[ii].EA[2] myEl.props[ii].alpha1[2] myEl.props[ii].rhoIyy[2] myEl.props[ii].rhoIzz[2] 0.0 0.0 myEl.props[ii].zcm[2] myEl.props[ii].ycm[2] 0.0 myEl.props[ii].a[2]], ' ')
        end
       end

    # Save Joint
    DelimitedFiles.open("$filename.jnt", "w") do io
           DelimitedFiles.writedlm(io, myjoint, '\t')
       end

    # Save Ort
    DelimitedFiles.open("$filename.ort", "w") do io
           DelimitedFiles.writedlm(io, [myort.elNum[:,1] myort.Psi_d myort.Theta_d myort.Twist_d myort.Length myort.Offset'], ' ')
       end

    # Save pBC
    DelimitedFiles.open("$filename.bc", "w") do io
           DelimitedFiles.writedlm(io, pBC, '\t')
       end

    # Save Blade
    # Used
    chord = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span),numadIn_bld.chord,mymesh.structuralSpanLocNorm[1,:])

    bldArray = zeros(length(mymesh.structuralSpanLocNorm),16)
    row = 1
    for ibld = 1:length(mymesh.structuralSpanLocNorm[:,1])
        for ispan = 1:length(mymesh.structuralSpanLocNorm[1,:])
            bldArray[row,1] = ibld
            bldArray[row,2] = mymesh.structuralSpanLocNorm[ibld,ispan]
            bldArray[row,3] = mymesh.structuralNodeNumbers[ibld,ispan]
            bldArray[row,4] = mymesh.structuralElNumbers[ibld,ispan]
            bldArray[row,14] = chord[ispan]
            bldArray[row,16] = 2*pi
            row += 1
        end
    end
    DelimitedFiles.open("$filename.bld", "w") do io
           DelimitedFiles.writedlm(io, bldArray, '\t')
       end

    return nothing
end


"""

    readMesh(filename)

Reads the mesh file and stores data in the mesh object.

input:
* `filename::string` string containing mesh path/to/filename.mesh

output:
* `mesh::OWENSFEA.Mesh` see OWENSFEA.Mesh
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

    mesh = OWENSFEA.Mesh(nodeNum,
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
* `BC::OWENSFEA.BC_struct`:   see OWENSFEA.BC_struct, object containing boundary condition data
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

    BC = OWENSFEA.BC_struct(numpBC,
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

    numBlades = Int(maximum(bladeNum))
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
* `el::OWENSFEA.El`:       see OWENSFEA.El    element data object
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

    sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, numElements)

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

        sectionPropsArray[i] = OWENSFEA.SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

    end
    close(fid) #close element file

    nodeNum = Int.(bladeData_struct.nodeNum)  #node number associated with blade section
    elNum = Int.(bladeData_struct.elementNum)    #element number associated with blade section
    bladeData = bladeData_struct.remaining  #blade data

    chord = zeros(maximum(Int,nodeNum),1)
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


    # println("EIyz, rhoIyz deactivated")

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
    el = OWENSFEA.El(sectionPropsArray,elLen,psi,theta,roll,rotationalEffects)

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
    println("Generator File to Generator properties not yet implemented")

    genprops = 0.0 #if generator property file does not exist, set object to null

    return genprops
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
