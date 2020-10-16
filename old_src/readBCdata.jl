mutable struct BC_struct
    numpBC
    pBC
    numsBC
    nummBC
    isConstrained
    map
    redVectorMap
end

function readBCdata(bcfilename,numNodes,numDofPerNode)
    #readBDdata  reads boundary condition file
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
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
