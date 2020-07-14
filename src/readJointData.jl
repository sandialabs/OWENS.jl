function readJointData(inputfile)
    #readJointData reads joint data file
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #   [joint] = readJointData(inputfile)
    #
    #   This function reads the joint data file and stores data in the joint
    #   object.
    #
    #      input:
    #      inputfile     = string containing joint data filename

    #      output:
    #      joint         = array containing joint data

    fid = open(inputfile) #open joint file
    file_length = 0
    while (!eof(fid))
        readline(fid)
        file_length = file_length+1
    end
    close(fid)

    joint = zeros(file_length,8)

    fid = open(inputfile) #open joint file
    count = 1
    while (!eof(fid))
        #read in nodal info associated with joint [joint #, master node #, slave node #, joint type]

        line = readline(fid)

        # Find where all of the delimiters are
        delimiter_idx = [0;collect.(Int,findall("	",line));length(line)+1]
        lineinfo = zeros(length(delimiter_idx)-1,1)
        # Extract the data from the beginning to the last delimiter
        for k = 2:length(delimiter_idx)
            lineinfo[k-1] = (parse(Float64,line[delimiter_idx[k-1][1]+1:delimiter_idx[k][1]-1]))
        end

        nodalInfo = lineinfo[1:4]

        if (isempty(nodalInfo))
            println("No joint data detected") #if no joint data detetcted return null joint object
            joint = []
            break
        end

        #reads in mass, stiffness, orientation of element attached to master joint
        massStiffOrt = lineinfo[5:8]
        temp = cat(nodalInfo',massStiffOrt',dims=2)
        joint[count,:] = temp        #store data in joint array
        count = count + 1

    end
    close(fid)    #close file

    return joint

end
