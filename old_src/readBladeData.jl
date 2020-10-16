

mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
end

function readBladeData(filename)
    #readBladeDAta reads blade data
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
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
