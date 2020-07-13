function [bladeData] = readBladeData(filename)
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

a = importCactusFile(filename,0,60,16,'	');

bladeNum = a(:,1);

numBlades = max(bladeNum);
#     numStruts = min(bladeNum);
#     if(numStruts>0)
#         numStruts = 0;
#     else
#         numStruts = abs(numStruts);
#     end

strutStartIndex = 0;
for i=1:length(bladeNum)
    if(isnan(a(i,end)))
        strutStartIndex = i;
        break;
    end
end



if(strutStartIndex~=0)
    #         strutDataBlock = a(strutStartIndex:end,:);
    #         [strutEntries, ~] = size(strutDataBlock);
    #         numNodesPerStrut = strutEntries/numStruts;
    #         numElPerStrut = numNodesPerStrut - 1;
else
    [temp,~]=size(a);
    strutStartIndex = temp + 1;
end

bladeDataBlock = a(1:strutStartIndex-1,:);

bladeData.numBlades = numBlades;  #assign data to bladeData object #TODO: Should not be loading this file in multiple times
bladeData.bladeNum = bladeDataBlock(:,1);
bladeData.h = bladeDataBlock(:,2);
bladeData.nodeNum = bladeDataBlock(:,3);
bladeData.elementNum = bladeDataBlock(:,4);
bladeData.remaining = bladeDataBlock(:,5:end);


end