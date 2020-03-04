function [bladeData] = readBladeData(filename)
%readBladeDAta reads blade data
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [bladeData] = readBladeData(filename)
%                    
%   This function reads blade data from file
%
%   input:
%   filename      = string containing file name for blade data file
%
%   output:
%   bladeData     = object containing blade data
try
    fid=fopen(filename,'r'); %open blade data file
    index = 1;
    while(~feof(fid)) %if blade data file exists
        bladeNumber = fscanf(fid,'%i',1);
        if(isempty(bladeNumber))
            break;
        end
        if(bladeNumber<0)
            break;
        end
        bladeNum(index,1) = bladeNumber;          %read blade number from file
        h(index,1) = fscanf(fid,'%f',1);          %read blade spanwise coordinates from file
        nodeNum(index,1) = fscanf(fid,'%i',1);    %read blade spanwise  coordinate node numbers from file
        elementNum(index,1) = fscanf(fid,'%i',1); %read blade section element numbers from file
        dump = fscanf(fid,'%f',12);               %read misc. data in blade file
        index = index + 1;
    end
    
    len=length(bladeNum);           %calculate total number of entries in blade file
    numBlades = bladeNum(len);      %extract total number of blades
    %assumes same # of blade nodes per blade
    nodesPerBlade = len/numBlades;  %calculate number of nodes per blade

    bladeData.numBlades = numBlades;  %assign data to bladeData object
    bladeData.bladeNum = bladeNum;
    bladeData.h = h;
    bladeData.nodeNum = nodeNum;
    bladeData.elementNum = elementNum;
catch
    bladeData = [];
end
end