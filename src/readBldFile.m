
function [structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers] = readBldFile(OWENSfile_root)
bldFn = [OWENSfile_root '.bld'];
%% READ IN BLD FILE
a = importCactusFile(bldFn,0,60,16,'	');

bladeNum = a(:,1);

numBlades = max(bladeNum);
% numStruts = min(bladeNum);
% if(numStruts>0)
% %     numStruts = 0;
% else
% %     numStruts = abs(numStruts);
% end

strutStartIndex = 0;
for i=1:length(bladeNum)
    if(isnan(a(i,end)))
        strutStartIndex = i;
        break;
    end
end



if(strutStartIndex~=0)
    %     strutDataBlock = a(strutStartIndex:end,:);
    %     [strutEntries, ~] = size(strutDataBlock);
    %     numNodesPerStrut = strutEntries/numStruts;
    %     numElPerStrut = numNodesPerStrut - 1;
else
    [temp,~]=size(a);
    strutStartIndex = temp + 1;
end

bladeDataBlock = a(1:strutStartIndex-1,:);
[bladeEntries, ~] = size(bladeDataBlock);
numNodesPerBlade = bladeEntries/numBlades;
% numElPerBlade = numNodesPerBlade - 1;

% [len,~]=size(bladeDataBlock);
structuralSpanLocNorm = zeros(numBlades,numNodesPerBlade);
structuralNodeNumbers = zeros(numBlades,numNodesPerBlade);
structuralElNumbers = zeros(numBlades,numNodesPerBlade);
for i=1:numBlades
    structuralSpanLocNorm(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,2)./bladeDataBlock(i*numNodesPerBlade,2);
    structuralNodeNumbers(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,3);
    structuralElNumbers(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,4);
end

bladeData.numBlades = numBlades;  %assign data to bladeData object
bladeData.bladeNum = bladeDataBlock(:,1);
bladeData.h = bladeDataBlock(:,2);
bladeData.nodeNum = bladeDataBlock(:,3);
bladeData.elementNum = bladeDataBlock(:,4);
bladeData.remaining = bladeDataBlock(:,5:end);
%%
end
