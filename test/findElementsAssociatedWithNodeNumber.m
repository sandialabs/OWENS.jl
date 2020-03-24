function [elList,localNode] = findElementsAssociatedWithNodeNumber(nodeNum,conn,jointData)
%findElementsAssociatedWithNodeNumber calculates reaction force at a node
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [elList,localNode] = findElementsAssociatedWithNodeNumber(nodeNum,...
%                        conn,jointData)
%
%   This function finds elements associated with a node number through
%   mesh connectivity or joint constraints
%
%   input:
%   nodeNum    = node number joint constraints are desired at
%   conn       = object containing mesh connectivity
%   jointData  = object containing joint informatino

%
%   output:
%   elList     = array containing a list of element numbers associated with
%                nodeNum
%   localNode  = array containing the local node number that correspond to
%                nodeNum in the list of associated elements

%search joint constraints
index = 1;
[numEl,~] = size(conn); %get number of elements in mesh
elList = zeros(1,1);
localNode = zeros(1,1);
if(~isempty(jointData))
    %first see if specified node is a slave node in a joint constraint
    res2 = find(ismember(jointData(:,3),nodeNum)); %search joint data slave nodes for node number
    %if it is, change it to the corresponding master node
    if(~isempty(res2))
        nodeNum = jointData((end),2);
        if length(jointData)>1
            error('Incorrect Joint Data and nodeNum, too many joints')
        end
    end
    
    res1 = find(ismember(jointData(:,2),nodeNum)); %search joint data master nodes for node number
    if ~isempty(res1)
        jointNodeNumbers = jointData(res1,3);
        
        elList = zeros(1,length(jointNodeNumbers)*numEl);
        localNode = zeros(length(jointNodeNumbers)*numEl);
        for j=1:length(jointNodeNumbers) %loop over joints
            for i=1:numEl
                res = ismember(conn(i,:),jointNodeNumbers(j)); %finds indices of nodeNum in connectivity of element i
                localNodeNumber = find(res); %finds the local node number of element i that corresponds to nodeNum
                if(localNodeNumber~=0) %assigns to an elementList and localNode list
                    elList(index) = i;
                    localNode(index) = localNodeNumber;
                    index = index + 1;
                end
            end
        end
    end
else
    error('empty jointData')
end


for i=1:numEl %loop over elements
    res = ismember(conn(i,:),nodeNum); %finds indices of nodeNum in connectivity of element i
    localNodeNumber = find(res); %finds the local node number of element i that corresponds to nodeNum
    if(localNodeNumber~=0) %assigns to an elementList and localNode list
        elList(index) = i;
        localNode(index) = localNodeNumber;
        index = index + 1;
    end
end


end
