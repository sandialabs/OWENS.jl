function [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
%readBDdata  reads boundary condition file
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
%
%   This function reads the boundray condition file and stores data in the
%   boundary condition object.
%
%      input:
%      bcfilename    = string containing boundary condition filename
%      numNodes      = number of nodes in structural model
%      numDofPerNode = number of degrees of freedom per node

%      output:
%      BC            = object containing boundary condition data

fid = fopen(bcfilename);       %open boundary condition file
numpBC = real(str2double(myfgetl(fid)));   %read in number of boundary conditions (displacement boundary conditions)
pBC = zeros(numpBC,3);         %initialize boundary conditions
for i=1:numpBC

    line = myfgetl(fid);

    % Find where all of the delimiters are
    %first two are boundary condition node number and local DOF number
    %third is boundary condition value (typically zero)
    delimiter_idx = [0.0,find(line == ' '),length(line)+1];
    % Extract the data from the beginning to the last delimiter
    for k = 2:length(delimiter_idx)
        pBC(i,k-1) = real(str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1)));
    end

end

totalNumDof = numNodes*numDofPerNode;

BC = struct('numpBC',zeros,'pBC',zeros,'numsBC',zeros,'nummBC',zeros,'isConstrained',zeros(totalNumDof,1));
BC.numpBC = numpBC;  %store boundary condition data  in boundayr condition object
BC.pBC = pBC;
BC.numsBC = 0;
BC.nummBC = 0;

fclose(fid);

%create a vector denoting constrained DOFs in the model (0 unconstrained, 1
%constrained)


%calculate constrained dof vector
isConstrained = zeros(totalNumDof,1);
constDof = (BC.pBC(:,1)-1)*numDofPerNode + BC.pBC(:,2);
index = 1;
for i=1:numNodes
    for j=1:numDofPerNode
        if(ismember((i-1)*numDofPerNode + j,constDof))
            isConstrained(index) = 1;
        end
        index = index + 1;
    end
end
BC.isConstrained = isConstrained;

end
