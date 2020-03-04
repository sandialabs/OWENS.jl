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
numpBC = fscanf(fid,'%i',1);   %read in number of boundary conditions (displacement boundary conditions)
pBC = zeros(numpBC,3);         %initialize boundary conditions
for i=1:numpBC
    
   dat1 = fscanf(fid,'%i',2); %boundary condition node number and local DOF number
   dat2 = fscanf(fid,'%f',1); %boundary condition value (typically zero)
   
   pBC(i,1:3) = [dat1(1) dat1(2) dat2]; %store boundary condition data in array
   
end

BC.numpBC = numpBC;  %store boundary condition data  in boundayr condition object
BC.pBC = pBC;
BC.numsBC = 0;
BC.nummBC = 0;

fclose(fid);

%create a vector denoting constrained DOFs in the model (0 unconstrained, 1
%constrained)

totalNumDof = numNodes*numDofPerNode;
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
