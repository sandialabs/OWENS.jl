function [bcMap] = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
%calculateBCMap   calculates a boundary condition map
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [bcMap] = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
%
%   This function creates a boundary condition map between full and reduced
%   dof listing as a result of constraints.
%
%      input:
%      numpBC            = number of boundary conditions
%      pBC               = array of boundary  condition data
%      numDofPerNode     = number of degrees of freedom per node
%      reducedDofList    = array of reduced DOF numbering
%
%      output:
%      elStorage         = map for boundary conditions between full and
%                          reduced dof list


constrainedDof = zeros(numpBC,1);
for i=1:numpBC
    constrainedDof(i) = (pBC(i,1)-1)*numDofPerNode + pBC(i,2);  %creates an array of constrained DOFs
end
constrainedDof = sort(constrainedDof);

reducedDOFCount = length(reducedDofList);

bcMap = zeros(reducedDOFCount,1);
index = 1;
for i=1:reducedDOFCount
    if(ismember(reducedDofList(i),constrainedDof))  %searches reduced DOF for constrained DOFs
        bcMap(i) = -1;
    else
        bcMap(i) = index;
        index = index + 1;
    end
end

end
