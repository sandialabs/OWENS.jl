function [Fnew] = applyBCModalVec(F,numpBC,bcMap)
#applyBCModal Applies boundary conditions to system for modal analysis
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
#                    
#   This function applies boundary conditions to a system matrix for modal
#   analysis
#
#      input:
#      K             = assembled global system matrix
#      BC            = struct of boundary condition information
#      numDofPerNode = number of degrees of freedom per node
 
#      output:
#      K             = global system matrix with boundary conditions
#      dofVector     = reduced DOF vector after imposing BCs

[numEq,dum]=size(F);
Fnew = zeros(numEq - numpBC,1);
indVec = zeros(numEq-numpBC,1);

index = 1;
for i=1:numEq
    if(bcMap(i) != -1)
        indVec(index) = i;
        index = index +1;
    end
end

#APPLY BCs FOR PRIMARY VARIABLE
if(numpBC > 0)
    Fnew = F(indVec);
else
    Fnew = F;
end

end




