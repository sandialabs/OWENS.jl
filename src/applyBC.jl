function [Kg,Fg] = applyBC(Kg,Fg,BC,numDofPerNode)
#applyBC Applies boundary conditions to system for static analysis
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [Kg,Fg] = applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode)
#
#   This function applies boundary conditions to the stiffness matrix and
#   load vector for a static analysis.
#
#      input:
#      Kg            = assembled global stiffness matrix
#      Fg            = assembled global load vector
#      BC            = struct of boundary condition information
#      u             = global displacement vector
#      iterationType = for nonlinear analysis, not used in BLAST
#      numDofPerNode = number of degrees of freedom per node

#      output:
#      Kg            = global stiffness matrix with boundary conditions
#      Fg            = global load vector with boundary condition


[numEq,~]=size(Kg);

#APPLY BCs FOR PRIMARY VARIABLE

if(BC.numpBC > 0)
    pBC = BC.pBC;
    [numpBC,~] = size(pBC);
    
    for i=1:numpBC
        nodeNumber = pBC(i,1);
        dofNumber = pBC(i,2);
        specVal = pBC(i,3);
        
        eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
        
        for j=1:numEq
            Kg(eqNumber,j) = 0.0;
            Fg(j) = Fg(j) - Kg(j,eqNumber)*specVal;
            Kg(j,eqNumber) = 0.0;
        end
        Fg(eqNumber) = specVal;
        Kg(eqNumber,eqNumber) = 1.0;
    end
end

#APPLY BCs FOR SECONDARY VARIABLE

# if(BC.numsBC > 0) # This does not appear to be used
#     sBC = BC.sBC;
#     [numsBC,~] = size(sBC);
#     
#     for i=1:numsBC
#         nodeNumber = sBC(i,1);
#         dofNumber = sBC(i,2);
#         specVal =  sBC(i,3);
#         
#         eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
#         
#         Fg(eqNumber) = Fg(eqNumber) + specVal;
#         
#     end
# end

end