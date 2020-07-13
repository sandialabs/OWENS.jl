function [sortedPhi] = sortEigOutput(eigValVec,eigVec,numModesSelected)
#sortEigOutput sorts output from eigensolve in ascending order
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [sortedPhi] = sortEigOutput(eigValVec,eigVec,numModesSelected)
#                    
#   This function reads the nodal terms file and stores data in the nodal
#   terms object.
#
#      input:
#      eigValVec        = vector of eigenvalues
#      eigVec           = matrix of eigenvectors(eigenvectors as columns)
#      numModesSelected = number of lower modes selected
#
#      output:
#      sortedPhi        = outputs an n x m sorted modal matrix of
#                         eigenvectors. n=number of DOFs, m = number of
#                         modes selected for modal matrix
#                   

[~,map] = sort(eigValVec); #sorts eigValVec in ascending order

sortedPhi = zeros(length(eigVec),numModesSelected); #initializes sortedPhiArray

for i=1:numModesSelected  #creates sortedPhi array with columns of eigenvectors corresponding to increasing eigenvalues
    j = map(i);
    sortedPhi(:,i) = eigVec(:,j);
end

end

