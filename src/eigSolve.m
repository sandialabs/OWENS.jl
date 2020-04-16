function [eigVec,eigVal] = eigSolve(M,C,K)%,numModes,flag)
%eigSolve   Calculates eigenvalues and vectors of structural dynamics rep.
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [eigVec,eigVal = eigSolve(M,C,K,numModes,flag)
%
%   This function calculates the eigenvalues and vectors of a structural
%   dynamic system.
%
%   input:
%   M         = system mass matrix
%   C         = system damping matrix
%   K         = system stiffness matrix
%   numModes  = number of lower system modes to extract
%
%   output:
%   eigVal    = diagonal matrix of eigenvalues
%   eigVec    = matrix of eigenvectors as columns
%   flag      = directs type of eigensolve
%                1 = all eigenvalues extracted
%                2 = subset of eigenvalues extracted

len=length(M);

% sysMat = zeros(828,828);
% sysMat(1:414,415:828) = eye(414);
% sysMat(415:828,1:414) = -M\K;
% sysMat(415:828,415:828) = -M\C;

sysMat = [zeros(len), eye(len); %constructs state space form (with mass matrix inverted)
    -M\K, -M\C];

[eigVec,eigVal] = eig(sysMat); %full eigenvalue solve.  Tried splitting the mass matrix out, was about equivalent in matlab, but nearly 2x the time in c++, also tried polyeig, also 2x the time.  EIGS not supported for C++ translation.


%%% if(flag==3)
%     sysMat=inv(M)*K;                      %eigenvalue solve on spring mass system only
%     [eigVec,eigVal] = eig(sysMat);
%     [eigVec] = sortEigOutput(diag(eigVal),eigVec,numModes); %eigenvalues/vectors sorted in ascending frequency before returning
%%% end


end
