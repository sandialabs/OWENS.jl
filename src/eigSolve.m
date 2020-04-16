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


%     if(meirovitchFlag)
%         disp('Solving with Meirovitch approach ...');
%         [eigVec,eigVal] = eigSolveMeirovitch(M,C,K,flag);
%     else
%         disp('Solving with standard state space approach ...');

len=length(M);
% eyelen = eye(len);
% zeroslen = zeros(len);

%%% if(flag == 1)
sysMat = [zeros(len), eye(len);   %constructs state space form (with mass matrix inverted)
    -M\K, -M\C];
%%% end

%%% if(flag == 2)
% sysMat2 = [zeroslen, eyelen;      %construct state space form and lets eigs invert mass matrix
%     -K, -C];
% sysMat2 = sparse(sysMat2);

% MMat2 = [eyelen,zeroslen;zeroslen,M];
% MMat2 = sparse(MMat2);
%%% end


%%% if(flag==1)
[eigVec,eigVal] = eig(sysMat);%,'vector');		  %full eigenvalue solve
% [eigVec0,eigVal0] = eig(sysMat2,MMat2,'vector');		  %full eigenvalue solve

[eigVec2,eigVal2] = polyeig(K,C,M);
[eigVal2,eigsortidx] = sort(eigVal2,'ComparisonMethod','abs');
eigVec2 = [eigVec2(:,eigsortidx);eigVec2(:,eigsortidx)];
eigVal2 = diag(eigVal2);


% residual0 = sysMat*eigVec0 - eigVec0*diag(eigVal0);
% residual0 = sysMat2*eigVec0 - MMat2*eigVec0*diag(eigVal0);
% fprintf('%e\n',real(residual0(1)))
% [sorted_eigVal,eigsortidx] = sort(eigVal0,'ComparisonMethod','abs');
% eigVec = eigVec0(:,eigsortidx);
% eigVal = diag(sorted_eigVal); % return to diag
% residual = sysMat*eigVec - eigVec*eigVal;
% output = max(max(real(residual)));
% fprintf('%e\n',output)
%%% end

%%% if(flag==2)
% [eigVec,eigVal] = eigs(sysMat2,MMat2,828,'SM');  %subest of modes for eigenvalue solve
% test = sysMat2*eigVec - MMat2*eigVec*eigVal;
% output = max(max(real(test)));
% fprintf('%e\n',output)
%%% end

%%% if(flag==3)
%     sysMat=inv(M)*K;                      %eigenvalue solve on spring mass system only
%     [eigVec,eigVal] = eig(sysMat);
%     [eigVec] = sortEigOutput(diag(eigVal),eigVec,numModes); %eigenvalues/vectors sorted in ascending frequency before returning
%%% end


end

% function [eigVec,eigVal] = eigSolveMeirovitch(M,C,K,flag)
%     	[len,dum] = size(M);
%
%     %Force symmetry and skewsymmetry
%     K=0.5*(K+K');
%     M=0.5*(M+M');
%     for i=1:len
%        for j=1:i
%            if(i~=j)
%            C(i,j) = -C(j,i);
%            end
%        end
%     end
%
%     Mstar = [M, zeros(len);
%              zeros(len), K];
%
%     Gstar = [C, K;-K,zeros(len)];
%
%     Rstar = chol(Mstar);
%     invRstar = inv(Rstar);
%
%     Kstar = Gstar'*(inv(Mstar))*Gstar;
%     Kstar = 0.5*(Kstar+Kstar');
%
% % 	sysMat = [zeros(len), eye(len);
% % 			  -(M^-1)*K, -(M^-1)*C];
%
%     sysMat = invRstar'*Kstar*invRstar;
%     %force symmetry of system matrix
%     sysMat = 0.5*(sysMat'+sysMat);
%
%     if(flag==1)
%         [eigVec,eigVal] = eig(sysMat);
%     end
%     if(flag==2)
%          sysMat=sparse(sysMat);
%         [eigVec,eigVal] = eigs(sysMat,20,'SM');
%     end
%
%
%     %sort and reconstruct typical vectors from meirovitch solve
%     [len,dum]=size(eigVal);
%     for i=1:len
%        valtemp(i)=eigVal(i,i);
%     end
%
%     [valtemp,map,posIndex] = bubbleSort(valtemp);
%
%     len = 20;
%     for i=1:len
%         vecNew(:,i)=eigVec(:,map(i));
%     end
%
%
%     for i=1:len
% %         eigVec(:,i) = inv(Rstar)*eigVec(:,i);
%        vecNew(:,i) = invRstar*vecNew(:,i);
%     end
%
%     eigVec=vecNew;
%
%     for i=1:len
%     eigValNew(i,i)=valtemp(i);
%     end
%     eigVal = eigValNew;
%
%     for i=1:len/2
%         index1=(i-1)*2+1;
%         index2=(i-1)*2+2;
%
%         vec1 = eigVec(:,index1);
%         vec2 = eigVec(:,index2);
%
%         temp1 = vec1 + 1i.*vec2;
%         temp2 = vec1 - 1i.*vec2;
%
%         eigVec(:,index1) = temp1;
%         eigVec(:,index2) = temp2;
%     end
% end
