function hand_verif_eig()

M = diag([3 1 3 1]);
C = [0.4 0 -0.3 0; 0 0 0 0; -0.3 0 0.5 -0.2; 0 0 -0.2 0.2];
K = [-7 2 4  0; 2 -4 2 0; 4 2 -9 3; 0 0 3 -3];

X_truth = [0.1828, 0.3421, 0.3989,-0.0621, 0.3890,-0.4143, 0.4575,-0.4563;...
    0.3530,-0.9296, 0.3330, 0.8571,-0.6366,-0.2717, 0.4981,-0.4985;...
    -0.5360,-0.0456,-0.1724,-0.3509,-0.3423, 0.1666, 0.5106,-0.5107;...
    0.7448, 0.1295,-0.8368, 0.3720, 0.5712, 0.8525, 0.5309,-0.5315];

e_truth =[-2.4498;-2.1536;-1.6248; 2.2279; 2.0364; 1.4752; 0.3353;-0.3466];

[e_truth,eigsortidx] = sort(e_truth,'ComparisonMethod','abs');
X_truth = X_truth(:,eigsortidx);
%Normalize Eigenvectors
X_truth = normalizeeigenvectors(X_truth);

len=length(M);
eyelen = eye(len);
zeroslen = zeros(len);

sysMat = [zeros(len), eye(len);   %constructs state space form (with mass matrix inverted)
    -M\K, -M\C];

sysMat2 = [zeroslen, eyelen;      %construct state space form and lets eigs invert mass matrix
    -K, -C];

MMat2 = [eyelen,zeroslen;zeroslen,M];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% split matrices eig %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[eigVec0,eigVal0] = eig(sysMat2,MMat2,'vector');		  %full eigenvalue solve
toc
[eigVal0,eigsortidx] = sort(eigVal0,'ComparisonMethod','abs');
eigVec0 = eigVec0(:,eigsortidx);

residual0 = sysMat2*eigVec0 - MMat2*eigVec0*diag(eigVal0);
fprintf('%s','split matrix eig residual: ');
fprintf('%e\n',max(max(abs(residual0))));

%Normalize Eigenvectors
eigVec0 = normalizeeigenvectors(eigVec0);

n_mismatch = checkequal(abs(eigVec0),abs(X_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvectors split matrix eig');
end


n_mismatch = checkequal(abs(eigVal0),abs(e_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvalues split matrix eig');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% single matrix eig %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[eigVec1,eigVal1] = eig(sysMat,'vector');		  %full eigenvalue solve
toc
[eigVal1,eigsortidx] = sort(eigVal1,'ComparisonMethod','abs');
eigVec1 = eigVec1(:,eigsortidx);

residual1 = sysMat*eigVec1 - eigVec1*diag(eigVal0);
fprintf('%s','single matrix eig residual: ');
fprintf('%e\n',max(max(abs(residual1))));

%Normalize Eigenvectors
eigVec1 = normalizeeigenvectors(eigVec1);

n_mismatch = checkequal(abs(eigVec1),abs(X_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvectors single matrix eig');
end


n_mismatch = checkequal(abs(eigVal1),abs(e_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvalues single matrix eig');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% polyeig %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[eigVec2,eigVal2] = polyeig(K,C,M);
toc
[eigVal2,eigsortidx] = sort(eigVal2,'ComparisonMethod','abs');
eigVec2 = eigVec2(:,eigsortidx);

residual2 = zeros(1,length(eigVal2));
for i = 1:length(eigVal2)
lambda = eigVal2(1);
x = eigVec2(:,1);
residual2(i) = real(max((M*lambda^2 + C*lambda + K)*x));
end

fprintf('%s','polyeig residual: ');
fprintf('%e\n',max(abs(residual2)));

%Normalize Eigenvectors
eigVec2 = normalizeeigenvectors(eigVec2);

n_mismatch = checkequal(abs(eigVec2),abs(X_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvectors polyeig');
end


n_mismatch = checkequal(abs(eigVal2),abs(e_truth));
if n_mismatch~=0
    fprintf('%s\n','MISMATCH: Eigenvalues polyeig');
end




% Plot Eigenvectors
figure()
plot(abs(X_truth(:,1)))
hold on
plot(abs(eigVec0(1:4,1)))
hold on
plot(abs(eigVec0(5:8,1)))
hold on
plot(abs(eigVec1(1:4,1)))
hold on
plot(abs(eigVec1(5:8,1)))
hold on
plot(abs(eigVec2(:,1)))
legend('Truth','Eigs Dual1','Eigs Dual12','Eigs single1','Eigs single2','Polyeig')
end

function [n_mismatch] = checkequal(first,second)
tol = 1e-3;
n_mismatch = 0;
if isequal(first,second)
    fprintf('%s\n', 'Elements are equal')
else
    for i = 1:length(second(:,1))
        for j = 1:length(second(1,:))
            if abs(abs(first(i,j))-abs(second(i,j)))>tol
                fprintf('%f\n',i);
                fprintf('%f\n',j);
                fprintf('%f\n',first(i,j));
                fprintf('%f\n',second(i,j));
                n_mismatch = n_mismatch+1;
            end
        end
    end
end

end

function [eigVec] = normalizeeigenvectors(eigVec)
for i=1:length(eigVec(1,:))
    eigVec(1:4,i) = eigVec(1:4,i)/max(abs(eigVec(1:4,i)));
end

if length(eigVec(:,1))>4
    for i=1:length(eigVec(1,:))
        eigVec(5:8,i) = eigVec(5:8,i)/max(abs(eigVec(5:8,i)));
    end
end
end
