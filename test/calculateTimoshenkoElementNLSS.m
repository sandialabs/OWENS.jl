function [output] = calculateTimoshenkoElementNLSS(input)
%calculateTimoshenkoElementNLSS performs selective nonlinear element calculations
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [output] = calculateTimoshenkoElementNLSS(input)
%
%   This function performs selective nonlinear element calculations. Only
%   stiffness matrix contributions are evaluate. No other calculations are
%   performed to facilitate efficiency gains.
%
%      input:
%      input      = object containing element input
%
%      output:
%      output     = object containing element data
%-------- assign input block ----------------
elementOrder   = input.elementOrder;
x              = input.x;
xloc           = input.xloc;
disp           = input.disp;
sectionProps   = input.sectionProps;
sweepAngle     = input.sweepAngle;
coneAngle      = input.coneAngle;
rollAngle      = input.rollAngle;

useDisp        = input.useDisp;
preStress      = input.preStress;
iterationType  = input.iterationType;

%--------------------------------------------
if(strcmp(input.analysisType,'M'))
    disp_iter=disp;
end


numGP = 1; %used reduced integration for nonlinear terms

%calculate quad points
[xi,weight] = getGP(numGP);

%Initialize element sub matrices and sub vectors
numNodesPerEl = length(x);

K22NL = zeros(numNodesPerEl);
K33NL = K22NL;
K12NL = K22NL;
K13NL = K22NL;
K23NL = K22NL;
%     SS33 = SS22;

%Convert frequencies from Hz to radians
%Sort displacement vector
%Written for 2 node element with 6 dof per node
twistAvg = rollAngle + 0.5*(sectionProps.twist(1) + sectionProps.twist(2));
lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0);

dispLocal = lambda*disp_iter(1:12)';

uNode = [dispLocal(1) dispLocal(7)];
vNode = [dispLocal(2) dispLocal(8)];
wNode = [dispLocal(3) dispLocal(9)];

%Integration loop
for i=1:numGP
    %Calculate shape functions at quad point i
    [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi(i),xloc);
    p_N1_x = p_N_x;
    p_N2_x = p_N_x;
    p_N3_x = p_N_x;
    
    integrationFactor = Jac * weight(i);
    
    %..... interpolate for value at quad point .....
    if(useDisp || preStress)
        EA   = interpolateVal(sectionProps.EA,N);
        
        uprime = interpolateVal(uNode,p_N1_x);
        vprime = interpolateVal(vNode,p_N2_x);
        wprime = interpolateVal(wNode,p_N3_x);
    end
    
    %nonlinear element calculations
    if(useDisp)
        [K12NL] = calculateElement1(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12NL);
        [K13NL] = calculateElement1(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13NL);
        [K22NL] = calculateElement1(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22NL);
        [K23NL] = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NL);
        [K33NL] = calculateElement1(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33NL);
        
        %K12NLhat = K12;
        %K13NLhat = K13;
        %nonlinear element tangent matrix component calculations
        % T_ij = K_ij + Khat_ij
        if(strcmp(iterationType,'NR'))
            [K22NLhat] = calculateElement1(EA*(uprime + vprime^2 + 0.5*wprime^2),integrationFactor,p_N2_x,p_N2_x,K22NLhat);
            [K33NLhat] = calculateElement1(EA*(uprime + wprime^2 + 0.5*vprime^2),integrationFactor,p_N3_x,p_N3_x,K33NLhat);
            [K23NLhat] = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NLhat);
        end
        
    end
    
    
end %END OF INTEGRATION LOOP

%---------------------------------------------

Ktemp = zeros(numNodesPerEl*6);
Ktemp(1:2,3:4) = K12NL;
Ktemp(1:2,5:6) = K13NL;
Ktemp(3:4,1:2) = 2*K12NL;
Ktemp(5:6,1:2) = 2*K13NL;
Ktemp(3:4,3:4) = K22NL;
Ktemp(5:6,5:6) = K33NL;
Ktemp(3:4,5:6) = K23NL;
Ktemp(5:6,3:4) = K23NL;

[Ke] = mapMatrixNonSym(Ktemp);

% transform matrices for sweep (currently hardcoded to 0 sweep)
% Note,a negative theta3, will sweep away from the direction of
% rotation
lambdaTran = lambda';
Ke = sparse(Ke);
lambda = sparse(lambda);
lambdaTran = sparse(lambdaTran);
Ke = lambdaTran*Ke*lambda;

if(strcmp(iterationType,'NR'))
    error('calcTimoElNLSS need some mods to be used with newton raphson');
end
%----- assign output block ----------------
output.Ke = full(Ke);

%------------------------------------------

end

function [valGP] = interpolateVal(valNode,N)
%This function interpolates a value using distinct values at valNode
%and the corresponding shape function N.
valGP = 0.0;
for i=1:length(N)
    valGP = valGP + N(i)*valNode(i);
end
end

%Element calculation functions---------------------------------

function [K] = calculateElement1(EA,integrationFactor,N1,N2,K)
%This function is a general routine to calculate an element matrix
len1 = length(N1);
len2 = length(N2);
for i=1:len1
    for j=1:len2
        K(i,j) = K(i,j) + EA*N1(i)*N2(j)*integrationFactor;
    end
end

end

function [Kel] = mapMatrixNonSym(Ktemp)
%----- function to form total stifness matrix and transform to desired
% DOF mapping

T = [1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

%map to FEA numbering
Kel = T'*Ktemp*T;

%declare map
% map = [1, 7, 2, 8, 3, 9,...
%       4, 10, 5, 11, 6, 12];
%
% %map to FEA numbering
% for i=1:a
%     I=map(i);
%     for j=1:a
%         J=map(j);
%         Kel(I,J) = Ktemp(i,j);
%     end
% end

end