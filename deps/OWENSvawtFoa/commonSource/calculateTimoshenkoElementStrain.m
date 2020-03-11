function [output] = calculateTimoshenkoElementStrain(input)
%calculateTimoshenkoElementNL performs nonlinear element calculations
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [output] = calculateTimoshenkoElementNL(input,elStorage)
%
%   This function performs nonlinear element calculations.
%
%      input:
%      input      = object containing element input
%      elStorage  = obect containing precalculated element data
%
%      output:
%      output     = object containing element data

%-------- assign input block ----------------
elementOrder   = input.elementOrder;
xloc           = input.xloc;
disp           = input.disp;
sectionProps   = input.sectionProps;
sweepAngle     = input.sweepAngle;
coneAngle      = input.coneAngle;
rollAngle      = input.rollAngle;
nlOn = input.nlOn;
%--------------------------------------------

numGP = 4;   %number of gauss points for full integration
%calculate quad points
[xi,~] = getGP(numGP);

p_disp_x = zeros(numGP,6);

%Initialize element sub matrices and sub vectors

%Sort displacement vector
%Written for 2 node element with 6 dof per node
twistAvg = rollAngle + 0.5*(sectionProps.twist(1) + sectionProps.twist(2));
lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0);

dispLocal = lambda*disp';    %'

uNode = [dispLocal(1) dispLocal(7)];
vNode = [dispLocal(2) dispLocal(8)];
wNode = [dispLocal(3) dispLocal(9)];
theta_xNode = [dispLocal(4)  dispLocal(10)];
theta_yNode = [dispLocal(5)  dispLocal(11)];
theta_zNode = [dispLocal(6)  dispLocal(12)];

%Integration loop
eps_xx_0 = zeros(1,numGP);
eps_xx_z = zeros(1,numGP);
eps_xx_y = zeros(1,numGP);
gam_xz_0 = zeros(1,numGP);
gam_xz_y = zeros(1,numGP);
gam_xy_0 = zeros(1,numGP);
gam_xy_z = zeros(1,numGP);
for i=1:numGP
    %Calculate shape functions at quad point i
    [N,p_N_x,~] = calculateShapeFunctions(elementOrder,xi(i),xloc);
    %N1 = N;
    %N2 = N;
    %N3 = N;
    %N4 = N;
    N5 = N;
    N6 = N;
    p_N1_x = p_N_x;
    p_N2_x = p_N_x;
    p_N3_x = p_N_x;
    p_N4_x = p_N_x;
    p_N5_x = p_N_x;
    p_N6_x = p_N_x;

    %calculate displacement derivatives at quad point i
    uprime = interpolateVal(uNode,p_N1_x);
    vprime = interpolateVal(vNode,p_N2_x);
    wprime = interpolateVal(wNode,p_N3_x);
    theta_x_prime = interpolateVal(theta_xNode,p_N4_x);
    theta_y_prime = interpolateVal(theta_yNode,p_N5_x);
    theta_y_gp = interpolateVal(theta_yNode,N5);
    theta_z_prime = interpolateVal(theta_zNode,p_N6_x);
    theta_z_gp = interpolateVal(theta_zNode,N6);
    p_disp_x(i,:) = [uprime, vprime, wprime, theta_x_prime, theta_y_prime, theta_z_prime];

    if(nlOn)
        eps_xx_0(i) = uprime + 0.5*(wprime^2 + vprime^2);
    else
        eps_xx_0(i) = uprime;
    end
    eps_xx_z(i) = theta_y_prime;
    eps_xx_y(i) = -theta_z_prime;
    gam_xz_0(i) =  theta_y_gp + wprime;
    gam_xz_y(i) =  theta_x_prime;
    gam_xy_0(i) = -theta_z_gp + vprime;
    gam_xy_z(i) =  -theta_x_prime;
end %END OF INTEGRATION LOOP

output.eps_xx_0 = eps_xx_0;
output.eps_xx_z = eps_xx_z;
output.eps_xx_y = eps_xx_y;
output.gam_xz_0 = gam_xz_0;
output.gam_xz_y = gam_xz_y;
output.gam_xy_0 = gam_xy_0;
output.gam_xy_z = gam_xy_z;
end


function [valGP] = interpolateVal(valNode,N)
%This function interpolates a value using distinct values at valNode
%and the corresponding shape function N.
valGP = 0.0;
for i=1:length(N)
    valGP = valGP + N(i)*valNode(i);
end
end
