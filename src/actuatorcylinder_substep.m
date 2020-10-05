
function [Rp, Tp, Zp, Mp,env] = actuatorcylinder_substep(turbines, env,us_param,step_AC,alpha_rad,cl_af,cd_af)
ntheta = turbines.ntheta;
% list comprehensions
centerX = zeros(1,length(turbines));
centerY = zeros(1,length(turbines));
radii = zeros(1,length(turbines));
for i = 1:length(turbines)
    centerX(i) = turbines(i).centerX;
    centerY(i) = turbines(i).centerY;
    radii(i) = turbines(i).R;
end

% assemble global matrices
[Ax, Ay, theta] = matrixAssemble(centerX, centerY, radii, ntheta);

% setup

wsave = us_param.awsave;
circ_step_num = floor((step_AC-1)/ntheta*turbines.B);
circular_step = step_AC-circ_step_num*ntheta/turbines.B; %env.circular_step;
V_wake_old = us_param.V_wake_old;
Vinf_nominal = env.Vinf_nominal;
tau = us_param.tau;


Vinf_used = zeros(1,turbines.B);
gustT = us_param.gusttime * Vinf_nominal / turbines.R;
dt = 1/(mean(turbines.omega) * ntheta / (2*pi));
dt_norm = dt*Vinf_nominal/turbines.R;
% IECGustFactor = 0;
%TODO: add check that ntheta is divisible by nblades and 2
% for step = 1:ntheta*env.N_Rev %
idx_sub = circular_step:ntheta/turbines.B:ntheta*2-ntheta/turbines.B+1+circular_step;
w0_sub = wsave(idx_sub);
idx = (i-1)*ntheta+1:i*ntheta;
% if mod(step-1,ntheta/turbines.B) == 0
%     circular_step = 1;
% end

ele_x = sin(double((step_AC:ntheta/turbines.B:step_AC+ntheta-ntheta/turbines.B+1)/ntheta*2*pi));
tr = double(step_AC)*double(dt_norm) - ele_x - double(us_param.gustX0);

for bld_i = 1:turbines.B
    if (tr(bld_i) >= 0) && (tr(bld_i)<=gustT)

        IECGustFactor = 1.0 - 0.37 * us_param.G_amp/Vinf_nominal * sin(3*pi*tr(bld_i)/gustT)  * (1.0 - cos(2*pi*tr(bld_i)/gustT));
        Vinf_used(bld_i) = Vinf_nominal*IECGustFactor;

    else
        Vinf_used(bld_i) = Vinf_nominal;
    end
end


env.Vinf(idx_sub(1:turbines.B)) = Vinf_used;
f_residual = @(wsub) frozen_residual(wsub,wsave,idx_sub, [Ax(idx, idx); Ay(idx, idx)], theta, 1.0, (turbines), env,alpha_rad,cl_af,cd_af);

options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
% tic
[w,~,info,~] = fsolve(f_residual,w0_sub,options);
%     toc

wnew = wsave;

wnew(idx_sub) = w;
u = wnew(idx);
v = wnew(ntheta + idx);
[~, ~, ~, ~, ~, ~, ~, ~, a_new] = radialforce(u, v, theta, turbines(i), env,alpha_rad,cl_af,cd_af);
%      constants, radius, average wake velocity

tau_near = tau(1)*turbines.R/V_wake_old;
tau_far = tau(2)*turbines.R/V_wake_old;

V_wake_old = V_wake_old*exp(double(-dt/tau_far))+(mean(Vinf_used)*(1-2*a_new))*(1-exp(double(-dt/tau_far)));

w_filtered = wsave*exp(double(-dt/tau_near))+wnew*(1-exp(double(-dt/tau_near)));

wsave = w_filtered;

% if info ~=1
%     fprintf(['fsolve terminated prematurely. info = ' num2str(info)]);
% end

idx = 1:ntheta;
u = w_filtered(idx);
v = w_filtered(ntheta + idx);
[~, ~, ~, ~, Rp_temp, Tp_temp, Zp_temp, ~, ~] = radialforce(u, v, theta, turbines(i), env,alpha_rad,cl_af,cd_af);

idx_sub_blade = circshift(idx_sub(1:turbines.B),circ_step_num);
Rp = Rp_temp(idx_sub_blade(1:turbines.B))';
Tp = Tp_temp(idx_sub_blade(1:turbines.B))';
Zp = Zp_temp(idx_sub_blade(1:turbines.B))';
Mp = zeros(1,length(Zp)); %TODO


% circular_step = circular_step + 1;

% env.circular_step = circular_step;
if env.steplast ~= step_AC
    us_param.awsave = wsave;
    env.V_wake_old = V_wake_old;
    env.steplast = step_AC;
end
env.idx_sub = idx_sub;

% end

% CP_cactus = [0.7126564;0.5164276;0.4486862;0.4202281;0.4053651;0.3960852;0.3898851;0.3852592;0.3817074;0.3790016;0.3860673;0.4094897;0.4422996;0.3769674;0.3803902;0.3773890;0.3753646;0.3737490;0.3725482;0.3717143];
%
% cactus_rev = linspace(1,length(CP_cactus),ntheta*19);
% interp_CP_cactus = interp1((1:length(CP_cactus)),CP_cactus,cactus_rev,'makima');
% useRev = 10;
% resid = sum(sqrt((interp_CP_cactus(1:ntheta*(useRev-1))-(CP(ntheta+1:ntheta*useRev))+tau(3)).^2));
%
% % residual = real(residual)+abs(imag(residual)^5);
%
% close all
% figure()
% scatter(step,Tp(1))
% hold on
% pause(0.001);
% plot((1:length(CP_cactus))*ntheta,CP_cactus,'k.-')
% hold on
% plot((1:length(cactus_rev))+ntheta,interp_CP_cactus,'k.-')
% pause(0.001);
% disp(tau);
% disp(resid);
end

% --- Influence Coefficients ---

% """
% applies for both Ay and Rx depending on which function ifunc(x, y, phi)
% is passed in
% """
function A = panelIntegration(xvec, yvec, thetavec, ifunc)

% initialize
nx = length(xvec);
ntheta = length(thetavec);
dtheta = thetavec(2) - thetavec(1);  % assumes equally spaced
A = zeros(nx, ntheta);

for i = 1:length(xvec)
    % redefine function so it has one parameter for use in quadgk
    if ifunc == 'Ayintegrand'
        integrand = @(phi) Ayintegrand(xvec(i), yvec(i), phi);
    else
        integrand = @(phi) Dxintegrand(xvec(i), yvec(i), phi);
    end

    for j = 1:length(thetavec)
        % an Adaptive Gauss-Kronrod quadrature integration.  Tried trapz but this was faster.

        [A(i, j), ~] = quadgk(integrand, thetavec(j)-dtheta/2.0, thetavec(j)+dtheta/2.0,'AbsTol',1e-12, 'RelTol', 1e-8, 'MaxIntervalCount',5000);
        %         A(i, j) = integral(integrand, thetavec(j)-dtheta/2.0, thetavec(j)+dtheta/2.0,'AbsTol',1e-10, 'RelTol', 1e-6);
        % x = linspace(thetavec(j)-dtheta/2.0, thetavec(j)+dtheta/2.0,10000);
        % A(i,j) = trapz(x,integrand); % integrate y w.r.t. x
    end

end

end


% """
% integrand used for computing Dx
% """
function output = Dxintegrand(x, y, phi)
v1 = x + sin(phi);
v2 = y - cos(phi);
% v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
output = (v1.*sin(phi) - v2.*cos(phi))./(2.*pi.*(v1.*v1 + v2.*v2));
end


% """
% integrand used for computing Ay
% """
function output = Ayintegrand(x, y, phi)
output = zeros(1,length(phi));

v1 = x + sin(phi);
v2 = y - cos(phi);
% occurs when integrating self, function symmetric around singularity, should integrate to zero
set_zero = and(abs(v1) < 1e-12 , abs(v2) < 1e-12);
output(set_zero) = 0.0;

output(~set_zero) = (v1.*cos(phi) + v2.*sin(phi))./(2*pi.*(v1.*v1 + v2.*v2));

end


function output = AyIJ(xvec, yvec, thetavec)
output = panelIntegration(xvec, yvec, thetavec, 'Ayintegrand');
end

function output = DxIJ(xvec, yvec, thetavec)
output = panelIntegration(xvec, yvec, thetavec, 'Dxintegrand');
end

function Wx = WxIJ(xvec, yvec, thetavec) %Used

% initialize
nx = length(xvec);
ntheta = length(thetavec);
dtheta = thetavec(2) - thetavec(1);  % assumes equally spaced
Wx = zeros(nx, ntheta);

for i = 1:length(xvec)
    if yvec(i) >= -1.0 && yvec(i) <= 1.0 && xvec(i) >= 0.0 && xvec(i)^2 + yvec(i)^2 >= 1.0
        % if yvec(i) >= -1.0 && yvec(i) <= 1.0 && (xvec(i) >= 0.0 || (xvec(i) >= -1 && xvec(i)^2 + yvec(i)^2 <= 1.0));
        thetak = acos(yvec(i));
        k = find((thetavec + dtheta/2) > thetak,1,'first');  % index of intersection
        Wx(i, k) = -1.0;
        Wx(i, ntheta-k+1) = 1.0;
    end
end
end

function Rx = DxII(thetavec) %Used

% initialize
ntheta = length(thetavec);
dtheta = thetavec(2) - thetavec(1);  % assumes equally spaced
Rx = dtheta/(4*pi)*ones(ntheta, ntheta);

for i = 1:length(thetavec)
    if i <= ntheta/2
        Rx(i, i) = (-1 + 1.0/ntheta)/2.0;
    else
        Rx(i, i) = (1 + 1.0/ntheta)/2.0;
    end
end

end

function Wx = WxII(thetavec) %Used

% initialize
ntheta = length(thetavec);
Wx = zeros(ntheta, ntheta);

for i = fix(ntheta/2)+1:ntheta
    Wx(i, ntheta+1-i) = -1;
end

end

function [theta,Dxself,Wxself,Ayself] = precomputeMatrices(ntheta) %Used

% precompute self influence matrices

% setup discretization (all the same, and uniformly spaced in theta)
dtheta = 2*pi/ntheta;
theta = dtheta/2:dtheta:2*pi;

Dxself = DxII(theta);
Wxself = WxII(theta);
Ayself = AyIJ(-sin(theta), cos(theta), theta);

save(['theta-' num2str(ntheta) '.mat'],'theta','Dxself','Wxself','Ayself')

end


function [Ax, Ay, theta] = matrixAssemble(centerX, centerY, radii, ntheta) %Used
%     """
%     centerX, centerY: array of x,y coordinates for centers of the VAWTs in the farm
%     radii: corresponding array of their radii
%     """


if ~isfile(['theta-' num2str(ntheta) '.mat'])
    [theta,Dxself,Wxself,Ayself] = precomputeMatrices(ntheta);
else
    precompMat = load(['theta-' num2str(ntheta) '.mat']);
    theta = precompMat.theta;
    Dxself = precompMat.Dxself;
    Wxself = precompMat.Wxself;
    Ayself = precompMat.Ayself;
end

% initialize global matrices
nturbines = length(radii);
Dx = zeros(nturbines*ntheta, nturbines*ntheta);
Wx = zeros(nturbines*ntheta, nturbines*ntheta);
Ay = zeros(nturbines*ntheta, nturbines*ntheta);

% iterate through turbines
for I = 1:length(radii)
    for J = 1:length(radii)

        % find normalized i locations relative to center of turbine J
        x = (centerX(I)-radii(I)*sin(theta) - centerX(J))/radii(J);
        y = (centerY(I)+radii(I)*cos(theta) - centerY(J))/radii(J);

        % self-influence is precomputed
        if I == J
            Dxsub = Dxself;
            Wxsub = Wxself;
            Aysub = Ayself;

            % pairs can be mapped for same radius
        elseif J < I && radii(I) == radii(J)

            % grab cross-diagonal I,J -> J,I matrix
            Dxsub = Dx((J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta);
            Aysub = Ay((J-1)*ntheta+1:J*ntheta, (I-1)*ntheta+1:I*ntheta);

            % mapping index for coefficients that are the same
            idx = [fix(ntheta/2)+1:ntheta, 1:fix(ntheta/2)];

            % directly map over
            Dxsub = Dxsub(idx, idx);
            Aysub = Aysub(idx, idx);

            % wake term must be recomptued
            Wxsub = WxIJ(x, y, theta);

            % % if VAWTs are very far apart we can approximate some of the influence coefficients
            % elseif approxfar && sqrt((centerX(I)-centerX(J))^2 + (centerY(I)-centerY(J))^2) > 10*radii(I)
            %     println("far apart")
            %     xc = (centerX(I) - centerX(J))/radii(J);
            %     yc = (centerY(I) - centerY(J))/radii(J);

            %     Rxsub = RxIJFar(xc, yc, theta);
            %     Wxsub = zeros(ntheta, ntheta)  % should have negligible wake impact;
            %     Aysub = AyIJFar(xc, yc, theta);

        else
            Dxsub = DxIJ(x, y, theta);
            Wxsub = WxIJ(x, y, theta);
            Aysub = AyIJ(x, y, theta);
        end

        % assemble into global matrix
        Dx((I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta) = Dxsub;
        Wx((I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta) = Wxsub;
        Ay((I-1)*ntheta+1:I*ntheta, (J-1)*ntheta+1:J*ntheta) = Aysub;

    end
end
Ax = Dx + Wx;

end
% ---------------------------


% ------- Force Coefficients ---------

function [q, ka, CT, CP, Rp, Tp, Zp, Qrev, a] = radialforce(uvec, vvec, thetavec, turbine, env,alpha_rad,cl_af,cd_af) %Used
% u, v, theta - arrays of size ntheta
% r, chord, twist, Vinf, Omega, rho, mu - scalars

% unpack
r = turbine.r;
chord = turbine.chord;
twist = turbine.twist;
delta = turbine.delta;
B = double(turbine.B);
Omega = turbine.omega;

Vinf = env.Vinf;
rho = env.rho;
%Velocities due to rotor deformation
V_vert = env.V_vert; %TODO: include vertical velocity effect
V_tang = env.V_tang;
V_rad = env.V_rad;
V_twist = env.V_twist; %TODO: include twist velocity effect along with dynamic stall model
% set the rotation direction
rotation = sign(Omega(1));

% velocity components and angles

Vn = V_rad + Vinf.*(1.0 + uvec).*sin(thetavec) - Vinf.*vvec.*cos(thetavec);
Vt = V_tang + rotation*(Vinf.*(1.0 + uvec).*cos(thetavec) + Vinf.*vvec.*sin(thetavec)) + abs(Omega).*r;
W = sqrt(Vn.^2 + Vt.^2);
phi = atan2(Vn, Vt);
alpha = phi - twist;
% Re = rho*W*chord/mu  % currently no Re dependence;

% airfoil
% cl = turbine.af.cl(alpha);
% cd = turbine.af.cd(alpha);
%     cl, cd = turbine.af(alpha);
cl = interp1(alpha_rad, cl_af, alpha,'makima');
cd = interp1(alpha_rad, cd_af, alpha,'makima');

% rotate force coefficients
cn = cl.*cos(phi) + cd.*sin(phi);
ct = cl.*sin(phi) - cd.*cos(phi);

% radial force
sigma = B*chord./r;
q = sigma./(4*pi).*cn.*(W./Vinf).^2;

% instantaneous forces
qdyn = 0.5*rho*W.^2;
Rp = -cn.*qdyn*chord;
Tp = ct.*qdyn.*chord./cos(delta);
Zp = -cn.*qdyn.*chord.*tan(delta);

% nonlinear correction factor
integrand = sigma./(4*pi).*(W./Vinf).^2 .* (cn.*sin(thetavec) - rotation*ct.*cos(thetavec)./cos(delta));
CT = pInt(thetavec, integrand);
if CT > 2.0
    a = 0.5*(1.0 + sqrt(1.0 + CT));
    ka = 1.0 / (a-1);

elseif CT > 0.96
    a = 1.0/7*(1 + 3.0*sqrt(7.0/2*CT - 3));
    ka = 18.0*a / (7*a^2 - 2*a + 4);

else
    a = 0.5*(1 - sqrt(1.0 - CT));
    ka = 1.0 / (1-a);
end

% power coefficient
H = 1.0;  % per unit height
Sref = 2*turbine.R*H;
Q = r.*Tp;
Qrev = B*pInt(thetavec, Q)/(2*pi);
P = Qrev * mean(abs(Omega));
CP = P / (0.5*rho*mean(Vinf)^3 * Sref);

end

% -----------------------------------------



% ------ Solve System --------------

function output = residual(w, A, theta, k, turbines, env) %Used

% setup
ntheta = length(theta);
nturbines = length(turbines);  %  int(length(w)/2/ntheta)
q = zeros(1,ntheta*nturbines);
ka = 0.0;

for i = 1:length(turbines)
    idx = (i-1)*ntheta+1:i*ntheta;

    u = w(idx);
    v = w(ntheta*nturbines + idx);

    [q(idx), ka, ~, ~, ~, ~, ~, ~]= radialforce(u, v, theta, turbines(i), env,alpha_rad,cl_af,cd_af);
end

if nturbines == 1  % if only one turbine use the k from the analysis;
    k = (ka);
end  % otherwise, use k that was input to this function

% reformat to multiply in correct locations
kmult = repmat(k, ntheta,1);
kmult = reshape(kmult,[ntheta*nturbines,1]);
kmult = [kmult; kmult];

output = (A*q').*kmult - w';
end

function output = frozen_residual(wsub,w,idx_sub, A, theta, k, turbines, env,alpha_rad,cl_af,cd_af) %Used
w(idx_sub) = wsub; %Reconstruct the full induced velocities matrix
% setup
ntheta = length(theta);
nturbines = length(turbines);  %  int(length(w)/2/ntheta)
q = zeros(1,ntheta*nturbines);
ka = 0.0;

for i = 1:length(turbines)
    idx = (i-1)*ntheta+1:i*ntheta;

    u = w(idx);
    v = w(ntheta*nturbines + idx);

    [q, ka, ~, ~, ~, ~, ~, ~]= radialforce(u, v, theta, turbines(i), env,alpha_rad,cl_af,cd_af);

end

if nturbines == 1  % if only one turbine use the k from the analysis;
    k = (ka);
end  % otherwise, use k that was input to this function

% reformat to multiply in correct locations
kmult = repmat(k, ntheta,1);
kmult = reshape(kmult,[ntheta*nturbines,1]);
kmult = [kmult; kmult];

disp(size(A))
disp(size(q))
disp(size(kmult))
disp(size(w))
output = (A*q(1,:)').*kmult - w';
end




% -----------------------------------------

% ---------- helper methods --------------

% trapezoidal integration
function integral = trapz(x, y)  % integrate y w.r.t. x %Used

integral = 0.0;
for i = 1:length(x)-1
    integral = integral + (x(i+1)-x(i))*0.5*(y(i) + y(i+1));
end

end

% integration for a periodic function where end points dont reach ends (uses trapezoidal method)
function integral = pInt(theta, f) %Used

integral = trapz(theta, f);

% add end points
dtheta = 2*theta(1);  % assumes equally spaced, starts at 0
integral = integral + dtheta * 0.5*(f(1) + f(end));

end
