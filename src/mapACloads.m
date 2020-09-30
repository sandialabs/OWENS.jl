function [Fg,ForceDof,env] = mapACloads(u_jLast,udot_j,Omega_j,t,PEy,QCy,NElem,NBlade,RefR,mesh,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,el,turbine3D,env,step_AC)

N_blade_nodes = length(structuralSpanLocNorm)+1;
% Initialize bladeForce
init_bladeForce = struct('N',zeros(1,NElem),'T',zeros(1,NElem),'M25',zeros(1,NElem));
bladeForce = repmat(init_bladeForce,1,NBlade);

for k = 1:length(turbine3D)
    %TODO: Incorporate deflections and changes in omega ->
    %r,twist,delta,omega all need to be vectors aligning with the ntheta
    %discretizations of the cylinder

    %TODO: ensure that the deflections aren't compounding after a
    %revolution.  They may be wrong

    %TODO: Verify units everywhere

    circ_step_num = floor((step_AC-1)/env{k}.ntheta*turbine3D{k}.B);
    circular_step = step_AC-circ_step_num*env{k}.ntheta/turbine3D{k}.B;
    idx_sub = circular_step:env{k}.ntheta/turbine3D{k}.B:env{k}.ntheta-env{k}.ntheta/turbine3D{k}.B+1+circular_step;

    %TODO: this is hard coded for 2 blades, need to simplify
    % Interpolate the deformations onto the aero model for the current step
    norm_disp_h = linspace(0,1,N_blade_nodes)';
    % 1 = Z deformation - not modeled in 2D AC method
    % 2 = Tangential deformation - no real effect on AC model
    offset = 3;
    turbine3D{k}.r(idx_sub(1)) = turbine3D{k}.r(idx_sub(1)) + interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    turbine3D{k}.r(idx_sub(2)) = turbine3D{k}.r(idx_sub(2)) + interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    offset = 4;
    turbine3D{k}.twist(idx_sub(1)) = turbine3D{k}.twist(idx_sub(1)) + interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    turbine3D{k}.twist(idx_sub(2)) = turbine3D{k}.twist(idx_sub(2)) + interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    offset = 5;
    turbine3D{k}.delta(idx_sub(1)) = turbine3D{k}.delta(idx_sub(1)) + interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    turbine3D{k}.delta(idx_sub(2)) = turbine3D{k}.delta(idx_sub(2)) + interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    % 6 = Sweep deformation, not modeled in AC method - assuming it is small so that it doesn't spill over into the next step/theta discretization

    turbine3D{k}.omega = Omega_j;

    % Interpolate deformation induced velocities onto the aero model for the most current step
    offset = 1;
    env{k}.V_vert(idx_sub(1)) = interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    env{k}.V_vert(idx_sub(2)) = interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    offset = 2;
    env{k}.V_tang(idx_sub(1)) = interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    env{k}.V_tang(idx_sub(2)) = interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    offset = 3;
    env{k}.V_rad(idx_sub(1)) = interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    env{k}.V_rad(idx_sub(2)) = interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    offset = 4;
    env{k}.V_twist(idx_sub(1)) = interp1(norm_disp_h,u_jLast(N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset),k/length(turbine3D));
    env{k}.V_twist(idx_sub(2)) = interp1(norm_disp_h,u_jLast(N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset),k/length(turbine3D));
    % 5 = Change in delta angle, not modeled in 2D AC method - second order
    % 6 = Change in sweep angle, not modeled in 2D AC method

    [Rp, Tp, Zp, Mp,env{k}] = actuatorcylinder_substep(turbine3D{k}, env{k}, step_AC);
    for j=1:NBlade
        bladeForce(j).N(k) = Rp(j);
        bladeForce(j).T(k) = -Tp(j);
        bladeForce(j).M25(k) = Mp(j);
    end
end



% scatter(t,bladeForce(1).T(floor(blade(j).NElem/2)));
% hold on
% pause(0.001)
%define these from params file
ft2m = 1 / 3.281;

%     RefAR = cactusGeom.RefAR*ft2m*ft2m;
RefR = RefR*ft2m;

spanLocNorm = zeros(NBlade,NElem);
for i=1:NBlade
    spanLocNorm(i,:) = PEy(1:NElem(1,1),1).*RefR(1,1)/(QCy(NElem(1,1)+1,1)*RefR(1,1));
end

%Initialize structuralLoad
init_structuralLoad = struct('N',zeros(1,length(structuralElNumbers)),'T',zeros(1,length(structuralElNumbers)),'M25',zeros(1,length(structuralElNumbers)));
structuralLoad = repmat(init_structuralLoad,1,NBlade);

for i=1:NBlade
    structuralLoad(i).N = linear_interp(spanLocNorm(i,:),bladeForce(i).N,structuralSpanLocNorm(i,:));
    structuralLoad(i).T = linear_interp(spanLocNorm(i,:),bladeForce(i).T,structuralSpanLocNorm(i,:));
    structuralLoad(i).M25 = linear_interp(spanLocNorm(i,:),bladeForce(i).M25,structuralSpanLocNorm(i,:));
end

[~,numNodesPerBlade] = size(structuralNodeNumbers);

%integrate over elements

%read element data in

numDofPerNode = 6;
%     [~,~,timeLen] = size(aeroDistLoadsArrayTime);
Fg = zeros(max(max(structuralNodeNumbers))*6,1);
for j = 1:NBlade
    for k = 1:numNodesPerBlade-1
        %get element data
        % orientation angle,xloc,sectionProps,element order]
        elNum = structuralElNumbers(j,k);
        %get dof map
        node1 = structuralNodeNumbers(j,k);
        node2 = structuralNodeNumbers(j,k+1);
        dofList = [(node1-1)*numDofPerNode+(1:6), (node2-1)*numDofPerNode+(1:6)];

        elInput.elementOrder = 1;
        elInput.x = [mesh.x(node1), mesh.x(node2)];
        elLength = sqrt((mesh.x(node2)-mesh.x(node1))^2 + (mesh.y(node2)-mesh.y(node1))^2 + (mesh.z(node2)-mesh.z(node1))^2);
        elInput.xloc = [0 elLength];
        elInput.sectionProps.twist = el.props{elNum}.twist;
        elInput.sweepAngle = el.psi(elNum);
        elInput.coneAngle = el.theta(elNum);
        elInput.rollAngle = el.roll(elNum);

        elInput.extDistF2Node =  [structuralLoad(j).T(k),   structuralLoad(j).T(k+1)];
        elInput.extDistF3Node = -[structuralLoad(j).N(k),   structuralLoad(j).N(k+1)];
        elInput.extDistF4Node = -[structuralLoad(j).M25(k), structuralLoad(j).M25(k+1)];

        [output] = calculateLoadVecFromDistForce(elInput);
        Fe = output.Fe;

        %asssembly
        for m = 1:length(dofList)
            Fg(dofList(m)) =  Fg(dofList(m))+Fe(m);
        end

    end
end

ForceDof = 1:length(Fg);

end
