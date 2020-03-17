function  [dispOut,FReaction_sp1] = structuralDynamicsTransientROM(model,mesh,el,dispData,Omega,OmegaDot,~,delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData)
%structuralDynamicsTransientROM perform transient analysis with ROM
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [dispOut,FReaction_sp1] = structuralDynamicsTransientROM(model,mesh,...
%                             el,dispData,Omega,OmegaDot,time,delta_t,...
%                             elStorage,rom,Fexternal,Fdof,CN2H,rbData)
%
%   This function performs transient structural dynamics analysis using a
%   reduced order model (ROM).
%
%   input:
%   model      = object containing model data
%   mesh       = object containing mesh data
%   el         = object containing element data
%   dispData   = object containing displacement data
%   Omega      = rotor speed (Hz)
%   OmegaDot   = rotor acceleratin (Hz)
%   time       = current simulation time
%   delta_t    = time step size
%   elStorage  = object containing stored element data
%   rom        = object containing reduced order model represnetation
%   Fexternal  = vector containing external force values
%   Fdof       = vector containing global DOF numbering associated with
%                external force values
%   CN2H       = transformation matrix from inertial frame to hub frame
%   rbData     = vector containing rigid body displacement, velocity, and
%                acceleration
%
%   output:
%   dispOut       = object containing displacement data at end of time step
%   FReaction_sp1 = vector containing reaction force at turbine base at
%                   end of time step

%-------- get model information -----------
numEl = mesh.numEl;
x = mesh.x;
% y = mesh.y;
% z = mesh.z;
conn = mesh.conn;
numNodes = length(x);
elementOrder = model.elementOrder;
BC = model.BC;
nlROM = model.nlOn;

numNodesPerEl = elementOrder + 1;
numDOFPerNode = 6;
totalNumDOF = numNodes * numDOFPerNode;
% [~,numReducedDOF]=size(model.jointTransform);
% nodalTerms = model.nodalTerms;
%-----------------------------------------

%------- intitialization -----------------
%-------------------------------------------
analysisType = 'M';
%------ newmark integration parameters ---------
alpha = 0.5;
gamma = 0.5;
beta = 0.5*gamma;

timeInt.delta_t = delta_t;
timeInt.a1 = alpha*delta_t;
timeInt.a2 = (1.0-alpha)*delta_t;
timeInt.a3 = 1.0/(beta*delta_t*delta_t);
timeInt.a4 = timeInt.a3*delta_t;
timeInt.a5 = 1.0/gamma-1.0;
timeInt.a6 = alpha/(beta*delta_t);
timeInt.a7 = alpha/beta - 1.0;
timeInt.a8 = delta_t*(alpha/gamma-1.0);
disp_s = dispData.displ_s;
%     dispdot_s = dispData.displdot_s;
%     dispddot_s = dispData.displddot_s;
%-----------------------------------------------

%initialize displacements, tolerance, uNorm, iteration count for nonlinear
%iteration
displ_iter = disp_s;
displ_last = disp_s;
iterationCount = 1;
maxIter = 50;
uNorm = 1.0e6;
tol = 1.0e-6;


while(iterationCount < maxIter && uNorm > tol)  %iteration loop
    Kg = zeros(totalNumDOF);  %initialize global stiffness and force vector
    Fg = zeros(totalNumDOF,1);
    %% Element Assembly Loop for NL Terms
    if(nlROM)
        for i=1:numEl
            %Calculate Ke and Fe for element i
            index = 1;                           %initialize element data
            elInput.analysisType = analysisType;
            elInput.elementOrder = elementOrder;
            elInput.modalFlag = true;
            elInput.timeInt = timeInt;
            elInput.xloc = [0.0 el.elLen(i)];
            elInput.sectionProps = el.props(i);
            elInput.sweepAngle = el.psi(i);
            elInput.coneAngle = el.theta(i);
            elInput.rollAngle = el.roll(i);
            elInput.aeroSweepAngle = 0.0;
            
            elx = zeros(1,numNodesPerEl);
            eldispiter = zeros(1,numNodesPerEl*numDOFPerNode);
            for j=1:numNodesPerEl
                
                elx(j) = x(conn(i,j));
                
                %get element nodal displacements at s and s-1 time step
                for k=1:numDOFPerNode
                    eldispiter(index) = displ_iter((conn(i,j)-1)*numDOFPerNode + k);
                    index = index + 1;
                end
            end
            
            elInput.disp= eldispiter;
            elInput.x = elx;
            elInput.omegaVec = zeros(3,1); %convert from platform frame to hub-frame
            elInput.omegaDotVec = zeros(3,1);
            elInput.Omega = 0.0;
            elInput.OmegaDot = 0.0;
            elInput.CN2H = CN2H;
            elInput.preStress = false;
            elInput.useDisp = 1;
            elInput.iterationType = 'DI';  %uses direct iteration
            
            %do element calculations
            [elOutput] = calculateTimoshenkoElementNLSS(elInput); %calculate nonlinear timoshenko element stiffness matrix
            
            [Kg] = assemblyMatrixOnly(elOutput.Ke,conn(i,:),numNodesPerEl,numDOFPerNode,Kg); %assemble nonlinear timoshenko element stiffness matrix
            %................................................
        end
    end
    %------- end element calculation and assembly ------------------
    %%
    %Apply external loads to structure
    for i=1:length(Fexternal)
        Fg(Fdof(i)) = Fg(Fdof(i)) + Fexternal(i);
    end
    
    %------ apply constraints on system -----------------------------------
    %%
    [Fg] = applyConstraintsVec(Fg,model.jointTransform);
    
    if(nlROM)
        [Kg] = applyConstraints(Kg,model.jointTransform);
    end
    %----------------------------------------------------------------------
    %%
    %Apply BCs to global system
    [Fg]     = applyBCModalVec(Fg,BC.numpBC,BC.map);
    
    eta_s = dispData.eta_s;
    etadot_s = dispData.etadot_s;
    etaddot_s = dispData.etaddot_s;
    
    Phi       = rom.Phi;
    Feta = (Phi')*Fg;   %transform global displacement vector to modal space
    
    if(nlROM)
        [Kgnl]   = applyBCModal(Kg,BC.numpBC,BC.map);
        Kgnlrom = Phi'*Kgnl*Phi;     %transform nonlinear stiffness matrix to modal space
    else
        Kgnlrom = zeros(length(rom.Mr)); %if nonlinear deactivated set Kgnlrom to zero matrix
    end
    
    
    %define omega_i and omegaDot_i and body accelerations
    omega_platform_s = rbData(4:6);
    omega_x = omega_platform_s(1);
    omega_y = omega_platform_s(2);
    omega_z = omega_platform_s(3) + Omega*2*pi;
    
    omega_platform_dot = rbData(7:9);
    omegaDot_x = omega_platform_dot(1);
    omegaDot_y = omega_platform_dot(2);
    omegaDot_z = omega_platform_dot(3) + OmegaDot*2*pi;
    
    a_x = rbData(1); %platform accelerations (in hub frame)
    a_y = rbData(2);
    a_z = rbData(3);
    
    if(model.gravityOn)
        g= 9.81; %gravitational acceleration m/s^2
    else
        g = 0.0;
    end
    
    a_x_n = 0.0; %accelerations in inertial frame
    a_y_n = 0.0;
    a_z_n = g;
    a_temp = CN2H*[a_x_n; a_y_n; a_z_n];
    
    a_x = a_x + a_temp(1);
    a_y = a_y + a_temp(2);
    a_z = a_z + a_temp(3);
    
    %calculate reduced order spin soft, cent force, body force,
    %circulatory, coriolis
    
    Seff = rom.SrOx2.*omega_x^2 + rom.SrOy2.*omega_y^2 + rom.SrOz2.*omega_z^2+...
        + rom.SrOxOy.*omega_x*omega_y + rom.SrOyOz.*omega_y*omega_z...
        + rom.SrOxOz.*omega_x*omega_z;
    
    FcentEff = rom.FrOx2.*omega_x^2 + rom.FrOy2.*omega_y^2 + rom.FrOz2.*omega_z^2+...
        + rom.FrOxOy.*omega_x*omega_y + rom.FrOyOz.*omega_y*omega_z...
        + rom.FrOxOz.*omega_x*omega_z + rom.FrOxdot*omegaDot_x ...
        + rom.FrOydot*omegaDot_y + rom.FrOzdot*omegaDot_z;
    
    FbodyEff = rom.FrAx.*a_x + rom.FrAy.*a_y + rom.FrAz.*a_z;
    
    Heff = rom.HrOx.*omegaDot_x + rom.HrOy.*omegaDot_y + rom.HrOz.*omegaDot_z;
    
    Geff = rom.GrOx*omega_x + rom.GrOy*omega_y + rom.GrOz*omega_z;
    
    %combine for effective reduced order stiffness, damping, force
    Keff = rom.Kr + Seff + Heff + Kgnlrom;
    Ceff = rom.Cr + Geff;
    Feff = Feta + FcentEff + FbodyEff;
    
    [eta_iter,etadot_iter,etaddot_iter] = timeIntegrateSubSystemEff(rom.Mr,Keff,Ceff,Feff,timeInt,eta_s,etadot_s,etaddot_s);
    
    %reconstruct constrained dof vector from boundary conditions
    dispVec = Phi*eta_iter;
    dispdotVec = Phi*etadot_iter;
    dispddotVec = Phi*etaddot_iter;
    
    %[dispVec,dispdotVec,dispddotVec]    = constructReducedDispVector(dispVec,dispdotVec,dispddotVec,numNodes,BC);
    [dispVec]    = constructReducedDispVectorSingle(dispVec,BC.redVectorMap);
    
    jointTransform = model.jointTransform;
    displ_iter = jointTransform*dispVec;
    
    if(nlROM)
        uNorm = calculateUnorm(displ_iter,displ_last);
        displ_last = displ_iter;
        iterationCount = iterationCount + 1;
    else
        uNorm = 0.0;
    end
    
end
if(iterationCount>=maxIter)
    error('Maximum iterations exceeded.');
else
    
end
%------ calculate reaction at turbine base ----------------------------
reactionNodeNumber = model.platformTurbineConnectionNodeNumber;
[FReaction_sp1] = calculateReactionForceAtNode(reactionNodeNumber,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H);
%----------------------------------------------------------------------
%Calculate strain
[elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ_iter,model.nlOn);
dispOut.elStrain = elStrain;
%%
[dispdotVec]    = constructReducedDispVectorSingle(dispdotVec,BC.redVectorMap);  % reconstruct reduced velocity vector with boundary conditions
[dispddotVec]    = constructReducedDispVectorSingle(dispddotVec,BC.redVectorMap); % reconstruct reduced acceleration vector with boundary conditions
displdot_iter = (jointTransform*dispdotVec);   %construct unconstrained velocity vector
displddot_iter = (jointTransform*dispddotVec);  %construct unconstrained acceleration vector

dispOut.displ_sp1 = displ_iter;         %store physical and modal displacement, velocity, and acceleratoin in dispOut object
dispOut.displdot_sp1 = displdot_iter;
dispOut.displddot_sp1 = displddot_iter;

dispOut.eta_sp1 = eta_iter;
dispOut.etadot_sp1 = etadot_iter;
dispOut.etaddot_sp1 = etaddot_iter;

end

function [unorm] = calculateUnorm(unew,uold)
%This function calculates a relative norm between two vectors: unew and
%uold
unorm = norm(unew-uold)/norm(unew);

end

function [Kg] = applyConstraints(Kg,transMatrix)
%This function transforms a matrix by the transformation matrix to
%enforce joint constraints
transMatrix = sparse(transMatrix);
Kg = sparse(Kg);
Kg = transMatrix'*Kg*transMatrix;
Kg = full(Kg);
end

function [Fg] = applyConstraintsVec(Fg,transMatrix)
%This function transforms a vector by the transformation matrix to
%enforce joint constraints
transMatrix = sparse(transMatrix);
Fg = transMatrix'*Fg;
end

function [vec1Red] = constructReducedDispVectorSingle(vec1,redVectorMap)
%This function reconstructs displacement vector to account for boundary
%conditions
len = length(redVectorMap);
vec1Red = zeros(len,1);

for i=1:len
    if(redVectorMap(i) == -1.0)
        vec1Red(i) = 0.0;
    else
        vec1Red(i) = vec1(redVectorMap(i));
    end
end

end
