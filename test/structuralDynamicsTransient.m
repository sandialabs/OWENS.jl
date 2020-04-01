function  [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,dispData,Omega,OmegaDot,~,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)
%structuralDynamicsTransient perform transient analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,...
%                             dispData,Omega,OmegaDot,time,delta_t,...
%                             elStorage,Fexternal,Fdof,CN2H,rbData)
%
%   This function performs transient structural dynamics analysis.
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
y = mesh.y;
z = mesh.z;
conn = mesh.conn;
numNodes = length(x);
elementOrder = model.elementOrder;
BC = model.BC;

numNodesPerEl = elementOrder + 1;
numDOFPerNode = 6;
totalNumDOF = numNodes * numDOFPerNode;
% [~,numReducedDOF]=size(model.jointTransform);
nodalTerms = model.nodalTerms;
nodalTermsCopy = nodalTerms;
%-----------------------------------------

%initialize displacements, tolerance, uNorm, iteration count for nonlinear
%iteration
unorm = 1e6;
tol = model.nlParams.tolerance;
maxIterations = model.nlParams.maxIterations;
iterationCount = 0;

elx=zeros(numNodesPerEl,1);
ely=zeros(numNodesPerEl,1);
elz=zeros(numNodesPerEl,2);
eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
eldisp_sm1 = zeros(1,numNodesPerEl*numDOFPerNode);
eldispdot = eldisp;
eldispddot = eldisp;
eldispiter = eldisp;
% if(model.nlOn)
%      iterationType = model.nlParams.iterationType;
% else
%      iterationType = 'LINEAR';
% end
iterationType = 'DI';
analysisType = model.analysisType;
timeInt = struct('delta_t',0.0,...
'a1',0.0,...
'a2',0.0,...
'a3',0.0,...
'a4',0.0,...
'a5',0.0,...
'a6',0.0,...
'a7',0.0,...
'a8',0.0);
if(strcmp(analysisType,'TNB'))
    %------ newmark integration parameters ---------
    alpha = 0.5;
    gamma = 0.5;
    beta = 0.5*gamma;

    timeInt.delta_t = delta_t;
    timeInt.a1 = alpha*delta_t;
    timeInt.a2 = (1.0-alpha)*delta_t;
    a3 = 1.0/(beta*delta_t*delta_t);
    timeInt.a3 = a3;
    timeInt.a4 = a3*delta_t;
    timeInt.a5 = 1.0/gamma-1.0;
    timeInt.a6 = alpha/(beta*delta_t);
    timeInt.a7 = alpha/beta - 1.0;
    timeInt.a8 = delta_t*(alpha/gamma-1.0);

    disp_s = dispData.displ_s;
    dispdot_s = dispData.displdot_s;
    dispddot_s = dispData.displddot_s;

    displddot_im1 = dispddot_s;
    displdot_im1 = dispdot_s;
    displ_im1 = disp_s;
elseif(strcmp(analysisType,'TD'))
    %------ dean integration parameters -------------
    alpha = 0.25;

    timeInt.delta_t = delta_t;
    timeInt.a1 = alpha*delta_t^2;
    timeInt.a2 = (1-2*alpha)*delta_t^2;
    timeInt.a3 = delta_t/2.0;
    timeInt.a4 = delta_t*delta_t;
    disp_s = dispData.displ_s;
    %     disp_sm1 = dispData.displ_sm1;
    %-------------------------------------------
else
    error('analysis type not supported, choose another')
end
%-----------------------------------------------

while(unorm>tol && iterationCount < maxIterations) %iteration loop
    %------- intitialization -----------------
    Kg = zeros(totalNumDOF,totalNumDOF); %initialize global stiffness and force vector
    Fg = zeros(totalNumDOF,1);

    nodalTerms = nodalTermsCopy;
    %-------------------------------------------

    %---- element  calculation and assembly ----------------------------------
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
        if(iterationCount == 0)
            elInput.firstIteration = true;
        else
            elInput.firstIteration = false;
        end

        for j=1:numNodesPerEl

            %get element cooridnates
            elx(j) = x(conn(i,j));
            ely(j) = y(conn(i,j));
            elz(j) = z(conn(i,j));

            %get element nodal displacements at s and s-1 time step
            for k=1:numDOFPerNode
                %                 if(strcmp(analysisType,'TD'))
                %                     eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k);
                %                     eldisp_sm1(index) = disp_sm1((conn(i,j)-1)*numDOFPerNode + k);
                %                     eldispiter(index) = displ_iter((conn(i,j)-1)*numDOFPerNode + k);
                %                 end
                if(strcmp(analysisType,'TNB'))
                    eldispiter(index) = displ_im1((conn(i,j)-1)*numDOFPerNode + k);
                    if(strcmp(iterationType,'NR'))
                        eldisp(index) = displ_im1((conn(i,j)-1)*numDOFPerNode + k);
                        eldispdot(index) = displdot_im1((conn(i,j)-1)*numDOFPerNode + k);
                        eldispddot(index) = displddot_im1((conn(i,j)-1)*numDOFPerNode + k);
                    elseif(strcmp(iterationType,'DI')||strcmp(iterationType,'LINEAR'))
                        eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k);
                        eldispdot(index) = dispdot_s((conn(i,j)-1)*numDOFPerNode + k);
                        eldispddot(index) = dispddot_s((conn(i,j)-1)*numDOFPerNode + k);
                    end
                end
                index = index + 1;
            end
        end

        %get concentrated terms associated with elemetn
        [massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(conn(i,:),model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad);

        elInput.concMass = massConc;
        elInput.concStiff = stiffConc;
        elInput.concLoad = loadConc;
        elInput.disp = eldisp;

        % specific to 'TD', but must be declared
        elInput.dispm1= eldisp_sm1;

        % specific to 'TNB' , but must be declared
        elInput.dispdot = eldispdot;
        elInput.dispddot = eldispddot;

        elInput.x = elx;
        elInput.y = ely;
        elInput.z = elz;
        elInput.accelVec = rbData(1:3);
        elInput.Omega = Omega;
        elInput.OmegaDot = OmegaDot;
        elInput.omegaVec = rbData(4:6);
        elInput.omegaDotVec = rbData(7:9);
        if(el.rotationalEffects(i))
            elInput.Omega = Omega;
            elInput.OmegaDot = OmegaDot;
        else
            elInput.Omega = 0.0;
            elInput.OmegaDot = 0.0;
        end

        elInput.displ_iter = eldispiter;
        elInput.useDisp = model.nlOn;
        elInput.preStress = false;
        elInput.iterationType = iterationType;
        elInput.freq = 0.0; %Is not used for this model type, but must be declared.
        elInput.aeroElasticOn = false;
        elInput.aeroForceOn = false;
        elInput.airDensity = model.airDensity;
        elInput.gravityOn = model.gravityOn;

        elInput.RayleighAlpha = model.RayleighAlpha;
        elInput.RayleighBeta = model.RayleighBeta;

        elInput.CN2H = CN2H;

        [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage(i)); %calculate timoshenko element

        [Kg,Fg] = assembly(elOutput.Ke,elOutput.Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); %assemble element stiffness matrix and force vector

        %         Erestotal = Erestotal + elOutput.Eres;
        %................................................
    end
    %------- end element calculation and assembly ------------------

    %%
    %----------------------------------------------------------------------

    %%
    %Apply external loads to structure
    for i=1:length(Fexternal)
        if(strcmp(analysisType,'TD'))
            Fg(Fdof(i)) = Fg(Fdof(i)) + Fexternal(i)*delta_t^2;
        end
        if(strcmp(analysisType,'TNB'))
            Fg(Fdof(i)) = Fg(Fdof(i)) + Fexternal(i);
        end
    end

    %------ apply constraints on system -----------------------------------
    [Kg] = applyConstraints(Kg,model.jointTransform);
    [Fg] = applyConstraintsVec(Fg,model.jointTransform);

    %----------------------------------------------------------------------
    %%

    %Apply BCs to global system
    [KgTotal,FgTotal] = applyBC(Kg,Fg,BC,numDOFPerNode);

    solution = KgTotal\FgTotal;  %solve for displacements
    solution = model.jointTransform*solution; %transform to full dof listing

    if(model.nlOn)  %calculate norm between current iteration and last iteration
        if(strcmp(iterationType,'NR'))
            [unorm] = calcUnorm(displ_im1+solution,displ_im1);
        else
            [unorm] = calcUnorm(solution,displ_im1);
        end
    else
        unorm = 0.0;
    end

    if(strcmp(iterationType,'NR'))
        %if newton raphson update u, udot, uddot at each iteration
        displ_im1 = displ_im1 + solution;
        cap_delta_displ = displ_im1 - dispData.displ_s;
        displddot_im1 = timeInt.a3*(cap_delta_displ) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s;
        displdot_im1  = -timeInt.a7*dispData.displdot_s -timeInt.a8*dispData.displddot_s + timeInt.a6*(cap_delta_displ);
    elseif(strcmp(iterationType,'DI')||strcmp(iterationType,'LINEAR'))
        displ_im1 = solution;
    else
        error('iteration type not supported, choose another')
    end

    iterationCount = iterationCount + 1;
end
%Calculate reaction at turbine base (hardwired to node number 1)
reactionNodeNumber = model.platformTurbineConnectionNodeNumber;
[FReaction] = calculateReactionForceAtNode(reactionNodeNumber,model,...
    mesh,el,elStorage,timeInt,dispData,displ_im1,rbData,Omega,OmegaDot,CN2H);

%Calculate strain
[elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ_im1,model.nlOn);
dispOut.elStrain = elStrain;
if(iterationCount>=maxIterations)
    error('Maximum iterations exceeded.');
end

FReaction_sp1 =FReaction;
displ_sp1 = displ_im1;
dispOut.displ_sp1 = displ_sp1;  %store displacement vector in dispOut

% Specific to TNB, but must be declared
displddot_sp1 = timeInt.a3*(displ_sp1-dispData.displ_s) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s; %store velocity vector in dispOut
dispOut.displddot_sp1 = displddot_sp1;
dispOut.displdot_sp1 = dispData.displdot_s + timeInt.a2*dispData.displddot_s + timeInt.a1*displddot_sp1;                    %store acceleration vector in dispOut



end


function [Kg] = applyConstraints(Kg,transMatrix)
%This function transforms a matrix by the transformation matrix to
%enforce joint constraints
Kg = transMatrix'*Kg*transMatrix;
end

function [Fg] = applyConstraintsVec(Fg,transMatrix)
%This function transforms a vector by the transformation matrix to
%enforce joint constraints
Fg = transMatrix'*Fg;
end

function [unorm] = calcUnorm(unew,uold)
%This function calculates a relative norm between two vectors: unew and
%uold
unorm = norm(unew-uold)/norm(unew);
end
