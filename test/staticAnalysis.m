function [displ,elStrain,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,Omega,~,elStorage)
%staticAnalysis performs static analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [displ,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,...
%                                    Omega,OmegaStart,elStorag
%
%   This function performs a static analysis and returns displacement
%   values and a flag denoting successful/unsuccessful analysis
%
%   input:
%   model          = object containing model information
%   mesh           = object containing mesh information
%   el             = object containing element information
%   displ          = displacement vector for use in pre-stressed analysis
%   Omega          = rotor speed (Hz)
%   OmegaStart     = rotor speed (Hz) from previous analysis if stepping
%                    through various rotor speeds, may be useful in load
%                    stepping
%   elStorage      = previously calculated element system matrices
%
%
%   output:
%   displ                    = vector of displacemetns
%   staticAnalysisSuccessful = boolean flag denoting successful static
%                              analysis


tic

numEl = mesh.numEl; %extract mesh information from mesh object
x = mesh.x;
y = mesh.y;
z = mesh.z;
conn = mesh.conn;
numNodes = length(x);

elementOrder = model.elementOrder; %extract element order from model

BC = model.BC; %extract boundary  condition object from model

numNodesPerEl = elementOrder + 1; %do initialization
numDOFPerNode = 6;
totalNumDOF = numNodes * numDOFPerNode;

nodalTerms = model.nodalTerms; %extract concentrated nodal terms from model
nodalTermsCopy = nodalTerms;   %extract extra copy of concentrated nodal terms from model

%load stepping paramters
loadStepCount = 1;
displPrev = displ; %copy of initial displacement vector
staticAnalysisComplete = false; %initialize to false
nlParams = model.nlParams;

if(nlParams.adaptiveLoadSteppingFlag)
    loadStepPrev = 0.0;
    loadStep = 1.0;
else
    loadStepPrev = 0.0;
    loadStep = nlParams.prescribedLoadStep(1);
    fprintf('Prescribed load step: %f\n',loadStep);
end

maxNumLoadSteps = nlParams.maxNumLoadSteps;
MAXIT = nlParams.maxIterations;
tolerance = nlParams.tolerance;
dispOld = displ;  %initialize dispOld, first iteration logic below for accurate calculations
%.........................................................................
staticAnalysisSuccessfulForLoadStep = false; %initialize variable
while(~staticAnalysisComplete && loadStepCount<maxNumLoadSteps)
    % staticAnalysisSuccessful = false; %initialize staticAnalysisSuccessful flag
    uNorm = 1.0; %initialize norm for convergence check
    
    iterationCount = 0; %initialize iteration count
    
    while(uNorm > tolerance && iterationCount < MAXIT) %iteration loop (convergence tolerance of 1.0e-6)
        Kg = zeros(totalNumDOF,totalNumDOF);   %initialize global stiffness matrix
        Fg = zeros(totalNumDOF,1);             %initialize global force vector
        nodalTerms = nodalTermsCopy;
        for i=1:numEl
            if(iterationCount >=1)
                gamm = 0.0;     %option for acceleration of iterative procedure (gamma = 0 or gamma=0.5 are typical)
                dispEval = dispOld.*gamm + displ.*(1-gamm);
            else
                dispEval = displ;
            end
            %Calculate Ke and Fe for element i
            index = 1;
            elInput.analysisType = 'M';  %define element input object flags and element properties from el object
            if(model.nlOn)
                elInput.iterationType = nlParams.iterationType;  %define nonlinear iteration type
            else
                elInput.iterationType = 'LINEAR';
            end
            eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
            elInput.displ_iter = eldisp;
            elInput.useDisp = model.nlOn;
            elInput.preStress = false;
            elInput.elementOrder = elementOrder;
            elInput.modalFlag = false;
            elInput.timeInt = struct('delta_t',0.0,...
                'a1',0.0,...
                'a2',0.0,...
                'a3',0.0,...
                'a4',0.0,...
                'a5',0.0,...
                'a6',0.0,...
                'a7',0.0,...
                'a8',0.0);
            elInput.xloc = [0.0 el.elLen(i)];
            elInput.sectionProps = el.props(i);
            elInput.sweepAngle = el.psi(i);
            elInput.coneAngle = el.theta(i);
            elInput.rollAngle = el.roll(i);
            elInput.aeroSweepAngle = 0.0;
            
            elx = zeros(numNodesPerEl,1); %initialize element coordinate list
            ely = elx;
            elz = elx;
            
            for j=1:numNodesPerEl       %define element coordinates and displacements associated with element
                elx(j) = x(conn(i,j));
                ely(j) = y(conn(i,j));
                elz(j) = z(conn(i,j));
                for k=1:numDOFPerNode
                    eldisp(index) = dispEval((conn(i,j)-1)*numDOFPerNode + k);
                    index = index + 1;
                end
            end
            
            %retrieve concentrated nodal terms associated with element
            [massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(conn(i,:),model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad);
            
            %assign concentrated nodal terms and coordinates to element input
            %object
            elInput.concMass = massConc;
            elInput.concStiff = stiffConc;
            elInput.concLoad = loadConc;
            elInput.disp = eldisp;
            elInput.x = elx;
            elInput.y = ely;
            elInput.z = elz;
            elInput.omegaVec = zeros(3,1);
            elInput.omegaDotVec = zeros(3,1);
            
            if(el.rotationalEffects(i)) %activate or deactivate rotational effects for element
                elInput.Omega = Omega;
                elInput.OmegaDot = 0.0;
            else
                elInput.Omega = 0.0;
                elInput.OmegaDot = 0.0;
            end
            
            %deactivate flutter type input
            if(~model.aeroElasticOn)
                elInput.freq = 0.0*2*pi;
                elInput.aeroElasticOn = 0; %turn off for static calculation
            end
            elInput.aeroForceOn = model.aeroForceOn;
            elInput.airDensity = model.airDensity;
            elInput.gravityOn = model.gravityOn;
            elInput.CN2H = eye(3);
            
            elInput.RayleighAlpha = model.RayleighAlpha; %not used for static analysis but in general the calculate element function will look for this data in the input struct
            elInput.RayleighBeta = model.RayleighBeta;
            
            elInput.loadStep = loadStep;
            elInput.dispDerivatives = false;
            [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage(i)); %do element calculation
            
            [Kg,Fg] = assembly(elOutput.Ke,elOutput.Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); %assemble element i into global system
        end
        
        [Fexternal, Fdof] = externalForcingStatic();  %get arbitrary external loads from externalForcingStatic() function
        for i=1:length(Fdof)
            Fg(Fdof(i)) =  Fg(Fdof(i)) + Fexternal(i)*loadStep; %modify assembled global load vector for external loads
        end
        
        %----------------------------------------------------------------------
        
        %apply general 6x6  mass, damping, and stiffness matrices to nodes
        [Kg,~,~] = applyGeneralConcentratedTerms(Kg,Kg,Kg,model.nodalTerms.concStiffGen,model.nodalTerms.concMassGen,model.nodalTerms.concDampGen);
        
        %APPLY BOUNDARY CONDITIONS
        [Kg] = applyConstraints(Kg,model.jointTransform); %modify global stiffness matrix for joint constraints using joint transform
        [Fg] = applyConstraintsVec(Fg,model.jointTransform); %modify global force vector for joint constraints using joint transform
        
        if(BC.numpBC==0)
            disp('*WARNING*: No boundary conditions detected. Fully fixing DOFs at Node 1 to faciliate static solve.');
            BC.numpBC = 6;
            BC.pBC = [1 1 0; 1 2 0; 1 3 0;1 4 0;1 5 0;1 6 0];
        end
        
        [Kg,Fg] = applyBC(Kg,Fg,BC,numDOFPerNode);  %apply boundary conditions to global stiffness matrix and force vector
        dispOld = displ;  %assign displacement vector from previous iteration
        
        if(strcmp(elInput.iterationType,'NR'))  %system solve, norm calculation for newton-raphson iteration
            delta_displ = Kg\Fg;
            delta_displ = model.jointTransform*delta_displ;
            displ = displ + delta_displ;
            uNorm = calculateNorm(displ-delta_displ,displ);
        elseif(strcmp(elInput.iterationType,'DI'))  %system solve, norm calculation for direct iteration
            displ_last = displ;
            displ = Kg\Fg;
            displ = model.jointTransform*displ;
            uNorm = calculateNorm(displ_last,displ);
            gamm = 0.5;
            displ = (1-gamm)*displ + gamm*displ_last;
        else                                        %system solve for linear case
            displ = Kg\Fg;
            displ = model.jointTransform*displ;
            uNorm = 0;
        end
        iterationCount = iterationCount +1;         %increment iteration count
    end
    
    loadStepCount = loadStepCount + 1; %increment load step count
    
    %update load step whether adaptive or prescribed
    [loadStep,loadStepPrev,displ,displPrev,staticAnalysisSuccessfulForLoadStep,staticAnalysisComplete] = updateLoadStep(iterationCount,nlParams,loadStep,loadStepPrev,loadStepCount,displPrev,displ);
    
end
staticAnalysisSuccessful = staticAnalysisSuccessfulForLoadStep;
t_static = toc;
disp('Elapsed time for static analysis(s):');
disp(t_static);

%     reactionNodeNumber = 1; %place holder for nodal reaction force
%     [FReaction] = calculateReactionForceAtNode(reactionNodeNumber,model,mesh,el,...
%         elStorage,[],[],displ,[],Omega,0,[]);

[elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ,model.nlOn);

end


function [uNorm] = calculateNorm(uPrev,u)
%Calculates a relative norm between two vectors: u and uPrev
len = length(u);
numerator = 0;
denom = 0;
for i=1:len
    numerator = numerator + (u(i)-uPrev(i))^2;
    denom = denom + u(i)^2;
end

uNorm = sqrt(numerator/denom);
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
