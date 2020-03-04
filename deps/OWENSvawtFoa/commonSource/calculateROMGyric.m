function  [rom] = calculateROMGyric(model,mesh,el,displ,omegaVec,omegaDotVec,accelVec,elStorage,rom0)
%calculateROMGyric Calculates reduced order model for a Gyric system
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [rom] = calculateROMGyric(model,mesh,el,displ,omegaVec,omegaDotVec,elStorage,rom0)
%                    
%   This function calculates a reduced order model with rotational/ rigid
%   body motion effects
%
%      input:
%      model        =  object containing model data
%      mesh         =  object containing mesh data
%      el           =  object containing elementdata
%      displ        =  displacement vector
%      omegaVec     =  vector of hub frame angular velocity components
%      omegaDotVec  =  vector of hub frame angular acceleration components
%      elStorage    =  object containing stored element data
%      rom0         =  object containing parked/conventional reduced order
%                      model

%      output:
%      rom          =  object containing reduced order model data

%initialization
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

Kg = zeros(totalNumDOF,totalNumDOF);
Cg = zeros(totalNumDOF,totalNumDOF);
Fg = zeros(totalNumDOF,1);

nodalTerms = model.nodalTerms;

    for i=1:numEl    %element loop
        %Calculate Ke and Fe for element i
        index = 1;                     %assign element i properties to elInput
        elInput.analysisType = 'RM0';  %set analysis type to initial ROM calculations
        elInput.useDisp = true;
        elInput.preStress = false;
        elInput.iterationType = 'LINEAR'; %linear analysis option
        elInput.elementOrder = elementOrder;
        elInput.modalFlag = true;
        elInput.xloc = [0.0 el.elLen(i)];
        elInput.sectionProps = el.props{i};
        elInput.sweepAngle = el.psi(i);
        elInput.coneAngle = el.theta(i);
        elInput.rollAngle = el.roll(i);
        elInput.aeroSweepAngle = 0.0;
        
        for j=1:numNodesPerEl %get element coordinates and displacements
            elx(j) = x(conn(i,j));
            ely(j) = y(conn(i,j));
            elz(j) = z(conn(i,j));
            for k=1:numDOFPerNode
                eldisp(index) = displ((conn(i,j)-1)*numDOFPerNode + k);
                index = index + 1;
            end
        end
        %get concentrated terms associated with element
        [massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(conn(i,:),model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad);
        
        elInput.concMass = massConc;
        elInput.concStiff = stiffConc;
        elInput.concLoad = loadConc;
        elInput.disp = eldisp;
        elInput.x = elx;
        elInput.y = ely;
        elInput.z = elz;
        Omega = 0.0;
        elInput.accelVec = accelVec;
        if(el.rotationalEffects(i))
            elInput.Omega = Omega;
            elInput.omegaVec = omegaVec;
            elInput.omegaDotVec = omegaDotVec;
            elInput.OmegaDot = 0;
        else
            elInput.Omega = 0.0;
            elInput.omegaVec = zeros(3,1);
            elInput.omegaDotVec = zeros(3,1);
            elInput.OmegaDot = 0.0;
        end
        
        elInput.aeroElasticOn = model.aeroElasticOn;
        elInput.aeroForceOn = false;
        elInput.gravityOn = model.gravityOn;
        
        elInput.RayleighAlpha = model.RayleighAlpha;
        elInput.RayleighBeta = model.RayleighBeta;
        
        if(model.aeroElasticOn)  %make assignments if aeroelastic analysis is active (needs check for consistency with ROM)
            elInput.freq = model.guessFreq*2.0*pi;
            elInput.airDensity = model.airDensity;
        end
        
        [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage{i});   %calculate element
        
        [Kg,Fg] = assembly(elOutput.Ke,elOutput.Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); %assemble element Kg and Fg 
        [Cg] = assemblyMatrixOnly(elOutput.Ce,conn(i,:),numNodesPerEl,numDOFPerNode,Cg); %assemble Mg
  
    end
    
    %----------------------------------------------------------------------
    %APPLY CONSTRAINT
    [Kg] = applyConstraints(Kg,model.jointTransform); %apply constraints to system matrices and force vector
    [Cg] = applyConstraints(Cg,model.jointTransform);
    [Fg] = applyConstraintsVec(Fg,model.jointTransform);
     
          
    %APPLY BOUNDARY CONDITIONS
    [KgTotal] = applyBCModal(Kg,model.BC.numpBC,BC.map); %apply boundary conditions to system matrices and force vector
    [CgTotal] = applyBCModal(Cg,model.BC.numpBC,BC.map);
    [FgTotal] = applyBCModalVec(Fg,model.BC.numpBC,BC.map);
    
    Phi = rom0.Phi;
    
    rom.Kr = Phi'*KgTotal*Phi - rom0.Kr; %removes contribution of structural stiffness to reduced K
    rom.Cr = Phi'*CgTotal*Phi - rom0.Cr; %removes contribution of structural damping to reduced C
    rom.Fr = Phi'*FgTotal - rom0.Fr;
      
end

function [Kg] = applyConstraints(Kg,transMatrix)
    %This function transforms a matrix by the transformation matrix to
    %enforce joint constraints
     Kg = transMatrix'*Kg*transMatrix;
end

function [U] = applyConstraintsVec(U,transMatrix)
    %This function transforms a vector by the transformation matrix to
    %enforce joint constraints
    U = transMatrix'*U;
end
