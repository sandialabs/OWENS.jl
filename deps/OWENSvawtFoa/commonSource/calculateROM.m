function  [rom] = calculateROM(model,mesh,el,displ,omegaVec,omegaDotVec,elStorage)
%calculateROM Calculates reduced order model (conventional system)
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [rom] = calculateROM(model,mesh,el,displ,omegaVec,omegaDotVec,elStorage)
%                    
%   This function calculates a reduced order model for a conventional
%   structural dynamics system (parked, non-rotating)
%
%      input:
%      model        =  object containing model data
%      mesh         =  object containing mesh data
%      el           =  object containing elementdata
%      displ        =  displacement vector
%      omegaVec     =  vector of hub frame angular velocity components
%      omegaDotVec  =  vector of hub frame angular acceleration  components
%      elStorage    =  object containing stored element data

%      output:
%      rom          =  object containing reduced order model data

%initializations
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
Mg = zeros(totalNumDOF,totalNumDOF);
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
%                 eldisp((j-1)*numDOFPerNode+k) = disp( (conn(i,j)-1)*numDOFPerNode + k);
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
        elInput.accelVec = zeros(3,1);
        elInput.gravityOn = model.gravityOn;
        
        elInput.RayleighAlpha = model.RayleighAlpha;
        elInput.RayleighBeta = model.RayleighBeta;
        
        if(el.rotationalEffects(i))
            elInput.Omega = Omega;
            elInput.omegaVec = omegaVec;
            elInput.omegaDotVec = omegaDotVec;
            elInput.OmegaDot = 0.0;
        else
            elInput.Omega = 0.0;
            elInput.omegaVec = omegaVec;
            elInput.omegaDotVec = omegaDotVec;
            elInput.OmegaDot = 0.0;
        end
        
        elInput.aeroElasticOn = model.aeroElasticOn;
        elInput.aeroForceOn = false;
        if(model.aeroElasticOn)   %make assignments if aeroelastic analysis is active (needs check for consistency with ROM)
            elInput.freq = model.guessFreq*2.0*pi;
            elInput.airDensity = model.airDensity;
        end
        
        [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage{i});   %calculate element
         
        [Kg,Fg] = assembly(elOutput.Ke,elOutput.Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); %assemble element Kg and Fg 
        [Mg] = assemblyMatrixOnly(elOutput.Me,conn(i,:),numNodesPerEl,numDOFPerNode,Mg); %assemble Mg
        [Cg] = assemblyMatrixOnly(elOutput.Ce,conn(i,:),numNodesPerEl,numDOFPerNode,Cg); %assemble Cg
      
    end


    %----------------------------------------------------------------------
    %APPLY CONSTRAINT
    [Kg] = applyConstraints(Kg,model.jointTransform);   %apply constraints to system matrices and force vector
    [Mg] = applyConstraints(Mg,model.jointTransform);
    [Cg] = applyConstraints(Cg,model.jointTransform);
    [Fg] = applyConstraintsVec(Fg,model.jointTransform);
   
           
    %APPLY BOUNDARY CONDITIONS
    [KgTotal] = applyBCModal(Kg,model.BC.numpBC,BC.map); %apply boundary conditions to system matrices and force vector
    [MgTotal] = applyBCModal(Mg,model.BC.numpBC,BC.map);
    [CgTotal] = applyBCModal(Cg,model.BC.numpBC,BC.map);
    [FgTotal] = applyBCModalVec(Fg,model.BC.numpBC,BC.map);
    
    if(Omega==0.0)
        solveFlag = 3;  %set solve flag to 3 for parked modal analysis (no Gyric effects: spring mass system)
    else
        solveFlag = 3;
    end
   [eigVec,eigVal] = eigSolve(MgTotal,CgTotal,KgTotal,...   %calculate eigenvalues and eigvenvectors of system
                      model.numModesForROM,solveFlag);
    len=length(KgTotal);                  
    Phi = real(eigVec);   %assign modal matrix to Phi
    rom.Phi = Phi;        
    rom.invPhi = pinv(Phi); %calculate psuedo inverse of Phi
    
    rom.Mr = Phi'*MgTotal*Phi;   %calculate ROM mass matrix
    rom.Kr = Phi'*KgTotal*Phi;   %calculate ROM stiffness matrix
    rom.Cr = Phi'*CgTotal*Phi;   %calculate ROM damping matrix
    rom.Fr = Phi'*FgTotal;       %calculate ROM force vector
    
    %apply modal damping
    dampRatio = 0.02;
    Cmodal = zeros(length(rom.Kr));
    for i=1:length(rom.Kr)
       Cmodal(i,i) = 2*dampRatio*sqrt(rom.Kr(i,i)*rom.Mr(i,i));
    end
    rom.CrModal = Cmodal;
      
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
