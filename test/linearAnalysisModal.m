function  [freq,damp,phase1,phase2,imagCompSign] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage)
%linearAnalysisModal performs modal analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freq,damp,phase1,phase2,imagCompSign] = linearAnalysisModal(model,mesh,
%                                                 el,displ,Omega,elStorage)
%
%   This function performs a modal analysis of a structural dynamics
%   system and returns freq, damping, and mode shapes
%
%   input:
%   model          = object containing model information
%   mesh           = object containing mesh information
%   el             = object containing element information
%   displ          = displacement vector for use in pre-stressed analysis
%   Omega          = rotor speed (Hz)
%   elStorage      = previously calculated element system matrices
%
%
%   output:
%   freq         = modal frequencies
%   damp         = modal damping ratios
%   phase1       = in phase mode shapes (real part of mode shape)
%   phase2       = out of phase mode shapes (imaginary part of mode shape)
%   imagCompSign = sign of imaginary component of eigenvalues
tic

numEl = mesh.numEl;  %extract mesh information from mesh object
x = mesh.x;
y = mesh.y;
z = mesh.z;
conn = mesh.conn;
numNodes = length(x);

elementOrder = model.elementOrder;  %extract element order from model

BC = model.BC; %extract boundary  condition object from model

numNodesPerEl = elementOrder + 1;  %do initialization
numDOFPerNode = 6;
totalNumDOF = numNodes * numDOFPerNode;

Kg = zeros(totalNumDOF,totalNumDOF);
Mg = zeros(totalNumDOF,totalNumDOF);
Cg = zeros(totalNumDOF,totalNumDOF);

nodalTerms = model.nodalTerms; %extract concentrated nodal terms from model


for i=1:numEl   %element loop
    eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
    %Calculate Ke and Fe for element i
    index = 1;
    elInput.analysisType = 'M';      %define element input object flags and element properties from el object
    elInput.iterationType = 'NONE';
    elInput.displ_iter = eldisp;
    elInput.useDisp = false;
    elInput.preStress = true;
    elInput.elementOrder = elementOrder;
    elInput.modalFlag = true;
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
            eldisp(index) = displ((conn(i,j)-1)*numDOFPerNode + k);
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
    
    elInput.aeroElasticOn = model.aeroElasticOn;   %set aeroelastic flag
    elInput.aeroForceOn = false;
    elInput.gravityOn = model.gravityOn;
    elInput.CN2H = eye(3);
    elInput.RayleighAlpha = model.RayleighAlpha;
    elInput.RayleighBeta = model.RayleighBeta;
    
    if(model.aeroElasticOn)
        elInput.freq = model.guessFreq*2.0*pi;     %set guess frequency if aeroelastic analysis
        elInput.airDensity = model.airDensity;     %set air density if aeroelastic analysis
    end
    
    [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage(i)); %do element calculation
    
    [Kg] = assemblyMatrixOnly(elOutput.Ke,conn(i,:),numNodesPerEl,numDOFPerNode,Kg); %assemble element into global stiffness matrix
    [Mg] = assemblyMatrixOnly(elOutput.Me,conn(i,:),numNodesPerEl,numDOFPerNode,Mg); %assemble element into global mass matrix
    [Cg] = assemblyMatrixOnly(elOutput.Ce,conn(i,:),numNodesPerEl,numDOFPerNode,Cg); %assemble element into global damping matrix
    
end

%apply general 6x6  mass, damping, and stiffness matrices to nodes
[Kg,Mg,Cg] = applyGeneralConcentratedTerms(Kg,Mg,Cg,model.nodalTerms.concStiffGen,model.nodalTerms.concMassGen,model.nodalTerms.concDampGen);

%----------------------------------------------------------------------
%APPLY CONSTRAINT
[Kg] = applyConstraints(Kg,model.jointTransform);  %modify global matrices for joint constraints using joint transform
[Mg] = applyConstraints(Mg,model.jointTransform);
[Cg] = applyConstraints(Cg,model.jointTransform);

%APPLY BOUNDARY CONDITIONS
[KgTotal] = applyBCModal(Kg,BC.numpBC,BC.map);     %apply boundary conditions to global matrices
[MgTotal] = applyBCModal(Mg,BC.numpBC,BC.map);
[CgTotal] = applyBCModal(Cg,BC.numpBC,BC.map);

if(Omega==0.0) %set eigensolver flag
    solveFlag = 2;
else
    solveFlag = 2;
end
[eigVec,eigVal] = eigSolve(MgTotal,CgTotal,KgTotal,... %eigensolve of global system
    model.numModesToExtract,solveFlag);
save eigVectors eigVec %save eigenvector for later use (if needed)


%extract frequency, damping, mode shapes from eigenvalues and vectors
[~,len] = size(eigVal);
% freq = zeros(len);
% damp = zeros(len);
% phase1 = zeros(len);
% phase2 = zeros(len);
% sortedModes = zeros(len);
% imagCompSign = zeros(len);
for i=1:len
    [freq(i),damp(i),phase1(:,:,i),phase2(:,:,i),sortedModes(:,:,i)] = extractFreqDamp(eigVal(i,i),eigVec(:,i),numDOFPerNode,model.jointTransform,model.reducedDOFList,model.BC,model.analysisType);
    imagCompSign(i) = sign(imag(eigVal(i,i)));
end

%write output
t_modal = toc;
disp('Elapsed time for modal analysis(s):');
disp(t_modal);

if(~strcmp(model.analysisType,'FA'))
    fidout=fopen(model.outFilename,'wt');
    [freq,damp,imagCompSign] = writeOutput(freq,damp,phase1,phase2,imagCompSign,fidout);
    fclose(fidout);
end

end

function [Kg] = applyConstraints(Kg,transMatrix)
%This function transforms a matrix by the transformation matrix to
%enforce joint constraints
transMatrix = sparse(transMatrix);
Kg = sparse(Kg);
Kg = transMatrix'*Kg*transMatrix;
Kg = full(Kg);
end
