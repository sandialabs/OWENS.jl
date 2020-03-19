function [Fpp] = elementPostProcess(elementNumber,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
%elementPostProcess post processes element for reaction force
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Fpp] = elementPostProcess(elementNumber,model,mesh,el,elStorage,....
%           timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
%
%   This function calculates the reaction force associated with an element.
%
%   input:
%   elementNumber  = node number joint constraints are desired at
%   model          = object containing model data
%   mesh           = object containing mesh data
%   elStorage      = object containing stored element data
%   el             = object containing element data
%   timeInt        = object containing time integration parameters
%   dispData       = object containing displacement data
%   displ_iter     = converged displacement solution
%   rbData         = vector containing rigid body displacement, velocity,
%                     and acceleration
%   Omega          = rotor speed (Hz)
%   OmegaDot       = rotor acceleratin (Hz)
%   CN2H           = transformation matrix from inertial frame to hub frame
%
%   output:
%   Fpp            = vector containing reaction force vector associated
%                    with element


%some initializations
elementOrder = model.elementOrder;
numNodesPerEl = elementOrder + 1;
numDOFPerNode = 6;
elx=zeros(numNodesPerEl,1);
ely=zeros(numNodesPerEl,1);
elz=zeros(numNodesPerEl,2);

eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
eldisp_sm1 = zeros(1,numNodesPerEl*numDOFPerNode);
eldispdot = eldisp;
eldispddot = eldisp;
eldispiter = eldisp;
disp_s = dispData.displ_s;
% not always used, but must be declared, specific to TNB and ROM
dispdot_s = dispData.displdot_s;
dispddot_s = dispData.displddot_s;
%unpack displacement information
analysisType = model.analysisType;

%unpack mesh information
x = mesh.x;
y = mesh.y;
z = mesh.z;
conn = mesh.conn;

%construct elInput for elementNumber
index = 1;
elInput.elementOrder = elementOrder;
elInput.modalFlag = true;
elInput.timeInt = timeInt;
elInput.xloc = [0.0 el.elLen(elementNumber)];
elInput.sectionProps = el.props(elementNumber);
elInput.sweepAngle = el.psi(elementNumber);
elInput.coneAngle = el.theta(elementNumber);
elInput.rollAngle = el.roll(elementNumber);
elInput.aeroSweepAngle = 0.0;
elInput.iterationType = 'DI';
elInput.useDisp = model.nlOn;
elInput.preStress = false;
elInput.aeroElasticOn = false;
elInput.aeroForceOn = false;

if(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM'))
    elInput.analysisType = 'TNB';
elseif(strcmp(analysisType,'S'))
    elInput.analysisType = 'M';
else
    error('analysisType not supported, choose another')
end

%unpack connectivity list, nodal terms, etc.
nodalTerms = model.nodalTerms;

for j=1:numNodesPerEl
    for k=1:numDOFPerNode
        %get element cooridnates
        elx(j) = x(conn(elementNumber,j));
        ely(j) = y(conn(elementNumber,j));
        elz(j) = z(conn(elementNumber,j));
        if(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM'))
            eldisp(index) = disp_s((conn(elementNumber,j)-1)*numDOFPerNode + k);
            eldispdot(index) = dispdot_s((conn(elementNumber,j)-1)*numDOFPerNode + k);
            eldispddot(index) = dispddot_s((conn(elementNumber,j)-1)*numDOFPerNode + k);
            eldispiter(index) = displ_iter((conn(elementNumber,j)-1)*numDOFPerNode + k);
        elseif(strcmp(analysisType,'S'))
            eldisp(index) = displ_iter((conn(elementNumber,j)-1)*numDOFPerNode + k);
        end
        index = index + 1;
    end
end

elInput.disp = eldisp;
elInput.dispdot = eldispdot;%specific to TNB and ROM
elInput.dispddot = eldispddot;%specific to TNB and ROM

if(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM'))
    elInput.displ_iter = eldispiter;
elseif(strcmp(analysisType,'S'))
    elInput.displ_iter = eldisp;
    eldispiter = eldisp;
else
    error('analysisType not supported, choose another')
end

%get concentrated terms associated with element
[massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(conn(elementNumber,:),model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad);

elInput.concMass = massConc;
elInput.concStiff = stiffConc;
elInput.concLoad = loadConc;
elInput.disp = eldisp;
elInput.dispm1 = eldisp_sm1;
elInput.x = elx;
elInput.y = ely;
elInput.z = elz;
elInput.gravityOn = model.gravityOn;

elInput.RayleighAlpha = model.RayleighAlpha;
elInput.RayleighBeta = model.RayleighBeta;

if(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM'))
    elInput.accelVec = rbData(1:3);
    elInput.omegaVec = rbData(4:6);
    elInput.omegaDotVec = rbData(7:9);
elseif(strcmp(analysisType,'S'))
    elInput.accelVec = zeros(3,1);
    elInput.omegaVec = zeros(3,1);
    elInput.omegaDotVec = zeros(3,1);
else
    error('analysisType not supported, choose another')
end

if(el.rotationalEffects(elementNumber))
    elInput.Omega = Omega;
    elInput.OmegaDot = OmegaDot;
else
    elInput.Omega = 0.0;
    elInput.OmegaDot = 0.0;
end

% Specific to TNB and ROM, but must be declared
elInput.CN2H = CN2H;
    

elInput.airDensity = model.airDensity;
elInput.freq = 0.0; %Is not used for this model type, but must be declared.
elInput.firstIteration = false;

%calculate element stiffness matrix and force vector
%(or effective stiffness matrix and force vector from time integration)
[elOutput] = calculateTimoshenkoElementNL(elInput,elStorage(elementNumber));

%post process for reaction force
FhatEl1PP = elOutput.Ke*eldispiter';
if(strcmp(analysisType,'TD'))
    denom = timeInt.a4;
elseif(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM')||strcmp(analysisType,'S'))
    denom = 1.0;
else
    error('analysisType not supported, choose another')
end
Fpp = (FhatEl1PP - elOutput.FhatLessConc)./denom;
%----------------------------------------------------------------------

end
