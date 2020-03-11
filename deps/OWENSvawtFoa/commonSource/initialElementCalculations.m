function [elStorage] = initialElementCalculations(model,el,mesh)
%initialElementCalculations  performs intitial element calculations
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [elStorage] = initialElementCalculations(model,el,mesh)
%
%   This function performs initial element calculation for use later in
%   analysis for efficiency gains.
%
%      input:
%      model               = object containing model information
%      el                  = object containing element information
%      mesh                = object containing mesh information
%
%      output:
%      elStorage           = object containing stored element data


%initial element calculation
numNodesPerEl = 2;

elStorage = cell(1,mesh.numEl);
for i=1:mesh.numEl
    %Calculate Ke and Fe for element i
    elInput.elementOrder = model.elementOrder; %assign elInput for element i
    elInput.modalFlag = true;
    elInput.xloc = [0.0 el.elLen(i)];
    elInput.sectionProps = el.props{i};
    elInput.sweepAngle = el.psi(i);
    elInput.coneAngle = el.theta(i);
    elInput.rollAngle = el.roll(i);
    elInput.aeroSweepAngle = 0.0;
    
    elx = zeros(1,numNodesPerEl);
    ely = zeros(1,numNodesPerEl);
    elz = zeros(1,numNodesPerEl);
    for j=1:numNodesPerEl
        
        %get element cooridnates
        elx(j) = mesh.x(mesh.conn(i,j));
        ely(j) = mesh.y(mesh.conn(i,j));
        elz(j) = mesh.z(mesh.conn(i,j));
        
    end
    
    %get concentrated terms associated with elemetn
    [massConc,~,~,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(mesh.conn(i,:),model.joint,model.nodalTerms.concMass,model.nodalTerms.concStiff,model.nodalTerms.concLoad);
    
    elInput.x = elx;
    elInput.y = ely;
    elInput.z = elz;
    elInput.concMassFlag = ~isempty(find(massConc,1));
    elInput.concMass = massConc; %only needed for structure mass props (not used in saved element matrices)
    
    elInput.Omega = 0;
    
    [elStorage{i}] = calculateTimoshenkoElementInitialRun(elInput); %initial element calculations for storage
end

end
