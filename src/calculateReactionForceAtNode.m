function [cummulativeForce] = calculateReactionForceAtNode(nodeNum,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
%calculateReactionForceAtNode calculates reaction force at a node
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [cummulativeForce] = calculateReactionForceAtNode(nodeNum,model,mesh,...
%    el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
%
%   This function calculates the reaction force at a node by post
%   processing all element associated with a node through connectivity or
%   joint constraints.
%
%   input:
%   nodeNum    = node number joint constraints are desired at
%   model      = object containing model data
%   mesh       = object containing mesh data
%   elStorage  = object containing stored element data
%   el         = object containing element data
%   timeInt    = object containing time integration parameters
%   dispData   = object containing displacement data
%   displ_iter = converged displacement solution
%   rbData     = vector containing rigid body displacement, velocity, and
%                acceleration
%   Omega      = rotor speed (Hz)
%   OmegaDot   = rotor acceleratin (Hz)
%   CN2H       = transformation matrix from inertial frame to hub frame
%
%   output:
%   cummulativeForce  = vector containing reaction force at nodeNum
conn = mesh.conn;  %get connectivity list
numDofPerNode = 6;

cummulativeForce = zeros(numDofPerNode,1); %initialize force at node

%find elements associated with nodeNum due to mesh connectivity or
%constraints
[elList,elListLocalNodeNumbers] = findElementsAssociatedWithNodeNumber(nodeNum,conn,model.joint);

%process elements for nodal reaction forces and compile to find total
%reaction force at specified node
for i=1:length(elList)
    [Fpp] = elementPostProcess(elList(i),model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H);
    localNode = elListLocalNodeNumbers(i);
    cummulativeForce = cummulativeForce + Fpp((localNode-1)*numDofPerNode+1:(localNode-1)*numDofPerNode+6);
end

end
