function [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage)
%calculateStructureMassProps   calculates mass properties of mesh
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage)
%
%   This function caclulates structural mass properties of the finite
%   element mesh (mass, moment of inertia, mass center) about the origin of
%   the mesh coordinate system.
%
%      input:
%      elStorage    = object containing arrays of stored element
%                     information
%
%      output:
%      structureMass       = mass of structure
%      structureMOI        = moment of inertia tensor of structgure
%      structureMassCenter = center of mass of structure

numElements = length(elStorage); %get number of elements

structureMass = 0; %initialize structure mass and moment of inertia
structureMOI = zeros(3);
temp = zeros(3,1);
for i=1:numElements %sum over elemetns contribution to mass and moment of inertia
    structureMass = structureMass + elStorage{i}.mel;
    structureMOI = structureMOI + elStorage{i}.moiel;
    temp = temp + elStorage{i}.xmel;
end

structureMassCenter = temp./structureMass; %calculate mass center

%modify moment of inertia to be about structure mass center
x = structureMassCenter(1);
y = structureMassCenter(2);
z = structureMassCenter(3);

structureMOI = structureMOI - structureMass*[(y^2+z^2), -x*y, -x*z;
    -x*y, (x^2+z^2),-y*z;
    -x*z,-y*z,(x^2+y^2)];

end