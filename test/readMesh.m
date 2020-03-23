function [mesh] = readMesh(filename)
%readMesh  reads mesh file and stores data in mesh object
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [mesh] = readMesh(filename)
%
%   This function reads the mesh file and stores data in the mesh object.
%
%      input:
%      filename      = string containing mesh filename

%      output:
%      mesh          = object containing mesh data

fid = fopen(filename,'r');   %open mesh file

% temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements
temp = getSplitLine(fid,'	');
numNodes = temp(1);
numEl = temp(2);

nodeNum = zeros(numNodes,1);
x = zeros(numNodes,1);
y = zeros(numNodes,1);
z = zeros(numNodes,1);

conn = zeros(numEl,2);
elNum = zeros(numEl,1);
temp = zeros(1,4);
    
for i=1:numNodes            % read in node number and node coordinates
    temp = getSplitLine(fid,'	');
    nodeNum(i) = temp(1);
    x(i) = temp(2);
    y(i) = temp(3);
    z(i) = temp(4);
end

for i=1:numEl               % read in element number and connectivity list
    temp = getSplitLine(fid,'	');
    elNum(i) = temp(1);
    
    conn(i,:) = temp(3:4);
end

fclose(fid);  %close mesh file

mesh.nodeNum = nodeNum; %store data in mesh object
mesh.numEl = numEl;
mesh.numNodes = numNodes;
mesh.x = x;
mesh.y = y;
mesh.z = z;
mesh.elNum = elNum;
mesh.conn = conn;

end
