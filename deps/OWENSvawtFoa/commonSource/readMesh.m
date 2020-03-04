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

temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements
numNodes = temp(1);
numEl = temp(2);

for i=1:numNodes            % read in node number and node coordinates
   nodeNum(i) = fscanf(fid,'%i',1);
   temp = fscanf(fid,'%f',3);
   x(i) = temp(1);
   y(i) = temp(2);
   z(i) = temp(3);
end

for i=1:numEl               % read in element number and connectivity list
   elNum(i) = fscanf(fid,'%i',1);
   dum = fscanf(fid,'%i',1);
   conn(i,:) = fscanf(fid,'%i',2);
   
   
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