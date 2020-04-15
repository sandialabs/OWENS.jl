function staticPlotter(meshFile,resfile,sf)
%staticPlotter  plots deformed mesh from static elasticity analysis
% **********************************************************************
% *                   Part of SNL VAWTGen                              *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   staticPlotter(meshFile,resfile,sf)
%                    
%   This functino plots the deformed mesh from a static elasticity
%   analysis.
%
%      input:
%      meshFile      = string containing mesh file name
%      resFile       = .mat file containing results of static analysis
%      sf            = scale factor for deformed mesh
%
%      output:      (NONE)
%
close all;
load(resfile);
[mesh,meshSeg] = readMeshVG(meshFile);

numNodes = length(displ)/6;
for i=1:numNodes
   u(i) = displ((i-1)*6 + 1); 
   v(i) = displ((i-1)*6 + 2); 
   w(i) = displ((i-1)*6 + 3); 
end


%add mode shape*scale factor + original components??
deformedMesh = mesh;
deformedMesh.x = deformedMesh.x + sf.*u;
deformedMesh.y = deformedMesh.y + sf.*v;
deformedMesh.z = deformedMesh.z + sf.*w;

%plot meshes

plotMesh(deformedMesh,'-r',meshSeg);
plotMesh(mesh,'-k',meshSeg);
grid on;

dispmag = sqrt(u.*u + v.*v + w.*w);
maxdisp = max(dispmag)

end