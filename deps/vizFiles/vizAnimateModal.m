function vizAnimateModal(meshFile,resultsFile,selectedMode,sf,outFileName)
%vizAnimateModal transient animation of modal analysis
% **********************************************************************
% *                   Part of SNL VAWTGen                              *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   vizAnimateModal(meshFile,resultsFile,selectedMode,sf,meshSeg,outFileName)
%                    
%   This function generates an animation of modal analysis results
%   for a structural dynamics simulation performed using the OWENS Toolkit.
%
%      input:
%      meshFile         = string containing mesh file name
%      resultsFile      = string containing results file name
%      sf               = scale factor for displacements
%      outFileName      = string containing filename for generated movie
%
%      output:          (NONE)

tic

close all; %close all open figures
[mesh,meshSeg] = readMeshVG(meshFile); %read mesh file

%read in output file
[modeShape,modeShapeOutPhase] = readResults(resultsFile,selectedMode,mesh.numNodes);

%set scale factor
modeShape = modeShape .* sf;
modeShapeOutPhase = modeShapeOutPhase.* sf;

% Create movie file with required parameters
  fps= 10;
  mov = VideoWriter(outFileName);
  mov.FrameRate = fps;
  mov.Quality = 50;
  open(mov);
    
[axdata] = calculateAxesBounds(mesh,modeShape,modeShapeOutPhase,1.5);

numframes=100;
for i=1:numframes
%     set(fig1,'NextPlot','replacechildren');
        delta_t=2.0;
        t = delta_t*i;
        freq=pi/50;
        timefac1=sin(freq.*t);
        timefac2=sin(freq.*t-pi/2);


%add mode shape*scale factor + original components??
deformedMesh1 = mesh;
deformedMesh1.x = mesh.x + timefac1.*modeShape(:,1)';
deformedMesh1.y = mesh.y + timefac1.*modeShape(:,2)';
deformedMesh1.z = mesh.z + timefac1.*modeShape(:,3)';

deformedMesh2 = mesh;
deformedMesh2.x = mesh.x + timefac2.*modeShapeOutPhase(:,1)';
deformedMesh2.y = mesh.y + timefac2.*modeShapeOutPhase(:,2)';
deformedMesh2.z = mesh.z + timefac2.*modeShapeOutPhase(:,3)';

deformedMesh3 = mesh;
deformedMesh3.x = mesh.x + timefac1.*modeShape(:,1)'+ timefac2.*modeShapeOutPhase(:,1)';
deformedMesh3.y = mesh.y + timefac1.*modeShape(:,2)'+timefac2.*modeShapeOutPhase(:,2)';
deformedMesh3.z = mesh.z + timefac1.*modeShape(:,3)'+timefac2.*modeShapeOutPhase(:,3)';


%plot mesh
plotMeshIso(deformedMesh3,'k-',axdata,meshSeg);

% put this plot in a movieframe
    movegui(gcf);
    F = getframe(gcf);
    writeVideo(mov,F);
end

% save movie
    close(mov);
toc
end

function [axdata] = calculateAxesBounds(mesh,modeShape,modeShapeOutPhase,axfac)
%This figure examines model coordinates, displacement, and scale factors
%and creates a bound for plot axes.

deformedMesh_1.x = mesh.x + 1.0.*modeShape(:,1)'+ 1.0.*modeShapeOutPhase(:,1)';
deformedMesh_1.y = mesh.y + 1.0.*modeShape(:,2)'+ 1.0.*modeShapeOutPhase(:,2)';
deformedMesh_1.z = mesh.z + 1.0.*modeShape(:,3)'+ 1.0.*modeShapeOutPhase(:,3)';

deformedMesh_2.x = mesh.x - 1.0.*modeShape(:,1)'- 1.0.*modeShapeOutPhase(:,1)';
deformedMesh_2.y = mesh.y - 1.0.*modeShape(:,2)'- 1.0.*modeShapeOutPhase(:,2)';
deformedMesh_2.z = mesh.z - 1.0.*modeShape(:,3)'- 1.0.*modeShapeOutPhase(:,3)';

deformedMesh2_1.x = mesh.x - 1.0.*modeShape(:,1)'+ 1.0.*modeShapeOutPhase(:,1)';
deformedMesh2_1.y = mesh.y - 1.0.*modeShape(:,2)'+ 1.0.*modeShapeOutPhase(:,2)';
deformedMesh2_1.z = mesh.z - 1.0.*modeShape(:,3)'+ 1.0.*modeShapeOutPhase(:,3)';

deformedMesh2_2.x = mesh.x + 1.0.*modeShape(:,1)'- 1.0.*modeShapeOutPhase(:,1)';
deformedMesh2_2.y = mesh.y + 1.0.*modeShape(:,2)'- 1.0.*modeShapeOutPhase(:,2)';
deformedMesh2_2.z = mesh.z + 1.0.*modeShape(:,3)'- 1.0.*modeShapeOutPhase(:,3)';

min1x = min(deformedMesh_1.x);
min1y = min(deformedMesh_1.y);
min1z = min(deformedMesh_1.z);
max1x = max(deformedMesh_1.x);
max1y = max(deformedMesh_1.y);
max1z = max(deformedMesh_1.z);

min2x = min(deformedMesh_2.x);
min2y = min(deformedMesh_2.y);
min2z = min(deformedMesh_2.z);
max2x = max(deformedMesh_2.x);
max2y = max(deformedMesh_2.y);
max2z = max(deformedMesh_2.z);

min3x = min(deformedMesh2_1.x);
min3y = min(deformedMesh2_1.y);
min3z = min(deformedMesh2_1.z);
max3x = max(deformedMesh2_1.x);
max3y = max(deformedMesh2_1.y);
max3z = max(deformedMesh2_1.z);

min4x = min(deformedMesh2_2.x);
min4y = min(deformedMesh2_2.y);
min4z = min(deformedMesh2_2.z);
max4x = max(deformedMesh2_2.x);
max4y = max(deformedMesh2_2.y);
max4z = max(deformedMesh2_2.z);

minx = min([min1x min2x min3x min4x]);
maxx = max([max1x max2x max3x max4x]);
miny = min([min1y min2y min3y min4y]);
maxy = max([max1y max2y max3y max4y]);
minz = min([min1z min2z min3z min4z]);
maxz = max([max1z max2z max3z max4z]);

axdata=[minx*axfac maxx*axfac miny*axfac maxy*axfac minz maxz];

end
