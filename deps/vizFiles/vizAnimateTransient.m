function vizAnimateTransient(meshFile,uData,sf,outFileName)
%vizAnimateTransient transient animation of transient simulation
% **********************************************************************
% *                   Part of SNL VAWTGen                              *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   vizAnimateTransient(meshFile,uData,sf,meshSeg,outFileName)
%                    
%   This function generates an animation of transient simulation results
%   for a structural dynamics simulation performed using the OWENS Toolkit.
%
%      input:
%      meshfile         = string containing mesh file name
%      uData            = array containing displacement history data
%      sf               = scale factor for displacements
%      outFileName      = string containing filename for generated movie
%
%      output:          (NONE)

tic
close all; %close all open figures;
[mesh,meshSeg] = readMeshVG(meshFile);  %read mesh file

% Create movie file with required parameters
  fps= 40;
  mov = VideoWriter(outFileName);
  mov.FrameRate = fps;
  mov.Quality = 50;
  open(mov);
    
[numDof,numframes] = size (uData);

% U = zeros(numframes,1);
% V = zeros(numframes,1);
% W = zeros(numframes,1);
% if(platformOn)
%     numDof = numDof - 6;
%     U = uData(numDof+1,:);
%     V = uData(numDof+2,:);
%     W = uData(numDof+3,:);
%     ROT1 = uData(numDof+4,:);
%     ROT2 = uData(numDof+5,:);
%     ROT3 = uData(numDof+6,:);
% end

%unpack u,v,w displacement data from uData array
numDofPerNode = 6;
for i=1:numDof/numDofPerNode
	u(i,:) = uData((i-1)*6 + 1,:);
	v(i,:) = uData((i-1)*6 + 2,:);
	w(i,:) = uData((i-1)*6 + 3,:);
end

%calculate bounds for axes window
[axdata] = calculateAxesBounds(mesh,u,v,w,sf,1.2);


for i=1:numframes %loop over number of frames

    %add mode shape*scale factor + original components??
        deformedMesh = mesh;
        deformedMesh.x = mesh.x + u(:,i)'.*sf; 
        deformedMesh.y = mesh.y + v(:,i)'.*sf;
        deformedMesh.z = mesh.z + w(:,i)'.*sf;

    %plot mesh
        plotMeshIso(deformedMesh,'k-',axdata,meshSeg);

    % put this plot in a movieframe
        movegui(gcf);
        F = getframe(gcf);
        writeVideo(mov,F);
end

% save movie
    close(mov);
toc
end


function [axdata] = calculateAxesBounds(mesh,u,v,w,sf,axfac)
%This figure examines model coordinates, displacement, and scale factors
%and creates a bound for plot axes.

%get bounds on axes
minx = 1e6;
miny = 1e6;
minz = 1e6;
maxx = -1e6;
maxy = -1e6;
maxz = -1e6;

[~,wid]=size(u);
for i=1:wid
    xd = mesh.x + u(:,i)'*sf;
    yd = mesh.y + v(:,i)'*sf;  
    zd = mesh.z + w(:,i)'*sf;
    
    if(maxx<max(xd))
        maxx=max(xd);
    end
    if(maxy<max(yd))
        maxy=max(yd);
    end
    if(maxz<max(zd))
        maxz=max(zd);
    end
    if(minx>min(xd))
        minx=min(xd);
    end
    if(miny>min(yd))
        miny=min(yd);
    end
    if(minz>min(zd))
        minz=min(zd);
    end
end


axdata=[minx*axfac maxx*axfac miny*axfac maxy*axfac minz maxz];

del = 1;
if(axdata(1) == axdata(2))
    axdata(1) = -del;
    axdata(2) = del;
end

if(axdata(3) == axdata(4))
    axdata(3) = -del;
    axdata(4) = del;
end

if(axdata(5) == axdata(6))
    axdata(5) = -del;
    axdata(6) = del;
end

end
