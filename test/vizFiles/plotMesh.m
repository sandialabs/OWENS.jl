function plotMesh(mesh,colorstring,meshSeg)
%plotMesh plots a mesh with various orthographic and isometric views
% **********************************************************************
% *                   Part of SNL VAWTGen                              *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%    plotMesh(mesh,colorstring,meshSeg)
%                    
%   This function plots a mesh in a 3D with isometric view.
%
%      input:
%      mesh            = object containing mesh data
%      colorstring     = string containing color for plotting
%      meshSeg         = array containing number of elements per mesh
%                        segment
%
%
%      output:          (NONE)

    conn = mesh.conn; %assign mesh connectivity and coordinates
    x=mesh.x;
    y=mesh.y;
    z=mesh.z;

    hold on;

    index=0;
    for i=1:length(meshSeg) %loop over mesh segments
       for j=1:meshSeg(i) %loop over elements in mesh segment i
            index=index+1; 
            n1=conn(index,1); %get node numbers for element j in mesh segment i
            n2=conn(index,2);

            %get nodal coordinates for element j in mesh segment i
            p1(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [x(n1), x(n2)];
            p2(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [y(n1), y(n2)];
            p3(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [z(n1), z(n2)];
       end

       %plot meshes in various figure windows
       hold on;
       subplot(2,2,1);
       plot3(p1,p2,p3,colorstring,'LineWidth',2);
       xlabel('h_1','FontSize',15); ylabel('h_2','FontSize',15); zlabel('h_3','FontSize',15);

       hold on;
       subplot(2,2,2);
       plot3(p1,p2,p3,colorstring,'LineWidth',2);
       xlabel('h_1','FontSize',15); ylabel('h_2','FontSize',15); zlabel('h_3','FontSize',15);

       hold on;
       subplot(2,2,3);
       plot3(p1,p2,p3,colorstring,'LineWidth',2);
       xlabel('h_1','FontSize',15); ylabel('h_2','FontSize',15); zlabel('h_3','FontSize',15);

       hold on;
       subplot(2,2,4);
       plot3(p1,p2,p3,colorstring,'LineWidth',2);
       xlabel('h_1','FontSize',15); ylabel('h_2','FontSize',15); zlabel('h_3','FontSize',15);

       clear p1 p2 p3;
   
    end
   %set orthographic views
   subplot(2,2,1);
   axis equal;
   view(0,0);
   
   subplot(2,2,2);
   axis  equal;
   view(0,90);
   
   subplot(2,2,3);
   axis equal;
   view(90,0);
   
   %set isoparametric views
   subplot(2,2,4);
   axis equal;
   view(3);

end