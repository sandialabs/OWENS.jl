function plotMeshIso(mesh,colorstring,axdata,meshSeg)
#plotMeshIso plots a mesh in a 3D with isometric view
# **********************************************************************
# *                   Part of SNL VAWTGen                              *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   plotMeshIso(mesh,colorstring,axdata,meshSeg)
#                    
#   This function plots a mesh in a 3D with isometric view.
#
#      input:
#      mesh            = object containing mesh data
#      colorstring     = string containing color for plotting
#      axdata          = axes bounds for 3D plot
#      meshSeg         = array containing number of elements per mesh
#                        segment
#
#
#      output:          (NONE)

conn = mesh.conn; #assign mesh connectivity and coordinates
x=mesh.x;
y=mesh.y;
z=mesh.z;

hold off;

index = 0;
for i=1:length(meshSeg)  #loop over mesh segments
    for j=1:meshSeg(i)  #loop over elements in mesh segment i
        
        index=index+1;
        n1=conn(index,1); #get node numbers for element j in mesh segment i
        n2=conn(index,2);
   
        #get nodal coordinates for element j in mesh segment i
        p1(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [x(n1), x(n2)]; 
        p2(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [y(n1), y(n2)];
        p3(:,((j-1)*2 + 1):((j-1)*2 + 2)) = [z(n1), z(n2)];

    end
    
    #plot element in 3D window
    plot3(p1,p2,p3,colorstring,'LineWidth',2);
    clear p1 p2 p3;
    hold on;
end


   #set axes 
   axis equal;
   axis(axdata);

   grid on;
   view(3);
   set(gcf,'Color',[1  1 1]);

end