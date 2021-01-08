function viz(meshFile,resultsFile,selectedMode,sf)
    #viz  visualizes mode shapes
    #   viz(meshFile,resultsFile,selectedMode,sf,meshSeg)
    #
    #   This function plots the mode shapes of a mode from a modal analysis
    #   performed using the OWENS toolkit.
    #
    #      input:
    #      meshFile         = string containing mesh file name
    #      resultsFile      = string containing results file name
    #      selectedMode     = integer denoting selected mode to plot
    #      sf               = scale factor for mode shape displacements
    #
    #      output:          (NONE)

    #     close all; #close all open figure windows #ble
    mesh = readMesh(meshFile) #read mesh file

    #read in output file
    freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90 = readResultsModalOut(resultsFile,mesh.numNodes)#ble: changed from readResults

    #add mode shape*scale factor + original components??
    deformedMesh = deepcopy(mesh)
    deformedMesh.x = deformedMesh.x .+ U_x_0[:,selectedMode] .* sf
    deformedMesh.y = deformedMesh.y .+ U_y_0[:,selectedMode] .* sf
    deformedMesh.z = deformedMesh.z .+ U_z_0[:,selectedMode] .* sf

    deformedMesh2 = deepcopy(mesh)
    deformedMesh2.x = deformedMesh2.x .+ U_x_90[:,selectedMode] .* sf
    deformedMesh2.y = deformedMesh2.y .+ U_y_90[:,selectedMode] .* sf
    deformedMesh2.z = deformedMesh2.z .+ U_z_90[:,selectedMode] .* sf


    # plot meshes
    fig = PyPlot.figure(selectedMode)
    ax1 = fig.add_subplot(221)#, projection="3d")
    ax2 = fig.add_subplot(222)#, projection="3d")
    ax3 = fig.add_subplot(223)#, projection="3d")
    ax4 = fig.add_subplot(224, projection="3d")

    plotMesh(deformedMesh,"-r",mesh.meshSeg,ax1,ax2,ax3,ax4)
    plotMesh(deformedMesh2,"--b",mesh.meshSeg,ax1,ax2,ax3,ax4)
    plotMesh(mesh,"-k",mesh.meshSeg,ax1,ax2,ax3,ax4)
    # ax1.view_init(0,0)
    # ax2.view_init(0,90)
    # ax3.view_init(90,0)
    ax4.view_init(45, 45)
    ax4.grid("off")
    # legend('in phase','out phase')
    # annotation('textbox', [0 0.9 1 0.1], 'String', ...
    # [strrep(resultsFile(1:end-16),'_','-') ' -- MODE: ' num2str(selectedMode) ' -- DOF: ' num2str((selectedMode-1)/2+1) ...
    # ' -- Freq: ' num2str(modalOut{selectedMode}.frequency) ' hz (' num2str(1/modalOut{selectedMode}.frequency,'#.3f') ' sec)'], ...
    # 'EdgeColor', 'none', 'HorizontalAlignment', 'center',...
    # 'FontSize',11)

end



function plotMesh(mesh,colorstring,meshSeg,ax1,ax2,ax3,ax4)
    #plotMesh plots a mesh with various orthographic and isometric views
    #    plotMesh(mesh,colorstring,meshSeg)
    #
    #   This function plots a mesh in a 3D with isometric view.
    #
    #      input:
    #      mesh            = object containing mesh data
    #      colorstring     = string containing color for plotting
    #      meshSeg         = array containing number of elements per mesh
    #                        segment
    #
    #
    #      output:          (NONE)

    conn = mesh.conn #assign mesh connectivity and coordinates
    x=mesh.x
    y=mesh.y
    z=mesh.z

    index=0
    for i=1:length(meshSeg) #loop over mesh segments
        p1 = zeros(meshSeg[i]*2)
        p2 = zeros(meshSeg[i]*2)
        p3 = zeros(meshSeg[i]*2)
        for j=1:meshSeg[i] #loop over elements in mesh segment i
            index=index+1
            n1=Int(conn[index,1]) #get node numbers for element j in mesh segment i
            n2=Int(conn[index,2])

            #get nodal coordinates for element j in mesh segment i
            p1[(j-1)*2 + 1] = x[n1]
            p1[(j-1)*2 + 2] = x[n2]
            p2[(j-1)*2 + 1] = y[n1]
            p2[(j-1)*2 + 2] = y[n2]
            p3[(j-1)*2 + 1] = z[n1]
            p3[(j-1)*2 + 2] = z[n2]
        end

        #plot meshes in various figure windows

        ax1.plot(p1,p3,colorstring)
        PyPlot.xlabel("h_1"); PyPlot.ylabel("h_3")#; PyPlot.zlabel("h_3")

        ax2.plot(p1,p2,colorstring)
        PyPlot.xlabel("h_1"); PyPlot.ylabel("h_2")#; PyPlot.zlabel("h_3")

        ax3.plot(p2,p3,colorstring)
        PyPlot.xlabel("h_2"); PyPlot.ylabel("h_3")#; PyPlot.zlabel("h_3")

        ax4.plot(p1,p2,p3,colorstring)
        PyPlot.xlabel("h_1"); PyPlot.ylabel("h_2"); PyPlot.zlabel("h_3")
    end

end
