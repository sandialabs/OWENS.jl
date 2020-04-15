function viz(meshFile,resultsFile,selectedMode,sf)
%viz  visualizes mode shapes
% **********************************************************************
% *                   Part of SNL VAWTGen                              *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   viz(meshFile,resultsFile,selectedMode,sf,meshSeg)
%
%   This function plots the mode shapes of a mode from a modal analysis
%   performed using the OWENS toolkit.
%
%      input:
%      meshFile         = string containing mesh file name
%      resultsFile      = string containing results file name
%      selectedMode     = integer denoting selected mode to plot
%      sf               = scale factor for mode shape displacements
%
%      output:          (NONE)

%     close all; %close all open figure windows %ble
figure %ble
[mesh,meshSeg] = readMeshVG(meshFile); %read mesh file

if 0 % original OWENS formulation
    %read in output file
    [modeShape,modeShapeOutPhase] = readResults(resultsFile,...
        selectedMode,mesh.numNodes);
    
    %set scale factor
    modeShape = modeShape .* sf;
    modeShapeOutPhase = modeShapeOutPhase.* sf;
    
    %add mode shape*scale factor + original components??
    deformedMesh = mesh;
    deformedMesh.x = deformedMesh.x + modeShape(:,1)';
    deformedMesh.y = deformedMesh.y + modeShape(:,2)';
    deformedMesh.z = deformedMesh.z + modeShape(:,3)';
    
    deformedMesh2 = mesh;
    deformedMesh2.x = deformedMesh2.x + modeShapeOutPhase(:,1)';
    deformedMesh2.y = deformedMesh2.y + modeShapeOutPhase(:,2)';
    deformedMesh2.z = deformedMesh2.z + modeShapeOutPhase(:,3)';
    
else % ble: changed formulation to plot frequency, etc    
    %read in output file
    modalOut = readResultsModalOut(resultsFile,mesh.numNodes);%ble: changed from readResults
    modeShape = modalOut{selectedMode}.InPhaseShape;
    modeShapeOutPhase = modalOut{selectedMode}.OutPhaseShape;
    
    %add mode shape*scale factor + original components??
    deformedMesh = mesh;
    deformedMesh.x = deformedMesh.x + modeShape.U_x' .* sf;
    deformedMesh.y = deformedMesh.y + modeShape.U_y' .* sf;
    deformedMesh.z = deformedMesh.z + modeShape.U_z' .* sf;
    
    deformedMesh2 = mesh;
    deformedMesh2.x = deformedMesh2.x + modeShapeOutPhase.U_x' .* sf;
    deformedMesh2.y = deformedMesh2.y + modeShapeOutPhase.U_y' .* sf;
    deformedMesh2.z = deformedMesh2.z + modeShapeOutPhase.U_z' .* sf;    
end

% plot meshes
plotMesh(deformedMesh,'-r',meshSeg);
plotMesh(deformedMesh2,'--b',meshSeg);
plotMesh(mesh,'-k',meshSeg);
% legend('in phase','out phase')
annotation('textbox', [0 0.9 1 0.1], 'String', ...
    [strrep(resultsFile(1:end-16),'_','-') ' -- MODE: ' num2str(selectedMode) ' -- DOF: ' num2str((selectedMode-1)/2+1) ...
    ' -- Freq: ' num2str(modalOut{selectedMode}.frequency) ' hz (' num2str(1/modalOut{selectedMode}.frequency,'%.3f') ' sec)'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center',...
    'FontSize',11)

end

