function [elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ,nlflag)
% calculate strains

elStrain = struct('eps_xx_0',cell(1,numEl),'eps_xx_z',cell(1,numEl),'eps_xx_y',...
    cell(1,numEl),'gam_xz_0',cell(1,numEl),'gam_xz_y',cell(1,numEl),'gam_xy_0',cell(1,numEl),'gam_xy_z',cell(1,numEl));

for i=1:numEl
    %Calculate Ke and Fe for element i
    index = 1;
    elInput.elementOrder = elementOrder;
    elInput.nlOn = nlflag;
    elInput.xloc = [0.0 el.elLen(i)];
    elInput.sectionProps = el.props(i);
    elInput.sweepAngle = el.psi(i);
    elInput.coneAngle = el.theta(i);
    elInput.rollAngle = el.roll(i);
    elInput.aeroSweepAngle = 0.0;

    eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
    for j=1:numNodesPerEl       %define element coordinates and displacements associated with element
        for k=1:numDOFPerNode
            eldisp(index) = displ((conn(i,j)-1)*numDOFPerNode + k);
            index = index + 1;
        end
    end

    elInput.disp = eldisp;
    temp = calculateTimoshenkoElementStrain(elInput);

    elStrain(i).eps_xx_0 = temp.eps_xx_0;
    elStrain(i).eps_xx_z = temp.eps_xx_z;
    elStrain(i).eps_xx_y = temp.eps_xx_y;
    elStrain(i).gam_xz_0 = temp.gam_xz_0;
    elStrain(i).gam_xz_y = temp.gam_xz_y;
    elStrain(i).gam_xy_0 = temp.gam_xy_0;
    elStrain(i).gam_xy_z = temp.gam_xy_z;

end
end
