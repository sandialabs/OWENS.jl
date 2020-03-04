function [elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ,nlflag)
% calculate strains
           for i=1:numEl
            %Calculate Ke and Fe for element i
            index = 1;
            elInput.elementOrder = elementOrder;
            elInput.nlOn = nlflag;
            elInput.xloc = [0.0 el.elLen(i)];
            elInput.sectionProps = el.props{i};
            elInput.sweepAngle = el.psi(i);
            elInput.coneAngle = el.theta(i);
            elInput.rollAngle = el.roll(i);
            elInput.aeroSweepAngle = 0.0;

            for j=1:numNodesPerEl       %define element coordinates and displacements associated with element
                for k=1:numDOFPerNode
                    eldisp(index) = displ((conn(i,j)-1)*numDOFPerNode + k);
                    index = index + 1;
                end
            end

            elInput.disp = eldisp;
            [elStrain(i)] = calculateTimoshenkoElementStrain(elInput);
            
           end
end
