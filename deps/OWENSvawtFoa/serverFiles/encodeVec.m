function [workedFlag,type,output] = encodeVec(vec)
% as communication needs to have a strict format this vector this function
% will get the informations out of the sended vector depending on what
% reason it was send for. 

% first checksum
if sum(vec(1:end-1))==vec(end)
    workedFlag=true;
    
    type=vec(1);  % stating which type of transmission
    
    switch type
        % DO NOT CHANGE LEFT SIDE
        case 1  % interpolation data
            output = vec(2:end-1);  % spanwise coordinates
        case 2  % displacement statement
            output.time = vec(2);  %simulation time
            output.omega = vec(3); % rotational speed [rad/s]
            output.psi   = vec(4); % aszimuthal angle between p1 and h1 coordinates system
            elementsOfVec = vec(5);  % elements in one vector
            output.nOfBlades     = vec(6);  % amount of blades
            output.deformation = reshape(vec(7:end-1),elementsOfVec,4,output.nOfBlades);
            % creates 3D matrix every plane represents one blade
            % Displacement is a matrix with 4 columns
            % 1st xlocal displacement (h1 direction)
            % 2nd ylocal displ (h2 direction)
            % 3rd zlocal displacements
            % 4th pitching of airfoil
            
        case 3  % force
            elementsOfVec = vec(2);  % elements in one vector
            output.nOfBlades     = vec(3);  % amount of blades
            output.forces = reshape(vec(4:end-1),elementsOfVec,6,output.nOfBlades);
            % creates 3D matrix every plane represents one blade
            % Forces is a matrix 3D with 3 columns
            % 1st forcex (h1 direction)
            % 2nd frocey (h2 direction)
            % 3rd forcez (h3 direction)
            % 4th momentx (around h1)
            % 5th momenty (around h2)
            % 6th momentz (around h3)
        case 4
            output.time = vec(2);  % termination of simulation
        otherwise
            workedFlag = false;
            
    end
   
else
    workedFlag = false;
end
return;
