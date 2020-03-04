function vec = decodeVec(type,varargin)
%as communication needs to have a strict format this vector this function
% will get the informations into the sended vector depending on what
% reason it was send for. 
vec(1) = type;  % first digit of vector

switch type
        case 1  % interpolation data
            vec = [vec varargin{1}];  % spanwise coordinates
        case 2  % displacement statement
             vec = [vec varargin{1}];  %simulation time
             vec = [vec varargin{2}]; % rotational speed [rad/s]
             vec = [vec varargin{3}]; % aszimuthal angle between p1 and h1 coordinates system
             vec = [vec varargin{4}]; % elements in one vector
             vec = [vec varargin{5}]; % amount of blades
             % per blade
             for bb = 1:varargin{5}
             vec = [vec varargin{5+bb}(:)']; % displacements of blade bb
            % Displacement is a matrix with 4 columns
            % 1st xlocal displacement (h1 direction)
            % 2nd ylocal displ (h2 direction)
            % 3rd zlocal displacements
            % 4th pitching of airfoil
             end
            
        case 3  % force
             vec = [vec varargin{1}];  % elements in one vector
             vec = [vec varargin{2}];  % amount of blades
             % per blade
             for bb = 1:varargin{2}
                 forceBlade = varargin{3}(:,:,bb);
                 vec = [vec forceBlade(:)']; % displacements of blade bb
            % Forces is a matrix 3D with 6 columns
            % 1st forcex (h1 direction)
            % 2nd frocey (h2 direction)
            % 3rd forcez (h3 direction)
            % 4th momentx (around h1)
            % 5th momenty (around h2)
            % 6th momentz (around h3)
             end
        case 4
            vec = [vec varargin{1}];
    otherwise
           warning('Server:encode','Type of information not recognized');          
end
vec = [vec sum(vec)];
return;


