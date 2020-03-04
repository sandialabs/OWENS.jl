% set up paths 

% add the main directory of OWENS
VAWT_Toolbox_path_main = 'C:\DesignCodes\OWENSvawtFoa\';
addpath(VAWT_Toolbox_path_main)
% add sub folders of the OWENS directory
VAWT_Toolbox_path_1 = [VAWT_Toolbox_path_main 'commonSource'];
addpath(VAWT_Toolbox_path_1)
VAWT_Toolbox_path_2 = [VAWT_Toolbox_path_main 'modalSource'];
addpath(VAWT_Toolbox_path_2)
VAWT_Toolbox_path_3 = [VAWT_Toolbox_path_main 'transientSource'];
addpath(VAWT_Toolbox_path_3)
VAWT_Toolbox_path_4 = [VAWT_Toolbox_path_main 'utilitySource'];
addpath(VAWT_Toolbox_path_4)
VAWT_Toolbox_path_5 = [VAWT_Toolbox_path_main 'serverFiles'];
addpath(VAWT_Toolbox_path_5)
VAWT_Toolbox_path_6 = [VAWT_Toolbox_path_main 'processingScripts'];
addpath(VAWT_Toolbox_path_6)

% add the main directory of VAWTgen
VAWT_Toolbox_path_main = 'C:\DesignCodes\vawtGen_031714\';
addpath(VAWT_Toolbox_path_main)
% add sub folders of the OWENS directory
VAWT_Toolbox_path_1 = [VAWT_Toolbox_path_main 'vizFiles'];
addpath(VAWT_Toolbox_path_1)

% add the matlab file exhange folder
MATLAB_Toolbox_path = '\\snl\home\blennis\MATLAB';
addpath(genpath(MATLAB_Toolbox_path))

% print to screen completion of path generation
disp('VAWT Tools Design Code path setup script complete.  (including OWENS and VAWTGen)')