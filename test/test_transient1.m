
function test_transient1()


fprintf('%s\n','Starting')
tic

%% ****************** COORDINATE SYSTEM DEFINITION *********************** %
% 1 - x -  surge (pointed downstream)
% 2 - y -  sway  (right hand rule)
% 3 - z -  heave (pointed upwards)
% 4 - Ox - roll
% 5 - Oy - pitch
% 6 - Oz - yaw
% *********************************************************************** %

% use this benchmark file
bmOwens = './input_files_test/_15mTower_transient_dvawt_c_2_lcdt';
% append this name to the end of the saved files
outFileExt = '_15mTowerExt_NOcentStiff';

% convert rotational stiffness to N-m/rad
StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6];
convRotStiff = [1 1 1 180/pi 180/pi 180/pi];
StiffDiag = StiffDiag_Nm_deg .* convRotStiff;

% filename root to save the created nodal file
platformProp = struct('fileRoot','./input_files_test/1_FourColumnSemi_2ndPass',...
    'MassDiag',[9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9],...
    'StiffDiag_Nm_deg',[1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6],...
    'StiffDiag',StiffDiag);

%% ************************************************************************
% perform the transient simulations using OWENS
% *************************************************************************


% define the filename saving convention
fname = [platformProp.fileRoot outFileExt];

% *********************************************************************
% perform operations for the nodal file generation
% *********************************************************************
MassVal = platformProp.MassDiag;
StiffVal = platformProp.StiffDiag;

nodes = [1 1];
cmkType = {'M6' 'K6'};
cmkValues = zeros(6,6,2);
for dd = 1:6
    % set up mass matrix
    cmkValues(dd,dd,1) = MassVal(dd);
    % set up stiffness matrix
    cmkValues(dd,dd,2) = StiffVal(dd);
end
writeOwensNDL(fname, nodes, cmkType, cmkValues)

% *********************************************************************
% perform operations for the aerodynamic forces file generation
% *********************************************************************
CACTUSfileRoot = './input_files_test/DVAWT_2B_LCDT';
OWENSfileRoot = bmOwens;

% *********************************************************************
% run a modal analysis of the platform design
% *********************************************************************
operatingRPM = 7.2; % rpm
Nrpm = 10;    % number of rpm stations
Nmodes = 40;  % number of modes to calculate
timeStep = 2e-3;
timeSim = 0.1;       % [sec]
n_t = timeSim/timeStep; % length of time vector
timeArray = [0 timeSim+1];
rpmArray  = [operatingRPM operatingRPM];
omegaArrayHz = rpmArray ./ 60;
% EXAMPLE: owens(inputFile,'TNB',timeStep, nlBool, 0)
owens([fname '.owens'],'TNB', timeStep, floor(timeSim/timeStep), false, 0, timeArray, omegaArrayHz)

fprintf('%s\n','Function Finished')
end



