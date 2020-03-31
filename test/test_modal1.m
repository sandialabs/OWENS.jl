function test_modal1()

fprintf('%s\n','Starting Modal')
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
% perform the modal simulations using OWENS
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
% read in the main owens file
% *********************************************************************
fid = fopen([bmOwens '.owens']);
owensMain = textscan(fid, '%s','whitespace','\r');
fclose(fid)
% change the reference for the nodal input file
owensMain{1}{5} = [fname '.ndl'];

% save the owens input file
fid = fopen([fname '.owens'],'w');
fprintf(fid,'%s\n',owensMain{1}{:});
fclose(fid);

% *********************************************************************
% run a modal analysis of the platform design
% *********************************************************************
maxRPM = 10;%7.2; % rpm
Nrpm = 10;    % number of rpm stations
Nmodes = 40;  % number of modes to calculate

owens([fname '.owens'],'M', 0.5*maxRPM*2*pi/60, false, Nmodes)


fprintf('%s\n','Function Finished')

end










