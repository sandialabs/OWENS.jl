
function test_transient1()

tol1 = 1e-12;
if (~isdeployed)
%     clear all
%     
%     
%     VAWT_Toolbox_path_main = '../deps/OWENSvawtFoa/';
%     addpath(VAWT_Toolbox_path_main)
%     % add sub folders of the OWENS directory
%     VAWT_Toolbox_path_1 = [VAWT_Toolbox_path_main 'commonSource'];
%     addpath(VAWT_Toolbox_path_1)
%     VAWT_Toolbox_path_2 = [VAWT_Toolbox_path_main 'modalSource'];
%     addpath(VAWT_Toolbox_path_2)
%     VAWT_Toolbox_path_3 = [VAWT_Toolbox_path_main 'transientSource'];
%     addpath(VAWT_Toolbox_path_3)
%     VAWT_Toolbox_path_4 = [VAWT_Toolbox_path_main 'utilitySource'];
%     addpath(VAWT_Toolbox_path_4)
%     %     VAWT_Toolbox_path_5 = [VAWT_Toolbox_path_main 'serverFiles'];
%     %     addpath(VAWT_Toolbox_path_5)
%     VAWT_Toolbox_path_6 = [VAWT_Toolbox_path_main 'processingScripts'];
%     addpath(VAWT_Toolbox_path_6)
%     
%     % add the main directory of VAWTgen
%     % add sub folders of the OWENS directory
%     VAWT_Toolbox_path_1 = [VAWT_Toolbox_path_main '../vizFiles'];
%     addpath(VAWT_Toolbox_path_1)
end

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

processAeroLoadsBLE(CACTUSfileRoot, OWENSfileRoot)

% *********************************************************************
% run a modal analysis of the platform design
% *********************************************************************
operatingRPM = 7.2; % rpm
Nrpm = 10;    % number of rpm stations
Nmodes = 40;  % number of modes to calculate
timeStep = 2e-3;
timeSim = 0.1;       % [sec]
timeArray = [0 timeSim+1];
rpmArray  = [operatingRPM operatingRPM];
omegaArrayHz = rpmArray ./ 60;
% EXAMPLE: owens(inputFile,'TNB',timeStep, nlBool, 0)
owens([fname '.owens'],'TNB', timeStep, floor(timeSim/timeStep), false, 0, timeArray, omegaArrayHz)


if (~isdeployed)
    
    OLD = ('1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_OLD.mat');
    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.mat');
    
    FileInfo = dir(NEW);
    
    if (datetime - FileInfo.date) > duration(0,1,0)
        error('Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.');
    end
    
    old = load(OLD);
    new = load(NEW);
    
    total_mismatch = 0;
    
    varnames = fieldnames(old);
    for i = 1:length(varnames)
        
        num_mismatch = 0;
        old_data = old.(varnames{i});
        new_data = new.(varnames{i});
        
        if ~isnumeric(old_data)
            subvarnames = fieldnames(old_data);
            for j = 1:length(subvarnames)
                
                sub_old_data = old_data.(subvarnames{j});
                sub_new_data = new_data.(subvarnames{j});
                for ii = 1:length(sub_old_data(:,1))
                    for jj = 1:length(sub_old_data(1,:))
                        if (abs(sub_old_data(ii,jj)-sub_new_data(ii,jj))>tol1)
                            msg = sprintf('%20s%i%2s%i%6s%8e%6s%8e', subvarnames{j}  , ii , ':' , jj ,' OLD: ' , sub_old_data(ii,jj) , ' NEW: ' , sub_new_data(ii,jj));
                            fprintf('%s\n',msg)
                            num_mismatch = num_mismatch + 1;
                        end
                    end
                end
            end
            
        elseif ~(min(min(ismembertol(old_data,new_data,tol1))))
            for ii = 1:length(old_data(:,1))
                for jj = 1:length(old_data(1,:))
                    if (abs(old_data(ii,jj)-new_data(ii,jj))>tol1)
                        msg = sprintf('%20s%i%2s%i%6s%8e%6s%8e',varnames{i} , ii , ':' , jj ,' OLD: ' , old_data(ii,jj) , ' NEW: ' , new_data(ii,jj));
                        fprintf('%s\n',msg)
                        num_mismatch = num_mismatch + 1;
                    end
                end
            end
            
        end
        
        if num_mismatch == 0
            fprintf('%s\n',['PASSED | ' , varnames{i}])
        else
            fprintf('%s\n',['FAILED | ', varnames{i}])
        end
        
        total_mismatch = total_mismatch + num_mismatch;
        
    end
    
    toc
    
    if total_mismatch == 0
        fprintf('%s\n','TESTS PASSED')
    else
        fprintf('%s\n','!!!TESTS FAILED!!!')
    end
    
end

end
