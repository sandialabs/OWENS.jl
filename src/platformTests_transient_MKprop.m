clear all
% clc

% snl_vawt

% cd('C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\transientAnalysis')

if (~isdeployed)
    VAWT_Toolbox_path_main = '../deps/OWENSvawtFoa/';
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
%     VAWT_Toolbox_path_5 = [VAWT_Toolbox_path_main 'serverFiles'];
%     addpath(VAWT_Toolbox_path_5)
    VAWT_Toolbox_path_6 = [VAWT_Toolbox_path_main 'processingScripts'];
    addpath(VAWT_Toolbox_path_6)

    % add the main directory of VAWTgen
    % add sub folders of the OWENS directory
    VAWT_Toolbox_path_1 = [VAWT_Toolbox_path_main '../vizFiles'];
    addpath(VAWT_Toolbox_path_1)
end

%% ****************** COORDINATE SYSTEM DEFINITION *********************** %
% 1 - x -  surge (pointed downstream)
% 2 - y -  sway  (right hand rule)
% 3 - z -  heave (pointed upwards)
% 4 - Ox - roll
% 5 - Oy - pitch
% 6 - Oz - yaw
% *********************************************************************** %

if 1
    % use this benchmark file
    bmOwens = './input_files/_15mTower_transient_dvawt_c_2_lcdt';
    % append this name to the end of the saved files
    outFileExt = '_15mTowerExt_NOcentStiff';
else
    % use this benchmark file
    bmOwens = '_10mTower_transient_dvawt_c_2_lcdt';
    % append this name to the end of the saved files
    outFileExt = '_10mTowerExt_NOcentStiff';
end

platformProp = {};
% filename root to save the created nodal file
% % platformProp{end+1}.fileRoot = '1_FourColumnSemi_1stPass';
% % platformProp{end}.MassDiag = [9.1883e6 9.1883e6 1.3305e7 2.946e9 2.979e9 1.9471e9]';
% % platformProp{end}.StiffDiag = [1.329e5 1.329e5 2.585e6 4.132e6 4.132e6 1.076e6]';
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/1_FourColumnSemi_2ndPass';
platformProp{end}.MassDiag = [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9];
platformProp{end}.StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/2_ClassicSpar';
platformProp{end}.MassDiag = [3.0605e7 3.0605e7 1.7914e7 5.2104e10 5.2137e10 6.5356e8];
platformProp{end}.StiffDiag_Nm_deg = [1.348e5 1.348e5 2.327e6 1.181e7 1.181e7 8.189e5];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/3_RingPontoon';
platformProp{end}.MassDiag = [8.73e6 8.73e6 3.3648e7 5.6784e9 5.7114e9 9.2874e8];
platformProp{end}.StiffDiag_Nm_deg = [1.329e5 1.329e5 6.599e6 3.977e6 3.977e6 1.076e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/4_CompactSemi';
platformProp{end}.MassDiag = [8.5581e6 8.5582e6 7.9592e6 1.7936e9 1.8266e9 1.3048e9];
platformProp{end}.StiffDiag_Nm_deg = [1.592e5 1.592e5 3.826e6 3.485e6 3.485e6 1.654e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/5_AdvancedSpar';
platformProp{end}.MassDiag = [3.2512e7 3.2512e7 2.6111e7 5.275e10 5.2783e10 6.1514e8];
platformProp{end}.StiffDiag_Nm_deg = [1.38e5 1.38e5 2.362e6 9.115e6 9.115e6 1.384e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = './input_files/6_MultiCellularTLP';
platformProp{end}.MassDiag = [7.3879e6 7.3879e6 1.3414e7 2.3365e9 2.3695e9 1.374e9];
platformProp{end}.StiffDiag_Nm_deg = [1.623e5 1.623e5 7.135e7 2.34e8 2.34e8 1.118e6];

% convert rotational stiffness to N-m/rad
for rr = 1:length(platformProp)
    convRotStiff = [1 1 1 180/pi 180/pi 180/pi];
    platformProp{rr}.StiffDiag = platformProp{rr}.StiffDiag_Nm_deg .* convRotStiff;
end


%% ************************************************************************
% perform the transient simulations using OWENS
% *************************************************************************
for rr = 1%:length(platformProp)

    % define the filename saving convention
    fname = [platformProp{rr}.fileRoot outFileExt];

    % *********************************************************************
    % perform operations for the nodal file generation
    % *********************************************************************
    MassVal = platformProp{rr}.MassDiag;
    StiffVal = platformProp{rr}.StiffDiag;

    nodes = [1 1];
    cmkType = {'M6' 'K6'};
    for dd = 1:6
        % set up mass matrix
        cmkValues{1}(dd,dd) = MassVal(dd);
        % set up stiffness matrix
        cmkValues{2}(dd,dd) = StiffVal(dd);
    end
    writeOwensNDL(fname, nodes, cmkType, cmkValues)

    % *********************************************************************
    % perform operations for the aerodynamic forces file generation
    % *********************************************************************
    CACTUSfileRoot = './input_files/DVAWT_2B_LCDT';
    OWENSfileRoot = bmOwens;
    outputAeroFileName = './output_files/aa_TESTaero';

    processAeroLoadsBLE(CACTUSfileRoot, OWENSfileRoot, outputAeroFileName)


    % *********************************************************************
    % read in the main owens file
    % *********************************************************************
    fid = fopen([bmOwens '.owens']);
    owensMain = textscan(fid, '%s','whitespace','\r');
    fclose(fid);

    % save the owens input file
    fid = fopen([fname '.owens'],'w');
    fprintf(fid,'%s\n',owensMain{1}{:});
    fclose(fid);

    % *********************************************************************
    % run a modal analysis of the platform design
    % *********************************************************************
    operatingRPM = 7.2; % rpm
    Nrpm = 10;    % number of rpm stations
    Nmodes = 40;  % number of modes to calculate
    timeStep = 2e-3;
    timeSim = 0.2;       % [sec]
    timeArray = [0 timeSim+1];
    rpmArray  = [operatingRPM operatingRPM];
    omegaArrayHz = rpmArray ./ 60;
    % EXAMPLE: owens(inputFile,'TNB',timeStep, nlBool, 0)
    owens([fname '.owens'],'TNB', timeStep, floor(timeSim/timeStep), false, 0, timeArray, omegaArrayHz)
end


return
%% Plotting scripts -- Save the transient video

if 0
    % use this benchmark file
    bmOwens = '_15mTower_transient_dvawt_c_2_lcdt';
    % append this name to the end of the saved files
    outFileExt = '_15mTowerExt_NOcentStiff';
    saveFolder = 'C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\images\'
else
    % use this benchmark file
    bmOwens = '_10mTower_transient_dvawt_c_2_lcdt';
    % append this name to the end of the saved files
    outFileExt = '_10mTowerExt_NOcentStiff';
    saveFolder = 'C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\images\10mTowerExt\'
end



for rr=1%:length(platformProp)

    outname = [platformProp{rr}.fileRoot outFileExt '.out'];
    out = load([outname(1:end-4) '.mat'])

    vizAnimateTransient([bmOwens '.mesh'], out.uHist, 5, outname(1:end-4))
    %             if savePlot % save the plot
    %                 export_fig([saveFolder outname(1:end-4) '_MODE' num2str(df) '.pdf'],'-p5')
    %                 close gcf
    %             else % flip through the plots visually
    %                 pause
    %             end


end


%% Plot transient output

if ~isempty(out)
    out = load('1_FourColumnSemi_2ndPass_10mTowerExt_NOcentStiff.mat')
end

i1 = 1;
i2 = round(length(out.t));

figure
xVar = 'aziHist';
yVar = 'FReactionHist'; FRvar = 6;
plot(out.(xVar)(i1:i2), out.(yVar)(i1:i2,FRvar))
grid on; xlim([0 4*pi])

figure
xVar = 't';
yVar = 'genPower';
plot(out.(xVar)(i1:i2), out.(yVar)(i1:i2))
