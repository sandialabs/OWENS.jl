clear all
clc

cd('C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\modalAnalysis')

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
    bmOwens = 'platformTest';
    % append this name to the end of the saved files
    outFileExt = '_15mTowerExt_NOcentStiff';
else
    % use this benchmark file
    bmOwens = 'platformTest_10mTowerExtension';
    % append this name to the end of the saved files
    outFileExt = '_10mTowerExt_NOcentStiff';
end

platformProp = {};
% filename root to save the created nodal file
% % platformProp{end+1}.fileRoot = '1_FourColumnSemi_1stPass';
% % platformProp{end}.MassDiag = [9.1883e6 9.1883e6 1.3305e7 2.946e9 2.979e9 1.9471e9]';
% % platformProp{end}.StiffDiag = [1.329e5 1.329e5 2.585e6 4.132e6 4.132e6 1.076e6]';
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '1_FourColumnSemi_2ndPass';
platformProp{end}.MassDiag = [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9];
platformProp{end}.StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '2_ClassicSpar';
platformProp{end}.MassDiag = [3.0605e7 3.0605e7 1.7914e7 5.2104e10 5.2137e10 6.5356e8];
platformProp{end}.StiffDiag_Nm_deg = [1.348e5 1.348e5 2.327e6 1.181e7 1.181e7 8.189e5];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '3_RingPontoon';
platformProp{end}.MassDiag = [8.73e6 8.73e6 3.3648e7 5.6784e9 5.7114e9 9.2874e8];
platformProp{end}.StiffDiag_Nm_deg = [1.329e5 1.329e5 6.599e6 3.977e6 3.977e6 1.076e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '4_CompactSemi';
platformProp{end}.MassDiag = [8.5581e6 8.5582e6 7.9592e6 1.7936e9 1.8266e9 1.3048e9];
platformProp{end}.StiffDiag_Nm_deg = [1.592e5 1.592e5 3.826e6 3.485e6 3.485e6 1.654e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '5_AdvancedSpar';
platformProp{end}.MassDiag = [3.2512e7 3.2512e7 2.6111e7 5.275e10 5.2783e10 6.1514e8];
platformProp{end}.StiffDiag_Nm_deg = [1.38e5 1.38e5 2.362e6 9.115e6 9.115e6 1.384e6];
% filename root to save the created nodal file
platformProp{end+1}.fileRoot = '6_MultiCellularTLP';
platformProp{end}.MassDiag = [7.3879e6 7.3879e6 1.3414e7 2.3365e9 2.3695e9 1.374e9];
platformProp{end}.StiffDiag_Nm_deg = [1.623e5 1.623e5 7.135e7 2.34e8 2.34e8 1.118e6];

% convert rotational stiffness to N-m/rad
for rr = 1:length(platformProp)
    convRotStiff = [1 1 1 180/pi 180/pi 180/pi];
    platformProp{rr}.StiffDiag = platformProp{rr}.StiffDiag_Nm_deg .* convRotStiff;
end


%% ************************************************************************
% perform the modal simulations using OWENS
% *************************************************************************
for rr = 1:length(platformProp)
    % *********************************************************************
    % perform operations for the nodal file generation
    % *********************************************************************
    MassVal = platformProp{rr}.MassDiag;
    StiffVal = platformProp{rr}.StiffDiag;
    
    fname = [platformProp{rr}.fileRoot outFileExt];
    
    nodes = [1 1];
    cmkType = {'M6' 'K6'};
    for dd = 1:6;
        % set up mass matrix
        cmkValues{1}(dd,dd) = MassVal(dd);
        % set up stiffness matrix
        cmkValues{2}(dd,dd) = StiffVal(dd);
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
    rotSpdArrayRPM = 0:.1:10; % rpm
    rotSpdArrayHz = rotSpdArrayRPM ./ 60;
    centStiff = true;      % centripetal stiffening
    Nmodes = 16;  % number of modes to calculate
    
%     [freqMat{rr}] = campDiagramGen([fname '.owens'], [fname '_CampDiag'], rotSpdArrayHz, centStiff, Nmodes)
end


return
%% Plotting scripts -- Plot mode shapes

if 1
    outFileExt = '_15mTowerExt_NOcentStiff';
    bmOwens = 'platformTest';
    saveFolder = 'C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\images\'
else
    outFileExt = '_10mTowerExt_NOcentStiff';
    bmOwens = 'platformTest_10mTowerExtension';
    saveFolder = 'C:\data\OffshoreVAWT\SESPhaseI_OWENS\dvawt\carbon\dvawt_c_2_lcdt\images\10mTowerExt\'
end

savePlot = 0;

Nmodes = 16;
NperRevLines = 4;
minRPM = 0;
maxRPM = 8;
ratedRPM = 7.2;

close all
for rr=1:length(platformProp)
    
    outname = [platformProp{rr}.fileRoot outFileExt '_CampDiag.mat'];
    
    campDiagPlotter(outname, Nmodes, NperRevLines, minRPM, maxRPM)

    title(strrep(platformProp{rr}.fileRoot,'_',' '))
    xL = get(gca,'XLim')
    yL = get(gca,'YLim')
    hold on; plot([ratedRPM ratedRPM], [0 0.7], 'k:','LineWidth',1.5)
    
    if savePlot % save the plot
        export_fig([saveFolder 'CAMPDIAG_' outname(1:end-4) '.pdf'],'-p5')
        close gcf
    else % flip through the plots visually
%         pause
    end    
end










return
%% Read in the output files for comparison of mode frequencies, shape, etc.

Ndof = 40;
for ii=1:Ndof
    rowNamesHdr{ii} = num2str(ii);
end

% --------------- first data set ----------------- %
outFileExt = '_10mTowerExt_NOcentStiff';

for rr = 1:length(platformProp)
    resultsFile = [platformProp{rr}.fileRoot outFileExt '.out'];
    twr10m.(platformProp{rr}.fileRoot(3:end)) = readResultsModalOut(resultsFile,Ndof)
end

structData = twr10m;
names = fieldnames(structData);
tblTower10m = table(); tblTower10mPeriod = table();
for nn=1:length(names)   
    alldata = [structData.(names{nn}){:}]
    tblTower10m.(names{nn}) = [alldata.frequency]';
    tblTower10mPeriod.(names{nn}) = 1./[alldata.frequency]';
end
tblTower10m.Properties.RowNames = rowNamesHdr
tblTower10mPeriod.Properties.RowNames = rowNamesHdr

% --------------- second data set ----------------- %
outFileExt = '_15mTowerExt_NOcentStiff';

for rr = 1:length(platformProp)
    resultsFile = [platformProp{rr}.fileRoot outFileExt '.out'];
    twr15m.(platformProp{rr}.fileRoot(3:end)) = readResultsModalOut(resultsFile,Ndof)
end

structData = twr15m;
names = fieldnames(structData);
tblTower15m = table(); tblTower15mPeriod = table();
for nn=1:length(names)   
    alldata = [structData.(names{nn}){:}]
    tblTower15m.(names{nn}) = [alldata.frequency]';
    tblTower15mPeriod.(names{nn}) = 1./[alldata.frequency]';
end
tblTower15m.Properties.RowNames = rowNamesHdr
tblTower15mPeriod.Properties.RowNames = rowNamesHdr


% --------------- write tables to an excel file ----------------- %
filename = 'FrequencyComparison';

% write the first data set
writetable(tblTower10m,filename,'FileType','spreadsheet','Sheet','10m Tower Extension',...
    'WriteRowNames',true)
writetable(tblTower10mPeriod,filename,'FileType','spreadsheet','Sheet','10m Tower Extension - Period',...
    'WriteRowNames',true)
% write the second data set
writetable(tblTower15m,filename,'FileType','spreadsheet','Sheet','15m Tower Extension',...
    'WriteRowNames',true)
writetable(tblTower15mPeriod,filename,'FileType','spreadsheet','Sheet','15m Tower Extension - Period',...
    'WriteRowNames',true)


%% Analytical check of added mass properties and frequency calculations

% Assumptions: 
% * CG of rotor is directly above the CG of the platform (dx = dy = 22.5+5, dz = 0)

Mr = 568800;         % kg
IrollCG = 8.08e8;    % kg-m^2
IpitchCG = 8.41e8;   % kg-m^2
IyawCG = 3.48e7;     % kg-m^2
dCG = 22.5+5;        % m, distance from platform properties to rotor CG

% define the mass additions to the platform dynamics from the rotor
RotorMassMatrix = [Mr Mr Mr (IrollCG + Mr*dCG^2) (IpitchCG + Mr*dCG^2) IyawCG]

for rr = 1:length(platformProp)
    % create the combined system (platform + rotor) properties matrices
    systemProp{rr}.fileRoot = platformProp{rr}.fileRoot;
    systemProp{rr}.MassDiag = platformProp{rr}.MassDiag + RotorMassMatrix;
    systemProp{rr}.StiffDiag = platformProp{rr}.StiffDiag;
    % calculate the system natural periods and frequencies
    systemProp{rr}.natFreq = 1/(2*pi)*sqrt(systemProp{rr}.StiffDiag ./...
        systemProp{rr}.MassDiag)
    systemProp{rr}.natPeriod = 1 ./ systemProp{rr}.natFreq;
end

names = fieldnames(structData)

tblFreq = table(); tblPeriod = table();
for nn = 1:length(names)
    tblFreq.(names{nn}) = systemProp{nn}.natFreq';
    tblPeriod.(names{nn}) = systemProp{nn}.natPeriod';  
    rowNamesHdr2{nn} = num2str(nn);
end
tblFreq.Properties.RowNames = rowNamesHdr2;
tblPeriod.Properties.RowNames = rowNamesHdr2;

% save the natural freqency and periods
filename = 'ApproximatedNaturalPeriods_CombinedSystem';
writetable(tblFreq,filename,'FileType','spreadsheet','Sheet','Natural Frequency',...
    'WriteRowNames',true)
writetable(tblPeriod,filename,'FileType','spreadsheet','Sheet','Natural Period',...
    'WriteRowNames',true)


%% Analytical static analysis check

% Inputs
MaxRotorTorque = 7.29e3 * 1000;         % N-m
CutoutThrust = 1043 * 1000;             % N
CenterPressure = 76+5;                % m, center of pressure + 5m tower extension
MaxRPM = 7.2;                           % rpm


static = table(); rowHdrStr = {};
for rr = 1:length(platformProp)
    rowHdrStr{1} = 'Surge [m]';
    static.(names{rr})(1,1) = CutoutThrust / platformProp{rr}.StiffDiag_Nm_deg(1);
    rowHdrStr{2} = 'Pitch [deg]';
    static.(names{rr})(2,1) = CutoutThrust * CenterPressure / platformProp{rr}.StiffDiag_Nm_deg(5);
    rowHdrStr{3} = 'Yaw [deg]';
    static.(names{rr})(3,1) = MaxRotorTorque / platformProp{rr}.StiffDiag_Nm_deg(6);
end

static.Properties.RowNames = rowHdrStr

% save the natural freqency and periods
filename = 'StaticMotions_calculations';
writetable(static,filename,'FileType','spreadsheet','Sheet','Calculated',...
    'WriteRowNames',true)












