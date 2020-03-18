function [hydromodel] = readPlatformFile(hydromodel,platformFlag,platfilename)
%readPlatformFile reads platform file
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [model] = readPlatformFile(model,platformFlag,platfilename)
%
%   This function reads the platform input file and stores data in the
%   model object.
%
%      input:
%      model          = string containing nodal terms filename
%      platformFlag   = integer containing a flag to activate or deactivate
%                       platform simulation
%      platfilename   = string containing platform filename
%
%      output:
%      model          = object containing model data

if(platformFlag == 1)
    hydromodel.hydroOn = true;
elseif(platformFlag == 0)
    hydromodel.hydroOn = false;
else
    error('Platform translation flag not recognized.');
end

hydromodel.platformTurbineConnectionNodeNumber = 1; %default to node number 1, update with platform file
if(platformFlag == 1)
    
    fid = fopen(platfilename);
    
    if(fid>0) %check if file was opened successfully
        %             model.activePlatformDofs = fscanf(fid,'%i',6)';         [dum]=fgetl(fid);
        %             model.initialCondPlatformDofs = fscanf(fid,'%f',6)';    [dum]=fgetl(fid);
        %             model.Plat_Drag_Flag = fscanf(fid,'%i',1);    [dum]=fgetl(fid);
        %             model.Plat_Moor_Flag = fscanf(fid,'%i',1);    [dum]=fgetl(fid);
        %             model.Plat_Grav_Flag = fscanf(fid,'%i',1);    [dum]=fgetl(fid);
        %             model.Plat_Plot_Flag = fscanf(fid,'%i',1);    [dum]=fgetl(fid);
        %             model.Plat_RadDamp_Flag = fscanf(fid,'%i',1); [dum]=fgetl(fid);
        %             model.platformTurbineConnectionNodeNumber = fscanf(fid,'%i',1); [dum]=fgetl(fid);
        %             model.platformTurbineYawInteraction = fscanf(fid,'%i',1); [dum]=fgetl(fid);
        %             model.platformServerPort = fscanf(fid,'%i',1); [dum]=fgetl(fid);
        %             model.platformClientPort = fscanf(fid,'%i',1); [dum]=fgetl(fid);
        error('PLATFORM NOT FULLY ENABLED')
        
    else
        error('Platform file could not be opened.');
    end
end
end