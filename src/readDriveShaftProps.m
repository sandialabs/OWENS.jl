function [model] =  readDriveShaftProps(model,driveShaftFlag,dsfilename)
%readDriveShaftProps  reads driveshaft properties from file
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [model] =  readDriveShaftProps(model,driveShaftFlag,dsfilename)
%
%   This function reads driveshaft properties from file.
%
%   input:
%   model          = model object containing model information
%   driveShaftFlag = boolean flag for activating drive shaft effects
%   dsfilename     = string containing drive shaft property file name
%
%   output:
%   model          = model object containing model information

if(driveShaftFlag) %if drive shaft active, open file and read properties
    %     fid = fopen(dsfilename);
    %     a = fscanf(fid,'%f',3);     [dum] = fgetl(fid);
    %     b = fscanf(fid,'%f',2);     [dum] = fgetl(fid);
    %
    %     model.driveTrainOn = true;   %set drive shaft to active
    %
    %     model.driveShaftProps.k = a(1); %assign drive shaft stiffness
    %     model.driveShaftProps.c = a(2); %assign drive shaft damping
    %     model.JgearBox =a(3);           %assign gearbox MOI
    %
    %     model.gearRatio = b(1);         %assign gear ratio
    %     model.gearBoxEfficiency = b(2); %assign gear box efficiency
    error(['DRIVESHAFT NOT FULLY ENABLED while loading: ' dsfilename])
    
else   %if drive train is deactivated
    model.driveTrainOn = false;          %set drive shaft unactive
    
    model.driveShaftProps.k = 0.0;       %set drive shat properties to 0
    model.driveShaftProps.c = 0.0;
    model.JgearBox =0.0;
    
    model.gearRatio = 1.0;             %set gear ratio and efficiency to 1
    model.gearBoxEfficiency = 1.0;
end

end