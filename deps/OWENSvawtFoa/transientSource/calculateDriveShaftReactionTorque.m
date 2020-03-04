function [torque] = calculateDriveShaftReactionTorque(driveShaftProps,thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)
%calculateDriveShaftReactionTorque calculates reaction torque of driveshaft
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [torque] = calculateDriveShaftReactionTorque(driveShaftProps,...
%                thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)
%                    
%   This function calculates reaction torque of driveshaft
%
%   input:
%   driveShaftProps      = object containing driveshaft properties
%   thetaRotor           = azimuth position of rotor/rotor shaft (rad)
%   thetaGB              = azimuth position of gearbox shaft (rad)
%   thetaDotRotor        = angular velocity of rotor/rotor shaft (rad/s)
%   thetaDotGB           = angular velocity of gearbox shaft (rad/s)
%
%   output:     
%   torque   = reaction torque of drive shaft


    k = driveShaftProps.k;  %drive shaft stiffness
    c = driveShaftProps.c;  %drive shaft damping

    torque = k*(thetaRotor-thetaGB) + c*(thetaDotRotor-thetaDotGB);

end