function [rotorAzimuth,rotorSpeed,rotorAcceleration] = getRotorPosSpeedAccelAtTime(t0,time,aziInit,delta_t)
%getRotorPosSpeedAccelAtTime uses user defined function to get rotor pos.
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [rotorAzimuth,rotorSpeed,rotorAcceleration] = getRotorPosSpeedAccelAtTime(t0,time,aziInit)
%
%   This function uses the user defined function rotorSpeedProfile() to get
%   the azimuth, speed, and acceleration of the rotor.
%
%   input:
%   t0      = time at which azimuth integration is beginning
%   time    = current time that position, velocity, and acceleration are
%             being requested
%   aziInit = initial rotor azimuth angle integration will begin at
%
%   output:
%   rotorAzimuth = azimuth position of rotor (rad) at time
%   rotorSpeed   = rotor speed (Hz) at time
%   rotorAcceleration = rotor acceleration (Hz/s) at time


rotorSpeed = userDefinedRotorSpeedProfile(time); %get rotgor speed at time

dt = 0.01;%some small delta t used in estimating rotor acceleration
if((time-dt) < 0)
    dt = delta_t/2;
end

omega_p1 = userDefinedRotorSpeedProfile(time+dt); %get rotor speed slightly before and after time
omega_m1 = userDefinedRotorSpeedProfile(time-dt);
%--------------------------------------------------------------------------

%estimate rotor acceleration with difference calculation
rotorAcceleration = diff([omega_m1,omega_p1])/(2*dt);

%calculate rotor azimuth using trapezoidal rule
[rotorAzimuth]= trapezoidalRule(aziInit,userDefinedRotorSpeedProfile(t0),rotorSpeed,time-t0);

end

%simple trapezoidal rule integration
function [aziEnd]= trapezoidalRule(aziInit,rotorSpeedStart,rotorSpeedEnd,dt)
aziEnd = aziInit + 0.5*dt*(rotorSpeedStart+rotorSpeedEnd)/(2*pi);
end