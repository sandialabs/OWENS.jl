function [Omega] = userDefinedRotorSpeedProfile(~)
%userDefinedRotorSpeedProfile  returns specified rotor speed
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Omega] = userDefinedRotorSpeedProfile(time)
%
%   This function returns an arbitrary time-varying rotor speed profile
%   specified by the user. Rotor speed must be specified in Hz.
%
%   input:
%   time           = time
%
%   output:
%   Omega          = rotor speed (Hz)

Omega = 0.5;

end