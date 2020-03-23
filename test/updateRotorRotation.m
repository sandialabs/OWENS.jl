function [azi_sp1,Omega_sp1,OmegaDot_sp1] = updateRotorRotation(Irotor,Crotor,Krotor,...
    shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s,...
    delta_t)
%updateRotorRotation updates rotor rotation
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
%   [azi_sp1,Omega_sp1,OmegaDot_sp1] = updateRotorRotation(Irotor,Crotor,Krotor,...
%                                  shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s,...
%                                  delta_t)
%
%   This function updates the rotor rotation given rotor properties and external
%   torques
%
%   input:
%   Irotor      = rotor inertia
%   Crotor      = arbitrary rotor damping
%   Krotor      = arbitrary rotor stiffness
%   shaftTorque = torque from external forces on rotor
%   genTorque   = torque from generator
%   azi_s       = rotor azimuth (rad) at beginning of time step
%   Omega_s     = rotor speed (Hz) at beginning of time step
%   OmegaDot_s  = rotor acceleration (Hz/s) at beginning of time step
%   delta_t     = time step
%
%   output:
%   azi_sp1       = rotor azimuth (rad) at end of time step
%   Omega_sp1     = rotor speed (Hz/s) at end of time step
%   OmegaDot_sp1  = rotor acceleration (Hz/s) at end of time step
%
Frotor = shaftTorque + genTorque; %calculate effective torque on rotor
Omega_s = Omega_s*2*pi; %conversion form Hz to rad/s, etc.
OmegaDot_s = OmegaDot_s*2*pi;
[azi_sp1,Omega_sp1,OmegaDot_sp1] = timeIntegrateSubSystem(Irotor,Krotor,Crotor,Frotor,... %time integrate using Newmark-Beta
    delta_t,azi_s,Omega_s,OmegaDot_s);

Omega_sp1 = Omega_sp1/(2*pi); %convert to Hz, etc.
OmegaDot_sp1 = OmegaDot_sp1/(2*pi);

end