function [rom,elStorage]=reducedOrderModel(model,mesh,el,displ)
%reducedOrderModel   constructs a reduced order model
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freq,damp]=modalExec(model,mesh,el,displ,Omega,OmegaStart,fid)
%                    
%   This function executes modal analysis.
%
%   input:
%   model          = object containing model information
%   mesh           = object containing mesh information
%   el             = object containing element information
%   displ          = displacement vector for use in pre-stressed analysis
%
%   output:
%   rom            = object containing a reduced order model
%   elStorage      = object containing stored element matrices

[elStorage] = initialElementCalculations(model,el,mesh); %performs initial element calculations

[rom0] = calculateROM(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage);  %Omega   = 0; %calculates system matrices for parked condition

%calculates system matrices for various acceleration, rotor speed, rotor
%accelration combinations
omx = 1; omy = 1; omz = 1;
omxdot = 1; omydot = 1; omzdot = 1;
a_x = 1; a_y = 1; a_z = 1;
[rom1] = calculateROMGyric(model,mesh,el,displ,[omx 0.0 0.0],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_x = 1; %OmegaDot_i = 0; accel _i = 0
[rom2] = calculateROMGyric(model,mesh,el,displ,[0.0 omy 0.0],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_y = 1; %OmegaDot_i = 0; accel_i = 0
[rom3] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 omz],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_z = 1; %OmegaDot_i = 0; accel_i = 0
[rom4] = calculateROMGyric(model,mesh,el,displ,[omx omy 0.0],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_x,y = 1 %OmegaDot_i = 0; accel_i = 0
[rom5] = calculateROMGyric(model,mesh,el,displ,[0.0 omy omz],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_y,z = 1; %OmegaDot_i = 0; accel_i = 0
[rom6] = calculateROMGyric(model,mesh,el,displ,[omx 0.0 omz],[0.0 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_x,z = 1; %OmegaDot_i = 0; accel_i = 0
[rom7] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 0.0 0.0],[a_x 0.0 0.0],elStorage,rom0);  %Omega_i = 0; %OmegaDot_i = 0; accel_x = 1
[rom8] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 0.0 0.0],[0.0 a_y 0.0],elStorage,rom0);  %Omega_i = 0; %OmegaDot_i = 0; accel_y = 1
[rom9] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 0.0 0.0],[0.0 0.0 a_z],elStorage,rom0);  %Omega_i = 0; %OmegaDot_i = 0; accel_z = 1
[rom10] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[omxdot 0.0 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_i = 0; %OmegaDot_1 = 0; accel_i=0
[rom11] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 omydot 0.0],[0.0 0.0 0.0],elStorage,rom0);  %Omega_i = 0; %OmegaDot_2 = 0; accel_i=0
[rom12] = calculateROMGyric(model,mesh,el,displ,[0.0 0.0 0.0],[0.0 0.0 omzdot],[0.0 0.0 0.0],elStorage,rom0);  %Omega_i = 0; %OmegaDot_3 = 0; accel_i=0

%reduced order structural stiffness, mass, and damping
rom.Kr = rom0.Kr;
rom.Mr = rom0.Mr;
rom.Cr = rom0.Cr + rom0.CrModal;
rom.Fr = rom0.Fr;

%reduced order transformation matrices
rom.Phi = rom0.Phi;
rom.invPhi = rom0.invPhi;

%reduced order spin softening coefficient matrices
rom.SrOx2 = rom1.Kr;
rom.SrOy2 = rom2.Kr;
rom.SrOz2 = rom3.Kr;
rom.SrOxOy = rom4.Kr - rom1.Kr - rom2.Kr; 
rom.SrOyOz = rom5.Kr - rom2.Kr - rom3.Kr;
rom.SrOxOz = rom6.Kr - rom1.Kr - rom3.Kr;

%reduced order centrifugal load coefficient matrices
rom.FrOx2 = rom1.Fr;
rom.FrOy2 = rom2.Fr;
rom.FrOz2 = rom3.Fr;
rom.FrOxOy = rom4.Fr - rom1.Fr - rom2.Fr; 
rom.FrOyOz = rom5.Fr - rom2.Fr - rom3.Fr;
rom.FrOxOz = rom6.Fr - rom1.Fr - rom3.Fr;
rom.FrOxdot = rom10.Fr;
rom.FrOydot = rom11.Fr;
rom.FrOzdot = rom12.Fr;

%need to calculate reduced order acceleration/body force coefeficients
rom.FrAx = rom7.Fr;
rom.FrAy = rom8.Fr;
rom.FrAz = rom9.Fr;

%reduced order gyric coefficient matrices
rom.GrOx = rom1.Cr;
rom.GrOy = rom2.Cr;
rom.GrOz = rom3.Cr;

%reduced order circulatory coefficient matrices
rom.HrOx = 0.5*rom1.Cr;
rom.HrOy = 0.5*rom2.Cr;
rom.HrOz = 0.5*rom3.Cr;

end