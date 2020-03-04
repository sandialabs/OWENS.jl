function [unp1,udotnp1,uddotnp1] = timeIntegrateSubSystemEff(M,K,C,F,timeInt,u,udot,uddot)
%timeIntegrateSubSystemEff integrates a system using Newmark-Beta method
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%                    
%   [unp1,udotnp1,uddotnp1] = timeIntegrateSubSystemEff(M,K,C,F,timeInt,u,udot,uddot)
%
%   %This function perform integration of a system using the Newmark-Beta
%   method(constant-average acceleration sceheme). The integration
%   parameters are calculated before hand and store in the timeInt object.
%
%   input:
%   M        = system mass matrix
%   K        = system sttiffness matrix
%   C        = system damping matrix
%   F        = system force vector
%   timeInt  = object containing time integraton parameters
%   u        = displacement at beginning of time step
%   udot     = velocity at beginning of time step
%   uddot    = acceleration at beginning of time step
%
%
%   output:
%   unp1        = displacement at end of time step
%   udotnp1     = velocity at end of time step
%   uddotnp1    = acceleration at end of time step
%   
%transient integration of sub system using Newmark Beta method
    a1 = timeInt.a1;
    a2 = timeInt.a2;
    a3 = timeInt.a3;
    a4 = timeInt.a4;
    a5 = timeInt.a5;
    a6 = timeInt.a6;
    a7 = timeInt.a7;
    a8 = timeInt.a8;
    A = a3*u + a4*udot + a5*uddot;
    B = a6*u + a7*udot + a8*uddot;
    
    Khat = K + a3.*M + a6.*C;
    Fhat = F + M*(A) + C*(B);
    
    unp1 = Khat\Fhat;
    
    uddotnp1 = a3*(unp1-u) - a4*udot - a5*uddot;
    udotnp1 =  udot + a2*uddot + a1*uddotnp1;
        
end
