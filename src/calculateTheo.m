function [Theo] = calculateTheo(k)
%calculateTheo Calculates Theodorsen function value for unsteady aero
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Theo] = calculateTheo(k)
%
%   This function accepts a reduced frequency value and calculates the
%   complex Theodorsen function value.
%
%   input:
%   k     = reduced frequency
%
%   output:
%   Theo  = Theodorsen  constant value

Theo = 1.0 - 0.165/(1-0.0455/k*1i) - 0.335/(1-0.3/k*1i);

if(isnan(k) || isinf(k))
    Theo = 0.0;
end
end

