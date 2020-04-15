function [genTorque] = simpleGenerator(generatorProps,genSpeed)
%simpleGenerator   calculates generator torque
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [genTorque] = simpleGenerator(generatorProps,genSpeed)
%
%   This function caclulates generator torque for simple induction
%   generator
%
%   input:
%   generatorProps = object containing generator properties
%   genSpeed       = generator speed (Hz)
%
%   output:
%   genTorque      = generator torque

%assign generator properties form generatorProps object
ratedTorque        = generatorProps.ratedTorque;
ratedGenSlipPerc   = generatorProps.ratedGenSlipPerc;
zeroTorqueGenSpeed = generatorProps.zeroTorqueGenSpeed;
pulloutRatio       = generatorProps.pulloutRatio;

%calculate rated generator speed
ratedGenSpeed = zeroTorqueGenSpeed*(1.0 + 0.01*ratedGenSlipPerc);

%calculate slope between lower and upper torque limits for simple induction
%generator
midSlope = (ratedTorque/(ratedGenSpeed-zeroTorqueGenSpeed));

%calculate lower and upper torque limits of generator
upperTorqueLimit =ratedTorque*pulloutRatio;
lowerTorqueLimit = -upperTorqueLimit;

%calculate upper and lower generator speeds at which linear torque vs.
%speed regin begins/ends
upperGenSpeed = zeroTorqueGenSpeed + upperTorqueLimit/midSlope;
lowerGenSpeed = zeroTorqueGenSpeed - upperTorqueLimit/midSlope;

%calculate generator torque
if(genSpeed<lowerGenSpeed)
    genTorque = lowerTorqueLimit;
elseif(genSpeed>upperGenSpeed)
    genTorque = upperTorqueLimit;
else
    genTorque = midSlope*(genSpeed-zeroTorqueGenSpeed);
end

end

