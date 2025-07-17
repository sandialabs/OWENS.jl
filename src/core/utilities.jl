"""
    simpleGenerator(generatorProps,genSpeed)

Caclulates generator torque for simple induction generator

#Input
* `generatorProps` object containing generator properties, see ?inputs
* `genSpeed::float`       generator speed (Hz)

#Output
* `genTorque::float`      generator torque
"""
function simpleGenerator(inputs,genSpeed)

    #assign generator properties form inputs object
    ratedTorque        = inputs.ratedTorque;
    ratedGenSlipPerc   = inputs.ratedGenSlipPerc;
    zeroTorqueGenSpeed = inputs.zeroTorqueGenSpeed;
    pulloutRatio       = inputs.pulloutRatio;

    #calculate rated generator speed
    ratedGenSpeed = zeroTorqueGenSpeed*(1.0 + 0.01*ratedGenSlipPerc);

    #calculate slope between lower and upper torque limits for simple induction generator
    midSlope = (ratedTorque/(ratedGenSpeed-zeroTorqueGenSpeed));

    #calculate lower and upper torque limits of generator
    upperTorqueLimit =ratedTorque*pulloutRatio;
    lowerTorqueLimit = -upperTorqueLimit;

    #calculate upper and lower generator speeds at which linear torque vs. speed region begins/ends
    upperGenSpeed = zeroTorqueGenSpeed + upperTorqueLimit/midSlope;
    lowerGenSpeed = zeroTorqueGenSpeed - upperTorqueLimit/midSlope;

    #calculate generator torque
    if genSpeed<0.0#lowerGenSpeed
        genTorque = 0.0#lowerTorqueLimit;
    elseif genSpeed>upperGenSpeed
        genTorque = upperTorqueLimit;
    else
        genTorque = midSlope*(genSpeed-zeroTorqueGenSpeed);
    end

    return genTorque

end


function safeakima(x,y,xpt;extrapolate=false)
    if minimum(xpt)<(minimum(x)-(abs(minimum(x))*0.1+1e-4)) || maximum(xpt)>(maximum(x)+abs(maximum(x))*0.1)
        msg="Extrapolating on akima spline results in undefined solutions minimum(xpt)<minimum(x) $(minimum(xpt))<$(minimum(x)) or maximum(xpt)>maximum(x) $(maximum(xpt))>$(maximum(x))"
        if !extrapolate
            throw(OverflowError(msg))
        else
            @warn msg
        end
    end
    return FLOWMath.akima(x,y,xpt)
end