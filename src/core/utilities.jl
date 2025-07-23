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

"""
    setInitConditions(initDisps, numNodes, numDOFPerNode)

Creates the formatted initial conditions array needed by OWENSFEA

#Input
* `initDisps`: an array of length numDOFPerNode specifying the initial displacement of each DOF
* `numNodes`: the number of nodes in the given mesh
* `numDOFPerNode`: the number of unconstrained degrees of freedom calculated in each node

#Output
* `initCond`: array containing initial conditions.
    initCond(i,1) node number for init cond i.
    initCond(i,2) local DOF number for init cond i.
    initCond(i,3) value for init cond i.
"""
function createInitCondArray(initDisps, numNodes, numDOFPerNode)
    if initDisps == zeros(length(initDisps))
        initCond = []
    else
        initCond = zeros(numNodes*numDOFPerNode, 3)
        for i = 1:length(initDisps)
            if initDisps[i] != 0.0
                dof_initCond = hcat(collect(1:numNodes), ones(Int,numNodes)*i, ones(numNodes)*initDisps[i])
                initCond[(i-1)*numNodes+1:i*numNodes,:] = dof_initCond
            end
        end
        initCond = initCond[vec(mapslices(col -> any(col .!= 0), initCond, dims = 2)), :] #removes rows of all zeros
    end

    return initCond
end

function setBCs(fixedDOFs, fixedNodes, numNodes, numDOFPerNode) #node, dof, bc
    if (fixedDOFs == []) && (fixedNodes == [])
        pBC = []
    else
        pBC = zeros(Int, numNodes*numDOFPerNode, 3)
        for i = 1:length(fixedNodes)
            pBC[(i-1)*numDOFPerNode+1:i*numDOFPerNode, :] = hcat(ones(numDOFPerNode)*fixedNodes[i], collect(1:numDOFPerNode), zeros(Int, numDOFPerNode) )
        end
        for i = 1:length(fixedDOFs)
            newNodes = setdiff(1:numNodes, fixedNodes) # this avoids duplicating nodes already counted for by fixedNodes
            numNewNodes = length(newNodes)
            dofBCs = hcat(newNodes, ones(Int,numNewNodes)*fixedDOFs[i], zeros(Int,numNewNodes))
            pBC[length(fixedNodes)*numDOFPerNode+(i-1)*numNewNodes+1:length(fixedNodes)*numDOFPerNode+i*numNewNodes, :] = dofBCs
        end
        pBC = pBC[vec(mapslices(col -> any(col .!= 0), pBC, dims = 2)), :] #removes extra rows (i.e. rows of all zeros)
    end

    return pBC

end