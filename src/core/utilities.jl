"""
    simpleGenerator(generatorProps,genSpeed)

Caclulates generator torque for simple induction generator

#Input
* `generatorProps` object containing generator properties, see ?inputs
* `genSpeed::float`       generator speed (Hz)

#Output
* `genTorque::float`      generator torque
"""
function simpleGenerator(inputs, genSpeed)

    #assign generator properties form inputs object
    ratedTorque = inputs.ratedTorque;
    ratedGenSlipPerc = inputs.ratedGenSlipPerc;
    zeroTorqueGenSpeed = inputs.zeroTorqueGenSpeed;
    pulloutRatio = inputs.pulloutRatio;

    #calculate rated generator speed
    ratedGenSpeed = zeroTorqueGenSpeed*(1.0 + 0.01*ratedGenSlipPerc);

    #calculate slope between lower and upper torque limits for simple induction generator
    midSlope = (ratedTorque/(ratedGenSpeed-zeroTorqueGenSpeed));

    #calculate lower and upper torque limits of generator
    upperTorqueLimit = ratedTorque*pulloutRatio;
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

"""
    hydrodynInputWithResolvedPotFile(hd_input_file, potflowfile)

Return a HydroDyn input filename whose `PotFile` entry is absolute.

HydroDyn resolves the `PotFile` line relative to the Julia process working
directory, not necessarily relative to the HydroDyn input file. OWENS examples
often launch from different directories, so this helper stages a temporary copy
of the input file with the `PotFile` token replaced by an absolute root. An
explicit `potflowfile` argument takes precedence. When `potflowfile` is missing,
the existing HydroDyn `PotFile` value is resolved from the input file itself.
"""
function hydrodynInputWithResolvedPotFile(hd_input_file, potflowfile)
    hd_input_file_string = string(hd_input_file)
    if isnothing(hd_input_file) ||
       lowercase(strip(hd_input_file_string)) == "none" ||
       !isfile(hd_input_file_string)
        return hd_input_file
    end

    lines = readlines(hd_input_file_string, keep = true)
    potfile_line_index = findfirst(
        line -> occursin(r"\bPotFile\b", line) && !startswith(strip(line), "!"),
        lines,
    )
    isnothing(potfile_line_index) && return hd_input_file

    potfile_line = lines[potfile_line_index]
    potfile_root =
        _resolvedHydroDynPotFileRoot(potfile_line, hd_input_file_string, potflowfile)
    isnothing(potfile_root) && return hd_input_file

    lines[potfile_line_index] = _hydrodynPotFileLineWithRoot(potfile_line, potfile_root)

    staged_file = tempname() * ".dat"
    open(staged_file, "w") do io
        foreach(line -> write(io, line), lines)
    end
    return staged_file
end

function _resolvedHydroDynPotFileRoot(potfile_line, hd_input_file_string, potflowfile)
    potflowfile_string = isnothing(potflowfile) ? "" : strip(string(potflowfile))
    if !isempty(potflowfile_string) &&
       !(lowercase(potflowfile_string) in ("none", "nothing"))
        return abspath(potflowfile_string)
    end

    potfile_root = _hydrodynPotFileRootFromLine(potfile_line)
    if isempty(potfile_root) || lowercase(potfile_root) in ("none", "nothing", "unused")
        return nothing
    end
    isabspath(potfile_root) && return abspath(potfile_root)

    hd_input_abspath = abspath(hd_input_file_string)
    candidates = unique(
        abspath.([
            potfile_root,
            joinpath(dirname(hd_input_abspath), potfile_root),
            joinpath(dirname(dirname(hd_input_abspath)), potfile_root),
        ]),
    )
    existing_root = findfirst(_hydrodynPotentialFlowRootExists, candidates)
    !isnothing(existing_root) && return candidates[existing_root]

    return abspath(joinpath(dirname(hd_input_abspath), potfile_root))
end

function _hydrodynPotFileRootFromLine(potfile_line)
    line_body = chomp(potfile_line)
    root_match = match(r"^\s*(?:\"([^\"]*)\"|(\S+))(?=\s+PotFile\b)", line_body)
    if isnothing(root_match)
        throw(ArgumentError("could not parse HydroDyn PotFile line: $(strip(line_body))"))
    end
    root =
        isnothing(root_match.captures[1]) ? root_match.captures[2] : root_match.captures[1]
    return strip(root)
end

function _hydrodynPotFileLineWithRoot(potfile_line, potfile_root)
    line_ending =
        endswith(potfile_line, "\r\n") ? "\r\n" : endswith(potfile_line, "\n") ? "\n" : ""
    line_body = chomp(potfile_line)
    line_match = match(r"^(\s*)(?:\"[^\"]*\"|\S+)(\s+PotFile\b.*)$", line_body)
    if isnothing(line_match)
        throw(ArgumentError("could not parse HydroDyn PotFile line: $(strip(line_body))"))
    end
    return string(
        line_match.captures[1],
        "\"$(potfile_root)\"",
        line_match.captures[2],
        line_ending,
    )
end

function _hydrodynPotentialFlowRootExists(potfile_root)
    return any(
        ext -> isfile(string(potfile_root, ext)),
        (".1", ".3", ".hst", ".12d", ".12s"),
    )
end

"""
    completedHistoryRanges(last_saved_index, numTS)

Return the state-history and per-step-history index ranges that have been
filled by an unsteady run.

State histories include the initial condition and every completed time step.
Per-step histories, such as strain arrays, have no initial-condition row and
therefore contain one fewer filled entry.
"""
function completedHistoryRanges(last_saved_index::Integer, numTS::Integer)
    numTS >= 1 || throw(ArgumentError("numTS must be positive"))
    1 <= last_saved_index <= numTS ||
        throw(ArgumentError("last_saved_index must be between 1 and numTS"))

    state_range = 1:Int(last_saved_index)
    step_range = 1:(Int(last_saved_index)-1)
    return (state = state_range, step = step_range)
end

"""
    endOpenFASTModules(inputs; openfast=OWENSOpenFASTWrappers, platform_initialized=inputs.platformActive, hd_initialized=platform_initialized, md_initialized=platform_initialized, ad_initialized=inputs.AD15On)

End initialized OpenFAST native modules without masking an upstream simulation
error.

Cleanup is attempted independently for HydroDyn, MoorDyn, and AeroDyn. Any
cleanup exceptions are returned as named tuples and, by default, also emitted as
warnings so a `finally` block can preserve the original failure.
"""
function endOpenFASTModules(
    inputs;
    openfast = OWENSOpenFASTWrappers,
    platform_initialized = inputs.platformActive,
    hd_initialized = platform_initialized,
    md_initialized = platform_initialized,
    ad_initialized = inputs.AD15On,
    warn_on_error = true,
)
    cleanup_errors = NamedTuple{(:label, :exception),Tuple{String,Any}}[]

    function cleanup(label, close_function)
        try
            close_function()
        catch err
            push!(cleanup_errors, (label = label, exception = err))
            if warn_on_error
                @warn "OpenFAST module cleanup failed" cleanup_module = label exception =
                    (err, catch_backtrace())
            end
        end
    end

    if hd_initialized
        cleanup("HydroDyn", openfast.HD_End)
    end
    if md_initialized
        cleanup("MoorDyn", openfast.MD_End)
    end
    if ad_initialized
        cleanup("AeroDyn", openfast.endTurb)
    end
    return cleanup_errors
end

function safeakima(x, y, xpt; extrapolate = false)
    if minimum(xpt)<(minimum(x)-(abs(minimum(x))*0.1+1e-4)) ||
       maximum(xpt)>(maximum(x)+abs(maximum(x))*0.1)
        msg="Extrapolating on akima spline results in undefined solutions minimum(xpt)<minimum(x) $(minimum(xpt))<$(minimum(x)) or maximum(xpt)>maximum(x) $(maximum(xpt))>$(maximum(x))"
        if !extrapolate
            throw(OverflowError(msg))
        else
            @warn msg
        end
    end
    return FLOWMath.akima(x, y, xpt)
end

"""
    createInitCondArray(initDisps, numNodes, numDOFPerNode)
    createInitCondArray(initDisps, nodeNumbers[, numDOFPerNode])

Creates the formatted initial conditions array needed by OWENSFEA

#Input
* `initDisps`: an array of length numDOFPerNode specifying the initial displacement of each DOF
* `numNodes`: the number of nodes in the given mesh. This preserves the legacy
  `1:numNodes` node numbering behavior.
* `nodeNumbers`: explicit mesh node numbers. Use this form when the mesh node IDs
  are not contiguous or do not start at 1.
* `numDOFPerNode`: the number of unconstrained degrees of freedom calculated in each node

#Output
* `initCond`: array containing initial conditions.
    initCond(i,1) node number for init cond i.
    initCond(i,2) local DOF number for init cond i.
    initCond(i,3) value for init cond i.
"""
function createInitCondArray(
    initDisps,
    numNodes::Integer,
    numDOFPerNode::Integer = length(initDisps),
)
    numNodes >= 0 || throw(ArgumentError("numNodes must be non-negative"))
    return createInitCondArray(initDisps, collect(1:numNodes), numDOFPerNode)
end

function createInitCondArray(
    initDisps,
    nodeNumbers::AbstractVector,
    numDOFPerNode::Integer = length(initDisps),
)
    numDOFPerNode > 0 || throw(ArgumentError("numDOFPerNode must be positive"))
    length(initDisps) <= numDOFPerNode ||
        throw(ArgumentError("initDisps length cannot exceed numDOFPerNode"))

    for nodeNumber in nodeNumbers
        nodeNumber isa Integer || throw(ArgumentError("nodeNumbers must contain integers"))
    end

    if all(iszero, initDisps)
        return []
    end

    activeDOFs = findall(!iszero, initDisps)
    initCond = zeros(length(activeDOFs)*length(nodeNumbers), 3)
    irow = 1
    for idof in activeDOFs
        for nodeNumber in nodeNumbers
            initCond[irow, :] = [nodeNumber, idof, initDisps[idof]]
            irow += 1
        end
    end

    return initCond
end

function setBCs(fixedDOFs, fixedNodes, numNodes, numDOFPerNode) #node, dof, bc
    if (fixedDOFs == []) && (fixedNodes == [])
        pBC = []
    else
        pBC = zeros(Int, numNodes*numDOFPerNode, 3)
        for i = 1:length(fixedNodes)
            pBC[((i-1)*numDOFPerNode+1):(i*numDOFPerNode), :] = hcat(
                ones(numDOFPerNode)*fixedNodes[i],
                collect(1:numDOFPerNode),
                zeros(Int, numDOFPerNode),
            )
        end
        for i = 1:length(fixedDOFs)
            newNodes = setdiff(1:numNodes, fixedNodes) # this avoids duplicating nodes already counted for by fixedNodes
            numNewNodes = length(newNodes)
            dofBCs =
                hcat(newNodes, ones(Int, numNewNodes)*fixedDOFs[i], zeros(Int, numNewNodes))
            pBC[
                (length(
                    fixedNodes,
                )*numDOFPerNode+(i-1)*numNewNodes+1):(length(
                    fixedNodes,
                )*numDOFPerNode+i*numNewNodes),
                :,
            ] = dofBCs
        end
        pBC = pBC[vec(mapslices(col -> any(col .!= 0), pBC, dims = 2)), :] #removes extra rows (i.e. rows of all zeros)
    end

    return pBC

end
