function [loadStep,loadStepPrev,displ,displPrev,staticAnalysisSuccessful,staticAnalysisComplete] = updateLoadStep(iterationCount,loadStepParams,loadStep,loadStepPrev,loadStepCount,displPrev,displ)
#updateLoadStep    updates load step for static nonlinear analysis
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [loadStep,loadStepPrev,displ,displCopy,staticAnalysisSuccessful,...
#    staticAnalysisComplete] = updateLoadStep(iterationCount,...
#    loadStepParams,loadStep,loadStepPrev,loadStepCount,displCopy,displ)
#
#   This function updates the load stepping parameter whether through means
#   of adaptive loadstepping or a specified load step profile.
#
#      input:
#      iterationCount      = number of iterations for current load step
#      loadStepParams      = struct containing load step parameters
#      loadStep            = load step value for current load step
#      loadStepPrev        = load step value for previous load st ep
#      loadStepCount       = number of load steps performed up to this
#                            point
#      displPrev           = converged displacement vector form previous
#                            load step
#      displ               = displacement vector at current load step
#
#      output:
#      loadStep            = new load step value
#      loadStepPrev        = load step value for previous load step
#      displ               = most up to date displacement vector in load
#                            stepping procedure
#      displPrev           = displacement vector at end of previous load
#                            step
#      staticAnalysisSuccessful = boolean flag, true if load
#                                 step was completed successfully
#      staticAnalysisComplete   = boolean flag, true if analysis is complete
#
#
staticAnalysisComplete = false; #initialize variable

if(loadStepParams.adaptiveLoadSteppingFlag) #for adaptive load stepping option
    #check if maximum number of load steps has been exceeded.
    if(loadStepCount > loadStepParams.maxNumLoadSteps)
        error('Maximum number of load steps exceeded. Exiting.');
    end
    #calculate new loadstep adaptively
    [loadStep,loadStepPrev,staticAnalysisSuccessful,staticAnalysisComplete] = adaptiveLoadStepping(iterationCount,loadStepParams,loadStep,loadStepPrev);
    
else #for prescribed load stepping option
    if(iterationCount<loadStepParams.maxIterations) #see if previous load step was successful
        staticAnalysisSuccessful= true;
        if(loadStep == 1.0) #if load step = 1.0, analysis is complete
            staticAnalysisComplete = true;
        else
            staticAnalysisComplete = false;
            loadStepPrev = loadStep;
            loadStep = loadStepParams.prescribedLoadStep(loadStepCount);
            fprintf('Prescribed load step: #f\n',loadStep);
        end
    else
        staticAnalysisSuccessful = false;
    end
    
end
if(staticAnalysisSuccessful)
    displPrev = displ; #update displacementPrev variable if previous load step was successful.
else
    displ = displPrev; #reset to displ vector to that at start of load step if load step was unsuccessful
    
    if(~loadStepParams.adaptiveLoadSteppingFlag) #if
        error('Maximum number of iterations exceeded in prescribed loadstep profile. Exiting.');
    end
end

end

function [loadStep,loadStepPrev,staticAnalysisSuccessful,staticAnalysisComplete] = adaptiveLoadStepping(iterationCount,loadStepParams,loadStep,loadStepPrev)
#This function performs updates a loadstep adaptively.

staticAnalysisComplete = false;
loadStepOrig = loadStep;

if(iterationCount>=loadStepParams.maxIterations)                   #check for exceeding max iterations
    msgId = 1; #corresponds to a message  saying load step was unsuccessful and is being reduced
    staticAnalysisSuccessful = false;
else
    if(loadStep == 1.0)
        msgId = 2; #corresponds to a message saying loadstep was successful and load stepping is finished
    else
        msgId = 3; #corresponds to a message saying loads tep was successful and analysis is proceeding to next load step
    end
    staticAnalysisSuccessful=true;
    if(loadStep == 1.0)
        staticAnalysisComplete = true;
    end
    loadStepPrev = loadStep; #update previous load step and end loadstep
    loadStep = 1.0;
end

loadStepCopy = loadStep; #make copy of load step for later checks
loadStep = 0.5*(loadStep + loadStepPrev); #update load step


if(abs(loadStep-loadStepPrev)< loadStepParams.minLoadStepDelta) #enforces delta load step is not below the minimum specified value
    if(loadStep<loadStepPrev)
        loadStep = loadStep - loadStepParams.minLoadStepDelta;
    else
        loadStep = loadStep + loadStepParams.minLoadStepDelta;
    end
end

if(loadStep < loadStepParams.minLoadStep ) #check that load step is not below the minimum specified load step
    loadStep = loadStepParams.minLoadStep;
    if(loadStepCopy == loadStep)
        error('Minimum load step reached. Exiting.');
    end
elseif(loadStep > 1.0) #if load step has extended beyond 1.0, set to 1.0
    loadStep = 1.0;
end

#print loadstepping message to command line
if(msgId == 1)
    fprintf('Max iterations exceeded for loadstep = #f.\t\t Reducing load step to #f.\n',loadStepOrig,loadStep);
elseif(msgId == 2)
    fprintf('Nonlinear iteration successful for loadstep = #f.\t Nonlinear static analysis complete.\n',loadStep);
elseif(msgId == 3)
    fprintf('Nonlinear iteration successful for loadstep = #f.\t Increasing loadstep size to #f.\n',loadStepOrig,loadStep);
end


end