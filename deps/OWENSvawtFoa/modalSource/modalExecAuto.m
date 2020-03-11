function [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
%modalExec  Executive function for automated flutter/modal analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
%
%   This function executes modal analysis.
%
%   input:
%   model          = object containing model information
%   mesh           = object containing mesh information
%   el             = object containing element information
%   displ          = displacement vector for use in pre-stressed analysis
%   Omega          = rotor speed (Hz)
%   OmegaStart     = rotor speed (Hz) from previous analysis if stepping
%                    through various rotor speeds, may be useful in load
%                    stepping
%
%   output:
%   freq           = modal frequencies of system
%   damp           = modal damping ratios of system

[elStorage] = initialElementCalculations(model,el,mesh); %performs initial element calculations

% [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage); %calculate mass properties of structure
convergedFreq = zeros(length(omegaArray),model.numModesToExtract);
convergedDamp = zeros(length(omegaArray),model.numModesToExtract);
for j=1:length(omegaArray) %loops over rotor speeds of interest
    Omega = omegaArray(j);
    
    %Do nonlinear iteration if needed
    if(model.spinUpOn)
        model.aeroElasticOn = false;
        model.aeroForceOn = true;
        [displ,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,Omega,OmegaStart,elStorage); %performs static analysis about specified operating condition
        save displnl displ
    else
        staticAnalysisSuccessful = true;
    end
    
    if(staticAnalysisSuccessful)
        tol = 1.0e-3;
        %get nominal frequencies without aeroelastic effects
        model.aeroElasticOn = true;
        modelCopy = model;
        modelCopy.aeroElasticOn = false;
        modelCopy.guessFreq = 0;
        
        if(j==1) %for the first rotor speed perform a modal analysis with initial guess frequency of zero to obtain frequency estimate
            [freqOrig,~,~,~,~] = linearAnalysisModal(modelCopy,mesh,el,displ,Omega,elStorage);
        else
            freqOrig = convergedFreq(j-1,:); %if not first rotor speed, use converged frequency from last rotor speed
        end
        
        
        for i=2:2:model.numModesToExtract %do p-k iteration for flutter analysis of modes of interest
            model.guessFreq = freqOrig(i);
            converged = false;
            while(~converged)
                [freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); %performs modal analysis during iteration
                
                
                if(abs(freq(i) - model.guessFreq)<tol)  %check if modal frequency is converged
                    converged = true;
                    %             i
                    convergedFreq(j,i) = freq(i);    %if converged, store and move on to next frequency
                    convergedDamp(j,i) = damp(i);
                else
                    model.guessFreq = 0.5*(freq(i) + model.guessFreq); %if not converged select another guess frequency
                end
            end
            
        end
    else
        error('Static analysis unsuccessful. Exiting');
    end
end

freq = convergedFreq;
damp = convergedDamp;

save flutterRun freq damp omegaArray %save frequency, damping, and rotor speed array to .mat file
end

