function [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
#modalExec  Executive function for automated flutter/modal analysis
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
#
#   This function executes modal analysis.
#
#   input:
#   model          = object containing model information
#   mesh           = object containing mesh information
#   el             = object containing element information
#   displ          = displacement vector for use in pre-stressed analysis
#   Omega          = rotor speed (Hz)
#   OmegaStart     = rotor speed (Hz) from previous analysis if stepping
#                    through various rotor speeds, may be useful in load
#                    stepping
#
#   output:
#   freq           = modal frequencies of system
#   damp           = modal damping ratios of system

[elStorage] = initialElementCalculations(model,el,mesh); #performs initial element calculations

# [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage); #calculate mass properties of structure
convergedFreq = zeros(length(omegaArray),model.numModesToExtract);
convergedDamp = zeros(length(omegaArray),model.numModesToExtract);
for j=1:length(omegaArray) #loops over rotor speeds of interest
    Omega = omegaArray(j);
    
    #Do nonlinear iteration if needed
    if(model.spinUpOn)
        model.aeroElasticOn = false;
        model.aeroForceOn = true;
        [displ,~,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,Omega,OmegaStart,elStorage); #performs static analysis about specified operating condition
        #         save displnl displ
    else
        staticAnalysisSuccessful = true;
    end
    
    if(staticAnalysisSuccessful)
        tol = 1.0e-3;
        #get nominal frequencies without aeroelastic effects
        model.aeroElasticOn = true;
        modelCopy = model;
        modelCopy.aeroElasticOn = false;
        modelCopy.guessFreq = 0;
        
        if(j==1) #for the first rotor speed perform a modal analysis with initial guess frequency of zero to obtain frequency estimate
            [freqOrig,~,~,~,~] = linearAnalysisModal(modelCopy,mesh,el,displ,Omega,elStorage);
        else
            freqOrig = convergedFreq(j-1,:); #if not first rotor speed, use converged frequency from last rotor speed
        end
        
        for i=2:2:model.numModesToExtract
            
            rootfun = @(guessfreq) f0_linAnal(guessfreq,i,model,mesh,el,displ,Omega,elStorage);
            
            model.guessFreq = freqOrig(i);
            
            if Omega~=0.0 #If zero, then the output freq is always the input freq and no need to iterate
                options = optimset('Display','iter','TolX',1e-3);
                #                 [finalfreq,fval,exitflag,output] = fzero(rootfun,freqOrig(i),options);
                [finalfreq,fval,exitflag,output] = fminbnd(rootfun,0.0,freqOrig(i)*2,options);
                model.guessFreq = finalfreq;
                disp(exitflag);
            end
            
            [freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); #performs modal analysis during iteration
            convergedFreq(j,i) = freq(i);    #if converged, store and move on to next frequency
            convergedDamp(j,i) = damp(i);
            
            #             model.guessFreq = freqOrig(i); #do p-k iteration for flutter analysis of modes of interest
            #             converged = false;
            #             while(~converged) #TODO: this is terribly inefficient to calculate all of the modes, but only converge on one at a time, need to use some sort of ND root finder.
            #                 [freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); #performs modal analysis during iteration
            #
            fprintf('#f\n',j);
            fprintf('#f\n',i);
            #                 if(abs(freq(i) - model.guessFreq)<tol)  #check if modal frequency is converged
            #                     converged = true;
            #                     #             i
            #
            #                     convergedFreq(j,i) = freq(i);    #if converged, store and move on to next frequency
            #                     convergedDamp(j,i) = damp(i);
            #                 else
            #                     model.guessFreq = 0.5*(freq(i) + model.guessFreq); #if not converged select another guess frequency
            #                 end
            #             end
            
        end
    else
        error('Static analysis unsuccessful. Exiting');
    end
end

freq = convergedFreq;
damp = convergedDamp;

# save flutterRun freq damp omegaArray #save frequency, damping, and rotor speed array to .mat file
fidout=fopen([model.outFilename(1:end-4) '_FLUTTER.out'],'w');

for i = 1:length(freq)
    fprintf(fidout,'#e,#e\n',freq(i),damp(i));
end
fclose(fidout);

end

function output = f0_linAnal(guessfreq,i,model,mesh,el,displ,Omega,elStorage)
model.guessFreq = guessfreq;
[freq,~,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage);
output = abs(freq(i)-guessfreq);
end


