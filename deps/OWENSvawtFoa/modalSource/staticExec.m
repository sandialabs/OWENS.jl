function staticExec(model,mesh,el,displ,Omega,OmegaStart)
%staticExec  Executive function for static analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   staticExec(model,mesh,el,displ,Omega,OmegaStart,fid)
%                    
%   This function executes static analysis.
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
%   output:        (NONE)

    [elStorage] = initialElementCalculations(model,el,mesh); %performs initial element calculations

    [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage); %calculate mass properties of structure

    %Do nonlinear iteration if needed
    model.aeroElasticOn = false;
    model.aeroForceOn = true;
    [displ,elStrain,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,Omega,OmegaStart,elStorage); %performs static analysis about specified operating condition
    
    if(staticAnalysisSuccessful)
        save(model.outFilename,'displ','elStrain'); %saves static analysis displacement output
    else
        error('Static analysis unsuccessful. Exiting');
    end
end

