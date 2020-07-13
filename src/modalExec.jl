function [freq,damp] = modalExec(model,mesh,el,displ,Omega,OmegaStart)
#modalExec  Executive function for modal analysis
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [freq,damp]=modalExec(model,mesh,el,displ,Omega,OmegaStart,fid)
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

#Do nonlinear iteration if needed
if(model.spinUpOn)
model.aeroElasticOn = false;
model.aeroForceOn = true;
[displ,~,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,Omega,OmegaStart,elStorage); #performs static analysis about specified operating condition
#save displnl displ #saves displacement vector from static analysis TODO: This doesn't appear to be used
else
    staticAnalysisSuccessful = true;
end

if(staticAnalysisSuccessful)
[freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); #performs modal analysis
else
    error('Static analysis unsuccessful. Exiting');
end
end

