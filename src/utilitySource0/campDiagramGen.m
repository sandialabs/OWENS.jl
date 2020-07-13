function [freq] = campDiagramGen(inputFileName,outputFileName,rotorSpeedArray,spinUpOn,numModes)
#campDiagramGen performs modal analysis at various rotor speeds
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [freq] = campDiagramGen(inputFileName,outputFileName,...
#                           rotorSpeedArray,spinUpOn,numModes)
#                    
#   This function peforms modal analysis at  various rotor speeds for
#   creating a campbell diagram.
#
#      input:
#      inputFileName      = string containing .owens filename  
#      outputFileName     = string containing output filename
#      rotorSpeedArray    = array of rotor speeds (Hz) to perform modal
#                           analysis at
#      spinUpOn           = activates a nonlinear caclulation at static
#                           equilibrium condition and employs pre-stress 
#                           effects in modal analysis
#      numModes           = number of lower system modes to extract
# 
#      output:
#      freq               = n x m array of modal frequencies (Hz) for m
#                           modes at n rotor speeds
    
    freq = zeros(length(rotorSpeedArray),numModes);

    for i=1:length(rotorSpeedArray)
        rotorSpeed = rotorSpeedArray(i);

        [freq(i,:),~] = owens(inputFileName,'M',rotorSpeed,spinUpOn,numModes);
    end

    save(outputFileName,'freq','rotorSpeedArray');
end