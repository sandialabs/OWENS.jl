function [output] = readResultsModalOut(resultsFile,numNodes)
#readResults reads modal analysis results file from the OWENS toolkit
# **********************************************************************
# *                   Part of SNL VAWTGen                              *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [phase1,phase2] = readResults(resultsFile,modeNum,numNodes)
#
#   This function reads modal analysis results file from the OWENS toolkit.
#
#      input:
#      resultsfile      = string containing results file name
#      modeNum          = integer denoting selected mode to plot
#      sf               = scale factor for mode shape displacements
#      meshSeg          = array denoting the number of elements per
#                         component (cleans up the plot)
#
#      output:
#      phase1           = array containing mode shape of in-phase mode
#                         shape component
#      phase2           = array containing mode shape of out-of-phase mode
#                         shape component

# determine the number of modes in the output file
fid = fopen(resultsFile,'r');
A = textscan(fid, '#s','delimiter','/n');
fclose(fid);
numModes = sum(strcmp(A{1},'0 deg Mode Shape:'));

# read in the modeshapes and frequencies
fid = fopen(resultsFile,'r'); #open results file for reading

for modecount = 1:numModes
    # read in the header for the mode shape
    headerInd = 0;
    while ~headerInd
        line = fgetl(fid);
        if contains(line, 'Frequency')
            freq = strrep(line,'Frequency','');
            freq = strrep(freq,':','');
            freq = strrep(freq,' ','');
            output{modecount}.frequency = real(str2double(freq));
        end
        if contains(line, 'Damping')
            damp = strrep(line,'Damping','');
            damp = strrep(damp,':','');
            damp = strrep(damp,' ','');
            output{modecount}.damping = real(str2double(damp));
            headerInd = 1;
            fgetl(fid); # read in the line declaring the mode shape orientation
        end
    end

    # read in the 0-deg shape
    tblHdr = textscan(fid, '#s #s #s #s #s #s',1);
    phase1 = textscan(fid, '#f #f #f #f #f #f',numNodes);   # real part of eigenvector
    output{modecount}.InPhaseShape = array2table(cell2mat(phase1),'VariableNames',[tblHdr{:}]);

    # read in the 90-deg shape
    headerInd = 0;
    while ~headerInd
        line = fgetl(fid);
        if contains(line, '90 deg Mode Shape')
            headerInd = 1;
        end
    end

    tblHdr = textscan(fid, '#s #s #s #s #s #s',1);
    phase2 = textscan(fid, '#f #f #f #f #f #f', numNodes);   # imaginary part of eigenvector
    output{modecount}.OutPhaseShape = array2table(cell2mat(phase2),'VariableNames',[tblHdr{:}]);
end

# close the results file
fclose(fid);

end
