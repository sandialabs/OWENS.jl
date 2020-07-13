function [phase1,phase2] = readResults(resultsFile,modeNum,numNodes)
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

fid = fopen(resultsFile,'r'); #open results file for reading

modecount = 1;

while(modecount <= modeNum && ~feof(fid))
    dum = fscanf(fid,'#s',17);            #read in mode shape data up to modeNum
    for j=1:numNodes
       temp1 = fscanf(fid,'#e',6);
       phase1(j,:,modecount) = temp1';
    end

    dum = fscanf(fid,'#s',10);
    for j=1:numNodes
       temp2 = fscanf(fid,'#e',6);
       phase2(j,:,modecount) = temp2';
    end
    modecount = modecount + 1;
end
fclose(fid); #close results file

if(modeNum > modecount) #error if modeNum is greater than available modes
   error('modeNum is greater than number of modes in results file. Exiting');
end

phase1 = phase1(:,:,modeNum); #assign modeNum mode shapes to phase1 and phase2
phase2 = phase2(:,:,modeNum);

end

