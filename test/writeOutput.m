function [freqSorted,dampSorted,imagCompSignSorted] = writeOutput(freq,damp,phase1,phase2,imagComponentSign,fid)
%writeOutput writes output from modal analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freqSorted,dampSorted,imagCompSignSorted] = writeOutput(freq,damp,...
%                                                phase1,phase2,...
%                                                imagComponentSign,fid)
%
%   This function writes an output file for modal analysis.
%
%      input:
%      freq               = array of modal frequencies
%      damp               = array of modal damping ratios
%      phase1             = array of in phase mode shapes
%      phase2             = array of out of phase mode shapes
%      imagComponentSign  = array of sign of imaginary components
%      fid                = file identifier for output
%
%      output:
%      freqSorted         = array of sorted(by frequency) modal frequencies
%      dampSorted         = array of sorted(by frequency) modal damping ratios
%      imagCompSignSorted = array of sorted(by frequency) of imaginarycomponentSign array

[l,~] = size(phase1(:,:,1)); %gets number of nodes for mode shape printing

[freq,map,posIndex] = bubbleSort(freq); %sorts frequency

dampSorted = zeros(1,length(freq)-1);
freqSorted = zeros(1,length(freq)-1);
imagCompSignSorted = zeros(1,length(freq)-1);

index = 1;
for i=posIndex:1:posIndex+(length(freq)-1) %prints mode frequency, damping and in/out of phase mode shapes
    fprintf(fid,'MODE # %0.0f \n\n',index);
    fprintf(fid,'Frequency: %e: \n',freq(i));
    fprintf(fid,'Damping %e: \n',damp(map(i)));
    fprintf(fid,'0 deg Mode Shape:\n');
    fprintf(fid,'U_x          U_y          U_z          theta_x     theta_y     theta_z \n');
    dampSorted(i) = damp(map(i));
    freqSorted(i) = freq(i);
    imagCompSignSorted(i) = imagComponentSign(map(i));
    
    for j=1:l
        fprintf(fid,'%8.6f \t',phase1(j,:,map(i)));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'90 deg Mode Shape:\n');
    fprintf(fid,'U_x          U_y          U_z          theta_x     theta_y     theta_z \n');
    for j=1:l
        fprintf(fid,'%8.6f \t',phase2(j,:,map(i)));
        fprintf(fid,'\n');
    end
    
    if(i<posIndex+(length(freq)-1))
        fprintf(fid,'\n\n');
    end
    
    index = index + 1;
    
end
end
