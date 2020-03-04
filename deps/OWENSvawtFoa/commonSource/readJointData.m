function [joint] = readJointData(inputfile)
%readJointData reads joint data file
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [joint] = readJointData(inputfile)
%                    
%   This function reads the joint data file and stores data in the joint
%   object.
%
%      input:
%      inputfile     = string containing joint data filename

%      output:
%      joint         = array containing joint data

fid = fopen(inputfile); %open joint file

count =1;
    while(~feof(fid)) 
        nodalInfo = fscanf(fid,'%i',4); %read in nodal info associated with joint [joint #, master node #, slave node #, joint type]
        if(isempty(nodalInfo))
            disp('No joint data detected'); %if no joint data detetcted return null joint object
            joint = [];
            break;
        end
        
        massStiffOrt = fscanf(fid,'%f',4); %reads in mass, stiffness, orientation of element attached to master joint
        temp = cat(2,nodalInfo',massStiffOrt');
        joint(count,:) = temp;        %store data in joint array
        count = count + 1;
        
    end
fclose(fid);    %close file

end

