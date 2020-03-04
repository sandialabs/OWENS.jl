function [nodalTerms] = readNodalTerms(filename)
%readNodalTerms reads concentrated nodal terms file
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [nodalTerms] = readNodalTerms(filename)
%                    
%   This function reads the nodal terms file and stores data in the nodal
%   terms object.
%
%      input:
%      filename      = string containing nodal terms filename
%
%      output:
%      nodalTerms    = object containing concentrated nodal data

    concStiff=[]; %initialize stiffness, load, and mass arrays to null
    concLoad=[];
    concMass=[];
    concStiffGen = [];
    concMassGen = [];
    concDampGen = [];
    fid=fopen(filename,'r');  %open nodal terms file
    index = 1;
    if(fid>0)  %if file is opened successfully
    while(~feof(fid))
        
         temp = strsplit(fgetl(fid));
        
        nodeNum(index,1) = str2num(temp{1});   %read in node number
        termType{index} = temp{2};  %read in concentrated term type (M=mass, K = stiffness, F = force)
        if(length(temp) == 4)
            localDof{index} = str2num(temp{3});
            val(index) = str2num(temp{4});
        else
            localDof{index} = [str2num(temp{3}) str2num(temp{4})];
            val(index) = str2num(temp{5});
        end
%         val(index,1) = fscanf(fid,'%f',1);       %read in value for concentrated nodal term
        index = index + 1;
    end
    
    fIndex = 1;    %counters for concentrated nodal force, stiffness, and mass
    kIndex = 1;
    mIndex = 1;
    mIndex2 = 1;
	kIndex2 = 1;
	cIndex2 = 1;
    for i=1:length(nodeNum)
       if(strcmpi(termType{i},'F'))    %read in concentrated load data
           concLoad(fIndex).nodeNum = nodeNum(i);
           concLoad(fIndex).dof = localDof{i};
           concLoad(fIndex).val = val(i);
           fIndex = fIndex + 1;
       end
       
       if(strcmpi(termType{i},'K'))    %read in concentrated stiffness data
           concStiff(kIndex).nodeNum = nodeNum(i);
           concStiff(kIndex).dof = localDof{i};
           concStiff(kIndex).val = val(i);
           kIndex = kIndex + 1; 
       end
       
       if(strcmpi(termType{i},'M'))   %read in concentrated mass data
          concMass(mIndex).nodeNum = nodeNum(i);
          concMass(mIndex).dof = localDof{i};
          concMass(mIndex).val = val(i);
          mIndex = mIndex + 1;
       end
	   
	   if(strcmpi(termType{i},'M6'))   %read in concentrated mass data
          concMassGen(mIndex2).nodeNum = nodeNum(i);
          temp = localDof{i};
          concMassGen(mIndex2).dof1 = temp(1);
		  concMassGen(mIndex2).dof2 = temp(2);
          concMassGen(mIndex2).val = val(i);
          mIndex2 = mIndex2 + 1;
       end
   
	   if(strcmpi(termType{i},'K6'))   %read in concentrated stiffness data
          concStiffGen(kIndex2).nodeNum = nodeNum(i);
           temp = localDof{i};
          concStiffGen(kIndex2).dof1 = temp(1);
		  concStiffGen(kIndex2).dof2 = temp(2);
          concStiffGen(kIndex2).val = val(i);
          kIndex2 = kIndex2 + 1;
       end
	   
	   if(strcmpi(termType{i},'C6'))   %read in concentrated damp data
          concDampGen(cIndex2).nodeNum = nodeNum(i);
          temp = localDof{i};
          concDampGen(cIndex2).dof1 = temp(1);
		  concDampGen(cIndex2).dof2 = temp(2);
          concDampGen(cIndex2).val = val(i);
          cIndex2 = cIndex2 + 1;
       end
    end
        fclose(fid);
    end

	%store concentrated nodal term data in nodalTerms object
    nodalTerms.concLoad = concLoad;
    nodalTerms.concStiff = concStiff;
    nodalTerms.concMass = concMass;
	nodalTerms.concStiffGen = concStiffGen;
	nodalTerms.concMassGen = concMassGen;
	nodalTerms.concDampGen = concDampGen;
end