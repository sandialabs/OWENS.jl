function [initCond] = readInitCond(filename)
%readInitCond reads initial conditions
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [initCond] = readInitCond(filename)
%                    
%   This function reads initial conditions from file
%
%   input:
%   filename      = string containing file name for initial conditions file
%
%   output:
%   initCond      = array containing initial conditions

    initCond =[]; %initialize intial condition to null
    fid = fopen(filename); %open initial  conditions file
    
    if(fid>0) %if file exists begin read in
        index = 1;
        while(~feof(fid))
            temp1 = fscanf(fid,'%i',2); %read node number and local DOF number for initial cond.
            temp2 = fscanf(fid,'%f',1); %read value for initial cond.
            
            %place node number, dof number and value into array
            initCond(index,1:3) = [temp1(1), temp1(2), temp2(1)]; 
            
            index = index + 1;
        end
    end
end