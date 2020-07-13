function [nlParams] = readNLParamsFile(inputfile)
#readNLParamsFile   reads file for nonlinear iteration/load stepping
# **********************************************************************
# *                   Part of the SNL OWENS toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [nlParams] = readNLParamsFile(inputfile)
#
#   This function reads a nonlinear parameter for nonlinear iteration and
#   load stepping.
#
#      input:
#      inputfile    = string containing filename of main input file
#
#      output:
#      nlParams     = struct containing nonlinear parameter information.

filePrefix = inputfile(1:end-6); #get main file prefix
nlFile = [filePrefix,'.nl'];     #create .nl file name string

fid = fopen(nlFile,'r'); #attempt to open file

if(fid ~= -1) #if file can be opened read it
    #         iterationType = fscanf(fid,'#c',2); fgetl(fid); #read iteration type 'NR' = Newton Raphson, 'DI' = Direct Iteration
    #         tolerance = fscanf(fid,'#f',1); fgetl(fid);
    #         maxIterations    = fscanf(fid,'#i',1); fgetl(fid); #read in maximum iterations allowed per time step
    #         temp  = fscanf(fid,'#i',1);
    #         if(temp == 0) # if temp = 0 adaptive load stepping, read in params
    #             fgetl(fid);
    #             adaptiveLoadSteppingFlag = true;
    #             maxNumLoadSteps     = fscanf(fid,'#i',1); fgetl(fid);
    #             minLoadStep      = fscanf(fid,'#f',1); fgetl(fid);
    #             minLoadStepDelta= fscanf(fid,'#f',1); fgetl(fid);
    #             fclose(fid); #close file
    #         elseif(temp>0) #if temp > 0, prescribed load stepping profile, read in params
    #             adaptiveLoadSteppingFlag = false;
    #             prescribedLoadStep = zeros(temp,1);
    #             for i=1:temp
    #                 prescribedLoadStep(i) = fscanf(fid,'#f',1);
    #             end
    #             fgetl(fid);
    #             maxNumLoadSteps = 1e6;
    #             fclose(fid); #close file
    #         else
    #             fclose(fid);
    #             error('Load stepping parameter in .nl file not recognized. Exiting.');
    #         end
    error('NLParams not fully implemented')

else #default parameters if file cannot be opened (doesnt exist)
    adaptiveLoadSteppingFlag = true; #uses adaptive load stepping by default
    iterationType    = 'NR';
    tolerance        = 1.0e-6;
    maxIterations    = 50;
    maxNumLoadSteps  = 20;
    minLoadStep      = 0.05;
    minLoadStepDelta = 0.05;
end

#assign parameters to nlParams file
nlParams.iterationType = iterationType;
nlParams.adaptiveLoadSteppingFlag = adaptiveLoadSteppingFlag;
nlParams.tolerance = tolerance;
if(adaptiveLoadSteppingFlag)
    nlParams.maxIterations = maxIterations ;
    nlParams.maxNumLoadSteps = maxNumLoadSteps;
    nlParams.minLoadStepDelta = minLoadStepDelta;
    nlParams.minLoadStep = minLoadStep;
    nlParams.prescribedLoadStep = 1; # not used but must be declared
else
    error('NLParams not fully implemented')
    #     nlParams.prescribedLoadStep = prescribedLoadStep;
    #     nlParams.maxIterations = maxIterations ;
    #     nlParams.maxNumLoadSteps = maxNumLoadSteps;
end


end
