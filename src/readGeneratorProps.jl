function [genprops] = readGeneratorProps(generatorfilename)
#readGeneratorProps reads generator properties from file
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [genprops] = readGeneratorProps(generatorfilename)
#
#   This function reads generator properties from file.
#
#   input:
#   generatorfilenanme  = string containing generator property file name
#
#   output:
#   genprops          = model object containing generator properties

fid = fopen(generatorfilename); #open generator property file
if(fid~=-1) #if file can be opened
    #         genprops.ratedTorque = fscanf(fid,'#f',1); #store rated torque
    #         dum = fgetl(fid);
    #         genprops.zeroTorqueGenSpeed = fscanf(fid,'#f',1); #store zero torque generator zpeed
    #         dum = fgetl(fid);
    #         genprops.pulloutRatio = fscanf(fid,'#f',1); #store pullout ratio
    #         dum = fgetl(fid);
    #         genprops.ratedGenSlipPerc= fscanf(fid,'#f',1); #store rated generator slip percentage
    #         dum = fgetl(fid);
    #
    #         fclose(fid); #close generator propery file
    genprops = 0.0;
    error('GENERATOR NOT FULLY ENABLED')
else
    genprops = 0.0; #if generator property file does not exist, set object to null
end
end