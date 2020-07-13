function [FAero,FAeroDof] = aeroModule(model,time,u,Omega,azi,numDOFPerNode,d_input_streamAero,d_output_streamAero)
#aeroModule provides interface to an external aerodynamics module
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#                    
#   [FAero,FAeroDof] = aeroModule(model,time,u,Omega,azi,numDOFPerNode,d_input_streamAero,d_output_streamAero)
#
#   #This function provides an interface to an external aerodynamics module. It
#   provides blade displacements, rotor speed, and rotor position and
#   receives aerodynamic forces.
#
#   input:
#   model               = object  containing model data
#   time                = time at which forces are requested
#   u                   = diplacement vector for overall configuration
#   Omega               = rotor speed (Hz)
#   azi                 = rotor azimuth position (radians)
#   numDOFPerNode       = number of degrees of freedom per node
#   d_input_streamAero  = socket input stream for receiving data from platform
#                         module
#   d_output_streamAero = socket input stream for sending data to platform
#                         module
#
#   output:
#   rbData          = vector containing rigid body data terms
#   Ywec_j          = current estimate of platform states
#   Accel_j         = current estimate of acceleration states

    useDisplacements = false; #if set to false displacements sent to aeromodule are zero, if true actual displacements of blades are sent to aero module

    bladeData = model.bladeData;
    bladeNodes = bladeData.nodeNum;
    numBlades = bladeData.numBlades;    #get number of blades
    numEntriesH = length(bladeData.h);    #get number of entries in h for blade
    numEntriesHPerBlade = numEntriesH/numBlades; #get number of points per blade 
    
    dBlade = zeros(numEntriesHPerBlade,4,numBlades); #eventually populate with structural motions
    
    if(useDisplacements)
        for iblade=1:numBlades
            for jblade=1:numEntriesHPerBlade
                    nodeIndex = (iblade - 1)*numEntriesHPerBlade + jblade;
                    bladeNode = bladeNodes(nodeIndex);
                    dBlade(jblade,1,iblade) = u((bladeNode-1)*numDOFPerNode + 1); #uses previous timestep for now since aero code can't iterate
                    dBlade(jblade,2,iblade) = u((bladeNode-1)*numDOFPerNode + 2);
                    dBlade(jblade,3,iblade) = u((bladeNode-1)*numDOFPerNode + 3);
            end
        end
    end
    
    
   #dBlade = dBlade.*0; #zero out displacement matrix for now (one way coupling)
 

    #azi_j is current estimate (or prescribed) rotor speed at "time"
    #Omega_j is current estimate (or prescribed) rotor speed at "time"
    
    ##generalize this function call
    dataToAero = decodeVec(2,time,Omega*2*pi,azi,numEntriesHPerBlade,numBlades,dBlade(:,:,1),dBlade(:,:,2)); # 2times dBlade displacement per blade
    serverSendVector(dataToAero,d_output_streamAero)
    
    dataFromAero = clientRecvVector(d_input_streamAero); # receive force back from forcing module
    [~,~,aeroOutput] = encodeVec(dataFromAero);
	
	index = 1;
	aeroF = aeroOutput.forces;
    
    FAero = zeros(numBlades*numEntriesHPerBlade*numDOFPerNode,1);
    FAeroDof = FAero;
    
	for iblade=1:numBlades
		for jblade=1:numEntriesHPerBlade
			for kblade=1:numDOFPerNode
				FAero(index) = aeroF(jblade,kblade,iblade);
                if(jblade == 1 || jblade == numEntriesHPerBlade)
                    FAero(index) = 0.0; 
                end
				nodeIndex = (iblade - 1)*numEntriesHPerBlade + jblade;
				FAeroDof(index) = (bladeNodes(nodeIndex)-1)*numDOFPerNode + kblade;
				index = index + 1;
            end
		end
	end
end

