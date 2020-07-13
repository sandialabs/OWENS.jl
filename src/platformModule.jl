function [rbData,Ywec_j,Accel_j] = platformModule(tSpan,q,CP2H,F1,F2,d_input_stream,d_output_stream)
#platformModule provides interface to an external platform module
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#                    
#   [rbData,Ywec_j,Accel_j] = platformModule(tSpan,q,CP2H,F1,F2,d_input_stream,d_output_stream)
#
#   #This function provides an interface to an external platform module. It
#   provides forcing on the platform and receives rigid body motion in
#   return.
#
#   input:
#   tspan           = time span for platform module to integrate over.
#   q               = stored intial condition vector for platform module
#   CP2H            = transformation matrix between platform and hub frame
#   F1              = External forcing vector on platform at beginning of
#                     time step
#   F2              = Estimated extenral forcing vector on platform at end
#                     of time step
#   d_input_stream  = socket input stream for receiving data from platform
#                     module
#   d_output_stream = socket input stream for sending data to platform
#                     module
#
#   output:
#   rbData          = vector containing rigid body data terms
#   Ywec_j          = current estimate of platform states
#   Accel_j         = current estimate of acceleration states

#======= hydrodynamics module =======================
#get rigid body motions a_i, theta_i, omega_i, omegaDot_i

#[Fext] is passed to WavEC each time step
Fext = zeros(1,12);
forceTransformation = [CP2H', zeros(3);zeros(3),CP2H'];
Fext(1:6)  = (forceTransformation*F1)';
Fext(7:12) = (forceTransformation*F2)';

# Initial Condition - this builds as the solution progresses.  The very
# first time, q(1,:) is defined above in the initial setup. q grabs the
# most recent solution from Y(:,:) and passes it to WavEC as the IC.

# Send data
serverSendVector(tSpan,d_output_stream) ; # send time info [1x2]
    
serverSendVector(Fext,d_output_stream) ; # send padded external force [1x12]
    
serverSendVector(q,d_output_stream) ; # Send padded initial condition [1 x s_pad_sol]

# After these three sends complete, WavEC solves and sends back data
    
# Recieve data
Accel_j = clientRecvVector(d_input_stream)'; #rigid body acceleration
    
Ywec_j = clientRecvVector(d_input_stream)'; #rigid body states

rbData = [CP2H*Accel_j(1:3)', CP2H*Ywec_j(end-2:end)', CP2H*Accel_j(4:6)']; #transformed rigid body data

end
