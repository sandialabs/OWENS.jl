function vec = clientRecvVector(input_stream)
%clientRecvVector receives a vector of doubles using a server socket connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   vec = clientRecvVector(d_input_stream)
%                    
%   This function receives a vector of doubles using a server socket connection.
%
%   input:
%   input_stream   = socket input object receiving messages

%   output:        
%   vec              = vector of doubles received through socket

    import java.net.Socket
    import java.io.*

    len = clientRecv(input_stream);
    vec = zeros(len,1);
    
    for i=1:len
        vec(i) = clientRecv(input_stream);
    end
end