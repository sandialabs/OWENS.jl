function serverSendVector(vec,output_stream)
%serverSendVector sends a double using a server socket connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   serverSendVector(vec,d_output_stream)
%                    
%   This function sends a vector of doubles using a server socket
%   connection.
%
%   input:
%   vec              = vector of doubles to send  to send
%   output_stream    = socket output object directing sent messages

%   output:        (NONE)

    import java.net.ServerSocket
    import java.io.*
    
    len = length(vec);
    serverSend(len,output_stream);
    
    for i=1:len
        serverSend(vec(i),output_stream);
    end
    
end
