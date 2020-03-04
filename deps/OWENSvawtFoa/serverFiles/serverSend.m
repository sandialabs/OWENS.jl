function serverSend(message,output_stream)
%serverSend sends a double using a server socket connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   serverSend(message,d_output_stream)
%                    
%   This function sends a double using a server socket connection.
%
%   input:
%   message          = double value to send
%   output_stream    = socket output object directing sent messages

%   output:        (NONE)

    import java.net.ServerSocket
    import java.io.*

    output_stream.writeDouble(message);
    output_stream.flush;
    
end