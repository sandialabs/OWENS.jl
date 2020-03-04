function terminateServer(server,socket,verboseFlag)
%terminateServer terminates a server connection for a socket interface
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   terminateServer(server_socket,output_socket,verboseFlag)
%                    
%   This function terminates a server connection.
%
%   input:
%   server          = server object
%   socket          = socket object
%   verboseFlag     = flag controlling printing output messages

%   output:        (NONE)

    server.close;
    socket.close;
    if(verboseFlag)
        disp('Sever has been terminated.');
    end
end

