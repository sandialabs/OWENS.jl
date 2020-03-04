function terminateClient(socket,verboseFlag)
%terminateClient terminates a client connection for a socket interface
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   terminateClient(input_socket,verboseFlag)
%                    
%   This function terminates a client connection.
%
%   input:
%   socket     = socket object
%   verboseFlag      = flag controlling printing output messages

%   output:        (NONE)

    socket.close;
    if(verboseFlag)
        disp('Client has been terminated.');
    end
end



