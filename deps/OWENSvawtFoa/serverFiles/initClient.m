function [data_input_stream,input_socket] = initClient(host,port,number_of_retries)
%initClient initializes a server connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [data_input_stream,input_socket] = initClient(host,port,number_of_retries)
%                    
%   This function initializes a client connection using sockets.
%
%   input:
%   host               = host address for server to connect client to
%   port               = port number for server on host to connect client
%                        to
%   number_of_retries  = number of attempts to connect to server
%
%   output:        
%   data_input_stream    = input stream object
%   input_socket         = input socket ojbect

    import java.net.Socket
    import java.io.*
    
    retry        = 0;
    input_socket = [];
    data_input_stream = [];

    while true

        retry = retry + 1;
        if ((number_of_retries > 0) && (retry > number_of_retries))
            fprintf(1, 'Too many connection attempts. Exiting.');
            break;
        end
        
        try

            input_socket = Socket(host, port);

            % get a buffered data input stream from the socket
            input_stream   = input_socket.getInputStream;
            data_input_stream = DataInputStream(input_stream);
            disp('Client has connected to server.');
            break;
            
        catch
            if ~isempty(input_socket) 
                input_socket.close;
            end


        end
    end    

end

