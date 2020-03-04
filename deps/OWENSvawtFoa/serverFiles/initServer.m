function [data_output_stream,server_socket,output_socket] = initServer(output_port,number_of_retries)
%initServer initializes a server connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [data_output_stream,server_socket,output_socket] = initServer(output_port,...
%                                                   number_of_retries)
%                    
%   This function initializes a server connection using sockets.
%
%   input:
%   output_port        = port number for client connection
%   number_of_retries  = number of attempts to connect to output port
%
%   output:        
%   data_output_stream    = double value received through socket
%   server_socket         = server socket object
%   output_socket         = output socket object

    import java.net.ServerSocket
    import java.io.*

    retry             = 0;

    server_socket  = [];
    output_socket  = [];
    data_output_stream = [];

    while true

        retry = retry + 1;

        try
            if ((number_of_retries > 0) && (retry > number_of_retries))
                fprintf(1, 'Too many connection attempts. Exiting.');
                break;
            end

            % wait for 1 second for client to connect server socket
            server_socket = ServerSocket(output_port);
            server_socket.setSoTimeout(1001);

            output_socket = server_socket.accept;

            output_stream   = output_socket.getOutputStream;
            data_output_stream = DataOutputStream(output_stream);
            disp('Server has connected to client.');
            break;
            
        catch
            if ~isempty(server_socket)
                server_socket.close
            end

            if ~isempty(output_socket)
                output_socket.close
            end
           
        end
    end

end

