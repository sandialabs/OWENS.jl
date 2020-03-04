function message = clientRecv(input_stream)
%clientRecv receives a double using a server socket connection
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   message = clientRecv(input_stream)
%                    
%   This function receives a double using a server socket connection.
%
%   input:
%   input_stream   = socket input object receiving messages

%   output:        
%   message          = double value received through socket

    import java.net.Socket
    import java.io.*

    message = input_stream.readDouble();
end