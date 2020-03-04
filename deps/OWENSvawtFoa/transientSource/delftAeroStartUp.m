function [d_output_stream,d_input_stream,server_socket,input_socket,output_socket] = delftAeroStartUp(host,port,output_port,model)
%delftAeroStartUp   starts up TU Delft aerodynamics module
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [d_output_stream,d_input_stream,server_socket,...
%    input_socket,output_socket] = delftAeroStartUp(host,port,...
%                                                   output_port,model)
%                    
%   This function starts up the TU Delft aerodyanmics module from the 
%   OWENS toolkit
%
%   input:
%   host            = host address for client connectino
%   port            = port number for client connection
%   output_port     = port number for server connection
%   model           = model object
%
%   output:
%   d_output_stream = output stream object for server socket connections
%   d_input_stream  = input stream object for client socket connections
%   server_socket   = server object for server socket connections
%   input_socket    = input object for client socket connections

    number_of_retries = 100000;
%......................................
%%
    import java.net.ServerSocket
    import java.net.Socket
    import java.io.*
    
    %launch force generating matlab instance
	disp('Launching TU Delft UMPM Aerodynamics Code');
    if (ispc())
        aeroLaunchString='matlab -noFigureWindows -r "cd(''C:\work\projects\delftCode\version_110313''); START;" &';
        system(aeroLaunchString); 
    else
        error('Aero code not implemented for non pc platform');
        system('matlab -logfile UMPM.log -r START &');
%     system('at -f RunMatlabBatch.sh now'); % runs code in the background
    end
pause(6);

[d_output_stream,server_socket,output_socket] = initServer(output_port,number_of_retries); %initializes server to send data to aero module
[d_input_stream,input_socket] = initClient(host,port,number_of_retries); %initializes client to receive data from aeromodule

%%
%unpack some blade data from model
bladeData = model.bladeData;
numBlades = bladeData.numBlades;    %get number of blades
numEntriesH = length(bladeData.h);    %get number of entries in h for blade
numEntriesHPerBlade = numEntriesH/numBlades; %get number of points per blade 
hBlade = bladeData.h(1:numEntriesHPerBlade); %define an array of spanwise points for a blade

tsend = decodeVec(1,hBlade');
serverSendVector(tsend,d_output_stream)
disp('Interpolation data sent to aero code')


end
