function [Ywec,rigidDof,Accel,s_stsp,d_output_stream,d_input_stream,server_socket,input_socket,output_socket] = waveECStartUp(host,port,output_port,model)
%wavECStartUp   starts up wave ec module
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Ywec,rigidDof,Accel,s_stsp,d_output_stream,d_input_stream,...
%    server_socket,input_socket,output_socket] = waveECStartUp(host,port,...
%                                                output_port,model)
%                    
%   This function starts up the wav ec module from the OWENS toolkit
%
%   input:
%   host            = host address for client connectino
%   port            = port number for client connection
%   output_port     = port number for server connection
%   model           = model object
%
%   output:
%   Ywec            = initialization of platform module states
%   rigidDof        = initialzation of rigid body degrees of freedom
%   Accel           = initialization of platform accelerations
%   s_stsp          = size of convolution state space in platform module
%   d_output_stream = output stream object for server socket connections
%   d_input_stream  = input stream object for client socket connections
%   server_socket   = server object for server socket connections
%   input_socket    = input object for client socket connections

    number_of_retries = 1000;
%......................................
%%
    import java.net.ServerSocket
    import java.net.Socket
    import java.io.*
    
    %sets system commands to launch wavEC software
if 0
    hydroCodeExec = 'hello';
    hydroCodeDirectory = 'C:\Users\blennis\Desktop\AEP';
    matlabPath = 'matlab.exe';
    sysCallString = [matlabPath,' -sd ',hydroCodeDirectory,' -r ',hydroCodeExec,' &'];
    hm = pwd; cd(hydroCodeDirectory)
        system(sysCallString);
cd(hm)
    keyboard
    
else
    disp('launching WAVEC');    %launch wavEC depending on pc or unix environment
    if(ispc());
        hydroCodeExec = 'TimeDModel_socket_V2';
        hydroCodeDirectory = 'C:\DesignCodes\WavEC2Wire_v1.0';
        matlabPath = 'matlab.exe';
        sysCallString = [matlabPath,' -sd ',hydroCodeDirectory,' -r ',hydroCodeExec,' &'];
        hm = pwd; cd(hydroCodeDirectory);
        system(sysCallString);
        cd(hm);
    end

    if(isunix())
        %     hydroLaunchScript = '/home/bcowens/work/OWENSsingle/launchPlatformCode';
        %     sysCallString = ['bash ',hydroLaunchScript,' &'];
        unix(sysCallString);
    end
%     keyboard
    pause(5);
end

[d_output_stream,server_socket,output_socket] = initServer(output_port,number_of_retries); %initializes server to send data to wavEC
[d_input_stream,input_socket] = initClient(host,port,number_of_retries); %initializes client to receive data from wavEC

%%
%
% Notes on Input Files
% The WAMIT files are stored in \FreqData\ in the run directory.  They
% should not need to be modified by you for your tests.  They will need to
% be updated as the platform changes.
%
% The input files to run WavEC are located in the parent directory
%[Ywec,rigigDof,Accel,s_stsp,d_output_stream,d_input_stream] = waveECStartUp(host,port,output_port)
% LAYOUT.log
% - Contains x,y location of device, not necessary to modify
%
% MOORINGS.log
% - Contains coefficients for mooring forces, disabled in the current code
% - Mooring coefficients are hard-coded in TimeDModel for now, will be
%   included in a mooring.m module later
%
% PTO.log
% - Contains information on the PTO system, disabled in the current code
%
% ExcINPUT.log
% - Used to control what type of sea state you want.  
% - To disable waves, use a Regular Wave (option 4) with Hs = 0m
%
% RadINPUT.log
% - Contains radiation information for formation of the convolution.
% - The first line is where you define active DOF (ex. [1 1 1 1 1 1 0 0])
%   There must be 8 entries (to match WAMIT analysis), so use the 
%   first 6 for your rigid body and leave the last two as 0
% - The simulation time is also specified here, but this is overwritten in
%   the code to be the value OWENS sends through the socket.

%% Initialization Parameters
% These need to be defined by OWENS

% Setup total time vector for simulation
% Required by WavEC to determine length of time to generate waves
%time_owens=[0:.1:50]; % Constant dt
%    time_owens=[0:.05:10,10.01:.01:20,20.05:.05:50]; % Variable dt

% Setup flags for WavEC [1]=On, [0] = Off
% These are used to control various aspects of WavEC from OWENS
% Note: mooring flag must be [1] to do surge/sway tests
Drag_Flag = model.Plat_Drag_Flag;  % Activate Drag time_owens(1) time_owens(end)damping
Moor_Flag = model.Plat_Moor_Flag;  % Activate Mooring
Grav_Flag = model.Plat_Grav_Flag;  % Include gravity?
Plot_Flag = model.Plat_Plot_Flag;  % Plot Results in WavEC
RadDamp_Flag = model.Plat_RadDamp_Flag; %Rad damping (default 1)

Flags = [Drag_Flag; Moor_Flag; Grav_Flag; Plot_Flag; RadDamp_Flag;0]; %compile list of flags

%send active platform DOFs
serverSendVector(model.activePlatformDofs,d_output_stream) ;

% send to WavEC
serverSendVector([0.0 model.numTS*model.delta_t],d_output_stream) ;
serverSendVector(Flags,d_output_stream) ;

%% Get size of convolution vector, #dof, and padded initial condition
s_stsp = clientRecvVector(d_input_stream); % Size of convolution st-sp
s_dof = clientRecvVector(d_input_stream);  % Num of active dof 
s_pad_sol = clientRecvVector(d_input_stream); % Size of padded sol vector

% Initial Conditions Y = [zeros(s_stsp) [6DOF Pos IC] [6DOF Vel IC]]
% Use this line to define your initial (t=0) conditions
% Initialize the first [1 x s_stsp] convolutions with 0
% The next 6DOF are surge, sway, heave, roll, pitch, yaw displacements
% The last 6DOF are surge, sway, heave, roll, pitch, yaw velocities
Ywec(1,:)=[zeros(1,s_stsp) model.initialCondPlatformDofs [0 0 0 0 0 0]]; % Initial Conditions [1 x s_pad_sol]            %problems with bounded rotations ... no resistance for now
rigidDof(1,:)=Ywec(s_stsp+1:s_stsp+6);
% Placeholder as first solution will be at t=dt
Accel(1,:)=zeros(1,6); % Placeholder for t=0 [1 x 6] 

end
