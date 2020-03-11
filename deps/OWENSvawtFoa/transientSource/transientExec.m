function transientExec(model,mesh,el)
%transientExec performs modular transient analysis
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
%   transientExec(model,mesh,el)
%
%   %This function is an executable function for transient analysis. It
%   provides the interface of various external module with transient
%   structural dynamics analysis capability.
%
%   input:
%   model       = object containing model data
%   mesh        = object containing mesh data
%   el          = object containing element data
%
%
%   output: (NONE)



%% activate platform module
%............... flags for module activation ....................
model.aeroOn                 = false;

%modularIteration
moduleIteration = true;


if(model.aeroOn)
    aeroSendPort = 3200;
    aeroReceivePort = 4200;
    [d_output_streamAero,d_input_streamAero,server_socketAero,input_socketAero,output_socketAero] = delftAeroStartUp('localhost',aeroReceivePort,aeroSendPort,model);
end
%................................................................

%% Rotor mode initialization
%..........................................................................
OmegaInitial = model.OmegaInit; %Initial rotor speed (Hz)

if(model.turbineStartup == 1) %forced start-up using generator as motor
    disp('Running in forced starting mode.');
    model.generatorOn = true;
    %     Omega = OmegaInitial;
    rotorSpeedForGenStart = 0.0;
elseif(model.turbineStartup == 2) %self-starting mode
    disp('Running in self-starting mode.');
    model.generatorOn = false;
    %     Omega = OmegaInitial;
    rotorSpeedForGenStart = model.OmegaGenStart; %Spec rotor speed for generator startup Hz
else
    disp('Running in specified rotor speed mode');
    model.generatorOn = false;
    %     Omega = OmegaInitial;
    rotorSpeedForGenStart = 1e6; %ensures generator always off for practical purposes
end
%..........................................................................
% if(model.turbineStartup ==1 || model.turbineStartup==2)
%     plotGenSpeedVsTorque([0:.01:5],model.generatorProps);
% end
%%

%% state initialization
numDOFPerNode = 6;
model.totalNumDof = mesh.numNodes*numDOFPerNode;
%......... specify initial conditions .......................
u_s = zeros(model.totalNumDof,1);
[u_s] = setInitialConditions(model.initCond,u_s,numDOFPerNode);
u_sm1 = u_s;
udot_s = u_s*0;
uddot_s = u_s*0;
%............................................................

numTS = model.numTS;       %define number of time steps
delta_t = model.delta_t;   %define time step size
uHist(:,1) = u_s;          %store initial condition

%initialize omega_platform, omega_platform_dot, omegaPlatHist
% omega_platform = zeros(3,1);
% omega_platform_dot = zeros(3,1);
% omegaPlatHist(:,1) = omega_platform;

t = zeros(1,numTS+1);
FReactionHist = zeros(numTS+1,6);
% strainHist(numTS+1) = struct();
aziHist = zeros(1,numTS+1);
OmegaHist = zeros(1,numTS+1);
OmegaDotHist = zeros(1,numTS+1);
gbHist = zeros(1,numTS+1);
gbDotHist = zeros(1,numTS+1);
gbDotDotHist = zeros(1,numTS+1);
%genTorque = zeros(1,numTS+1);
genTorque = zeros(1,numTS+1);
genPower = zeros(1,numTS+1);
torqueDriveShaft = zeros(1,numTS+1);
Ywec = zeros(1,numTS+1);
rigidDof = zeros(1,numTS+1);

t(1) = 0.0; %initialize various states and variables
gb_s = 0;
gbDot_s = 0;
gbDotDot_s = 0;
azi_s = 0;
Omega_s = OmegaInitial;
OmegaDot_s = 0;
genTorque_s = 0;
torqueDriveShaft_s = 0;
% azi_sm1 = -Omega*delta_t*2*pi;
aziHist(1) = azi_s;
OmegaHist(1) = Omega_s;
OmegaDotHist(1) = OmegaDot_s;
FReactionsm1 = zeros(6,1);
FReactionHist(1,:) = FReactionsm1;
FReaction_j = FReactionsm1;
gbHist(1) = gb_s;
gbDotHist(1) = gbDot_s;
gbDotDotHist(1) = gbDotDot_s;
genTorque(1) = genTorque_s;
torqueDriveShaft(1) = torqueDriveShaft_s;
%%
try
    toc
catch
    tic
end

%% structural dynamics initialization
%..........................................................................
if(strcmp(model.analysisType,'ROM')) %initialize reduced order model
    %calculate constrained dof vector
    numDofPerNode = 6;
    isConstrained = zeros(model.totalNumDof,1);
    constDof = (model.BC.pBC(:,1)-1)*numDofPerNode + model.BC.pBC(:,2);
    index = 1;
    for i=1:mesh.numNodes
        for j=1:numDofPerNode
            if(ismember((i-1)*numDofPerNode + j,constDof))
                isConstrained(index) = 1;
            end
            index = index + 1;
        end
    end
    model.BC.isConstrained = isConstrained;


    [rom,elStorage]=reducedOrderModel(model,mesh,el,u_s); %construct reduced order model

    %set up inital values in modal space
    jointTransformTrans = model.jointTransform'; %'
    u_sRed = jointTransformTrans*u_s(1:end);
    udot_sRed = jointTransformTrans*udot_s(1:end);
    uddot_sRed = jointTransformTrans*uddot_s(1:end);

    BC = model.BC;
    [u_s2] = applyBCModalVec(u_sRed,BC.numpBC,BC.map);
    [udot_s2] = applyBCModalVec(udot_sRed,BC.numpBC,BC.map);
    [uddot_s2] = applyBCModalVec(uddot_sRed,BC.numpBC,BC.map);

    invPhi = rom.invPhi;

    eta_s     = invPhi*u_s2;
    etadot_s  = invPhi*udot_s2;
    etaddot_s = invPhi*uddot_s2;
else
    [elStorage] = initialElementCalculations(model,el,mesh); %perform initial element calculations for conventional structural dynamics analysis
end

%calculate structural/platform moi
[~,structureMOI,~]=calculateStructureMassProps(elStorage);
%..........................................................................

%% Main Loop - iterate for a solution at each time step, i
for i=1:numTS

    %     i %TODO add verbose printing
    if(mod(i,100)==0) %print command that displays progress of time stepping
        fprintf('%s\n',['Iteration: ' i])
    end

    %% check for specified rotor speed at t(i) + delta_t
    model.omegaControl = false;
    if(model.turbineStartup == 0)
        model.omegaControl = true;
        if(model.usingRotorSpeedFunction) %use user specified rotor speed profile function
            [~,omegaCurrent,~] = getRotorPosSpeedAccelAtTime(t(i),t(i)+delta_t,0.0,delta_t);
            Omega_s = omegaCurrent;
        else %use discreteized rotor speed profile function
            [omegaCurrent,OmegaDotCurrent,terminateSimulation] = omegaSpecCheck(t(i)+delta_t,model.tocp,model.Omegaocp,delta_t);
            if(terminateSimulation)
                break;
            end
            Omega_s = omegaCurrent;
            OmegaDot_s = OmegaDotCurrent;
        end
    end
    %%

    %% initialize "j" Gauss-Sidel iteration
    u_j=u_s;
    azi_j = azi_s;
    Omega_j = Omega_s;
    OmegaDot_j = OmegaDot_s;
    gb_j = gb_s;
    gbDot_j = gbDot_s;
    gbDotDot_j = gbDotDot_s;
    genTorque_j = genTorque_s;

    needsAeroCalcAtThisTimestep = true;

    if(model.hydroOn)  %initialize  platform module related variables
        Ywec_j = Ywec(i,:);
        Ywec_jLast = Ywec_j;
        % 		Accel_j = Accel;
        % 		Accel_jLast = Accel;
    end

    TOL = 1e-8;  %gauss-seidel iteration tolerance for various modules
    MAXITER = 50; %max iteration for various modules
    numIterations = 1; uNorm = 1e6; platNorm = 1e6; aziNorm = 1e6; gbNorm = 1e6; %initialize norms for various module states
    %%

    while((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER)) %module gauss-seidel iteration loop

        rbData = zeros(1,9);
        %calculate CP2H (platform frame to hub frame transformation matrix)
        CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1];

        %.... inertial frame to platform transformation matrix ...........
        if(model.hydroOn)
            CN2P=calculateLambdaSlim(Ywec_j(s_stsp+6),Ywec_j(s_stsp+5),Ywec_j(s_stsp+4));
        else
            CN2P=eye(3);
        end
        %.........................................

        CN2H = CP2H*CN2P;

%         %% evaluate platform module
%         %%====================================================
%         if(model.hydroOn)
%             % 	Accel_jLast= Accel_j;
%             Ywec_jLast = Ywec_j;
%             if(model.platformTurbineYawInteraction == 0)
%                 FReaction0 = [-FReactionHist(i,1:5)'; 0.0]; %'
%                 FReaction1 =  [-FReaction_j(1:5); 0.0];
%             elseif(model.platformTurbineYawInteraction == 1)
%                 FReaction0 = (-FReactionHist(i,1:6)');
%                 FReaction1 =  (-FReaction_j(1:6));
%             elseif(model.platformTurbineYawInteraction == 2)
%                 FReaction0 = [-FReactionHist(i,1:5)'; genTorque_s];
%                 FReaction1 =  [-FReaction_j(1:5); genTorque_j];
%             else
%                 error('PlatformTurbineYawInteraction flag not recognized.');
%             end
%             [rbData,Ywec_j,~] = platformModule([t(i) t(i)+delta_t],Ywec(i,:),CP2H,FReaction0,FReaction1,d_input_streamPlatform,d_output_streamPlatform);
%         end
%         %====================================================
        %%

        %% evaluate generator module
        %===== generator module ===========================
        genTorque_j = 0;
        if (model.generatorOn)
            if(model.driveTrainOn)
                if(model.useGeneratorFunction)
                    [genTorqueHSS0] = userDefinedGenerator(gbDot_j*model.gearRatio);
                else
                    [genTorqueHSS0] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio);
                end
            else
                if(model.useGeneratorFunction)
                    [genTorqueHSS0] = userDefinedGenerator(Omega_j);
                else
                    [genTorqueHSS0] = simpleGenerator(model.generatorProps,Omega_j);
                end
            end
            %should eventually account for Omega = gbDot*gearRatio here...
            genTorque_j = genTorqueHSS0*model.gearRatio*model.gearBoxEfficiency; %calculate generator torque on LSS side
            %         genTorqueAppliedToTurbineRotor0 = -genTorque0;
            %         genTorqueAppliedToPlatform0 = genTorqueHSS0;
        end
        %==================================================
        %%

        %% evaluate drivetrain module
        % %===== drivetrain module ==========================
        torqueDriveShaft_j = genTorque_j;
        gb_jLast = gb_j;
        if(~model.omegaControl)
            if(model.driveTrainOn)
                [torqueDriveShaft_j] = calculateDriveShaftReactionTorque(model.driveShaftProps,...
                    azi_j,gb_j,Omega_j*2*pi,gbDot_j*2*pi);

                [gb_j,gbDot_j,gbDotDot_j] = updateRotorRotation(model.JgearBox,0,0,...
                    -genTorque_j,torqueDriveShaft_j,...
                    gb_s,gbDot_s,gbDotDot_s,delta_t);
            else
                gb_j = azi_j;
                gbDot_j = Omega_j;
                gbDotDot_j = OmegaDot_j;
            end
        else
            gb_j = azi_j;
            gbDot_j = omegaCurrent*2*pi;
            gbDotDot_j = 0;
        end
        % %==================================================

        %% rotor speed update
        %===== update rotor speed =========================
        azi_jLast = azi_j;
        if(model.omegaControl)
            if(model.usingRotorSpeedFunction)
                [azi_j,Omega_j,OmegaDot_j] = getRotorPosSpeedAccelAtTime(t(i),t(i)+delta_t,azi_s,delta_t);
            else
                Omega_j = Omega_s;
                OmegaDot_j = OmegaDot_s;
                azi_j = azi_s + Omega_j*delta_t*2*pi;
            end
        end
        if(~model.omegaControl)
            Crotor = 0;
            Krotor = 0;
            [azi_j,Omega_j,OmegaDot_j] = updateRotorRotation(structureMOI(3,3),Crotor,Krotor,...
                -FReaction_j(6),-torqueDriveShaft_j,...
                azi_s,Omega_s,OmegaDot_s,delta_t);
        end
        %===================================================
        %%

        %% evaluate aerodynamic module (CACTUS ONE-WAY)
        %%======= aerodynamics module ======================
        if(model.aeroLoadsOn) % TODO: this is odd since we load in the aero loads elsewhere
%             disp('using cactus aero loads...');
%             %----- calculate aerodynamic loads --------------------------------
%             % this is place holder
%             [FAero] = getAeroLoads(model.bladeData,model.aeroloadfile,t(i),el.props,model.totalNumDof);
%             % these loads will need to account for a turbine with a
%             % platform with different orientation than the inertial system.
%             FAeroDof = (1:length(u_s));
%             %------------------------------------------------------------------
        else
            FAero = [];
            FAeroDof = [];
        end
        %==================================================
        %% evaluate aerodynamic module (TU DELFT)
        %%======= aerodynamics module ======================
        if(model.aeroOn && needsAeroCalcAtThisTimestep)

            [FAero,FAeroDof] = aeroModule(model,t(i) + delta_t,u_j,Omega_j,azi_j,numDOFPerNode,d_input_streamAero,d_output_streamAero);

            %set aero forces flag
            needsAeroCalcAtThisTimestep = false;

        end
        %==================================================
        %%

        %% compile external forcing on rotor
        %compile forces to supply to structural dynamics solver
        [Fexternal, Fdof] = externalForcing(t(i)+delta_t);
        Fexternal = [Fexternal; FAero];
        Fdof = [Fdof, FAeroDof];


        %% evaluate structural dynamics
        %call structural dynamics solver
        if(strcmp(model.analysisType,'TD'))  %initialization of structural dynamics displacements, velocities, accelerations, etc.
            dispData.displ_s = u_s;
            dispData.displ_sm1 = u_sm1;
        end

        if(strcmp(model.analysisType,'TNB'))
            dispData.displ_s = u_s;
            dispData.displdot_s = udot_s;
            dispData.displddot_s = uddot_s;
        end

        if(strcmp(model.analysisType,'ROM'))
            dispData.displ_s = u_s;
            dispData.displdot_s = udot_s;
            dispData.displddot_s = uddot_s;

            dispData.eta_s     = eta_s;
            dispData.etadot_s  = etadot_s;
            dispData.etaddot_s = etaddot_s;
        end

        if(strcmp(model.analysisType,'ROM'))
            % evalulate structural dynamics using reduced order model
            [dispOut,FReaction_j] = structuralDynamicsTransientROM(model,mesh,el,dispData,Omega_j,OmegaDot_j,t(i),delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData);
        else
            % evalulate structural dynamics using conventional representation
            [dispOut,FReaction_j] = structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t(i),delta_t,elStorage,Fexternal,Fdof,CN2H,rbData);
        end
        %update last iteration displacement vector
        u_jLast = u_j;
        u_j = dispOut.displ_sp1;             %update current estimates of velocity, acceleration
        if(strcmp(model.analysisType,'TNB'))
            udot_j  = dispOut.displdot_sp1;
            uddot_j = dispOut.displddot_sp1;
        end

        if(strcmp(model.analysisType,'ROM'))
            udot_j  = dispOut.displdot_sp1;
            uddot_j = dispOut.displddot_sp1;

            eta_j = dispOut.eta_sp1;
            etadot_j = dispOut.etadot_sp1;
            etaddot_j = dispOut.etaddot_sp1;
        end
        %%

        %% calculate norms
        uNorm = norm(u_j-u_jLast)/norm(u_j);            %structural dynamics displacement iteration norm
        aziNorm = norm(azi_j - azi_jLast)/norm(azi_j);  %rotor azimuth iteration norm

        if(model.hydroOn)
            platNorm = norm(Ywec_j-Ywec_jLast)/norm(Ywec_j); %platform module states iteration norm
        else
            platNorm = 0.0;
        end

        if(model.driveTrainOn)
            gbNorm = norm(gb_j - gb_jLast)/norm(gb_j); %gearbox states iteration norm
        else
            gbNorm = 0.0;
        end

        if(moduleIteration == false)
            break;
        end

        numIterations = numIterations + 1;
    end %end iteration while loop

    %% calculate converged generator torque/power
    if(~isempty(model.generatorProps))
        if(model.generatorOn || (model.turbineStartup==0))
            [genTorquePlot] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio);
        else
            genTorquePlot = 0;
        end
    else
        genTorquePlot = 0;
    end
    [genPowerPlot] = genTorquePlot*(gbDot_j*2*pi)*model.gearRatio;


    %% update timestepping variables and other states, store in history arrays
    if(strcmp(model.analysisType,'TD'))
        u_sm1 = u_s;
        u_s = u_j;
    end
    if(strcmp(model.analysisType,'TNB'))
        u_s = u_j;
        udot_s = udot_j;
        uddot_s = uddot_j;
    end
    if(strcmp(model.analysisType,'ROM'))
        u_s = u_j;
        udot_s = udot_j;
        uddot_s = uddot_j;

        eta_s = eta_j;
        etadot_s = etadot_j;
        etaddot_s = etaddot_j;
    end

    uHist(:,i+1) = u_s;
    FReactionHist(i+1,:) = FReaction_j;
    strainHist(:,i) = dispOut.elStrain;
%     strainHist(:,i).eps_xx_0 = temp.eps_xx_0;
%     strainHist(:,i).eps_xx_z = temp.eps_xx_z;
%     strainHist(:,i).eps_xx_y = temp.eps_xx_y;
%     strainHist(:,i).gam_xz_0 = temp.gam_xz_0;
%     strainHist(:,i).gam_xz_y = temp.gam_xz_y;
%     strainHist(:,i).gam_xy_0 = temp.gam_xy_0;
%     strainHist(:,i).gam_xy_z = temp.gam_xy_z;
    t(i+1) = t(i) + delta_t;

    azi_s = azi_j;
    Omega_s = Omega_j;
    OmegaDot_s = OmegaDot_j;

    genTorque_s = genTorque_j;
    torqueDriveShaft_s = torqueDriveShaft_j;

    aziHist(i+1) = azi_s;
    OmegaHist(i+1) = Omega_s;
    OmegaDotHist(i+1) = OmegaDot_s;

    gb_s = gb_j;
    gbDot_s = gbDot_j;
    gbDotDot_s = gbDotDot_j;

    gbHist(i+1) = gb_s;
    gbDotHist(i+1) = gbDot_s;
    gbDotDotHist(i+1) = gbDotDot_s;

    %genTorque(i+1) = genTorque_s;
    genTorque(i+1) = genTorquePlot;
    genPower(i+1) = genPowerPlot;
    torqueDriveShaft(i+1) = torqueDriveShaft_s;

    if(model.hydroOn)
        error('Hydro Model not fully implemented');
        %         Ywec(i+1,:) = Ywec_j;
        %         rigidDof(i+1,:)=Ywec(i+1,s_stsp+1:s_stsp+6);

    else
        rigidDof(i) = 0;
    end

    FReactionHist(i+1,:) = FReaction_j;
    %%

    %% check rotor speed for generator operation
    if(Omega_s>= rotorSpeedForGenStart)
        model.generatorOn = true;
    else
        model.generatorOn = false;
    end
    %%

end %end timestep loop

%% kill platform module process
if(model.hydroOn)
    serverSendVector(-1.0,d_output_streamPlatform);
    terminateServer(server_socketPlatform,output_socketPlatform,1);
    terminateClient(input_socketPlatform,1);
end

%% kill aerodynamic module process
if(model.aeroOn)
    serverSendVector(decodeVec(4,-1.0),d_output_streamAero) ;
    terminateClient(input_socketAero,1); %close down client connection to forcing module
    terminateServer(server_socketAero,output_socketAero,1); %close down server on this side
end

%%
%toc
% save aeroOutputArray
%save simulation data in .mat file
save(model.outFilename,'t','uHist','aziHist','OmegaHist','OmegaDotHist','gbHist','gbDotHist','gbDotDotHist','FReactionHist','rigidDof','genTorque','genPower','torqueDriveShaft','strainHist');

end

function [OmegaCurrent,OmegaDotCurrent,terminateSimulation] = omegaSpecCheck(tCurrent,tocp,Omegaocp,delta_t)

if(tocp(length(tocp))<tCurrent)
    disp('Simulation time is greater than that specified in control points for prescribed rotor speed profile.');
    disp('Terminating simulation.');
    terminateSimulation = true;
    OmegaCurrent = [];
    OmegaDotCurrent = [];
else
    OmegaCurrent = interp1(tocp,Omegaocp,tCurrent); %interpolated discreteized profile for current omega

    %calculate current rotor acceleration
    dt = delta_t/2.0;
    omega_m1 = interp1(tocp,Omegaocp,tCurrent-dt);
    omega_p1 = interp1(tocp,Omegaocp,tCurrent+dt);

    OmegaDotCurrent = diff([omega_m1,omega_p1])/(dt*2);

    terminateSimulation = false;
    if(isnan(OmegaCurrent) || isnan(OmegaDotCurrent))
        error('Omega calcualted a NaN. Exiting.');
    end
end

end
