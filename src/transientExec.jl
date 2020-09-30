
include("externalForcing.jl")
include("setInitialConditions.jl")
include("call_structuralDynamicsTransient.jl")

mutable struct ElStrain
    eps_xx_0
    eps_xx_z
    eps_xx_y
    gam_xz_0
    gam_xz_y
    gam_xy_0
    gam_xy_z
end

mutable struct DispOut
    elStrain
    displ_sp1
    displddot_sp1
    displdot_sp1
end

mutable struct DispData
    displ_s
    displdot_s
    displddot_s
end

mutable struct Slice
    r
    RefR
    chord
    twist
    delta
    B
    alpha_rad
    cl
    cd
    Omega
    centerX
    centerY
end

mutable struct Env
    rho
    mu
    G_amp
    gusttime
    gustX0
    N_Rev
    idx_sub
    wsave
    Vinf_nominal
    V_vert
    V_tang
    V_rad
    V_twist
    Vinf
    V_wake_old
    steplast
    ntheta
    tau
end

function transientExec(model,mesh,el)
    #transientExec performs modular transient analysis
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #
    #   transientExec(model,mesh,el)
    #
    #   #This function is an executable function for transient analysis. It
    #   provides the interface of various external module with transient
    #   structural dynamics analysis capability.
    #
    #   input:
    #   model       = object containing model data
    #   mesh        = object containing mesh data
    #   el          = object containing element data
    #
    #
    #   output: (NONE)



    ## activate platform module
    #............... flags for module activation ....................
    aeroOn = false
    #modularIteration
    moduleIteration = true
    # Get AeroLoads
    # mat"$aeroLoads = processAeroLoadsBLE($model.aeroloadfile, $model.owensfile)"
    aeroLoadsFile_root = model.aeroloadfile[1:end-16] #cut off the _ElementData.csv
    OWENSfile_root = model.owensfile[1:end-6] #cut off the .owens

    QCx=[0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00]
    QCy=[0.00000e+00,2.44680e-01,4.89360e-01,7.34040e-01,9.78720e-01,1.22340e+00,1.46808e+00,1.71276e+00,1.95744e+00,2.20212e+00,2.44680e+00]
    QCz=[0.00000e+00,3.76945e-01,6.77076e-01,8.73458e-01,9.77684e-01,1.00000e+00,9.44658e-01,8.13441e-01,6.09216e-01,3.25362e-01,0.00000e+00]
    CtoR=[1.11093e-01,1.11093e-01,8.20390e-02,6.35940e-02,5.68145e-02,5.55467e-02,5.88008e-02,6.82860e-02,9.11773e-02,1.11093e-01,1.11093e-01]
    RefR=177.2022
    NBlade = 2
    PEy=[1.22340e-01,3.67020e-01,6.11700e-01,8.56380e-01,1.10106e+00,1.34574e+00,1.59042e+00,1.83510e+00,2.07978e+00,2.32446e+00]
    NElem = 10

    mat"[$structuralSpanLocNorm,$structuralNodeNumbers,$structuralElNumbers] = readBldFile($OWENSfile_root)"

    step_AC = 0
    ntheta = 36#1/(model.OmegaInit*model.delta_t)
    Vinf = 25 #TODO
    rho = 1.225 #TODO
    mu = 1.7894e-5 #TODO

    H = [] #For plotting turbine
    plotTurbine = false

    mat"[$alpha_rad,$cl,$cd] = readaerodyn('airfoils/NACA_0015_RE3E5.dat')"
    ft2m = 1 / 3.281
    RefR = RefR*ft2m

    xyz = zeros(length(QCz),3)
    xyz[:,1] = RefR*QCz
    xyz[:,2] = RefR*QCx
    xyz[:,3] = RefR*QCy

    n_slices = length(QCz)-1

    # h_ends = xyz(:,3)
    # h_frac = (h_ends(2:end) - h_ends(1:end-1))/h_ends(end)
    # h = (h_ends(2:end) + h_ends(1:end-1))/2

    delta_xs = xyz[2:end,1] - xyz[1:end-1,1]
    delta_zs = xyz[2:end,3] - xyz[1:end-1,3]

    # element_planf_A = sqrt(delta_xs.^2+delta_zs.^2)*chord
    # element_planf_L = sqrt(delta_xs.^2+delta_zs.^2)

    delta = atan.(delta_xs./delta_zs)

    r = (xyz[2:end,1]+xyz[1:end-1,1])/2
    twist = ones(n_slices)*0*pi/180
    chord = RefR*CtoR
    # Single Slice
    slice = Slice(ones(ntheta),
        RefR,
        0.0,
        zeros(1,ntheta),
        zeros(1,ntheta),
        Float64(NBlade),
        alpha_rad,
        cl,
        cd,
        ones(1,ntheta)*model.OmegaInit*2*pi,
        0.0,
        0.0)

    # turbine built from bottom up
    turbine3D = fill(slice, n_slices)
    for i = 1:n_slices
        turbine3D[i].r = ones(1,ntheta)*r[i]
        turbine3D[i].chord = chord[i]
        turbine3D[i].twist = ones(1,ntheta)*twist[i]
        turbine3D[i].delta = ones(1,ntheta)*delta[i]
    end

    env1 = Env(rho,
        mu,
        0.0, #m/s gust
        0.8, #sec
        13.0, #radaii offset back is starting point for gust
        20,
        zeros(1,NBlade*2),
        zeros(1,ntheta*2),
        Vinf,
        zeros(1,ntheta),
        zeros(1,ntheta),
        zeros(1,ntheta),
        zeros(1,ntheta),
        ones(ntheta)*Vinf,
        Vinf,
        0,
        ntheta,
        [0.3025, 2.9500])

    env = fill(env1, n_slices)
    # Declare Variable Type, are set later
    udot_j = 0.0
    uddot_j = 0.0
    torqueDriveShaft_j = 0.0

    elStrain = fill(ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)

    dispOut = DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492))
    #................................................................

    ## Rotor mode initialization
    #..........................................................................
    OmegaInitial = model.OmegaInit #Initial rotor speed (Hz)

    if (model.turbineStartup == 1) #forced start-up using generator as motor
        println("Running in forced starting mode.")
        model.generatorOn = true  #TODO: clean this redundant/conflicting logic up
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = 0.0
    elseif (model.turbineStartup == 2) #self-starting mode
        println("Running in self-starting mode.")
        model.generatorOn = false
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = model.OmegaGenStart #Spec rotor speed for generator startup Hz
    else
        println("Running in specified rotor speed mode")
        model.generatorOn = false
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = 1e6 #ensures generator always off for practical purposes
    end
    #..........................................................................
    # if (model.turbineStartup ==1 || model.turbineStartup==2)
    #     plotGenSpeedVsTorque([0:.01:5],model.generatorProps)
    # end
    ##
    ## state initialization
    numDOFPerNode = 6
    model.totalNumDof = mesh.numNodes*numDOFPerNode
    #......... specify initial conditions .......................
    u_s = zeros(model.totalNumDof,1)
    u_s = setInitialConditions(model.initCond,u_s,numDOFPerNode)
    u_sm1 = u_s
    udot_s = u_s*0
    uddot_s = u_s*0
    #............................................................

    numTS = Int(model.numTS)       #define number of time steps
    delta_t = model.delta_t   #define time step size
    uHist = zeros(length(u_s),numTS+1)
    uHist[:,1] = u_s          #store initial condition

    #initialize omega_platform, omega_platform_dot, omegaPlatHist
    # omega_platform = zeros(3,1)
    # omega_platform_dot = zeros(3,1)
    # omegaPlatHist(:,1) = omega_platform

    t = zeros(numTS+1)
    FReactionHist = zeros(numTS+1,6)
    # strainHist(numTS+1) = struct()
    strainHist = fill(ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)),mesh.numEl,numTS)

    aziHist = zeros(numTS+1)
    OmegaHist = zeros(numTS+1)
    OmegaDotHist = zeros(numTS+1)
    gbHist = zeros(numTS+1)
    gbDotHist = zeros(numTS+1)
    gbDotDotHist = zeros(numTS+1)
    genTorque = zeros(numTS+1)
    genPower = zeros(numTS+1)
    torqueDriveShaft = zeros(numTS+1)
    Ywec = zeros(numTS+1,1)
    rigidDof = zeros(numTS+1)

    t[1] = 0.0 #initialize various states and variables
    gb_s = 0
    gbDot_s = 0
    gbDotDot_s = 0
    azi_s = 0
    Omega_s = OmegaInitial
    OmegaDot_s = 0
    genTorque_s = 0
    torqueDriveShaft_s = 0
    # azi_sm1 = -Omega*delta_t*2*pi
    aziHist[1] = azi_s
    OmegaHist[1] = Omega_s
    OmegaDotHist[1] = OmegaDot_s
    FReactionsm1 = zeros(6)
    FReactionHist[1,:] = FReactionsm1
    FReaction_j = FReactionsm1
    gbHist[1] = gb_s
    gbDotHist[1] = gbDot_s
    gbDotDotHist[1] = gbDotDot_s
    genTorque[1] = genTorque_s
    torqueDriveShaft[1] = torqueDriveShaft_s
    ##

    ## structural dynamics initialization
    #..........................................................................
    if (occursin("ROM",model.analysisType)) #initialize reduced order model
        error("ROM not fully implemented")
        #     #calculate constrained dof vector
        #     numDofPerNode = 6
        #     isConstrained = zeros(model.totalNumDof,1)
        #     constDof = (model.BC.pBC(:,1)-1)*numDofPerNode + model.BC.pBC(:,2)
        #     index = 1
        #     for i=1:mesh.numNodes
        #         for j=1:numDofPerNode
        #             if (ismember((i-1)*numDofPerNode + j,constDof))
        #                 isConstrained(index) = 1
        #             end
        #             index = index + 1
        #         end
        #     end
        #     model.BC.isConstrained = isConstrained
        #
        #
        #     [rom,elStorage]=reducedOrderModel(model,mesh,el,u_s) #construct reduced order model
        #
        #     #set up inital values in modal space
        #     jointTransformTrans = model.jointTransform' #'
        #     u_sRed = jointTransformTrans*u_s(1:end)
        #     udot_sRed = jointTransformTrans*udot_s(1:end)
        #     uddot_sRed = jointTransformTrans*uddot_s(1:end)
        #
        #     BC = model.BC
        #     [u_s2] = applyBCModalVec(u_sRed,BC.numpBC,BC.map)
        #     [udot_s2] = applyBCModalVec(udot_sRed,BC.numpBC,BC.map)
        #     [uddot_s2] = applyBCModalVec(uddot_sRed,BC.numpBC,BC.map)
        #
        #     invPhi = rom.invPhi
        #
        #     eta_s     = invPhi*u_s2
        #     etadot_s  = invPhi*udot_s2
        #     etaddot_s = invPhi*uddot_s2
    else
        mat"$elStorage = initialElementCalculations($model,$el,$mesh)" #perform initial element calculations for conventional structural dynamics analysis
    end

    #calculate structural/platform moi
    mat"[$unused1,$structureMOI,$unused2]=calculateStructureMassProps($elStorage)"
    #..........................................................................

    ## Main Loop - iterate for a solution at each time step, i
    for i=1:numTS

        #     i #TODO add verbose printing
        if (mod(i,100)==0) #print command that displays progress of time stepping
            println("Iteration: $i")
        end

        ## check for specified rotor speed at t[i] + delta_t
        model.omegaControl = false
        if (model.turbineStartup == 0)
            model.omegaControl = true
            if (model.usingRotorSpeedFunction) #use user specified rotor speed profile function
                _,omegaCurrent,_ = getRotorPosSpeedAccelAtTime(t[i],t[i]+delta_t,0.0,delta_t)
                Omega_s = omegaCurrent
            else #use discreteized rotor speed profile function
                omegaCurrent,OmegaDotCurrent,terminateSimulation = omegaSpecCheck(t[i]+delta_t,model.tocp,model.Omegaocp,delta_t)
                if (terminateSimulation)
                    break
                end
                Omega_s = omegaCurrent
                OmegaDot_s = OmegaDotCurrent
            end
        else
            omegaCurrent = 0.0
        end
        ##

        ## initialize "j" Gauss-Sidel iteration
        u_j=u_s
        azi_j = azi_s
        Omega_j = Omega_s
        OmegaDot_j = OmegaDot_s
        gb_j = gb_s
        gbDot_j = gbDot_s
        gbDotDot_j = gbDotDot_s
        genTorque_j = genTorque_s

        needsAeroCalcAtThisTimestep = true

        #initialize  platform module related variables only used if (model.hydroOn)
        Ywec_j = Ywec[i,:]
        Ywec_jLast = Ywec_j
        # 		Accel_j = Accel
        # 		Accel_jLast = Accel

        TOL = 1e-8  #gauss-seidel iteration tolerance for various modules
        MAXITER = 50 #max iteration for various modules
        numIterations = 1
        uNorm = 1e6
        platNorm = 1e6
        aziNorm = 1e6
        gbNorm = 1e6 #initialize norms for various module states
        ##

        while ((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER)) #module gauss-seidel iteration loop

            rbData = zeros(9)
            #calculate CP2H (platform frame to hub frame transformation matrix)
            CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]

            #.... inertial frame to platform transformation matrix ...........
            if (model.hydroOn)
                CN2P=1.0*LinearAlgebra.I(3) # temporary declaration: original CN2P below
                #             CN2P=calculateLambdaSlim(Ywec_j(s_stsp+6),Ywec_j(s_stsp+5),Ywec_j(s_stsp+4))
            else
                CN2P=1.0*LinearAlgebra.I(3)
            end
            #.........................................

            CN2H = CP2H*CN2P

            #         ## evaluate platform module
            #         ##-------------------------------------
            #         if (model.hydroOn)
            #             # 	Accel_jLast= Accel_j
            #             Ywec_jLast = Ywec_j;
            #             if (model.platformTurbineYawInteraction == 0)
            #                 FReaction0 = [-FReactionHist(i,1:5)'; 0.0]; #'
            #                 FReaction1 =  [-FReaction_j(1:5); 0.0];
            #             elseif (model.platformTurbineYawInteraction == 1)
            #                 FReaction0 = (-FReactionHist(i,1:6)');
            #                 FReaction1 =  (-FReaction_j(1:6));
            #             elseif (model.platformTurbineYawInteraction == 2)
            #                 FReaction0 = [-FReactionHist(i,1:5)'; genTorque_s];
            #                 FReaction1 =  [-FReaction_j(1:5); genTorque_j];
            #             else
            #                 error('PlatformTurbineYawInteraction flag not recognized.');
            #             end
            #             [rbData,Ywec_j,_] = platformModule([t[i] t[i]+delta_t],Ywec(i,:),CP2H,FReaction0,FReaction1,d_input_streamPlatform,d_output_streamPlatform);
            #         end
            #         #-------------------------------------
            ##

            ## evaluate generator module
            #----- generator module ---------------------------
            genTorque_j = 0
            if (model.generatorOn)
                if (model.driveTrainOn)
                    if (model.useGeneratorFunction)
                        genTorqueHSS0 = userDefinedGenerator(gbDot_j*model.gearRatio)
                    else
                        error("simpleGenerator not fully implemented")#[genTorqueHSS0] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio)
                    end
                else
                    if (model.useGeneratorFunction)
                        genTorqueHSS0 = userDefinedGenerator(Omega_j)
                    else
                        error("simpleGenerator not fully implemented")#[genTorqueHSS0] = simpleGenerator(model.generatorProps,Omega_j)
                    end
                end
                #should eventually account for Omega = gbDot*gearRatio here...
                genTorque_j = genTorqueHSS0*model.gearRatio*model.gearBoxEfficiency #calculate generator torque on LSS side
                #         genTorqueAppliedToTurbineRotor0 = -genTorque0
                #         genTorqueAppliedToPlatform0 = genTorqueHSS0
            end
            #-------------------------------------       ##

            ## evaluate drivetrain module
            # #------ drivetrain module ---------------------------------
            torqueDriveShaft_j = genTorque_j
            gb_jLast = gb_j
            if (!model.omegaControl)
                if (model.driveTrainOn)
                    torqueDriveShaft_j = calculateDriveShaftReactionTorque(model.driveShaftProps,
                    azi_j,gb_j,Omega_j*2*pi,gbDot_j*2*pi)

                    gb_j,gbDot_j,gbDotDot_j = updateRotorRotation(model.JgearBox,0,0,
                    -genTorque_j,torqueDriveShaft_j,gb_s,gbDot_s,gbDotDot_s,delta_t)
                else
                    gb_j = azi_j
                    gbDot_j = Omega_j
                    gbDotDot_j = OmegaDot_j
                end
            else
                gb_j = azi_j
                gbDot_j = omegaCurrent*2*pi
                gbDotDot_j = 0
            end
            #-------------------------------------        ## rotor speed update
            #------ update rotor speed ---------------------------------
            azi_jLast = azi_j
            if (model.omegaControl)
                if (model.usingRotorSpeedFunction)
                    azi_j,Omega_j,OmegaDot_j = getRotorPosSpeedAccelAtTime(t[i],t[i]+delta_t,azi_s,delta_t)
                else
                    Omega_j = Omega_s
                    OmegaDot_j = OmegaDot_s
                    azi_j = azi_s + Omega_j*delta_t*2*pi
                end
            end
            if (!model.omegaControl)
                Crotor = 0
                Krotor = 0
                azi_j,Omega_j,OmegaDot_j = updateRotorRotation(structureMOI[3,3],Crotor,Krotor,
                -FReaction_j[6],-torqueDriveShaft_j,
                azi_s,Omega_s,OmegaDot_s,delta_t)
            end
            #-------------------------------------        ##

            ## evaluate aerodynamic module (CACTUS ONE-WAY)
            ##------== aerodynamics module ------------------------==
            if (model.aeroLoadsOn) # TODO: this is odd since we load in the aero loads elsewhere
                #             println('using cactus aero loads...')
                #             #----- calculate aerodynamic loads --------------------------------
                #             # this is place holder
                #             [FAero] = getAeroLoads(model.bladeData,model.aeroloadfile,t[i],el.props,model.totalNumDof)
                #             # these loads will need to account for a turbine with a
                #             # platform with different orientation than the inertial system.
                #             FAeroDof = (1:length(u_s))
                #             #------------------------------------------------------------------
            else
                FAero = []
                FAeroDof = []
            end
            #-------------------------------------       ## evaluate aerodynamic module (TU DELFT)
            ##------== aerodynamics module ------------------------==
            if (aeroOn && needsAeroCalcAtThisTimestep)

                #             [FAero,FAeroDof] = aeroModule(model,t[i] + delta_t,u_j,Omega_j,azi_j,numDOFPerNode,d_input_streamAero,d_output_streamAero)
                #
                #             #set aero forces flag
                #             needsAeroCalcAtThisTimestep = false

            end
            #-------------------------------------       ##

            ## compile external forcing on rotor
            #compile forces to supply to structural dynamics solver
            # Fexternal_sub, Fdof_sub = externalForcing(t[i]+delta_t,aeroLoads)

            step_AC = ceil(azi_j/(2*pi/ntheta)) #current_rot_angle/angle_per_step = current step (rounded since the structural ntheta is several thousand, whereas the actuator cylinder really can't handle that many)
            t_used = t[i]
            # function runme()
            #     println("here")
            #     println("here")
            #     println("here")
                # println("here")
            mat"[$Fexternal_sub, $Fdof_sub, $env2] = mapACloads($u_j,$udot_j,$Omega_j,$t_used,$PEy,$QCy,$NElem,$NBlade,$RefR,$mesh,$structuralSpanLocNorm,$structuralNodeNumbers,$structuralElNumbers,$el,$turbine3D,$env,$step_AC)"
            # end
            # Juno.@enter runme()

            if isempty(FAeroDof)
                Fdof = Fdof_sub
                Fexternal = Fexternal_sub
            else
                Fdof = [Fdof_sub, FAeroDof]
                Fexternal = [Fexternal_sub; FAero]
            end


            ## evaluate structural dynamics
            #call structural dynamics solver
            #initialization of structural dynamics displacements, velocities, accelerations, etc.
            #             dispData.displ_sm1 = u_sm1 #Not even used
            dispData = DispData(u_s,udot_s,uddot_s)

            #         if (occursin('ROM',model.analysisType))
            #             dispData.displ_s = u_s
            #             dispData.displdot_s = udot_s
            #             dispData.displddot_s = uddot_s
            #
            #             dispData.eta_s     = eta_s
            #             dispData.etadot_s  = etadot_s
            #             dispData.etaddot_s = etaddot_s
            #         end

            if (occursin("ROM",model.analysisType))
                error("ROM not fully implemented")
                #             # evalulate structural dynamics using reduced order model
                #             [dispOut,FReaction_j] = structuralDynamicsTransientROM(model,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData)
            else
                # evalulate structural dynamics using conventional representation
                t_in = t[i]
                # mat"[$elStrain,$dispOut,$FReaction_j] = structuralDynamicsTransient($model,$mesh,$el,$dispData,$Omega_j,$OmegaDot_j,$t_in,$delta_t,$elStorage,$Fexternal,$Fdof,$CN2H,$rbData)"
                # Juno.@enter call_structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t_in,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)
                elStrain,dispOut,FReaction_j = call_structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t_in,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData) #TODO: figure out how to pass structures
            end
            #update last iteration displacement vector
            u_jLast = u_j
            u_j = dispOut.displ_sp1             #update current estimates of velocity, acceleration
            # Only used for TNB, but must be declared
            udot_j  = dispOut.displdot_sp1
            uddot_j = dispOut.displddot_sp1


            if (occursin("ROM",model.analysisType))
                #             udot_j  = dispOut.displdot_sp1
                #             uddot_j = dispOut.displddot_sp1
                #
                #             eta_j = dispOut.eta_sp1
                #             etadot_j = dispOut.etadot_sp1
                #             etaddot_j = dispOut.etaddot_sp1
                error("ROM not fully implemented")
            end
            ##

            ## calculate norms
            uNorm = LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
            aziNorm = LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm

            if (model.hydroOn)
                platNorm = LinearAlgebra.norm(Ywec_j-Ywec_jLast)/LinearAlgebra.norm(Ywec_j) #platform module states iteration norm
            else
                platNorm = 0.0
            end

            if (model.driveTrainOn)
                gbNorm = LinearAlgebra.norm(gb_j - gb_jLast)/LinearAlgebra.norm(gb_j) #gearbox states iteration norm
            else
                gbNorm = 0.0
            end

            if (moduleIteration == false)
                break
            end

            numIterations = numIterations + 1
        end #end iteration while loop

        ## calculate converged generator torque/power
        if (model.useGeneratorFunction)
            if (model.generatorOn || (model.turbineStartup==0))
                println("simpleGenerator not fully implemented")#[genTorquePlot] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio)
                genTorquePlot = 0
            else
                genTorquePlot = 0
            end
        else
            genTorquePlot = 0
        end
        genPowerPlot = genTorquePlot*(gbDot_j*2*pi)*model.gearRatio


        ## update timestepping variables and other states, store in history arrays
        if (occursin("TD",model.analysisType))
            u_sm1 = u_s
            u_s = u_j
        end
        if (occursin("TNB",model.analysisType))
            u_s = u_j
            udot_s = udot_j
            uddot_s = uddot_j
        end
        if (occursin("ROM",model.analysisType))
            #         u_s = u_j
            #         udot_s = udot_j
            #         uddot_s = uddot_j
            #
            #         eta_s = eta_j
            #         etadot_s = etadot_j
            #         etaddot_s = etaddot_j
            error("ROM not fully implemented")
        end

        uHist[:,i+1] = u_s
        FReactionHist[i+1,:] = FReaction_j
        for ii = 1:length(elStrain)
            strainHist[ii,i] = ElStrain(elStrain[ii].eps_xx_0,
                elStrain[ii].eps_xx_z,
                elStrain[ii].eps_xx_y,
                elStrain[ii].gam_xz_0,
                elStrain[ii].gam_xz_y,
                elStrain[ii].gam_xy_0,
                elStrain[ii].gam_xy_z)
        end
        t[i+1] = t[i] + delta_t

        azi_s = azi_j
        Omega_s = Omega_j
        OmegaDot_s = OmegaDot_j

        genTorque_s = genTorque_j
        torqueDriveShaft_s = torqueDriveShaft_j

        aziHist[i+1] = azi_s
        OmegaHist[i+1] = Omega_s
        OmegaDotHist[i+1] = OmegaDot_s

        gb_s = gb_j
        gbDot_s = gbDot_j
        gbDotDot_s = gbDotDot_j

        gbHist[i+1] = gb_s
        gbDotHist[i+1] = gbDot_s
        gbDotDotHist[i+1] = gbDotDot_s

        #genTorque[i+1] = genTorque_s
        genTorque[i+1] = genTorquePlot
        genPower[i+1] = genPowerPlot
        torqueDriveShaft[i+1] = torqueDriveShaft_s

        if (model.hydroOn)
            error("Hydro Model not fully implemented")
            #         Ywec(i+1,:) = Ywec_j
            #         rigidDof(i+1,:)=Ywec(i+1,s_stsp+1:s_stsp+6)

        else
            rigidDof[i] = 0
        end

        FReactionHist[i+1,:] = FReaction_j
        ##

        ## check rotor speed for generator operation
        if (Omega_s>= rotorSpeedForGenStart)
            model.generatorOn = true
        else
            model.generatorOn = false
        end
        ##

    end #end timestep loop

    ## kill platform module process
    if (model.hydroOn)
        #     serverSendVector(-1.0,d_output_streamPlatform)
        #     terminateServer(server_socketPlatform,output_socketPlatform,1)
        #     terminateClient(input_socketPlatform,1)
    end

    ## kill aerodynamic module process
    if (aeroOn)
        #     serverSendVector(decodeVec(4,-1.0),d_output_streamAero)
        #     terminateClient(input_socketAero,1) #close down client connection to forcing module
        #     terminateServer(server_socketAero,output_socketAero,1) #close down server on this side
    end

    ##
    #toc
    # save aeroOutputArray
    #save simulation data in .mat file
    # save(model.outFilename,'t','uHist','aziHist','OmegaHist','OmegaDotHist','gbHist','gbDotHist','gbDotDotHist','FReactionHist','rigidDof','genTorque','genPower','torqueDriveShaft','strainHist')
    # fprintf('#s\n','Output Saving Not Currently Implemented')

    #Writefile
    writefile = true
    if !writefile
        println("NOT WRITING Verification File")
    else
        println("WRITING Verification File")
        mat"write_verification($model,$t,$aziHist,$OmegaHist,$OmegaDotHist,$gbHist,$gbDotHist,$gbDotDotHist,$FReactionHist,$rigidDof,$genTorque,$genPower,$torqueDriveShaft,$uHist,$strainHist)"
    end

end

function omegaSpecCheck(tCurrent,tocp,Omegaocp,delta_t)

    if (tocp[length(tocp)]<tCurrent)
        println("Simulation time is greater than that specified in control points for prescribed rotor speed profile.")
        println("Terminating simulation.")
        terminateSimulation = true
        OmegaCurrent = 0.0
        OmegaDotCurrent = 0.0
    else
        OmegaCurrent = FLOWMath.linear(tocp,Omegaocp,tCurrent) #interpolated discreteized profile for current omega

        #calculate current rotor acceleration
        dt = delta_t/2.0
        omega_m1 = FLOWMath.linear(tocp,Omegaocp,tCurrent-dt)
        omega_p1 = FLOWMath.linear(tocp,Omegaocp,tCurrent+dt)

        OmegaDotCurrent = diff([omega_m1,omega_p1])/(dt*2)
        OmegaDotCurrent = OmegaDotCurrent[1]

        terminateSimulation = false
        if (isnan(OmegaCurrent) || isnan(OmegaDotCurrent))
            error("Omega calcualted a NaN. Exiting.")
        end
    end
    return OmegaCurrent,OmegaDotCurrent,terminateSimulation
end
