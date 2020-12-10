
include("call_structuralDynamicsTransient.jl")

#These are from the matlab version of the actuator cylinder
# mutable struct Slice
#     r
#     RefR
#     chord
#     twist
#     delta
#     B
#     alpha_rad
#     cl
#     cd
#     Omega
#     centerX
#     centerY
# end
#
# mutable struct Env
#     rho
#     mu
#     G_amp
#     gusttime
#     gustX0
#     N_Rev
#     idx_sub
#     wsave
#     Vinf_nominal
#     V_vert
#     V_tang
#     V_rad
#     V_twist
#     Vinf
#     V_wake_old
#     steplast
#     ntheta
#     tau
# end

function Unsteady(model,mesh,el)
    #Unsteady performs modular transient analysis
    #
    #   Unsteady(model,mesh,el)
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
    aeroOn = false #TODO: clean this up
    CACTUS = true #TODO: not hardcoded
    println("CACTUS AERO: $CACTUS")
    #modularIteration
    moduleIteration = true
    aeroLoadsFile_root = model.aeroloadfile[1:end-16] #cut off the _ElementData.csv
    OWENSfile_root = model.owensfile[1:end-6] #cut off the .owens
    # Get AeroLoads
    if CACTUS
        # mat"$aeroLoads = processAeroLoadsBLE($model.aeroloadfile, $model.owensfile)"


        d1 = string(aeroLoadsFile_root, ".geom")
        d2 = string(aeroLoadsFile_root, "_ElementData.csv")
        d3 = string(OWENSfile_root, ".bld")
        d4 = string(OWENSfile_root, ".el")
        d5 = string(OWENSfile_root, ".ort")
        d6 = string(OWENSfile_root, ".mesh")

        # mat"[$aerotimeArray1,$aeroForceValHist1,$aeroForceDof1] = mapCactusLoadsFile($d1,$d2,$d3,$d4,$d5,$d6)"
        aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = mapCactusLoadsFile(d1,d2,d3,d4,d5,d6)

    else
        #TODO: not hard coded
        QCx=[0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00]
        QCy=[0.00000e+00,2.44680e-01,4.89360e-01,7.34040e-01,9.78720e-01,1.22340e+00,1.46808e+00,1.71276e+00,1.95744e+00,2.20212e+00,2.44680e+00]
        QCz=[0.00000e+00,3.76945e-01,6.77076e-01,8.73458e-01,9.77684e-01,1.00000e+00,9.44658e-01,8.13441e-01,6.09216e-01,3.25362e-01,0.00000e+00]
        CtoR=[1.11093e-01,1.11093e-01,8.20390e-02,6.35940e-02,5.68145e-02,5.55467e-02,5.88008e-02,6.82860e-02,9.11773e-02,1.11093e-01,1.11093e-01]
        RefR=177.2022
        NBlade = 2
        PEy=[1.22340e-01,3.67020e-01,6.11700e-01,8.56380e-01,1.10106e+00,1.34574e+00,1.59042e+00,1.83510e+00,2.07978e+00,2.32446e+00]
        NElem = 10

        bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers = readBladeData(string(OWENSfile_root, ".bld"))
        # mat"[$structuralSpanLocNorm,$structuralNodeNumbers,$structuralElNumbers] = readBldFile($OWENSfile_root)"

        step_AC = 0
        ntheta = 36#1/(model.OmegaInit*model.delta_t)
        Vinf = 25.0 #TODO
        rho = 1.225 #TODO
        mu = 1.7894e-5 #TODO

        H = [] #For plotting turbine
        plotTurbine = false

        # mat"[$alpha_rad,$cl_af,$cd_af] = readaerodyn('airfoils/NACA_0015_RE3E5.dat')"
        path,_ = splitdir(@__FILE__)
        af = VAWTAero.readaerodyn("$path/../test/airfoils/NACA_0015_RE3E5.dat") #TODO: make this path smarter
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
        # slice = Slice(ones(ntheta),
        #     RefR,
        #     0.0,
        #     zeros(1,ntheta),
        #     zeros(1,ntheta),
        #     Float64(NBlade),
        #     alpha_rad,
        #     cl,
        #     cd,
        #     ones(1,ntheta)*model.OmegaInit*2*pi,
        #     0.0,
        #     0.0)

        slice = VAWTAero.Turbine(RefR,
        zeros(ntheta),
        zeros(1),
        zeros(ntheta),
        zeros(ntheta),
        zeros(ntheta),
        2,
        1.0,
        af,
        ntheta,
        false,
        0.0,
        0.0)

        # turbine built from bottom up
        turbine3D = fill(slice, n_slices)
        for i = 1:n_slices
            turbine3D[i].r[:] .= ones(ntheta)*r[i]
            turbine3D[i].chord[:] .= chord[i]
            turbine3D[i].twist[:] .= ones(ntheta)*twist[i]
            turbine3D[i].delta[:] .= ones(ntheta)*delta[i]
        end

        # env1 = Env(rho,
        #     mu,
        #     0.0, #m/s gust
        #     0.8, #sec
        #     13.0, #radaii offset back is starting point for gust
        #     20,
        #     zeros(1,NBlade*2),
        #     zeros(1,ntheta*2),
        #     Vinf,
        #     zeros(1,ntheta),
        #     zeros(1,ntheta),
        #     zeros(1,ntheta),
        #     zeros(1,ntheta),
        #     ones(ntheta)*Vinf,
        #     Vinf,
        #     0,
        #     ntheta,
        #     [0.3025, 2.9500])

        env1 = VAWTAero.Environment(rho,
        mu,
        ones(ntheta).*Vinf,
        "AC", #TODO: use this to unify the DMS and AC steady methods, like in the unsteady method
        "None",
        false, #TODO: generalize this name since it is now being used as an option on the solution method for both AC and DMS, not just DMS anymore
        Vinf,
        zeros(ntheta))

        env = fill(env1, n_slices)

        if env1.AModel == "AC"
            N_aw = ntheta*2 #Number of either "a" (DMS induction factor) or w (u and v induction factors)
        elseif env1.AModel == "DMS"
            N_aw = ntheta
        else
            error("Aeromodel not recognized, choose AC or DMS")
        end

        us_param = VAWTAero.UnsteadyParams(true,
        0,
        0.0,
        0.0,
        0.0,
        zeros(Int,1),
        0.0,
        [0.3,3.0],
        zeros(Int,2*2),
        zeros(N_aw),
        zeros(1))
    end


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
    # strainHist = fill(ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)),mesh.numEl,numTS)
    eps_xx_0_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_z_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_z_hist = zeros(4,mesh.numEl,numTS)

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
        # mat"$elStorage = initialElementCalculations($model,$el,$mesh)" #perform initial element calculations for conventional structural dynamics analysis
        elStorage = initialElementCalculations(model,el,mesh) #perform initial element calculations for conventional structural dynamics analysis
    end

    #calculate structural/platform moi
    _,structureMOI,_=calculateStructureMassProps(elStorage)
    #..........................................................................

    ## Main Loop - iterate for a solution at each time step, i
    for i=1:numTS

        #     i #TODO add verbose printing
        # if (mod(i,100)==0) #print command that displays progress of time stepping
        println("Time Step: $i")
        # end

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

        TOL = 1e-5  #gauss-seidel iteration tolerance for various modules
        MAXITER = 50 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        platNorm = 1e5
        aziNorm = 1e5
        gbNorm = 1e5 #initialize norms for various module states
        ##

        while ((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER)) #module gauss-seidel iteration loop
            println("$(numIterations)   uNorm: $(uNorm)    platNorm: $(platNorm)    aziNorm: $(aziNorm)    gbNorm: $(gbNorm)")
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

            # Freq Domain of Hydro does give stiffness and damping
            # Should be consistent with how OrcaFlex/Hydrodyn would work for modularity

            # Possible best way: Kevin run simple case to produce example reaction forces at bases for time history,
            # then Ryan writes code that calculates new motion and sitffness/damping

            # Other way: OWENS solves for motions and hydro solves for reaction forces

            # Kevin look at wavec2wire (dissertation might help)
            # rigid body motion should not be calculated within hydro

            # Near term coupling is easier in time domain

            # Assignments
            # - Kevin look at wavec2wire, OrcaFlex, hydrodyn, mass and stiffness approx
            # - Kevin give Ryan example of mass and stiffness approx how it's being used/read in
            # - Owens gives motions, hydro gives mass/stiffness, owens recalculates motions, and reiterate

            # Hydro Coupling
            #1) modify the mesh to include a platform cg offset #option
            #2) modify the element properties so that this element is rigid
            #3) ensure that this element is non-rotating - deflected solution is solved in one step (perhaps just give rotated stiffenss, forces, etc, so the node is actually rotating, as long as the mass isn't being cyntrifical)
            #4) apply the hydro F, C, K, M to the nodal term (we have calm sea, as well as waves and currents). Mooring will also apply forces.
            #5) use flags to turn degrees of freedom off

            # Unit Testing
            # 1) Turn all dof off and the solution should be the same
            # 2) allow heave and the turbine should go up and down
            # 3) no aero, but mass offset should cause the turbine to tilt slightly

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
            println("starting Aero Loads")
            if CACTUS
                Fexternal_sub, Fdof_sub = externalForcing(t[i]+delta_t,aerotimeArray,aeroForceValHist,aeroForceDof)
            else
                step_AC = ceil(azi_j/(2*pi/ntheta)) #current_rot_angle/angle_per_step = current step (rounded since the structural ntheta is several thousand, whereas the actuator cylinder really can't handle that many)
                t_used = t[i]
                # function runme()
                #     println("here")
                #     println("here")
                #     println("here")
                # println("here")
                # Juno.@enter mapACloads(u_j,udot_j,Omega_j,t_used,PEy,QCy,NElem,NBlade,RefR,mesh,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,el,turbine3D,env,step_AC,us_param)

                Fexternal_sub, Fdof_sub, env = mapACloads(u_j,udot_j,Omega_j,t_used,PEy,QCy,NElem,NBlade,RefR,mesh,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,el,turbine3D,env,step_AC,us_param)

                # end
                # Juno.@enter runme()
            end

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
            println("starting struct dyn")
            if (occursin("ROM",model.analysisType))
                error("ROM not fully implemented")
                #             # evalulate structural dynamics using reduced order model
                #             [dispOut,FReaction_j] = structuralDynamicsTransientROM(model,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,rom,Fexternal,Fdof,CN2H,rbData)
            else
                # evalulate structural dynamics using conventional representation
                t_in = t[i]
                # mat"[$elStrain,$dispOut,$FReaction_j] = structuralDynamicsTransient($model,$mesh,$el,$dispData,$Omega_j,$OmegaDot_j,$t_in,$delta_t,$elStorage,$Fexternal,$Fdof,$CN2H,$rbData)"
                # Juno.@enter call_structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t_in,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)
                start = time()
                elStrain,dispOut,FReaction_j = call_structuralDynamicsTransient(model,mesh,el,dispData,Omega_j,OmegaDot_j,t_in,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData) #TODO: figure out how to pass structures
                println("$(time()-start)")
            end
            #update last iteration displacement vector
            u_jLast = u_j
            u_j = dispOut.displ_sp1             #update current estimates of velocity, acceleration
            # Only used for TNB, but must be declared
            udot_j  = dispOut.displdot_sp1
            uddot_j = dispOut.displddot_sp1

            println("ran structdyn")
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
            if numIterations==MAXITER
                @warn "Maximum Iterations Met"
            end
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
            eps_xx_0_hist[:,ii,i] = elStrain[ii].eps_xx_0
            eps_xx_z_hist[:,ii,i] = elStrain[ii].eps_xx_z
            eps_xx_y_hist[:,ii,i] = elStrain[ii].eps_xx_y
            gam_xz_0_hist[:,ii,i] = elStrain[ii].gam_xz_0
            gam_xz_y_hist[:,ii,i] = elStrain[ii].gam_xz_y
            gam_xy_0_hist[:,ii,i] = elStrain[ii].gam_xy_0
            gam_xy_z_hist[:,ii,i] = elStrain[ii].gam_xy_z
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
        # mat"write_verification($model,$t,$aziHist,$OmegaHist,$OmegaDotHist,$gbHist,$gbDotHist,$gbDotDotHist,$FReactionHist,$rigidDof,$genTorque,$genPower,$torqueDriveShaft,$uHist,$strainHist)"

        filename = string(model.outFilename[1:end-3], "h5")
        HDF5.h5open(filename, "w") do file
            # HDF5.write(file,"model",model)
            HDF5.write(file,"t",t)
            HDF5.write(file,"aziHist",aziHist)
            HDF5.write(file,"OmegaHist",OmegaHist)
            HDF5.write(file,"OmegaDotHist",OmegaDotHist)
            HDF5.write(file,"gbHist",gbHist)
            HDF5.write(file,"gbDotHist",gbDotHist)
            HDF5.write(file,"gbDotDotHist",gbDotDotHist)
            HDF5.write(file,"FReactionHist",FReactionHist)
            HDF5.write(file,"rigidDof",rigidDof)
            HDF5.write(file,"genTorque",genTorque)
            HDF5.write(file,"genPower",genPower)
            HDF5.write(file,"torqueDriveShaft",torqueDriveShaft)
            HDF5.write(file,"uHist",uHist)
            HDF5.write(file,"eps_xx_0_hist",eps_xx_0_hist)
            HDF5.write(file,"eps_xx_z_hist",eps_xx_z_hist)
            HDF5.write(file,"eps_xx_y_hist",eps_xx_y_hist)
            HDF5.write(file,"gam_xz_0_hist",gam_xz_0_hist)
            HDF5.write(file,"gam_xz_y_hist",gam_xz_y_hist)
            HDF5.write(file,"gam_xy_0_hist",gam_xy_0_hist)
            HDF5.write(file,"gam_xy_z_hist",gam_xy_z_hist)
        end

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

function userDefinedGenerator(input)
    return 0.0 #this is what was in the original file...
end

function getRotorPosSpeedAccelAtTime(t0,time,aziInit,delta_t)
    #getRotorPosSpeedAccelAtTime uses user defined function to get rotor pos.
    #   [rotorAzimuth,rotorSpeed,rotorAcceleration] = getRotorPosSpeedAccelAtTime(t0,time,aziInit)
    #
    #   This function uses the user defined function rotorSpeedProfile() to get
    #   the azimuth, speed, and acceleration of the rotor.
    #
    #   input:
    #   t0      = time at which azimuth integration is beginning
    #   time    = current time that position, velocity, and acceleration are
    #             being requested
    #   aziInit = initial rotor azimuth angle integration will begin at
    #
    #   output:
    #   rotorAzimuth = azimuth position of rotor (rad) at time
    #   rotorSpeed   = rotor speed (Hz) at time
    #   rotorAcceleration = rotor acceleration (Hz/s) at time


    rotorSpeed = userDefinedRotorSpeedProfile(time) #get rotgor speed at time

    dt = 0.01#some small delta t used in estimating rotor acceleration
    if ((time-dt) < 0)
        dt = delta_t/2
    end

    omega_p1 = userDefinedRotorSpeedProfile(time+dt) #get rotor speed slightly before and after time
    omega_m1 = userDefinedRotorSpeedProfile(time-dt)
    #--------------------------------------------------------------------------

    #estimate rotor acceleration with difference calculation
    rotorAcceleration = diff([omega_m1,omega_p1])/(2*dt)

    #calculate rotor azimuth using trapezoidal rule
    rotorAzimuth = trapezoidalRule(aziInit,userDefinedRotorSpeedProfile[t0],rotorSpeed,time-t0)

    return rotorAzimuth,rotorSpeed,rotorAcceleration
end

#simple trapezoidal rule integration
function trapezoidalRule(aziInit,rotorSpeedStart,rotorSpeedEnd,dt)
    return (aziInit + 0.5*dt*(rotorSpeedStart+rotorSpeedEnd)/(2*pi))
end

function userDefinedRotorSpeedProfile(time)
    return 0.5 #this is what was originally in the file...
end

function setInitialConditions(initCond,u,numDOFPerNode)
    #setInitialConditions sets initial conditions
    #   [u] =  setInitialConditions(initCond,u,numDOFPerNode)
    #
    #   This function reads initial conditions from file
    #
    #   input:
    #   initCond      = array containing initial conditions
    #                     initCond(i,1) = node number for init cond i
    #                     initCond(i,2) = local DOF number for init cond i
    #                     initCond(i,3) = value for init cond i
    #   u             = displacement vector
    #   numDOFPerNode = number of degrees of freedom per node
    #
    #   output:
    #    u             = displacement vector modified for initial conditions

    len=size(initCond) #get number of specified initial conditions
    #unspecified initial conditions are assumed to
    #be zero

    for i=1:len[1] #loop over initial conditions
        if (initCond[i,2]>numDOFPerNode) #error check
            error("setInitalConditios:: DOF greater than numDOFPerNode")
        end
        index = (initCond[i,1]-1)*numDOFPerNode + initCond[i,2] #calculate global DOF number for initial condition
        u[index] = initCond[i,3] #specify initial condition using global DOF number
    end

    return u

end

function externalForcing(time,timeArray,ForceValHist,ForceDof)

    #owens externalForcing function for the OWENS toolkit
    #   [Fexternal, Fdof] = externalForcing(time,aeroLoads)
    #
    #   This function specifies external forcing for a transient analysis.
    #   Fexternal is a vector of loads and Fdof is a corresponding vector of
    #   degrees of freedom the concentrated loads in Fexternal correspond to.
    #   The input time allows for arbitrary time varying loads
    #   The global degree of freedom number corresponding with the local degree
    #   of freedom of a node may be calculated by:
    #   globalDOFNumber = (nodeNumber-1)*6 + localDOFnumber
    #   The localDOFnumber may range from 1 to 6 such that 1 corresponds to a
    #   force in "x direction" of the co-rotating hub frame. 2 and 3
    #   corresponds to a force in the "y" and "z directions" respectively. 4,
    #   5, and 6 correspond to a moment about the "x", "y", and "z" directions
    #   respectively.

    #
    #      input:
    #      time         = simulation time
    #
    #      output:
    #      Fexternal     = vector of external loads (forces/moments)
    #      Fdof          = vector of corresponding DOF numbers to apply loads to


    #     if(time < 0.2)
    #         Fexternal = 1e6
    #         Fdof = 20*6+1
    #     else
    #         Fexternal = []
    #         Fdof = []
    #     end

    #temp = load('aeroLoads.mat')
    # timeArray = aeroLoads["timeArray"]
    # ForceValHist = aeroLoads["ForceValHist"]
    # ForceDof = aeroLoads["ForceDof"]
    Fexternal = zeros(length(ForceDof))

    for i = 1:length(ForceDof)
        Fexternal[i] = FLOWMath.linear(timeArray,ForceValHist[i,:],time)
    end
    Fdof = ForceDof

    return Fexternal, Fdof

end

function calculateDriveShaftReactionTorque(driveShaftProps,thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)
    #calculateDriveShaftReactionTorque calculates reaction torque of driveshaft
    #   [torque] = calculateDriveShaftReactionTorque(driveShaftProps,...
    #                thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)
    #
    #   This function calculates reaction torque of driveshaft
    #
    #   input:
    #   driveShaftProps      = object containing driveshaft properties
    #   thetaRotor           = azimuth position of rotor/rotor shaft (rad)
    #   thetaGB              = azimuth position of gearbox shaft (rad)
    #   thetaDotRotor        = angular velocity of rotor/rotor shaft (rad/s)
    #   thetaDotGB           = angular velocity of gearbox shaft (rad/s)
    #
    #   output:
    #   torque   = reaction torque of drive shaft


    k = driveShaftProps.k  #drive shaft stiffness
    c = driveShaftProps.c  #drive shaft damping

    return k*(thetaRotor-thetaGB) + c*(thetaDotRotor-thetaDotGB)

end

function updateRotorRotation(Irotor,Crotor,Krotor,
    shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s,
    delta_t)
    #updateRotorRotation updates rotor rotation
    #
    #   [azi_sp1,Omega_sp1,OmegaDot_sp1] = updateRotorRotation(Irotor,Crotor,Krotor,...
    #                                  shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s,...
    #                                  delta_t)
    #
    #   This function updates the rotor rotation given rotor properties and external
    #   torques
    #
    #   input:
    #   Irotor      = rotor inertia
    #   Crotor      = arbitrary rotor damping
    #   Krotor      = arbitrary rotor stiffness
    #   shaftTorque = torque from external forces on rotor
    #   genTorque   = torque from generator
    #   azi_s       = rotor azimuth (rad) at beginning of time step
    #   Omega_s     = rotor speed (Hz) at beginning of time step
    #   OmegaDot_s  = rotor acceleration (Hz/s) at beginning of time step
    #   delta_t     = time step
    #
    #   output:
    #   azi_sp1       = rotor azimuth (rad) at end of time step
    #   Omega_sp1     = rotor speed (Hz/s) at end of time step
    #   OmegaDot_sp1  = rotor acceleration (Hz/s) at end of time step
    #
    Frotor = shaftTorque + genTorque #calculate effective torque on rotor
    Omega_s = Omega_s*2*pi #conversion form Hz to rad/s, etc.
    OmegaDot_s = OmegaDot_s*2*pi
    azi_sp1,Omega_sp1,OmegaDot_sp1 = timeIntegrateSubSystem(Irotor,Krotor,Crotor,Frotor, #time integrate using Newmark-Beta
    delta_t,azi_s,Omega_s,OmegaDot_s)

    Omega_sp1 = Omega_sp1/(2*pi) #convert to Hz, etc.
    OmegaDot_sp1 = OmegaDot_sp1/(2*pi)

    return azi_sp1,Omega_sp1,OmegaDot_sp1

end


function timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)
    #timeIntegrateSubSystem integrates a system using Newmark-Beta method
    #
    #   [unp1,udotnp1,uddotnp1] = timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)
    #
    #   #This function perform integration of a system using the Newmark-Beta
    #   method(constant-average acceleration sceheme).
    #
    #   input:
    #   M        = system mass matrix
    #   K        = system sttiffness matrix
    #   C        = system damping matrix
    #   F        = system force vector
    #   delta_t  = time step
    #   u        = displacement at beginning of time step
    #   udot     = velocity at beginning of time step
    #   uddot    = acceleration at beginning of time step
    #
    #
    #   output:
    #   unp1        = displacement at end of time step
    #   udotnp1     = velocity at end of time step
    #   uddotnp1    = acceleration at end of time step
    #
    alpha = 0.5 #constant avg accel scheme
    gamma = 0.5
    beta = 0.5*gamma

    a1 = alpha*delta_t
    a2 = (1.0-alpha)*delta_t
    a3 = 1.0/(beta*delta_t*delta_t)
    a4 = a3*delta_t
    a5 = 1.0/gamma-1.0
    a6 = alpha/(beta*delta_t)
    a7 = alpha/beta - 1.0
    a8 = delta_t*(alpha/gamma-1.0)

    A = a3*u + a4*udot + a5*uddot
    B = a6*u + a7*udot + a8*uddot

    Khat = K + a3.*M + a6.*C
    Fhat = F + M*(A') + C*(B')

    unp1 = Khat\Fhat

    uddotnp1 = a3*(unp1-u) - a4*udot - a5*uddot
    udotnp1 =  udot + a2*uddot + a1*uddotnp1

    return unp1,udotnp1,uddotnp1

end
