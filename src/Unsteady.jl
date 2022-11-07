"""

Unsteady(model,topModel,mesh,el,aero;getLinearizedMatrices=false)

Executable function for transient analysis. Provides the interface of various
    external module with transient structural dynamics analysis capability.

    # Input
    * `inputs::Model`: see ?Model
    * `topModel::FEAModel`: see ?GyricFEA.FEAModel
    * `mesh::Mesh`: see ?GyricFEA.Mesh
    * `el::El`: see ?GyricFEA.El
    * `bin::Bin`: see ?Bin
    * `aero::function`: Fexternal, Fdof = aero(t) where Fexternal is the force on each affected mesh dof and Fdof is the corresponding DOFs affected
    * `getLinearizedMatrices::Bool`: Flag to save the linearized matrices
    * `elStorage::ElStorage.ElStorage`: Optional object containing stored element matrices
    * `u_s::Array{<:float}`: Optional warm start of top deflections, of length Nnodes x Ndof


    # Output
    * `t`: time array
    * `aziHist`: azimuthal history array
    * `OmegaHist`: rotational speed array history
    * `OmegaDotHist`: rotational acceleration array history
    * `gbHist`: gearbox position history array
    * `gbDotHist`: gearbox velocity history array
    * `gbDotDotHist`: gearbox acceleration history array
    * `FReactionHist`: Base reaction 6dof forces history
    * `rigidDof`:
    * `genTorque`: generator torque history
    * `genPower`: generator power history
    * `torqueDriveShaft`: driveshaft torque history
    * `uHist`: mesh displacement history for each dof
    * `epsilon_x_hist`: strain history for eps_xx_0 for each dof
    * `epsilon_y_hist`: strain history for eps_xx_z for each dof
    * `epsilon_z_hist`: strain history for eps_xx_y for each dof
    * `kappa_x_hist`: strain history for gam_xz_0 for each dof
    * `kappa_y_hist`: strain history for gam_xz_y for each dof
    * `kappa_z_hist`: strain history for gam_xy_0 for each dof
    """
function Unsteady(inputs;topModel=nothing,topMesh=nothing,topEl=nothing,
    aeroVals=nothing,aeroDOFs=nothing,aero=nothing,deformAero=nothing,
    bottomModel=nothing,bottomMesh=nothing,bottomEl=nothing,bin=nothing,
    getLinearizedMatrices=false,
    system=nothing,assembly=nothing, #TODO: should we initialize them in here? Unify the interface for ease?
    topElStorage = nothing,bottomElStorage = nothing, u_s = nothing)

    #..........................................................................
    #                             INITIALIZATION
    #..........................................................................

    if (!inputs.topsideOn) && (!inputs.hydroOn)
        error("No structure is being simulated!")
    end

    ## General
    delta_t = inputs.delta_t
    numTS = Int(inputs.numTS)
    numDOFPerNode = 6
    CN2H = LinearAlgebra.I(3) # hub and inertial frames initialize as copies
    t = range(0, length=numTS, step=delta_t)
    integrator = 0.0 #for generator control algorithm
    integrator_j = 0.0
    topDispOut = [] #TODO: better way to control scope

    # g = [0.0, 0.0, -9.80665]
    uHist,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FReactionHist,FTwrBsHist,aziHist,OmegaHist,OmegaDotHist,gbHist,
    gbDotHist,gbDotDotHist,genTorque,genPower,torqueDriveShaft,uHist_prp,
    FPtfmHist,FHydroHist,FMooringHist = allocate_general(inputs,topMesh,numDOFPerNode,numTS,assembly)

    if inputs.topsideOn
        # Allocate memory for topside
        u_s,udot_s,uddot_s,u_sm1,topDispData1,topDispData2,topElStrain,gb_s,
        gbDot_s,gbDotDot_s,azi_s,Omega_s,OmegaDot_s,genTorque_s,
        torqueDriveShaft_s,topFexternal,topFexternal_hist = allocate_topside(inputs,topMesh,topEl,topModel,numDOFPerNode,u_s)
    end

    ## Hydrodynamics/mooring module initialization and coupling variables
    if inputs.hydroOn
        bottom_totalNumDOF,u_s_ptfm_n,udot_s_ptfm_n,uddot_s_ptfm_n,u_sm1_ptfm_n,bottomDispData,prpDOFs,u_s_prp_n,
        udot_s_prp_n,uddot_s_prp_n,jac,numMooringLines,FHydro_n,FMooring_n,outVals,
        mooringTensions = allocate_bottom(t,numTS,delta_t,inputs,bottomMesh,bottomEl,bottomModel,bin,numDOFPerNode)
    end

    ## Rotor mode initialization
    rotorSpeedForGenStart = initialize_generator!(inputs)

    ## Structural dynamics initialization
    if isnothing(topElStorage) && inputs.topsideOn
        topElStorage = GyricFEA.initialElementCalculations(topModel,topEl,topMesh) #perform initial element calculations for conventional structural dynamics analysis
    end
    if isnothing(bottomElStorage) && inputs.hydroOn
        bottomElStorage = GyricFEA.initialElementCalculations(bottomModel,bottomEl,bottomMesh) #perform initial element calculations for conventional structural dynamics analysis
    end
    if inputs.analysisType=="ROM"
        if inputs.topsideOn
            top_rom,topJointTransformTrans,u_sRed,udot_sRed,uddot_sRed,
            topBC,u_s2,udot_s2,uddot_s2,top_invPhi,eta_s,etadot_s,
            etaddot_s = initialize_ROM(topElStorage,topModel,topMesh,topEl,u_s,udot_s,uddot_s)

            topDispData1.eta_s = eta_s
            topDispData1.etadot_s = etadot_s
            topDispData1.etaddot_s = etaddot_s
            topDispData2.eta_s = eta_s
            topDispData2.etadot_s = etadot_s
            topDispData2.etaddot_s = etaddot_s
        end

        if inputs.hydroOn

            bottom_rom,bottomJointTransformTrans,u_sRed_ptfm_n,
            udot_sRed_ptfm_n,uddot_sRed_ptfm_n,bottomBC,u_s2_ptfm_n,udot_s2_ptfm_n,
            uddot_s2_ptfm_n,bottom_invPhi,eta_s_ptfm_n,etadot_s_ptfm_n,
            etaddot_s_ptfm_n = initialize_ROM(bottomElStorage,bottomModel,bottomMesh,bottomEl,u_s_ptfm_n,udot_s_ptfm_n,uddot_s_ptfm_n)

            bottomDispData.eta_s = eta_s_ptfm_n
            bottomDispData.etadot_s = etadot_s_ptfm_n
            bottomDispData.etaddot_s = etaddot_s_ptfm_n
        end
    end

    if inputs.topsideOn
        topsideMass, topsideMOI, topsideCG = GyricFEA.calculateStructureMassProps(topElStorage)
        topModel.jointTransform, topModel.reducedDOFList = GyricFEA.createJointTransform(topModel.joint,topMesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints

        if inputs.hydroOn
            hydro_topside_nodal_coupling!(bottomModel,bottomMesh,topsideMass,topModel,topsideCG,topsideMOI,numDOFPerNode)
        end # if inputs.hydroOn
    end # if inputs.topsideOn

    if inputs.topsideOn
        uHist[1,:] = u_s          #store initial condition
        aziHist[1] = azi_s
        OmegaHist[1] = Omega_s
        OmegaDotHist[1] = OmegaDot_s
        FReactionsm1 = zeros(Real, 6)
        FReactionHist[1,:] = FReactionsm1
        topFReaction_j = FReactionsm1
        # topWeight = [0.0, 0.0, topsideMass*-9.80665, 0.0, 0.0, 0.0] #TODO: propogate gravity, or remove since this isn't used
        gbHist[1] = gb_s
        gbDotHist[1] = gbDot_s
        gbDotDotHist[1] = gbDotDot_s
        genTorque[1] = genTorque_s
        torqueDriveShaft[1] = torqueDriveShaft_s
    end

    if inputs.hydroOn
        uHist_prp[1,:] = u_s_prp_n
    end

    #..................................................................
    #                          INITIAL SOLVE
    #..................................................................

    ## Evaluate mooring and hydrodynamics at t=0 based on initial conditions
    if inputs.hydroOn
    # function initial_solve_hydro(inputs,bottom_totalNumDOF,numDOFPerNode,FHydro_n,FMooring_n,t,delta_t)
        # Initial coupled bottomside solve using reaction force from topside
        bottomFexternal = zeros(6)
        bottomFDOFs = collect(bottom_totalNumDOF-numDOFPerNode+1:bottom_totalNumDOF)
        FPtfm_n = FHydro_n + FMooring_n

        if inputs.analysisType=="ROM"
            _, bottomCoupledDisps, FHydro_n, FMooring_n, outVals, jac = OWENS_HD_Coupled_Solve(t[1], delta_t, true, jac, numDOFPerNode, prpDOFs, FPtfm_n,
            bottomDispData, bottomModel, bottomMesh, bottomEl, bottomElStorage,
            bottomFexternal, bottomFDOFs, CN2H, bottom_rom)
        else
            _, bottomCoupledDisps, FHydro_n, FMooring_n, outVals, jac = OWENS_HD_Coupled_Solve(t[1], delta_t, true, jac, numDOFPerNode, prpDOFs, FPtfm_n,
            bottomDispData, bottomModel, bottomMesh, bottomEl, bottomElStorage,
            bottomFexternal, bottomFDOFs, CN2H)
        end

        # Update DispData with the new accelerations ONLY to account for added mass from HydroDyn
        # uddot_s_ptfm_h = bottomCoupledDisps.displddot_sp1
        uddot_s_ptfm_n = bottomCoupledDisps.displddot_sp1
        uddot_s_prp_n = uddot_s_ptfm_n[prpDOFs]
        # uddot_s_prp_n = frame_convert(uddot_s_ptfm_h[prpDOFs], CH2N)
        # u_sp1_prp_predState = frame_convert(bottomCoupledDisps.displ_sp1[prpDOFs], CH2N) # this is kind of like ED%x in FAST, which is advanced using an ABM4 method.
        # Since we don't have a 4th order integration method, we do this as a stopgap before motions are actually updated in the first time marching structural solve.
        # TODO: is this even needed? is the normal u_sp1_ptfm prediction accurate enough?
        #       what if we're including correction steps?
        u_sp1_prp_predState = bottomCoupledDisps.displ_sp1[prpDOFs]
        # bottomDispData = GyricFEA.DispData(u_s_ptfm_h, udot_s_ptfm_h, uddot_s_ptfm_h, u_sm1_ptfm_h)
        bottomDispData = GyricFEA.DispData(u_s_ptfm_n, udot_s_ptfm_n, uddot_s_ptfm_n, u_sm1_ptfm_n)
        FPtfm_n = FHydro_n + FMooring_n

        FPtfmHist[1,:] = FPtfm_n
        FHydroHist[1,:] = FHydro_n
        FMooringHist[1,:] = FMooring_n

        ## Copy values of the inital outputs to arrays for interpolation/extrapolation
        recent_u_prp = repeat(u_s_prp_n, 1, inputs.interpOrder+1)
        recent_udot_prp = repeat(udot_s_prp_n, 1, inputs.interpOrder+1)
        recent_uddot_prp = repeat(uddot_s_prp_n, 1, inputs.interpOrder+1)
        recent_FPtfm = repeat(FPtfm_n, 1, inputs.interpOrder+1)
        recent_times = collect(range(-delta_t*inputs.interpOrder, 0.0, inputs.interpOrder+1))
    end

    #..........................................................................
    #                                MAIN LOOP
    #..........................................................................

    ### Iterate for a solution at t+dt
    for i=1:numTS-1 # we compute for the next time step, so the last step of our desired time series is computed in the second to last numTS value
        println(i)
        ## Print current simulation time to terminal
        if isinteger(t[i])
            now = Int(t[i])
            if now == 1
                println("\nSimulation Time: $now second")
            else
                println("\nSimulation Time: $now seconds")
            end
        end

        ## Check for specified rotor speed at t+dt #TODO: fix this so that it can be probably accounted for in RK4
        if inputs.topsideOn
            if (inputs.turbineStartup == 0)
                inputs.omegaControl = true #TODO: are we setting this back?
                if (inputs.usingRotorSpeedFunction) #use user specified rotor speed profile function
                    _,omegaCurrent,_ = getRotorPosSpeedAccelAtTime(t[i],t[i+1],0.0,delta_t)
                    Omega_s = omegaCurrent
                else #use discreteized rotor speed profile function
                    omegaCurrent,OmegaDotCurrent,terminateSimulation = omegaSpecCheck(t[i+1],inputs.tocp,inputs.Omegaocp,delta_t)
                    if (terminateSimulation)
                        break
                    end
                    Omega_s = omegaCurrent
                    OmegaDot_s = OmegaDotCurrent
                end
            else
                omegaCurrent = 0.0
            end

            ## Initialize "j" Gauss-Seidel iteration variables
            u_j=u_s
            udot_j=udot_s
            uddot_j=uddot_s

            azi_j = azi_s
            Omega_j = Omega_s
            OmegaDot_j = OmegaDot_s
            gb_j = gb_s
            gbDot_j = gbDot_s
            gbDotDot_j = gbDotDot_s
            genTorque_j = genTorque_s
            torqueDriveShaft_j = torqueDriveShaft_s
        end

        ## Extrapolate platform motions at t+dt to send to HydroDyn/MoorDyn #TODO: use Adams-Bashforth method if i > 4?
        if inputs.hydroOn
            top_grav_setting = copy(topModel.gravityOn)
            topModel.gravityOn = false
            u_s_prp_n = extrap_pred_vals(recent_u_prp, recent_times, t[i+1], inputs.interpOrder)
            udot_s_prp_n = extrap_pred_vals(recent_udot_prp, recent_times, t[i+1], inputs.interpOrder)
            uddot_s_prp_n = extrap_pred_vals(recent_uddot_prp, recent_times, t[i+1], inputs.interpOrder)
            u_s_prp_predState = copy(u_sp1_prp_predState)
            FPtfm_n = extrap_pred_vals(recent_FPtfm, recent_times, t[i+1], inputs.interpOrder)
            FHydro_n = copy(FHydro_n) # this doesn't have to be used outside the while loop outside of tracking purposes
            FMooring_n = copy(FMooring_n)
        end

        if inputs.topsideOn
            #TODO: put these in the model
            TOL = 1e-4  #gauss-seidel iteration tolerance for various modules
            MAXITER = 300 #max iteration for various modules
            numIterations = 1
            uNorm = 1e5
            aziNorm = 1e5
            platNorm = 0.0 #TODO: remove this?
            gbNorm = 0.0 #initialize norms for various module states
            integrator = integrator_j #don't compound integrator within the convergence loop.

            if inputs.analysisType=="GX"
                # systemout = deepcopy(system)
                strainGX = zeros(3,length(assembly.elements))
                curvGX = zeros(3,length(assembly.elements))
            end

            ## Gauss-Seidel predictor-corrector loop
            while ((uNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))
                # println("Iteration $numIterations")

                #------------------
                # GENERATOR MODULE
                #------------------
                genTorque_j = 0
                if inputs.generatorOn
                        if inputs.useGeneratorFunction
                            specifiedOmega,_,_ = omegaSpecCheck(t[i]+delta_t,inputs.tocp,inputs.Omegaocp,delta_t)
                            newVinf = FLOWMath.akima(inputs.tocp_Vinf,inputs.Vinfocp,t[i])
                            genTorqueHSS0,integrator_j,controlnamecurrent = userDefinedGenerator(newVinf,t[i],azi_j,Omega_j,OmegaHist[i],OmegaDot_j,OmegaDotHist[i],delta_t,integrator,specifiedOmega) #;operPhase
                        else
                            genTorqueHSS0 = simpleGenerator(inputs,Omega_j)
                        end
                    #should eventually account for Omega = gbDot*gearRatio here...
                    genTorque_j = genTorqueHSS0*inputs.gearRatio*inputs.gearBoxEfficiency #calculate generator torque on LSS side
                    #         genTorqueAppliedToTurbineRotor0 = -genTorque0
                    #         genTorqueAppliedToPlatform0 = genTorqueHSS0
                end

                #-------------------
                # DRIVETRAIN MODULE
                #-------------------
                torqueDriveShaft_j = genTorque_j
                gb_jLast = gb_j
                if (!inputs.omegaControl)
                    if (inputs.driveTrainOn)
                        torqueDriveShaft_j = calculateDriveShaftReactionTorque(inputs.driveShaftProps,
                        azi_j,gb_j,Omega_j*2*pi,gbDot_j*2*pi)

                        gb_j,gbDot_j,gbDotDot_j = updateRotorRotation(inputs.JgearBox,0,0,
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

                # Update rotor speed
                azi_jLast = azi_j
                if inputs.omegaControl
                    if (inputs.usingRotorSpeedFunction)
                        azi_j,Omega_j,OmegaDot_j = getRotorPosSpeedAccelAtTime(t[i],t[i+1],azi_s,delta_t)
                    else
                        Omega_j = Omega_s
                        OmegaDot_j = OmegaDot_s
                        azi_j = azi_s + Omega_j*delta_t*2*pi
                    end
                elseif !inputs.omegaControl
                    Crotor = 0
                    Krotor = 0
                    azi_j,Omega_j,OmegaDot_j = updateRotorRotation(topsideMOI[3,3],Crotor,Krotor,
                    -topFReaction_j[6],-torqueDriveShaft_j,
                    azi_s,Omega_s,OmegaDot_s,delta_t)
                else
                    error("omega control option not correctly specified")
                end

                #---------------------
                # AERODYNAMICS MODULE
                #---------------------
                # Calculate new aerodynamic loading

                # Update reference frame transformation and convert aerodynamic loads to hub reference frame
                if inputs.hydroOn
                    CN2H = calcHubRotMat(u_s_prp_predState[4:6], azi_j)
                    CN2H_no_azi = calcHubRotMat(u_s_prp_predState[4:6], 0.0)
                    # CN2H = LinearAlgebra.I(3)
                    uddot_s_prp_h = frame_convert(uddot_s_prp_n, CN2H)
                    uddot_s_prp_h[3] = -1*uddot_s_prp_h[3]
                    udot_s_prp_h = frame_convert(udot_s_prp_n, CN2H)
                    rbData = vcat(uddot_s_prp_h[1:3], udot_s_prp_h[4:6], uddot_s_prp_h[4:6])
                else
                    CN2H = calcHubRotMat(zeros(3), azi_j)
                    CN2H_no_azi = calcHubRotMat(zeros(3), 0.0)
                    rbData = zeros(Real, 9)
                end
                CH2N = LinearAlgebra.transpose(CN2H)


                #################################################################
                if !isnothing(aero)
                    if inputs.aeroLoadsOn > 0 #0 off, 1 one way, 1.5 one way with deformation from last timestep, 2 two way
                        runaero = true
                        if (inputs.aeroLoadsOn==1 || inputs.aeroLoadsOn==1.5) && numIterations!=1
                            runaero = false
                        end
                        if runaero
                            aeroVals,aeroDOFs,Xp,Yp,Zp,z3Dnorm = run_aero_with_deform(aero,deformAero,topMesh,topEl,u_j,inputs,numIterations,t[i],azi_j,Omega_j)
                        end
                    end
                end
                #################################################################
                if inputs.aeroLoadsOn > 0
                    if isnothing(aeroVals)
                        error("aeroVals must be specified if OWENS.Inputs.aeroLoadsOn")
                    elseif isnothing(aeroDOFs)
                        error("aeroDOFs must be specified if OWENS.Inputs.aeroLoadsOn")
                    end

                    if length(size(aeroVals))==1 || size(aeroVals)[2]==1 #i.e. the standard aero force input as a long array
                        # Fill in forces and dofs if they were specified not in full arrays TODO: make this more efficient
                        full_aeroVals = zeros(topMesh.numNodes*6)
                        for i_idx = 1:length(aeroDOFs)
                            full_aeroVals[Int(aeroDOFs[i_idx])] = aeroVals[i_idx]
                        end
                        aeroDOFs = collect(1:topMesh.numNodes*6)
                        for iter_i = 1:floor(Int,length(full_aeroVals)/6)
                            topFexternal[6*(iter_i-1)+1:6*(iter_i-1)+6] = frame_convert(full_aeroVals[6*(iter_i-1)+1:6*(iter_i-1)+6], CN2H_no_azi)
                        end
                    else # the other aero input as a 2D array
                        topFexternal = frame_convert(aeroVals[i+1,:], CN2H)
                    end
                else
                    topFexternal = zeros(numDOFPerNode)
                    aeroDOFs = copy(topFexternal).*0.0
                end

                #------------------------------------
                # TOPSIDE STRUCTURAL MODULE
                #------------------------------------
                # println(Float64.(rbData))
                if inputs.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                    topElStrain, topDispOut, topFReaction_j = GyricFEA.structuralDynamicsTransientROM(topModel,topMesh,topEl,topDispData1,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,top_rom,topFexternal,Int.(aeroDOFs),CN2H,rbData)
                elseif inputs.analysisType=="GX"
                    topElStrain, topDispOut, topFReaction_j,systemout  = structuralDynamicsTransientGX(topModel,topMesh,Xp,Yp,Zp,z3Dnorm,system,assembly,t,Omega_j,OmegaDot_j,delta_t,numIterations,i,strainGX,curvGX)
                else # evalulate structural dynamics using conventional representation
                    topElStrain, topDispOut, topFReaction_j = GyricFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topDispData1,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,topFexternal,Int.(aeroDOFs),CN2H,rbData;predef = topModel.nlParams.predef)
                end

                u_jLast = copy(u_j)
                u_j = topDispOut.displ_sp1
                udot_j = topDispOut.displdot_sp1
                uddot_j = topDispOut.displddot_sp1

                ## calculate norms
                uNorm = LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
                aziNorm = LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm
                gbNorm = LinearAlgebra.norm(gb_j - gb_jLast)/LinearAlgebra.norm(gb_j) #gearbox states iteration norm if it is off, the norm will be zero

                numIterations = numIterations + 1

                if inputs.analysisType=="GX" && (!(uNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                    system = deepcopy(systemout)
                end

                # Strain stiffening, save at the end of the simulation, at the last while loop iteration, mutates elStorage
                if i==numTS-1 && inputs.analysisType=="TNB" && topModel.nlParams.predef=="update" && (!(uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                    GyricFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topDispData2,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,topFexternal,Int.(aeroDOFs),CN2H,rbData;predef = topModel.nlParams.predef)
                end

                #TODO: verbosity
                # println("$(numIterations)   uNorm: $(uNorm) ")

                if numIterations==MAXITER
                    @warn "Maximum Iterations Met Breaking Iteration Loop"
                    break
                end

            end #end iteration while loop

            if inputs.analysisType=="ROM"
                eta_j = topDispOut.eta_sp1 #eta_j
                etadot_j = topDispOut.etadot_sp1 #etadot_j
                etaddot_j = topDispOut.etaddot_sp1 #etaddot_j
                topDispData1 = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1, eta_j, etadot_j, etaddot_j)
            else
                topDispData1 = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1)
            end


        end # if inputs.topsideOn

        #------------------------------------
        # COUPLED BOTTOMSIDE STRUCTURAL/HYDRO/MOORING MODULES
        #------------------------------------
        if inputs.hydroOn
            if inputs.topsideOn
                bottomFexternal = frame_convert(-1*topFReaction_j, LinearAlgebra.transpose(CN2H)) # in hub frame already
                # bottomFexternal = zeros(6)
            end

            ## Evaluate hydro-structural dynamics
            # FAST updates HD/MD using t+dt inputs extrapolated from previous time steps, NOT from the new ElastoDyn motions
            OpenFASTWrappers.HD_UpdateStates(t[i], t[i+1], u_s_prp_n, udot_s_prp_n, uddot_s_prp_n)
            if inputs.interpOrder == 1
                OpenFASTWrappers.MD_UpdateStates(0, t[i], t[i+1], u_s_prp_n, udot_s_prp_n, uddot_s_prp_n)
            elseif inputs.interpOrder == 2
                OpenFASTWrappers.MD_UpdateStates(t[i]-delta_t, t[i], t[i+1], u_s_prp_n, udot_s_prp_n, uddot_s_prp_n)
            end

            if inputs.analysisType=="ROM"
                bottomUncoupledDisps, bottomCoupledDisps, FHydro_n, FMooring_n, outVals, jac = OWENS_HD_Coupled_Solve(t[i+1], delta_t, false, jac, numDOFPerNode, prpDOFs, FPtfm_n,
                bottomDispData, bottomModel, bottomMesh, bottomEl, bottomElStorage,
                bottomFexternal, bottomFDOFs, CN2H, bottom_rom)
            else
                bottomUncoupledDisps, bottomCoupledDisps, FHydro_n, FMooring_n, outVals, jac = OWENS_HD_Coupled_Solve(t[i+1], delta_t, false, jac, numDOFPerNode, prpDOFs, FPtfm_n,
                bottomDispData, bottomModel, bottomMesh, bottomEl, bottomElStorage,
                bottomFexternal, bottomFDOFs, CN2H)
            end

            # update displacement and velocity estimates based on inputs at t+dt
            u_s_ptfm_n = bottomUncoupledDisps.displ_sp1
            u_s_prp_n = u_s_ptfm_n[prpDOFs]
            udot_s_ptfm_n = bottomUncoupledDisps.displdot_sp1
            udot_s_prp_n = udot_s_ptfm_n[prpDOFs]

            # update current acceleration estimates from the coupled solve to account for added mass
            uddot_s_ptfm_n = bottomCoupledDisps.displddot_sp1
            uddot_s_prp_n = uddot_s_ptfm_n[prpDOFs]

            # update displacement and velocity predictions for next time step
            u_sp1_prp_predState = bottomCoupledDisps.displ_sp1[prpDOFs]

            FPtfm_n = FHydro_n + FMooring_n

            u_sm1_ptfm_n = copy(u_s_ptfm_n)
            if inputs.analysisType=="ROM"
                eta_s_ptfm_n = bottomUncoupledDisps.eta_sp1 #eta_j
                etadot_s_ptfm_n = bottomUncoupledDisps.etadot_sp1 #etadot_j
                etaddot_s_ptfm_n = bottomCoupledDisps.etaddot_sp1 #etaddot_j
                bottomDispData = GyricFEA.DispData(u_s_ptfm_n, udot_s_ptfm_n, uddot_s_ptfm_n, u_sm1_ptfm_n, eta_s_ptfm_n, etadot_s_ptfm_n, etaddot_s_ptfm_n)
            else
                bottomDispData = GyricFEA.DispData(u_s_ptfm_n, udot_s_ptfm_n, uddot_s_ptfm_n, u_sm1_ptfm_n)
            end

            #------------------------------------
            # TOPSIDE STRUCTURAL MODULE W/ PLATFORM RIGID BODY MOTIONS
            #------------------------------------
            platNorm = 0.0 #TODO: remove this?
            if inputs.topsideOn
                numIterations = 1
                uNorm = 1e5
                aziNorm = 1e5
                gbNorm = 0.0 #initialize norms for various module states #TODO: this while loop for the aero side could be turned into a function, just determine what the hydro coupled adds and make it optional
                while ((uNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))

                    topModel.gravityOn = top_grav_setting
                    CN2H = calcHubRotMat(u_s_prp_predState[4:6], azi_j)
                    CN2H_no_azi = calcHubRotMat(u_s_prp_predState[4:6], 0.0)
                    # CN2P = calcHubRotMat(u_s_prp_predState[4:6], 0.0)
                    uddot_s_prp_h = frame_convert(uddot_s_prp_n, CN2H)
                    uddot_s_prp_h[3] = -1*uddot_s_prp_h[3]
                    udot_s_prp_h = frame_convert(udot_s_prp_n, CN2H)
                    rbData = vcat(uddot_s_prp_h[1:3], udot_s_prp_h[4:6], uddot_s_prp_h[4:6])

                    if !isnothing(aero)
                        if inputs.aeroLoadsOn > 0 #0 off, 1 one way, 1.5 one way with deformation from last timestep, 2 two way
                            runaero = true
                            if (inputs.aeroLoadsOn==1 || inputs.aeroLoadsOn==1.5) && numIterations!=1
                                runaero = false
                            end
                            if runaero
                                aeroVals,aeroDOFs,Xp,Yp,Zp,z3Dnorm = run_aero_with_deform(aero,deformAero,topMesh,topEl,u_j,inputs,numIterations,t[i],azi_j,Omega_j)
                            end
                        end
                    end

                    if inputs.aeroLoadsOn > 0
                        if length(size(aeroVals))==1 || size(aeroVals)[2]==1 #i.e. the standard aero force input as a long array
                            # Fill in forces and dofs if they were specified not in full arrays TODO: make this more efficient
                            full_aeroVals = zeros(topMesh.numNodes*6)
                            for i_idx = 1:length(aeroDOFs)
                                full_aeroVals[aeroDOFs[i_idx]] = aeroVals[i_idx]
                            end
                            aeroDOFs = collect(1:topMesh.numNodes*6)
                            for iter_i = 1:floor(Int,length(aeroVals)/6)
                                topFexternal[6*(iter_i-1)+1:6*(iter_i-1)+6] = frame_convert(full_aeroVals[6*(iter_i-1)+1:6*(iter_i-1)+6], CN2H_no_azi)
                            end
                        else
                            topFexternal = frame_convert(aeroVals[i+1,:], CN2H)
                        end
                    end

                    if inputs.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                        topElStrain, topDispOut, topFReaction_j = GyricFEA.structuralDynamicsTransientROM(topModel,topMesh,topEl,topDispData2,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,top_rom,topFexternal,Int.(aeroDOFs),CN2H,rbData)
                    elseif inputs.analysisType=="GX"
                        topElStrain, topDispOut, topFReaction_j,systemout  = structuralDynamicsTransientGX(topModel,topmesh,Xp,Yp,Zp,z3Dnorm,system,assembly,t,Omega_j,OmegaDot_j,delta_t,numIterations,i,strainGX,curvGX)
                    else # evalulate structural dynamics using conventional representation
                        topElStrain, topDispOut, topFReaction_j = GyricFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topDispData2,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,topFexternal,Int.(aeroDOFs),CN2H,rbData;predef = topModel.nlParams.predef)
                    end

                    u_jLast = copy(u_j)
                    u_j = topDispOut.displ_sp1
                    udot_j = topDispOut.displdot_sp1
                    uddot_j = topDispOut.displddot_sp1

                    if inputs.analysisType=="GX" && (!(uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                        system = deepcopy(systemout)
                    end

                    # Strain stiffening, save at the end of the simulation, at the last while loop iteration, mutates elStorage
                    if i==numTS-1 && inputs.analysisType=="TNB" && topModel.nlParams.predef=="update" && (!(uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                        GyricFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topDispData2,Omega_s,OmegaDot_s,t[i+1],delta_t,topElStorage,topFexternal,Int.(aeroDOFs),CN2H,rbData;predef = topModel.nlParams.predef)
                    end

                    ## calculate norms
                    uNorm = 0.0#LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
                    aziNorm = 0.0#LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm
                    gbNorm = 0.0#LinearAlgebra.norm(gb_j - gb_jLast)/LinearAlgebra.norm(gb_j) #gearbox states iteration norm if it is off, the norm will be zero

                    numIterations = numIterations + 1
                    if numIterations==MAXITER
                        @warn "Maximum Iterations Met Breaking Iteration Loop"
                        break
                    end

                end # iteration while loop

                if inputs.analysisType=="ROM"
                    eta_j = topDispOut.eta_sp1 #eta_j
                    etadot_j = topDispOut.etadot_sp1 #etadot_j
                    etaddot_j = topDispOut.etaddot_sp1 #etaddot_j
                    topDispData2 = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1, eta_j, etadot_j, etaddot_j)
                else
                    topDispData2 = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1)
                end

            end # if inputs.topsideOn

        end # if inputs.hydroOn

        ## update timestepping variables and other states, store in history arrays
        if inputs.topsideOn
            ## calculate converged generator torque/power
            genTorquePlot = 0
            if (inputs.useGeneratorFunction)
                if (inputs.generatorOn || (inputs.turbineStartup==0))
                    genTorqueHSS0 = simpleGenerator(inputs,Omega_j)
                end
            end
            genPowerPlot = genTorquePlot*(gbDot_j*2*pi)*inputs.gearRatio

            u_sm1 = copy(u_s)
            u_s = u_j
            udot_s = udot_j
            uddot_s = uddot_j

            azi_s = azi_j
            Omega_s = Omega_j
            OmegaDot_s = OmegaDot_j

            genTorque_s = genTorque_j
            torqueDriveShaft_s = torqueDriveShaft_j

            gb_s = gb_j
            gbDot_s = gbDot_j
            gbDotDot_s = gbDotDot_j

            uHist[i+1,:] = u_s
            FReactionHist[i+1,:] = topFReaction_j
            topFexternal_hist[i+1,1:length(topFexternal)] = topFexternal
            # FTwrBsHist[i+1,:] = -topFReaction_j  + topFWeight_j

            aziHist[i+1] = azi_s
            OmegaHist[i+1] = Omega_s
            OmegaDotHist[i+1] = OmegaDot_s

            gbHist[i+1] = gb_s
            gbDotHist[i+1] = gbDot_s
            gbDotDotHist[i+1] = gbDotDot_s

            #genTorque[i+1] = genTorque_s
            genTorque[i+1] = genTorquePlot
            genPower[i+1] = genPowerPlot
            torqueDriveShaft[i+1] = torqueDriveShaft_s

            if inputs.analysisType=="GX"
                strainGX = topElStrain[1]
                curvGX = topElStrain[2]
                for ii = 1:length(strainGX[1,:])
                    epsilon_x_hist[:,ii,i] .= strainGX[1,ii]
                    kappa_y_hist[:,ii,i] .= curvGX[2,ii]
                    kappa_z_hist[:,ii,i] .= curvGX[3,ii]
                    epsilon_z_hist[:,ii,i] .= strainGX[3,ii]
                    kappa_x_hist[:,ii,i] .= curvGX[1,ii]
                    epsilon_y_hist[:,ii,i] .= strainGX[2,ii]
                end
            else
                for ii = 1:length(topElStrain)
                    epsilon_x_hist[:,ii,i] = topElStrain[ii].epsilon_x
                    kappa_y_hist[:,ii,i] = topElStrain[ii].kappa_y
                    kappa_z_hist[:,ii,i] = topElStrain[ii].kappa_z
                    epsilon_z_hist[:,ii,i] = topElStrain[ii].epsilon_z
                    kappa_x_hist[:,ii,i] = topElStrain[ii].kappa_x
                    epsilon_y_hist[:,ii,i] = topElStrain[ii].epsilon_y
                end
            end
            ## check rotor speed for generator operation
            if Omega_s >= rotorSpeedForGenStart
                inputs.generatorOn = true
            else
                inputs.generatorOn = false
            end

        end

        if inputs.hydroOn
            uHist_prp[i+1,:] = u_s_prp_n
            FPtfmHist[i+1,:] = FPtfm_n
            FHydroHist[i+1,:] = FHydro_n
            FMooringHist[i+1, :] = FMooring_n

            # Shift up window of extrapolation vectors
            recent_u_prp = hcat(recent_u_prp[:, 2:end], u_s_prp_n)
            recent_udot_prp = hcat(recent_udot_prp[:, 2:end], udot_s_prp_n)
            recent_uddot_prp = hcat(recent_uddot_prp[:, 2:end], uddot_s_prp_n)
            recent_FPtfm = hcat(recent_FPtfm[:, 2:end], FPtfm_n)
            recent_times = vcat(recent_times[2:end], t[i+1])
        end

    end #end timestep loop

    println("Simulation Complete.")

    # End FAST module links
    if inputs.hydroOn
        OpenFASTWrappers.HD_End()
        OpenFASTWrappers.MD_End()
    end

    outputData(inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist)

    return t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,topFexternal_hist
end

function structuralDynamicsTransientGX(topModel,mesh,Xp,Yp,Zp,z3Dnorm,system,assembly,t,Omega_j,OmegaDot_j,delta_t,numIterations,i,strainGX,curvGX)
    distributed_loads = Dict()

    height = maximum(mesh.z)
    for jbld = 1:length(mesh.structuralElNumbers[:,1])
        XpGXspl1 = FLOWMath.Akima(z3Dnorm.*height,Xp[jbld,end,:])#,GXz/maximum(GXz))
        YpGXspl1 = FLOWMath.Akima(z3Dnorm.*height,Yp[jbld,end,:])#,GXz/maximum(GXz))
        ZpGXspl1 = FLOWMath.Akima(z3Dnorm.*height,Zp[jbld,end,:])#,GXz/maximum(GXz))

        for ipt = Int.(mesh.structuralNodeNumbers[jbld,1]:mesh.structuralNodeNumbers[jbld,end-1])
            iel = findfirst(x->x==ipt,assembly.start)
            s1 = assembly.points[assembly.start[ipt]][3]
            s2 = assembly.points[assembly.stop[ipt]][3]
            if !isempty(iel)
                distributed_loads[iel] = GXBeam.DistributedLoads(assembly,iel;s1,s2,fx = (s) -> XpGXspl1(s),
                    fy = (s) -> YpGXspl1(s))#, fz = (s) -> ZpGXspl1(s))
            else
                println("Empty at $ipt")
            end
        end
    end
    #TODO: pull in parametrically from inputs.BC info
    prescribed_conditions = Dict(
        # fixed base
        1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # fixed top, but free to rotate around z-axis
        # 50 => GXBeam.PrescribedConditions(Fx = 1e4*sin(20*t)),
        23 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0),
        )

    linear_velocity = [0.0,0.0,0.0]
    angular_velocity = [0.0,0.0,Omega_j*2*pi]
    linear_acceleration = [0.0,0.0,0.0]
    angular_acceleration = [0.0,0.0,OmegaDot_j*2*pi]

    tvec = [t[i],t[i]+delta_t]

    if eltype(topModel.gravityOn) == Bool && topModel.gravityOn == true
        gravity = [0.0,0.0,-9.81]
    elseif eltype(topModel.gravityOn) == Bool && topModel.gravityOn == false
        gravity = [0.0,0.0,0.0]
    end

    if eltype(topModel.gravityOn) == Float64
        gravity = feamodel.gravityOn
    end


    reset_state = false
    initialize = false

    if i == 1 && numIterations == 1

        GXBeam.steady_state_analysis!(system,assembly; prescribed_conditions, #TODO: just gravity and rotation for the first solve? Could add on forces in next solve below
        gravity,angular_velocity,linear=false,reset_state=true)
    end

    systemout, history, converged = GXBeam.time_domain_analysis!(deepcopy(system),assembly, tvec;
    reset_state,initialize,linear_velocity,angular_velocity,linear_acceleration,
    angular_acceleration,prescribed_conditions,gravity,linear=false)#!topModel.nlOn)

    if !converged
        println("GX failed to converge")
    end
    # elStrain
    state = GXBeam.AssemblyState(systemout, assembly;
    prescribed_conditions)

    for iel = 1:length(state.elements)
        strainGX[:,iel] = GXBeam.element_strain(assembly.elements[iel],state.elements[iel].Fi,state.elements[iel].Mi)
        curvGX[:,iel] = GXBeam.element_curvature(assembly.elements[iel],state.elements[iel].Fi,state.elements[iel].Mi)
    end

    # disp
    disp_sp1 = zeros(length(history[end].points)*6)
    idx = 1
    for ipt = 1:length(history[end].points)
        for iu = 1:3
            disp_sp1[idx] = history[end].points[ipt].u[iu]
            idx += 1
        end
        for itheta = 1:3
            disp_sp1[idx] = history[end].points[ipt].theta[itheta]
            idx += 1
        end
    end
    dispOut = GyricFEA.DispOut(nothing, disp_sp1,copy(disp_sp1).*0.0,copy(disp_sp1).*0.0)
    FReaction_j = [-history[end].points[1].F[1];-history[end].points[1].F[2];-history[end].points[1].F[3];
    -history[end].points[1].M[1];-history[end].points[1].M[2];-history[end].points[1].M[3]]
    return (strainGX,curvGX), dispOut, FReaction_j,systemout
end

function run_aero_with_deform(aero,deformAero,mesh,el,u_j,inputs,numIterations,t_i,azi_j,Omega_j)

    if inputs.tocp_Vinf == -1
        newVinf = -1
    else
        newVinf = FLOWMath.akima(inputs.tocp_Vinf,inputs.Vinfocp,t_i)
    end

    if inputs.aeroLoadsOn > 1
        # Transform Global Displacements to Local
        numDofPerNode = 6

        u_j_local = zeros(Int(max(maximum(mesh.structuralNodeNumbers))*6))
        for jbld = 1:length(mesh.structuralElNumbers[:,1])
            for kel = 1:length(mesh.structuralElNumbers[1,:])-1
                # orientation angle,xloc,sectionProps,element order]
                elNum = Int(mesh.structuralElNumbers[jbld,kel])
                #get dof map
                node1 = Int(mesh.structuralNodeNumbers[jbld,kel])
                node2 = Int(mesh.structuralNodeNumbers[jbld,kel+1])
                dofList = [(node1-1)*numDofPerNode.+(1:6);(node2-1)*numDofPerNode.+(1:6)]

                globaldisp = u_j[dofList]

                twist = el.props[elNum].twist
                sweepAngle = el.psi[elNum]
                coneAngle = el.theta[elNum]
                rollAngle = el.roll[elNum]

                twistAvg = rollAngle + 0.5*(twist[1] + twist[2])
                lambda = GyricFEA.calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)
                localdisp = lambda*globaldisp

                #asssembly
                for m = 1:length(dofList)
                    u_j_local[dofList[m]] =  u_j_local[dofList[m]]+localdisp[m]
                end

            end
        end

        disp_x = [u_j[i] for i = 1:6:length(u_j)]
        disp_y = [u_j[i] for i = 2:6:length(u_j)]
        disp_z = [u_j[i] for i = 3:6:length(u_j)]
        disp_twist = [u_j_local[i] for i = 4:6:length(u_j_local)]
        # disp_twist2 = [u_j_local[i] for i = 5:6:length(u_j_local)]
        # disp_twist3 = [u_j_local[i] for i = 6:6:length(u_j_local)]

        bld_x = copy(mesh.structuralNodeNumbers[:,1:end-1]).*0.0
        bld_y = copy(mesh.structuralNodeNumbers[:,1:end-1]).*0.0
        bld_z = copy(mesh.structuralNodeNumbers[:,1:end-1]).*0.0
        bld_twist = copy(mesh.structuralNodeNumbers[:,1:end-1]).*0.0

        for jbld = 1:length(mesh.structuralNodeNumbers[:,1])
            bld_indices = Int.(mesh.structuralNodeNumbers[jbld,1:end-1])
            bld_x[jbld,:] = mesh.x[bld_indices]+disp_x[bld_indices]
            bld_y[jbld,:] = mesh.y[bld_indices]+disp_y[bld_indices]
            bld_z[jbld,:] = mesh.z[bld_indices]+disp_z[bld_indices]
            # flatten blade x,y
            bld_x[jbld,:] = sqrt.(bld_x[jbld,:].^2 .+bld_y[jbld,:].^2) #TODO: a better way via the blade offset azimuth?
            bld_twist[jbld,:] = -disp_twist[bld_indices] #the bending displacements are in radians
            # The local structural FOR follows right hand rule, so at the bottom of the blade, the x-vector is pointing outward, and a positive
            # rotation about x would make the blade twist into the turbine.  In AC and DMS, if we are looking at say Andrew's 2016 paper, fig 10,
            # the blade has looped up and is pointing at us, so positive twist would increase the aoa.  In the AC and DMS equations, aoa is decreased by twist
            # so, we should negate
            # PyPlot.figure(111)
            # PyPlot.plot(mesh.x[bld_indices]+disp_x[bld_indices],mesh.z[bld_indices]+disp_z[bld_indices].*1,label="disp")
            # PyPlot.plot(bld_twist[jbld,:]*180/pi*1,mesh.z[bld_indices]+disp_z[bld_indices].*1,label="twist1")
            # PyPlot.plot(disp_twist2[bld_indices]*180/pi*1,mesh.z[bld_indices]+disp_z[bld_indices].*1,label="twist2")
            # PyPlot.plot(disp_twist3[bld_indices]*180/pi*1,mesh.z[bld_indices]+disp_z[bld_indices].*1,label="twist3")
            # PyPlot.legend()
            # sleep(.1)
        end

    else
        bld_x = -1
        bld_z = -1
        bld_twist = -1
    end

    # println("Calling Aero $(Omega_j*60) RPM $newVinf Vinf")
    deformAero(azi_j;newOmega=Omega_j*2*pi,newVinf,bld_x,bld_z,bld_twist) #TODO: implement deformation induced velocities
    aeroVals,aeroDOFs,Xp,Yp,Zp,z3Dnorm = aero(t_i,azi_j)
    return aeroVals,aeroDOFs,Xp,Yp,Zp,z3Dnorm
end

function outputData(inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist)

    #Writefile
    if inputs.outFilename!="none"
        println("WRITING Verification File to $(inputs.outFilename)")

        filename = string(inputs.outFilename[1:end-3], "h5")
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
            HDF5.write(file,"genTorque",genTorque)
            HDF5.write(file,"genPower",genPower)
            HDF5.write(file,"torqueDriveShaft",torqueDriveShaft)
            HDF5.write(file,"uHist",uHist)
            HDF5.write(file,"uHist_prp",uHist_prp)
            HDF5.write(file,"epsilon_x_hist",epsilon_x_hist)
            HDF5.write(file,"epsilon_y_hist",epsilon_y_hist)
            HDF5.write(file,"epsilon_z_hist",epsilon_z_hist)
            HDF5.write(file,"kappa_x_hist",kappa_x_hist)
            HDF5.write(file,"kappa_y_hist",kappa_y_hist)
            HDF5.write(file,"kappa_z_hist",kappa_z_hist)
        end
    end
end

function allocate_topside(inputs,topMesh,topEl,topModel,numDOFPerNode,u_s)
    if isnothing(topModel)
        error("topMesh must be specified if OWENS.Inputs.topsideOn")
    elseif isnothing(topMesh)
        error("topMesh must be specified if OWENS.Inputs.topsideOn")
    elseif isnothing(topEl)
        error("topEl must be specified if OWENS.Inputs.topsideOn")
    end

    top_totalNumDOF = topMesh.numNodes*numDOFPerNode
    if isnothing(u_s)
        u_s = zeros(top_totalNumDOF)
    end
    u_s = GyricFEA.setInitialConditions(topModel.initCond, u_s, numDOFPerNode)
    udot_s = copy(u_s).*0.0
    uddot_s = copy(u_s).*0.0
    u_sm1 = copy(u_s)

    topDispData1 = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1)
    topDispData2 = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1)
    topElStrain = fill(GyricFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), topMesh.numEl)

    gb_s = 0
    gbDot_s = 0
    gbDotDot_s = 0
    azi_s = 0
    Omega_s = copy(inputs.OmegaInit)
    OmegaDot_s = 0
    genTorque_s = 0
    torqueDriveShaft_s = 0
    topFexternal = zeros(top_totalNumDOF)
    topFexternal_hist = zeros(Int(inputs.numTS),top_totalNumDOF)
    return u_s,udot_s,uddot_s,u_sm1,topDispData1,topDispData2,topElStrain,gb_s,gbDot_s,gbDotDot_s,azi_s,Omega_s,OmegaDot_s,genTorque_s,torqueDriveShaft_s,topFexternal,topFexternal_hist
end


function allocate_bottom(t,numTS,delta_t,inputs,bottomMesh,bottomEl,bottomModel,bin,numDOFPerNode)
    if isnothing(bottomModel)
        error("bottomMesh must be specified if OWENS.Inputs.hydroOn")
    elseif isnothing(bottomMesh)
        error("bottomMesh must be specified if OWENS.Inputs.hydroOn")
    elseif isnothing(bottomEl)
        error("bottomEl must be specified if OWENS.Inputs.hydroOn")
    elseif isnothing(bin)
        error("bin must be specified if OWENS.Inputs.hydroOn")
    end

    bottom_totalNumDOF = bottomMesh.numNodes*numDOFPerNode
    u_s_ptfm_n = zeros(Real, bottom_totalNumDOF)
    u_s_ptfm_n = GyricFEA.setInitialConditions(bottomModel.initCond, u_s_ptfm_n, numDOFPerNode)
    udot_s_ptfm_n = copy(u_s_ptfm_n).*0.0
    uddot_s_ptfm_n = copy(u_s_ptfm_n).*0.0
    u_sm1_ptfm_n = copy(u_s_ptfm_n)
    bottomDispData = GyricFEA.DispData(u_s_ptfm_n, udot_s_ptfm_n, uddot_s_ptfm_n, u_sm1_ptfm_n)

    prpDOFs = collect(1:6) #TODO: add this to bottomModel
    u_s_prp_n = Vector(u_s_ptfm_n[prpDOFs]) # the hub and global reference frames start as the same, so we don't need to frame_convert here
    udot_s_prp_n = Vector(udot_s_ptfm_n[prpDOFs])
    uddot_s_prp_n = Vector(uddot_s_ptfm_n[prpDOFs])

    jac = Array{Real}(undef, numDOFPerNode*2, numDOFPerNode*2)
    numMooringLines = 3

    if inputs.outFilename == "none"
        hd_outFilename = "hydrodyn_temp.out"
    else
        hd_outFilename = inputs.outFilename
    end

    FHydro_n = zeros(Float32, numDOFPerNode) #Vector{Real}(undef, numDOFPerNode)
    FMooring_n = zeros(Float32, numDOFPerNode) #Vector{Real}(undef, numDOFPerNode)
    outVals = Vector{Float32}(undef, numDOFPerNode+1) # Rigid body displacement in 6DOF + wave elevation
    mooringTensions = Vector{Float32}(undef, numMooringLines*2) # Fairlead + anchor tension for each line

    OpenFASTWrappers.HD_Init(bin.hydrodynLibPath, hd_outFilename, hd_input_file=inputs.hd_input_file, PotFile=inputs.potflowfile, t_initial=t[1], dt=delta_t, t_max=t[1]+(numTS-1)*delta_t, interp_order=inputs.interpOrder)
    OpenFASTWrappers.MD_Init(bin.moordynLibPath, md_input_file=inputs.md_input_file, init_ptfm_pos=u_s_prp_n, interp_order=inputs.interpOrder, WtrDpth=200)

    return bottom_totalNumDOF,u_s_ptfm_n,udot_s_ptfm_n,uddot_s_ptfm_n,u_sm1_ptfm_n,bottomDispData,prpDOFs,u_s_prp_n,udot_s_prp_n,uddot_s_prp_n,jac,numMooringLines,FHydro_n,FMooring_n,outVals,mooringTensions
end

"""

OWENS_HD_Coupled_Solve(time, dt, calcJacobian, jac, numDOFPerNode, ptfm_dofs, FPtfm_old, FMooring_new,
dispIn, topModel, mesh, el, elStorage, Omega, OmegaDot,
other_Fexternal, other_Fdof, CN2H, u_ptfm_n, udot_ptfm_n, uddot_ptfm_n, rom=0)

Internal, performs the input-output solve procedure between OWENS and HydroDyn. This is required
due to the close interdependency between the rigid body acceleration of the platform (calculated in
OWENS/GyricFEA) and the hydrodynamic forcing due to the platform added mass (calculated in HydroDyn).
This is basically a mirror of the `ED_HD_InputOutputSolve` subroutine from OpenFAST, only replacing
ElastoDyn with GyricFEA.

# Input
* `time::float`: current simulation time
* `dt::float`: simulation time step
* `calcJacobian::bool`: flag on whether or not to calculate the force/acceleration Jacobians (needs to be true at initialization)
* `jac::Array{<:float}`: existing array containing the platform hydro force and acceleration Jacobians
* `numDOFPerNode::int`: number of degrees of freedom per node (typically 6)
* `ptfm_dofs::Vector{<:int}`: indices of the nodal vector representing the platform (typically 1 through 6)
* `FPtfm_old::Vector{<:float}`: previously calculated total platform loads in the inertial frame (N)
* `FMooring_new::Vector{<:float}`: mooring loads calculated in the current time step/correction in the inertial frame (N)
* `dispIn::DispData`: see ?GyricFEA.DispData
* `topModel::FEAModel`: see ?GyricFEA.FEAModel
* `mesh::Mesh`: see ?GyricFEA.Mesh
* `el::El`: see ?GyricFEA.El,
* `elStorage::ElStorage`: see ?GyricFEA.elStorage
* `Omega::float`: rotor speed (Hz)
* `OmegaDot::float`: rotor acceleration (Hz/s)
* `other_Fexternal::Vector{<:float}`: other external loads acting on the structure besides at the platform (typically aero loads) (N)
* `other_Fdof::Vector{<:float}`: vector of nodes with indices corresponding to the loads given in other_Fexternal
* `CN2H::Array{<:float}`: rotation matrix defining the transformation from the inertial to the hub reference frames
* `u_ptfm_n::Vector{<:float}`: the platform position in the inertial reference frame at the current simulation time (m)
* `udot_ptfm_n::Vector{<:float}`: the platform velocity in the inertial reference frame at the current simulation time (m/s)
* `uddot_ptfm_n::Vector{<:float}`: the platform acceleration in the inertial reference frame at the current simulation time (m/s/s)
* `rom::int/ROM`: see ?GyricFEA.ROM. Instead defaults to 0 if the reduced order model is not used

# Output
* `elStrain_out`: object containing strains on each element returned by GyricFEA
* `dispOut`: object containing motions on each element returned by GyricFEA
* `FReaction_out`: vector of the total loads acting on the specified nodes
* `FHydro_out`: vector of the hydrodynamic loads acting on the platform returned by HydroDyn
* `out_vals`: vector containing the output parameters specified at the end of the HydroDyn input file
* `jac_out`: the platform hydro forces and accelerations Jacobian calculated internally if calcJacobian=true (passthrough if calcJacobian=false)
"""

function OWENS_HD_Coupled_Solve(time, dt, calcJacobian, jac, numDOFPerNode, prpDOFS, FPtfm_old,
    dispIn, topModel, mesh, el, elStorage,
    other_Fexternal, other_Fdof, CN2H, rom=0)

    # !!! Make sure structuralDynamicsTransient and MD_CalcOutput have run before this so you can send dispOut and FMooring_new here !!!
    # !!! other_Fexternal is assumed to already be in the hub reference frame !!!
    #
    # allocate new variables
    FHydro_2 = Vector{Float32}(undef, numDOFPerNode)
    FMooring_new = Vector{Float32}(undef, numDOFPerNode)
    FHydro_new = Vector{Float32}(undef, numDOFPerNode)
    outVals = Vector{Float32}(undef, numDOFPerNode+1)
    mooringTensions = Vector{Float32}(undef, 6) # Fairlead + anchor tension for each line TODO: don't hardcode this
    u = Vector{Float64}(undef, numDOFPerNode*2)

    FMultiplier = 1e6

    # Calculate outputs at the current time, based on inputs at the current time
    total_Fexternal = vcat(other_Fexternal, FPtfm_old) # platform loads are old inputs from last time step/correction
    total_Fdof = vcat(other_Fdof, prpDOFS)
    if rom != 0
        _ ,dispsUncoupled, _ = GyricFEA.structuralDynamicsTransientROM(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,rom,total_Fexternal,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
    else
        _ ,dispsUncoupled, _ = GyricFEA.structuralDynamicsTransient(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,total_Fexternal,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
    end
    uddot_prp_n = dispsUncoupled.displddot_sp1[prpDOFS]
    if time > 0
        udot_prp_n = dispsUncoupled.displdot_sp1[prpDOFS]
        u_prp_n = dispsUncoupled.displ_sp1[prpDOFS]
    else
        udot_prp_n = zeros(Real, numDOFPerNode)
        u_prp_n = zeros(Real, numDOFPerNode)
    end
    FHydro_2[:], _ = OpenFASTWrappers.HD_CalcOutput(time, u_prp_n, udot_prp_n, uddot_prp_n, FHydro_2, outVals)
    FMooring_new[:], _ = OpenFASTWrappers.MD_CalcOutput(time, u_prp_n, udot_prp_n, uddot_prp_n, FMooring_new, mooringTensions)

    # Calculate the residual
    # For consistency, everything in the u vector is in the inertial reference frame
    u[1:numDOFPerNode] = FPtfm_old / FMultiplier
    u[numDOFPerNode+1:numDOFPerNode*2] = uddot_prp_n
    residual = calcHydroResidual(uddot_prp_n, FHydro_2, FMooring_new, u, FMultiplier)  # note that residual[7:12] will always be zero.
    # Ordinarily, the accelerations would account for differences in motions once aerodynamics are factored in.
    # However, since the aerodynamics aren't included here in this meshing approach, we can remove a run of structuralDynamicsTransient by hardcoding it like this.

    # Calculate the Jacobian. Since there's no single function associated with the motions/forces, we will manually
    # perturb each degree of freedom one at a time and recalculate the outputs so we can get a full gradient.
    if calcJacobian

        for dof = collect(1:numDOFPerNode) # forces
            FPtfm_perturb = copy(FPtfm_old)
            u_perturb = copy(u)
            FPtfm_perturb[dof] += 1E6
            u_perturb[dof] += 1
            total_Fexternal_perturb = vcat(other_Fexternal, FPtfm_perturb)
            if !isa(rom, Number)
                _ ,dispOut_perturb, _ = GyricFEA.structuralDynamicsTransientROM(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,rom,total_Fexternal_perturb,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
            else
                _, dispOut_perturb, _ = GyricFEA.structuralDynamicsTransient(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,total_Fexternal_perturb,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
            end
            uddot_ptfm_perturb_n = dispOut_perturb.displddot_sp1[prpDOFS]
            residual_perturb = calcHydroResidual(uddot_ptfm_perturb_n, FHydro_2, FMooring_new, u_perturb, FMultiplier)
            jac[:,dof] = residual_perturb - residual
        end

        for dof = collect(1:numDOFPerNode) # accelerations
            uddot_prp_perturb = copy(uddot_prp_n)
            u_perturb = copy(u)
            uddot_prp_perturb[dof] += 1
            u_perturb[dof+numDOFPerNode] += 1
            FHydro_perturb = Vector{Float32}(undef, numDOFPerNode) # this is reset each time, otherwise HydroDyn returns garbage values after the first iteration
            FHydro_perturb[:], _ = OpenFASTWrappers.HD_CalcOutput(time, u_prp_n, udot_prp_n, uddot_prp_perturb, FHydro_perturb, outVals)
            residual_perturb = calcHydroResidual(uddot_prp_n, FHydro_perturb, FMooring_new, u_perturb, FMultiplier)
            jac[:,dof+numDOFPerNode] = residual_perturb - residual
        end

    end  # if calcJacobian

    # Solve for delta_u: jac*delta_u = -residual
    delta_u = jac\-residual
    jac_out = jac

    # Update inputs #TODO: maybe update u and udot for HD as well using Newmark approach?
    # del = 0.5
    # alp = 0.167
    u += delta_u
    FPtfm_rev = FPtfm_old + delta_u[1:numDOFPerNode]*FMultiplier
    uddot_prp_n_rev = uddot_prp_n + delta_u[numDOFPerNode+1:end]
    # udot_prp_n_rev = udot_prp_n + (1-del)*dt*uddot_prp_n + del*dt*uddot_prp_n_rev
    # u_prp_n_rev = u_prp_n + udot_prp_n*dt + (0.5 - alp)*dt^2*uddot_prp_n + alp*dt^2*uddot_prp_n_rev

    total_Fexternal_rev = vcat(other_Fexternal, FPtfm_rev)

    # Rerun OWENS and HydroDyn with updated inputs
    if rom != 0
        _, dispsCoupled, _ = GyricFEA.structuralDynamicsTransientROM(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,rom,total_Fexternal_rev,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
    else
        _, dispsCoupled, _ = GyricFEA.structuralDynamicsTransient(topModel,mesh,el,dispIn,0.0,0.0,time,dt,elStorage,total_Fexternal_rev,Int.(total_Fdof),LinearAlgebra.I(3),zeros(Real, 9))
    end
    FHydro_new[:], outVals[:] = OpenFASTWrappers.HD_CalcOutput(time, u_prp_n, udot_prp_n, uddot_prp_n_rev, FHydro_new, outVals)


    return dispsUncoupled, dispsCoupled, FHydro_new, FMooring_new, outVals, jac_out

end

function initialize_generator!(inputs)
    if (inputs.turbineStartup == 1) #forced start-up using generator as motor
        println("Running in forced starting mode.")
        inputs.generatorOn = true  #TODO: clean this redundant/conflicting logic up
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = 0.0
    elseif (inputs.turbineStartup == 2) #self-starting mode
        println("Running in self-starting mode.")
        inputs.generatorOn = false
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = inputs.OmegaGenStart #Spec rotor speed for generator startup Hz
    else
        println("Running in specified rotor speed mode")
        inputs.generatorOn = false
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = 1e6 #ensures generator always off for practical purposes
    end
    return rotorSpeedForGenStart
end

function initialize_ROM(IElStorage,IModel,IMesh,IEl,u_s,udot_s,uddot_s)
    I_rom = GyricFEA.reducedOrderModel(IElStorage,IModel,IMesh,IEl,u_s) #construct reduced order model

    #set up inital values in modal space
    IjointTransformTrans = IModel.jointTransform'
    u_sRed = IjointTransformTrans*u_s
    udot_sRed = IjointTransformTrans*udot_s
    uddot_sRed = IjointTransformTrans*uddot_s

    IBC = IModel.BC
    u_s2 = GyricFEA.applyBCModalVec(u_sRed,IBC.numpBC,IBC.map)
    udot_s2 = GyricFEA.applyBCModalVec(udot_sRed,IBC.numpBC,IBC.map)
    uddot_s2 = GyricFEA.applyBCModalVec(uddot_sRed,IBC.numpBC,IBC.map)

    I_invPhi = I_rom.invPhi

    eta_s      = I_invPhi*u_s2
    etadot_s   = I_invPhi*udot_s2
    etaddot_s  = I_invPhi*uddot_s2
    return I_rom,IjointTransformTrans,u_sRed,udot_sRed,uddot_sRed,IBC,u_s2,udot_s2,uddot_s2,I_invPhi,eta_s,etadot_s,etaddot_s
end

function hydro_topside_nodal_coupling!(bottomModel,bottomMesh,topsideMass,topModel,topsideCG,topsideMOI,numDOFPerNode)
    ## Implement the mass matrix of the topside as a concentrated mass on the bottom side
    topsideMassMat = [ #TODO: I'm sure there's a more efficient way of doing this
    topsideMass                0.0                        0.0                        0.0                        topsideMass*topsideCG[3]  -topsideMass*topsideCG[2]
    0.0                        topsideMass                0.0                       -topsideMass*topsideCG[3]   0.0                        topsideMass*topsideCG[1]
    0.0                        0.0                        topsideMass                topsideMass*topsideCG[2]  -topsideMass*topsideCG[1]   0.0
    0.0                       -topsideMass*topsideCG[3]   topsideMass*topsideCG[2]   topsideMOI[1,1]            topsideMOI[1,2]            topsideMOI[1,3]
    -topsideMass*topsideCG[3]   0.0                       -topsideMass*topsideCG[1]   topsideMOI[2,1]            topsideMOI[2,2]            topsideMOI[2,3]
    -topsideMass*topsideCG[2]   topsideMass*topsideCG[1]   0.0                        topsideMOI[3,1]            topsideMOI[3,2]            topsideMOI[3,3]
    # topsideMass                0.0                        0.0                        0.0                        0.0                        0.0
    # 0.0                        topsideMass                0.0                        0.0                        0.0                        0.0
    # 0.0                        0.0                        topsideMass                0.0                        0.0                        0.0
    # 0.0                        0.0                        0.0                        topsideMOI[1,1]            topsideMOI[1,2]            topsideMOI[1,3]
    # 0.0                        0.0                        0.0                        topsideMOI[2,1]            topsideMOI[2,2]            topsideMOI[2,3]
    # 0.0                        0.0                        0.0                        topsideMOI[3,1]            topsideMOI[3,2]            topsideMOI[3,3]
    ]

    nodalinputdata = Array{Any}(undef, 36,5)

    idx = 0
    for ii = 1:6
        for jj = 1:6
            idx += 1
            nodalinputdata[idx,:] = [length(bottomMesh.z) "M6" ii jj topsideMassMat[ii,jj]]
        end
    end

    topsideConcTerms = GyricFEA.applyConcentratedTerms(bottomMesh.numNodes, numDOFPerNode, data=nodalinputdata, jointData=topModel.joint,applyTop2Bot=true)

    bottomModel.nodalTerms.concMass += topsideConcTerms.concMass
end

function allocate_general(inputs,topMesh,numDOFPerNode,numTS,assembly)

    ## History array initialization
    if inputs.topsideOn
        uHist = zeros(Real, numTS, topMesh.numNodes*numDOFPerNode)
        if inputs.analysisType == "GX"
            nel = length(assembly.start) #TODO: these should be the same.
        else
            nel = topMesh.numEl
        end
        epsilon_x_hist = zeros(Real, 4,nel,numTS)
        epsilon_y_hist = zeros(Real, 4,nel,numTS)
        epsilon_z_hist = zeros(Real, 4,nel,numTS)
        kappa_x_hist = zeros(Real, 4,nel,numTS)
        kappa_y_hist = zeros(Real, 4,nel,numTS)
        kappa_z_hist = zeros(Real, 4,nel,numTS)
    else
        uHist = zeros(Real, numTS, numDOFPerNode)
        epsilon_x_hist = zeros(Real, 4,1, numTS)
        epsilon_y_hist = zeros(Real, 4,1, numTS)
        epsilon_z_hist = zeros(Real, 4,1, numTS)
        kappa_x_hist = zeros(Real, 4,1, numTS)
        kappa_y_hist = zeros(Real, 4,1, numTS)
        kappa_z_hist = zeros(Real, 4,1, numTS)
    end

    FReactionHist = zeros(Real, numTS,6) #TODO: these can go with the general init
    FTwrBsHist = zeros(Real, numTS, 6)
    aziHist = zeros(Real, numTS)
    OmegaHist = zeros(Real, numTS)
    OmegaDotHist = zeros(Real, numTS)
    gbHist = zeros(Real, numTS)
    gbDotHist = zeros(Real, numTS)
    gbDotDotHist = zeros(Real, numTS)
    genTorque = zeros(Real, numTS)
    genPower = zeros(Real, numTS)
    torqueDriveShaft = zeros(Real, numTS)

    uHist_prp = zeros(Real, numTS,numDOFPerNode)
    FPtfmHist = zeros(Real, numTS,numDOFPerNode)
    FHydroHist = zeros(Real, numTS,numDOFPerNode)
    FMooringHist = zeros(Real, numTS,numDOFPerNode)

    return uHist,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FReactionHist,FTwrBsHist,aziHist,OmegaHist,OmegaDotHist,gbHist,
    gbDotHist,gbDotDotHist,genTorque,genPower,torqueDriveShaft,uHist_prp,
    FPtfmHist,FHydroHist,FMooringHist
end
