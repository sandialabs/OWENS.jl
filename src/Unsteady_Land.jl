mutable struct TopData
    delta_t
    numTS
    numDOFPerNode
    CN2H
    t
    integrator
    integrator_j
    topDispOut
    uHist
    epsilon_x_hist
    epsilon_y_hist
    epsilon_z_hist
    kappa_x_hist
    kappa_y_hist
    kappa_z_hist
    FReactionHist
    FTwrBsHist
    aziHist
    OmegaHist
    OmegaDotHist
    gbHist
    gbDotHist
    gbDotDotHist
    genTorque
    genPower
    torqueDriveShaft
    uHist_prp
    FPtfmHist
    FHydroHist
    FMooringHist
    rbData
    rbDataHist
    u_s
    udot_s
    uddot_s
    u_sm1
    topDispData1
    topDispData2
    topElStrain
    gb_s
    gbDot_s
    gbDotDot_s
    azi_s
    Omega_s
    OmegaDot_s
    genTorque_s
    torqueDriveShaft_s
    topFexternal
    topFexternal_hist
    rotorSpeedForGenStart
    top_rom
    topJointTransformTrans
    u_sRed
    udot_sRed
    uddot_sRed
    topBC
    u_s2
    udot_s2
    uddot_s2
    top_invPhi
    eta_s
    etadot_s
    etaddot_s
    topsideMass
    topsideMOI
    topsideCG
    u_j
    udot_j
    uddot_j
    azi_j
    Omega_j
    OmegaDot_j
    gb_j
    gbDot_j
    gbDotDot_j
    genTorque_j
    torqueDriveShaft_j
    FReactionsm1
    topFReaction_j
end
"""

Unsteady(model,topModel,mesh,el,aero;getLinearizedMatrices=false)

Executable function for transient analysis. Provides the interface of various
    external module with transient structural dynamics analysis capability.

    # Input
    * `inputs::Model`: see ?Model
    * `topModel::FEAModel`: see ?OWENSFEA.FEAModel
    * `mesh::Mesh`: see ?OWENSFEA.Mesh
    * `el::El`: see ?OWENSFEA.El
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
    * `FReactionHist`: Nodal reaction 6dof forces history
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
function Unsteady_Land(inputs;topModel=nothing,topMesh=nothing,topEl=nothing,
    aeroVals=nothing,aeroDOFs=nothing,aero=nothing,deformAero=nothing,
    bottomModel=nothing,bottomMesh=nothing,bottomEl=nothing,bin=nothing,
    getLinearizedMatrices=false, verbosity=0,
    system=nothing,assembly=nothing, #TODO: should we initialize them in here? Unify the interface for ease?
    topElStorage = nothing,bottomElStorage = nothing, u_s = nothing, meshcontrolfunction = nothing,userDefinedGenerator=nothing,turbsimfile=nothing)

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
    FPtfmHist,FHydroHist,FMooringHist,rbData,rbDataHist = allocate_general(inputs,topModel,topMesh,numDOFPerNode,numTS,assembly)

    # Allocate memory for topside
    u_s,udot_s,uddot_s,u_sm1,topDispData1,topDispData2,topElStrain,gb_s,
    gbDot_s,gbDotDot_s,azi_s,Omega_s,OmegaDot_s,genTorque_s,
    torqueDriveShaft_s,topFexternal,topFexternal_hist = allocate_topside(inputs,topMesh,topEl,topModel,numDOFPerNode,u_s,assembly)

    ## Rotor mode initialization
    rotorSpeedForGenStart = initialize_generator!(inputs)

    ## Structural dynamics initialization
    if isnothing(topElStorage)
        topElStorage = OWENSFEA.initialElementCalculations(topModel,topEl,topMesh) #perform initial element calculations for conventional structural dynamics analysis
    end

    top_rom,topJointTransformTrans,u_sRed,udot_sRed,uddot_sRed,
    topBC,u_s2,udot_s2,uddot_s2,top_invPhi,eta_s,etadot_s,
    etaddot_s = initialize_ROM(deepcopy(topElStorage),deepcopy(topModel),deepcopy(topMesh),deepcopy(topEl),deepcopy(u_s),deepcopy(udot_s),deepcopy(uddot_s))

    topDispData1.eta_s = eta_s
    topDispData1.etadot_s = etadot_s
    topDispData1.etaddot_s = etaddot_s
    topDispData2.eta_s = eta_s
    topDispData2.etadot_s = etadot_s
    topDispData2.etaddot_s = etaddot_s


    topsideMass, topsideMOI, topsideCG = OWENSFEA.calculateStructureMassProps(topElStorage)

    # TODO: clean this up
    u_j = 0.0
    udot_j = 0.0
    uddot_j = 0.0
    azi_j = 0.0
    Omega_j = 0.0
    OmegaDot_j = 0.0
    gb_j = 0.0
    gbDot_j = 0.0
    gbDotDot_j = 0.0
    genTorque_j = 0.0
    torqueDriveShaft_j = 0.0
    FReactionsm1 = 0.0
    topFReaction_j = 0.0
    # Package up into data struct
    topdata = TopData(delta_t,numTS,numDOFPerNode,CN2H,t,integrator,integrator_j,topDispOut,
    uHist,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FReactionHist,FTwrBsHist,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,
    genTorque,genPower,torqueDriveShaft,uHist_prp,FPtfmHist,FHydroHist,FMooringHist,rbData,
    rbDataHist,u_s,udot_s,uddot_s,u_sm1,topDispData1,topDispData2,topElStrain,gb_s,gbDot_s,
    gbDotDot_s,azi_s,Omega_s,OmegaDot_s,genTorque_s,torqueDriveShaft_s,topFexternal,topFexternal_hist,
    rotorSpeedForGenStart,top_rom,topJointTransformTrans,u_sRed,udot_sRed,uddot_sRed,topBC,u_s2,udot_s2,
    uddot_s2,top_invPhi,eta_s,etadot_s,etaddot_s,topsideMass,topsideMOI,topsideCG,u_j,udot_j,uddot_j,azi_j,
    Omega_j,OmegaDot_j,gb_j,gbDot_j,gbDotDot_j,genTorque_j,torqueDriveShaft_j,FReactionsm1,topFReaction_j)
    
    topModel.jointTransform, topModel.reducedDOFList = OWENSFEA.createJointTransform(topModel.joint,topMesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints


    topdata.uHist[1,:] = topdata.u_s          #store initial condition
    topdata.aziHist[1] = topdata.azi_s
    topdata.OmegaHist[1] = topdata.Omega_s
    topdata.OmegaDotHist[1] = topdata.OmegaDot_s
    topdata.FReactionsm1 = zeros(length(topdata.u_s))
    topdata.FReactionHist[1,:] = topdata.FReactionsm1
    topdata.topFReaction_j = topdata.FReactionsm1
    # topWeight = [0.0, 0.0, topsideMass*-9.80665, 0.0, 0.0, 0.0] #TODO: propogate gravity, or remove since this isn't used
    topdata.gbHist[1] = gb_s
    topdata.gbDotHist[1] = gbDot_s
    topdata.gbDotDotHist[1] = gbDotDot_s
    topdata.rbDataHist[1,:] = zeros(9)
    topdata.genTorque[1] = genTorque_s
    topdata.torqueDriveShaft[1] = torqueDriveShaft_s


    #..........................................................................
    #                                MAIN LOOP
    #..........................................................................

    ### Iterate for a solution at t+dt
    i=0
    timeconverged = false
    while (i<numTS-1) && timeconverged == false # we compute for the next time step, so the last step of our desired time series is computed in the second to last numTS value
        i += 1
        # println(i)
        ## Print current simulation time to terminal at each .1 second
        if isapprox((t[i]*10)%1,0;atol=5e-2)
            now = round(t[i];digits=1)
            if now == 1
                println("\nSimulation Time: $(now) second of $((numTS-1)*delta_t) seconds")
            else
                println("\nSimulation Time: $(now) seconds of $((numTS-1)*delta_t) seconds")
            end
        end

        ## Check for specified rotor speed at t+dt #TODO: fix this so that it can be probably accounted for in RK4
        if (inputs.turbineStartup == 0)
            inputs.omegaControl = true #TODO: are we setting this back?
            if (inputs.usingRotorSpeedFunction) #use user specified rotor speed profile function
                _,omegaCurrent,_ = getRotorPosSpeedAccelAtTime(t[i],t[i+1],0.0,delta_t)
                topdata.Omega_s = omegaCurrent
            else #use discreteized rotor speed profile function
                omegaCurrent,OmegaDotCurrent,terminateSimulation = omegaSpecCheck(t[i+1],inputs.tocp,inputs.Omegaocp,delta_t)
                if (terminateSimulation)
                    break
                end
                topdata.Omega_s = omegaCurrent
                topdata.OmegaDot_s = OmegaDotCurrent
            end
        else
            omegaCurrent = 0.0
        end

        ## Initialize "j" Gauss-Seidel iteration variables
        topdata.u_j = topdata.u_s
        topdata.udot_j = topdata.udot_s
        topdata.uddot_j = topdata.uddot_s

        topdata.azi_j = topdata.azi_s
        topdata.Omega_j = topdata.Omega_s
        topdata.OmegaDot_j = topdata.OmegaDot_s
        topdata.gb_j = topdata.gb_s
        topdata.gbDot_j = topdata.gbDot_s
        topdata.gbDotDot_j = topdata.gbDotDot_s
        topdata.genTorque_j = topdata.genTorque_s
        topdata.torqueDriveShaft_j = topdata.torqueDriveShaft_s

        #TODO: put these in the model
        TOL = inputs.iteration_parameters.TOL#1e-4  #gauss-seidel iteration tolerance for various modules
        MAXITER = inputs.iteration_parameters.MAXITER#2 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        aziNorm = 1e5
        platNorm = 0.0 #TODO: remove this?
        gbNorm = 0.0 #initialize norms for various module states
        topdata.integrator = topdata.integrator_j #don't compound integrator within the convergence loop.

        if inputs.analysisType=="GX"
            # systemout = deepcopy(system)
            strainGX = zeros(3,length(assembly.elements))
            curvGX = zeros(3,length(assembly.elements))
        end

        newVinf = 0.0 #TODO

        ## Gauss-Seidel predictor-corrector loop
        while ((uNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))
            # println("Iteration $numIterations")

            #------------------
            # GENERATOR MODULE
            #------------------
            topdata.genTorque_j = 0

            if inputs.generatorOn
                    if inputs.useGeneratorFunction
                        specifiedOmega,_,_ = omegaSpecCheck(t[i]+topdata.delta_t,inputs.tocp,inputs.Omegaocp,topdata.delta_t)
                        newVinf = FLOWMath.akima(inputs.tocp_Vinf,inputs.Vinfocp,t[i]) #TODO: ifw sampling of same file as aerodyn
                        if isnothing(userDefinedGenerator)
                            genTorqueHSS0,topdata.integrator_j,controlnamecurrent = internaluserDefinedGenerator(newVinf,t[i],topdata.azi_j,topdata.Omega_j,topdata.OmegaHist[i],topdata.OmegaDot_j,topdata.OmegaDotHist[i],topdata.delta_t,topdata.integrator,specifiedOmega) #;operPhase
                        else
                            if !isnothing(turbsimfile) #&& inputs.AD15On
                                velocity = OWENSOpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)],t[i])
                                newVinf = velocity[1]
                            end
                            genTorqueHSS0 = userDefinedGenerator(t[i],topdata.Omega_j*60,newVinf) #;operPhase
                        end
                    else
                        genTorqueHSS0 = simpleGenerator(inputs,topdata.Omega_j)
                    end
                
                #should eventually account for Omega = gbDot*gearRatio here...
                topdata.genTorque_j = genTorqueHSS0*inputs.gearRatio*inputs.gearBoxEfficiency #calculate generator torque on LSS side
                
                #         genTorqueAppliedToTurbineRotor0 = -genTorque0
                #         genTorqueAppliedToPlatform0 = genTorqueHSS0
            else
                if !isnothing(turbsimfile) && verbosity >=1#&& inputs.AD15On
                    velocity = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)],t[i])
                    newVinf = velocity[1]
                end
            end
            
            #-------------------
            # DRIVETRAIN MODULE
            #-------------------
            topdata.torqueDriveShaft_j = topdata.genTorque_j
            gb_jLast = topdata.gb_j
            if (!inputs.omegaControl)
                if (inputs.driveTrainOn)
                    topdata.torqueDriveShaft_j = calculateDriveShaftReactionTorque(inputs.driveShaftProps,
                    topdata.azi_j,topdata.gb_j,topdata.Omega_j*2*pi,topdata.gbDot_j*2*pi)

                    topdata.gb_j,topdata.gbDot_j,topdata.gbDotDot_j = updateRotorRotation(inputs.JgearBox,0,0,
                    -topdata.genTorque_j,topdata.torqueDriveShaft_j,topdata.gb_s,topdata.gbDot_s,topdata.gbDotDot_s,topdata.delta_t)
                else
                    topdata.gb_j = topdata.azi_j
                    topdata.gbDot_j = topdata.Omega_j
                    topdata.gbDotDot_j = topdata.OmegaDot_j
                end
            else
                topdata.gb_j = topdata.azi_j
                topdata.gbDot_j = omegaCurrent*2*pi
                topdata.gbDotDot_j = 0
            end

            # Update rotor speed
            azi_jLast = topdata.azi_j
            if inputs.omegaControl
                if (inputs.usingRotorSpeedFunction)
                    topdata.azi_j,topdata.Omega_j,topdata.OmegaDot_j = getRotorPosSpeedAccelAtTime(t[i],t[i+1],azi_s,delta_t)
                else
                    topdata.Omega_j = topdata.Omega_s
                    topdata.OmegaDot_j = topdata.OmegaDot_s
                    topdata.azi_j = topdata.azi_s + topdata.Omega_j*topdata.delta_t*2*pi
                end
            elseif !inputs.omegaControl
                Crotor = 0
                Krotor = 0
                topdata.azi_j,topdata.Omega_j,topdata.OmegaDot_j = updateRotorRotation(topdata.topsideMOI[3,3],Crotor,Krotor,
                -topdata.topFReaction_j[6],-topdata.torqueDriveShaft_j,
                topdata.azi_s,topdata.Omega_s,topdata.OmegaDot_s,topdata.delta_t)
            else
                error("omega control option not correctly specified")
            end

            #---------------------
            # AERODYNAMICS MODULE
            #---------------------
            # Calculate new aerodynamic loading

            # Update reference frame transformation and convert aerodynamic loads to hub reference frame

            topdata.CN2H = calcHubRotMat(zeros(3), azi_j)
            CN2H_no_azi = calcHubRotMat(zeros(3), 0.0)
            rbData = zeros(9)

            CH2N = LinearAlgebra.transpose(CN2H)


            #################################################################
            if !isnothing(aero)
                if inputs.aeroLoadsOn > 0 #0 off, 1 one way, 1.5 one way with deformation from last timestep, 2 two way
                    runaero = true
                    if (inputs.aeroLoadsOn==1 || inputs.aeroLoadsOn==1.5) && numIterations!=1
                        runaero = false
                    end
                    
                    if runaero
                        if inputs.AD15On
                            aeroVals,aeroDOFs = run_aero_with_deformAD15(aero,deformAero,[topMesh],[topEl],[topdata],[inputs],t[i])
                            aeroVals = aeroVals[1]
                            aeroDOFs = aeroDOFs[1]
                        else
                            aeroVals,aeroDOFs = run_aero_with_deform(aero,deformAero,topMesh,topEl,topdata.u_j,inputs,numIterations,t[i],topdata.azi_j,topdata.Omega_j)
                        end
                    end
                end
            end

            #################################################################
            if isnan(maximum(aeroVals))
                @warn "Nan detected in aero forces"
            end
            if inputs.aeroLoadsOn > 0
                if isnothing(aeroVals)
                    error("aeroVals must be specified if OWENS.Inputs.aeroLoadsOn")
                elseif isnothing(aeroDOFs)
                    error("aeroDOFs must be specified if OWENS.Inputs.aeroLoadsOn")
                end

                if inputs.AD15On
                    # AD15 is in global frame, so no frame conversion???
                    topdata.topFexternal = aeroVals
                else
                    if length(size(aeroVals))==1 || size(aeroVals)[2]==1 #i.e. the standard aero force input as a long array
                        # Fill in forces and dofs if they were specified not in full arrays TODO: make this more efficient
                        full_aeroVals = zeros(topMesh.numNodes*6)
                        for i_idx = 1:length(aeroDOFs)
                            full_aeroVals[Int(aeroDOFs[i_idx])] = aeroVals[i_idx]
                        end
                        aeroDOFs = collect(1:topMesh.numNodes*6)
                        for iter_i = 1:floor(Int,length(full_aeroVals)/6)
                            topdata.topFexternal[6*(iter_i-1)+1:6*(iter_i-1)+6] = frame_convert(full_aeroVals[6*(iter_i-1)+1:6*(iter_i-1)+6], CN2H_no_azi)
                        end
                    else # the other aero input as a 2D array
                        topdata.topFexternal = frame_convert(aeroVals[i+1,:], CN2H)
                    end
                end
            else
                topdata.topFexternal = zeros(numDOFPerNode)
                aeroDOFs = copy(topdata.topFexternal).*0.0
            end
            aeroVals = topdata.topFexternal

            if meshcontrolfunction !== nothing
                # add to the loads based on the inputs, TODO: CN2H
                meshforces, meshdofs, timeconverged = meshcontrolfunction(topMesh,u_j,t[i])
                for idx_main in aeroDOFs
                    for (idx,meshdof_idx) in enumerate(meshdofs)
                        if idx_main == meshdof_idx
                            topdata.topFexternal[idx_main] += meshforces[idx]
                        end
                    end
                end
            end

            #------------------------------------
            # TOPSIDE STRUCTURAL MODULE
            #------------------------------------
            # println(Float64.(rbData))
            if inputs.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j = OWENSFEA.structuralDynamicsTransientROM(topModel,topMesh,topEl,topdata.topDispData1,topdata.Omega_s,topdata.OmegaDot_s,t[i+1],topdata.delta_t,topElStorage,topdata.top_rom,topdata.topFexternal,Int.(aeroDOFs),topdata.CN2H,topdata.rbData)
            elseif inputs.analysisType=="GX"                                                                                
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j,systemout  = structuralDynamicsTransientGX(topModel,topMesh,topdata.topFexternal,Int.(aeroDOFs),system,assembly,t,topdata.Omega_j,topdata.OmegaDot_j,topdata.delta_t,numIterations,i,strainGX,curvGX)
            else # evalulate structural dynamics using conventional representation
                topdata.topElStrain, topdata.topDispOut, topdata.topFReaction_j = OWENSFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topdata.topDispData1,topdata.Omega_s,topdata.OmegaDot_s,t[i+1],topdata.delta_t,topElStorage,topdata.topFexternal,Int.(aeroDOFs),topdata.CN2H,topdata.rbData;predef = topModel.nlParams.predef)
            end

            u_jLast = copy(topdata.u_j)
            topdata.u_j = topdata.topDispOut.displ_sp1
            topdata.udot_j = topdata.topDispOut.displdot_sp1
            topdata.uddot_j = topdata.topDispOut.displddot_sp1

            ## calculate norms
            uNorm = LinearAlgebra.norm(topdata.u_j .- u_jLast)/LinearAlgebra.norm(topdata.u_j)            #structural dynamics displacement iteration norm
            aziNorm = LinearAlgebra.norm(topdata.azi_j .- azi_jLast)/LinearAlgebra.norm(topdata.azi_j)  #rotor azimuth iteration norm
            if inputs.generatorOn
                gbNorm = LinearAlgebra.norm(topdata.gb_j .- gb_jLast)/LinearAlgebra.norm(topdata.gb_j) #gearbox states iteration norm if it is off, the norm will be zero
            end

            numIterations = numIterations + 1

            if inputs.analysisType=="GX" && (!(uNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                system = deepcopy(systemout)
            end

            # Strain stiffening, save at the end of the simulation, at the last while loop iteration, mutates elStorage
            if (i==numTS-1 || timeconverged == true) && inputs.analysisType=="TNB" && topModel.nlParams.predef=="update" && (!(uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) || (numIterations >= MAXITER))
                OWENSFEA.structuralDynamicsTransient(topModel,topMesh,topEl,topdata.topDispData2,topdata.Omega_s,topdata.OmegaDot_s,t[i+1],topdata.delta_t,topdata.topElStorage,topdata.topFexternal,Int.(aeroDOFs),topdata.CN2H,topdata.rbData;predef = topModel.nlParams.predef)
            end

            if verbosity>4
                println("$(numIterations) uNorm: $(uNorm) aziNorm: $(aziNorm) gbNorm: $(gbNorm)")
            end

            if numIterations==MAXITER
                if inputs.iteration_parameters.iterwarnings
                    @warn "Maximum Iterations Met Breaking Iteration Loop"
                end
                break
            end

        end #end iteration while loop

        if verbosity >=1
            println("Gen Torque: $(topdata.genTorque_j)")
            println("RPM: $(topdata.Omega_j*60)")
            println("Vinf: $(newVinf)")
            # velocitymid = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)/2],t[i])
            # velocityquarter = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,maximum(topMesh.z)/4],t[i])
            # println("Velocity mid: $(velocitymid[1])")
            # println("Velocity quarter: $(velocityquarter[1])")
        end

        if inputs.analysisType=="ROM"
            eta_j = topdata.topDispOut.eta_sp1 #eta_j
            etadot_j = topdata.topDispOut.etadot_sp1 #etadot_j
            etaddot_j = topdata.topDispOut.etaddot_sp1 #etaddot_j
            topdata.topDispData1 = OWENSFEA.DispData(topdata.u_j, topdata.udot_j, topdata.uddot_j, topdata.u_sm1, eta_j, etadot_j, etaddot_j)
        else
            topdata.topDispData1 = OWENSFEA.DispData(topdata.u_j, topdata.udot_j, topdata.uddot_j, topdata.u_sm1)
        end

        if isnan(maximum(u_j))
            @warn "Nan detected in displacements"
            break
        end   

        ## update timestepping variables and other states, store in history arrays
        ## calculate converged generator torque/power

        genPowerPlot = topdata.genTorque_s*(gbDot_j*2*pi)*inputs.gearRatio

        topdata.u_sm1 = copy(topdata.u_s)
        topdata.u_s = topdata.u_j
        topdata.udot_s = topdata.udot_j
        topdata.uddot_s = topdata.uddot_j

        topdata.azi_s = topdata.azi_j
        topdata.Omega_s = topdata.Omega_j
        topdata.OmegaDot_s = topdata.OmegaDot_j

        topdata.genTorque_s = topdata.genTorque_j
        topdata.torqueDriveShaft_s = topdata.torqueDriveShaft_j

        topdata.gb_s = topdata.gb_j
        topdata.gbDot_s = topdata.gbDot_j
        topdata.gbDotDot_s = topdata.gbDotDot_j

        topdata.uHist[i+1,:] = topdata.u_s
        topdata.FReactionHist[i+1,:] = topdata.topFReaction_j
        topdata.topFexternal_hist[i+1,1:length(topdata.topFexternal)] = topdata.topFexternal
        # FTwrBsHist[i+1,:] = -topFReaction_j  + topFWeight_j

        topdata.aziHist[i+1] = topdata.azi_s
        topdata.OmegaHist[i+1] = topdata.Omega_s
        topdata.OmegaDotHist[i+1] = topdata.OmegaDot_s

        topdata.gbHist[i+1] = topdata.gb_s
        topdata.gbDotHist[i+1] = topdata.gbDot_s
        topdata.gbDotDotHist[i+1] = topdata.gbDotDot_s

        #genTorque[i+1] = genTorque_s
        topdata.genTorque[i+1] = topdata.genTorque_s
        topdata.genPower[i+1] = genPowerPlot
        topdata.torqueDriveShaft[i+1] = topdata.torqueDriveShaft_s

        if inputs.analysisType=="GX"
            strainGX = topdata.topElStrain[1]
            curvGX = topdata.topElStrain[2]
            for ii = 1:length(strainGX[1,:])
                topdata.epsilon_x_hist[:,ii,i] .= strainGX[1,ii]
                topdata.kappa_y_hist[:,ii,i] .= curvGX[2,ii]
                topdata.kappa_z_hist[:,ii,i] .= curvGX[3,ii]
                topdata.epsilon_z_hist[:,ii,i] .= strainGX[3,ii]
                topdata.kappa_x_hist[:,ii,i] .= curvGX[1,ii]
                topdata.epsilon_y_hist[:,ii,i] .= strainGX[2,ii]
            end
        else
            for ii = 1:length(topdata.topElStrain)
                topdata.epsilon_x_hist[:,ii,i] = topdata.topElStrain[ii].epsilon_x
                topdata.kappa_y_hist[:,ii,i] = topdata.topElStrain[ii].kappa_y
                topdata.kappa_z_hist[:,ii,i] = topdata.topElStrain[ii].kappa_z
                topdata.epsilon_z_hist[:,ii,i] = topdata.topElStrain[ii].epsilon_z
                topdata.kappa_x_hist[:,ii,i] = topdata.topElStrain[ii].kappa_x
                topdata.epsilon_y_hist[:,ii,i] = topdata.topElStrain[ii].epsilon_y
            end
        end
        ## check rotor speed for generator operation
        if topdata.Omega_s >= topdata.rotorSpeedForGenStart
            inputs.generatorOn = true
        else
            inputs.generatorOn = false
        end

    end #end timestep loop

    println("Simulation Complete.")

    outputData(topMesh,inputs,t[1:i],topdata.aziHist[1:i],topdata.OmegaHist[1:i],topdata.OmegaDotHist[1:i],topdata.gbHist[1:i],topdata.gbDotHist[1:i],topdata.gbDotDotHist[1:i],
    topdata.FReactionHist[1:i,:],topdata.genTorque[1:i],topdata.genPower[1:i],topdata.torqueDriveShaft[1:i],topdata.uHist[1:i,:],topdata.uHist_prp[1:i,:],
    topdata.epsilon_x_hist[:,:,1:i],topdata.epsilon_y_hist[:,:,1:i],topdata.epsilon_z_hist[:,:,1:i],topdata.kappa_x_hist[:,:,1:i],topdata.kappa_y_hist[:,:,1:i],
    topdata.kappa_z_hist[:,:,1:i])
    
    return t[1:i], topdata.aziHist[1:i],topdata.OmegaHist[1:i],topdata.OmegaDotHist[1:i],topdata.gbHist[1:i],topdata.gbDotHist[1:i],topdata.gbDotDotHist[1:i],
    topdata.FReactionHist[1:i,:],topdata.FTwrBsHist[1:i,:],topdata.genTorque[1:i],topdata.genPower[1:i],topdata.torqueDriveShaft[1:i],topdata.uHist[1:i,:],
    topdata.uHist_prp[1:i,:],topdata.epsilon_x_hist[:,:,1:i],topdata.epsilon_y_hist[:,:,1:i],topdata.epsilon_z_hist[:,:,1:i],topdata.kappa_x_hist[:,:,1:i],
    topdata.kappa_y_hist[:,:,1:i],topdata.kappa_z_hist[:,:,1:i],topdata.FPtfmHist[1:i,:],topdata.FHydroHist[1:i,:],topdata.FMooringHist[1:i,:],
    topdata.topFexternal_hist[1:i,:],topdata.rbDataHist[1:i,:]
end
