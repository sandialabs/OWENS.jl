"""
Internal, gets specified rotor speed at time
"""
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

"""
Internal, generator definintion
"""
function userDefinedGenerator(input)
    return 0.0 #this is what was in the original file...
end

"""

    getRotorPosSpeedAccelAtTime(t0,time,aziInit)

Uses the user defined function rotorSpeedProfile() to get
the azimuth, speed, and acceleration of the rotor.

#Input
* `t0`      time at which azimuth integration is beginning
* `time`    current time that position, velocity, and acceleration are being requested
* `aziInit` initial rotor azimuth angle integration will begin at

#Output
* `rotorAzimuth` azimuth position of rotor (rad) at time
* `rotorSpeed`   rotor speed (Hz) at time
* `rotorAcceleration` rotor acceleration (Hz/s) at time
"""
function getRotorPosSpeedAccelAtTime(t0,time,aziInit,delta_t)

    rotorSpeed = userDefinedRotorSpeedProfile(time) #get rotor speed at time

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

"""
Internal, simple trapezoidal rule integration
"""
function trapezoidalRule(aziInit,rotorSpeedStart,rotorSpeedEnd,dt)
    return (aziInit + 0.5*dt*(rotorSpeedStart+rotorSpeedEnd)/(2*pi))
end

# 3x3 transformation matrix between coordinate systems (with small angle assumption)
function transMat(theta1, theta2, theta3)
    theta11      = theta1 * theta1
    theta22      = theta2 * theta2
    theta33      = theta3 * theta3

    sqrdSum      = theta11 + theta22 + theta33
    sqrt1sqrdSum = sqrt( 1.0 + sqrdSum )
    comDenom     = sqrdSum * sqrt1sqrdSum

    theta12S     = theta1 * theta2 * ( sqrt1sqrdSum - 1.0 )
    theta13S     = theta1 * theta3 * ( sqrt1sqrdSum - 1.0 )
    theta23S     = theta2 * theta3 * ( sqrt1sqrdSum - 1.0 )


    # Define the transformation matrix:
    transMat = Array{Float64}(undef, 3,3)
    if comDenom == 0.0  # All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity

        transMat = LinearAlgebra.I(3)

    else  # At least one angle is nonzero

        transMat[1,1] = ( theta11*sqrt1sqrdSum + theta22              + theta33              ) / comDenom
        transMat[2,2] = ( theta11              + theta22*sqrt1sqrdSum + theta33              ) / comDenom
        transMat[3,3] = ( theta11              + theta22              + theta33*sqrt1sqrdSum ) / comDenom
        transMat[1,2] = (  theta3*sqrdSum + theta12S ) / comDenom
        transMat[2,1] = ( -theta3*sqrdSum + theta12S ) / comDenom
        transMat[1,3] = ( -theta2*sqrdSum + theta13S ) / comDenom
        transMat[3,1] = (  theta2*sqrdSum + theta13S ) / comDenom
        transMat[2,3] = (  theta1*sqrdSum + theta23S ) / comDenom
        transMat[3,2] = ( -theta1*sqrdSum + theta23S ) / comDenom

    end

    return transMat

end
"""
Internal, unused, userDefinedRotorSpeedProfile
"""
function userDefinedRotorSpeedProfile(time)
    return 0.5 #this is what was originally in the file...
end

"""

    externalForcing(time,timeArray,ForceValHist,ForceDof)

Internal, linear time interpolation on the input forces (ForceValHist) for each specified ForceDof

This function specifies external forcing for a transient analysis.
Fexternal is a vector of loads and Fdof is a corresponding vector of
degrees of freedom the concentrated loads in Fexternal correspond to.
The input time allows for arbitrary time varying loads
The global degree of freedom number corresponding with the local degree
of freedom of a node may be calculated by:
globalDOFNumber = (nodeNumber-1)*6 + localDOFnumber
The localDOFnumber may range from 1 to 6 such that 1 corresponds to a
force in "x direction" of the co-rotating hub frame. 2 and 3
corresponds to a force in the "y" and "z directions" respectively. 4,
5, and 6 correspond to a moment about the "x", "y", and "z" directions
respectively.


#Input
* `time`: Current time
* `timeArray`: time associated with ForceValHist
* `ForceValHist`: Forces for each time for each Dof
* `ForceDof`: Dofs within ForceValHist

#Output
* `Fexternal`: vector of external loads (forces/moments)
* `Fdof`:      vector of corresponding DOF numbers to apply loads to
"""
function externalForcing(time,timeArray,ForceValHist,ForceDof)

    Fexternal = zeros(length(ForceDof))

    for i = 1:length(ForceDof)
        Fexternal[i] = FLOWMath.linear(timeArray,ForceValHist[i,:],time)
    end
    Fdof = ForceDof

    return Fexternal, Fdof
end

"""

    calculateDriveShaftReactionTorque(driveShaftProps,thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)

Internal, calculates reaction torque of driveshaft

#Input
* `driveShaftProps`:    object containing driveshaft properties
* `thetaRotor`:         azimuth position of rotor/rotor shaft (rad)
* `thetaGB`:            azimuth position of gearbox shaft (rad)
* `thetaDotRotor`:      angular velocity of rotor/rotor shaft (rad/s)
* `thetaDotGB`:         angular velocity of gearbox shaft (rad/s)

#Output
* `torque`: reaction torque of drive shaft
"""
function calculateDriveShaftReactionTorque(driveShaftProps,thetaRotor,thetaGB,thetaDotRotor,thetaDotGB)

    k = driveShaftProps.k  #drive shaft stiffness
    c = driveShaftProps.c  #drive shaft damping

    return k*(thetaRotor-thetaGB) + c*(thetaDotRotor-thetaDotGB)

end

"""
updateRotorRotation updates rotor rotation

    updateRotorRotation(Irotor,Crotor,Krotor,shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s,delta_t)

Internal, updates the rotor rotation given rotor properties and external torques

#Input
* `Irotor`:      rotor inertia
* `Crotor`:      arbitrary rotor damping
* `Krotor`:      arbitrary rotor stiffness
* `shaftTorque`: torque from external forces on rotor
* `genTorque`:   torque from generator
* `azi_s`:       rotor azimuth (rad) at beginning of time step
* `Omega_s`:     rotor speed (Hz) at beginning of time step
* `OmegaDot_s`:  rotor acceleration (Hz/s) at beginning of time step
* `delta_t`:     time step

#Output
* `azi_sp1`:       rotor azimuth (rad) at end of time step
* `Omega_sp1`:     rotor speed (Hz/s) at end of time step
* `OmegaDot_sp1`:  rotor acceleration (Hz/s) at end of time step

"""
function updateRotorRotation(Irotor,Crotor,Krotor,shaftTorque,genTorque,azi_s,Omega_s,OmegaDot_s, delta_t)

    Frotor = shaftTorque + genTorque #calculate effective torque on rotor
    Omega_s = Omega_s*2*pi #conversion form Hz to rad/s, etc.
    OmegaDot_s = OmegaDot_s*2*pi
    azi_sp1,Omega_sp1,OmegaDot_sp1 = timeIntegrateSubSystem(Irotor,Krotor,Crotor,Frotor, #time integrate using Newmark-Beta
    delta_t,azi_s,Omega_s,OmegaDot_s)

    Omega_sp1 = Omega_sp1/(2*pi) #convert to Hz, etc.
    OmegaDot_sp1 = OmegaDot_sp1/(2*pi)

    return azi_sp1,Omega_sp1,OmegaDot_sp1

end

"""

    timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)

Internal, performs integration of a system using the Newmark-Beta method(constant-average acceleration sceheme).

#Input
* `M`:       system mass matrix
* `K`:       system sttiffness matrix
* `C`:       system damping matrix
* `F`:       system force vector
* `delta_t`: time step
* `u`:       displacement at beginning of time step
* `udot`:    velocity at beginning of time step
* `uddot`:   acceleration at beginning of time step


#Output
* `unp1`:       displacement at end of time step
* `udotnp1`:    velocity at end of time step
* `uddotnp1`:   acceleration at end of time step

"""
function timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)

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


function calc_hydro_residual(new_accels, new_hydro_frcs, md_frc, u, frc_multiplier)

    new_frcs = md_frc + new_hydro_frcs

    u_resid = Vector{Float32}(undef, length(new_frcs) + length(new_accels))
    u_resid[1:length(new_frcs)] = u[1:length(new_frcs)] - new_frcs/frc_multiplier
    u_resid[length(new_frcs)+1:end] = u[length(new_frcs)+1:end] - new_accels

    return u_resid
end


function OWENS_HD_Coupling(time, dt, calcJacobian, jac, numDOFPerNode, ptfm_dofs, frc_hydro_in,
    frc_mooring_in, dispData, feamodel, mesh, el, elStorage, Omega, OmegaDot,
    other_Fexternal, other_Fdof, CN2H)

    # !!! Make sure structuralDynamicsTransient and MD_CalcOutput have run before this so you can send dispOut and frc_mooring_in here !!!
    # 
    # allocate new variables
    frc_hydro2 = zero(frc_hydro_in)
    if calcJacobian
        frc_hydro_perturb_hd = zero(frc_hydro_in)
    end
    frc_hydro_out = zero(frc_hydro_in) # TODO: all of the frc_hydro's may be able to be combined for efficiency
    out_vals = zeros(numDOFPerNode+1)
    u = Vector{Float32}(undef, numDOFPerNode*2)
    

    frc_multiplier = 1e6

    dispIn = deepcopy(dispData) # TODO: does this need to be copied, or can we just use the original dispData?
    u_ptfm = dispIn.displ_s[ptfm_dofs]
    udot_ptfm = dispIn.displdot_s[ptfm_dofs]
    uddot_ptfm = dispIn.displddot_s[ptfm_dofs]
    
    u[1:numDOFPerNode] = frc_hydro_in / frc_multiplier
    u[numDOFPerNode+1:numDOFPerNode*2] = uddot_ptfm

    # Calculate outputs at the current time, based on inputs at the current time
    total_Fexternal = [other_Fexternal; frc_hydro_in + frc_mooring_in]
    total_Fdof = [other_Fdof; ptfm_dofs]
    _ ,dispOut2, _ = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal,Int.(total_Fdof),CN2H,zeros(9))
    frc_hydro2[:], _ = VAWTHydro.HD_CalcOutput(time, u_ptfm, udot_ptfm, uddot_ptfm, frc_hydro2, out_vals)

   # Calculate the residual
   uddot_ptfm2 = dispOut2.displddot_sp1[ptfm_dofs]
   residual = calc_hydro_residual(uddot_ptfm2, frc_hydro2, frc_mooring_in, u, frc_multiplier)

   # Calculate the Jacobian. Since there's no single function associated with the motions/forces, we will manually
   # perturb each degree of freedom one at a time and recalculate the outputs so we can get a full gradient.
    if calcJacobian

        for dof = collect(1:numDOFPerNode) # forces
            frc_hydro_perturb_owens = deepcopy(frc_hydro_in)
            u_perturb = copy(u)
            frc_hydro_perturb_owens[dof] += 1E6
            u_perturb[dof] += 1
            total_Fexternal_perturb = [other_Fexternal; frc_hydro_perturb_owens + frc_mooring_in]
            _, dispOut_perturb_owens, _ = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal_perturb,Int.(total_Fdof),CN2H,zeros(9))
            uddot_ptfm_perturb = dispOut_perturb_owens.displddot_sp1[ptfm_dofs]
            residual_perturb = calc_hydro_residual(uddot_ptfm_perturb, frc_hydro2, frc_mooring_in, u_perturb, frc_multiplier)
            jac[:,dof] = residual_perturb - residual
        end

        for dof = collect(1:numDOFPerNode) # accelerations
            uddot_ptfm_perturb_hd = deepcopy(uddot_ptfm)
            u_perturb = copy(u)
            uddot_ptfm_perturb_hd[dof] += 1
            u_perturb[dof+numDOFPerNode] += 1
            frc_hydro_perturb_hd[:], _ = VAWTHydro.HD_CalcOutput(time, u_ptfm, udot_ptfm, uddot_ptfm_perturb_hd, frc_hydro_perturb_hd, out_vals)
            residual_perturb = calc_hydro_residual(uddot_ptfm2, frc_hydro_perturb_hd, frc_mooring_in, u_perturb, frc_multiplier)
            jac[:,dof+numDOFPerNode] = residual_perturb - residual
        end

   end  # if calcJacobian

   # Solve for delta_u: jac*delta_u = -residual
   delta_u = jac\-residual
   jac_out = jac

   # Update inputs
   frc_hydro_rev = frc_hydro_in + delta_u[1:numDOFPerNode]*frc_multiplier
   total_Fexternal = [other_Fexternal; frc_hydro_rev + frc_mooring_in]
   dispData_rev = deepcopy(dispIn)
   dispData_rev.displddot_s[ptfm_dofs] = dispIn.displddot_s[ptfm_dofs] + delta_u[numDOFPerNode+1:numDOFPerNode*2]
   uddot_ptfm = dispData_rev.displddot_s[ptfm_dofs]

   # Rerun OWENS and HydroDyn with updated inputs
   elStrain_out,dispOut,FReaction_out = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal,Int.(total_Fdof),CN2H,zeros(9)) # TODO: should we use dispData_rev instead of dispIn here?
   frc_hydro_out[:], out_vals[:] = VAWTHydro.HD_CalcOutput(time, u_ptfm, udot_ptfm, uddot_ptfm, frc_hydro_out, out_vals)

   return elStrain_out, dispOut, FReaction_out, frc_hydro_out, out_vals, jac_out

end


"""

    Unsteady(model,feamodel,mesh,el,aero;getLinearizedMatrices=false)

Executable function for transient analysis. Provides the interface of various
external module with transient structural dynamics analysis capability.

# Input
* `model::Model`: see ?Model
* `feamodel::FEAModel`: see ?GyricFEA.FEAModel
* `mesh::Mesh`: see ?GyricFEA.Mesh
* `el::El`: see ?GyricFEA.El
* `aero::function`: Fexternal, Fdof = aero(t) where Fexternal is the force on each affected mesh dof and Fdof is the corresponding DOFs affected
* `getLinearizedMatrices::Bool`: Flag to save the linearized matrices


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
* `eps_xx_0_hist`: strain history for eps_xx_0 for each dof
* `eps_xx_z_hist`: strain history for eps_xx_z for each dof
* `eps_xx_y_hist`: strain history for eps_xx_y for each dof
* `gam_xz_0_hist`: strain history for gam_xz_0 for each dof
* `gam_xz_y_hist`: strain history for gam_xz_y for each dof
* `gam_xy_0_hist`: strain history for gam_xy_0 for each dof
* `gam_xy_z_hist`: strain history for gam_xy_z for each dof
"""
function Unsteady(model,feamodel,mesh,el,bin,aero,deformAero;getLinearizedMatrices=false)

    # Declare Variable Type, are set later
    udot_j = 0.0
    uddot_j = 0.0
    torqueDriveShaft_j = 0.0

    elStrain = fill(GyricFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)

    dispOut = GyricFEA.DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492)) # I don't think this needs to be initialized here
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

    ## state initialization
    numDOFPerNode = 6
    totalNumDof = mesh.numNodes*numDOFPerNode

    #......... specify initial conditions .......................
    u_s = zeros(totalNumDof)
    u_s = GyricFEA.setInitialConditions(feamodel.initCond,u_s,numDOFPerNode)
    u_sm1 = copy(u_s)
    udot_s = zero(u_s)
    uddot_s = zero(u_s)
    dispData = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1)
    #............................................................

    numTS = Int(model.numTS)       #define number of time steps
    delta_t = model.delta_t   #define time step size
    uHist = zeros(numTS, length(u_s))
    uHist[1,:] = u_s          #store initial condition

    t = range(0, length=numTS, step=delta_t)

    FReactionHist = zeros(numTS,6)

    eps_xx_0_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_z_hist = zeros(4,mesh.numEl,numTS)
    eps_xx_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xz_y_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_0_hist = zeros(4,mesh.numEl,numTS)
    gam_xy_z_hist = zeros(4,mesh.numEl,numTS)

    aziHist = zeros(numTS)
    OmegaHist = zeros(numTS)
    OmegaDotHist = zeros(numTS)
    gbHist = zeros(numTS)
    gbDotHist = zeros(numTS)
    gbDotDotHist = zeros(numTS)
    genTorque = zeros(numTS)
    genPower = zeros(numTS)
    torqueDriveShaft = zeros(numTS)

    # initialize various states and variables
    gb_s = 0
    gbDot_s = 0
    gbDotDot_s = 0
    azi_s = 0
    Omega_s = OmegaInitial
    OmegaDot_s = 0
    genTorque_s = 0
    torqueDriveShaft_s = 0

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

    ## structural dynamics initialization
    #..........................................................................
    if model.analysisType=="ROM" #initialize reduced order model
        #calculate constrained dof vector
        #TODO: This is already done, remove this redundant code
        # numDofPerNode = 6
        # isConstrained = zeros(totalNumDof)
        # constDof = (feamodel.BC.pBC[:,1]-1)*numDofPerNode + feamodel.BC.pBC[:,2]
        # index = 1
        # for i=1:mesh.numNodes
        #     for j=1:numDofPerNode
        #         if ((i-1)*numDofPerNode + j in constDof)
        #             isConstrained[index] = 1
        #         end
        #         index = index + 1
        #     end
        # end
        # feamodel.BC.isConstrained = isConstrained

        rom,elStorage = GyricFEA.reducedOrderModel(feamodel,mesh,el,u_s) #construct reduced order model

        #set up inital values in modal space
        jointTransformTrans = feamodel.jointTransform' #'
        u_sRed = jointTransformTrans*u_s
        udot_sRed = jointTransformTrans*udot_s
        uddot_sRed = jointTransformTrans*uddot_s

        BC = feamodel.BC
        u_s2 = GyricFEA.applyBCModalVec(u_sRed,BC.numpBC,BC.map)
        udot_s2 = GyricFEA.applyBCModalVec(udot_sRed,BC.numpBC,BC.map)
        uddot_s2 = GyricFEA.applyBCModalVec(uddot_sRed,BC.numpBC,BC.map)

        invPhi = rom.invPhi

        eta_s     = invPhi*u_s2
        etadot_s  = invPhi*udot_s2
        etaddot_s = invPhi*uddot_s2

    else
        elStorage = GyricFEA.initialElementCalculations(feamodel,el,mesh) #perform initial element calculations for conventional structural dynamics analysis
    end

    #calculate structural/platform moi
    _,structureMOI,_ = GyricFEA.calculateStructureMassProps(elStorage)
    #..........................................................................

    feamodel.jointTransform, feamodel.reducedDOFList = GyricFEA.createJointTransform(feamodel.joint,mesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints

    # Initialize hydro libraries and working variables
    if model.hydroOn
        jac = Array{Float64}(undef, numDOFPerNode*2, numDOFPerNode*2)
        numMooringLines = 3

        if model.outFilename == "none"
            hd_outFilename = "hydrodyn_temp.out"
        else
            hd_outFilename = model.outFilename
        end

        ptfm_dofs = numDOFPerNode+1:numDOFPerNode*2 # platform BCs are modeled at the second node (bottom node is fixed)
        u_s_ptfm = Vector(u_s[ptfm_dofs])
        udot_s_ptfm = Vector(udot_s[ptfm_dofs])
        uddot_s_ptfm = Vector(uddot_s[ptfm_dofs])

        uHist_ptfm = zeros(numTS,length(u_s_ptfm))
        uHist_ptfm[1,:] = u_s_ptfm

        VAWTHydro.HD_Init(bin.hydrodynLibPath, hd_outFilename, hd_input_file=model.hd_input_file, PotFile=model.potflowfile, t_initial=t[1], dt=delta_t, t_max=t[1]+(numTS-1)*delta_t)
        VAWTHydro.MD_Init(bin.moordynLibPath, md_input_file=model.md_input_file, init_ptfm_pos=u_s_ptfm, interp_order=model.interpOrder)
        
        frc_hydro_n = zeros(Float32, numDOFPerNode)
        frc_hydro_h = zeros(Float32, numDOFPerNode)
        frc_mooring_n = zeros(Float32, numDOFPerNode)
        frc_mooring_h = zeros(Float32, numDOFPerNode)
        out_vals = zeros(Float32, numDOFPerNode+1) # Rigid body displacement in 6DOF + wave elevation
        mooring_tensions = zeros(Float32, numMooringLines*2) # Fairlead + anchor tension for each line
        
        # calculate initial HydroDyn/MoorDyn states
        # TODO: transform these to the inertial frame, since u_s is in the hub frame
        #       (u_s and azi_s are initialized to zero, so it doesn't matter right now, but it will if we add nonzero initial conditions)
        frc_mooring_n[:], mooring_tensions[:] = VAWTHydro.MD_CalcOutput(t[1], u_s_ptfm, udot_s_ptfm, uddot_s_ptfm, frc_mooring_n, mooring_tensions)
        elStrain, dispOut, FReaction, frc_hydro_n, out_vals, jac = OWENS_HD_Coupling(t[1], delta_t, true, jac, numDOFPerNode, ptfm_dofs, frc_hydro_n,
                                                                                    frc_mooring_n, dispData, feamodel, mesh, el, elStorage, Omega_s, OmegaDot_s,
                                                                                    0.0, 1, zeros(3,3))
        

        # Spd = wave.resource.jonswap_spectrum(f=model.plat_model.hydro.freq, Tp=6, Hs=1)
    end

    ptfm_roll = 0 #placeholder for numerical stability while debugging #TODO reverse azi_j to get platform motions in inertial frame
    ptfm_pitch = 0
    ptfm_yaw = 0

    ## Main Loop - iterate for a solution at each time step, i
    for i=1:numTS-1 # we compute for the next time step, so the last step of our desired time series is computed in the second to last numTS value

        if isinteger(t[i])
            now = Int(t[i])
            if now == 1
                println("\nSimulation Time: $now second")
            else
                println("\nSimulation Time: $now seconds")
            end
        end

        ## check for specified rotor speed at t+dt
        if (model.turbineStartup == 0)
            model.omegaControl = true #TODO: are we setting this back?
            if (model.usingRotorSpeedFunction) #use user specified rotor speed profile function
                _,omegaCurrent,_ = getRotorPosSpeedAccelAtTime(t[i],t[i+1],0.0,delta_t)
                Omega_s = omegaCurrent
            else #use discreteized rotor speed profile function
                omegaCurrent,OmegaDotCurrent,terminateSimulation = omegaSpecCheck(t[i+1],model.tocp,model.Omegaocp,delta_t)
                if (terminateSimulation)
                    break
                end
                Omega_s = omegaCurrent
                OmegaDot_s = OmegaDotCurrent
            end
        else
            omegaCurrent = 0.0
        end

        ## initialize "j" Gauss-Sidel iteration
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

        #initialize  platform module related variables only used if (model.hydroOn)
        if model.hydroOn
            u_j_ptfm = Vector(u_j[ptfm_dofs])
            udot_j_ptfm = Vector(udot_j[ptfm_dofs])
            uddot_j_ptfm = Vector(uddot_j[ptfm_dofs])
        end

        #TODO: put these in the model
        TOL = 1e-4  #gauss-seidel iteration tolerance for various modules
        MAXITER = 300 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        platNorm = 0.0
        aziNorm = 1e5
        gbNorm = 0.0 #initialize norms for various module states

        # Assignments
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
        Fexternal = 0.0 #TODO: do this right, especially if there are only hydro forces
        Fdof = 1

        while ((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER)) #module gauss-seidel iteration loop
            println("Iteration $numIterations")

            rbData = zeros(9)

            #calculate CP2H (platform frame to hub frame transformation matrix)
            CP2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]

            #.... inertial frame to platform transformation matrix ...........
            if (model.hydroOn)
                CN2P = transMat(ptfm_roll, ptfm_pitch, ptfm_yaw)
            else
                CN2P=1.0*LinearAlgebra.I(3)            
            end
            CN2H = CN2P*CP2H
            iCN2H = inv(CN2H)

            #.........................................
            ## evaluate generator module
            #----- generator module ---------------------------
            genTorque_j = 0
            if model.generatorOn
                if model.driveTrainOn
                    if model.useGeneratorFunction
                        genTorqueHSS0 = userDefinedGenerator(gbDot_j*model.gearRatio)
                    else
                        error("simpleGenerator not fully implemented")#[genTorqueHSS0] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio)
                    end
                else
                    if model.useGeneratorFunction
                        genTorqueHSS0 = userDefinedGenerator(Omega_j)
                    else
                        genTorqueHSS0 = simpleGenerator(model,Omega_j)
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

            #------ update rotor speed ---------------------------------
            azi_jLast = azi_j
            if model.omegaControl
                if (model.usingRotorSpeedFunction)
                    azi_j,Omega_j,OmegaDot_j = getRotorPosSpeedAccelAtTime(t[i],t[i+1],azi_s,delta_t)
                else
                    Omega_j = Omega_s
                    OmegaDot_j = OmegaDot_s
                    azi_j = azi_s + Omega_j*delta_t*2*pi
                end
            elseif !model.omegaControl
                Crotor = 0
                Krotor = 0
                azi_j,Omega_j,OmegaDot_j = updateRotorRotation(structureMOI[3,3],Crotor,Krotor,
                -FReaction_j[6],-torqueDriveShaft_j,
                azi_s,Omega_s,OmegaDot_s,delta_t)
            else
                error("omega control option not correctly specified")
            end

            #initialization of structural dynamics displacements, velocities, accelerations, etc.
            

            # --------------------------------------------------------
            #                      Solve Procedure
            # --------------------------------------------------------

            ## Calculate aerodynamic forcing
            deformAero(Omega_j*2*pi)
            # Fexternal, Fdof = aero(t[i]) #TODO: implement turbine deformation and deformation induced velocities

            # Evaluate structural dynamics based on aerodynamic forcing and outputs from time t.
            if model.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,rom,Fexternal,Int.(Fdof),CN2H,rbData)
            else # evalulate structural dynamics using conventional representation
                # println(Fexternal)
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,Fexternal,Int.(Fdof),CN2H,rbData)
            end
                  
            # update current motion estimates and structures
            u_jLast = u_j
            u_jLast_ptfm = u_jLast[ptfm_dofs]
            u_j = dispOut.displ_sp1
            u_j_ptfm = u_j[ptfm_dofs]
            udot_j  = dispOut.displdot_sp1
            udot_j_ptfm = udot_j[ptfm_dofs]
            uddot_j = dispOut.displddot_sp1
            uddot_j_ptfm = uddot_j[ptfm_dofs]

            if model.analysisType=="ROM"
                dispData_j = GyricFEA.DispData(u_j,udot_j,uddot_j,u_sm1,eta_s,etadot_s,etaddot_s) #TODO: fix the j/s mixing in this line
            else
                dispData_j = GyricFEA.DispData(u_j,udot_j,uddot_j,u_sm1)
            end

            ## Evaluate hydro-structural dynamics
            if model.hydroOn

                # transform platform motions from hub frame to inertial frame
                u_j_ptfm[1:3] = u_j_ptfm[1]*iCN2H[1,:] + u_j_ptfm[2]*iCN2H[2,:] + u_j_ptfm[3]*iCN2H[3,:]
                u_j_ptfm[4:6] = u_j_ptfm[4]*iCN2H[1,:] + u_j_ptfm[5]*iCN2H[2,:] + u_j_ptfm[6]*iCN2H[3,:]
                udot_j_ptfm[1:3] = udot_j_ptfm[1]*iCN2H[1,:] + udot_j_ptfm[2]*iCN2H[2,:] + udot_j_ptfm[3]*iCN2H[3,:]
                udot_j_ptfm[4:6] = udot_j_ptfm[4]*iCN2H[1,:] + udot_j_ptfm[5]*iCN2H[2,:] + udot_j_ptfm[6]*iCN2H[3,:]
                uddot_j_ptfm[1:3] = uddot_j_ptfm[1]*iCN2H[1,:] + uddot_j_ptfm[2]*iCN2H[2,:] + uddot_j_ptfm[3]*iCN2H[3,:]
                uddot_j_ptfm[4:6] = uddot_j_ptfm[4]*iCN2H[1,:] + uddot_j_ptfm[5]*iCN2H[2,:] + uddot_j_ptfm[6]*iCN2H[3,:]
                
                println(u_j_ptfm)
                println(udot_j_ptfm)
                println(uddot_j_ptfm)
                
                VAWTHydro.HD_UpdateStates(t[i], t[i+1], u_j_ptfm, udot_j_ptfm, uddot_j_ptfm)
                if model.interpOrder == 1
                    VAWTHydro.MD_UpdateStates(0, t[i], t[i+1], u_j_ptfm, udot_j_ptfm, uddot_j_ptfm)
                elseif model.interpOrder == 2
                    VAWTHydro.MD_UpdateStates(t[i]-delta_t, t[i], t[i+1], u_j_ptfm, udot_j_ptfm, uddot_j_ptfm)
                end

                frc_mooring_n[:], mooring_tensions[:] = VAWTHydro.MD_CalcOutput(t[i+1], u_j_ptfm, udot_j_ptfm, uddot_j_ptfm, frc_mooring_n, mooring_tensions)
                

                elStrain, dispOut_j, FReaction_j, frc_hydro_n, out_vals, jac = OWENS_HD_Coupling(t[1], delta_t, false, jac, numDOFPerNode, ptfm_dofs, frc_hydro_n,
                                                                                    frc_mooring_n, dispData_j, feamodel, mesh, el, elStorage, Omega_j, OmegaDot_j,
                                                                                    Fexternal, Fdof, zeros(3,3))


                println(frc_hydro_n)
                println(frc_mooring_n)

                # transform forces/moments calculated in inertial reference frame back to hub reference frame for the structural solve
                frc_hydro_h[1:3] = frc_hydro_n[1]*CN2H[1,:] + frc_hydro_n[2]*CN2H[2,:] + frc_hydro_n[3]*CN2H[3,:]
                frc_hydro_h[4:6] = frc_hydro_n[4]*CN2H[1,:] + frc_hydro_n[5]*CN2H[2,:] + frc_hydro_n[6]*CN2H[3,:]
                frc_mooring_h[1:3] = frc_mooring_n[1]*CN2H[1,:] + frc_mooring_n[2]*CN2H[2,:] + frc_mooring_n[3]*CN2H[3,:]
                frc_mooring_h[4:6] = frc_mooring_n[4]*CN2H[1,:] + frc_mooring_n[5]*CN2H[2,:] + frc_mooring_n[6]*CN2H[3,:]
            
                # update current motion estimates and structures
                u_jLast = u_j
                u_jLast_ptfm = u_jLast[ptfm_dofs]
                u_j = dispOut_j.displ_sp1
                u_j_ptfm = u_j[ptfm_dofs]
                udot_j  = dispOut_j.displdot_sp1
                udot_j_ptfm = udot_j[ptfm_dofs]
                uddot_j = dispOut_j.displddot_sp1
                uddot_j_ptfm = uddot_j[ptfm_dofs]
            end


            ## calculate norms
            uNorm = 0.0 #LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
            aziNorm = 0.0 #LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm
            platNorm = 0.0 #LinearAlgebra.norm(u_j_ptfm - u_jLast_ptfm)/LinearAlgebra.norm(u_j_ptfm) #platform module states iteration norm if it is off, the norm will be zero
            gbNorm = 0.0 #LinearAlgebra.norm(gb_j - gb_jLast)/LinearAlgebra.norm(gb_j) #gearbox states iteration norm if it is off, the norm will be zero

            numIterations = numIterations + 1
            if numIterations==MAXITER
                @warn "Maximum Iterations Met Breaking Iteration Loop"
                break
            end

        end #end iteration while loop

        ## calculate converged generator torque/power
        genTorquePlot = 0
        if (model.useGeneratorFunction)
            if (model.generatorOn || (model.turbineStartup==0))
                println("simpleGenerator not fully implemented")#[genTorquePlot] = simpleGenerator(model.generatorProps,gbDot_j*model.gearRatio)
            end
        end
        genPowerPlot = genTorquePlot*(gbDot_j*2*pi)*model.gearRatio

        ## update timestepping variables and other states, store in history arrays

        u_s = u_j
        udot_s = udot_j
        uddot_s = uddot_j

        if model.analysisType=="ROM"
            eta_s = dispOut.eta_sp1 #eta_j
            etadot_s = dispOut.etadot_sp1 #etadot_j
            etaddot_s = dispOut.etaddot_sp1 #etaddot_j
        end

        uHist[i+1,:] = u_s
        uHist_ptfm[i+1,:] = u_j_ptfm
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

        ## check rotor speed for generator operation
        if Omega_s >= rotorSpeedForGenStart
            model.generatorOn = true
        else
            model.generatorOn = false
        end

    end #end timestep loop

    # End FAST module links
    if model.hydroOn
        VAWTHydro.HD_End()
        VAWTHydro.MD_End()
    end

    #Writefile
    if model.outFilename=="none"
        println("NOT WRITING Verification File")
    else
        println("WRITING Verification File")

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
            HDF5.write(file,"genTorque",genTorque)
            HDF5.write(file,"genPower",genPower)
            HDF5.write(file,"torqueDriveShaft",torqueDriveShaft)
            HDF5.write(file,"uHist",uHist)
            HDF5.write(file,"uHist_ptfm",uHist_ptfm)
            HDF5.write(file,"eps_xx_0_hist",eps_xx_0_hist)
            HDF5.write(file,"eps_xx_z_hist",eps_xx_z_hist)
            HDF5.write(file,"eps_xx_y_hist",eps_xx_y_hist)
            HDF5.write(file,"gam_xz_0_hist",gam_xz_0_hist)
            HDF5.write(file,"gam_xz_y_hist",gam_xz_y_hist)
            HDF5.write(file,"gam_xy_0_hist",gam_xy_0_hist)
            HDF5.write(file,"gam_xy_z_hist",gam_xy_z_hist)
        end

    end
    return t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_ptfm,eps_xx_0_hist,eps_xx_z_hist,eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist
end
