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

"""
    transMat(theta1, theta2, theta3)

Internal, computes the 3x3 transformation matrix for given input rotations. The generated matrix
is the closest orthonormal matrix to the Bernoulli-Euler transformation matrix from beam theory,
which assumes small rotations. A full description of this matrix is found in the
"FASTCoordinateSystems.doc" document by Jason Jonkman.
"""
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

"""
    extrap_pred_vals(curr_vals, ts, t_out, interp_order)

Internal, calculates predicted values at t_out based on previous values at earlier times.

# Input
* `curr_vals::Vector{<:float}`: input values
* `ts::Vector{<:float}`: time values corresponding to curr_vals. Left to right, times go from earliest to most recent.
* `t_out::float`: time for values to be extrapolated tocp
* `interp_order::int`: order of the spline fit for the extrapolation, 0 flat, 1 linear, 2 quadratic

# Output
* `pred_vals`: extrapolated values at t_out
"""
function extrap_pred_vals(curr_vals, ts, t_out, interp_order)

    pred_vals = zeros(length(curr_vals[:,1]))
    for i = 1:size(curr_vals, 1)
        spl = Dierckx.Spline1D(ts, curr_vals[i,:], k=interp_order, bc="extrapolate")
        pred_vals[i] = Dierckx.evaluate(spl, t_out)
    end

    return pred_vals

end

"""
    frame_convert(init_frame_vals, trans_mat)

Internal, transfers 6 DOFs element-wise to a new reference frame

# Input
* `init_frame_vals::Vector{<:float}`: Values in 6 degrees of freedom in the initial reference frame
* `trans_mat::Array{<:float}`: Transformation matrix to the output reference frame

# Output
* `out_frame_vals`: Values in 6 degrees of freedom in the output reference frame
"""
function frame_convert(init_frame_vals, trans_mat)

    out_frame_vals = zero(init_frame_vals)
    out_frame_vals[1:3] = init_frame_vals[1]*trans_mat[1,:] + init_frame_vals[2]*trans_mat[2,:] + init_frame_vals[3]*trans_mat[3,:]
    out_frame_vals[4:6] = init_frame_vals[4]*trans_mat[1,:] + init_frame_vals[5]*trans_mat[2,:] + init_frame_vals[6]*trans_mat[3,:]       
    
    return out_frame_vals

end

"""

    calc_hydro_residual(new_accels, new_hydro_frcs, md_frc, u, frc_multiplier)

Internal, finds the residual of the platform accelerations/forces when adding in new values
from HydroDyn, MoorDyn, and GyricFEA.

# Input
* `new_accels::Vector{<:float}`: the new platform accelerations in the inertial reference frame from GyricFEA  (m/s/s)
* `new_hydro_frcs::Vector{<:float}`: the new platform hydro forces in the inertial reference frame from HydroDyn (N)
* `md_frc::Vector{<:float}`: the mooring forces at the platform in the inertial reference frame from MoorDyn (N)
* `u::Vector{<:float}`: input vector containing platform loads and accelerations prior to the new additions
* `frc_multiplier::Number`: scaling factor used to make the loads/accelerations in the same magnitude for the Jacobians

# Output
* `u_resid`: vector containing residual between the input u and the new load/accelerations
"""
function calc_hydro_residual(new_accels, new_hydro_frcs, md_frc, u, frc_multiplier)

    new_frcs = md_frc + new_hydro_frcs

    u_resid = Vector{Float32}(undef, length(new_frcs) + length(new_accels))
    u_resid[1:length(new_frcs)] = u[1:length(new_frcs)] - new_frcs/frc_multiplier
    u_resid[length(new_frcs)+1:end] = u[length(new_frcs)+1:end] - new_accels

    return u_resid
end

"""

    OWENS_HD_Coupled_Solve(time, dt, calcJacobian, jac, numDOFPerNode, ptfm_dofs, ptfm_mass, rZT, rYT, frc_hydro_in,
        frc_mooring_in, dispIn, feamodel, mesh, el, elStorage, Omega, OmegaDot,
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
* `ptfm_mass::float`: mass of the floating platform (not including the ballast) (kg)
* `rZT::Vector{<:float}`: see ?Model.ptfmref2bs
* `rYT::Vector{<:float}`: see ?Model.ptfmcom2bs
* `frc_hydro_in::Vector{<:float}`: previously calculated hydrodynamic forces in the inertial frame (N)
* `frc_mooring_in::Vector{<:float}`: previously calculated mooring forces in the inertial frame (N)
* `dispIn::DispData`: see ?GyricFEA.DispData
* `feamodel::FEAModel`: see ?GyricFEA.FEAModel
* `mesh::Mesh`: see ?GyricFEA.Mesh
* `el::El`: see ?GyricFEA.El,
* `elStorage::ElStorage`: see ?GyricFEA.elStorage
* `Omega::float`: rotor speed (Hz)
* `OmegaDot::float`: rotor acceleration (Hz/s)
* `other_Fexternal::Vector{<:float}`: other external loads acting on the structure besides the hydro and mooring loads (typically aero loads) (N)
* `other_Fdof::Vector{<:float}`: vector of nodes with indices corresponding to the loads given in other_Fexternal`
* `CN2H::Array{<:float}`: rotation matrix defining the transformation from the inertial to the hub reference frames
* `u_ptfm_n::Vector{<:float}`: the platform position in the inertial reference frame at the current simulation time (m)
* `udot_ptfm_n::Vector{<:float}`: the platform velocity in the inertial reference frame at the current simulation time (m/s)
* `uddot_ptfm_n::Vector{<:float}`: the platform acceleration in the inertial reference frame at the current simulation time (m/s/s)
* `rom::int/ROM`: see ?GyricFEA.ROM. Instead defaults to 0 if the reduced order model is not used

# Output
* `elStrain_out`: object containing strains on each element returned by GyricFEA
* `dispOut`: object containing motions on each element returned by GyricFEA
* `FReaction_out`: vector of the total loads acting on the specified nodes
* `frc_hydro_out`: vector of the hydrodynamic loads acting on the platform returned by HydroDyn
* `out_vals`: vector containing the output parameters specified at the end of the HydroDyn input file
* `jac_out`: the platform hydro forces and accelerations Jacobian calculated internally if calcJacobian=true (passthrough if calcJacobian=false)
"""

function OWENS_HD_Coupled_Solve(time, dt, calcJacobian, jac, numDOFPerNode, ptfm_dofs, frc_ptfm_old, frc_mooring_new,
    dispIn, feamodel, mesh, el, elStorage, Omega, OmegaDot,
    other_Fexternal, other_Fdof, CN2H, u_ptfm_n, udot_ptfm_n, uddot_ptfm_n, rom=0)

    # !!! Make sure structuralDynamicsTransient and MD_CalcOutput have run before this so you can send dispOut and frc_mooring_new here !!!
    # !!! Hydro/mooring force inputs/outputs are in the inertial reference frame, platform motion inputs/outputs are in the hub reference frame !!!
    # 
    # allocate new variables
    frc_hydro2 = Vector{Float32}(undef, numDOFPerNode)
    if calcJacobian
        
    end
    frc_hydro_out = Vector{Float32}(undef, numDOFPerNode) # TODO: all of the frc_hydro's may be able to be combined for efficiency
    out_vals = Vector{Float32}(undef, numDOFPerNode+1)
    u = Vector{Float32}(undef, numDOFPerNode*2)
    
    frc_multiplier = 1e6

    # u_ptfm_n = frame_convert(dispIn.displ_s[ptfm_dofs], LinearAlgebra.inv(CN2H))
    # udot_ptfm_n = frame_convert(dispIn.displdot_s[ptfm_dofs], LinearAlgebra.inv(CN2H))
    # uddot_ptfm_n = frame_convert(dispIn.displddot_s[ptfm_dofs], LinearAlgebra.inv(CN2H))
    
    # For consistency, everything in the u vector is in the inertial reference frame
    u[1:numDOFPerNode] = frc_ptfm_old / frc_multiplier
    u[numDOFPerNode+1:numDOFPerNode*2] = uddot_ptfm_n

    # Calculate outputs at the current time, based on inputs at the current time
    total_Fexternal = [other_Fexternal; frame_convert(frc_ptfm_old, CN2H)] # platform loads are old inputs from last time step/correction
    total_Fdof = [other_Fdof; ptfm_dofs]
    if rom != 0
        _ ,dispOut2, _ = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,rom,total_Fexternal,Int.(total_Fdof),CN2H,zeros(9))
    else
        _ ,dispOut2, _ = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal,Int.(total_Fdof),CN2H,zeros(9))
    end
    frc_hydro2[:], _ = VAWTHydro.HD_CalcOutput(time, u_ptfm_n, udot_ptfm_n, uddot_ptfm_n, frc_hydro2, out_vals)

   # Calculate the residual
   uddot_ptfm_n2 = frame_convert(dispOut2.displddot_sp1[ptfm_dofs], LinearAlgebra.inv(CN2H))
   residual = calc_hydro_residual(uddot_ptfm_n2, frc_hydro2, frc_mooring_new, u, frc_multiplier)

   # Calculate the Jacobian. Since there's no single function associated with the motions/forces, we will manually
   # perturb each degree of freedom one at a time and recalculate the outputs so we can get a full gradient.
    if calcJacobian

        for dof = collect(1:numDOFPerNode) # forces
            frc_ptfm_perturb = copy(frc_ptfm_old)
            u_perturb = copy(u)
            frc_ptfm_perturb[dof] += 1E6
            u_perturb[dof] += 1
            total_Fexternal_perturb = [other_Fexternal; frame_convert(frc_ptfm_perturb, CN2H)]
            if !isa(rom, Number)
                _ ,dispOut_perturb, _ = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,rom,total_Fexternal_perturb,Int.(total_Fdof),CN2H,zeros(9))
            else
                _, dispOut_perturb, _ = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal_perturb,Int.(total_Fdof),CN2H,zeros(9))
            end
            uddot_ptfm_perturb_n = frame_convert(dispOut_perturb.displddot_sp1[ptfm_dofs], LinearAlgebra.inv(CN2H))
            residual_perturb = calc_hydro_residual(uddot_ptfm_perturb_n, frc_hydro2, frc_mooring_new, u_perturb, frc_multiplier)
            jac[:,dof] = residual_perturb - residual
        end

        for dof = collect(1:numDOFPerNode) # accelerations
            uddot_ptfm_perturb = copy(uddot_ptfm_n)
            u_perturb = copy(u)
            uddot_ptfm_perturb[dof] += 1
            u_perturb[dof+numDOFPerNode] += 1
            frc_hydro_perturb = Vector{Float32}(undef, numDOFPerNode) # this is reset each time, otherwise HydroDyn returns garbage values after the first iteration
            frc_hydro_perturb[:], _ = VAWTHydro.HD_CalcOutput(time, u_ptfm_n, udot_ptfm_n, uddot_ptfm_perturb, frc_hydro_perturb, out_vals)
            residual_perturb = calc_hydro_residual(uddot_ptfm_n2, frc_hydro_perturb, frc_mooring_new, u_perturb, frc_multiplier)
            jac[:,dof+numDOFPerNode] = residual_perturb - residual
        end
 
    end  # if calcJacobian

    # Solve for delta_u: jac*delta_u = -residual
    delta_u = jac\-residual
    jac_out = jac

    # Update inputs
    frc_ptfm_rev = frc_ptfm_old + delta_u[1:numDOFPerNode]*frc_multiplier
    uddot_ptfm_n_rev = uddot_ptfm_n + delta_u[numDOFPerNode+1:end]

    total_Fexternal_rev = [other_Fexternal; frame_convert(frc_ptfm_rev, CN2H)]
    # dispData_rev = deepcopy(dispIn)
    # dispData_rev.displddot_s[ptfm_dofs] = frame_convert(dispIn.displddot_s[ptfm_dofs], LinearAlgebra.inv(CN2H)) + delta_u[numDOFPerNode+1:numDOFPerNode*2]
    # uddot_ptfm_n = dispData_rev.displddot_s[ptfm_dofs]

    # Rerun OWENS and HydroDyn with updated inputs
    if rom != 0
        elStrain_out,dispOut,FReaction_out = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,rom,total_Fexternal_rev,Int.(total_Fdof),CN2H,zeros(9))
    else
        elStrain_out,dispOut,FReaction_out = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispIn,Omega,OmegaDot,time,dt,elStorage,total_Fexternal_rev,Int.(total_Fdof),CN2H,zeros(9)) # TODO: should we use dispData_rev instead of dispIn here?
    end
    frc_hydro_out[:], out_vals[:] = VAWTHydro.HD_CalcOutput(time, u_ptfm_n, udot_ptfm_n, uddot_ptfm_n_rev, frc_hydro_out, out_vals)


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
* `bin::Bin`: see ?Bin
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

    #..........................................................................
    #                             INITIALIZATION
    #..........................................................................

    ## General
    delta_t = model.delta_t
    numTS = Int(model.numTS)
    t = range(0, length=numTS, step=delta_t)
    numDOFPerNode = 6
    totalNumDof = mesh.numNodes*numDOFPerNode
    CN2H = LinearAlgebra.I(3) # hub and inertial frames initialize as copies

    ## Initial conditions
    u_s = zeros(totalNumDof)
    u_s = GyricFEA.setInitialConditions(feamodel.initCond,u_s,numDOFPerNode)
    u_sm1 = copy(u_s)
    udot_s = zero(u_s)
    uddot_s = zero(u_s)
    rbData = zeros(9)
    dispData = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1)

    gb_s = 0
    gbDot_s = 0
    gbDotDot_s = 0
    azi_s = 0
    Omega_s = copy(model.OmegaInit)
    OmegaDot_s = 0
    genTorque_s = 0
    torqueDriveShaft_s = 0

    Fexternal = 0.0
    Fdof = 1

    ## Hydrodynamics/mooring module initialization and coupling variables
    if model.hydroOn

        ptfm_dofs = 1:numDOFPerNode
        ptfm_mass = feamodel.nodalTerms.concMass[1].val #TODO: make it so the platform concentrated mass terms doesn't have to be the first concentrated term defined

        u_s_ptfm = Vector(u_s[ptfm_dofs])
        udot_s_ptfm = Vector(udot_s[ptfm_dofs])
        uddot_s_ptfm = Vector(uddot_s[ptfm_dofs])

        jac = Array{Float64}(undef, numDOFPerNode*2, numDOFPerNode*2)
        numMooringLines = 3

        if model.outFilename == "none"
            hd_outFilename = "hydrodyn_temp.out"
        else
            hd_outFilename = model.outFilename
        end
        
        frc_hydro_n = zeros(numDOFPerNode) #Vector{Float32}(undef, numDOFPerNode)
        frc_mooring_n = zeros(numDOFPerNode) #Vector{Float32}(undef, numDOFPerNode)
        out_vals = Vector{Float32}(undef, numDOFPerNode+1) # Rigid body displacement in 6DOF + wave elevation
        mooring_tensions = Vector{Float32}(undef, numMooringLines*2) # Fairlead + anchor tension for each line

        VAWTHydro.HD_Init(bin.hydrodynLibPath, hd_outFilename, hd_input_file=model.hd_input_file, PotFile=model.potflowfile, t_initial=t[1], dt=delta_t, t_max=t[1]+(numTS-1)*delta_t)
        VAWTHydro.MD_Init(bin.moordynLibPath, md_input_file=model.md_input_file, init_ptfm_pos=u_s_ptfm, interp_order=model.interpOrder)
    end

    ## Rotor mode initialization

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

    ## Structural dynamics initialization
    if model.analysisType=="ROM"
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

    _,structureMOI,_ = GyricFEA.calculateStructureMassProps(elStorage)

    feamodel.jointTransform, feamodel.reducedDOFList = GyricFEA.createJointTransform(feamodel.joint,mesh.numNodes,6) #creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
    
    ## History array initialization
    uHist = zeros(numTS, length(u_s))
    uHist[1,:] = u_s          #store initial condition

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

    if model.hydroOn
        uHist_ptfm = zeros(numTS,length(u_s_ptfm))
        uHist_ptfm[1,:] = u_s_ptfm
        FHydroHist = zeros(numTS,6)
    end

    #..........................................................................
    #                          INITIAL COUPLED SOLVE
    #..........................................................................

    ## TODO: I think structuralDynamicsTransient *might* need to be called here if we want to mirror ElastoDyn-HydroDyn coupling perfectly, but since we already initialize displacements
    ##       without running it, I think it's good unless we see instability.

    # Calculate aerodynamic forcing
    deformAero(Omega_s*2*pi)
    # frc_aero_h, aero_dofs = aero(t[i]) #TODO: implement turbine deformation and deformation induced velocities
    frc_aero_h = 0.0
    aero_dofs = 1

    ## Evaluate mooring and hydrodynamics at t=0 based on initial conditions
    if model.hydroOn
        moms_ptfm2bs_n = vcat(zeros(3), # 3 force DOFs (no additions here, only moments below)
        LinearAlgebra.cross(model.ptfmref2bs, (frc_hydro_n[1:3] + frc_mooring_n[1:3])) - # add the moments about the platform node due to the hydrodynamic forces at the platform reference point
        ptfm_mass * LinearAlgebra.cross(model.ptfmcom2bs, [0,0,9.81])) # add the moments about the platform node due to gravitational effects at the platform center of mass (gravity acts in the negative direction)
        frc_ptfm_n = frc_hydro_n + frc_mooring_n + moms_ptfm2bs_n

        frc_mooring_n = Vector{Float32}(undef, numDOFPerNode) # this is reset here, otherwise MoorDyn returns garbage values in MD_CalcOutput
        frc_mooring_n[:], mooring_tensions[:] = VAWTHydro.MD_CalcOutput(t[1], u_s_ptfm, udot_s_ptfm, uddot_s_ptfm, frc_mooring_n, mooring_tensions)
        if model.analysisType=="ROM"
            elStrain, dispOut, FReaction, frc_hydro_n, out_vals, jac = OWENS_HD_Coupled_Solve(t[1], delta_t, true, jac, numDOFPerNode, ptfm_dofs, frc_ptfm_n, frc_mooring_n,
                                                                                    dispData, feamodel, mesh, el, elStorage, Omega_s, OmegaDot_s,
                                                                                    frc_aero_h, aero_dofs, CN2H, u_s_ptfm, udot_s_ptfm, uddot_s_ptfm, rom)
        else
            elStrain, dispOut, FReaction, frc_hydro_n, out_vals, jac = OWENS_HD_Coupled_Solve(t[1], delta_t, true, jac, numDOFPerNode, ptfm_dofs, frc_ptfm_n, frc_mooring_n,
                                                                                    dispData, feamodel, mesh, el, elStorage, Omega_s, OmegaDot_s,
                                                                                    frc_aero_h, aero_dofs, CN2H, u_s_ptfm, udot_s_ptfm, uddot_s_ptfm)
        end

        # Don't need to do a frame conversion here since the inertial and hub frames are aligned at t=0
        u_sm1 = copy(u_s) # not sure if this is necessary (or correct), since the initial u_s isn't technically at t=-dt
        u_s = dispOut.displ_sp1
        udot_s = dispOut.displdot_sp1
        uddot_s = dispOut.displddot_sp1
        dispData = GyricFEA.DispData(u_s,udot_s,uddot_s,u_sm1)
        
        moms_ptfm2bs_n_s = vcat(zeros(3), # 3 force DOFs (no additions here, only moments below)
        LinearAlgebra.cross(model.ptfmref2bs, (frc_hydro_n[1:3] + frc_mooring_n[1:3])) - # add the moments about the platform node due to the hydrodynamic forces at the platform reference point
        ptfm_mass * LinearAlgebra.cross(model.ptfmcom2bs, [0,0,9.81])) # add the moments about the platform node due to gravitational effects at the platform center of mass (gravity acts in the negative direction)
        frc_ptfm_n_s = frc_hydro_n + frc_mooring_n + moms_ptfm2bs_n_s

        FHydroHist[1,:] = frc_hydro_n

    ## Copy values of the inital outputs to arrays for interpolation/extrapolation
        recent_u_ptfm = repeat(u_s[ptfm_dofs], 1, model.interpOrder+1)
        recent_udot_ptfm = repeat(udot_s[ptfm_dofs], 1, model.interpOrder+1)
        recent_uddot_ptfm = repeat(uddot_s[ptfm_dofs], 1, model.interpOrder+1)
        recent_frc_ptfm = repeat(frc_ptfm_n_s, 1, model.interpOrder+1)
        recent_times = collect(range(-delta_t*model.interpOrder, 0.0, model.interpOrder+1))
    end

    #..........................................................................
    #                                MAIN LOOP
    #..........................................................................

    ### Iterate for a solution at t+dt
    for i=1:numTS-1 # we compute for the next time step, so the last step of our desired time series is computed in the second to last numTS value

        ## Print current simulation time to terminal
        if isinteger(t[i])
            now = Int(t[i])
            if now == 1
                println("\nSimulation Time: $now second")
            else
                println("\nSimulation Time: $now seconds")
            end
        end

        ## Check for specified rotor speed at t+dt
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

        ## Guess platform motions and forces at t+dt to get proper rotations for reference frame conversions (and to make hydro coupling stable)
        if model.hydroOn
            u_s_ptfm = extrap_pred_vals(recent_u_ptfm, recent_times, t[i+1], model.interpOrder)
            udot_s_ptfm = extrap_pred_vals(recent_udot_ptfm, recent_times, t[i+1], model.interpOrder)
            uddot_s_ptfm = extrap_pred_vals(recent_uddot_ptfm, recent_times, t[i+1], model.interpOrder)
            frc_ptfm_n_s = extrap_pred_vals(recent_frc_ptfm, recent_times, t[i+1], model.interpOrder)
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

        if model.hydroOn
            u_j_ptfm_n = u_s_ptfm
            udot_j_ptfm_n = udot_s_ptfm
            uddot_j_ptfm_n = uddot_s_ptfm
            frc_ptfm_n_j = frc_ptfm_n_s
            frc_hydro_n_j = frc_hydro_n # this doesn't have to be used outside the while loop outside of tracking purposes
        end

        #TODO: put these in the model
        TOL = 1e-4  #gauss-seidel iteration tolerance for various modules
        MAXITER = 300 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        platNorm = 0.0
        aziNorm = 1e5
        gbNorm = 0.0 #initialize norms for various module states

        ## Gauss-Seidel predictor-corrector loop
        while ((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))
            # println("Iteration $numIterations")

            #------------------
            # GENERATOR MODULE
            #------------------
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

            #-------------------
            # DRIVETRAIN MODULE
            #-------------------
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

            # Update rotor speed
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

            #---------------------
            # AERODYNAMICS MODULE
            #---------------------

            # Calculate aerodynamic forcing
            deformAero(Omega_j*2*pi)
            # frc_aero_h, aero_dofs = aero(t[i]) #TODO: implement turbine deformation and deformation induced velocities
            frc_aero_h_j = 0.0
            aero_dofs = 1

            #------------------------------------
            # STRUCTURAL/HYDRO/MOORING MODULES
            #------------------------------------

            # Update reference frame transformation and convert external forces to hub reference frame (note aero forces are already assumed to be in the hub frame)
            CT2H = [cos(azi_j) sin(azi_j) 0;-sin(azi_j) cos(azi_j) 0;0 0 1]
            if (model.hydroOn)
                CP2T = transMat(u_j[model.towertop_dofs[4]]-u_j_ptfm_n[4], -(u_j[model.towertop_dofs[5]]-u_j_ptfm_n[5]), u_j[model.towertop_dofs[6]]-u_j_ptfm_n[6])
                CN2P = transMat(u_j_ptfm_n[4], -u_j_ptfm_n[5], u_j_ptfm_n[6])
            else
                CP2T = transMat(u_j[model.towertop_dofs[4]], -u_j[model.towertop_dofs[5]], u_j[model.towertop_dofs[6]])
                CN2P=1.0*LinearAlgebra.I(3)
            end
            CN2H = CN2P*CP2T*CT2H
            CH2N = inv(CN2H)

            # Evaluate structural dynamics based on values from last iteration (or extrapolated values from last time step if it's the first iteration) and new aerodynamic loads
            Fexternal = [frc_aero_h_j; frame_convert(frc_ptfm_n_j, CN2H)] # apply both aero and previous platform loads (in the hub reference frame!) here, since it is outside the coupling routine
            Fdof = [aero_dofs; ptfm_dofs]
            if model.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i+1],delta_t,elStorage,rom,Fexternal,Int.(Fdof),CN2H,rbData)
            else # evalulate structural dynamics using conventional representation
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i+1],delta_t,elStorage,Fexternal,Int.(Fdof),CN2H,rbData)
            end
                  
            # update current motion estimates and structures
            u_j = dispOut.displ_sp1
            udot_j  = dispOut.displdot_sp1
            uddot_j = dispOut.displddot_sp1

            ## Evaluate hydro-structural dynamics
            if model.hydroOn

                # dispData_j = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1)

                # transform platform motions from hub frame to inertial frame
                u_j_ptfm_n2 = frame_convert(u_j[ptfm_dofs], CH2N)
                udot_j_ptfm_n2 = frame_convert(udot_j[ptfm_dofs], CH2N)
                uddot_j_ptfm_n2 = frame_convert(uddot_j[ptfm_dofs], CH2N)
                
                # FAST updates HD/MD using inputs from last time step/correction instead of transferring the new ElastoDyn outputs, which is why I'm using u*_j_ptfm_n instead of u*_j_ptfm_n2 here
                VAWTHydro.HD_UpdateStates(t[i], t[i+1], u_j_ptfm_n, udot_j_ptfm_n, uddot_j_ptfm_n)
                if model.interpOrder == 1
                    VAWTHydro.MD_UpdateStates(0, t[i], t[i+1], u_j_ptfm_n, udot_j_ptfm_n, uddot_j_ptfm_n)
                elseif model.interpOrder == 2
                    VAWTHydro.MD_UpdateStates(t[i]-delta_t, t[i], t[i+1], u_j_ptfm_n, udot_j_ptfm_n, uddot_j_ptfm_n)
                end

                # Calculate new mooring forces on the platform since MoorDyn doesn't use acceleration inputs (despite the input arguments). Note the new GyricFEA motion outputs are now used as inputs.
                frc_mooring_n_j = Vector{Float32}(undef, numDOFPerNode) # this is reset here, otherwise MoorDyn returns garbage values in MD_CalcOutput
                frc_mooring_n_j[:], mooring_tensions[:] = VAWTHydro.MD_CalcOutput(t[i+1], u_j_ptfm_n2, udot_j_ptfm_n2, uddot_j_ptfm_n2, frc_mooring_n_j, mooring_tensions)

                if model.analysisType=="ROM"
                    elStrain, dispOut_j, FReaction_j, frc_hydro_n_j, out_vals, jac = OWENS_HD_Coupled_Solve(t[i+1], delta_t, false, jac, numDOFPerNode, ptfm_dofs, frc_ptfm_n_j, frc_mooring_n_j,
                                                                                            dispData, feamodel, mesh, el, elStorage, Omega_j, OmegaDot_j,
                                                                                            frc_aero_h_j, aero_dofs, CN2H, u_j_ptfm_n2, udot_j_ptfm_n2, uddot_j_ptfm_n2, rom)
                else
                    elStrain, dispOut_j, FReaction_j, frc_hydro_n_j, out_vals, jac = OWENS_HD_Coupled_Solve(t[i+1], delta_t, false, jac, numDOFPerNode, ptfm_dofs, frc_ptfm_n_j, frc_mooring_n_j,
                                                                                            dispData, feamodel, mesh, el, elStorage, Omega_s, OmegaDot_j,
                                                                                            frc_aero_h_j, aero_dofs, CN2H, u_j_ptfm_n2, udot_j_ptfm_n2, uddot_j_ptfm_n2)
                end
            
                # update current motion estimates and structures
                u_jLast = copy(u_j)
                u_j = dispOut_j.displ_sp1
                udot_j  = dispOut_j.displdot_sp1
                uddot_j = dispOut_j.displddot_sp1
                u_jLast_ptfm_n = frame_convert(u_jLast[ptfm_dofs], CH2N)
                u_j_ptfm_n = frame_convert(u_j[ptfm_dofs], CH2N)
                udot_j_ptfm_n = frame_convert(udot_j[ptfm_dofs], CH2N)
                uddot_j_ptfm_n = frame_convert(uddot_j[ptfm_dofs], CH2N)

                moms_ptfm2bs_n_j = vcat(zeros(3), # 3 force DOFs (no additions here, only moments below)
                LinearAlgebra.cross(model.ptfmref2bs, (frc_hydro_n_j[1:3] + frc_mooring_n_j[1:3])) - # add the moments about the platform node due to the hydrodynamic forces at the platform reference point
                ptfm_mass * LinearAlgebra.cross(model.ptfmcom2bs, [0,0,9.81])) # add the moments about the platform node due to gravitational effects at the platform center of mass (gravity acts in the negative direction)
                frc_ptfm_n_j = frc_hydro_n_j + frc_mooring_n_j + moms_ptfm2bs_n_j

                # dispData_j = GyricFEA.DispData(u_j, udot_j, uddot_j, u_sm1)

                # println(u_j_ptfm_n)
                # println(udot_j_ptfm_n)
                # println(uddot_j_ptfm_n)
                # println(frc_hydro_n_j)
                # println(frc_mooring_n_j)
            end

            ## calculate norms
            uNorm = 0.0 #LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
            aziNorm = 0.0 #LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm
            platNorm = 0.0 #LinearAlgebra.norm(u_j_ptfm_n - u_jLast_ptfm_n)/LinearAlgebra.norm(u_j_ptfm_n) #platform module states iteration norm if it is off, the norm will be zero
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
        u_sm1 = copy(u_s)
        u_s = u_j
        udot_s = udot_j
        uddot_s = uddot_j

        if model.analysisType=="ROM"
            eta_s = dispOut.eta_sp1 #eta_j
            etadot_s = dispOut.etadot_sp1 #etadot_j
            etaddot_s = dispOut.etaddot_sp1 #etaddot_j
            dispData = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1, eta_s, etadot_s, etaddot_s)
        else
            dispData = GyricFEA.DispData(u_s, udot_s, uddot_s, u_sm1)
        end

        if model.hydroOn
            u_s_ptfm = u_j_ptfm_n
            udot_s_ptfm = udot_j_ptfm_n
            uddot_s_ptfm = uddot_j_ptfm_n
            frc_ptfm_n_s = frc_ptfm_n_j

            uHist_ptfm[i+1,:] = u_s_ptfm
            FHydroHist[i+1,:] = frc_hydro_n_j

            # Shift up window of extrapolation vectors
            recent_u_ptfm = hcat(recent_u_ptfm[:, 2:end], u_s_ptfm)
            recent_udot_ptfm = hcat(recent_udot_ptfm[:, 2:end], udot_s_ptfm)
            recent_uddot_ptfm = hcat(recent_uddot_ptfm[:, 2:end], uddot_s_ptfm)
            recent_frc_ptfm = hcat(recent_frc_ptfm[:, 2:end], frc_ptfm_n_s)
            recent_times = vcat(recent_times[2:end], t[i+1])
        end
        
        uHist[i+1,:] = u_s
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

    println("Simulation Complete.")

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
    return t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_ptfm,eps_xx_0_hist,eps_xx_z_hist,eps_xx_y_hist,gam_xz_0_hist,gam_xz_y_hist,gam_xy_0_hist,gam_xy_z_hist,FHydroHist
end
