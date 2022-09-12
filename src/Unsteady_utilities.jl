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
        spl = FLOWMath.Akima(tocp,Omegaocp)
        OmegaCurrent = spl(tCurrent)
        OmegaDotCurrent = FLOWMath.derivative(spl,tCurrent)

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

# If startup

# If normal operation

# If shutdown

function userDefinedGenerator(newVinf,t,gb_j,omega,omegalast,omegadot,omegadotlast,dt,integrator,omegasetpoint)
    # omega is in hz
    omega_RPM = omega*60
    omegalast_RPM = omegalast*60

    omegadot_RPM = omegadot*60
    omegadotlast_RPM = omegadotlast*60

    # if omegasetpoint*60<=0.99 && isapprox(0.0,omega_RPM;atol=1.0)
    #     controlnamecurrent = "off"
    #     println("off setpoint $(omegasetpoint*60) RPM $(omega_RPM)")
    #     controllerQ = 0
    # elseif isapprox(omegasetpoint,omega;atol=omegasetpoint*0.07)#operPhase == "normal" #Level controller
        controlnamecurrent = "normal"
        # println(" ")
        # println("normal setpoint $(omegasetpoint*60) RPM $(omega*60)")
        Kpfactor = 1.0#newVinf^2 / 10.0^2
        omega_RPM0 = omegasetpoint*60 #33.92871
        Kp = 10.62992345720471*Kpfactor
        Ki = 5.2553876053628725*Kpfactor
        Kd = 0.0
        Q0 = -320.1533668398164*Kpfactor
        integrator = integrator + (omega_RPM - omega_RPM0)*dt
        deriv = (omega_RPM-omegalast_RPM)/dt
        controllerQ = Q0 + Kp*(omega_RPM) + Kd*deriv + Ki*integrator
        # if t>40.0
        #     controllerQ = 145.0
        # end
        # if omega_RPM<0.0 || t>75.0
        #     controllerQ = 20*(omega_RPM)
        # end
        # println("$t $controllerQ")
    #     # # Synchronous Generator: 17m generator from SAND-78-0577 x10 and scaled up for rotation rate
    #     # k = 57.6 #N-m-s/rad
    #     # D = 1.27e3 #N-m/rad
    #     # ws = omega_RPM0 / 60 * 2 * pi
    #     # scale = 10 #x larger than the 17m
    #     # gbMultiplier = 52.94
    #     # controllerQ = (k*(gb_j-ws*t) + D*(omega*2*pi-ws))/1000 * gbMultiplier * scale
    #
    #     # println("Integrator $(integrator[1]) Ki*Int: $(Ki*integrator[1])")
    #     # println("Kp*omegadot_RPM $(Kp*omegadot_RPM)")
    #     # println("controllerQ $(controllerQ)")
    # elseif omegasetpoint>omega#operPhase == "startup" #Rate controller
    #     controlnamecurrent = "startup"
    #     println(" ")
    #     println("startup setpoint $(omegasetpoint*60) RPMdot $(omegadot_RPM)")
    #     omegadot_RPM0 = 0.1076
    #     Kp = 0.0
    #     Kd = 610.8381799970701*0 #TODO: derivative based on GB position
    #     Ki = 37.79082483023065
    #     Q0 = 39.74991266082707
    #     integrator = integrator + (omegadot_RPM - omegadot_RPM0)*dt
    #     deriv = (omegadot_RPM-omegadotlast_RPM)/dt
    #     QKp = Kp*(omegadot_RPM)
    #
    #     controllerQ = Q0 + QKp + Kd*deriv + Ki*integrator
    #
    #     if controllerQ>100.0
    #         controllerQ = 100.0
    #     elseif controllerQ<-100.0
    #         controllerQ = -100.0
    #     end
    #
    #     println("Integrator $(integrator[1]) Ki*Int: $(Ki*integrator[1])")
    #     println("Kp*omegadot_RPM $(Kp*omegadot_RPM)")
    #     println("controllerQ $(controllerQ)")
    #
    # elseif omegasetpoint<omega#operPhase == "alarmstop" #Rate controller
    #     controlnamecurrent = "alarmstop"
    #     println(" ")
    #     println("alarmstop setpoint $(omegasetpoint*60) RPMdot $(omegadot_RPM)")
    #     omegadot_RPM0 = -1.04474
    #     Kp = -5.592014753661118
    #     Kd = -0.6525496964068226*0 #TODO: derivative based on GB position
    #     Ki = -3.2809838749728284
    #     Q0 = 123.94659715178108
    #     integrator = integrator + (omegadot_RPM - omegadot_RPM0)*dt
    #     deriv = (omegadot_RPM-omegadotlast_RPM)/dt
    #     controllerQ = Q0 + Kp*(omegadot_RPM) + Kd*deriv + Ki*integrator
    #     println("Integrator $(integrator[1]) Ki*Int: $(Ki*integrator[1])")
    #     println("Kp*omegadot_RPM $(Kp*omegadot_RPM)")
    #     println("controllerQ $(controllerQ)")
    # end
    #
    # if omega_RPM>41.0
    #     controllerQ = 150.0*2
    # end

    return controllerQ*1000,integrator,controlnamecurrent
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
    transMat = Array{Float32}(undef, 3,3)
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

    Fexternal = zeros(Float32, length(ForceDof))

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
    Omega_s = Omega_s*2*pi #conversion from Hz to rad/s
    OmegaDot_s = OmegaDot_s*2*pi
    azi_sp1,Omega_sp1,OmegaDot_sp1,Fhat = timeIntegrateSubSystem(Irotor,Krotor,Crotor,Frotor, #time integrate using Newmark-Beta
    delta_t,azi_s,Omega_s,OmegaDot_s)

    Omega_sp1 = Omega_sp1/(2*pi) #convert to Hz, etc.
    OmegaDot_sp1 = OmegaDot_sp1/(2*pi)

    return azi_sp1,Omega_sp1,OmegaDot_sp1,Frotor

end

"""

    timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)

Internal, performs integration of a system using the Newmark-Beta method (constant-average acceleration sceheme).

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

    return unp1,udotnp1,uddotnp1,Fhat

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

    pred_vals = zeros(Float32, length(curr_vals[:,1]))

    # reduce times to relative values (fixes extrapolation problems when t gets large)
    t = ts .- ts[3]
    tout = t_out - ts[3]

    if interp_order == 0
        pred_vals = copy(curr_vals[:,end])
    elseif interp_order == 1
        k = tout / t[end-1]
        for i = 1:size(curr_vals, 1)
            pred_vals[i] = curr_vals[i,end] + (curr_vals[i,end-1] - curr_vals[i,end])*k
        end
    elseif interp_order == 2
        k = tout / (t[end-1]*t[end-2]*(t[end-1]-t[end-2]))
        for i = 1:size(curr_vals, 1)
            pred_vals[i] =   curr_vals[i,end] +
            (t[end-2]^2 * (curr_vals[i,end] - curr_vals[i,end-1]) + t[end-1]^2 * (-curr_vals[i,end] + curr_vals[i,end-2]))*k +
            ((t[end-1]-t[end-2])*curr_vals[i,end] + t[end-2]*curr_vals[i,end-1] - t[end-1]*curr_vals[i,end-2])*k*tout
        end
    else
        error("interp_order must equal 0, 1, or 2")
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
    out_frame_vals[1:3] = trans_mat * init_frame_vals[1:3]
    out_frame_vals[4:6] = trans_mat * init_frame_vals[4:6]

    return out_frame_vals

end

function calcHubRotMat(ptfmRot, azi_j)

    CN2P = transMat(ptfmRot[1], ptfmRot[2], ptfmRot[3])
    CP2H = [cos(azi_j) sin(azi_j) 0; -sin(azi_j) cos(azi_j) 0;0 0 1]

    CN2H = CN2P*CP2H

    return CN2H

end

"""

    calc_hydro_residual(new_accels, new_hydro_frcs, md_frc, u, FMultiplier)

Internal, finds the residual of the platform accelerations/forces when adding in new values
from HydroDyn, MoorDyn, and GyricFEA.

# Input
* `new_accels::Vector{<:float}`: the new platform accelerations in the inertial reference frame from GyricFEA  (m/s/s)
* `new_hydro_frcs::Vector{<:float}`: the new platform hydro forces in the inertial reference frame from HydroDyn (N)
* `md_frc::Vector{<:float}`: the mooring forces at the platform in the inertial reference frame from MoorDyn (N)
* `u::Vector{<:float}`: input vector containing platform loads and accelerations prior to the new additions
* `FMultiplier::Number`: scaling factor used to make the loads/accelerations in the same magnitude for the Jacobians

# Output
* `u_resid`: vector containing residual between the input u and the new load/accelerations
"""
function calcHydroResidual(accels_new, FHydro_new, FMooring_new, u, FMultiplier)

    Fnew = FMooring_new + FHydro_new

    resid = Vector{Float32}(undef, length(Fnew) + length(accels_new))
    resid[1:length(Fnew)] = u[1:length(Fnew)] - Fnew/FMultiplier
    resid[length(Fnew)+1:end] = u[length(Fnew)+1:end] - accels_new

    return resid
end
