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
* `deformAero::function`: deformAero(input) in work, currently just updates the backend aero state variables based on current RPM
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
* `epsilon_x_hist`: strain history for epsilon_x for each element at the 4 quad points
* `kappa_y_hist`: strain history for kappa_y for each element at the 4 quad points
* `kappa_z_hist`: strain history for kappa_z for each element at the 4 quad points
* `epsilon_z_hist`: strain history for epsilon_z for each element at the 4 quad points
* `kappa_x_hist`: strain history for kappa_x for each element at the 4 quad points
* `epsilon_y_hist`: strain history for epsilon_y for each element at the 4 quad points
"""
function Unsteady(model,feamodel,mesh,el,aero,deformAero;getLinearizedMatrices=false,assembly=nothing,system=nothing)

    # Declare Variable Type, are set later
    udot_j = 0.0
    uddot_j = 0.0
    torqueDriveShaft_j = 0.0
    Fhat = 0.0

    elStrain = fill(GyricFEA.ElStrain(zeros(4),zeros(4),zeros(4),zeros(4),zeros(4),zeros(4)), mesh.numEl)

    dispOut = GyricFEA.DispOut(elStrain,zeros(492,492),zeros(492,492),zeros(492,492))
    #................................................................

    ## Rotor mode initialization
    #..........................................................................
    OmegaInitial = model.OmegaInit #Initial rotor speed (Hz)

    if model.turbineStartup == 1 || model.turbineStartup == "forced" #forced start-up using generator as motor
        println("Running in forced starting mode.")
        model.generatorOn = true  #TODO: clean this redundant/conflicting logic up
        #     Omega = OmegaInitial
        rotorSpeedForGenStart = 0.0
    elseif model.turbineStartup == 2 || model.turbineStartup == "self" #self-starting mode
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
    Omega_jlast = 0.0
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

    if model.analysisType == "GX"
        npt = 1
        nel = length(assembly.start) #TODO: these should be the same.
    else
        npt = 4
        nel = mesh.numEl
    end
    epsilon_x_hist = zeros(npt,nel,numTS)
    kappa_y_hist = zeros(npt,nel,numTS)
    kappa_z_hist = zeros(npt,nel,numTS)
    epsilon_z_hist = zeros(npt,nel,numTS)
    kappa_x_hist = zeros(npt,nel,numTS)
    epsilon_y_hist = zeros(npt,nel,numTS)

    aziHist = zeros(numTS+1)
    OmegaHist = zeros(numTS+1)
    OmegaDotHist = zeros(numTS+1)
    FhatHist = zeros(numTS+1)
    gbHist = zeros(numTS+1)
    gbDotHist = zeros(numTS+1)
    gbDotDotHist = zeros(numTS+1)
    genTorque = zeros(numTS+1)
    genPower = zeros(numTS+1)
    torqueDriveShaft = zeros(numTS+1)
    Ywec = zeros(numTS+1,1)
    rigidDof = zeros(numTS+1) #TODO: figure this out

    t[1] = 0.0 #initialize various states and variables
    gb_s = 0
    gbDot_s = OmegaInitial
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

    if model.hydroOn
        # Spd = wave.resource.jonswap_spectrum(f=model.plat_model.hydro.freq, Tp=6, Hs=1)
    end

    Fexternal = 0.0 #TODO: do this right, especially if there are only hydro forces
    struct_X = 0.0
    struct_Y = 0.0
    struct_Z = 0.0
    Fdof = 0.0
    Fexternal_global = 0.0
    z3Dnorm = 0.0
    integrator = zeros(1) #for generator control algorithm
    ## Main Loop - iterate for a solution at each time step, i
    for i=1:numTS

        #     i #TODO add verbose printing
        # if (mod(i,100)==0) #print command that displays progress of time stepping
        println("Time Step: $i, time $(t[i])")
        # end

        ## check for specified rotor speed at t[i] + delta_t
        if (model.turbineStartup == 0)
            model.omegaControl = true #TODO: are we setting this back?
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
            # omegaCurrent = 0.0 TODO: figure this out
        end

        ## initialize "j" Gauss-Sidel iteration
        u_j=copy(u_s)
        azi_j = azi_s
        Omega_j = Omega_s
        OmegaDot_j = OmegaDot_s
        gb_j = gb_s
        gbDot_j = gbDot_s
        gbDotDot_j = gbDotDot_s
        genTorque_j = genTorque_s

        #initialize  platform module related variables only used if (model.hydroOn)
        Ywec_j = Ywec[i,:]
        Ywec_jLast = copy(Ywec_j)

        #TODO: put these in the model
        TOL = 1e-3  #gauss-seidel iteration tolerance for various modules
        MAXITER = 11 #max iteration for various modules
        numIterations = 1
        uNorm = 1e5
        platNorm = 0.0
        aziNorm = 1e5
        gbNorm = 0.0 #initialize norms for various module states

        ## evaluate platform module
        ##-------------------------------------
        if model.hydroOn
            # ds = model.plat_model.get_waveExcitation(Spd, time=[t[i],t[i]+delta_t], seed=1)
            # wave6dof_F_M = ds."fexc".data[1,:]
        end
        #-------------------------------------

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
        if model.analysisType=="GX"
            # systemout = deepcopy(system)
            strainGX = zeros(3,length(assembly.elements))
            curvGX = zeros(3,length(assembly.elements))
        end

        while ((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER)) #module gauss-seidel iteration loop
            # println("$(numIterations)   uNorm: $(uNorm)    platNorm: $(platNorm)    aziNorm: $(aziNorm)    gbNorm: $(gbNorm)")
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

            ## evaluate generator module
            #----- generator module ---------------------------
            genTorque_j = 0
            if model.generatorOn
                if model.useGeneratorFunction
                    specifiedOmega,_,_ = omegaSpecCheck(t[i]+delta_t,model.tocp,model.Omegaocp,delta_t)
                    genTorqueHSS0 = userDefinedGenerator(Omega_j,OmegaHist[i],delta_t,integrator,specifiedOmega)
                else
                    genTorqueHSS0 = simpleGenerator(model,Omega_j)
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
                    #
                    # if model.JgearBox==0.0
                    #     @error "model.JgearBox cannot be 0 if modeling the drivetrain. The SNL17m's was 243.0 n-s^2-m"
                    # end
                    # model.JgearBox = structureMOI[3,3]
                    println("here")
                    println(structureMOI[3,3])
                    println(-genTorque_j)
                    println(torqueDriveShaft_j)
                    gb_j,gbDot_j,gbDotDot_j,Fhat = updateRotorRotation(model.JgearBox,0,0,
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
                    azi_j,Omega_j,OmegaDot_j = getRotorPosSpeedAccelAtTime(t[i],t[i]+delta_t,azi_s,delta_t)
                else
                    Omega_j = Omega_s
                    OmegaDot_j = OmegaDot_s
                    azi_j = azi_s + Omega_j*delta_t*2*pi
                end
            elseif !model.omegaControl
                Crotor = 0
                Krotor = 0
                azi_j,Omega_j,OmegaDot_j,Fhat = updateRotorRotation(structureMOI[3,3],Crotor,Krotor,
                -FReaction_j[6],-torqueDriveShaft_j,
                azi_s,Omega_s,OmegaDot_s,delta_t)
            else
                error("omega control option not correctly specified")
            end

            ## compile external forcing on rotor
            #compile forces to supply to structural dynamics solver

            if model.aeroLoadsOn>0 #0 off, 1 one way, 2 two way
                runaero = true
                if model.aeroLoadsOn==1 && numIterations!=1
                    runaero = false
                end
                if runaero
                    if model.tocp_Vinf == -1
                        newVinf = -1
                    else
                        newVinf = FLOWMath.akima(model.tocp_Vinf,model.Vinfocp,t[i])
                    end

                    if model.aeroLoadsOn==2
                        # Transform Local Displacements to Global
                        numDofPerNode = 6
                        #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
                        u_j_global = zeros(Int(max(maximum(mesh.structuralNodeNumbers))*6))
                        for jbld = 1:length(mesh.structuralElNumbers[:,1])
                            for kel = 1:length(mesh.structuralElNumbers[1,:])-1
                                # orientation angle,xloc,sectionProps,element order]
                                elNum = Int(mesh.structuralElNumbers[jbld,kel])
                                #get dof map
                                node1 = Int(mesh.structuralNodeNumbers[jbld,kel])
                                node2 = Int(mesh.structuralNodeNumbers[jbld,kel+1])
                                dofList = [(node1-1)*numDofPerNode.+(1:6);(node2-1)*numDofPerNode.+(1:6)]

                                localdisp = u_j[dofList]
                                x = [mesh.x[node1], mesh.x[node2]]
                                y = [mesh.y[node1], mesh.y[node2]]
                                z = [mesh.z[node1], mesh.z[node2]]

                                twist = el.props[elNum].twist
                                sweepAngle = el.psi[elNum]
                                coneAngle = el.theta[elNum]
                                rollAngle = el.roll[elNum]

                                twistAvg = rollAngle + 0.5*(twist[1] + twist[2])
                                lambda = GyricFEA.calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)
                                globaldisp = inv(lambda')*localdisp

                                #asssembly
                                for m = 1:length(dofList)
                                    u_j_global[dofList[m]] =  u_j_global[dofList[m]]+globaldisp[m]
                                end

                            end
                        end

                        disp_x = [u_j_global[i] for i = 1:6:length(u_j_global)]
                        disp_y = [u_j_global[i] for i = 2:6:length(u_j_global)]
                        disp_z = [u_j_global[i] for i = 3:6:length(u_j_global)]
                        disp_twist = [u_j[i] for i = 4:6:length(u_j)]

                        bld_x = zero(mesh.structuralElNumbers[:,1:end-1])
                        bld_y = zero(mesh.structuralElNumbers[:,1:end-1])
                        bld_z = zero(mesh.structuralElNumbers[:,1:end-1])
                        bld_twist = zero(mesh.structuralElNumbers[:,1:end-1])

                        for jbld = 1:length(mesh.structuralElNumbers[:,1])
                            bld_indices = Int.(mesh.structuralElNumbers[jbld,1:end-1])
                            bld_x[jbld,:] = mesh.x[bld_indices]+disp_x[bld_indices]
                            bld_y[jbld,:] = mesh.y[bld_indices]+disp_y[bld_indices]
                            bld_z[jbld,:] = mesh.z[bld_indices]+disp_z[bld_indices]
                            # flatten blade x,y
                            bld_x[jbld,:] = sqrt.(bld_x[jbld,:].^2 .+bld_y[jbld,:].^2)
                            bld_twist[jbld,:] = disp_twist[bld_indices]
                        end
                    else
                        bld_x = -1
                        bld_z = -1
                        bld_twist = -1
                    end

                    println("Calling Aero $(Omega_j*60) RPM $newVinf Vinf")
                    deformAero(azi_j;newOmega=Omega_j*2*pi,newVinf,bld_x,bld_z,bld_twist) #TODO: implement deformation induced velocities
                    Fexternal, Fdof,Xp,Yp,Zp,z3Dnorm,Fexternal_global = aero(t[i],azi_j)
                end
            end

            if model.hydroOn
                Fdof = [Fdof; Int.(Fdof[:,1]); collect(1:6)] #TODO: tie into ndof per node
                Fexternal = [Fexternal; wave6dof_F_M]
            end

            ## evaluate structural dynamics
            #initialization of structural dynamics displacements, velocities, accelerations, etc.

            if model.analysisType=="ROM"
                dispData = GyricFEA.DispData(u_s,udot_s,uddot_s,u_sm1,eta_s,etadot_s,etaddot_s)
            else
                dispData = GyricFEA.DispData(u_s,udot_s,uddot_s,u_sm1)
            end

            if model.analysisType=="ROM" # evalulate structural dynamics using reduced order model
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransientROM(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,rom,Fexternal,Int.(Fdof),CN2H,rbData)
            elseif model.analysisType=="GX"

                # function getForces(Fexternal_global,Fdof,node,dof)
                #
                #     gdof = (node-1)*6+dof
                #     idxGdof = findfirst(x->x==gdof,Fdof[:])
                #
                #     if !isnothing(idxGdof)
                #         # println("1 $(size(Fexternal_global))")
                #         # println("2 $(idxGdof)")
                #         return Fexternal_global[idxGdof]
                #     else
                #         return nothing
                #     end
                # end
                #
                # # create dictionary of prescribed conditions
                #
                # prescribed_conditions = Dict()
                # for inode = 1:length(assembly.start)
                #
                #     prescribed_conditions[assembly.start[inode]] = GXBeam.PrescribedConditions(Fx = getForces(Fexternal_global,Fdof,assembly.start[inode],1),
                #     Fy = getForces(Fexternal_global,Fdof,assembly.start[inode],2), Fz = getForces(Fexternal_global,Fdof,assembly.start[inode],3),
                #     Mx = getForces(Fexternal_global,Fdof,assembly.start[inode],4), My = getForces(Fexternal_global,Fdof,assembly.start[inode],5),
                #     Mz = getForces(Fexternal_global,Fdof,assembly.start[inode],6))
                # end
                # prescribed_conditions[1] = GXBeam.PrescribedConditions(ux=0, uy=0, uz=0,theta_x=0, theta_y=0, theta_z=0)
                # top_idx = Int(feamodel.joint[7,2])
                # prescribed_conditions[top_idx] = GXBeam.PrescribedConditions(ux=0, uy=0)
                # #TODO: parametric way to pull in pBC

                distributed_loads = Dict()
                # B1
                # GXz = [assembly.points[ipt][3] for ipt = 24:46]
                for jbld = 1:length(mesh.structuralElNumbers[:,1])
                    XpGXspl1 = FLOWMath.Akima(z3Dnorm,Xp[jbld,:,1])#,GXz/maximum(GXz))
                    YpGXspl1 = FLOWMath.Akima(z3Dnorm,Yp[jbld,:,1])#,GXz/maximum(GXz))
                    ZpGXspl1 = FLOWMath.Akima(z3Dnorm,Zp[jbld,:,1])#,GXz/maximum(GXz))

                    for ipt = mesh.structuralNodeNumbers[jbld,1]:mesh.structuralNodeNumbers[jbld,end-1] #TODO: el or node?
                        iel = findfirst(x->x==ipt,assembly.start)
                        if !isempty(iel)
                            distributed_loads[ipt] = GXBeam.DistributedLoads(assembly,iel;fx = (s) -> XpGXspl1(s),
                                fy = (s) -> YpGXspl1(s), fz = (s) -> ZpGXspl1(s))
                        else
                            println("Empty at $ipt")
                        end
                    end
                end

                # create dictionary of prescribed conditions

                # prescribed_conditions = Dict()
                # # for ibld = 1:length(mesh.structuralNodeNumbers[:,1])
                # #     for (inode,nodenum) in enumerate(Int.(mesh.structuralNodeNumbers[ibld,:]))
                # #         L = assembly.elements[nodenum].L
                # #         prescribed_conditions[nodenum] = GXBeam.PrescribedConditions(Fx = -struct_X[ibld,1,inode]*L,
                # #         Fy = -struct_Y[ibld,1,inode]*L, Fz = struct_Z[ibld,1,inode]*L)
                # #     end
                # # end
                # prescribed_conditions[1] = GXBeam.PrescribedConditions(ux=0, uy=0, uz=0,theta_x=0, theta_y=0, theta_z=0)
                # top_idx = Int(feamodel.joint[7,2])
                # prescribed_conditions[top_idx] = GXBeam.PrescribedConditions(ux=0, uy=0)
                # #TODO: parametric way to pull in pBC
                #
                # # Create Dictionary of Distributed Loads
                # distributed_loads = Dict()
                # for ibld = 1:length(mesh.structuralNodeNumbers[:,1])
                #     for (iel,elnum) in enumerate(Int.(mesh.structuralNodeNumbers[ibld,:]))
                #         distributed_loads[elnum] = GXBeam.DistributedLoads(assembly,iel;fx = (s) -> -struct_X[ibld,1,iel],
                #         fy = (s) -> -struct_Y[ibld,1,iel], fz = (s) -> struct_Z[ibld,1,iel])
                #     end
                # end


                prescribed_conditions = (t) -> begin
                    Dict(
                    # fixed base
                    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
                    # fixed top, but free to rotate around z-axis
                    # 50 => GXBeam.PrescribedConditions(Fx = 1e4*sin(20*t)),
                    23 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0),
                    )
                end

                if i == 1 && numIterations == 1
                    reset_state = false
                    initialize = false
                else
                    reset_state = false
                    initialize = false
                end
                linear_velocity = [0.0,0.0,0.0]
                angular_velocity = [0.0,0.0,Omega_j*2*pi]
                linear_acceleration = [0.0,0.0,0.0]
                angular_acceleration = [0.0,0.0,OmegaDot_j*2*pi]

                tvec = [t[i],t[i]+delta_t]

                gravity = [0.0,0.0,-9.81]

                systemout, history, converged = GXBeam.time_domain_analysis!(deepcopy(system),assembly, tvec;
                reset_state,initialize,linear_velocity,angular_velocity,linear_acceleration,
                angular_acceleration,prescribed_conditions,distributed_loads,gravity,linear=false)#!feamodel.nlOn)

                if !converged
                    println("GX Didn't Converge")
                end
                # elStrain
                state = GXBeam.AssemblyState(systemout, assembly;
                prescribed_conditions)

                for iel = 1:length(state.elements)
                    strainGX[:,iel] = GXBeam.element_strain(assembly.elements[iel],state.elements[iel].F,state.elements[iel].M)
                    curvGX[:,iel] = GXBeam.element_curvature(assembly.elements[iel],state.elements[iel].F,state.elements[iel].M)
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
                dispOut = GyricFEA.DispOut(nothing, disp_sp1,zero(disp_sp1),zero(disp_sp1))
                FReaction_j = [-history[end].points[1].F[1];-history[end].points[1].F[2];-history[end].points[1].F[3];
                -history[end].points[1].M[1];-history[end].points[1].M[2];-history[end].points[1].M[3]]
            else # evalulate structural dynamics using conventional representation
                elStrain,dispOut,FReaction_j = GyricFEA.structuralDynamicsTransient(feamodel,mesh,el,dispData,Omega_j,OmegaDot_j,t[i],delta_t,elStorage,Fexternal,Int.(Fdof),CN2H,rbData)
            end

            #update last iteration displacement vector
            u_jLast = u_j
            u_j = dispOut.displ_sp1             #update current estimates of velocity, acceleration
            # Only used for TNB, but must be declared
            udot_j  = dispOut.displdot_sp1
            uddot_j = dispOut.displddot_sp1

            ## calculate norms
            uNorm = LinearAlgebra.norm(u_j-u_jLast)/LinearAlgebra.norm(u_j)            #structural dynamics displacement iteration norm
            aziNorm = LinearAlgebra.norm(azi_j - azi_jLast)/LinearAlgebra.norm(azi_j)  #rotor azimuth iteration norm
            platNorm = LinearAlgebra.norm(Ywec_j-Ywec_jLast)/LinearAlgebra.norm(Ywec_j) #platform module states iteration norm if it is off, the norm will be zero
            gbNorm = LinearAlgebra.norm(gb_j - gb_jLast)/LinearAlgebra.norm(gb_j) #gearbox states iteration norm if it is off, the norm will be zero
            Omega_jlast = Omega_j

            if model.analysisType=="GX" && !((uNorm > TOL || platNorm > TOL || aziNorm > TOL || gbNorm > TOL) && (numIterations < MAXITER))
                system = deepcopy(systemout)
            end
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
                genTorquePlot = simpleGenerator(model,gbDot_j*model.gearRatio)
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

        uHist[1:length(u_s),i+1] = u_s #TODO: resolve this for GX
        FReactionHist[i+1,:] = FReaction_j
        FhatHist[i+1] = Fhat
        if model.analysisType=="GX"
            for ii = 1:length(strainGX[1,:])
                epsilon_x_hist[:,ii,i] .= strainGX[1,ii]
                kappa_y_hist[:,ii,i] .= curvGX[2,ii]
                kappa_z_hist[:,ii,i] .= curvGX[3,ii]
                epsilon_z_hist[:,ii,i] .= strainGX[3,ii]
                kappa_x_hist[:,ii,i] .= curvGX[1,ii]
                epsilon_y_hist[:,ii,i] .= strainGX[2,ii]
            end
        else
            for ii = 1:length(elStrain)
                epsilon_x_hist[:,ii,i] = elStrain[ii].epsilon_x
                kappa_y_hist[:,ii,i] = elStrain[ii].kappa_y
                kappa_z_hist[:,ii,i] = elStrain[ii].kappa_z
                epsilon_z_hist[:,ii,i] = elStrain[ii].epsilon_z
                kappa_x_hist[:,ii,i] = elStrain[ii].kappa_x
                epsilon_y_hist[:,ii,i] = elStrain[ii].epsilon_y
            end
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

        genTorque[i+1] = genTorque_s
        # genTorque[i+1] = genTorquePlot
        genPower[i+1] = genPowerPlot
        torqueDriveShaft[i+1] = torqueDriveShaft_s

        # Ywec(i+1,:) = Ywec_j #TODO figure out platform deflection, etc
        # rigidDof(i+1,:)=Ywec(i+1,s_stsp+1:s_stsp+6)

        ## check rotor speed for generator operation
        if Omega_s >= rotorSpeedForGenStart
            model.generatorOn = true
        else
            model.generatorOn = false
        end

    end #end timestep loop

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
            HDF5.write(file,"rigidDof",rigidDof)
            HDF5.write(file,"genTorque",genTorque)
            HDF5.write(file,"genPower",genPower)
            HDF5.write(file,"torqueDriveShaft",torqueDriveShaft)
            HDF5.write(file,"uHist",uHist)
            HDF5.write(file,"eps_xx_0_hist",epsilon_x_hist)
            HDF5.write(file,"eps_xx_z_hist",kappa_y_hist)
            HDF5.write(file,"eps_xx_y_hist",kappa_z_hist)
            HDF5.write(file,"gam_xz_0_hist",epsilon_z_hist)
            HDF5.write(file,"gam_xz_y_hist",kappa_x_hist)
            HDF5.write(file,"gam_xy_0_hist",epsilon_y_hist)
            HDF5.write(file,"gam_xy_z_hist",-kappa_x_hist)
        end

    end
    return t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,rigidDof,genTorque,genPower,torqueDriveShaft,uHist,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,kappa_x_hist,epsilon_y_hist,-kappa_x_hist,FhatHist
end

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

function userDefinedGenerator(omega,omegalast,dt,integrator,omega0)
    # omega is in hz
    omega_RPM = omega*60
    omega_RPM0 = omega0*60 # 33.92871
    Kp = 11.1599
    Ki = 5.2693
    Kd = 3.3076
    Q0 = -345.7731*0
    integrator[1] = integrator[1] + (omega_RPM - omega_RPM0)*dt
    deriv = (omega-omegalast)/dt
    controllerQ = Q0 + Kp*(omega_RPM-omega_RPM0) + Kd*deriv + Ki*integrator[1]
    return controllerQ*1000
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

    # Fexternal = zeros(length(ForceDof))
    #
    # for i = 1:length(ForceDof)
    #     Fexternal[i] = FLOWMath.linear(timeArray,ForceValHist[i,:],time)
    # end
    # Fdof = ForceDof

    if (time < 0.5)
        Fexternal = [1e6,1e6,1e6];
        Fdof = [24*6-5,24*6-4,24*6];
    else
        Fexternal = [];
        Fdof = [];
    end

    return Fexternal, Fdof,0,0,0,0,0
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
