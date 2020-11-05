function mapACloads(u_jLast,udot_j,Omega_j,t,PEy,QCy,NElem,NBlade,RefR,mesh,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,el,turbine3D,env,step_AC,us_param)

    N_blade_nodes = length(structuralSpanLocNorm[1,:])+1
    # Initialize bladeForce
    bladeForce_N = zeros(NBlade,NElem)
    bladeForce_T = zeros(NBlade,NElem)
    bladeForce_M25 = zeros(NBlade,NElem)

    for k = 1:length(turbine3D)
        #TODO: Incorporate deflections and changes in omega ->
        #r,twist,delta,omega all need to be vectors aligning with the ntheta
        #discretizations of the cylinder

        #TODO: ensure that the deflections aren't compounding after a
        #revolution.  They may be wrong

        #TODO: Verify units everywhere

        circ_step_num = floor((step_AC-1)/turbine3D[k].ntheta*turbine3D[k].B)
        circular_step = step_AC-circ_step_num*turbine3D[k].ntheta/turbine3D[k].B
        idx_sub = Int.(collect(circular_step:turbine3D[k].ntheta/turbine3D[k].B:turbine3D[k].ntheta-turbine3D[k].ntheta/turbine3D[k].B+1+circular_step))

        #TODO: this is hard coded for 2 blades, need to simplify
        # Interpolate the deformations onto the aero model for the current step
        norm_disp_h = LinRange(0,1,N_blade_nodes)
        # 1 = Z deformation - not modeled in 2D AC method
        # 2 = Tangential deformation - no real effect on AC model
        offset = 3
        turbine3D[k].r[idx_sub[1]] = turbine3D[k].r[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        turbine3D[k].r[idx_sub[2]] = turbine3D[k].r[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        offset = 4
        turbine3D[k].twist[idx_sub[1]] = turbine3D[k].twist[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        turbine3D[k].twist[idx_sub[2]] = turbine3D[k].twist[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        offset = 5
        turbine3D[k].delta[idx_sub[1]] = turbine3D[k].delta[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        turbine3D[k].delta[idx_sub[2]] = turbine3D[k].delta[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        # 6 = Sweep deformation, not modeled in AC method - assuming it is small so that it doesn't spill over into the next step/theta discretization

        turbine3D[k].omega[:] .= Omega_j

        # Interpolate deformation induced velocities onto the aero model for the most current step
        offset = 1
        env[k].V_vert[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        env[k].V_vert[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        offset = 2
        env[k].V_tang[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        env[k].V_tang[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        offset = 3
        env[k].V_rad[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        env[k].V_rad[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        offset = 4
        env[k].V_twist[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
        env[k].V_twist[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
        # 5 = Change in delta angle, not modeled in 2D AC method
        # 6 = Change in sweep angle, not modeled in 2D AC method

        # envin = env[k]
        # turbine_in = turbine3D[k]
        # mat"[$Rp, $Tp, $Zp, $Mp, $envout] = actuatorcylinder_substep($turbine_in, $envin,$us_param, $step_AC,$alpha_rad,$cl_af,$cd_af)"
        # env[k].rho = envout["rho"]
        # env[k].mu = envout["mu"]
        # env[k].G_amp = envout["G_amp"]
        # env[k].gusttime = envout["gusttime"]
        # env[k].gustX0 = envout["gustX0"]
        # env[k].N_Rev = envout["N_Rev"]
        # env[k].idx_sub = envout["idx_sub"]
        # env[k].wsave = envout["wsave"]
        # env[k].Vinf_nominal = envout["Vinf_nominal"]
        # env[k].V_vert = envout["V_vert"]
        # env[k].V_tang = envout["V_tang"]
        # env[k].V_rad = envout["V_rad"]
        # env[k].V_twist = envout["V_twist"]
        # env[k].Vinf = envout["Vinf"]
        # env[k].V_wake_old = envout["V_wake_old"]
        # env[k].steplast = envout["steplast"]
        # turbine3D[k].ntheta = envout["ntheta"]
        # env[k].tau = envout["tau"]

        Q, Rp, Tp, Zp, Vinf_used, alpha, cl, cd, Vloc, Re = VAWTAero.Unsteady_Step(turbine3D[k],env[k],us_param,step_AC)
        Mp = zeros(length(Zp)) #TODO: fix moment CALCS in VAWTAero



        for j=1:NBlade
            delta = turbine3D[k].delta[idx_sub[j]]
            bladeForce_N[j,k] = -Rp[j]*cos(delta) + -Zp[j]*sin(delta)
            bladeForce_T[j,k] = -Tp[j] #TODO: fix RPI's difficulty to converge when running in reverse (may have to change RPI indexing)
            bladeForce_M25[j,k] = Mp[j]
        end
    end



    # scatter(t,bladeForce[1].T[floor(blade[j].NElem/2)])
    # hold on
    # pause[0.001)
    #define these from params file
    ft2m = 1 / 3.281

    #     RefAR = cactusGeom.RefAR*ft2m*ft2m
    RefR = RefR*ft2m

    spanLocNorm = zeros(NBlade,NElem)
    for i=1:NBlade
        spanLocNorm[i,:] = PEy[1:NElem[1,1],1].*RefR[1,1]/(QCy[NElem[1,1]+1,1]*RefR[1,1])
    end

    #Initialize structuralLoad

    structuralLoad_N = zeros(NBlade,length(structuralElNumbers[1,:]))
    structuralLoad_T = zeros(NBlade,length(structuralElNumbers[1,:]))
    structuralLoad_M25 = zeros(NBlade,length(structuralElNumbers[1,:]))

    for i=1:NBlade
        structuralLoad_N[i,:] = FLOWMath.linear(spanLocNorm[i,:],bladeForce_N[i,:],structuralSpanLocNorm[i,:])
        structuralLoad_T[i,:] = FLOWMath.linear(spanLocNorm[i,:],bladeForce_T[i,:],structuralSpanLocNorm[i,:])
        structuralLoad_M25[i,:]= FLOWMath.linear(spanLocNorm[i,:],bladeForce_M25[i,:],structuralSpanLocNorm[i,:])
    end

    _,numNodesPerBlade = size(structuralNodeNumbers)

    #integrate over elements

    #read element data in

    numDofPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6))
    for j = 1:NBlade
        for k = 1:numNodesPerBlade-1
            #get element data
            # orientation angle,xloc,sectionProps,element order]
            elNum = Int(structuralElNumbers[j,k])
            #get dof map
            node1 = Int(structuralNodeNumbers[j,k])
            node2 = Int(structuralNodeNumbers[j,k+1])
            dofList = [(node1-1)*numDofPerNode.+(1:6), (node2-1)*numDofPerNode.+(1:6)]

            elementOrder = 1
            x = [mesh.x[node1], mesh.x[node2]]
            elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
            xloc = [0 elLength]
            twist = el.props[elNum].twist
            sweepAngle = el.psi[elNum]
            coneAngle = el.theta[elNum]
            rollAngle = el.roll[elNum]

            extDistF2Node =  [structuralLoad_T[j,k],   structuralLoad_T[j,k+1]]
            extDistF3Node = -[structuralLoad_N[j,k],   structuralLoad_N[j,k+1]]
            extDistF4Node = -[structuralLoad_M25[j,k], structuralLoad_M25[j,k+1]]

            mat"[$Fe] = calculateLoadVecFromDistForce($elementOrder,$x,$xloc,$twist,$sweepAngle,$coneAngle,$rollAngle,$extDistF2Node,$extDistF3Node,$extDistF4Node)"

            #asssembly
            for m = 1:length(dofList)
                Fg[dofList[m]] =  Fg[dofList[m]].+Fe[m]
            end

        end
    end

    ForceDof = Float64.(1:length(Fg))

    return Fg,ForceDof,env

end


function calculateStructureMassProps(elStorage)
    #calculateStructureMassProps   calculates mass properties of mesh
    #   [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage)
    #
    #   This function caclulates structural mass properties of the finite
    #   element mesh (mass, moment of inertia, mass center) about the origin of
    #   the mesh coordinate system.
    #
    #      input:
    #      elStorage    = object containing arrays of stored element
    #                     information
    #
    #      output:
    #      structureMass       = mass of structure
    #      structureMOI        = moment of inertia tensor of structgure
    #      structureMassCenter = center of mass of structure

    numElements = length(elStorage) #get number of elements

    structureMass = 0.0 #initialize structure mass and moment of inertia
    structureMOI = zeros(3,3)
    temp = zeros(3,1)
    for i=1:numElements #sum over elemetns contribution to mass and moment of inertia
        structureMass = structureMass + elStorage[i].mel
        structureMOI = structureMOI + elStorage[i].moiel
        temp = temp + elStorage[i].xmel
    end

    structureMassCenter = temp./structureMass #calculate mass center

    #modify moment of inertia to be about structure mass center
    x = structureMassCenter[1]
    y = structureMassCenter[2]
    z = structureMassCenter[3]

    structureMOI = structureMOI - structureMass*[(y^2+z^2) -x*y -x*z
    -x*y (x^2+z^2) -y*z
    -x*z -y*z (x^2+y^2)]


    return structureMass,structureMOI,structureMassCenter

end

function calculateLambda(theta1,theta2,theta3)
    #calculateLambda Calculates transformation matrix from element to hub frame
    #   [lambda] = calculateLambda(theta1,theta2,theta3 )
    #
    #   This function calculates a transformation matrix to transform the
    #   element degree of freedom vector (12 DOFs) from the hub frame to
    #   the element frame. The transformation matrix is constructed via the
    #   direction cosine matrices of a 3-2-1 Euler rotation sequence.
    #
    #      input:
    #      theta1        = angle (rad) of rotation for 1st rotation
    #                      of 3-2-1 sequence
    #      theta2        = angle (rad) of rotation for 2nd rotation
    #                      of 3-2-1 sequence
    #      theta3        = angle (rad) of rotation for 3rd rotation
    #                      of 3-2-1 sequence

    #      output:
    #      lambda        = 12 x 12 transformation matrix

    # dcm that is created is [dcm] = [M1(theta3)][M2(theta2)][M3(theta1)]

    ct1 = cos(theta1); st1=sin(theta1);
    ct2 = cos(theta2); st2=sin(theta2);
    ct3 = cos(theta3); st3=sin(theta3);

    fac1 = st3*st2
    fac2 = ct3*st2
    dcm = [ct2*ct1           ct2*st1          -st2
    fac1*ct1-ct3*st1  fac1*st1+ct3*ct1  st3*ct2
    fac2*ct1+st3*st1  fac2*st1-st3*ct1  ct3*ct2]


    lambda = zeros(12,12)
    lambda[1:3,1:3] = dcm
    lambda[4:6,4:6] = dcm
    lambda[7:9,7:9] = dcm
    lambda[10:12,10:12] = dcm

    return lambda

end

function createTda(jointType,slaveNodeNum,masterNodeNum,psi,theta,joint)
    #This function creates a constraint transformation matrix for a single
    #joint. Tda is this matrix, dDOF contains a listing of dependent global
    #DOFs associated with this joint, and aDOF contains a listing of active
    #global DOFs associated with this joint.

    if (jointType == 4 && (abs(abs(theta)-90)<1.0e-3 || (abs(abs(theta)-270)<1.0e-3) ))
        theta = 0.0
        jointType = 3
    end

    #calculate transformation matrix from hub frame to joint frame
    Lambda = calculateLambda(psi*pi/180.0,theta*pi/180.0,0.0)

    #Tda is a local mapping of dependent DOFs to active DOFs at a node
    # u_d = Tda * u_a
    # such that u_d is a list of local dependent slave DOFs at a jont and
    # u_a is a list of local dependent slave DOFs at a joint.


    if (jointType == 0) #for weld/fixed joint type
        activeDof0 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof0 = [1 2 3 4 5 6] #slave DOF list at joint

        #determine local active DOFs associated with slave node
        slaveActiveDof0 = determineActiveDofsFromSlaveNode(slaveDof0,6)

        Rda0 = -1.0*LinearAlgebra.I(6) #from constraint equation for fixed joint
        Rdd0 = 1.0*LinearAlgebra.I(6)

        Tda,dDOF,aDOF = getNodeMaps(Rdd0,Rda0,masterNodeNum,slaveNodeNum,slaveDof0,activeDof0,slaveActiveDof0)

    elseif (jointType == 1) #for pinned joint type
        activeDof1 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof1 = [1 2 3] #slave DOF list at joint

        #determine local active DOFs associated with slave node
        slaveActiveDof1 = determineActiveDofsFromSlaveNode(slaveDof1,6)

        Rda1 = -[1.0*LinearAlgebra.I(3), zeros(3,6)] #from constraint equation for pinned joint
        Rdd1 = 1.0*LinearAlgebra.I(3)

        Tda,dDOF,aDOF = getNodeMaps(Rdd1,Rda1,masterNodeNum,slaveNodeNum,slaveDof1,activeDof1,slaveActiveDof1)

    elseif (jointType == 2)     #hinge axis along localy "2" frame of joint

        if ((abs(abs(psi)-90))<1.0e-3 || (abs(abs(psi)-270))<1.0e-3)
            activeDof2 = [1 2 3 4 5 6]
            slaveDof2  = [1 2 3 5 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof2 = determineActiveDofsFromSlaveNode(slaveDof2,6)

            globalConstraintEqMatrix2 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3) zeros(3,3)
            zeros(2,3) -[0 1 0;0 0 1] zeros(2,3) [0 1 0;0 0 1]]

        else
            activeDof2 = [1 2 3 4 5 6]
            slaveDof2 = [1 2 3 4 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof2 = determineActiveDofsFromSlaveNode(slaveDof2,6)

            localConstraintEqMatrix2 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 0 1] zeros(2,3) [1 0 0;0 0 1]]
            globalConstraintEqMatrix2 = localConstraintEqMatrix2*Lambda
        end
        #extract Rda from globalConstraintEqMatrix2
        Rda2 = zeros(length(globalConstraintEqMatrix2[:,1]),length(activeDof2))
        for i=1:length(activeDof2)
            ind = activeDof2[i]
            Rda2[:,i] = globalConstraintEqMatrix2[:,ind]
        end

        for i=1:length(slaveActiveDof2)
            ind = slaveActiveDof2[i]+6
            Rda2[:,i] = globalConstraintEqMatrix2[:,ind]
        end

        #extract Rdd from globalConstraintEqMatrix2
        Rdd2 = zeros(length(globalConstraintEqMatrix2[:,1]),length(slaveDof2))
        for i=1:length(slaveDof2)
            ind = slaveDof2[i]+6
            Rdd2[:,i] = globalConstraintEqMatrix2[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd2,Rda2,masterNodeNum,slaveNodeNum,slaveDof2,activeDof2,slaveActiveDof2)

    elseif (jointType == 3)     #hinge axis along local "1" frame of joint

        if ((abs(abs(theta)-90))<1.0e-3 || (abs(abs(theta)-270))<1.0e-3)
            activeDof3 = [1 2 3 4 5 6]
            slaveDof3  = [1 2 3 4 5]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            globalConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 1 0] zeros(2,3) [1 0 0;0 1 0]]

        elseif ((abs(abs(psi)-90))<1.0e-3 || (abs(abs(psi)-270))<1.0e-3)

            activeDof3 = [1 2 3 4 5 6]
            slaveDof3  = [1 2 3 4 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            globalConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[1 0 0;0 0 1] zeros(2,3) [1 0 0;0 0 1]]

        else

            activeDof3 = [1 2 3 4 5 6]
            slaveDof3 = [1 2 3 5 6]

            #determine local active DOFs associated with slave node
            slaveActiveDof3 = determineActiveDofsFromSlaveNode(slaveDof3,6)

            localConstraintEqMatrix3 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
            zeros(2,3) -[0 1 0;0 0 1] zeros(2,3) [0 1 0;0 0 1]]
            globalConstraintEqMatrix3 = localConstraintEqMatrix3*Lambda

        end

        #extract Rda from globalConstraintEqMatrix3
        Rda3 = zeros(length(globalConstraintEqMatrix3[:,1]),length(activeDof3))
        for i=1:length(activeDof3)
            ind = activeDof3[i]
            Rda3[:,i] = globalConstraintEqMatrix3[:,ind]
        end

        for i=1:length(slaveActiveDof3)
            ind = slaveActiveDof3[i]+6
            Rda3[:,i] = globalConstraintEqMatrix3[:,ind]
        end

        #extract Rdd from globalConstraintEqMatrix3
        Rdd3 = zeros(length(globalConstraintEqMatrix3[:,1]),length(slaveDof3))
        for i=1:length(slaveDof3)
            ind = slaveDof3[i]+6
            Rdd3[:,i] = globalConstraintEqMatrix3[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd3,Rda3,masterNodeNum,slaveNodeNum,slaveDof3,activeDof3,slaveActiveDof3)

    elseif (jointType == 4)     #hinge axis along local "3" frame of joint

        activeDof4 = [1 2 3 4 5 6]
        slaveDof4  = [1 2 3 4 5]

        #determine local active DOFs associated with slave node
        slaveActiveDof4 = determineActiveDofsFromSlaveNode(slaveDof4,6)

        localConstraintEqMatrix4 = [-1.0*LinearAlgebra.I(3) zeros(3,3) 1.0*LinearAlgebra.I(3)  zeros(3,3)
        zeros(2,3) -[1 0 0;0 1 0] zeros(2,3) [1 0 0;0 1 0]]
        globalConstraintEqMatrix4 = localConstraintEqMatrix4*Lambda

        #extract Rda from globalConstraintEqMatrix4
        Rda4 = zeros(length(globalConstraintEqMatrix4[:,1]),length(activeDof4))
        for i=1:length(activeDof4)
            ind = activeDof4[i]
            Rda4[:,i] = globalConstraintEqMatrix4[:,ind]
        end

        for i=1:length(slaveActiveDof4)
            ind = slaveActiveDof4[i]+6
            Rda4[:,i] = globalConstraintEqMatrix4[:,ind]
        end

        #extract Rdd from globalConstraintEqMatrix4
        Rdd4 = zeros(length(globalConstraintEqMatrix4[:,1]),length(slaveDof4))
        for i=1:length(slaveDof4)
            ind = slaveDof4[i]+6
            Rdd4[:,i] = globalConstraintEqMatrix4[:,ind]
        end

        Tda,dDOF,aDOF = getNodeMaps(Rdd4,Rda4,masterNodeNum,slaveNodeNum,slaveDof4,activeDof4,slaveActiveDof4)


    elseif (jointType == 5) #rigid bar constraint type
        activeDof5 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof5 = [1 2 3 4 5 6] #slave DOF list at joint
        #determine local active DOFs associated with slave node
        slaveActiveDof5 = determineActiveDofsFromSlaveNode(slaveDof5,6)

        Rdd5 = 1.0*LinearAlgebra.I(6);  #need to define lx,ly,lz, from mesh level
        lx = joint[5]
        ly = joint[6]
        lz = joint[7]

        Rda5 = -1.0*LinearAlgebra.I(6)
        Rda5[1:3,4:6] = [0 -lz ly;lz 0 -lx;-ly lx 0]

        Tda,dDOF,aDOF = getNodeMaps(Rdd5,Rda5,masterNodeNum,slaveNodeNum,slaveDof5,activeDof5,slaveActiveDof5);

    else
        error("Correct jointType not specified, should be 1, 2, 3, 4, or 5")
    end
    return Tda,dDOF,aDOF
end

function getNodeMaps(Rdd,Rda,masterNodeNum,slaveNodeNum,slaveDof,activeDof,slaveActiveDof)

    if (abs(LinearAlgebra.det(Rdd)) < 1.0e-3)
        error("Singular joint transformation matrix. Exiting")
    end

    Tda = -Rdd\Rda #calculate Tda #TODO

    numSlaveDOFs = length(slaveDof) #get number of joint DOFs for this joint
    numActiveDOFsFromMasterNode = length(activeDof) #get number of active DOFs for this joint

    dDOF = zeros(Int,numSlaveDOFs) #initialize arrays
    aMap = zeros(Int,numActiveDOFsFromMasterNode,1)

    for i=1:numSlaveDOFs
        #get global DOF numbers of slave DOFs for this joint
        dDOF[i,1] = (slaveNodeNum-1)*6 + slaveDof[i]
    end

    for i=1:numActiveDOFsFromMasterNode
        #get global DOF numbers of active DOFs for this joint from master nodes
        aMap[i,1] = (masterNodeNum-1)*6 + activeDof[i]
    end

    #determine global active DOFs associated with slave node
    aMap2 =zeros(length(slaveActiveDof))
    for i=1:length(slaveActiveDof)
        aMap2[i] = (slaveNodeNum-1)*6 + slaveActiveDof[i]
    end

    if (!isempty(aMap2)) #create overall map of active DOFs associated with this joint
        aDOF = [aMap;aMap2]
    else
        aDOF = aMap
    end
    return Tda,dDOF,aDOF
end

function determineActiveDofsFromSlaveNode(slaveDof,numDofPerNode)
    #This function determines the local master DOF associated with a local slave DOF.
    # Get size
    count = 1;
    for i=1:numDofPerNode #loop over number of DOF per node
        if !in(i,slaveDof) #if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node
            #         slaveNodeActiveDof(count) = i;
            count = count + 1;

        end
    end

    if count>1
        count = 1
        slaveNodeActiveDof = zeros(count)
        for i=1:numDofPerNode #loop over number of DOF per node
            if !in(i,slaveDof) #if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node
                slaveNodeActiveDof[count] = i
                count = count + 1

            end
        end
    else
        slaveNodeActiveDof = []
    end

    return slaveNodeActiveDof
end

function extractdaInfo(joint,numNodes,numDofPerNode)
    #This function gets the total number of DOFs in the model, active
    #number of DOFs in the model, and a list of slave DOFs that will be
    #eliminated by joint constraints.

    adNumDof = numNodes*numDofPerNode; #total number of DOFs (active and dependent)

    numJoints=size(joint)[1] #get number of joints

    #Get Count
    dependentCount = 0
    count = 1
    for i=1:numJoints #loop over number of joints
        if (joint[i,4]==0 || joint[i,4]==5) #for a "fixed/weld" joint type or rigid bar constraint
            for j=1:6
                count = count + 1;
            end

        end

        if (joint[i,4]==1) #for a "pinned" joint type
            for j=1:3
                count = count + 1;
            end
        end

        if (joint[i,4]==2) #for a single axis hinge joint along a local "2" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

        if (joint[i,4]==3) #for a single axis hinge =joint along a local "1" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

        if (joint[i,4]==4) #for a single axis hinge joint along a local "3" axis of a joint
            for j=1:5
                count = count + 1;
            end
        end

    end

    slaveDof = zeros(count)
    count = 1
    for i=1:numJoints #loop over number of joints
        if (joint[i,4]==0 || joint[i,4]==5) #for a "fixed/weld" joint type or rigid bar constraint
            con = [1 2 3 4 5 6] #all DOFs of a slave node are constrained
            dependentCount = dependentCount + 6 #increment number of dependent DOFs
            for j=1:6
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end

        end

        if (joint[i,4]==1) #for a "pinned" joint type
            con = [1 2 3] #only translational (first 3) DOFs of a slave node are  constrained
            dependentCount = dependentCount + 3 #increment number of dependent DOFs
            for j=1:3
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==2) #for a single axis hinge joint along a local "2" axis of a joint
            if ((abs(abs(joint[i,7])-90))<1.0e-3 || (abs(abs(joint[i,7])-270))<1.0e-3)
                con=[1 2 3 5 6]
            else
                con=[1 2 3 4 6] #all but 5th DOF of a  slave node are constrained
            end
            dependentCount = dependentCount + 5 #increment number of dependent DOFs
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==3) #for a single axis hinge =joint along a local "1" axis of a joint
            if ((abs(abs(joint[i,8])-90))<1.0e-3 || (abs(abs(joint[i,8])-270))<1.0e-3)
                con = [1 2 3 4 5]
            elseif ((abs(abs(joint[i,7])-90))<1.0e-3 || (abs(abs(joint[i,7])-270))<1.0e-3)
                con = [1 2 3 4 6]
            else
                con = [1 2 3 5 6] #all but the 4th DOF of a slave node are constrained
            end
            dependentCount = dependentCount + 5 #increment number of dependent DOFs
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

        if (joint[i,4]==4) #for a single axis hinge joint along a local "3" axis of a joint
            if ((abs(abs(joint[i,8])-90))<1.0e-3 || (abs(abs(joint[i,8])-270))<1.0e-3)
                con = [1 2 3 5 6]
                if ((abs(abs(joint(i,7))-90))<1.0e-3 || (abs(abs(joint(i,7))-270))<1.0e-3)
                    con = [1 2 3 4 6]
                end
            else
                con = [1 2 3 4 5]
            end
            dependentCount = dependentCount + 5 #all but the 6th DOF of a slave node are constrained
            for j=1:5
                slaveDof[count] = numDofPerNode*(joint[i,3]-1) + con[j] #assign slave DOFs (joint(i,3) is the slave node number associated with this joint
                count = count + 1
            end
        end

    end

    aNumDof = adNumDof - dependentCount #calculate number of active DOFs in the model

    return adNumDof,aNumDof,slaveDof

end

function createJointTransform(joint,numNodes,numDofPerNode)
    #createJointTransform   Creates transformation matrix for joint constaints
    #   [jointTransform,reducedDOF] = createJointTransform(joint,numNodes,
    #                                                      numDofPerNode)
    #
    #   This function calculates the eigenvalues and vectors of a structural
    #   dynamic system.
    #
    #   input:
    #   joint         = object containing joint data
    #   numModes      = number of nodes in mesh
    #   numDofPerNode = number of degrees of freedom per node
    #
    #   output:
    #   jointTransform = joint transformation matrix
    #   reducedDOF     = map of original DOF numbering to reduced DOF numbering

    numJoints=size(joint)[1]  #get number of joints in model

    #extract number of active DOFs, number of dependent DOFs, slave DOF numbers
    adNumDof,aNumDof,slaveDof = extractdaInfo(joint,numNodes,numDofPerNode)

    #initialize joint transformation matrix
    jointTransform = zeros(adNumDof,aNumDof)

    #form reduced DOF vector which maps original DOF numbering to reduced DOF
    #numbering

    #Get Count
    count = 1
    for i=1:numNodes*numDofPerNode #loop over total number of DOFs in model
        if !in(i,slaveDof)
            #         reducedDOF(count) = i #if DOF is NOT a slave DOF include it in reducedDOF
            count = count + 1
        end
    end

    reducedDOF = zeros(Int,count-1)
    count = 1
    for i=1:numNodes*numDofPerNode #loop over total number of DOFs in model
        if !in(i,slaveDof)
            reducedDOF[count] = i #if DOF is NOT a slave DOF include it in reducedDOF
            count = count + 1
        end
    end

    #create identity portion of transformation matrix (This is done by Craig,
    #but here the original DOF ordering is retained
    for i=1:aNumDof #loop over number of active DOFs
        jointTransform[reducedDOF[i],i] = 1.0 #mapping of active DOFs in full DOF list to reduced DOF list
    end

    #impose Tda portion of identity matrix and map to appropriate locations

    for i=1:numJoints # loop of number of joints in the model
        jointType = joint[i,4] #get joint type
        slaveNodeNum = joint[i,3] #get slave node number associated with joint
        masterNodeNum = joint[i,2] #get master node number associated with joint
        psi = joint[i,7] #get psi orientation angle associated with joint
        theta = joint[i,8] #get theta orientation angle associated with joint

        #Tda is a local transform between dependent and active DOFs for nodes
        #associated with a particular joint, dDOF is a listing of dependent
        #global DOFs associated with this joint, aDOF is a listing of
        #active global DOFs associated with this joint.
        Tda,dDOF,aDOF =  createTda(jointType,slaveNodeNum,masterNodeNum,psi,theta,joint[i,:])

        for m=1:length(aDOF) #loop over global active DOFs associated with joint
            for k = 1:length(dDOF) #loop over global dependent DOFs associated with joint
                entry=findall(x->x==aDOF[m],reducedDOF)[1]  #determine reduced DOF associated with active DOF from original DOF listing
                jointTransform[dDOF[k],entry] = Tda[k,m]  #map local joint transformation matrix (Tda) to entries in global transformation matrix (jointTransform)
            end
        end
    end

    return jointTransform, reducedDOF
end


function calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained)
    #This function searches over all DOFs in a structural model and
    #determines and returns "dofVector" containing only unconstrained DOFs

    #loop over all DOFs in the model checking if constrained by BC or not
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                #             dofVector(index) = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    dofVector = zeros(Int,index)
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if (isConstrained[(i-1)*numDofPerNode + j]) == 0
                dofVector[index] = (i-1)*numDofPerNode + j #DOF vector only contains unconstrained DOFs
                index = index + 1
            end
        end
    end

    return dofVector
end

function constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,BC)
    #This function creates a map of unconstrained DOFs between a full
    #listing and reduced listing (aftger constraints have been applied)

    bcdoflist=zeros(Int, BC.numpBC)

    #create a listing of constrained DOFs from boundary condition file
    for i=1:BC.numpBC
        bcnodenum = BC.pBC[i,1]
        bcdofnum = BC.pBC[i,2]
        bcdoflist[i] = (bcnodenum-1)*numDofPerNode + bcdofnum
    end

    dofList = calculateReducedDOFVector(numNodes,numDofPerNode,BC.isConstrained) #calculate a reduced (unconstrained) DOF vector

    redVectorMap = zeros(numReducedDof)

    for i=1:numReducedDof

        if (i in bcdoflist)              #creates a map of unconstrained reduced DOFs
            redVectorMap[i] = -1.0
        else
            index = findall(x->x==i,dofList)[1]
            redVectorMap[i] = index
        end

    end
    return redVectorMap
end



function calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
    #calculateBCMap   calculates a boundary condition map
    #   [bcMap] = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
    #
    #   This function creates a boundary condition map between full and reduced
    #   dof listing as a result of constraints.
    #
    #      input:
    #      numpBC            = number of boundary conditions
    #      pBC               = array of boundary  condition data
    #      numDofPerNode     = number of degrees of freedom per node
    #      reducedDofList    = array of reduced DOF numbering
    #
    #      output:
    #      elStorage         = map for boundary conditions between full and
    #                          reduced dof list


    constrainedDof = zeros(numpBC)
    for i=1:numpBC
        constrainedDof[i] = (pBC[i,1]-1)*numDofPerNode + pBC[i,2]  #creates an array of constrained DOFs
    end
    constrainedDof = sort(constrainedDof)

    reducedDOFCount = length(reducedDofList)

    bcMap = zeros(reducedDOFCount)
    index = 1
    for i=1:reducedDOFCount
        if reducedDofList[i] in constrainedDof  #searches reduced DOF for constrained DOFs
            bcMap[i] = -1
        else
            bcMap[i] = index
            index = index + 1
        end
    end

    return bcMap

end


function calculateShapeFunctions(elementOrder,xi,x)
    #calculateShapeFunctions Calculates Lagrange shape functions
    #   [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi,x)
    #
    #   This function calculates the Lagrange shape function, shape
    #   function derivative, and Jacobian to map between the local element
    #   domain and physical length of the element. The shape function
    #   derivative is defined with respect to the physical length domain. The
    #   shape functions may be linear or quadratic in order.
    #
    #      input:
    #      elementOrder = order of element: 1 linear, 2 quadratic
    #      xi           = guass point values to evaluate shape functions at
    #      x            = nodal coordinates in physical length domain
    #
    #      output:
    #      N            = shape function value at specified gauss points
    #      p_N_x        = shape function derivative w.r.t physical length
    #                     domain at specified gauss points
    #      Jac          = Jacobian for mat between local element domain and
    #                     physical length domain.

    # N shape function
    # p_N_xi partial derivative of shape function w.r.t. xi

    #Linear interpolation functions
    N = zeros(elementOrder+1)
    p_N_xi = zeros(elementOrder+1)
    if elementOrder == 1
        N[1] = 0.5*(1.0 - xi)
        N[2] = 0.5*(1.0 + xi)

        p_N_xi[1] = -0.5
        p_N_xi[2] = 0.5
    end

    #Quadratic interpolation functions
    if elementOrder == 2
        N[1] = 0.5*(xi-1.0)*xi
        N[2] = 1.0-xi^2
        N[3] = 0.5*(xi+1.0)*xi

        p_N_xi[1] = xi - 0.5
        p_N_xi[2] = -2.0*xi
        p_N_xi[3] = xi + 0.5
    end

    numNodesPerEl = length(N)
    Jac=0.0
    for i=1:numNodesPerEl
        Jac = Jac + p_N_xi[i]*x[i]
    end

    p_N_x = zeros(numNodesPerEl)
    for i=1:numNodesPerEl
        p_N_x[i] = p_N_xi[i]/Jac
    end
    return N,p_N_x,Jac
end


function interpolateVal(valNode,N)
    valGP = 0.0
    for i=1:length(N)
        valGP = valGP + N[i]*valNode[i]
    end
    return valGP
end

#Element calculation functions---------------------------------

function calculateElement1(EA,integrationFactor,N1,N2,K)
    #This function is a general routine to calculate an element matrix
    len1 = length(N1)
    len2 = length(N2)
    for i=1:len1
        for j=1:len2
            K[i,j] = K[i,j] + EA*N1[i]*N2[j]*integrationFactor
        end
    end
    return K
end

# function [F] = calculateVec1(f,integrationFactor,N,F)
# #This function is a general routine to calculate an element vector
#     len=length(N)
#     for i=1:len
#         F(i) = F(i) + f*N(i)*integrationFactor
#     end
#
# end

function getGP(numGP)
    #getGP Defines Gauss point information for numerical integration
    #   [xi,weight] = getGP(numGP)
    #
    #   This function defines gauss point coordinates in a local element frame
    #   and the associated weights for Gaussian quadrature numerical
    #   integration.
    #
    #      input:
    #      numGP        = number of quad points used for integration
    #
    #      output:
    #      xi           = list of quad point coordinates in local element frame
    #      weight       = associated weights for quad point coordinate

    #define Gauss integration points
    xi = zeros(numGP)
    weight = zeros(numGP)
    if (numGP == 1)
        xi[1] = 0
        weight[1] = 2.0
    elseif (numGP == 2)
        xi[1] = -sqrt(1/3)
        xi[2] = sqrt(1/3)
        weight[1] = 1.0
        weight[2] = 1.0
    elseif (numGP == 3)
        xi[1] = -sqrt(3.0/5.0)
        xi[2] = 0.0
        xi[3] = sqrt(3.0/5.0)
        weight[1] = 5.0/9.0
        weight[2] = 8.0/9.0
        weight[3] = 5.0/9.0
    elseif (numGP == 4)
        xi[1] = sqrt((3.0-2*sqrt(6.0/5.0))/7.0)
        xi[2] = -sqrt((3.0-2*sqrt(6.0/5.0))/7.0)
        xi[3] = sqrt((3.0+2*sqrt(6.0/5.0))/7.0)
        xi[4] = -sqrt((3.0+2*sqrt(6.0/5.0))/7.0)

        weight[1] = (18+sqrt(30))/36.0
        weight[2] = (18+sqrt(30))/36.0
        weight[3] = (18-sqrt(30))/36.0
        weight[4] = (18-sqrt(30))/36.0
    end
    return xi, weight
end

function calculateElementMass(rhoA,rhoIyy,rhoIzz,rhoIyz,rhoJ,ycm,zcm,x,y,z,integrationFactor,M,Itens,xm)
    #This function calculates element mass properties.
    M = M + rhoA*integrationFactor
    y = y + ycm
    z = z + zcm

    #Total MOI = parallel axis theorem + local MOI
    Itens = Itens + rhoA.*integrationFactor.*[(y^2+z^2)  -x*y  -x*z
    -x*y  (x^2+z^2) -y*z
    -x*z -y*z (x^2+y^2)] + integrationFactor.*[rhoJ 0 0
                                                0 rhoIyy rhoIyz
                                                0 rhoIyz rhoIzz]

    xm[1] =  xm[1] + x*rhoA*integrationFactor
    xm[2] =  xm[2] + y*rhoA*integrationFactor
    xm[3] =  xm[3] + z*rhoA*integrationFactor

    return M,Itens,xm
end


function ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,nodalStiffnessTerms,nodalLoads)
    #ConcMassAssociatedWithElement gets concentrated terms associated w/ el
    #   [mass,stiff,load,modJoint,modNodalMassTerms,...
    #    modNodalStiffnessTerms,modNodalLoads] = ...
    #     ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,...
    #     nodalStiffnessTerms,nodalLoads)
    #
    #   This function compiles concentrated mass, stiffness, and load
    #   associated with a node from both ndl and joint files. The mod*
    #   variables are passed back with these terms removed to prevent
    #   duplicate application of shared nodal terms between elements
    #
    #      input:
    #      conn                = connectivity list for element
    #      joint               = joint array for nodal terms
    #      nodalMassTerms      = listing of concentrated nodal mass terms
    #      nodalStiffnessTerms = listing of concentrated nodal stiffness terms
    #      nodalLoads          = listing of concentrated nodal loads terms
    #
    #
    #      output:
    #      mass                = array of concentrated mass associated with element
    #      stiff               = array of concentrated stiffness associated with
    #                            element
    #      load                = array of concentrated loads associated with element
    #      modJoint            = modified joint object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modNodalMassTerms   = modified nodal mass object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modalStiffnessTerms = modified nodal stiffness object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modNodalLoads       = modified nodal loads object removing nodal terms that
    #                            have/will be applied to the element calculations

    node1 = conn[1] #define node #1 and node #2
    node2 = conn[2]

    mass1=0  #initialize concentrated mass amd moi for nodes
    mass2=0
    moix1=0
    moiy1=0
    moiz1=0
    moix2=0
    moiy2=0
    moiz2=0

    stiff1x=0 #initialize concentrated stifness for nodes
    stiff2x=0
    stiff1y=0
    stiff2y=0
    stiff1z=0
    stiff2z=0
    stiff1mx=0
    stiff2mx=0
    stiff1my=0
    stiff2my=0
    stiff1mz=0
    stiff2mz=0

    f1x = 0   #initialize concentrated loads/moments
    f2x = 0
    f1y = 0
    f2y = 0
    f1z = 0
    f2z = 0
    m1x =0
    m2x =0
    m1y =0
    m2y =0
    m1z =0
    m2z =0

    modJoint = joint                         #create copies of joint, and nodal mass, stiffness, loads arrays
    modNodalMassTerms = nodalMassTerms
    modNodalStiffnessTerms = nodalStiffnessTerms
    modNodalLoads = nodalLoads

    numJoints,_=size(joint)    #get number of joints in model

    if numJoints > 0
        node1flag=joint[:,2].==node1  #see if nodes are associated with a joint constraint as a master node
        node2flag=joint[:,2].==node2
    else
        node1flag = false
        node2flag = false
        mass1 = 0.0
        mass2 = 0.0
    end

    for i=1:numJoints           #if nodes are associated with joint constraint, use (if any) mass and stiffness specification from the joint file
        if node1flag[i]==1
            mass1 = mass1+joint[i,5]
            #             stiff1x = stiff1x + joint[i,6]
            #             stiff1y = stiff1y + joint[i,6]
            #             stiff1z = stiff1z + joint[i,6]
            #             stiff1mx = stiff1mx + joint[i,6]
            #             stiff1my = stiff1my + joint[i,6]
            #             stiff1mz = stiff1mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
        if node2flag[i]==1
            mass2 = mass2+joint[i,5]
            #             stiff2x = stiff2x + joint[i,6]
            #             stiff2y = stiff2y + joint[i,6]
            #             stiff2z = stiff2z + joint[i,6]
            #             stiff2mx = stiff2mx + joint[i,6]
            #             stiff2my = stiff2my + joint[i,6]
            #             stiff2mz = stiff2mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
    end

    #apply concentrated mass/stiffness from NDL file

    for i=1:length(nodalMassTerms)   #if node is specified in nodal mass terms file add to mass properties for this node
        node1flagM=nodalMassTerms[i].nodeNum.==node1
        node2flagM=nodalMassTerms[i].nodeNum.==node2
        if node1flagM==1
            if nodalMassTerms[i].dof == 1
                mass1 = mass1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 4
                moix1 = moix1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy1 = moiy1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz1 = moiz1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
        if node2flagM==1
            mass2 = mass2+nodalMassTerms[i].val
            modNodalMassTerms[i].val = 0.0

            if nodalMassTerms[i].dof == 4
                moix2 = moix2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy2 = moiy2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz2 = moiz2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
    end



    for i=1:length(nodalStiffnessTerms)     #if node is specified in nodal stiffness terms file add to stiffness properties for this node
        node1flagK=nodalStiffnessTerms[i].nodeNum.==node1
        node2flagK=nodalStiffnessTerms[i].nodeNum.==node2
        if node1flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff1x = stiff1x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff1y = stiff1y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff1z = stiff1z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff1mx = stiff1mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff1my = stiff1my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff1mz = stiff1mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
        if node2flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff2x = stiff2x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff2y = stiff2y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff2z = stiff2z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff2mx = stiff2mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff2my = stiff2my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff2mz = stiff2mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
    end

    for i=1:length(nodalLoads)  #if node is specified in nodal forces terms file add to concentrated force for this node
        node1flagF=nodalLoads[i].nodeNum.==node1
        node2flagF=nodalLoads[i].nodeNum.==node2
        if node1flagF==1
            if nodalLoads[i].dof==1
                f1x = f1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f1y = f1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f1z = f1z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m1x = m1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m1y = m1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m1z = m1z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
        if node2flagF==1
            if nodalLoads[i].dof==1
                f2x = f2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f2y = f2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f2z = f2z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m2x = m2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m2y = m2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m2z = m2z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
    end


    #compile nodal concentrated terms into mass, stiffness, and load arrays
    mass = [mass1 mass2
    moix1 moix2
    moiy1 moiy2
    moiz1 moiz2]

    stiff = [stiff1x stiff2x
    stiff1y stiff2y
    stiff1z stiff2z
    stiff1mx stiff2mx
    stiff1my stiff2my
    stiff1mz stiff2mz]

    load = [f1x f2x
    f1y f2y
    f1z f2z
    m1x m2x
    m1y m2y
    m1z m2z]

    return mass,stiff,load,modJoint,modNodalMassTerms,modNodalStiffnessTerms,modNodalLoads

end
