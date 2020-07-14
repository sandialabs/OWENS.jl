include("calculateLambda.jl")

mutable struct JntTransform
    jointTransform
    reducedDOF
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

    elseif(jointType == 1) #for pinned joint type
        activeDof1 = [1 2 3 4 5 6] #active DOF list at joint
        slaveDof1 = [1 2 3] #slave DOF list at joint

        #determine local active DOFs associated with slave node
        slaveActiveDof1 = determineActiveDofsFromSlaveNode(slaveDof1,6)

        Rda1 = -[1.0*LinearAlgebra.I(3), zeros(3,6)] #from constraint equation for pinned joint
        Rdd1 = 1.0*LinearAlgebra.I(3)

        Tda,dDOF,aDOF = getNodeMaps(Rdd1,Rda1,masterNodeNum,slaveNodeNum,slaveDof1,activeDof1,slaveActiveDof1)

    elseif (jointType == 2)     #hinge axis along localy "2" frame of joint

        if((abs(abs(psi)-90))<1.0e-3 || (abs(abs(psi)-270))<1.0e-3)
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

    if(!isempty(aMap2)) #create overall map of active DOFs associated with this joint
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
                if((abs(abs(joint(i,7))-90))<1.0e-3 || (abs(abs(joint(i,7))-270))<1.0e-3)
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
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
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

    return JntTransform(jointTransform,reducedDOF)
end
