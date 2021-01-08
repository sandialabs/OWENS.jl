#TODO: write file output to match OWENS original files
#TODO: especially write output for NUMAD excel input
function create_mesh(;Ht = 15.0, #tower height before blades attach
    Hb = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nstrut = 2,
    strut_mout_ratio = 0.1, #distance from top/bottom
    ntelem = 20, #tower elements
    nbelem = 20, #blade elements
    nselem = 2,  #strut elements
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1)) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

    v_space_t = (Ht+Hb)/(ntelem)
    v_space_b = (Hb)/(nbelem)

    bshapex = R .* bshapex./maximum(bshapex)
    bshapez = Hb .* bshapez./maximum(bshapez)

    ####################################
    ##------------Tower--------------##
    ####################################
    function discretize_z(Ht,Hb,strut_mout_ratio,nelem; offset=0)
        mesh_z_inner = collect(LinRange(0,Ht+Hb,nelem+1))
        # Insert Bottom Blade Mount Point
        idx_bot_tower_bld_inner = 1
        if offset ==0 #but only if it isn't being used for the blade
            # insert!(mesh_z_inner,idx_bot_tower_bld_inner,Ht)
            mesh_z_inner = sort([mesh_z_inner;Ht])
            idx_bot_tower_bld_inner = findall(x -> x == Ht,mesh_z_inner)[1]
        end

        # Insert Bottom Strut Mount Point
        strut_z = Ht+strut_mout_ratio*Hb
        mesh_z_inner = sort([mesh_z_inner;strut_z])
        idx_bot_tower_strut_inner = findall(x -> x == strut_z,mesh_z_inner)[1]

        # Insert Top Strut Mount Point
        strut_z = Ht+(1-strut_mout_ratio)*Hb
        mesh_z_inner = sort([mesh_z_inner;strut_z])
        idx_top_tower_strut_inner = findall(x -> x == strut_z,mesh_z_inner)[1]

        # Don't Insert Top Blade Mount Point since it is already there at the exact spot
        idx_top_tower_blade_inner = length(mesh_z_inner)
        return mesh_z_inner, idx_bot_tower_bld_inner+offset, idx_bot_tower_strut_inner+offset, idx_top_tower_strut_inner+offset, idx_top_tower_blade_inner+offset
    end
    mesh_z, idx_bot_tower_bld, idx_bot_tower_strut, idx_top_tower_strut, idx_top_tower_blade = discretize_z(Ht,Hb,strut_mout_ratio,ntelem)
    # Create the x and y components
    mesh_x = zero(mesh_z)
    mesh_y = zero(mesh_z)

    # intra-tower connectivity
    conn = zeros(length(mesh_z)-1,2)
    conn[:,1] = collect(1:length(mesh_z)-1)
    conn[:,2] = collect(2:length(mesh_z))

    ####################################
    ##------------Blades--------------##
    ####################################

    #connection points on tower already exist, as do the strut connection points if we reuse the tower z discretization as the blade discretization
    bld_Z, idx_bot_lbld_tower, idx_bot_lbld_strut, idx_top_lbld_strut, idx_top_lbld_tower = discretize_z(0.0,Hb,strut_mout_ratio,nbelem; offset = length(mesh_z))

    idx_bot_lbld_strut += 1
    idx_top_lbld_strut += 1

    idx_bot_rbld_tower = idx_bot_lbld_tower + nbelem+1+nstrut
    idx_bot_rbld_strut = idx_bot_lbld_strut + nbelem+1+nstrut
    idx_top_rbld_strut = idx_top_lbld_strut + nbelem+1+nstrut
    idx_top_rbld_tower = idx_top_lbld_tower + nbelem+1+nstrut

    if bshapex == zeros(nbelem+1)
        bld_X = R.*(1.0.-4.0.*(bld_Z/Hb.-.5).^2)
    else
        # Ensure the blade shape conforms to the turbine height and radius specs
        bld_X = FLOWMath.akima(bshapez,bshapex,bld_Z)
    end
    bld_Y = zero(bld_X)
    bld_Z .+= Ht

    #TODO: arbitrary number of blades rotated about center
    mesh_z = [mesh_z;bld_Z;bld_Z]
    mesh_x = [mesh_x;-bld_X;bld_X]
    mesh_y = [mesh_y;bld_Y;bld_Y]


    # intra-blade connectivity and blade-tower connections
    # Left Blade
    conn_b = zeros(length(bld_Z)-1,2)
    conn_b[:,1] = collect(idx_bot_lbld_tower:1:idx_top_lbld_tower-1)
    conn_b[:,2] = collect(idx_bot_lbld_tower+1:1:idx_top_lbld_tower)
    conn = [conn;conn_b]

    # Right Blade
    conn_b[:,1] = collect(idx_bot_rbld_tower:1:idx_top_rbld_tower-1)
    conn_b[:,2] = collect(idx_bot_rbld_tower+1:1:idx_top_rbld_tower)
    conn = [conn;conn_b]


    ####################################
    ##------------Struts--------------##
    ####################################

    bot_strut_z = strut_mout_ratio*Hb+Ht
    top_strut_z = (1-strut_mout_ratio)*Hb+Ht
    if bshapex == zeros(nbelem+1)
        bot_strut_x = R.*(1.0.-4.0.*((bot_strut_z-Ht)/Hb.-.5).^2)
        top_strut_x = R.*(1.0.-4.0.*((top_strut_z-Ht)/Hb.-.5).^2)
    else
        bot_strut_x = FLOWMath.akima(bshapez,bshapex,(bot_strut_z-Ht))
        top_strut_x = FLOWMath.akima(bshapez,bshapex,(top_strut_z-Ht))
    end

    bot_strut_X = collect(LinRange(0.0,bot_strut_x,nselem+1))
    top_strut_X = collect(LinRange(0.0,top_strut_x,nselem+1))

    idx_bot_lstrut_tower = length(mesh_z)+1
    idx_bot_lstrut_blade = idx_bot_lstrut_tower + 2
    idx_bot_rstrut_tower = idx_bot_lstrut_blade + 1
    idx_bot_rstrut_blade = idx_bot_rstrut_tower + 2
    idx_top_lstrut_tower = idx_bot_rstrut_blade + 1
    idx_top_lstrut_blade = idx_top_lstrut_tower + 2
    idx_top_rstrut_tower = idx_top_lstrut_blade + 1
    idx_top_rstrut_blade = idx_top_rstrut_tower + 2
    mesh_z = [mesh_z;ones(nselem+1).*bot_strut_z ; ones(nselem+1).*bot_strut_z ; ones(nselem+1).*top_strut_z ; ones(nselem+1).*top_strut_z]
    mesh_x = [mesh_x;-bot_strut_X;bot_strut_X;-top_strut_X;top_strut_X]
    mesh_y = [mesh_y;zeros((nselem+1)*4)]

    conn_s = zeros(8,2)

    conn_s[:,1] = [idx_bot_lstrut_tower,idx_bot_lstrut_tower+1,
    idx_bot_rstrut_tower,idx_bot_rstrut_tower+1,
    idx_top_lstrut_tower,idx_top_lstrut_tower+1,
    idx_top_rstrut_tower,idx_top_rstrut_tower+1]

    conn_s[:,2] = [idx_bot_lstrut_tower+1,idx_bot_lstrut_blade,
    idx_bot_rstrut_tower+1,idx_bot_rstrut_blade,
    idx_top_lstrut_tower+1,idx_top_lstrut_blade,
    idx_top_rstrut_tower+1,idx_top_rstrut_blade]

    conn = [conn;conn_s]
    #space out the mesh numerically to avoid numerical issues
    for i = 1:length(mesh_z)-1
        if isapprox(mesh_z[i],mesh_z[i+1];atol = 1e-6)
            mesh_z[i+1] -= 1e-4
        end
    end



    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:idx_top_tower_blade] .= 1 #Tower
    # meshtype[idx_bot_lbld_tower:idx_top_rbld_tower] .= 0 #Blades
    meshtype[idx_bot_lstrut_tower:end] .= 2 #Struts

    #TODO: this is hard coded for two blades, make arbitrary
    meshSeg = zeros(1+2+nstrut*2) #tower, two blades, and two struts that support the two blades

    meshSeg[1] = ntelem+nstrut+1
    meshSeg[2:3] .= nbelem+nstrut
    meshSeg[4:end] .= nselem

    mymesh = Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,conn,meshtype,meshSeg)

    ####################################
    ##----------Joint Matrix----------##
    ####################################

    #Connect Tower Bottom to Blades
    jointconn = [idx_bot_tower_bld idx_bot_lbld_tower
    idx_bot_tower_bld idx_bot_rbld_tower
    #Connect Tower Bottom to Struts
    idx_bot_tower_strut idx_bot_lstrut_tower
    idx_bot_tower_strut idx_bot_rstrut_tower
    #Connect Tower Top to Struts
    idx_top_tower_strut idx_top_lstrut_tower
    idx_top_tower_strut idx_top_rstrut_tower
    #Connect Tower Top to Blades
    idx_top_tower_blade idx_top_lbld_tower
    idx_top_tower_blade idx_top_rbld_tower
    #Connect Bottom Left Blade to Strut
    idx_bot_lbld_strut idx_bot_lstrut_blade
    #Connect Top Left Blade to Strut
    idx_top_lbld_strut idx_top_lstrut_blade
    #Connect Bottom Right Blade to Strut
    idx_bot_rbld_strut idx_bot_rstrut_blade
    #Connect Top Right Blade to Strut
    idx_top_rbld_strut idx_top_rstrut_blade]

    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh)
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,2],ort.elNum) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2]-1,ort.elNum) #TODO: we get away with this since the elements are increasing and there are no two point objects (and this is only a problem for the top of the tower connecting to the blade tops)
        end
        if length(elnum_of_joint)==0
            elnum_of_joint = findall(x->x==jointconn[jnt,1]-1,ort.elNum)
        end
        Psi_d_joint[jnt] = ort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = ort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

    #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint) zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    return mymesh, ort, myjoint
end

function calculateElementOrientation(mesh)
    #calculateElementOrientation calculates the orienation of elements in mesh
    #   [elOr] = calculateElementOrientation(mesh)
    #
    #   This function calculates the orientation of elements in a mesh.
    #
    #      input:
    #      mesh       = object containing mesh data
    #
    #      output:
    #      elOr       = object containing element orientation data

    numEl = mesh.numEl #get number of elements
    Psi_d=zeros(numEl) #initialize Psi, Theta, Twist, and Offset Arrays
    Theta_d=zeros(numEl)
    twist_d=zeros(numEl)
    Offset=zeros(3,numEl)    #offset is the hub frame coordinate of node 1 of the element
    elNum=zeros(numEl) #initialize element number array


    #calculate "mesh centroid"
    meshCentroid = [Statistics.mean(mesh.x) Statistics.mean(mesh.y) Statistics.mean(mesh.z)] #calculate a geometric centroid using all nodal coordinates
    lenv = zeros(numEl)
    for i = 1:numEl #loop over elements
        n1 = Int(mesh.conn[i,1]) #n1 := node number for node 1 of element i
        n2 = Int(mesh.conn[i,2]) #n2 := node number for node 2 of element i

        p1 = [mesh.x[n1] mesh.y[n1] mesh.z[n1]] #nodal coordinates of n1
        p2 = [mesh.x[n2] mesh.y[n2] mesh.z[n2]] #nodal coordinates of n2
        Offset[:,i] = p1 #set offset as position of n1

        v=p2-p1 #define vector from p1 to p2
        lenv[i] = LinearAlgebra.norm(v) #calculate element lengtt

        Psi_d[i],Theta_d[i] = calculatePsiTheta(v) #calculate elment Psi and Theta angles for orientation
        elNum[i] = mesh.conn[i,1] #get elemetn number

        nVec1,nVec2,nVec3 = rigidBodyRotation(0,0,1,[Psi_d[i],Theta_d[i]],[3,2]) #tranform a local normal "flapwise" vector in the element frame to the hub frame
        nVec = [nVec1 nVec2 nVec3]

        #for consistency, force the "flapwise" normal vector of an element to be
        #away from the machine

        # Mesh Type: 0-blade 1-tower 2-strut
        if mesh.type[i]==2
            refVector = [0;0;1]
        elseif mesh.type[i]==1
            refVector = [1;0;0]
        else
            refVector = p1-meshCentroid
        end

        refVector = refVector./LinearAlgebra.norm(refVector)
        dotTest = LinearAlgebra.dot(nVec,refVector)

        if dotTest<0 && abs(dotTest)>1.0e-4 #if vectors are not more or less aligned the normal vector is pointed inwards, towards the turbine (twist_d by 180 degrees)
            twist_d[i] = 180.0
        elseif abs(dotTest)<1.0e-4
            twist_dtemp = 90.0
            nVec1,nVec2,nVec3 = rigidBodyRotation(0,0,1,[Psi_d[i],Theta_d[i],twist_dtemp],[3,2,1])
            nVec = [nVec1 nVec2 nVec3]
            dotTest2 = LinearAlgebra.dot(nVec,refVector)
            if dotTest2<0 && abs(dotTest2)>1.0e-4 #if vectors are not more or less aligned the normal vector is pointed inwards, towards the turbine (twist_d by 180 degrees)
                twist_d[i] = twist_dtemp+180.0
                if abs(abs(twist_d[i])-270.0) < 1.0e-3
                    twist_d[i] = -90.0
                end
            else
                twist_d[i] = twist_dtemp
            end
        else  #the normal vector is pointed outwards, away from the turbine (no twist_d necessary)
            twist_d[i] = 0.0
        end

    end

    #assign data to element orientation (Ort) object
    return Ort(Psi_d,Theta_d,twist_d,lenv,elNum,Offset)
end

function createGeneralTransformationMatrix(angleArray,axisArray)
    #createGeneralTransformationMatrix  calculates general transformation matrix
    #   [dcmTotal] = createGeneralTransformationMatrix(angleArray,axisArray)
    #
    #   This function calculates the transformation matrix assocaited with a
    #   general Euler rotation sequence.
    #
    #      input:
    #      angleArray      = array of angles for Euler rotation sequence
    #      axisArray       = array of axis of rotatoins for Euler rotation
    #                        sequences
    #
    #      output:
    #      dcmTotal        = transformation matrix of specified euler rotation
    #                        sequence

    numRotations = length(angleArray) #get number of rotation to perform
    dcmArray = zeros(3,3,numRotations) #initialize individual rotation direction cosine matrix arrays

    for i=1:numRotations #calculate individual single rotatio direction cosine matrices
        dcmArray[:,:,i] = createSingleRotationDCM(angleArray[i],axisArray[i])
    end

    dcmTotal = dcmArray[:,:,1] #initialize dcmTotal as first rotation

    #multiply consecutive rotation sequency direction cosine matrices to arrive at overall transformation matrix
    for i=2:1:numRotations
        dcmTotal = dcmArray[:,:,i]*dcmTotal
    end

    return dcmTotal

end

function createSingleRotationDCM(angleDeg,axisNum)
    #This function creates a direction cosine matrix (dcm) associated
    #with a rotation of angleDeg about axisNum.

    angleRad = angleDeg*pi/180.0 #convert angle to radians

    if axisNum == 1 #direction cosine matrix about 1 axis
        dcm = [1.0 0.0 0.0
        0.0 cos(angleRad) sin(angleRad)
        0.0 -sin(angleRad) cos(angleRad)]
    elseif axisNum == 2 #direction cosine matrix about 2 axis
        dcm = [cos(angleRad) 0.0 -sin(angleRad)
        0.0 1.0 0.0
        sin(angleRad) 0.0 cos(angleRad)]
    elseif axisNum == 3 #direction cosine matrix about 3 axis
        dcm = [cos(angleRad) sin(angleRad) 0.0
        -sin(angleRad) cos(angleRad) 0.0
        0.0 0.0 1.0]
    else  #error catch
        error("Error: createSingleRotationDCM. Axis number must be 1, 2, or 3.")
    end

    return dcm

end

function rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)
    #rigidBodyRotation rotates a vector through a rotation sequence
    #   [H1,H2,H3] = rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)
    #
    #   This function performs a coordinate transformation from a local
    #   body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
    #   rotation sequence
    #
    #      input:
    #      B1        = array containing body frame 1 coordinates of points to be
    #                  mapped to the hub frame
    #      B2        = array containing body frame 2 coordinates of points to be
    #                  mapped to the hub frame
    #      B3        = array containing body frame 3 coordinates of points to be
    #                  mapped to the hub frame
    #     AngleArray = Array of angles for Euler rotation sequence
    #     AxisArray  = Array of axes for Euler rotation sequence
    #
    #      output:
    #      H1        = array containg hub frame 1 coordinates of points mapped
    #                  to the hub frame from body frame
    #      H2        = array containg hub frame 2 coordinates of points mapped
    #                  to the hub frame from body frame
    #      H3        = array containg hub frame 3 coordinates of points mapped
    #                  to the hub frame from body frame

    #This function performs a coordinate transformation from a local
    #body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
    #rotation sequence

    #That is CHtoB = [M3(SweepAngle)][M2(Theta)][M3(Psi)];

    #calculate coordinate transformation matrix from element frame to
    #hub frame (CBtoH)
    dcm = createGeneralTransformationMatrix(AngleArray,AxisArray)
    C = dcm'

    #transform body coordinatized vector to be coordinatized in the hub
    #frame
    H1 = C[1,1].*B1 + C[1,2].* B2 + C[1,3].*B3
    H2 = C[2,1].*B1 + C[2,2].* B2 + C[2,3].*B3
    H3 = C[3,1].*B1 + C[3,2].* B2 + C[3,3].*B3

    return H1,H2,H3
end

function calculatePsiTheta(v)
    #calculatePsiTheta calculates the orienation of a single element
    #   [Psi,Theta] = calculatePsiTheta(v)
    #
    #   This function calculates the orientation of a single element. A local
    #   element frame is related to a hub frame through a transformation matrix
    #   CHtoE (transforming a vector from an element frame E to a global frame
    #   H) such that CHtoE = [M2(Theta)]*[M3(Psi)]. Here [M2( )] is a direction
    #   cosine matrix about a 2 axis and [M3( )] is a direction cosine matrix
    #   about a 3 axis.
    #
    #      input:
    #      v          = vector from node 1 to node 2 of an element
    #
    #      output:
    #      Psi        = "3" angle for element orientation (deg)
    #      Theta      = "2" angle for element orientation (deg)
    #                   see above for definition

    v = v./LinearAlgebra.norm(v) #normalize vector by its length
    Psi_d = atan(v[2],v[1])*180.0/pi #calculate sweep angle, convert to deg
    Theta_d = -asin(v[3])*180.0/pi #calculate theta angle, convert to deg

    return Psi_d,Theta_d
end

function makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)
    #readBDdata  reads boundary condition file
    #   [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
    #
    #   This function reads the boundray condition file and stores data in the
    #   boundary condition object.
    #
    #      input:
    #      bcfilename    = string containing boundary condition filename
    #      numNodes      = number of nodes in structural model
    #      numDofPerNode = number of degrees of freedom per node

    #      output:
    #      BC            = object containing boundary condition data

    totalNumDof = numNodes*numDofPerNode

    numsBC = 0
    nummBC = 0

    #create a vector denoting constrained DOFs in the model (0 unconstrained, 1
    #constrained)

    #calculate constrained dof vector
    isConstrained = zeros(totalNumDof,1)
    constDof = (pBC[:,1].-1)*numDofPerNode + pBC[:,2]
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if ((i-1)*numDofPerNode + j in constDof)
                isConstrained[index] = 1
            end
            index = index + 1
        end
    end
    numpBC = length(pBC[:,1])

    map = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDOFList)
    numReducedDof = length(jointTransform[1,:])
    redVectorMap = constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,numpBC,pBC,isConstrained) #create a map between reduced and full DOF lists

    BC = BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    map,
    redVectorMap)

    return BC

end
