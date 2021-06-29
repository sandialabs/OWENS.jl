#TODO: write file output to match OWENS original files
#TODO: especially write output for NUMAD excel input
function create_mesh(;Ht = 15.0, #tower height before blades attach
    Hb = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nstrut = 2,
    nblade = 2, #TODO: this isn't really used downstream, it's hard coded for two, need to change
    strut_mout_ratio = 0.1, #distance from top/bottom
    ntelem = 20, #tower elements
    nbelem = 20, #blade elements
    nselem = 2,  #strut elements
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1)) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

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
    ####################################
    ##------------Tower--------------##
    ####################################
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
        bshapex = R .* bshapex./maximum(bshapex)
        bshapez = Hb .* bshapez./maximum(bshapez)
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

    #########################
    # .bld equivalent
    #########################

    #TODO: this is hard coded for two blades, make arbitrary
    meshSeg = zeros(1+2+nstrut*2) #tower, two blades, and two struts that support the two blades

    meshSeg[1] = ntelem+nstrut+1 # the +1 is the blade mount point
    meshSeg[2:3] .= nbelem+nstrut
    meshSeg[4:end] .= nselem

    structuralSpanLocNorm = zeros(nblade,nbelem+nstrut+1)
    structuralNodeNumbers = zeros(nblade,nbelem+nstrut+1)
    structuralElNumbers = zeros(nblade,nbelem+nstrut+1)

    # for iblade = 1:nblade
    span_len = sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-Ht).^2.0)
    structuralSpanLocNorm[1,:] = span_len./maximum(span_len)
    structuralSpanLocNorm[2,:] = structuralSpanLocNorm[1,:]

    structuralNodeNumbers[1,:] = collect(idx_bot_lbld_tower:idx_top_lbld_tower)
    structuralNodeNumbers[2,:] = collect(idx_bot_rbld_tower:idx_top_rbld_tower)

    structuralElNumbers[1,:] = structuralNodeNumbers[1,:].-1
    structuralElNumbers[1,end] = -1 #TODO: figure out why this is in the original OWENS setup and if it is used
    structuralElNumbers[2,:] = structuralNodeNumbers[2,:].-2
    structuralElNumbers[1,end] = -1
    # end

    mymesh = Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,conn,meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

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

function create_arcus_mesh(;
    Ht = 15.0, #tower height before blades attach
    Hb = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nblade = 3,
    ntelem = 4, #tower elements
    nbelem = 30, #blade elements
    ncelem = 10,  #cable elements
    c_mount_ratio = 0.05, #fraction of blade where the cable mounts
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1)) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

    ##################################
    #             _
    #           /| |\
    #          | | | )
    #           \|_|/
    #             |
    # Wires as tower tension support, for now are mounted to the blades like a bow
    ####################################
    ###------------Tower--------------##
    ####################################
    mesh_z = collect(LinRange(0,Ht,ntelem+1))
    # Create the x and y components
    mesh_x = zero(mesh_z)
    mesh_y = zero(mesh_z)

    t_topidx = length(mesh_z)

    # intra-tower connectivity
    conn = zeros(length(mesh_z)-1,2)
    conn[:,1] = collect(1:length(mesh_z)-1)
    conn[:,2] = collect(2:length(mesh_z))

    #####################################
    ###------------Blades--------------##
    #####################################

    #connection points on tower are simply the top of the tower connecting to the bottom of the blades
    bld_Z = collect(LinRange(0.0,Hb,nbelem+1))

    # The cable connections must be added though
    if maximum(Hb*c_mount_ratio .== bld_Z) # if we are at exactly an existing node, then offset our mount point
        c_mount_ratio += 1e-6
    end
    bld_Z = sort([bld_Z;Hb*c_mount_ratio])

    if bshapex == zeros(nbelem+1)
        bshapex = R .* bshapex./maximum(bshapex)
        bshapez = Hb .* bshapez./maximum(bshapez)
        bld_Y = R.*(1.0.-4.0.*(bld_Z/Hb.-.5).^2)
    else
        # Ensure the blade shape conforms to the turbine height and radius specs
        bld_Y = FLOWMath.akima(bshapez,bshapex,bld_Z)
    end
    bld_X = zero(bld_Y)

    bld_Z .+= Ht

    b_Z = []
    b_X = []
    b_Y = []
    # Now using standard VAWT convention, blade 1 is zero degrees at top dead center, or North/Y+
    # and they are offset counter clockwise
    b_topidx = zeros(Int,nblade)
    b_botidx = zeros(Int,nblade) .+ length(mesh_z)
    conn_b = zeros(length(bld_Z)-1,2)
    for ibld = 1:nblade
        myangle = (ibld-1)*2.0*pi/nblade
        b_Z = [b_Z;bld_Z]
        b_X = [b_X;bld_Y.*sin(myangle).+bld_X.*cos(myangle)]
        b_Y = [b_Y;bld_X.*sin(myangle).+bld_Y.*cos(myangle)]

        # Element joint indices
        b_botidx[ibld] = length(mesh_z)+1 + length(bld_Z)*(ibld-1)
        b_topidx[ibld] = length(mesh_z)+1 + length(bld_Z)*ibld-1

        # Intraconnectivity
        conn_b[:,1] = collect(b_botidx[ibld]:1:b_topidx[ibld]-1)
        conn_b[:,2] = collect(b_botidx[ibld]+1:1:b_topidx[ibld])
        conn = [conn;conn_b]
    end

    # Add to the mesh
    mesh_z = [mesh_z;b_Z]
    mesh_x = [mesh_x;b_X]
    mesh_y = [mesh_y;b_Y]

    # pick out the blade to cable mounting indices
    b2c_botidx = findall(x->x==Hb*c_mount_ratio+Ht,mesh_z)[1:nblade]


    #####################################
    ###------------Cables--------------##
    #####################################

    #Connect from the bottom to the top
    # For each blade, find the mounting location and draw a line to the top center confluence of blades
    c2b_botidx = zeros(Int,nblade)
    c2b_topidx = zeros(Int,nblade)
    conn_c = zeros(ncelem,2)
    for ibld = 1:nblade
        cstartidx = b2c_botidx[ibld]
        cxstart = mesh_x[cstartidx]
        cystart = mesh_y[cstartidx]
        czstart = mesh_z[cstartidx]

        # x y z end are at the top of the tower
        cxend = 0.0
        cyend = 0.0
        czend = Ht + Hb

        # Now draw the lines
        c_x = collect(LinRange(cxstart,cxend,ncelem+1))
        c_y = collect(LinRange(cystart,cyend,ncelem+1))
        c_z = collect(LinRange(czstart,czend,ncelem+1))

        # joint connections
        c2b_botidx[ibld] = length(mesh_z)+1
        c2b_topidx[ibld] = c2b_botidx[ibld]+length(c_z)-1

        # and add to the mesh
        mesh_x = [mesh_x;c_x]
        mesh_y = [mesh_y;c_y]
        mesh_z = [mesh_z;c_z]

        # Intraconnectivity
        conn_c[:,1] = collect(c2b_botidx[ibld]:1:c2b_topidx[ibld]-1)
        conn_c[:,2] = collect(c2b_botidx[ibld]+1:1:c2b_topidx[ibld])
        conn = [conn;conn_c]

    end

    #######################################
    ###   Cleanup/Derived parameters    ###
    #######################################

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
    meshtype[1:t_topidx] .= 1 #Tower
    # meshtype[idx_bot_lbld_tower:idx_top_rbld_tower] .= 0 #Blades
    meshtype[c2b_botidx[1]-1:end] .= 2 #Struts

    #########################
    # .bld equivalent
    #########################

    ncable = nblade

    # For a single blade
    meshSeg = zeros(1+nblade+ncable) #tower, blades, and cables

    meshSeg[1] = ntelem
    meshSeg[2:nblade+1] .= nbelem+1 # +1 for the cable mount
    meshSeg[nblade+2:end] .= ncelem

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z)) # +1 for the cable mount
    structuralNodeNumbers = zeros(nblade,length(bld_Z))
    structuralElNumbers = zeros(nblade,length(bld_Z))

    for iblade = 1:nblade

        # Normalized Span
        span_len = sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-Ht).^2.0)
        structuralSpanLocNorm[iblade,:] = span_len./maximum(span_len)

        # Node Numbers
        structuralNodeNumbers[iblade,:] = collect(b_botidx[iblade]:b_topidx[iblade])

        # Element Numbers
        structuralElNumbers[iblade,:] = structuralNodeNumbers[iblade,:].-iblade
        structuralElNumbers[iblade,end] = -1 #TODO: figure out why this is in the original OWENS setup and if it is used
    end

    mymesh = Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,conn,meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ######################################
    ####----------Joint Matrix----------##
    ######################################

    # Connect Tower Top to Blades bottom, then each cable to each blade bottom
    # Then each cable to each blade top, then the latter two blade tops to the first

    jointconn = zeros(Int,nblade*2-1+ncable*2,2)
    for ibld = 1:nblade
        # connect tower to blades
        jointconn[ibld,:] = [t_topidx b_botidx[ibld]]

        # connect cables to blades bases
        jointconn[ibld+nblade,:] = [b2c_botidx[ibld] c2b_botidx[ibld]]

        # connect cables to blades tops
        jointconn[ibld+nblade*2,:] = [b_topidx[ibld] c2b_topidx[ibld]]

        #Connect first Blade to all other blades
        if ibld>1
            jointconn[ibld-1+nblade*3,:] = [b_topidx[1] b_topidx[ibld]]
        end
    end

    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh)
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,2],ort.elNum) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2]-1,ort.elNum) #TODO: we get away with this since the elements are increasing and there are no two point objects
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

    map = GyricFEA.calculateBCMap(numpBC,pBC,numDofPerNode,reducedDOFList)
    numReducedDof = length(jointTransform[1,:])
    redVectorMap = GyricFEA.constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,numpBC,pBC,isConstrained) #create a map between reduced and full DOF lists

    BC = BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    map,
    redVectorMap)

    return BC

end



# NuMad_geom_xlscsv_file = "$module_path/../test/data/NuMad_Geom_SNL_5MW_D_Carbon_LCDT.csv"

function getPreCompOutput(numadIn;plyprops = plyproperties())
    # numadIn = readNuMadGeomCSV(NuMad_geom_xlscsv_file)

    n_stations = length(numadIn.span)


    mat = Array{Array{Composites.Material,1}}(undef,n_stations)
    lam = Array{Composites.Laminate,1}(undef,n_stations)
    precompinput = Array{PreComp.Input,1}(undef,n_stations)
    precompoutput = Array{PreComp.Output,1}(undef,n_stations)

    normalchord = numadIn.chord #TODO: if any sweep occurs the normal chord will need to be calculated since it isn't the regular chord anymore

    sloc = numadIn.span
    twist_d = zeros(length(numadIn.twist_d))
    for ii = 1:length(numadIn.twist_d) #TODO: real interpolation in the reading file
        twist_d[ii] = numadIn.twist_d[ii]
    end

    twistrate_d = PreComp.tw_rate(n_stations,sloc[1:n_stations],twist_d)
    leloc = numadIn.xoffset

    for i_station = 1:n_stations

        ######################################
        ######## Airfoil Shape Input #########
        ######################################

        #TODO: use actual airfoil at each section, ensure square matrix via interpolation so it can be stored as a 2D matrix
        af_xy = DelimitedFiles.readdlm("$(module_path)/../test/airfoils/$(numadIn.airfoil[i_station]).csv",',',Float64,skipstart = 0)

        # Normalize the surface points and make sure that they start at the trailing edge and loop around starting on the bottom side #TODO: add a check for both of these
        #TODO: simplify this since there is circshift happening on these points later on
        if af_xy[2,2]>af_xy[end-1,2]
            xaf1 = reverse(af_xy[:,1])./maximum(af_xy[:,1])
            yaf1 = reverse(af_xy[:,2])./maximum(af_xy[:,1])
        else
            xaf1 = af_xy[:,1]./maximum(af_xy[:,1])
            yaf1 = af_xy[:,2]./maximum(af_xy[:,1])
        end

        xaf = zeros(n_stations,length(xaf1))
        yaf = zeros(n_stations,length(yaf1))
        for i = 1:n_stations
            xaf[i,:] = xaf1[:] #break links since using in multiple areas
            yaf[i,:] = yaf1[:]
        end
        # find leading edge
        lei = argmin(abs.(xaf[i_station,:]))
        # shift so leading edge is first
        xpc = circshift(xaf[i_station,:],-(lei-1))
        ypc = circshift(yaf[i_station,:],-(lei-1))

        # assemble input
        # precompinput[i_station],mat[i_station],lam[i_station] = layup(normalchord[i_station],twist_d[i_station],twistrate_d[i_station],xpc,ypc,lam_t[i_station,:],usedmaterials,webloc,plyprops,leloc[i_station],orientation)

        ##################################
        ######## Materials Input #########
        ##################################
        mat_single = []
        usedmaterials = numadIn.stack_mat_types #["highmodulus_uni","highmodulus_weave","SNL_foam","highmodulus_weave","SNL_foam","highmodulus_weave","SNL_foam","highmodulus_weave"] #TODO: hook this up to the numad materials
        e1 = zeros(length(usedmaterials))
        e2 = zeros(length(usedmaterials))
        g12 = zeros(length(usedmaterials))
        anu12 = zeros(length(usedmaterials))
        density = zeros(length(usedmaterials))
        ply_thickness = zeros(length(usedmaterials))

        for i_mat = 1:length(usedmaterials)
            matnames = plyprops.names
            idx = usedmaterials[i_mat] #findall(matnames -> matnames == usedmaterials[i_mat],matnames) #TODO: determine best way to modify workflow to use names instead of a blind index that potentially might not align with the materials
            material = plyprops.plies[idx] #[idx[1]]
            push!(mat_single,material)
            e1[i_mat] = material.e1
            e2[i_mat] = material.e2
            g12[i_mat] = material.g12
            anu12[i_mat] = material.nu12
            density[i_mat] = material.rho
            ply_thickness[i_mat] = material.t
        end

        mat[i_station] = mat_single

        ##################################
        ############# Upper ##############
        ##################################
        xsec_nodeU = Float64.(numadIn.segments[i_station,numadIn.segments[i_station,:].>=0.0])


        seg_idxU = (numadIn.segments[i_station,:].>0.0)[2:end]
        n_laminaU = zeros(Int,sum(seg_idxU))
        # idx_le = argmin(abs.(numadIn.segments[i_station,:]))-1
        idx = 1
        for seg_idx = 1:numadIn.n_segments
            if seg_idxU[seg_idx] == true
                n_laminaU[idx] = length(numadIn.skin_seq[i_station,seg_idx].seq)
                idx += 1
            end
        end

        n_pliesU = zeros(Int,sum(n_laminaU))
        mat_lamU = zeros(Int,sum(n_laminaU))
        t_lamU = zeros(sum(n_laminaU)) #TODO: hook this into the optimization parameters and or the material properties
        tht_lamU = zeros(sum(n_laminaU)) #TODO: same with this
        idx = 1
        for seg_idx = 1:numadIn.n_segments
            if seg_idxU[seg_idx] == true # make sure we are using the upper segment
                for seq_idx = 1:length(numadIn.skin_seq[i_station,seg_idx].seq)
                    mat_idx = numadIn.skin_seq[i_station,seg_idx].seq[seq_idx]
                    mat_lamU[idx] = mat_idx
                    n_pliesU[idx] = numadIn.stack_layers[i_station,mat_idx]
                    t_lamU[idx] = ply_thickness[mat_idx] #n_pliesU[idx]*ply_thickness[mat_idx]
                    idx += 1
                end
            end
        end

        ##################################
        ############# Lower ##############
        ##################################
        xsec_nodeL = Float64.(abs.(reverse(numadIn.segments[i_station,numadIn.segments[i_station,:].<=0.0]))) #TODO: fix types and verify positive increasing is correct for precomp even on the bottom

        seg_idxL = (numadIn.segments[i_station,:].<=0.0)[2:end]
        n_laminaL = zeros(Int,sum(seg_idxL))
        idx = 1
        for seg_idx = 1:numadIn.n_segments
            if seg_idxL[seg_idx] == true
                n_laminaL[idx] = length(numadIn.skin_seq[i_station,seg_idx].seq)
                idx += 1
            end
        end

        n_pliesL = zeros(Int,sum(n_laminaL))
        mat_lamL = zeros(Int,sum(n_laminaL))
        t_lamL = zeros(sum(n_laminaL)) #TODO: hook this into the optimization parameters and or the material properties
        tht_lamL = zeros(sum(n_laminaL)) #TODO: same with this
        idx = 1
        for seg_idx = 1:numadIn.n_segments
            if seg_idxL[seg_idx] == true
                for seq_idx = 1:length(numadIn.skin_seq[i_station,seg_idx].seq)
                    mat_idx = numadIn.skin_seq[i_station,seg_idx].seq[seq_idx]
                    mat_lamL[idx] = mat_idx
                    n_pliesL[idx] = numadIn.stack_layers[i_station,mat_idx]
                    t_lamL[idx] = ply_thickness[mat_idx]#n_pliesL[idx]*ply_thickness[mat_idx]
                    idx += 1
                end
            end
        end

        ##################################
        ############# Web(s) #############
        ##################################

        loc_web = zeros(numadIn.n_web)
        n_laminaW = zeros(Int,numadIn.n_web)
        # println("You must define shear webs at each spanwise station, just set the ply thicknesses to zero if not desired")
        for web_idx = 1:numadIn.n_web
            idx_loc_web= numadIn.web_dp[i_station,web_idx].seq[1]+1
            loc_web[web_idx] = abs(numadIn.segments[i_station,idx_loc_web])
            n_laminaW[web_idx] = length(numadIn.web_seq[i_station,web_idx].seq)
        end

        # Now ensure that there aren't any airfoil points already where the webs are located
        xpc_filtered = [] #TODO: make this more efficient
        ypc_filtered = [] #TODO: make this more efficient
        for i_af = 1:length(xpc)
            alreadyPushed = false
            if length(loc_web)>=1
                for j_web = 1:length(loc_web)
                    if !isapprox(xpc[i_af],loc_web[j_web],atol = 1e-4) && alreadyPushed == false
                        push!(xpc_filtered,xpc[i_af])
                        push!(ypc_filtered,ypc[i_af])
                        alreadyPushed = true
                    end
                end
            else
                push!(xpc_filtered,xpc[i_af])
                push!(ypc_filtered,ypc[i_af])
            end
        end
        xpc_filtered = Float64.(xpc_filtered)
        ypc_filtered = Float64.(ypc_filtered)

        n_pliesW = zeros(Int,sum(n_laminaW))
        mat_lamW = zeros(Int,sum(n_laminaW))
        t_lamW = zeros(sum(n_laminaW)) #TODO: hook this into the optimization parameters and or the material properties
        tht_lamW = zeros(sum(n_laminaW)) #TODO: same with this
        idx = 1
        for web_idx = 1:numadIn.n_web
            for seq_idx = 1:length(numadIn.web_seq[i_station,web_idx].seq)
                mat_idx = numadIn.web_seq[i_station,web_idx].seq[seq_idx]
                mat_lamW[idx] = mat_idx
                n_pliesW[idx] = numadIn.stack_layers[i_station,mat_idx]
                t_lamW[idx] = ply_thickness[mat_idx] #n_pliesW[idx]*ply_thickness[mat_idx]
                idx += 1
            end
        end

        ########################################
        ## Create the Precomp Input Structure ##
        ########################################
        precompinput[i_station] = PreComp.Input(
        normalchord[i_station],
        -twist_d[i_station],-twistrate_d[i_station],
        leloc[i_station],xpc_filtered,ypc_filtered,
        e1,e2,g12,anu12,density,
        xsec_nodeU,n_laminaU,n_pliesU,t_lamU,tht_lamU,mat_lamU,
        xsec_nodeL,n_laminaL,n_pliesL,t_lamL,tht_lamL,mat_lamL,
        loc_web,n_laminaW,n_pliesW,t_lamW,tht_lamW,mat_lamW)

        # calculate composite properties: stiffness, mass, etc
        precompoutput[i_station] = PreComp.properties(precompinput[i_station])

    end

    return precompoutput,precompinput
end

function getSectPropsFromPreComp(usedUnitSpan,numadIn,precompoutput)
    # usedUnitSpan is node positions, as is numadIn.span, and the precomp calculations
    # create spline of the precomp output to be used with the specified span array
    len_pc = length(precompoutput)
    ei_flap = zeros(len_pc)
    ei_lag = zeros(len_pc)
    gj = zeros(len_pc)
    ea = zeros(len_pc)
    s_fl = zeros(len_pc)
    s_af = zeros(len_pc)
    s_al = zeros(len_pc)
    s_ft = zeros(len_pc)
    s_lt = zeros(len_pc)
    s_at = zeros(len_pc)
    x_sc = zeros(len_pc)
    y_sc = zeros(len_pc)
    x_tc = zeros(len_pc)
    y_tc = zeros(len_pc)
    mass = zeros(len_pc)
    flap_iner = zeros(len_pc)
    lag_iner = zeros(len_pc)
    tw_iner_d = zeros(len_pc)
    x_cm = zeros(len_pc)
    y_cm = zeros(len_pc)

    # extract the values from the precomp outputs
    for i_pc = 1:len_pc
        ei_flap[i_pc] = precompoutput[i_pc].ei_flap
        ei_lag[i_pc] = precompoutput[i_pc].ei_lag
        gj[i_pc] = precompoutput[i_pc].gj
        ea[i_pc] = precompoutput[i_pc].ea
        s_fl[i_pc] = precompoutput[i_pc].s_fl
        s_af[i_pc] = precompoutput[i_pc].s_af
        s_al[i_pc] = precompoutput[i_pc].s_al
        s_ft[i_pc] = precompoutput[i_pc].s_ft
        s_lt[i_pc] = precompoutput[i_pc].s_lt
        s_at[i_pc] = precompoutput[i_pc].s_at
        x_sc[i_pc] = precompoutput[i_pc].x_sc
        y_sc[i_pc] = precompoutput[i_pc].y_sc
        x_tc[i_pc] = precompoutput[i_pc].x_tc
        y_tc[i_pc] = precompoutput[i_pc].y_tc
        mass[i_pc] = precompoutput[i_pc].mass
        flap_iner[i_pc] = precompoutput[i_pc].flap_iner
        lag_iner[i_pc] = precompoutput[i_pc].lag_iner
        tw_iner_d[i_pc] = precompoutput[i_pc].tw_iner_d
        x_cm[i_pc] = precompoutput[i_pc].x_cm
        y_cm[i_pc] = precompoutput[i_pc].y_cm
    end

    # Now create the splines and sample them at the used span
    origUnitSpan = numadIn.span./numadIn.span[end]
    usedUnitSpan = usedUnitSpan./maximum(usedUnitSpan)
    ei_flap_used = FLOWMath.akima(origUnitSpan,ei_flap,usedUnitSpan)
    ei_lag_used = FLOWMath.akima(origUnitSpan,ei_lag,usedUnitSpan)
    gj_used = FLOWMath.akima(origUnitSpan,gj,usedUnitSpan)
    ea_used = FLOWMath.akima(origUnitSpan,ea,usedUnitSpan)
    s_fl_used = FLOWMath.akima(origUnitSpan,s_fl,usedUnitSpan)
    s_af_used = FLOWMath.akima(origUnitSpan,s_af,usedUnitSpan)
    s_al_used = FLOWMath.akima(origUnitSpan,s_al,usedUnitSpan)
    s_ft_used = FLOWMath.akima(origUnitSpan,s_ft,usedUnitSpan)
    s_lt_used = FLOWMath.akima(origUnitSpan,s_lt,usedUnitSpan)
    s_at_used = FLOWMath.akima(origUnitSpan,s_at,usedUnitSpan)
    x_sc_used = FLOWMath.akima(origUnitSpan,x_sc,usedUnitSpan)
    y_sc_used = FLOWMath.akima(origUnitSpan,y_sc,usedUnitSpan)
    x_tc_used = FLOWMath.akima(origUnitSpan,x_tc,usedUnitSpan)
    y_tc_used = FLOWMath.akima(origUnitSpan,y_tc,usedUnitSpan)
    mass_used = FLOWMath.akima(origUnitSpan,mass,usedUnitSpan)
    flap_iner_used = FLOWMath.akima(origUnitSpan,flap_iner,usedUnitSpan)
    lag_iner_used = FLOWMath.akima(origUnitSpan,lag_iner,usedUnitSpan)
    tw_iner_d_used = FLOWMath.akima(origUnitSpan,tw_iner_d,usedUnitSpan)
    x_cm_used = FLOWMath.akima(origUnitSpan,x_cm,usedUnitSpan)
    y_cm_used = FLOWMath.akima(origUnitSpan,y_cm,usedUnitSpan)

    ac_used = FLOWMath.akima(origUnitSpan,numadIn.aerocenter,usedUnitSpan)
    twist_d_used = FLOWMath.akima(origUnitSpan,numadIn.twist_d,usedUnitSpan)
    chord_used = FLOWMath.akima(origUnitSpan,numadIn.chord,usedUnitSpan)

    sectionPropsArray = Array{OWENS.SectionPropsArray, 1}(undef, length(usedUnitSpan)-1)

    for i=1:length(usedUnitSpan)-1

        #structural properties
        ac = -([ac_used[i], ac_used[i+1]].-0.5)
        twist_d=[twist_d_used[i], twist_d_used[i+1]] # indegrees #TODO: update all angles to be in radians unless explicitely indicated
        rhoA = [mass_used[i], mass_used[i+1]]
        EIyy = [ei_flap_used[i], ei_flap_used[i+1]]
        EIzz = [ei_lag_used[i], ei_lag_used[i+1]]
        if (minimum(abs.(EIyy .- EIzz)) < 1.0e-3)
            EIzz = EIzz.*1.0001
        end
        GJ = [gj_used[i], gj_used[i+1]]
        EA = [ea_used[i], ea_used[i+1]]

        rhoIyy = [flap_iner_used[i], flap_iner_used[i+1]]
        rhoIzz = [lag_iner_used[i], lag_iner_used[i+1]]
        rhoJ = [flap_iner_used[i]+lag_iner_used[i+1], flap_iner_used[i]+lag_iner_used[i+1]]
        zcm = [x_cm_used[i], x_cm_used[i+1]]
        ycm = [y_cm_used[i], y_cm_used[i+1]]
        a = [y_tc_used[i], y_tc_used[i+1]]

        #coupling factors
        EIyz = [0.0, 0.0]
        alpha1 = [0.0, 0.0] #This is always 0 in the element file, and it is unclear what it is used for since I can't find it being used in the code
        alpha2 = [0.0, 0.0]
        alpha3 = [0.0, 0.0]
        alpha4 = [0.0, 0.0]
        alpha5 = [0.0, 0.0]
        alpha6 = [0.0, 0.0]
        rhoIyz = [0.0, 0.0]
        b = [0.0, 0.0]
        a0 = [2*pi, 2*pi] #TODO: shouldn't the lift slope for a cylinder be 0? So this should depend on the airfoil used.
        aeroCenterOffset = [0.0, 0.0]

        #TODO: not all of the precomp data is used, need to include it for a more accurate solution
        sectionPropsArray[i] = SectionPropsArray(ac,twist_d,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

    end

    #TODO: why are we doing this post-calc, why not during???
    for i=1:length(usedUnitSpan)-1

            sectionPropsArray[i].b = 0.5.*[chord_used[i], chord_used[i+1]] #element semi chord
            # sectionPropsArray[i].a0 = [bladeData[i,12], bladeData[i+1,12]]         #element lift curve slope (needed for flutter analysis) TODO: enable coupling between actual airfoil lift slope

            #convert "a" to semichord fraction aft of halfchord
            sectionPropsArray[i].a = (sectionPropsArray[i].a .+ 0.25*2*sectionPropsArray[i].b .- sectionPropsArray[i].b)./sectionPropsArray[i].b

            #convert "ac" to semichord fraction foreward of halfchord TODO: why are we doing it this way???
            sectionPropsArray[i].ac = sectionPropsArray[i].ac.*2

            #physical aero center offset from elastic axis
            sectionPropsArray[i].aeroCenterOffset = sectionPropsArray[i].ac .* sectionPropsArray[i].b .- sectionPropsArray[i].a
    end

    println("EIyz, rhoIyz deactivated") #TODO: why is this, especially when I believe precomp calculates them
    return sectionPropsArray

end
