
# mutable struct Component
#     name
#     x
#     y
#     z
#     conn
#     joint_connection_to
# end


function create_mesh(;Htwr_base = 15.0,
    Htwr_blds = 147.148-15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    AD15hubR = 2.0,
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    nselem = 4,
    nguyelem = 20,
    Nguy_sets = 3,
    guyanchor_radius = 10.0, #meters
    guy_twr_mountpoint = [0.5,0.95], # This puts struts at top and bottom, as a fraction of the blade position
    strut_twr_mountpoint = [0.125,0.5,0.95], # This puts struts at top and bottom, as a fraction of the blade position
    strut_bld_mountpoint = [0.25,0.5,0.75], # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    bshapey = zeros(nbelem+1), # but magnitude for this is relevant
    angularOffset = 0.0, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = false,
    verbosity = 0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr = true,
    isHAWT=false)

    Nguy_vert = length(guy_twr_mountpoint)
    Nstrut = length(strut_bld_mountpoint)

    if length(strut_bld_mountpoint) != length(strut_twr_mountpoint)
        error("strut_twr_mountpoint must be the same length as strut_bld_mountpoint")
    end

    ##################################
    #             _
    #           /_|_\
    #          |  |  )
    #           \-|-/
    ####################################
    ###------------Tower--------------##
    ####################################
    mesh_z = collect(LinRange(0,Htwr_blds+Htwr_base,ntelem+1))

    # Insert mount point base
    mesh_z = sort([mesh_z;Htwr_base])
    t_botidx = findall(x->isapprox(x,Htwr_base,atol=1e-5),mesh_z)[1]#:nblade]

    t2s_idx = zeros(Int, Nstrut)
    # Insert strut mount points
    for istrut = 1:Nstrut
        if maximum(Hbld*strut_twr_mountpoint[istrut]+Htwr_base .== mesh_z) # if we are at exactly an existing node, then offset our mount point
            @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of tower elements to space out the mesh"
            strut_twr_mountpoint[istrut] += 1e-6
        end
        mesh_z = sort([mesh_z;Hbld*strut_twr_mountpoint[istrut]+Htwr_base])

        # pick out the strut mounting indices
        t2s_idx[istrut] = findall(x->isapprox(x,Hbld*strut_twr_mountpoint[istrut]+Htwr_base,atol=1e-5*Hbld),mesh_z)[1]
    end

    if Nguy_sets>0
        t2g_idx = zeros(Int, Nguy_vert)
        # Insert guy mount points
        for iguy = 1:Nguy_vert
            if maximum(Hbld*guy_twr_mountpoint[iguy]+Htwr_base .== mesh_z) # if we are at exactly an existing node, then offset our mount point
                @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of tower elements to space out the mesh"
                guy_twr_mountpoint[iguy] += 1e-6
            end
            mesh_z = sort([mesh_z;Hbld*guy_twr_mountpoint[iguy]+Htwr_base])

            # pick out the strut mounting indices
            t2g_idx[iguy] = findall(x->isapprox(x,Hbld*guy_twr_mountpoint[iguy]+Htwr_base,atol=1e-5*Hbld),mesh_z)[1]
        end
    end

    # Create the x and y components of same size as mesh_z now that the strut mount points are inserted
    mesh_x = zero(mesh_z)
    mesh_y = zero(mesh_z)

    # pick out the tower top index
    t_topidx = length(mesh_z)

    # intra-tower connectivity
    conn = zeros(length(mesh_z)-1,2)
    conn[:,1] = collect(1:length(mesh_z)-1)
    conn[:,2] = collect(2:length(mesh_z))

    #####################################
    ###------------Blades--------------##
    #####################################

    #connection points on tower are simply the bottom of the tower offset connecting to the bottom of the blades, which is jointed if Darrieus in the joint matrix below
    bld_Z = collect(LinRange(0.0,Hbld,nbelem+1))

    # Insert bottom strut mount point
    for istrut = 1:Nstrut
        if maximum(Hbld*strut_bld_mountpoint[istrut] .== bld_Z) # if we are at exactly an existing node, then offset our mount point
            @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of elements to space out the mesh"
            strut_bld_mountpoint[istrut] -= 1e-6
        end
        bld_Z = sort([bld_Z;Hbld*strut_bld_mountpoint[istrut]])
    end

    if bshapex == zeros(nbelem+1)
        bld_Y = R.*(1.0.-4.0.*(bld_Z/Hbld.-.5).^2)
    else
        # Ensure the blade shape conforms to the turbine height and radius specs
        bshapex = R .* bshapex./maximum(bshapex)
        bshapez = Hbld .* bshapez./maximum(bshapez)
        bld_Y = OWENS.safeakima(bshapez,bshapex,bld_Z)
    end

    if bshapey == zeros(nbelem+1)
        bld_X = zero(bld_Y)
    else
        bld_X = OWENS.safeakima(bshapez,bshapey,bld_Z)
    end

    # AeroDyn Compatability
    AD15bldNdIdxRng = zeros(Int64,0,2)

    bld_Z .+= Htwr_base

    b_Z = []
    b_X = []
    b_Y = []
    # Now using standard VAWT convention, blade 1 is zero degrees at top dead center, or North/Y+
    # and they are offset counter clockwise
    b_topidx = zeros(Int,nblade)
    b_botidx = zeros(Int,nblade) .+ length(mesh_z)
    conn_b = zeros(length(bld_Z)-1,2)
    for ibld = 1:nblade
        myangle = (ibld-1)*2.0*pi/nblade + angularOffset
        b_Z = [b_Z;bld_Z]
        b_X = [b_X;bld_Y.*sin(myangle).+bld_X.*cos(myangle)]
        b_Y = [b_Y;-bld_X.*sin(myangle).+bld_Y.*cos(myangle)]

        # Element joint indices
        b_botidx[ibld] = length(mesh_z)+1 + length(bld_Z)*(ibld-1)
        b_topidx[ibld] = length(mesh_z)+1 + length(bld_Z)*ibld-1

        # Intraconnectivity
        conn_b[:,1] = collect(b_botidx[ibld]:1:b_topidx[ibld]-1)
        conn_b[:,2] = collect(b_botidx[ibld]+1:1:b_topidx[ibld])
        conn = [conn;conn_b]

        if AD15_ccw #Clockwise, the blades roots are at the top, trailing edge is always positive y
            AD15bldNdIdxRng = [AD15bldNdIdxRng; b_topidx[ibld] b_botidx[ibld]]    # top of blade is root 
        elseif !(AD15_ccw) #Clockwise, the blades roots are at the bottom
            AD15bldNdIdxRng = [AD15bldNdIdxRng; b_botidx[ibld] b_topidx[ibld]]    # bottom of blade is root
        end
    end

    # pick out the strut mounting indices
    b2s_idx = zeros(Int,nblade,Nstrut)
    for istrut = 1:Nstrut
        b2s_idx[:,istrut] = findall(x->x==Hbld*strut_bld_mountpoint[istrut]+Htwr_base,b_Z)[1:nblade] .+ length(mesh_z)
    end

    # Add to the mesh
    mesh_z = [mesh_z;b_Z]
    mesh_x = [mesh_x;b_X]
    mesh_y = [mesh_y;b_Y]

    #####################################
    ###------------Struts--------------##
    #####################################

    function createstrut(sstartidx,sendidx,mesh_x,mesh_y,mesh_z,conn,AD15hubR,nel;sxend=nothing,syend=nothing,szend=nothing)

        sxstart = mesh_x[sstartidx]
        systart = mesh_y[sstartidx]
        szstart = mesh_z[sstartidx]

        if isnothing(sxend)
            sxend = mesh_x[sendidx]
            syend = mesh_y[sendidx]
            szend = mesh_z[sendidx]
        end

        # Now draw the lines
        s_x = collect(LinRange(sxstart,sxend,nel+1))
        s_y = collect(LinRange(systart,syend,nel+1))
        s_z = collect(LinRange(szstart,szend,nel+1))

        hubIdx = 1
        if AD15hubR > 1e-6
            lenXY = sqrt((sxend - sxstart)^2 + (syend - systart)^2)   # strut length in XY
            minR2 = lenXY 
            for i = 1:nel+1  # step through to find closest point to hub radius on x-y plane
                R2 = AD15hubR - sqrt((s_x[i] - sxstart)^2 + (s_y[i] - systart)^2)
                if abs(R2) < abs(minR2)
                    hubIdx = i
                    minR2 = R2
                end
            end
            R_temp = minR2

            s_x[hubIdx] = s_x[hubIdx] + R_temp/lenXY*(sxend-sxstart)
            s_y[hubIdx] = s_y[hubIdx] + R_temp/lenXY*(syend-systart)
            s_z[hubIdx] = s_z[hubIdx] + R_temp/lenXY*(szend-szstart)

            if verbosity>0
                println("Hub crossing at idx $hubIdx and radially at $R_temp with AD15 hub radius of $AD15hubR")
                println("Moving strut point from [$(s_x[hubIdx]),$(s_y[hubIdx]),$(s_z[hubIdx])] to [$(s_x[hubIdx]),$(s_y[hubIdx]),$(s_z[hubIdx])]")
            end

        end

        hubIdx += length(mesh_z)

        # joint connections
        s2b_idx_internal = length(mesh_z)+1
        s2t_idx_internal = s2b_idx_internal+length(s_z)-1

        # and add to the mesh
        mesh_x = [mesh_x;s_x]
        mesh_y = [mesh_y;s_y]
        mesh_z = [mesh_z;s_z]

        # Intraconnectivity
        conn_s = zeros(nel,2)
        conn_s[:,1] = collect(s2b_idx_internal:1:s2t_idx_internal-1)
        conn_s[:,2] = collect(s2b_idx_internal+1:1:s2t_idx_internal)
        conn = [conn;conn_s]

        return s2b_idx_internal,s2t_idx_internal,mesh_x,mesh_y,mesh_z,conn,hubIdx
    end

    #Connect from the tower to the blades
    # For each blade, find the mounting location and draw a line
    s2b_idx = zeros(Int,nblade,Nstrut)
    s2t_idx = zeros(Int,nblade,Nstrut)

    # Bottom Struts
    for istrut = 1:Nstrut
        for ibld = 1:nblade
            s2t_idx[ibld,istrut],s2b_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,hubIsectIdx = createstrut(t2s_idx[istrut],b2s_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,AD15hubR,nselem)
            
            AD15bldNdIdxRng = [AD15bldNdIdxRng; hubIsectIdx  s2b_idx[ibld,istrut]] #AD15 struts always start at hub regardless of rotation, but watch out for airfoil orientation!
        end
    end

    if Nguy_sets>0
        #Connect from the tower to the ground
        # For each guy set, find the mounting location and draw a line
        g2ground_idx = zeros(Int,Nguy_sets,Nguy_vert)
        g2t_idx = zeros(Int,Nguy_sets,Nguy_vert)
        # Guy Wires
        for iguy_vert = 1:Nguy_vert
            for iguy_set = 1:Nguy_sets
                myangle = (iguy_set-1)*2.0*pi/Nguy_sets + angularOffset
                
                g2t_idx[iguy_set,iguy_vert],g2ground_idx[iguy_set,iguy_vert],mesh_x,mesh_y,mesh_z,conn,_ = createstrut(t2g_idx[iguy_vert],nothing,mesh_x,mesh_y,mesh_z,conn,0.0,nguyelem;
                sxend = guyanchor_radius.*sin(myangle),
                syend = guyanchor_radius.*cos(myangle),
                szend=0.0)
            end
        end
    end
    #######################################
    ###   Cleanup/Derived parameters    ###
    #######################################

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower, treat struts like blades
    meshtype = zeros(Int,numEl)

    # Find elnum associated with t_topidx
    topel_idx = findall(x->x==t_topidx,conn[:,2])
    meshtype[1:topel_idx[1]] .= 1 #Tower

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,1+nblade+Nstrut*nblade+Nguy_sets*Nguy_vert) #tower, blades, and struts

    meshSeg[1] = topel_idx[1]
    meshSeg[2:nblade+1] .= nbelem+Nstrut #+Nstrut for strut mount points
    meshSeg[nblade+2:end] .= nselem

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z))
    structuralNodeNumbers = zeros(nblade,length(bld_Z))
    structuralElNumbers = zeros(nblade,length(bld_Z))

    for iblade = 1:nblade

        # Normalized Span
        span_len = bld_Z.-Htwr_base#sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-Htwr_base).^2.0)
        structuralSpanLocNorm[iblade,:] = span_len./maximum(span_len)

        # Node Numbers
        structuralNodeNumbers[iblade,:] = collect(b_botidx[iblade]:b_topidx[iblade])

        # Element Numbers
        structuralElNumbers[iblade,:] = structuralNodeNumbers[iblade,:].-iblade
        structuralElNumbers[iblade,end] = -1 #TODO: figure out why this is in the original OWENS setup and if it is used
    end

    mymesh = OWENSFEA.Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,Int.(conn),meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ######################################
    ####----------Joint Matrix----------##
    ######################################

    # Connect Tower Top to Blades bottom, then each cable to each blade bottom
    # Then each cable to each blade top, then the latter two blade tops to the first

    njoint = nblade*2+Nstrut*nblade*2+Nguy_sets*Nguy_vert #create the full array and then pull out zeros below if H-VAWT where the blades aren't connected to the tower.
    jointconn = zeros(Int,njoint,2)
    jointtype = zeros(njoint)
    ijoint = 0
    for ibld = 1:nblade
        if connectBldTips2Twr || isHAWT
            # connect tower to blades
            ijoint +=1
            jointconn[ijoint,:] = [t_botidx b_botidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect tower to each strut
            ijoint +=1
            jointconn[ijoint,:] = [t2s_idx[istrut] s2t_idx[ibld,istrut]]
        end

        if connectBldTips2Twr && !isHAWT
            # connect tower to blades tops
            ijoint +=1
            jointconn[ijoint,:] = [t_topidx b_topidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect strut to blade bottom
            ijoint +=1
            jointconn[ijoint,:] = [s2b_idx[ibld,istrut] b2s_idx[ibld,istrut]]
        end
    end

    for iguy_set = 1:Nguy_sets
        for iguy_vert = 1:Nguy_vert
            # connect strut to blade bottom
            ijoint +=1
            jointconn[ijoint,:] = [t2g_idx[iguy_vert] g2t_idx[iguy_set,iguy_vert]]
        end
    end


    # Reduce the matrix based on if the blades got connected or not, throwing out all the zero rows
    bitlogic = jointconn[:,1] .!= 0.0
    jointconn = jointconn[bitlogic,:]
    jointtype = jointtype[bitlogic]

   

    njoint = length(jointconn[:,1]) # reset the number of joints
    myort = OWENS.calculateElementOrientation(mymesh)
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,1],conn[:,1]) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,1],conn[:,2])
        end
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2],conn[:,2])
        end
        Psi_d_joint[jnt] = myort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = myort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

    #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn jointtype zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    # Blade and strut starting and ending node and element numbers
    AD15bldElIdxRng = zeros(Int64,0,2)
    for i = 1:size(AD15bldNdIdxRng,1)
        if AD15bldNdIdxRng[i,2] > AD15bldNdIdxRng[i,1]  # ascending order
            idx1 = findfirst(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,1])
            idx2 = findfirst(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
        else    # upside down oriented blade
            idx1 = findlast(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,2])
            idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,1])
        end

        if isnothing(idx2)
            idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
        end
        AD15bldElIdxRng = [AD15bldElIdxRng; idx1 idx2]
    end
    custom_mesh_outputs = (g2ground_idx)
    return mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs
end