# Meshing Rules of Thumb:
# Always start from center and go outwards
# Blades start from bottom and go up
# Joints start from center and go outwards

#TODO: write file output to match OWENS original files
#TODO: especially write output for NUMAD excel input

# mesh for a single continuous beam
# Rotating beam case from Finite Element Solution of Nonlinear Intrinsic Equations for Curved Composite Beams by Hodges, Shang, and Cesnic
# ______
#       \
function mesh_beam(;L1 = 31.5, #first section of beam length
    L2 = 6.0, #second section of beam length
    Nelem1 = 13,
    Nelem2 = 3,
    angleD = 45.0, # angle of second section of beam relative to first (0 is straight)
    zeroOffset = 2.5,
    vertical=true,
    x_shape = nothing,
    y_shape = nothing,
    z_shape = nothing) #offset from 0 before first beam begins

    angle = angleD/360*2*pi

    if isnothing(x_shape)
        # First Section
        mesh_x1 = collect(LinRange(zeroOffset,zeroOffset+L1,Nelem1+1))
        mesh_y1 = zero(mesh_x1)
        mesh_z1 = zero(mesh_x1)

        # Angled Section
        mesh_x2 = zeroOffset+L1.+cos(angle).*collect(LinRange(0.0,L2,Nelem2+1))
        mesh_y2 = -sin(angle).*collect(LinRange(0.0,L2,Nelem2+1))
        mesh_z2 = zero(mesh_x2)
    else

        # First Section
        mesh_x1 = collect(LinRange(zeroOffset,zeroOffset+L1,Nelem1+1))
        mesh_y1 = safeakima(x_shape,y_shape,mesh_x1)
        mesh_z1 = safeakima(x_shape,z_shape,mesh_x1)

        # Angled Section
        mesh_x2 = zeroOffset+L1.+collect(LinRange(0.0,L2,Nelem2+1))
        mesh_y2 = safeakima(x_shape,y_shape,mesh_x2)
        mesh_z2 = safeakima(x_shape,z_shape,mesh_x2)
    end

    # intra-beam connectivity
    conn1 = zeros(length(mesh_z1)-1,2)
    conn1[:,1] = collect(1:length(mesh_z1)-1)
    conn1[:,2] = collect(2:length(mesh_z1))
    conn2 = zeros(length(mesh_z2)-1,2)
    conn2[:,1] = length(mesh_z1).+collect(1:length(mesh_z2)-1)
    conn2[:,2] = length(mesh_z1).+collect(2:length(mesh_z2))

    conn = [conn1;conn2]

    if vertical
        mesh_z = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_x = [mesh_z1;mesh_z2]
    else
        mesh_x = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_z = [mesh_z1;mesh_z2]
    end

    # #space out the mesh numerically to avoid numerical issues
    # for i = 1:length(mesh_x)-1
    #     if isapprox(mesh_x[i],mesh_x[i+1];atol = 1e-6)
    #         mesh_x[i+1] -= 1e-4
    #     end
    # end

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:end] .= 0 #Tower

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,2)

    meshSeg[1] = Nelem1
    meshSeg[2] = Nelem2

    # Not used for the beam case
    structuralSpanLocNorm = []
    structuralNodeNumbers = []
    structuralElNumbers = []
    # end

    mymesh = OWENSFEA.Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,Int.(conn),meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ####################################
    ##----------Joint Matrix----------##
    ####################################

    #Connect L1 tip to L2 base
    jointconn = [length(mesh_x1) length(mesh_x1)+1]

    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh) #TODO: consider getting rid of ort struct for simplification since it isn't used hardly at all
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

# mesh for a single continuous beam
# Rotating beam case from Finite Element Solution of Nonlinear Intrinsic Equations for Curved Composite Beams by Hodges, Shang, and Cesnic
# ______
#       \
function mesh_beam_centered(;L1 = 6.0, #first section of beam length
    L2 = -6.0, #second section of beam length
    Nelem1 = 2,
    Nelem2 = 2,
    angleD = 0.0, # angle of second section of beam relative to first (0 is straight)
    zeroOffset = 0.0,
    shaftOffset = 0.25,
    vertical=false,
    x_shape = nothing,
    y_shape = nothing,
    z_shape = nothing) #offset from 0 before first beam begins

    angle = angleD/360*2*pi

    if isnothing(x_shape)
        # First Section
        mesh_x1 = collect(LinRange(zeroOffset,zeroOffset+L1,Nelem1+1))
        mesh_y1 = zero(mesh_x1)
        mesh_z1 = zero(mesh_x1)

        # Angled Section
        mesh_x2 = zeroOffset.+collect(LinRange(0.0,L2,Nelem2+1))
        mesh_y2 = -sin(angle).*collect(LinRange(0.0,L2,Nelem2+1))
        mesh_z2 = zero(mesh_x2)
    else

        # First Section
        mesh_x1 = collect(LinRange(zeroOffset,zeroOffset+L1,Nelem1+1))
        mesh_y1 = safeakima(x_shape,y_shape,mesh_x1)
        mesh_z1 = safeakima(x_shape,z_shape,mesh_x1)

        # Angled Section
        mesh_x2 = zeroOffset+L1.+collect(LinRange(0.0,L2,Nelem2+1))
        mesh_y2 = safeakima(x_shape,y_shape,mesh_x2)
        mesh_z2 = safeakima(x_shape,z_shape,mesh_x2)
    end

    # Push a shaft onto the mesh

    mesh_x1 = [0.0;mesh_x1]
    mesh_y1 = [0.0;mesh_y1]
    mesh_z1 = [0.0;mesh_z1]

    # intra-beam connectivity
    conn1 = zeros(length(mesh_z1)-1,2)
    conn1[:,1] = collect(1:length(mesh_z1)-1)
    conn1[:,2] = collect(2:length(mesh_z1))
    conn2 = zeros(length(mesh_z2)-1,2)
    conn2[:,1] = length(mesh_z1).+collect(1:length(mesh_z2)-1)
    conn2[:,2] = length(mesh_z1).+collect(2:length(mesh_z2))

    conn = [conn1;conn2]

    if vertical
        mesh_z = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_x = [mesh_z1;mesh_z2]
        mesh_x[2:end] .+= shaftOffset
    else
        mesh_x = [mesh_x1;mesh_x2]
        mesh_y = [mesh_y1;mesh_y2]
        mesh_z = [mesh_z1;mesh_z2]
        mesh_z[2:end] .+= shaftOffset
    end

    # #space out the mesh numerically to avoid numerical issues
    # for i = 1:length(mesh_x)-1
    #     if isapprox(mesh_x[i],mesh_x[i+1];atol = 1e-6)
    #         mesh_x[i+1] -= 1e-4
    #     end
    # end

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = ones(Int,numEl) # tower
    meshtype[2:end] .= 0 #Blade

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,3)

    meshSeg[1] = 1 # shaft
    meshSeg[2] = Nelem1
    meshSeg[3] = Nelem2

    # used for this beam case
    structuralSpanLocNorm = collect([LinRange(0,1,Nelem1+1);;LinRange(0,-1,Nelem1+1)]')
    structuralNodeNumbers = collect([1:Nelem1+1 Nelem1+2:length(mesh_z)-1]') .+ 1
    structuralElNumbers = collect([1:Nelem1 Nelem1+1:Nelem1+Nelem2]') .+ 1
    # end

    mymesh = OWENSFEA.Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,Int.(conn),meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ####################################
    ##----------Joint Matrix----------##
    ####################################

    #Connect L2 base to L1 base
    jointconn = [2 length(mesh_x1)+1]

    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh) #TODO: consider getting rid of ort struct for simplification since it isn't used hardly at all
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
    myjoint = [Float64.(1:1:njoint) Int.(jointconn) zeros(njoint) zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    return mymesh, ort, myjoint
end

"""
mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = create_mesh_struts(;Htwr_base = 15.0,
Htwr_blds = 147.148-15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    AD15hubR = 2.0,
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    nselem = 4,
    strut_twr_mountpoint = [0.01,0.5,0.9],
    strut_bld_mountpoint = [0.01,0.5,0.9],
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    bshapey = zeros(nbelem+1), # but magnitude for this is relevant
    angularOffset = 0.0, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = false,
    verbosity=0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr =true)

Standard Mesh Matching 5MW, 34m configurations

#Inputs
* `Htwr_base::float`: height of tower before blades attach (m)
* `Htwr_blds::float`: height of the tower once the blades attach (m)
* `Hbld::float`: blade height (m)
* `R::float`: bade radius (m)
* `nblade::int`: number of blades
* `ntelem::int`: number of tower elements
* `nbelem::int`: number of blade elements
* `nselem::int`: number of strut elements
* `strut_twr_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the tower # This puts struts at top and bottom, as a fraction of the blade position
* `strut_bld_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the blade # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
* `bshapex::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapez::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapey::Array{<:float}`: Blade shape, magnitude IS relevant #TODO: resolve this
* `angularOffset::float`: (rad) angular offset of mesh generation, typically used to match CACTUS input.  Value of 0 puts blade 1 at the "north" position and the others populate counterclockwise when looking down
* `AD15_ccw::boolean`: Use AD15 convention of VAWT counter-clockwise with blade root at top (blade points down)
* `AD15hubR::float`: AD15 has a hub radius, so the struts do not go all the way to the center of the axis of rotation, while the structural mesh does.
# `verbosity::int`: 0 nothing, 1 basic, 2 lots: amount of printed information 
* `connectBldTips2Twr::book`: True for Darrieus style, false for H-VAWT, but the blade shapes should be appropriate

#Outputs
* `mymesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `myort::OWENSFEA.Ort`: see ?OWENSFEA.Ort
* `myjoint:Array{<:float}`: see ?OWENSFEA.FEAModel.joint
* `AD15bldNdIdxRng`: indices for start and end of all blades for AD15 (includes struts).  Note that strut start nodes may be inside the strut (strut connects to tower, AD15 blade connects to hub wich is a few nodes away from tower)
* `AD15bldElIdxRng`: range of elements for start and end of all AD15 blades (includes struts)
"""
function create_mesh_struts(;Htwr_base = 15.0,
    Htwr_blds = 147.148-15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    AD15hubR = 2.0,
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    nselem = 4,
    strut_twr_mountpoint = [0.125,0.5,0.95], # This puts struts at top and bottom, as a fraction of the blade position
    strut_bld_mountpoint = [0.25,0.5,0.75], # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    bshapey = zeros(nbelem+1), # but magnitude for this is relevant
    angularOffset = 0.0, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = false,
    verbosity = 0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr = true)

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
            @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of elements to space out the mesh"
            strut_twr_mountpoint[istrut] += 1e-6
        end
        mesh_z = sort([mesh_z;Hbld*strut_twr_mountpoint[istrut]+Htwr_base])

        # pick out the strut mounting indices
        t2s_idx[istrut] = findall(x->isapprox(x,Hbld*strut_twr_mountpoint[istrut]+Htwr_base,atol=1e-5*Hbld),mesh_z)[1]
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
        bld_Y = safeakima(bshapez,bshapex,bld_Z)
    end

    if bshapey == zeros(nbelem+1)
        bld_X = zero(bld_Y)
    else
        bld_X = safeakima(bshapez,bshapey,bld_Z)
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

    function createstrut(sstartidx,sendidx,mesh_x,mesh_y,mesh_z,conn,AD15hubR)

        sxstart = mesh_x[sstartidx]
        systart = mesh_y[sstartidx]
        szstart = mesh_z[sstartidx]

        sxend = mesh_x[sendidx]
        syend = mesh_y[sendidx]
        szend = mesh_z[sendidx]

        # Now draw the lines
        s_x = collect(LinRange(sxstart,sxend,nselem+1))
        s_y = collect(LinRange(systart,syend,nselem+1))
        s_z = collect(LinRange(szstart,szend,nselem+1))

        hubIdx = 1
        if AD15hubR > 1e-6
            lenXY = sqrt((sxend - sxstart)^2 + (syend - systart)^2)   # strut length in XY
            minR2 = lenXY 
            for i = 1:nselem+1  # step through to find closest point to hub radius on x-y plane
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
        conn_s = zeros(nselem,2)
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
            s2t_idx[ibld,istrut],s2b_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,hubIsectIdx = createstrut(t2s_idx[istrut],b2s_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,AD15hubR)
            
            AD15bldNdIdxRng = [AD15bldNdIdxRng; hubIsectIdx  s2b_idx[ibld,istrut]] #AD15 struts always start at hub regardless of rotation, but watch out for airfoil orientation!
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

    meshSeg = zeros(Int,1+nblade+Nstrut*nblade) #tower, blades, and struts

    if connectBldTips2Twr == false
        meshSeg[1] = ntelem+Nstrut+1
    else
        meshSeg[1] = ntelem+Nstrut*nblade-nblade+1 #connects at top node
    end
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

    njoint = nblade*2+Nstrut*nblade*2 #create the full array and then pull out zeros below if H-VAWT where the blades aren't connected to the tower.
    jointconn = zeros(Int,njoint,2)
    jointtype = zeros(njoint)
    for ibld = 1:nblade
        if connectBldTips2Twr
            # connect tower to blades
            jointconn[ibld,:] = [t_botidx b_botidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect tower to each strut
            jointconn[ibld+nblade*istrut,:] = [t2s_idx[istrut] s2t_idx[ibld,istrut]]
        end

        if connectBldTips2Twr
            # connect tower to blades tops
            jointconn[ibld+nblade*(Nstrut+1),:] = [t_topidx b_topidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect strut to blade bottom
            jointconn[ibld+nblade*(Nstrut+1+istrut),:] = [s2b_idx[ibld,istrut] b2s_idx[ibld,istrut]]
        end
    end

    # Reduce the matrix based on if the blades got connected or not, throwing out all the zero rows
    bitlogic = jointconn[:,1] .!= 0.0
    jointconn = jointconn[bitlogic,:]
    jointtype = jointtype[bitlogic]

   

    njoint = length(jointconn[:,1]) # reset the number of joints
    myort = calculateElementOrientation(mymesh)
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
    return mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng
end

"""
create_arcus_mesh(;Htwr_base = 15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    ncelem = 4,
    strut_mountpoint = 0.01, # This puts struts at top and bottom
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    joint_type = 0,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0)


ARCUS mesh configuration: no tower between blades, no struts, but cables from top center attaching to specified blade mount point at base

#Inputs
* `Htwr_base::float`: height of tower before blades attach (m)
* `Hbld::float`: blade height (m)
* `R::float`: bade radius (m)
* `nblade::int`: number of blades
* `ntelem::int`: number of tower elements
* `nbelem::int`: number of blade elements
* `ncelem::int`: number of strut elements
* `c_mount_ratio::float`: factor of blade height where the struts attach on both top and bottom
* `bshapex::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapez::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `joint_type::int`: 0 is fixed, 1 is about x-axis, 2 is y-axis, etc
* `cables_connected_to_blade_base::bool`: = true,
* `angularOffset::float`: (rad) angular offset of mesh generation, typically used to match CACTUS input.  Value of 0 puts blade 1 at the "north" position and the others populate counterclockwise when looking down

#Outputs
* `mymesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `ort::OWENSFEA.Ort`: see ?OWENSFEA.Ort
* `myjoint:Array{<:float}`: see ?OWENSFEA.FEAModel.joint

"""
function create_arcus_mesh(;
    Htwr_base = 15.0, #tower height before blades attach
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nblade = 3,
    ntelem = 4, #tower elements
    nbelem = 30, #blade elements
    ncelem = 10,  #cable elements
    c_mount_ratio = 0.05, #fraction of blade where the cable mounts
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    joint_type = 0,
    AD15_ccw = false,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0)

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
    mesh_z = collect(LinRange(0,Htwr_base,ntelem+1))
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
    bld_Z = collect(LinRange(0.0,Hbld,nbelem+1))

    # The cable connections must be added though
    if maximum(Hbld*c_mount_ratio .== bld_Z) # if we are at exactly an existing node, then offset our mount point
        c_mount_ratio += 1e-6
    end
    bld_Z = sort([bld_Z;Hbld*c_mount_ratio])

    if bshapex == zeros(nbelem+1)
        bshapex = R .* bshapex./maximum(bshapex)
        bshapez = Hbld .* bshapez./maximum(bshapez)
        bld_Y = R.*(1.0.-4.0.*(bld_Z/Hbld.-.5).^2)
    else
        # Ensure the blade shape conforms to the turbine height and radius specs
        bld_Y = safeakima(bshapez,bshapex,bld_Z)
    end
    bld_X = zero(bld_Y)

    bld_Z .+= Htwr_base

    # AeroDyn Compatability
    AD15bldNdIdxRng = zeros(Int64,0,2)

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
        b_Y = [b_Y;bld_X.*sin(myangle).+bld_Y.*cos(myangle)]

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

    # Add to the mesh
    mesh_z = [mesh_z;b_Z]
    mesh_x = [mesh_x;b_X]
    mesh_y = [mesh_y;b_Y]

    # pick out the blade to cable mounting indices
    b2c_botidx = findall(x->x==Hbld*c_mount_ratio+Htwr_base,mesh_z)[1:nblade]


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
        czend = Htwr_base + Hbld

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

    # #space out the mesh numerically to avoid numerical issues
    # for i = 1:length(mesh_z)-1
    #     if isapprox(mesh_z[i],mesh_z[i+1];atol = 1e-6)
    #         mesh_z[i+1] -= 1e-4
    #     end
    # end

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
    meshSeg = zeros(Int,1+nblade+ncable) #tower, blades, and cables

    meshSeg[1] = ntelem
    meshSeg[2:nblade+1] .= nbelem+1 # +1 for the cable mount
    meshSeg[nblade+2:end] .= ncelem

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z)) # +1 for the cable mount
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

    if cables_connected_to_blade_base
        jointconn = zeros(Int,nblade*2-1+ncable*2,2)
        for ibld = 1:nblade
            # connect tower to blades
            jointconn[ibld,:] = [t_topidx b_botidx[ibld]]

            # connect cables to blades bases
            jointconn[ibld+nblade,:] = [b2c_botidx[ibld] c2b_botidx[ibld]]

            # connect cables to blades tops - connect to master blade, same one that the other blades connect to
            jointconn[ibld+nblade*2,:] = [b_topidx[1] c2b_topidx[ibld]]

            #Connect first Blade to all other blades
            if ibld>1
                jointconn[ibld-1+nblade*3,:] = [b_topidx[1] b_topidx[ibld]]
            end
        end
    else
        jointconn = zeros(Int,nblade*2-1+ncable,2)
        for ibld = 1:nblade
            # connect tower to blades
            jointconn[ibld,:] = [t_topidx b_botidx[ibld]]

            # # connect cables to blades bases
            # jointconn[ibld+nblade,:] = [b2c_botidx[ibld] c2b_botidx[ibld]]

            # connect cables to blades tops - connect to master blade, same one that the other blades connect to
            jointconn[ibld+nblade,:] = [b_topidx[1] c2b_topidx[ibld]]

            #Connect first Blade to all other blades
            if ibld>1
                jointconn[ibld-1+nblade*2,:] = [b_topidx[1] b_topidx[ibld]]
            end
        end
    end



    njoint = length(jointconn[:,1])
    ort = calculateElementOrientation(mymesh)
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
        Psi_d_joint[jnt] = ort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = ort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

    #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint).+joint_type zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]


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

    return mymesh, ort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng
end

function calculateElementOrientation2(mesh)

    # Note on gimbal lock:
    #   when calculating a (roll -> pitch -> yaw) rotation sequence (twist -> theta -> psi) on a vertical element, it is ambiguous
    #   if the roll or yaw (twist or psi) should be set.  It is therefore up to the designer to pick which is used.  For gimbal lock
    #   the yaw (psi) is currently defined as zero by choice and the roll (twist) set to non-zero.  It could be equivalently
    #   calculated as the roll (twist) is zero, the yaw (psi) is non-zero.  The latter scenario may actually be simpler for
    #   calculating DCM's when coupling to other codes.

    numEl = mesh.numEl #get number of elements
    Psi_d=zeros(numEl) #initialize Psi, Theta, Twist, and Offset Arrays
    Theta_d=zeros(numEl)
    twist_d=zeros(numEl)
    twist_d2=zeros(numEl)
    Offset=zeros(3,numEl)    #offset is the hub frame coordinate of node 1 of the element
    vsave=zeros(numEl,3)    #offset is the hub frame coordinate of node 1 of the element
    elNum=zeros(numEl,2) #initialize element number array


    #calculate "mesh centroid"
    meshCentroid = [0.0 0.0 Statistics.mean(mesh.z)] #calculate a geometric centroid using all nodal coordinates
    lenv = zeros(numEl)
    for i = 1:numEl #loop over elements

        n1 = Int(mesh.conn[i,1]) #n1 := node number for node 1 of element i
        n2 = Int(mesh.conn[i,2]) #n2 := node number for node 2 of element i

        p1 = [mesh.x[n1] mesh.y[n1] mesh.z[n1]] #nodal coordinates of n1
        p2 = [mesh.x[n2] mesh.y[n2] mesh.z[n2]] #nodal coordinates of n2
        Offset[:,i] = p1 #set offset as position of n1

        v=p2-p1 #define vector from p1 to p2
        vsave[i,:] = v./LinearAlgebra.norm(v)
        lenv[i] = LinearAlgebra.norm(v) #calculate element lengtt

        Psi_d[i],Theta_d[i] = calculatePsiTheta(v) #calculate elment Psi and Theta angles for orientation
        elNum[i,:] = mesh.conn[i,:] #get node number map

        if isapprox(v[1],0.0;atol=1e-5) && isapprox(v[2],0.0;atol=1e-5) #If we are vertical, correct for gimbal lock
            # Calculate the twist for the blades such that the normal vector points away from the turbine center of rotation.
            # get vector from center of rotation to element.
            if mesh.type[i]==0 || mesh.type[i]==2 #Mesh Type: 0-blade 1-tower 2-strut
                v2 = [0 0 0]-(p1+p2)/2
                twist_d[i], _ = calculatePsiTheta(v2)
                if v[1]*cosd(Psi_d[i])+v[2]*sind(Psi_d[i]) < 0.0 # if normal vector is pointing towards center of turbine
                    twist_d[i] += 180
                    twist_d[i] = twist_d[i]%360.0 #correct if greater than 360 degrees
                end
            elseif mesh.type[i]==1
                twist_d[i] = 0.0        # This assumes a radially symmetric tower such that twist does not matter
            # elseif mesh.type[i]==2
            #     twist_d[i] = 0.0        # This assumption implies that the strut is symmetric with no pitch angle
            end
        else #otherwise, proceed as before
            nVec1,nVec2,nVec3 = rigidBodyRotation(0,0,1,[Psi_d[i],Theta_d[i]],[3,2]) #tranform a local normal "flapwise" vector in the element frame to the hub frame
            nVec = [nVec1 nVec2 nVec3]

            # Mesh Type: 0-blade 1-tower 2-strut
            if mesh.type[i]==2
                refVector = meshCentroid-p1
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

    end

    # for i = 1:length(twist_d)
    #     println("$i $(twist_d[i]) $(twist_d2[i]) vec $(vsave[i,1]) $(vsave[i,2]) $(vsave[i,3]) psi $(Psi_d[i]) theta $(Theta_d[i])")#-Gimbal[i])")
    # end


    #assign data to element orientation (Ort) object
    return OWENSFEA.Ort(Psi_d,Theta_d,twist_d,lenv,elNum,Offset)
end

"""

    calculateElementOrientation(mesh)

Calculates the orientation of elements in a mesh.

#Input
* `mesh::OWENSFEA.Mesh`  see ?OWENSFEA.Mesh object containing mesh data

#Output
* `elOr::OWENSFEA.Ort`  see ?OWENSFEA.Ort object containing element orientation data
"""
function calculateElementOrientation(mesh)

    # Note on gimbal lock:
    #   when calculating a (roll -> pitch -> yaw) rotation sequence (twist -> theta -> psi) on a vertical element, it is ambiguous
    #   if the roll or yaw (twist or psi) should be set.  It is therefore up to the designer to pick which is used.  For gimbal lock
    #   the yaw (psi) is currently defined as zero by choice and the roll (twist) set to non-zero.  It could be equivalently
    #   calculated as the roll (twist) is zero, the yaw (psi) is non-zero.  The latter scenario may actually be simpler for
    #   calculating DCM's when coupling to other codes.

    numEl = mesh.numEl #get number of elements
    Psi_d=zeros(numEl) #initialize Psi, Theta, Twist, and Offset Arrays
    Theta_d=zeros(numEl)
    twist_d=zeros(numEl)
    Offset=zeros(3,numEl)    #offset is the hub frame coordinate of node 1 of the element
    elNum=zeros(numEl,2) #initialize element number array

    lenv = zeros(numEl)
    for i = 1:numEl #loop over elements

        n1 = Int(mesh.conn[i,1]) #n1 := node number for node 1 of element i
        n2 = Int(mesh.conn[i,2]) #n2 := node number for node 2 of element i

        p1 = [mesh.x[n1] mesh.y[n1] mesh.z[n1]] #nodal coordinates of n1
        p2 = [mesh.x[n2] mesh.y[n2] mesh.z[n2]] #nodal coordinates of n2
        Offset[:,i] = p1 #set offset as position of n1
        elpos = (p1.+p2)./2

        Psi_d[i] = -atand(elpos[1],elpos[2]).-90.0 #global yaw position, this calculation is agnostic to vertical elements, keep in mind that top dead center is 0 degrees yaw

        if mesh.type[i] == 4 # treat tangentially aligned mesh components different than radially aligned
            Psi_d[i] -= 90.0
            twist_d[i] += 90.0
        end

        # Now with the global yaw position know, get the node points in a consistent frame of reference to calculate the delta, or the slope of the element
        p1[1],p1[2],p1[3] = rigidBodyRotation(p1[1],p1[2],p1[3],[-Psi_d[i]],[3])
        p2[1],p2[2],p2[3] = rigidBodyRotation(p2[1],p2[2],p2[3],[-Psi_d[i]],[3])
        
        v=p2-p1
        v[abs.(v).<1e-7] .= 0.0 #zero out close to zero differences

        Theta_d[i]  = atand(v[1],v[3]).-90.0

        lenv[i] = LinearAlgebra.norm(v) #calculate element length

        elNum[i,:] = mesh.conn[i,:] #get node number map

    end

    #assign data to element orientation (Ort) object
    # yaw, blade slope, angle of attack
    return OWENSFEA.Ort(Psi_d,Theta_d,twist_d,lenv,elNum,Offset)
end

"""

    createGeneralTransformationMatrix(angleArray,axisArray)

Calculates the transformation matrix assocaited with a general Euler rotation sequence.

#Input
* `angleArray`:      = array of angles for Euler rotation sequence
* `axisArray`:       = array of axis of rotatoins for Euler rotation

#Output
* `dcmTotal`:        = transformation matrix of specified euler rotation sequence
"""
function createGeneralTransformationMatrix(angleArray,axisArray)

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

"""
Creates a direction cosine matrix (dcm) associated with a rotation of angleDeg about axisNum.
"""
function createSingleRotationDCM(angleDeg,axisNum)

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

"""

    rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)

Performs a coordinate transformation from a local
body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
rotation sequence

#Input
* `B1`:         array containing body frame 1 coordinates of points to be mapped to the hub frame
* `B2`:         array containing body frame 2 coordinates of points to be mapped to the hub frame
* `B3`:         array containing body frame 3 coordinates of points to be mapped to the hub frame
* `AngleArray`:  Array of angles for Euler rotation sequence
* `AxisArray`:   Array of axes for Euler rotation sequence

#Output
* `H1`:         array containg hub frame 1 coordinates of points mapped to the hub frame from body frame
* `H2`:         array containg hub frame 2 coordinates of points mapped to the hub frame from body frame
* `H3`:         array containg hub frame 3 coordinates of points mapped to the hub frame from body frame

That is CHtoB = [M3(SweepAngle)][M2(Theta)][M3(Psi)];
"""
function rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)

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

"""

    calculatePsiTheta(v)

Calculates the orientation of a single element. A local
element frame is related to a hub frame through a transformation matrix
CHtoE (transforming a vector from an element frame E to a global frame
H) such that CHtoE = [M2(Theta)]*[M3(Psi)]. Here [M2( )] is a direction
cosine matrix about a 2 axis and [M3( )] is a direction cosine matrix
about a 3 axis.

#Input
* `v`:      vector from node 1 to node 2 of an element

#Output
* `Psi`:    "3" angle for element orientation (deg)
* `Theta`:  "2" angle for element orientation (deg)
"""
function calculatePsiTheta(v)

    v = v./LinearAlgebra.norm(v) #normalize vector by its length #NOTE: eps() offset is for automatic gradient propogation
    # avoid gimbal lock
    if isapprox(v[2], 0.0; atol=eps(), rtol=0) && isapprox(v[1], 0.0; atol=eps(), rtol=0)
        Psi_d = 0.0     # Define as zero due to gimbal lock state (we can instead set twist)
    else
        Psi_d = atan(v[2],v[1])*180.0/pi #calculate sweep angle, convert to deg
    end
    Theta_d = -asin(v[3]-eps())*180.0/pi # atan(v[3],v[1])*180.0/pi #calculate theta angle, convert to deg

    return Psi_d,Theta_d
end

"""

    getOWENSPreCompOutput(numadIn;yscale=1.0,plyprops = plyproperties())

Takes numad formatted inputs for composite layup and material properties and runs them through OWENSPreComp

#Inputs
* `numadIn::NuMad`: see ?NuMad inputs
* `yscale::float`: airfoil thickness scaling
* `plyprops::plyproperties()`: see ?plyproperties for input material properties

#Outputs
* `precompoutput::OWENSPreComp.Output`: see ?OWENSPreComp.Input
* `precompinput::OWENSPreComp.Input`: see ?OWENSPreComp.properties
* `lam_U::Composites.Laminate`: laminate stacks used for post processing, size (n_stations, n_upper chorwise stations), upper surface, see ?Composites.Laminate
* `lam_L::Composites.Laminate`: laminate stacks used for post processing, size (n_stations, n_lower chorwise stations), lower surface, see ?Composites.Laminate
* `lam_W::Composites.Laminate`: laminate stacks used for post processing, size (n_stations, n_webs), shear webs, see ?Composites.Laminate
"""
function getOWENSPreCompOutput(numadIn;yscale=1.0,plyprops = plyproperties())
    # Assume each station has the same number of segments and webs (some of which may be of zero width/thickness)

    n_stations = length(numadIn.span)

    if length(yscale) == 1
        yscale = ones(n_stations).*yscale
    end

    precompinput = Array{OWENSPreComp.Input,1}(undef,n_stations)
    precompoutput = Array{OWENSPreComp.Output,1}(undef,n_stations)
    lam_U = Array{Composites.Laminate,2}(undef,n_stations,sum((numadIn.segments[1,:].>0.0)[2:end]))
    lam_L = Array{Composites.Laminate,2}(undef,n_stations,sum((numadIn.segments[1,:].<=0.0)[2:end]))
    lam_W = Array{Composites.Laminate,2}(undef,n_stations,numadIn.n_web)

    normalchord = numadIn.chord #TODO: if any sweep occurs the normal chord will need to be calculated since it isn't the regular chord anymore

    sloc = numadIn.span
    twist_d = zeros(length(numadIn.twist_d))
    for ii = 1:length(numadIn.twist_d) #TODO: real interpolation in the reading file
        twist_d[ii] = numadIn.twist_d[ii]
    end

    twistrate_d = OWENSPreComp.tw_rate(n_stations,sloc[1:n_stations],twist_d)
    leloc = numadIn.xoffset

    for i_station = 1:n_stations

        ######################################
        ######## Airfoil Shape Input #########
        ######################################

        #TODO: ensure square matrix via interpolation so it can be stored as a 2D matrix
        af_xy = DelimitedFiles.readdlm("$(numadIn.airfoil[i_station]).csv",',',Float64,skipstart = 0)

        # Normalize the surface points and make sure that they start at the trailing edge and loop around starting on the bottom side #TODO: add a check for both of these
        #TODO: simplify this since there is circshift happening on these points later on
        if af_xy[2,2]>af_xy[end-1,2]
            xaf = reverse(af_xy[:,1])./maximum(af_xy[:,1])
            yaf = reverse(af_xy[:,2])./maximum(af_xy[:,1]).*yscale[i_station]
        else
            xaf = af_xy[:,1]./maximum(af_xy[:,1])
            yaf = af_xy[:,2]./maximum(af_xy[:,1]).*yscale[i_station]
        end

        # find leading edge
        lei = argmin(abs.(xaf))
        # shift so leading edge is first
        xpc = circshift(xaf,-(lei-1))
        ypc = circshift(yaf,-(lei-1))

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

        n_pliesU = zeros(sum(n_laminaU))
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

        n_pliesL = zeros(sum(n_laminaL))
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

        n_pliesW = zeros(sum(n_laminaW))
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
        precompinput[i_station] = OWENSPreComp.Input(
        normalchord[i_station],
        -twist_d[i_station],-twistrate_d[i_station],
        leloc[i_station],xpc_filtered,ypc_filtered,
        e1,e2,g12,anu12,density,
        xsec_nodeU,n_laminaU,n_pliesU,t_lamU,tht_lamU,mat_lamU,
        xsec_nodeL,n_laminaL,n_pliesL,t_lamL,tht_lamL,mat_lamL,
        loc_web,n_laminaW,n_pliesW,t_lamW,tht_lamW,mat_lamW)

        # calculate composite properties: stiffness, mass, etc
        precompoutput[i_station] = OWENSPreComp.properties(precompinput[i_station])
        # Store all lamina at each segment for each station

        # Compensate for non-integer number of plies by adjusting thickness and setting number of plies to int
        t_lamU = (t_lamU.*n_pliesU)./floor.(n_pliesU)
        t_lamL = (t_lamL.*n_pliesL)./floor.(n_pliesL)
        t_lamW = (t_lamW.*n_pliesW)./floor.(n_pliesW)
        
        # Handle case when floor goes to 0
        t_lamU[isinf.(t_lamU)] .= eps()
        t_lamL[isinf.(t_lamL)] .= eps()
        t_lamW[isinf.(t_lamW)] .= eps()

        start_idx = 1
        for i_lam = 1:length(n_laminaU)
            end_idx = sum(n_laminaU[1:i_lam])
            lam_U[i_station,i_lam] = Composites.Laminate(mat_lamU[start_idx:end_idx],floor.(Int,n_pliesU[start_idx:end_idx]),t_lamU[start_idx:end_idx],tht_lamU[start_idx:end_idx])
            start_idx = end_idx+1
        end
        start_idx = 1
        for i_lam = 1:length(n_laminaL)
            end_idx = sum(n_laminaL[1:i_lam])
            lam_L[i_station,i_lam] = Composites.Laminate(mat_lamL[start_idx:end_idx],floor.(Int,n_pliesL[start_idx:end_idx]),t_lamL[start_idx:end_idx],tht_lamL[start_idx:end_idx])
            start_idx = end_idx+1
        end
        start_idx = 1
        for i_lam = 1:length(n_laminaW)
            end_idx = sum(n_laminaW[1:i_lam])
            lam_W[i_station,i_lam] = Composites.Laminate(mat_lamW[start_idx:end_idx],floor.(Int,n_pliesW[start_idx:end_idx]),t_lamW[start_idx:end_idx],tht_lamW[start_idx:end_idx])
            start_idx = end_idx+1
        end
    end

    return precompoutput,precompinput,lam_U,lam_L,lam_W
end

"""

    getSectPropsFromOWENSPreComp(usedUnitSpan,numadIn,precompoutput;GX=false,precompinputs=nothing)

Arranges the precomp output into the sectional properties required by OWENSFEA

#Inputs
* `usedUnitSpan::Array{<:float}`: Array specifying the relative (0-1) mesh z locations where the sectional properties are actually called for (interpolation)
* `numadIn::NuMAD`: see ?NuMAD
* `precompoutput::OWENSPreComp.properties`: see ?OWENSPreComp.properties
* `GX::bool`: optional specifies if GX outputs should be output
* `precompinputs`: optional

#Outputs
* `sectionPropsArray::SectionPropsArray`: see ?OWENSFEA.SectionPropsArray, if !GX bool
* `stiff`: if GX bool
* `mass`: if GX bool
stiff, mass

"""
function getSectPropsFromOWENSPreComp(usedUnitSpan,numadIn,precompoutput;GX=false,precompinputs=nothing)
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
        gj[i_pc] = precompoutput[i_pc].gj.*5.1
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
    ei_flap_used = safeakima(origUnitSpan,ei_flap,usedUnitSpan)
    ei_lag_used = safeakima(origUnitSpan,ei_lag,usedUnitSpan)
    gj_used = safeakima(origUnitSpan,gj,usedUnitSpan)
    ea_used = safeakima(origUnitSpan,ea,usedUnitSpan)
    s_fl_used = safeakima(origUnitSpan,s_fl,usedUnitSpan)
    s_af_used = safeakima(origUnitSpan,s_af,usedUnitSpan)
    s_al_used = safeakima(origUnitSpan,s_al,usedUnitSpan)
    s_ft_used = safeakima(origUnitSpan,s_ft,usedUnitSpan)
    s_lt_used = safeakima(origUnitSpan,s_lt,usedUnitSpan)
    s_at_used = safeakima(origUnitSpan,s_at,usedUnitSpan)
    x_sc_used = safeakima(origUnitSpan,x_sc,usedUnitSpan)
    y_sc_used = safeakima(origUnitSpan,y_sc,usedUnitSpan)
    x_tc_used = safeakima(origUnitSpan,x_tc,usedUnitSpan)
    y_tc_used = safeakima(origUnitSpan,y_tc,usedUnitSpan)
    mass_used = safeakima(origUnitSpan,mass,usedUnitSpan)
    flap_iner_used = safeakima(origUnitSpan,flap_iner,usedUnitSpan)
    lag_iner_used = safeakima(origUnitSpan,lag_iner,usedUnitSpan)
    tw_iner_d_used = safeakima(origUnitSpan,tw_iner_d,usedUnitSpan)
    x_cm_used = safeakima(origUnitSpan,x_cm,usedUnitSpan).*0.0
    y_cm_used = safeakima(origUnitSpan,y_cm,usedUnitSpan).*0.0

    ac_used = safeakima(origUnitSpan,numadIn.aerocenter,usedUnitSpan)
    twist_d_used = safeakima(origUnitSpan,numadIn.twist_d,usedUnitSpan)
    chord_used = safeakima(origUnitSpan,numadIn.chord,usedUnitSpan)

    if GX

        stiff = Array{Array{Float64,2}, 1}(undef, length(usedUnitSpan)-1)
        mass = Array{Array{Float64,2}, 1}(undef, length(usedUnitSpan)-1)

        for i=1:length(usedUnitSpan)-1
            GA = ea_used[i]/2.6*5/6
            stiff[i]= [ea_used[i] 0.0 0.0 s_at_used[i] s_af_used[i] s_al_used[i]
                            0.0 GA 0.0 0.0 0.0 0.0
                            0.0 0.0 GA 0.0 0.0 0.0
                            s_at_used[i] 0.0 0.0 gj_used[i] s_ft_used[i] s_lt_used[i]
                            s_af_used[i] 0.0 0.0 s_ft_used[i] ei_flap_used[i] s_fl_used[i]
                            s_al_used[i] 0.0 0.0 s_lt_used[i] s_fl_used[i] ei_lag_used[i]]

            ux3 = (mass_used[i]*y_cm_used[i])
            ux2 = (mass_used[i]*x_cm_used[i])
            mass[i] = [mass_used[i] 0.0 0.0 0.0 ux3 -ux2
                          0.0 mass_used[i] 0.0 -ux3 0.0 0.0
                          0.0 0.0 mass_used[i] ux2 0.0 0.0
                          0.0 -ux3 ux2 flap_iner_used[i]+lag_iner_used[i] 0.0 0.0
                          ux3 0.0 0.0 0.0 flap_iner_used[i] -tw_iner_d_used[i]
                          -ux2 0.0 0.0 0.0 -tw_iner_d_used[i] lag_iner_used[i]]
        end

        return stiff, mass
    end


    sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, length(usedUnitSpan)-1)

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
        rhoJ = [flap_iner_used[i]+lag_iner_used[i], flap_iner_used[i+1]+lag_iner_used[i+1]]
        zcm = [0.0,0.0]# [x_cm_used[i], x_cm_used[i+1]]
        ycm = [0.0,0.0]# [y_cm_used[i], y_cm_used[i+1]]
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
        sectionPropsArray[i] = OWENSFEA.SectionPropsArray(ac,twist_d,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

    end

    if !isnothing(precompinputs)
        # Airfoil data
        nspl = 100
        # Spline the airfoil data to a common size
        myxafpc = zeros(length(precompinputs),nspl*2+1)
        myyafpc = zeros(length(precompinputs),nspl*2+1)
        myzafpc = numadIn.span./maximum(numadIn.span)
        mychord = zeros(length(precompinputs))
        for ipci = 1:length(precompinputs)
            mychord[ipci] = precompinputs[ipci].chord
            xaf = precompinputs[ipci].xnode .* mychord[ipci]
            yaf = precompinputs[ipci].ynode .* mychord[ipci]

            LE_pt = findfirst(x->x==mychord[ipci],xaf)

            xaf_top = xaf[1:LE_pt]
            yaf_top = yaf[1:LE_pt]

            xaf_bot = xaf[LE_pt:end]
            yaf_bot = yaf[LE_pt:end]

            myxpts_top = LinRange(xaf_top[1],xaf_top[end],nspl)
            myxpts_bot = LinRange(xaf_bot[1],xaf_bot[end],nspl)

            myypts_top = safeakima(xaf_top,yaf_top,myxpts_top)
            myypts_bot = safeakima(reverse(xaf_bot),reverse(yaf_bot),reverse(myxpts_bot))

            myxafpc[ipci,:] = [myxpts_top;myxpts_bot;myxpts_top[1]]
            myyafpc[ipci,:] = [myypts_top;reverse(myypts_bot);myypts_top[1]]
            # PyPlot.figure()
            # PyPlot.plot(xaf,yaf,"k",label="orig")
            # PyPlot.plot(xaf_top,yaf_top,"r",label="top")
            # PyPlot.plot(xaf_bot,yaf_bot,"b",label="bot")
            # PyPlot.plot(myxpts_top,myypts_top,"r.",label="mytop")
            # PyPlot.plot(myxpts_bot,myypts_bot,"b.",label="mybot")
            # PyPlot.plot(myxaf,myyaf,"k+-",label="myaf")
        end

        # Spline the airfoil data to align with the mesh
        myxaf = zeros(length(usedUnitSpan)-1,nspl*2+1)
        myyaf = zeros(length(usedUnitSpan)-1,nspl*2+1)
        myzaf = cumsum(diff(usedUnitSpan))
        myzaf = myzaf .-= myzaf[1]
        # chord = safeakima(myzafpc,mychord,myzaf)

        for iline = 1:nspl*2+1
            myxaf[:,iline] = safeakima(myzafpc,myxafpc[:,iline],myzaf)
            myyaf[:,iline] = safeakima(myzafpc,myyafpc[:,iline],myzaf)
        end
    else
        myxaf = nothing
        myyaf = nothing
    end

    for i=1:length(usedUnitSpan)-1

            sectionPropsArray[i].b = 0.5.*[chord_used[i], chord_used[i+1]] #element semi chord
            # sectionPropsArray[i].a0 = [bladeData[i,12], bladeData[i+1,12]]         #element lift curve slope (needed for flutter analysis) TODO: enable coupling between actual airfoil lift slope

            #convert "a" to semichord fraction aft of halfchord
            sectionPropsArray[i].a = (sectionPropsArray[i].a .+ 0.25*2*sectionPropsArray[i].b .- sectionPropsArray[i].b)./sectionPropsArray[i].b

            #convert "ac" to semichord fraction foreward of halfchord TODO: why are we doing it this way???
            sectionPropsArray[i].ac = sectionPropsArray[i].ac.*2

            #physical aero center offset from elastic axis
            sectionPropsArray[i].aeroCenterOffset = sectionPropsArray[i].ac .* sectionPropsArray[i].b .- sectionPropsArray[i].a

            #airfoil coordinates, if supplied
            if !isnothing(precompinputs)
                sectionPropsArray[i].xaf = myxaf[i,:]
                sectionPropsArray[i].yaf = myyaf[i,:]
            end
    end

    # println("EIyz, rhoIyz deactivated") #TODO: why is this, especially when I believe precomp calculates them
    return sectionPropsArray

end


"""
create_hawt_mesh(;Ht = 15.0,
    Hb = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    ncelem = 4,
    strut_mountpoint = 0.01, # This puts struts at top and bottom
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    joint_type = 0,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0)


ARCUS mesh configuration: no tower between blades, no struts, but cables from top center attaching to specified blade mount point at base

#Inputs
* `Ht::float`: height of tower before blades attach (m)
* `Hb::float`: blade height (m)
* `R::float`: bade radius (m)
* `nblade::int`: number of blades
* `ntelem::int`: number of tower elements
* `nbelem::int`: number of blade elements
* `ncelem::int`: number of strut elements
* `c_mount_ratio::float`: factor of blade height where the struts attach on both top and bottom
* `bshapex::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapez::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `joint_type::int`: 0 is fixed, 1 is about x-axis, 2 is y-axis, etc
* `cables_connected_to_blade_base::bool`: = true,
* `angularOffset::float`: (rad) angular offset of mesh generation, typically used to match CACTUS input.  Value of 0 puts blade 1 at the "north" position and the others populate counterclockwise when looking down

#Outputs
* `mymesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `ort::OWENSFEA.Ort`: see ?OWENSFEA.Ort
* `myjoint:Array{<:float}`: see ?OWENSFEA.FEAModel.joint

"""
function create_hawt_mesh(;
    hub_depth = 15.0, #Hub Beam Depth
    tip_precone = 1.0, #blade precone
    R = 54.014, # m bade radius
    AD15hubR = 7.0,
    nblade = 3,
    ntelem = 4, #tower elements
    nbelem = 30, #blade elements
    ncelem = 10,  #cable elements
    c_mount_ratio = 0.05, #fraction of blade where the cable mounts
    bshapex = LinRange(0.0,1.0,nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    joint_type = 0,
    angularOffset = 0.0,
    AD15_ccw = false)

    ##################################
    #            
    #           |   _|_
    #          / \   |
    #           
    #           
    # 
    ####################################
    ###------------Tower--------------##
    ####################################
    mesh_z = collect(LinRange(0,hub_depth,ntelem+1))
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
    bld_Y = collect(LinRange(0.0,R,nbelem+1))

    bshapex = R .* bshapex./maximum(bshapex)
    if maximum(bshapez) != 0.0 #TODO more robust
        bshapez = tip_precone .* bshapez./maximum(bshapez)
    end
    bld_Z = safeakima(bshapex,bshapez,bld_Y)

    bld_X = zero(bld_Z)

    bld_Z .+= hub_depth

    b_Z = []
    b_X = []
    b_Y = []
    AD15bldNdIdxRng = zeros(Int64,0,2)
    # Now using standard VAWT convention, blade 1 is zero degrees at top dead center, or North/Y+
    # and they are offset counter clockwise
    b_topidx = zeros(Int,nblade)
    b_botidx = zeros(Int,nblade) .+ length(mesh_z)
    conn_b = zeros(length(bld_Z)-1,2)
    for ibld = 1:nblade
        myangle = (ibld-1)*2.0*pi/nblade + angularOffset

        # find the closest point to hub radius and move it to exactly match
        b_xstart = bld_X[1] #mesh_x[b_startidx]
        b_ystart = bld_Y[1] #mesh_y[b_startidx]
        b_zstart = bld_Z[1] #mesh_z[b_startidx]

        b_xend = bld_X[end] #mesh_x[b_endidx]
        b_yend = bld_Y[end] #mesh_y[b_endidx]
        b_zend = bld_Z[end] #mesh_z[b_endidx]
        hubIdx=1
        if AD15hubR > 1e-6
            lenXY = sqrt((b_xend - b_xstart)^2 + (b_yend - b_ystart)^2)   # strut length in XY
            minR2 = lenXY 
            for i = 1:nbelem+1  # step through to find closest point to hub radius on x-y plane
                R2 = AD15hubR - sqrt((bld_X[i] - b_xstart)^2 + (bld_Y[i] - b_ystart)^2)
                if abs(R2) < abs(minR2)
                    hubIdx = i
                    minR2 = R2
                end
            end
            R_temp = minR2
            println("Hub crossing at idx $hubIdx at $R_temp with hub radius of $AD15hubR")
            print("Moving strut point from [$(bld_X[hubIdx]),$(bld_Y[hubIdx]),$(bld_Z[hubIdx])] to ")
            bld_X[hubIdx] = bld_X[hubIdx] + R_temp/lenXY*(b_xend-b_xstart)
            bld_Y[hubIdx] = bld_Y[hubIdx] + R_temp/lenXY*(b_yend-b_ystart)
            bld_Z[hubIdx] = bld_Z[hubIdx] + R_temp/lenXY*(b_zend-b_zstart)
            print("[$(bld_X[hubIdx]),$(bld_Y[hubIdx]),$(bld_Z[hubIdx])]\n")
        end
        # set index of where this point is in the mesh
        

        b_Z = [b_Z;bld_Z]
        b_X = [b_X;bld_Y.*sin(myangle).+bld_X.*cos(myangle)]
        b_Y = [b_Y;bld_X.*sin(myangle).+bld_Y.*cos(myangle)]

        hubIdx = length(mesh_z) + length(bld_Z)*(ibld-1) + hubIdx

        # Element joint indices
        b_botidx[ibld] = length(mesh_z)+1 + length(bld_Z)*(ibld-1)
        b_topidx[ibld] = length(mesh_z)+1 + length(bld_Z)*ibld-1

        # Intraconnectivity
        conn_b[:,1] = collect(b_botidx[ibld]:1:b_topidx[ibld]-1)
        conn_b[:,2] = collect(b_botidx[ibld]+1:1:b_topidx[ibld])
        conn = [conn;conn_b]

        if AD15_ccw
            AD15bldNdIdxRng = [AD15bldNdIdxRng; hubIdx b_topidx[ibld]]    # top of blade is root
        else
            AD15bldNdIdxRng = [AD15bldNdIdxRng; b_topidx[ibld] hubIdx]    # bottom of blade is root
        end
    end

    # Add to the mesh
    mesh_z = [mesh_z;b_Z]
    mesh_x = [mesh_x;b_X]
    mesh_y = [mesh_y;b_Y]

    #######################################
    ###   Cleanup/Derived parameters    ###
    #######################################

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:t_topidx] .= 1 #Tower
    # meshtype[idx_bot_lbld_tower:idx_top_rbld_tower] .= 0 #Blades

    #########################
    # .bld equivalent
    #########################

    # For a single blade
    meshSeg = zeros(1+nblade) #tower, blades, and cables

    meshSeg[1] = ntelem
    meshSeg[2:nblade+1] .= nbelem

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z))
    structuralNodeNumbers = zeros(nblade,length(bld_Z))
    structuralElNumbers = zeros(nblade,length(bld_Z))

    for iblade = 1:nblade

        # Normalized Span
        span_len = sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-hub_depth).^2.0)
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

    jointconn = zeros(Int,nblade,2)
    for ibld = 1:nblade
        # connect tower to blades
        jointconn[ibld,:] = [t_topidx b_botidx[ibld]]
    end

    njoint = length(jointconn[:,1])
    ort = OWENS.calculateElementOrientation(mymesh)
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
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint).+joint_type zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    AD15bldElIdxRng = zeros(Int64,0,2)
    for i = 1:size(AD15bldNdIdxRng,1)
        if AD15bldNdIdxRng[i,2] > AD15bldNdIdxRng[i,1]  # ascending order
            idx1 = findfirst(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,1])
            idx2 = findfirst(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
        else    # upside down oriented blade
            idx1 = findlast(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,2])
            idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,1])
        end
        AD15bldElIdxRng = [AD15bldElIdxRng; idx1 idx2]
        #println("Idx: $(AD15bldNdIdxRng[i,:])    idx: $([idx1 idx2])")
    end
    return mymesh, ort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng
end


"""
create_hawt_biwing_mesh(;Ht = 15.0,
    Hb = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    ncelem = 4,
    strut_mountpoint = 0.01, # This puts struts at top and bottom
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    joint_type = 0,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0)

#Inputs
* `Ht::float`: height of tower before blades attach (m)
* `Hb::float`: blade height (m)
* `R::float`: bade radius (m)
* `nblade::int`: number of blades
* `ntelem::int`: number of tower elements
* `nbelem::int`: number of blade elements
* `ncelem::int`: number of strut elements
* `c_mount_ratio::float`: factor of blade height where the struts attach on both top and bottom
* `bshapex::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapez::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `joint_type::int`: 0 is fixed, 1 is about x-axis, 2 is y-axis, etc
* `cables_connected_to_blade_base::bool`: = true,
* `angularOffset::float`: (rad) angular offset of mesh generation, typically used to match CACTUS input.  Value of 0 puts blade 1 at the "north" position and the others populate counterclockwise when looking down

#Outputs
* `mymesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `ort::OWENSFEA.Ort`: see ?OWENSFEA.Ort
* `myjoint:Array{<:float}`: see ?OWENSFEA.FEAModel.joint

"""
function create_hawt_biwing_mesh(;
    nblade = 3,
    hub_depth = 15.0, #Hub Beam Depth
    R_root = 10.0, # m biwing radius
    R_biwing = 30.0, # outer radius
    R_tip = 54.014, # outer radius
    AD15hubR = 7.0,
    ntelem = 4, #tower elements
    nbelem_root = 30, #biwing elements, for each 
    nbelem_biwing = 30, #tip elements
    nbelem_tip = 30, #tip elements
    bshapex_root = LinRange(0.0,R_root,nbelem_root+1), #Blade shape, magnitude is relevant
    bshapez_root = zeros(nbelem_root+1), #Blade shape, magnitude is relevant
    bshapex_biwing_U = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_U = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_biwing_L = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_L = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_tip = LinRange(R_biwing,R_tip,nbelem_tip+1), #Blade shape, magnitude is relevant
    bshapez_tip = zeros(nbelem_tip+1), #Blade shape, magnitude is relevant
    joint_type = 0,
    angularOffset = 0.0,
    AD15_ccw = false)

    ##################################
    #            
    #         -<>-o-<>-
    #           
    ####################################
    ###------------Tower--------------##
    ####################################
    mesh_z = collect(LinRange(0,hub_depth,ntelem+1))
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

    # Root
    bld_root_Y = collect(LinRange(0.0,R_root,nbelem_root+1))
    bld_root_Z = safeakima(bshapex_root,bshapez_root,bld_root_Y)    
    bld_root_X = zero(bld_root_Y)

    # Biwing upper
    bld_biwing_U_Y = collect(LinRange(R_root,R_biwing,nbelem_biwing+1))
    bld_biwing_U_Z = safeakima(bshapex_biwing_U,bshapez_biwing_U,bld_biwing_U_Y)
    bld_biwing_U_X = zero(bld_biwing_U_Y)

    # Biwing lower
    bld_biwing_L_Y = collect(LinRange(R_root,R_biwing,nbelem_biwing+1))
    bld_biwing_L_Z = safeakima(bshapex_biwing_L,bshapez_biwing_L,bld_biwing_L_Y)
    bld_biwing_L_X = zero(bld_biwing_L_Y)

    # Tip
    bld_tip_Y = collect(LinRange(R_biwing,R_tip,nbelem_tip+1))
    bld_tip_Z = safeakima(bshapex_tip,bshapez_tip,bld_tip_Y)
    bld_tip_X = zero(bld_tip_Y)

    # Iterate around the hub

    bldsec_Z = [bld_root_Z,bld_biwing_U_Z,bld_biwing_L_Z,bld_tip_Z]
    bld_Z = [bld_root_Z;bld_biwing_U_Z;bld_biwing_L_Z;bld_tip_Z] .+ hub_depth
    N_sec = length(bldsec_Z)

    bldsec_X = [bld_root_X,bld_biwing_U_X,bld_biwing_L_X,bld_tip_X]
    bld_X = [bld_root_X;bld_biwing_U_X;bld_biwing_L_X;bld_tip_X]

    bldsec_Y = [bld_root_Y,bld_biwing_U_Y,bld_biwing_L_Y,bld_tip_Y]
    bld_Y = [bld_root_Y;bld_biwing_U_Y;bld_biwing_L_Y;bld_tip_Y]

    AD15bldNdIdxRng = zeros(Int64,0,2)

    # Now using standard VAWT convention, blade 1 is zero degrees at top dead center, or North/Y+
    # and they are offset counter clockwise
    b_sectopidx = zeros(Int,nblade,length(bldsec_Z))
    b_secbotidx = zeros(Int,nblade,length(bldsec_Z)) .+ length(mesh_z)
    startlenmeshz = length(mesh_z)
    for ibld = 1:nblade
        myangle = (ibld-1)*2.0*pi/nblade + angularOffset

        if ibld != 1
            startlenmeshz = b_sectopidx[ibld-1,end]
        end
        for isec = 1:N_sec
            mesh_z = [mesh_z;bldsec_Z[isec] .+ hub_depth]
            mesh_x = [mesh_x;bldsec_Y[isec].*sin(myangle).+bldsec_X[isec].*cos(myangle)]
            mesh_y = [mesh_y;bldsec_X[isec].*sin(myangle).+bldsec_Y[isec].*cos(myangle)]

            # Element joint indices #TODO: here is where the error is, was originally length of the local blade section, not the whole mesh
            if isec == 1
                b_secbotidx[ibld,isec] = startlenmeshz+1 #+ length(bldsec_Z[isec])*(ibld-1)
                b_sectopidx[ibld,isec] = startlenmeshz+1 + length(bldsec_Z[isec])-1#*ibld-1
            else
                b_secbotidx[ibld,isec] = b_sectopidx[ibld,isec-1]+1 #+ length(bldsec_Z[isec])*(ibld-1)
                b_sectopidx[ibld,isec] = b_sectopidx[ibld,isec-1]+1 + length(bldsec_Z[isec])-1#*ibld-1
            end

            # hubIdx = length(mesh_z) + length(bld_Z)*(ibld-1) + hubIdx

            # Intraconnectivity
            conn_b = zeros(length(bldsec_Z[isec])-1,2)
            conn_b[:,1] = collect(b_secbotidx[ibld,isec]:1:b_sectopidx[ibld,isec]-1)
            conn_b[:,2] = collect(b_secbotidx[ibld,isec]+1:1:b_sectopidx[ibld,isec])
            conn = [conn;conn_b]
        end

        #  # find the closest point to hub radius and move it to exactly match
        #  b_xstart = bld_X[1] #mesh_x[b_startidx]
        #  b_ystart = bld_Y[1] #mesh_y[b_startidx]
        #  b_zstart = bld_Z[1] #mesh_z[b_startidx]
 
        #  b_xend = bld_X[end] #mesh_x[b_endidx]
        #  b_yend = bld_Y[end] #mesh_y[b_endidx]
        #  b_zend = bld_Z[end] #mesh_z[b_endidx]
        #  hubIdx=1
        #  if AD15hubR > 1e-6
        #      lenXY = sqrt((b_xend - b_xstart)^2 + (b_yend - b_ystart)^2)   # strut length in XY
        #      minR2 = lenXY 
        #      for i = 1:nbelem+1  # step through to find closest point to hub radius on x-y plane
        #          R2 = AD15hubR - sqrt((bld_X[i] - b_xstart)^2 + (bld_Y[i] - b_ystart)^2)
        #          if abs(R2) < abs(minR2)
        #              hubIdx = i
        #              minR2 = R2
        #          end
        #      end
        #      R_temp = minR2
        #      println("Hub crossing at idx $hubIdx at $R_temp with hub radius of $AD15hubR")
        #      print("Moving strut point from [$(bld_X[hubIdx]),$(bld_Y[hubIdx]),$(bld_Z[hubIdx])] to ")
        #      bld_X[hubIdx] = bld_X[hubIdx] + R_temp/lenXY*(b_xend-b_xstart)
        #      bld_Y[hubIdx] = bld_Y[hubIdx] + R_temp/lenXY*(b_yend-b_ystart)
        #      bld_Z[hubIdx] = bld_Z[hubIdx] + R_temp/lenXY*(b_zend-b_zstart)
        #      print("[$(bld_X[hubIdx]),$(bld_Y[hubIdx]),$(bld_Z[hubIdx])]\n")
        #  end
        #  # set index of where this point is in the mesh
         

        # if AD15_ccw
        #     AD15bldNdIdxRng = [AD15bldNdIdxRng; hubIdx b_topidx[ibld]]    # top of blade is root
        # else
        #     AD15bldNdIdxRng = [AD15bldNdIdxRng; b_topidx[ibld] hubIdx]    # bottom of blade is root
        # end
    end

    #######################################
    ###   Cleanup/Derived parameters    ###
    #######################################

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:t_topidx] .= 1 #Tower
    # meshtype[idx_bot_lbld_tower:idx_top_rbld_tower] .= 0 #Blades

    #########################
    # .bld equivalent
    #########################
    #TODO: fix this when mapping aero
    # For a single blade
    meshSeg = zeros(1+nblade*4) #tower, blades, and cables

    meshSeg[1] = ntelem
    meshSeg[2:4:end] .= nbelem_root
    meshSeg[3:4:end] .= nbelem_biwing
    meshSeg[4:4:end] .= nbelem_biwing
    meshSeg[5:4:end] .= nbelem_tip

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z))
    structuralNodeNumbers = zeros(nblade,length(bld_Z))
    structuralElNumbers = zeros(nblade,length(bld_Z))

    for iblade = 1:nblade

        # Normalized Span
        span_len = sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-hub_depth).^2.0)
        structuralSpanLocNorm[iblade,:] = span_len./maximum(span_len)

        # Node Numbers
        structuralNodeNumbers[iblade,:] = collect(b_secbotidx[iblade,1]:b_sectopidx[iblade,4])

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

    jointconn = zeros(Int,nblade*5,2)
    for ibld = 1:nblade
        # connect tower to blades
        jointconn[ibld,:] = [t_topidx b_secbotidx[ibld,1]]
        
        # Connect The blade root to the upper biwing
        jointconn[ibld+nblade,:] = [b_sectopidx[ibld,1] b_secbotidx[ibld,2]]

        # Connect The blade root to the lower biwing
        jointconn[ibld+nblade*2,:] = [b_sectopidx[ibld,1] b_secbotidx[ibld,3]]

        # Connect The upper biwing to the tip
        jointconn[ibld+nblade*3,:] = [b_sectopidx[ibld,2] b_secbotidx[ibld,4]]

        # Connect The lower biwing to the tip
        jointconn[ibld+nblade*4,:] = [b_sectopidx[ibld,3] b_secbotidx[ibld,4]]

    end

    # PyPlot.figure()
    # PyPlot.plot(mymesh.x,mymesh.z,"b.")
    #  for myi = 1:length(mymesh.x)
    #      PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
    #      PyPlot.draw()
    #      sleep(0.1)
    #  end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.axis("equal")

    # PyPlot.figure()
    # PyPlot.scatter3D(mymesh.x,mymesh.y,mymesh.z,"b.")
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.zlabel("z")
    # PyPlot.axis("equal")

    njoint = length(jointconn[:,1])
    ort = OWENS.calculateElementOrientation(mymesh)
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
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint).+joint_type zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    AD15bldElIdxRng = zeros(Int64,0,2)
    # for i = 1:size(AD15bldNdIdxRng,1)
    #     if AD15bldNdIdxRng[i,2] > AD15bldNdIdxRng[i,1]  # ascending order
    #         idx1 = findfirst(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,1])
    #         idx2 = findfirst(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
    #     else    # upside down oriented blade
    #         idx1 = findlast(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,2])
    #         idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,1])
    #     end
    #     AD15bldElIdxRng = [AD15bldElIdxRng; idx1 idx2]
    #     #println("Idx: $(AD15bldNdIdxRng[i,:])    idx: $([idx1 idx2])")
    # end

    return mymesh, ort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng
end
