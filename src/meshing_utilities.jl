
struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
end

function create_mesh(;Ht = 15.0, #tower height before blades attach
    Hb = 147.148-15.0, #blade height
    R = 55.0, # m bade radius
    nstrut = 2,
    strut_mout_ratio = 0.1, #distance from top/bottom
    ntelem = 20, #tower elements
    nbelem = 20, #blade elements
    nselem = 2) #strut elements

    v_space_t = (Ht+Hb)/(ntelem)
    v_space_b = (Hb)/(nbelem)

    ####################################
    ##------------Tower--------------##
    ####################################
    function discretize_z(Ht,Hb,nelem; offset=0)
        mesh_z_inner = collect(LinRange(0,Ht+Hb,nelem+1))
        # Insert Bottom Blade Mount Point
        bot_tower_bld_idx = round(Int,Ht/v_space_t+1.0+length(mesh_z_inner)-nelem)-1
        if offset ==0 #but only if it isn't being used for the blade
            insert!(mesh_z_inner,bot_tower_bld_idx,Ht)
            bot_tower_bld_idx+=1
        end

        # Insert Bottom Strut Mount Point
        strut_z = Ht+strut_mout_ratio*Hb
        bot_tower_strut_idx = round(Int,strut_z/v_space_t+length(mesh_z_inner)-nelem)
        insert!(mesh_z_inner,bot_tower_strut_idx,strut_z)

        # Insert Top Strut Mount Point
        strut_z = Ht+(1-strut_mout_ratio)*Hb
        extra_offset = 1
        if offset != 0
            extra_offset = 2
        end
        top_tower_strut_idx = round(Int,strut_z/v_space_t+extra_offset+length(mesh_z_inner)-nelem)
        insert!(mesh_z_inner,top_tower_strut_idx,strut_z)

        # Don't Insert Top Blade Mount Point since it is already there at the exact spot
        top_tower_blade_idx = length(mesh_z_inner)
        return mesh_z_inner, bot_tower_bld_idx+offset, bot_tower_strut_idx+offset, top_tower_strut_idx+offset, top_tower_blade_idx+offset
    end
    mesh_z, bot_tower_bld_idx, bot_tower_strut_idx, top_tower_strut_idx, top_tower_blade_idx = discretize_z(Ht,Hb,ntelem)
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
    bld_Z, bot_lbld_tower_idx, bot_lbld_strut_idx, top_lbld_strut_idx, top_lbld_tower_idx = discretize_z(0.0,Hb,nbelem; offset = length(mesh_z))

    bot_rbld_tower_idx = bot_lbld_tower_idx + nbelem+1+nstrut
    bot_rbld_strut_idx = bot_lbld_strut_idx + nbelem+1+nstrut
    top_rbld_strut_idx = top_lbld_strut_idx + nbelem+1+nstrut
    top_rbld_tower_idx = top_lbld_tower_idx + nbelem+1+nstrut
    # bld_Z = mesh_z[bot_tower_bld_idx:end].-Ht
    bld_X = R.*(1.0.-4.0.*(bld_Z/Hb.-.5).^2)
    bld_Y = zero(bld_X)
    bld_Z .+= Ht

    #TODO: arbitrary number of blades rotated about center
    mesh_z = [mesh_z;bld_Z;bld_Z]
    mesh_x = [mesh_x;-bld_X;bld_X]
    mesh_y = [mesh_y;bld_Y;bld_Y]


    # intra-blade connectivity and blade-tower connections
    # Left Blade
    conn_b = zeros(length(bld_Z)-1,2)
    conn_b[:,1] = collect(bot_lbld_tower_idx:1:top_lbld_tower_idx-1)
    conn_b[:,2] = collect(bot_lbld_tower_idx+1:1:top_lbld_tower_idx)
    conn = [conn;conn_b]

    # Right Blade
    conn_b[:,1] = collect(bot_rbld_tower_idx:1:top_rbld_tower_idx-1)
    conn_b[:,2] = collect(bot_rbld_tower_idx+1:1:top_rbld_tower_idx)
    conn = [conn;conn_b]


    ####################################
    ##------------Struts--------------##
    ####################################

    bot_strut_z = strut_mout_ratio*Hb+Ht
    top_strut_z = (1-strut_mout_ratio)*Hb+Ht
    bot_strut_x = R.*(1.0.-4.0.*((bot_strut_z-Ht)/Hb.-.5).^2)
    top_strut_x = R.*(1.0.-4.0.*((top_strut_z-Ht)/Hb.-.5).^2)

    bot_strut_X = collect(LinRange(0.0,bot_strut_x,nselem+1))
    top_strut_X = collect(LinRange(0.0,top_strut_x,nselem+1))

    bot_lstrut_tower_idx = length(mesh_z)+1
    bot_lstrut_blade_idx = bot_lstrut_tower_idx + 2
    bot_rstrut_tower_idx = bot_lstrut_blade_idx + 1
    bot_rstrut_blade_idx = bot_rstrut_tower_idx + 2
    top_lstrut_tower_idx = bot_rstrut_blade_idx + 1
    top_lstrut_blade_idx = top_lstrut_tower_idx + 2
    top_rstrut_tower_idx = top_lstrut_blade_idx + 1
    top_rstrut_blade_idx = top_rstrut_tower_idx + 2
    mesh_z = [mesh_z;ones(nselem+1).*bot_strut_z ; ones(nselem+1).*bot_strut_z ; ones(nselem+1).*top_strut_z ; ones(nselem+1).*top_strut_z]
    mesh_x = [mesh_x;-bot_strut_X;bot_strut_X;-top_strut_X;top_strut_X]
    mesh_y = [mesh_y;zeros((nselem+1)*4)]

    conn_s = zeros(8,2)

    conn_s[:,1] = [bot_lstrut_tower_idx,bot_lstrut_tower_idx+1,
    bot_rstrut_tower_idx,bot_rstrut_tower_idx+1,
    top_lstrut_tower_idx,top_lstrut_tower_idx+1,
    top_rstrut_tower_idx,top_rstrut_tower_idx+1]

    conn_s[:,2] = [bot_lstrut_tower_idx+1,bot_lstrut_blade_idx,
    bot_rstrut_tower_idx+1,bot_rstrut_blade_idx,
    top_lstrut_tower_idx+1,top_lstrut_blade_idx,
    top_rstrut_tower_idx+1,top_rstrut_blade_idx]

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
    mymesh = Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,conn)

    return mymesh
end

function readMesh(filename)
#readMesh  reads mesh file and stores data in mesh object
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [mesh] = readMesh(filename)
#
#   This function reads the mesh file and stores data in the mesh object.
#
#      input:
#      filename      = string containing mesh filename

#      output:
#      mesh          = object containing mesh data

fid = open(filename,"r")   #open mesh file

# temp = fscanf(fid,'#i',2)   #read in number of nodes and number of elements
line = readline(fid)
temp = split(line)

numNodes = parse(Int,temp[1])
numEl = parse(Int,temp[2])

nodeNum = zeros(numNodes,1)
x = zeros(numNodes,1)
y = zeros(numNodes,1)
z = zeros(numNodes,1)

conn = zeros(numEl,2)
elNum = zeros(numEl,1)

for i=1:numNodes            # read in node number and node coordinates
    line = readline(fid)
    temp = split(line)
    nodeNum[i] = parse(Float64,temp[1])
    x[i] = parse(Float64,temp[2])
    y[i] = parse(Float64,temp[3])
    z[i] = parse(Float64,temp[4])
end

for i=1:numEl               # read in element number and connectivity list
    line = readline(fid)
    temp = split(line)
    elNum[i] = parse(Float64,temp[1])

    conn[i,1] = parse(Float64,temp[3])
    conn[i,2] = parse(Float64,temp[4])
end

close(fid)  #close mesh file

mesh = Mesh(nodeNum,
numEl,
numNodes,
x,
y,
z,
elNum,
conn)

return mesh

end
