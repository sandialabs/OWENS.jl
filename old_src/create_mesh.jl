using PyPlot
close("all")
include("readMesh.jl")
mesh = readMesh("old_src/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh")

# blade

Ht = 15.0 #tower height before blades attach
Hb = 147.148-15.0 #blade height
R = 55.0 # m bade radius
Nstrut = 2
strut_mout_ratio = 0.1 #distance from top/bottom

# create points
nbelem = 22 #blade elements
v_space = Hb/nbelem
ntelem = 3 #tower elements

# Tower

#                   tower_x             up to bottom strut
tower = collect(LinRange(0,Ht,ntelem))
base_idx = ntelem

start = tower[end]+v_space
stop = strut_mout_ratio*Hb+Ht
tower2strut = collect(LinRange(start,stop,max(round(Int,(stop-start)/v_space),2)))
b_tower_strut_idx = base_idx+max(round(Int,(stop-start)/v_space),2)

start = tower2strut[end]+v_space
stop = (1-strut_mout_ratio)*Hb+Ht
strut2strut = collect(LinRange(start,stop,max(round(Int,(stop-start)/v_space),2)))
t_tower_strut_idx = b_tower_strut_idx+max(round(Int,(stop-start)/v_space),2)

start = strut2strut[end]+v_space
stop = Hb+Ht
strut2top = collect(LinRange(start,stop,max(round(Int,(stop-start)/v_space),2)))
top_idx = t_tower_strut_idx + max(round(Int,(stop-start)/v_space),2)

mesh_z = [tower;tower2strut;strut2strut;strut2top]
mesh_x = zeros(length(mesh_z))
mesh_y = zeros(length(mesh_z))

# intra-tower connectivity
conn = zeros(length(mesh_z)-1,2)
conn[:,1] = collect(1:length(mesh_z)-1)
conn[:,2] = collect(2:length(mesh_z))

# Blades

#connection points on tower already exist, as do the strut connection points if we reuse the tower z discretization as the blade discretization
bld_Z = mesh_z[ntelem+1:end-1].-Ht
bld_X = R.*(1.0.-4.0.*(bld_Z/Hb.-.5).^2)
bld_Y = zero(bld_X)
bld_Z .+= Ht

#TODO: arbitrary number of blades rotated about centers
l_blade_start_idx = top_idx+1
r_blade_start_idx = l_blade_start_idx+length(bld_Z)
mesh_z = [mesh_z;bld_Z;bld_Z]
mesh_x = [mesh_x;-bld_X;bld_X]
mesh_y = [mesh_y;bld_Y;bld_Y]
l_blade_strut_bot_idx = b_tower_strut_idx-ntelem+top_idx
l_blade_strut_top_idx = t_tower_strut_idx-ntelem+top_idx
r_blade_strut_bot_idx = b_tower_strut_idx-ntelem+r_blade_start_idx-1
r_blade_strut_top_idx = t_tower_strut_idx-ntelem+r_blade_start_idx-1
# intra-blade connectivity and blade-tower connections
# Left Blade
conn_b = zeros(length(bld_Z)+2-1,2)
#              connect l blade bottom, lblade, connect l top blade, connect r blade bottom, rblade, connect rblade top
conn_b[:,1] = [base_idx;collect(l_blade_start_idx:1:l_blade_start_idx+length(bld_Z)-1)]
conn_b[:,2] = [collect(l_blade_start_idx:1:l_blade_start_idx+length(bld_Z)-1);top_idx]
conn = [conn;conn_b]

# Right Blade
conn_b[:,1] = [base_idx;collect(r_blade_start_idx:1:r_blade_start_idx+length(bld_Z)-1)]
conn_b[:,2] = [collect(r_blade_start_idx:1:r_blade_start_idx+length(bld_Z)-1);top_idx]
conn = [conn;conn_b]

#Connection Within Blades
# conn_b = zeros(,2)
# conn_b[:,1] =
# conn_b[:,2] =

# Struts
bot_strut_z = strut_mout_ratio*Hb+Ht
top_strut_z = (1-strut_mout_ratio)*Hb+Ht
bot_strut_x = R.*(1.0.-4.0.*((bot_strut_z-Ht)/Hb.-.5).^2)/2
top_strut_x = R.*(1.0.-4.0.*((top_strut_z-Ht)/Hb.-.5).^2)/2
bot_left_strut_idx = length(mesh_z)+1
bot_right_strut_idx = length(mesh_z)+2
top_left_strut_idx = length(mesh_z)+3
top_right_strut_idx = length(mesh_z)+4
mesh_z = [mesh_z;bot_strut_z;bot_strut_z;top_strut_z;top_strut_z]
mesh_x = [mesh_x;-bot_strut_x;bot_strut_x;-top_strut_x;top_strut_x]
mesh_y = [mesh_y;0.0;0.0;0.0;0.0]

# Struts only have one point inbetween, so we just connect them to the tower and blades
conn_s = zeros(2,2)
# Lower Left Strut
            #point on tower, strut point, point on blade
conn_s[:,1] = [b_tower_strut_idx;bot_left_strut_idx]
conn_s[:,2] = [bot_left_strut_idx;l_blade_strut_bot_idx]
conn = [conn;conn_s]

# Lower Right Strut
            #point on tower, strut point, point on blade
conn_s[:,1] = [b_tower_strut_idx;bot_right_strut_idx]
conn_s[:,2] = [bot_right_strut_idx;r_blade_strut_bot_idx]
conn = [conn;conn_s]

# Top Left Strut
            #point on tower, strut point, point on blade
conn_s[:,1] = [t_tower_strut_idx;top_left_strut_idx]
conn_s[:,2] = [top_left_strut_idx;l_blade_strut_top_idx]
conn = [conn;conn_s]

# Top Right Strut
            #point on tower, strut point, point on blade
conn_s[:,1] = [t_tower_strut_idx;top_right_strut_idx]
conn_s[:,2] = [top_right_strut_idx;r_blade_strut_top_idx]
conn = [conn;conn_s]


figure()
plot(mesh_x,mesh_z,".-")
plot(mesh.x,mesh.z,".-")
axis("equal")

figure()
for i = 1:length(conn[:,1])
    println("$(conn[i,1]) $(conn[i,2])")
    plot([mesh_x[Int(conn[i,1])];mesh_x[Int(conn[i,2])]],[mesh_z[Int(conn[i,1])];mesh_z[Int(conn[i,2])]],".-")
    axis("equal")
    pause(0.5)
end

# figure()
# for i = 1:length(mesh.x)-1
#     plot(mesh.x[i:i+1],mesh.z[i:i+1],".-")
#     # pause(0.9)
# end
# plot(shapeX,shapeZ,".-")
# axis("equal")
