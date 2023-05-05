using GXBeam, LinearAlgebra, FLOWMath, ModelGen
path = splitdir(@__FILE__)[1]

# --- Tower Definition --- #

tower_height = 41.9
tower_stiffness = Diagonal([8.4404e9, 2.7053e9, 2.7053e9, 7.0428e9, 9.4072e9, 9.4072e9])
tower_mass = Diagonal([330.755, 330.755, 330.755, 737.282, 368.641, 368.641])

# --- Blade Definition --- #

# geometry
blade_xyz = [
     0.0        0.0   0.0;
     2.26837    0.0   1.257;
     3.63183    0.0   2.095;
     6.76113    0.0   4.19;
     9.55882    0.0   6.285;
    11.8976     0.0   8.38;
    13.7306     0.0  10.475;
    15.1409     0.0  12.57;
    16.1228     0.0  14.665;
    16.7334     0.0  16.76;
    17.0133     0.0  18.855;
    17.0987     0.0  20.95;
    16.9615     0.0  23.045;
    16.5139     0.0  25.14;
    15.7435     0.0  27.235;
    14.5458     0.0  29.33;
    12.9287     0.0  31.425;
    10.8901     0.0  33.52;
     8.547      0.0  35.615;
     5.93739    0.0  37.71;
     3.05842    0.0  39.805;
     1.87486    0.0  40.643;
     0.0        0.0  41.9
]

# section boundaries (z-coordinate)
blade_transition = [0.0, 5.8, 11.1, 29.0, 34.7, 41.9]

# root section properties

blade_stiffness1 = [
    3.74563e9  0.0        0.0        0.0        0.0        1.14382e8;
    0.0        1.20052e9  0.0        0.0        0.0        0.0;
    0.0        0.0        1.20052e9  0.0        0.0        0.0;
    0.0        0.0        0.0        1.87992e7  0.0        0.0;
    0.0        0.0        0.0        0.0        2.24336e7  0.0;
    1.14382e8  0.0        0.0        0.0        0.0        4.09242e8;
]

blade_mass1 = [
    146.781     0.0       0.0      0.0     28.7783     0.0;
      0.0     146.781     0.0    -28.7783   0.0        0.0;
      0.0       0.0     146.781    0.0      0.0        0.0;
      0.0     -28.7783    0.0     16.7793   0.0        0.0;
     28.7783    0.0       0.0      0.0      0.879112  -0.0;
      0.0       0.0       0.0      0.0     -0.0       15.9002;
]

# transition section properties

blade_stiffness2 = [
    2.22783e9  0.0        0.0        0.0        0.0        4.20422e7;
    0.0        7.14048e8  0.0        0.0        0.0        0.0;
    0.0        0.0        7.14048e8  0.0        0.0        0.0;
    0.0        0.0        0.0        6.55493e6  0.0        0.0;
    0.0        0.0        0.0        0.0        7.35548e6  0.0;
    4.20422e7  0.0        0.0        0.0        0.0        1.84227e8;
]

blade_mass2 = [
    87.3025    0.0      0.0       0.0      16.7316     0.0;
     0.0      87.3025   0.0     -16.7316    0.0        0.0;
     0.0       0.0     87.3025   -0.0       0.0        0.0;
     0.0     -16.7316  -0.0       7.47649   0.0        0.0;
    16.7316    0.0      0.0       0.0       0.288241  -0.0;
     0.0       0.0      0.0       0.0      -0.0        7.18825;
]

# center section properties

blade_stiffness3 = [
    1.76888e9  0.0        0.0        0.0        0.0        2.34071e7;
    0.0        5.66947e8  0.0        0.0        0.0        0.0;
    0.0        0.0        5.66947e8  0.0        0.0        0.0;
    0.0        0.0        0.0        4.00804e6  0.0        0.0;
    0.0        0.0        0.0        0.0        4.34302e6  0.0;
    2.34071e7  0.0        0.0        0.0        0.0        1.09341e8;
]

blade_mass3 = [
    69.3173    0.0      0.0       0.0      11.5831     0.0;
     0.0      69.3173   0.0     -11.5831    0.0        0.0;
     0.0       0.0     69.3173   -0.0       0.0        0.0;
     0.0     -11.5831  -0.0       4.44282   0.0        0.0;
    11.5831    0.0      0.0       0.0       0.170191  -0.0;
     0.0       0.0      0.0       0.0      -0.0        4.27263;
]

# --- Strut Definition --- #

strut_locations = [1.257, 40.643]
strut_stiffness = blade_stiffness1
strut_mass = blade_mass1

# --- Define Assembly --- #

# Tower

# number of tower sections
nt = 22

# tower points
x = zeros(nt+1)
y = zeros(nt+1)
z = vcat(0, range(strut_locations[1], strut_locations[2]; length=nt-1), tower_height)
pt_t = [[x[i],y[i],z[i]] for i = 1:nt+1]

# tower frame of reference
Twist_d_el = -90.0
# apply the twist rotation, which is about the x (1) axis
DCM_roll = [1.0 0.0 0.0
0.0 cosd(Twist_d_el) -sind(Twist_d_el)
0.0 sind(Twist_d_el) cosd(Twist_d_el)]
frame_t = fill([0 0 -1; 0 1 0; 1 0 0]*DCM_roll, nt)

# tower stiffness
stiff_t = fill(tower_stiffness, nt)

# tower mass
mass_t = fill(tower_mass, nt)

# Blades

# number of blade sections
nbr = 4 # root
nbt = 3 # transition
nbc = 8 # center
nb = 2*nbr + 2*nbt + nbc # total number of blade sections

# interpolation parameter coordinates
new_z = vcat(0.0,
    range(strut_locations[1], 5.8, length=nbr)[1:end-1],
    range(5.8, 11.1, length=nbt+1)[1:end-1],
    range(11.1, 29.0, length=nbc+1)[1:end-1],
    range(29.0, 34.7, length=nbt+1)[1:end-1],
    range(34.7, strut_locations[2], length=nbr),
    tower_height)

# blade points
x = FLOWMath.akima(blade_xyz[:,3], blade_xyz[:,1], new_z)
y = zero(new_z)
z = new_z
pt_bl = [[-x[i],y[i],z[i]] for i = 1:nb+1] # left blade
pt_br = [[x[i],y[i],z[i]] for i = 1:nb+1] # right blade
# left blade frame of reference
frame_bl = Vector{Matrix{Float64}}(undef, nb)
for i = 1:nb
    r = pt_bl[i+1] - pt_bl[i]
    n = norm(r)
    s = r[3]/n
    c = r[1]/n
    frame_bl[i] = [c 0 -s; 0 1 0; s 0 c]
end

# right blade frame of reference
frame_br = Vector{Matrix{Float64}}(undef, nb)
for i = 1:nb
    r = pt_br[i+1] - pt_br[i]
    n = norm(r)
    s = r[3]/n
    c = r[1]/n

    Twist_d_el = 180.0
    # apply the twist rotation, which is about the x (1) axis
    DCM_roll = [1.0 0.0 0.0
    0.0 cosd(Twist_d_el) -sind(Twist_d_el)
    0.0 sind(Twist_d_el) cosd(Twist_d_el)]

    # # apply theta rotation, which is the tilt angle, or about the y (2) axis in global
    # DCM_pitch = [cosd(Theta_d_el) 0.0 sind(Theta_d_el)
    #     0.0 1.0 0.0
    #     -sind(Theta_d_el) 0.0 cosd(Theta_d_el)]

    # # apply Psi rotation, which is about Z (3) axis in global
    # DCM_yaw = [cosd(Psi_d_el) -sind(Psi_d_el) 0.0
    #     sind(Psi_d_el) cosd(Psi_d_el) 0.0
    #     0.0 0.0 1.0]

    frame_br[i] = [c 0 -s; 0 1 0; s 0 c]*DCM_roll
end


# blade stiffness
stiff_b = vcat(
    fill(blade_stiffness1, nbr),
    fill(blade_stiffness2, nbt),
    fill(blade_stiffness3, nbc),
    fill(blade_stiffness2, nbt),
    fill(blade_stiffness1, nbr)
)

# blade mass
mass_b = vcat(
    fill(blade_mass1, nbr),
    fill(blade_mass2, nbt),
    fill(blade_mass3, nbc),
    fill(blade_mass2, nbt),
    fill(blade_mass1, nbr)
)

# Struts

# number of strut sections per strut
ns = 3

# lower left strut points
x = range(0.0, pt_bl[2][1]; length=ns+1)
y = zeros(ns+1)
z = fill(strut_locations[1], ns+1)
pt_s1 = [[x[i],y[i],z[i]] for i = 1:ns+1]

# lower right strut points
x = range(0.0, pt_br[2][1]; length=ns+1)
y = zeros(ns+1)
z = fill(strut_locations[1], ns+1)
pt_s2 = [[x[i],y[i],z[i]] for i = 1:ns+1]

# upper left strut points
x = range(0.0, pt_bl[end-1][1]; length=ns+1)
y = zeros(ns+1)
z = fill(strut_locations[2], ns+1)
pt_s3 = [[x[i],y[i],z[i]] for i = 1:ns+1]

# upper right strut points
x = range(0.0, pt_br[end-1][1]; length=ns+1)
y = zeros(ns+1)
z = fill(strut_locations[2], ns+1)
pt_s4 = [[x[i],y[i],z[i]] for i = 1:ns+1]

# strut frame of reference
frame_sR = fill([1 0 0; 0 -1 0; 0 0 -1], ns)
frame_sL = fill([-1 0 0; 0 1 0; 0 0 -1], ns)

# strut stiffness
stiff_s = fill(strut_stiffness, ns)

# strut mass
mass_s = fill(strut_mass, ns)

# Combine Tower, Blades, and Struts

# combine points
points = vcat(pt_t, pt_bl, pt_br, pt_s1, pt_s2, pt_s3, pt_s4)

# define element connectivity
istart = cumsum([1, nt+1, nb+1, nb+1, ns+1, ns+1, ns+1])
istop = cumsum([nt+1, nb+1, nb+1, ns+1, ns+1, ns+1, ns+1])
start = vcat([istart[i]:istop[i]-1 for i = 1:length(istart)]...)
stop = vcat([istart[i]+1:istop[i] for i = 1:length(istart)]...)

# use zero-length elements as joints

nj = 12 # number of joints

joints = [
    istart[1]     istart[2]; # tower - bottom of left blade
    istart[1]     istart[3]; # tower - bottom of right blade
    istart[1]+1   istart[4]; # tower - lower left strut
    istart[1]+1   istart[5]; # tower - lower right strut
    istop[1]-1    istart[6]; # tower - upper left strut
    istop[1]-1    istart[7]; # tower - upper right strut
    istop[2]      istop[1]; # top of left blade - tower
    istop[3]      istop[1]; # top of right blade - tower
    istop[4]      istart[2]+1; # lower left strut - left blade
    istop[5]      istart[3]+1; # lower right strut - right blade
    istop[6]      istop[2]-1; # upper left strut - left blade
    istop[7]      istop[3]-1; # upper right strut - right blade
]

frame_j = fill([1 0 0; 0 1 0; 0 0 1], nj)

stiff_j = fill(zeros(6,6), nj) # will be modeled as infinitely stiff

mass_j = fill(zeros(6,6), nj)

# add joint connectivity
start = vcat(start, joints[:,1])
stop = vcat(stop, joints[:,2])

# combine frames
frames = vcat(frame_t, frame_bl, frame_br, frame_sL, frame_sR, frame_sL, frame_sR, frame_j)

# combine stiffness
stiffness = vcat(stiff_t, stiff_b, stiff_b, stiff_s, stiff_s, stiff_s, stiff_s, stiff_j)

# combine mass
mass = vcat(mass_t, mass_b, mass_b, mass_s, mass_s, mass_s, mass_s, mass_j)

xplot = [point[1] for point in points]
yplot = [point[2] for point in points]
zplot = [point[3] for point in points]

function rotate_normal(i_el, points, start, stop, frames;vec=[0,1,0.0],normal_len=1)

    # Map from element to node
    nodenum1 = start[i_el]
    nodenum2 = stop[i_el]

    # Use a line and rotate it about the angles, different starting vectors show different angles.
    myvec = vec.*normal_len

    xplot1 = (points[nodenum1][1]+points[nodenum2][1])/2
    yplot1 = (points[nodenum1][2]+points[nodenum2][2])/2
    zplot1 = (points[nodenum1][3]+points[nodenum2][3])/2

    # Offset the myvector by the location of the element
    myvec = frames[i_el]*myvec .+ [xplot1,yplot1,zplot1]
    x_el_plot1 = [xplot1,myvec[1]]
    y_el_plot1 = [yplot1,myvec[2]]
    z_el_plot1 = [zplot1,myvec[3]]
    return x_el_plot1, y_el_plot1, z_el_plot1
end
import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(15, 15))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

PyPlot.figure()
PyPlot.scatter3D(xplot,yplot,zplot,color=plot_cycle[1])

 # Add the orientation vectors, ort is on elements
 for i_el = 1:length(start)
    x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[10,0,0.0])
    PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[2])
end

for i_el = 1:length(start)
    x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[0,5,0.0])
    PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[3])
end

for i_el = 1:length(start)
    x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[0,0,5.0])
    PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[4])
end

PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[2],label="X-norm")
PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[3],label="Y-norm")
PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[4],label="Z-norm")
PyPlot.legend()

PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")
PyPlot.axis("equal")

# create assembly
assembly = Assembly(points, start, stop;
    frames=frames,
    stiffness=stiffness,
    mass=mass)

# --- Define Prescribed Conditions --- #

# create dictionary of prescribed conditions
prescribed_conditions = Dict(
    # fixed base
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    # fixed top, but free to rotate around z-axis
    istop[1] => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0),
)

# --- Perform Analysis --- #

# revolutions per minute
rpm = 0:1:40

# gravity vector
gravity = [0, 0, -9.81]

# number of modes
nmode = 10

# number of eigenvalues
nev = 2*nmode

# initialize system storage
system = DynamicSystem(assembly)


println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","Fx","Fy","Fz","Mx","My","Mz"]
t = [0,1]

# map el props to points using con
userPointData = zeros(length(userPointNames),length(t),length(points))
EA_points = zeros(length(points))
EIyy_points = zeros(length(points))
EIzz_points = zeros(length(points))

# Time-invariant data
for iel = 1:length(start)
    node1 = start[iel]
    node2 = stop[iel]
    EA_points[node1] = stiffness[iel][1,1]
    EA_points[node2] = stiffness[iel][1,1]
    EIyy_points[node1] = stiffness[iel][2,2]
    EIyy_points[node2] = stiffness[iel][2,2]
    EIzz_points[node1] = stiffness[iel][3,3]
    EIzz_points[node2] = stiffness[iel][3,3]
end

# fill in the big matrix
for it = 1:length(t)

    userPointData[1,it,:] = EA_points
    userPointData[2,it,:] = EIyy_points
    userPointData[3,it,:] = EIzz_points
    # userPointData[4,it,:] = FReactionHist[it,1:6:end]
    # userPointData[5,it,:] = FReactionHist[it,2:6:end]
    # userPointData[6,it,:] = FReactionHist[it,3:6:end]
    # userPointData[7,it,:] = FReactionHist[it,4:6:end]
    # userPointData[8,it,:] = FReactionHist[it,5:6:end]
    # userPointData[9,it,:] = FReactionHist[it,6:6:end]
end

azi=[0.0,pi/8]#./aziHist*1e-6
uHist = [zeros(length(points)*6) zeros(length(points)*6)]'
saveName = "$path/vtk/pleasework"

af = [100.0001 0
100 0.221
95 1.412
90 2.534
80 4.591
70 6.412
60 7.986
50 9.265
40 10.156
30 10.504
25 10.397
20 10.04
15 9.354
10 8.195
7.5 7.35
5 6.221
2.5 4.576
1.25 3.315
0 0
1.25 -3.315
2.5 -4.576
5 -6.221
7.5 -7.35
10 -8.195
15 -9.354
20 -10.04
25 -10.397
30 -10.504
40 -10.156
50 -9.265
60 -7.986
70 -6.412
80 -4.591
90 -2.534
95 -1.412
100 -0.221
100.0001 0]./100.0

sections = zeros(3,length(af[:,1]),length(assembly.points))

for ipt = 1:length(assembly.points)

    elNum = findfirst(x->x==ipt,start) # Get element number
    if isnothing(elNum)
        elNum = findfirst(x->x==ipt,stop) # Get element number
    end

    myaf = frames[elNum]*[zero(af[:,1]) (af[:,1].-maximum(af[:,1])/2) af[:,2]]'
    sections[1,:,ipt] = myaf[1,:]
    sections[2,:,ipt] = myaf[2,:]
    sections[3,:,ipt] = myaf[3,:]

end

ModelGen.gyricFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)


# # storage for results
# freq = zeros(length(rpm), nmode)

# # perform an analysis for each rotation rate
# for (i,rpm) in enumerate(rpm)

#     global system, Up

#     # set turbine rotation
#     angular_velocity = [0, 0, rpm*(2*pi)/60]

#     # eigenvalues and (right) eigenvectors
#     system, λ, V, converged = eigenvalue_analysis!(system, assembly;
#         prescribed_conditions = prescribed_conditions,
#         angular_velocity = angular_velocity,
#         gravity = gravity,
#         nev = nev
#         )

#     # check convergence
#     @assert converged

#     if i > 1
#         # construct correlation matrix
#         C = Up*system.M*V

#         # correlate eigenmodes
#         perm, corruption = correlate_eigenmodes(C)

#         # re-arrange eigenvalues
#         λ = λ[perm]

#         # update left eigenvector matrix
#         Up = left_eigenvectors(system, λ, V)
#         Up = Up[perm,:]
#     else
#         # update left eigenvector matrix
#         Up = left_eigenvectors(system, λ, V)
#     end

#     # save frequencies
#     freq[i,:] = [imag(λ[k])/(2*pi) for k = 1:2:nev]

# end