import MAT
using Test
import DelimitedFiles
import ModelGen
import FLOWMath
import GyricFEA
using GXBeam
import OWENS
path = splitdir(@__FILE__)[1]


#Put in one place so its not repeated for all of the analyses
include("$(path)/34mSetup.jl")

system, assembly, sections, frames, points, start, stop, stiff, mass = ModelGen.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld;damp_coef=0.01,VTKmeshfilename="$path/vtk/SNL34m")
println("setup complete")
# function rotate_normal(i_el, points, start, stop, frames;vec=[0,1,0.0],normal_len=1)

#     # Map from element to node
#     nodenum1 = start[i_el]
#     nodenum2 = stop[i_el]

#     # Use a line and rotate it about the angles, different starting vectors show different angles.
#     myvec = vec.*normal_len

#     xplot1 = (points[nodenum1][1]+points[nodenum2][1])/2
#     yplot1 = (points[nodenum1][2]+points[nodenum2][2])/2
#     zplot1 = (points[nodenum1][3]+points[nodenum2][3])/2

#     # Offset the myvector by the location of the element
#     myvec = frames[i_el]*myvec .+ [xplot1,yplot1,zplot1]
#     x_el_plot1 = [xplot1,myvec[1]]
#     y_el_plot1 = [yplot1,myvec[2]]
#     z_el_plot1 = [zplot1,myvec[3]]
#     return x_el_plot1, y_el_plot1, z_el_plot1
# end
# import PyPlot
# PyPlot.close("all")
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(15, 15))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=4.0)
# PyPlot.rc("legend", frameon=true)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# xplot = [point[1] for point in points]
# yplot = [point[2] for point in points]
# zplot = [point[3] for point in points]

# PyPlot.figure()
# PyPlot.scatter3D(xplot,yplot,zplot,color=plot_cycle[1])

#  # Add the orientation vectors, ort is on elements
#  for i_el = 1:length(start)
#     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[10,0,0.0])
#     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[2])
# end

# for i_el = 1:length(start)
#     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[0,5,0.0])
#     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[3])
# end

# for i_el = 1:length(start)
#     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, points, start, stop, frames;vec=[0,0,5.0])
#     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[4])
# end

# PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[2],label="X-norm")
# PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[3],label="Y-norm")
# PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[4],label="Z-norm")
# PyPlot.legend()

# PyPlot.xlabel("x")
# PyPlot.ylabel("y")
# PyPlot.zlabel("z")
# PyPlot.axis("equal")

println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","Fx","Fy","Fz","Mx","My","Mz"]
t = [0,1]

# map el props to points using con
userPointData = zeros(length(userPointNames),length(t),mymesh.numNodes)
EA_points = zeros(mymesh.numNodes)
EIyy_points = zeros(mymesh.numNodes)
EIzz_points = zeros(mymesh.numNodes)

# Time-invariant data
for iel = 1:length(myel.props)
    nodes = mymesh.conn[iel,:]
    EA_points[Int.(nodes)] = myel.props[iel].EA
    EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
    EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
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
uHist = [zeros(mymesh.numNodes*6) zeros(mymesh.numNodes*6)]'
saveName = "$path/vtk/pleasework"
ModelGen.gyricFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)


# node, dof, bc
top_idx = Int(myjoint[7,2])
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
top_idx 1 0
top_idx 2 0
top_idx 3 0]

# filename = "$(path)/data/newmesh_34m"
# ModelGen.saveOWENSfiles(filename,mymesh,myort,myjoint,myel,pBC,numadIn_bld)


# if testModal
##############################################
# Modal Test
#############################################
displ = zeros(mymesh.numNodes*6)
numModes = 16

mymodel = GyricFEA.FEAModel(;analysisType = "M",
        # outFilename = "$path/data/outplat.out",
        joint = myjoint,
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        gravityOn=true,
        nlOn = true,
        tolerance = 1e-6,
        spinUpOn = true,
        iterationType = "NR",
        numNodes = mymesh.numNodes,
        numModes)  # number of modes to calculate)



# Campbell Diagram generation
rotSpdArrayRPM = 0.0:5.0:40 # rpm
rotorSpeedArrayHZ = rotSpdArrayRPM ./ 60.0
centStiff = true      # centripetal stiffening
NperRevLines = 8
freq = zeros(length(rotorSpeedArrayHZ),numModes)
for i=1:length(rotorSpeedArrayHZ)
    println("$i of $(length(rotorSpeedArrayHZ))")
    rotorSpeed = rotorSpeedArrayHZ[i]
    local Omega = rotorSpeed
    global displ
    OmegaStart = Omega
    
    freqtemp,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90,displ=OWENS.Modal(mymodel,mymesh,myel;displ,Omega,OmegaStart,returnDynMatrices=false)
    
    freq[i,:] = freqtemp[1:numModes]
end

elapsedtime = time() - starttime
println(freq[1,:])


########################################
############ GXBeam ####################
########################################

starttime2 = time()

# create dictionary of prescribed conditions
prescribed_conditions = Dict(
    # fixed base
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    # fixed top, but free to rotate around z-axis
    top_idx => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0),
)

# --- Perform Analysis --- #
# revolutions per minute
rpm = 0:5:40

# gravity vector
gravity = [0, 0, -9.81]

# number of modes
nmode = 10

# number of eigenvalues
nev = 2*nmode

# storage for results
freq2 = zeros(length(rpm), nmode)
λ = []
eigenstates =[]
# perform an analysis for each rotation rate
for (i,rpm) in enumerate(rpm)

    global system, Up, λ, eigenstates

    # set turbine rotation
    angular_velocity = [0, 0, rpm*(2*pi)/60]

    # eigenvalues and (right) eigenvectors
    system, λ, V, converged = eigenvalue_analysis!(system, assembly;
        prescribed_conditions = prescribed_conditions,
        angular_velocity = angular_velocity,
        gravity = gravity,
        nev = nev
        )

    # check convergence
    @assert converged

    if i > 1
        # construct correlation matrix
        C = Up*system.M*V

        # correlate eigenmodes
        perm, corruption = correlate_eigenmodes(C)

        # re-arrange eigenvalues
        λ = λ[perm]

        # update left eigenvector matrix
        Up = left_eigenvectors(system, λ, V)
        Up = Up[perm,:]
    else
        # update left eigenvector matrix
        Up = left_eigenvectors(system, λ, V)
    end

    # save frequencies
    freq2[i,:] = [imag(λ[k])/(2*pi) for k = 1:2:nev]

    state = AssemblyState(system, assembly;
    prescribed_conditions = prescribed_conditions)
    eigenstates = [AssemblyState(V[:,k],system, assembly;
    prescribed_conditions = prescribed_conditions) for k = 1:nev]

end

elapsedtime2 = time() - starttime2

println("OWENS: $elapsedtime")
println("GX: $elapsedtime2")

state = AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)

write_vtk("$path/vtk/Campbell5_34m", assembly, state, 
    λ[5], eigenstates[5]; sections,mode_scaling = 500.0)

    
freqOWENS = [freq[:,i] for i=1:4:numModes-2]
freqGX = [freq2[:,i] for i=1:2:numModes-6-2]

for imode = 1:length(freqOWENS)
    for ifreq = 1:length(freqOWENS[imode])
        println("imode $imode, ifreq $ifreq")
        atol = freqGX[imode][ifreq]*0.06
        @test isapprox(freqOWENS[imode][ifreq],freqGX[imode][ifreq];atol)
    end
end

# import PyPlot
# PyPlot.close("all")
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(4.5, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=4.0)
# PyPlot.rc("legend", frameon=true)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# PyPlot.figure()
# for i=1:1:numModes-2
#        PyPlot.plot(rotSpdArrayRPM,freq[:,i],color=plot_cycle[1],"b-") #plot mode i at various rotor speeds
# end

# # Plot 34m experimental data
# SNL34_flap = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_Flatwise.csv",',',skipstart = 0)
# SNL34_lead = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_LeadLag.csv",',',skipstart = 0)

# PyPlot.plot(SNL34_flap[:,1],SNL34_flap[:,2],"k.",label="Flapwise Gauges")
# PyPlot.plot(SNL34_lead[:,1],SNL34_lead[:,2],"kx",label="Lead-Lag Gauges")
# PyPlot.plot(0,0,"k--",label="Tower-Mode")

# SNL34_1F = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1F.csv",',',skipstart = 0)
# SNL34_1BE = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1BE.csv",',',skipstart = 0)
# SNL34_2FA = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_2FA.csv",',',skipstart = 0)
# SNL34_2FS = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_2FS.csv",',',skipstart = 0)
# SNL34_1TO = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_1TO.csv",',',skipstart = 0)
# SNL34_3F = DelimitedFiles.readdlm("$(path)/data/SAND-91-2228_Data/4.5_3F.csv",',',skipstart = 0)

# PyPlot.plot(SNL34_1F[:,1],SNL34_1F[:,2],"k-")
# PyPlot.plot(SNL34_1BE[:,1],SNL34_1BE[:,2],"k-")
# PyPlot.plot(SNL34_2FA[:,1],SNL34_2FA[:,2],"k-")
# PyPlot.plot(SNL34_2FS[:,1],SNL34_2FS[:,2],"k-")
# PyPlot.plot(SNL34_1TO[:,1],SNL34_1TO[:,2],"k--")
# PyPlot.plot(SNL34_3F[:,1],SNL34_3F[:,2],"k-")

# #plot per rev lines
# for i=1:NperRevLines
#     linex=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5]
#     liney=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5].*i./60.0
#     PyPlot.plot(linex,liney,"--k",linewidth=0.5)
#     PyPlot.annotate("$i P",xy=(0.95*linex[2],liney[2]+.05+(i-1)*.01))
# end
# PyPlot.grid()
# PyPlot.xlabel("Rotor Speed (RPM)")
# PyPlot.ylabel("Frequency (Hz)")
# PyPlot.plot(0,0,"k-",label="Experimental")
# PyPlot.plot(0,0,color=plot_cycle[1],"-",label="OWENS")
# PyPlot.legend()
# # PyPlot.ylim([0.0,0.8])
# # PyPlot.savefig("$(path)/../figs/34mCampbell.pdf",transparent = true)

# # Add to figure
# for i=1:2:numModes-6-2
#        PyPlot.plot(rpm,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
# end
# PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
# PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
# PyPlot.ylim([0,6.01])
# # PyPlot.savefig("$(path)/../figs/34mCampbellWGX.pdf",transparent = true)
