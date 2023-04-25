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
println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = ModelGen.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld)


starttime2 = time()

prescribed_conditions = Dict(
# fixed base
1 => PrescribedConditions(ux=0, uy=0, uz=0,theta_x=0, theta_y=0, theta_z=0),
# fixed top
top_idx => PrescribedConditions(ux=0, uy=0, uz=0),
)

gravity = [0, 0, -9.81]

rotorRPM = collect(0.0:5:40.0)

freq2 = zeros(length(rotorRPM),numModes)
λ = []
eigenstates =[]
for (i,rpm) in enumerate(rotorRPM)

    w0 = [0, 0, rpm*(2*pi)/60]

    # eigenvalues and (right) eigenvectors
    global system, Up, λ, eigenstates
    nev = numModes*2

    system, λ, V, converged = eigenvalue_analysis!(system, assembly;
    prescribed_conditions = prescribed_conditions,
    angular_velocity = w0,
    linear = false,
    reset_state = true,
    gravity = gravity,
    nev)

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
#        PyPlot.plot(rotorRPM,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
# end
# PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
# PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
# PyPlot.ylim([0,6.01])
# # PyPlot.savefig("$(path)/../figs/34mCampbellWGX.pdf",transparent = true)
