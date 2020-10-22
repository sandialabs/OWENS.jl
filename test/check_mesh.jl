# using PyPlot
# close("all")
using Test
import DelimitedFiles
import OWENS
path,_ = splitdir(@__FILE__)
# include("$(path)/../src/OWENS.jl")
mesh = OWENS.readMesh("$(path)/data/unit_test_5MW.mesh")
joint = DelimitedFiles.readdlm("$(path)/data/unit_test_5MW.jnt",'\t',skipstart = 0)

#TODO: ort file, nodal file, and blade file

# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL5MW_bld_z = [15.0, 21.61004296, 28.20951408, 28.2148, 34.81955704, 41.4296, 48.03964296, 54.63911408, 61.24915704, 67.8592, 74.46924296, 81.06871408, 87.67875704, 94.2888, 100.89884296, 107.49831408, 114.10835704, 120.7184, 127.32844296, 133.92791408, 133.9332, 140.53795704, 147.148].-15.0
SNL5MW_bld_x = -[0.0, -10.201, -20.361, -20.368290684, -29.478, -36.575, -42.579, -47.177, -50.555, -52.809, -53.953, -54.014, -53.031, -51.024, -47.979, -43.942, -38.768, -32.91, -25.587, -17.587, -17.580079568, -8.933, 8.0917312607e-15]

mymesh,myort,myjoint = OWENS.create_mesh(;Ht = 15.0, #tower height before blades attach
                                    Hb = 147.148-15.0, #blade height
                                    R = 54.014, # m bade radius
                                    nstrut = 2,
                                    strut_mout_ratio = 0.1, #distance from top/bottom
                                    ntelem = 20, #tower elements
                                    nbelem = 20, #blade elements
                                    nselem = 2,  #strut elements
                                    bshapex=SNL5MW_bld_x,
                                    bshapez=SNL5MW_bld_z) #use defaults


# println("Element Point Psi_d Theta_d")
# for i = 1:length(el.elNum)
#     println("$(i) $(el.elNum[i]) $(el.Psi_d[i]) $(el.Theta_d[i])")
# end
#
# for i = 1:length(mesh.conn[:,1])
#     println(mesh.conn[i,2])
# end

# Test
tol = 0.1

# Mesh
@test isapprox(mesh.nodeNum,mymesh.nodeNum;atol=tol)
@test isapprox(mesh.numEl,mymesh.numEl;atol=tol)
@test isapprox(mesh.numNodes,mymesh.numNodes;atol=tol)
for i = 1:length(mesh.x)
    @test isapprox(mesh.x[i],mymesh.x[i];atol=tol) # It looks like the original blade shape of the SNL5MW VAWT was done by hand...
    @test isapprox(mesh.y[i],mymesh.y[i];atol=tol)
    @test isapprox(mesh.z[i],mymesh.z[i];atol=tol)
end
@test isapprox(mesh.elNum,mymesh.elNum;atol=tol)
for i = 1:length(mesh.conn[:,1])
    # println(i)
    # println("$(mesh.conn[i,:]) $(mymesh.conn[i,:]) $(isapprox(mesh.conn[i,:],mymesh.conn[i,:]))")
    @test isapprox(mesh.conn[i,:],mymesh.conn[i,:];atol=tol)
end

# Joints
for i = 1:length(joint[:,1])
    for j = 1:length(joint[1,:])
        if isapprox(abs(joint[i,j]),180.0,atol=1.0) && isapprox(abs(myjoint[i,j]),180.0,atol=1.0) #180 and -180 are the same
            @test isapprox(abs(joint[i,j]),abs(myjoint[i,j]);atol=0.3)
        else
            @test isapprox(joint[i,j],myjoint[i,j];atol=0.3) #Within 0.5 degrees since the 5MW blade shape was done by hand
        end
    end
end

# figure()
# plot(mymesh.x,mymesh.z,"k.-",markersize=10.0)
# plot(mesh.x,mesh.z,"b.")
# axis("equal")
#
# figure()
# plot(LinRange(0,1,length(mymesh.z)),mesh.z[1:length(mymesh.z)]-mymesh.z,"k.-")
# ylabel("z")
#
# figure()
# plot(LinRange(0,1,length(mymesh.z)),mesh.x[1:length(mymesh.x)]-mymesh.x,"k.-")
# ylabel("x")
#
# figure()
# plot(LinRange(0,1,length(mymesh.conn[:,2])),mesh.conn[:,2].-mymesh.conn[:,2],"k.-")
# ylabel("connection")
