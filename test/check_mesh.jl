# using PyPlot
# close("all")
using Test
path,_ = splitdir(@__FILE__)
include("$(path)/../old_src/readMesh.jl")
include("$(path)/../old_src/create_mesh.jl")
mesh = readMesh("$(path)/../old_src/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.mesh")

mymesh = create_mesh() #use defaults

# Test
tol = 1e-6

@test isapprox(mesh.nodeNum,mymesh.nodeNum;atol=tol)
@test isapprox(mesh.numEl,mymesh.numEl;atol=tol)
@test isapprox(mesh.numNodes,mymesh.numNodes;atol=tol)
for i = 1:length(mesh.x)
    @test isapprox(mesh.x[i],mymesh.x[i];atol=3.0) # It looks like the original blade shape of the SNL5MW VAWT was done by hand...
    @test isapprox(mesh.y[i],mymesh.y[i];atol=0.3)
    @test isapprox(mesh.z[i],mymesh.z[i];atol=0.3)
end
@test isapprox(mesh.elNum,mymesh.elNum;atol=tol)
@test isapprox(mesh.conn,mymesh.conn;atol=tol)

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
