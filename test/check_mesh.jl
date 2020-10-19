# using PyPlot
# close("all")
using Test
import OWENS
path,_ = splitdir(@__FILE__)
# include("$(path)/../src/OWENS.jl")
mesh = OWENS.readMesh("$(path)/data/unit_test_5MW.mesh")

mymesh = OWENS.create_mesh() #use defaults

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
