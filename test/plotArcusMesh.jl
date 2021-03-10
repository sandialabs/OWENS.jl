using PyPlot
close("all")

const module_path = splitdir(@__FILE__)[1]

include("$module_path/../src/OWENS.jl")

println("Start")
println("here")

mymesh, myort, myjoint = OWENS.create_arcus_mesh(;nblade=4,c_mount_ratio = 0.2)
figure()
plot3D(mymesh.x,mymesh.y,mymesh.z, "b.")
xlim([-maximum(mymesh.z)/2,maximum(mymesh.z)/2])
ylim([-maximum(mymesh.z)/2,maximum(mymesh.z)/2])
