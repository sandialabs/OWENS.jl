using PyPlot
close("all")

const module_path = splitdir(@__FILE__)[1]

include("$module_path/../src/OWENS.jl")

println("Start")
println("here")

mymesh, myort, myjoint = OWENS.create_arcus_mesh()
figure()
plot3D(mymesh.x,mymesh.y,mymesh.z, "b.")
