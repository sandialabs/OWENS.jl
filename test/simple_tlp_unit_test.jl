filepath = splitdir(@__FILE__)[1]
include("$(filepath)/../src/OWENS.jl")
# using OWENS
import HDF5

p = OWENS.hydro.tlp_platform(r_spar=2,draft=30, height=5, width=2, length=10, num=3, ofst=1)
p.make_mesh(mshRefFactor=1,
    # clcurv=360/200,
    write_file=false)
freq,T,Zi,rao = p.run_hydro()
#TODO: use numpy to avoid copying data

#write to file
filename = "$filepath/data/tpl_simple_unit_data.h5"
# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"freq_old",freq)
#     HDF5.write(file,"T_old",T)
#     HDF5.write(file,"Zi_old",Zi)
#     HDF5.write(file,"rao_old",rao)
# end

freq_old = HDF5.h5read(filename,"freq_old")
T_old = HDF5.h5read(filename,"T_old")
Zi_old = HDF5.h5read(filename,"Zi_old")
rao_old = HDF5.h5read(filename,"rao_old")

tol = 1e-6
@test isapprox(freq,freq_old,atol = tol)
@test isapprox(T,T_old,atol = tol)
@test isapprox(Zi,Zi_old,atol = tol)
@test isapprox(rao,rao_old,atol = tol)
