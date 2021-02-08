import PyPlot
PyPlot.close("all")
using Test
import DelimitedFiles
# import OWENS
path,_ = splitdir(@__FILE__)
include("$(path)/../src/OWENS.jl")
mesh = OWENS.readMesh("$(path)/data/unit_test_5MW.mesh")
joint = DelimitedFiles.readdlm("$(path)/data/unit_test_5MW.jnt",'\t',skipstart = 0)

plots = false

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

if plots

    PyPlot.figure()
    PyPlot.plot(mymesh.x,mymesh.z,"k.-",markersize=10.0)
    PyPlot.plot(mesh.x,mesh.z,"b.")
    PyPlot.axis("equal")

    PyPlot.figure()
    PyPlot.plot(mymesh.y,mymesh.z,"k.-",markersize=10.0)
    PyPlot.plot(mesh.y,mesh.z,"b.")
    PyPlot.legend(["mymesh","originalmesh"])
    # PyPlot.axis("equal")

    PyPlot.figure()
    PyPlot.plot(LinRange(0,1,length(mymesh.z)),mesh.z[1:length(mymesh.z)]-mymesh.z,"k.-")
    PyPlot.ylabel("z")

    PyPlot.figure()
    PyPlot.plot(LinRange(0,1,length(mymesh.x)),mesh.x[1:length(mymesh.x)]-mymesh.x,"k.-")
    PyPlot.ylabel("x")

    PyPlot.figure()
    PyPlot.plot(LinRange(0,1,length(mymesh.y)),mesh.y[1:length(mymesh.y)]-mymesh.y,"k.-")
    PyPlot.ylabel("y")

    PyPlot.figure()
    PyPlot.plot(LinRange(0,1,length(mymesh.conn[:,2])),mesh.conn[:,2].-mymesh.conn[:,2],"k.-")
    PyPlot.ylabel("connection")
end


##################################
### Blade Sectional Properties ###
##################################


bladeData,_,_,_ = OWENS.readBladeData("$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.bld") #reads overall blade data file
el = OWENS.readElementData(mymesh.numEl,"$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.el","$(path)/data/input_files_test/_15mTower_transient_dvawt_c_2_lcdt.ort",bladeData) #read element data file (also reads orientation and blade data file associated with elements)

NuMad_props_xlscsv_file = "$path/data/NuMAD_Props_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
xlsprops_blade = DelimitedFiles.readdlm(NuMad_props_xlscsv_file,',',skipstart = 0)

span_loc = Float64.(xlsprops_blade[3:23,1])
chord = Float64.(xlsprops_blade[3:23,2])
tw_aero = Float64.(xlsprops_blade[3:23,3])
ei_flap = Float64.(xlsprops_blade[3:23,4])
ei_lag = Float64.(xlsprops_blade[3:23,5])
gj = Float64.(xlsprops_blade[3:23,6])
ea = Float64.(xlsprops_blade[3:23,7])
s_fl = Float64.(xlsprops_blade[3:23,8])
s_af = Float64.(xlsprops_blade[3:23,9])
s_al = Float64.(xlsprops_blade[3:23,10])
s_ft = Float64.(xlsprops_blade[3:23,11])
s_lt = Float64.(xlsprops_blade[3:23,12])
s_at = Float64.(xlsprops_blade[3:23,13])
x_sc = Float64.(xlsprops_blade[3:23,14])
y_sc = Float64.(xlsprops_blade[3:23,15])
x_tc = Float64.(xlsprops_blade[3:23,16])
y_tc = Float64.(xlsprops_blade[3:23,17])
mass = Float64.(xlsprops_blade[3:23,18])
flap_iner = Float64.(xlsprops_blade[3:23,19])
lag_iner = Float64.(xlsprops_blade[3:23,20])
tw_iner = Float64.(xlsprops_blade[3:23,21])
x_cm = Float64.(xlsprops_blade[3:23,22])
y_cm = Float64.(xlsprops_blade[3:23,23])

#Blades

# Geometry
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

# Materials
NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# Precomp Outputs Restructured into owens format
bld_precompoutput,bld_precompinput = OWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(span_loc,numadIn_bld,bld_precompoutput)

if plots
    PyPlot.figure()
    toplot0 = zeros(length(span_loc))
    for ii = 1:length(span_loc)-1
        toplot0[ii] = sectionPropsArray_bld[ii].rhoA[1]
    end
    toplot0[end] = sectionPropsArray_bld[end].rhoA[2]

    PyPlot.plot(span_loc,toplot0,"r.-")
    PyPlot.plot(span_loc,mass,"b.-")
end


for ii = 1:length(sectionPropsArray_bld)
    tol = 0.05 # 5% error
    tol2 = 0.16 # %
    if !(isapprox(tw_aero[ii],sectionPropsArray_bld[ii].twist[1],atol=tol))
        println("tw_aero failed at $tol tolerance $(isapprox(tw_aero[ii],sectionPropsArray_bld[ii].twist[1],atol=tol)), increasing tolerance")
        @test isapprox(tw_aero[ii],sectionPropsArray_bld[ii].twist[1],atol=tol2)
    end
    if !(abs((ei_flap[ii]-sectionPropsArray_bld[ii].EIyy[1])/ei_flap[ii])<tol)
        println("ei_flap failed at $tol tolerance $(abs((ei_flap[ii]-sectionPropsArray_bld[ii].EIyy[1])/ei_flap[ii])), increasing tolerance")
        @test abs((ei_flap[ii]-sectionPropsArray_bld[ii].EIyy[1])/ei_flap[ii])<tol2
    end
    if !(abs((ei_lag[ii]-sectionPropsArray_bld[ii].EIzz[1])/ei_lag[ii])<tol)
        println("ei_lag failed at $tol tolerance $(abs((ei_lag[ii]-sectionPropsArray_bld[ii].EIzz[1])/ei_lag[ii])), increasing tolerance")
        @test abs((ei_lag[ii]-sectionPropsArray_bld[ii].EIzz[1])/ei_lag[ii])<tol2
    end
    if !(abs((gj[ii]-sectionPropsArray_bld[ii].GJ[1])/gj[ii])<tol)
        println("gj failed at $tol tolerance $(abs((gj[ii]-sectionPropsArray_bld[ii].GJ[1])/gj[ii])), increasing tolerance")
        @test abs((gj[ii]-sectionPropsArray_bld[ii].GJ[1])/gj[ii])<tol2
    end
    if !(abs((ea[ii]-sectionPropsArray_bld[ii].EA[1])/ea[ii])<tol)
        println("ea failed at $tol tolerance $(abs((ea[ii]-sectionPropsArray_bld[ii].EA[1])/ea[ii])), increasing tolerance")
        @test abs((ea[ii]-sectionPropsArray_bld[ii].EA[1])/ea[ii])<tol2
    end
    if !(abs((mass[ii]-sectionPropsArray_bld[ii].rhoA[1])/mass[ii])<tol)
        println("mass failed at $tol tolerance $(abs((mass[ii]-sectionPropsArray_bld[ii].rhoA[1])/mass[ii])), increasing tolerance")
        @test abs((mass[ii]-sectionPropsArray_bld[ii].rhoA[1])/mass[ii])<tol2
    end
    if !(abs((flap_iner[ii]-sectionPropsArray_bld[ii].rhoIyy[1])/flap_iner[ii])<tol)
        println("flap_iner failed at $tol tolerance $(abs((flap_iner[ii]-sectionPropsArray_bld[ii].rhoIyy[1])/flap_iner[ii])), increasing tolerance")
        @test abs((flap_iner[ii]-sectionPropsArray_bld[ii].rhoIyy[1])/flap_iner[ii])<tol2
    end
    if !(abs((lag_iner[ii]-sectionPropsArray_bld[ii].rhoIzz[1])/lag_iner[ii])<tol)
        println("lag_iner failed at $tol tolerance $(abs((lag_iner[ii]-sectionPropsArray_bld[ii].rhoIzz[1])/lag_iner[ii])), increasing tolerance")
        @test abs((lag_iner[ii]-sectionPropsArray_bld[ii].rhoIzz[1])/lag_iner[ii])<tol2
    end
    if !(abs((x_cm[ii]-sectionPropsArray_bld[ii].zcm[1])/x_cm[ii])<tol)
        println("x_cm failed at $tol tolerance $(abs((x_cm[ii]-sectionPropsArray_bld[ii].zcm[1])/x_cm[ii])), increasing tolerance")
        @test abs((x_cm[ii]-sectionPropsArray_bld[ii].zcm[1])/x_cm[ii])<tol2*2
    end
    if !(abs((y_cm[ii]-sectionPropsArray_bld[ii].ycm[1])/y_cm[ii])<tol)
        println("y_cm failed at $tol tolerance $(abs((y_cm[ii]-sectionPropsArray_bld[ii].ycm[1])/y_cm[ii])), increasing tolerance")
        @test abs((y_cm[ii]-sectionPropsArray_bld[ii].ycm[1])/y_cm[ii])<tol2
    end

end
