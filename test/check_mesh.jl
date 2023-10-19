using Test
import DelimitedFiles
import GyricFEA
import OWENS
import FLOWMath
path,_ = splitdir(@__FILE__)
# include("$path/../src/OWENS.jl")

mesh = OWENS.readMesh("$(path)/data/unit_test_5MW.mesh")
joint = DelimitedFiles.readdlm("$(path)/data/unit_test_5MW.jnt",'\t',skipstart = 0)

# Use the SNL5MW as the baseline check
#The 15 is subtracted off at the end of the line
SNL5MW_bld_z = [15.0, 21.61004296, 28.20951408, 28.2148, 34.81955704, 41.4296, 48.03964296, 54.63911408, 61.24915704, 67.8592, 74.46924296, 81.06871408, 87.67875704, 94.2888, 100.89884296, 107.49831408, 114.10835704, 120.7184, 127.32844296, 133.92791408, 133.9332, 140.53795704, 147.148].-15.0
SNL5MW_bld_x = -[0.0, -10.201, -20.361, -20.368290684, -29.478, -36.575, -42.579, -47.177, -50.555, -52.809, -53.953, -54.014, -53.031, -51.024, -47.979, -43.942, -38.768, -32.91, -25.587, -17.587, -17.580079568, -8.933, 8.0917312607e-15]

mymesh,myort,myjoint = OWENS.create_mesh_struts(;Ht=15.0,
Hb = 147.148-15.0, #blade height
R = 54.014, # m bade radius
nblade = 2,
ntelem = 20, #tower elements
nbelem = 20, #blade elements
nselem = 2,
strut_mountpointtop = 0.1,
strut_mountpointbot = 0.1,
bshapex = SNL5MW_bld_x, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
bshapez = SNL5MW_bld_z,
angularOffset = -pi/2) #Blade shape, magnitude is irrelevant, scaled based on height and radius above

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

# # Joints TODO: redo test since joints have been rewritten
# jointminormismatch = 0
# for i = 1:length(joint[:,1])
#     for j = 1:length(joint[1,:])
#         if isapprox(abs(joint[i,j]),180.0,atol=1.0) && isapprox(abs(myjoint[i,j]),180.0,atol=1.0) #180 and -180 are the same
#             @test isapprox(abs(joint[i,j]),abs(myjoint[i,j]);atol=1.1)
#         elseif isapprox(abs(joint[i,j]),abs(myjoint[i,j]);atol=1.1)
#             @test isapprox(abs(joint[i,j]),abs(myjoint[i,j]);atol=1.1) #Within 0.5 degrees since the 5MW blade shape was done by hand
#         else
#             global jointminormismatch += 1
#         end
#     end
# end

# @test jointminormismatch<7

# import PyPlot
# PyPlot.close("all")
# PyPlot.pygui(true)

# PyPlot.figure()
# for idot = 1:length(mymesh.x)

#     PyPlot.plot(mymesh.x[idot],mymesh.z[idot],"ko-",markersize=10.0)
#     PyPlot.plot(mesh.x[idot],mesh.z[idot],"b.-")
#     PyPlot.axis("equal")
#     sleep(0.2)
# end

# PyPlot.figure()
# PyPlot.plot(mymesh.y,mymesh.z,"k.-",markersize=10.0)
# PyPlot.plot(mesh.y,mesh.z,"b.")
# PyPlot.legend(["mymesh","originalmesh"])
# # PyPlot.axis("equal")

# PyPlot.figure()
# PyPlot.plot(LinRange(0,1,length(mymesh.z)),mesh.z[1:length(mymesh.z)]-mymesh.z,"k.-")
# PyPlot.ylabel("z")

# PyPlot.figure()
# PyPlot.plot(LinRange(0,1,length(mymesh.x)),mesh.x[1:length(mymesh.x)]-mymesh.x,"k.-")
# PyPlot.ylabel("x")

# PyPlot.figure()
# PyPlot.plot(LinRange(0,1,length(mymesh.y)),mesh.y[1:length(mymesh.y)]-mymesh.y,"k.-")
# PyPlot.ylabel("y")

# PyPlot.figure()
# PyPlot.plot(LinRange(0,1,length(mymesh.conn[:,2])),mesh.conn[:,2].-mymesh.conn[:,2],"k.-")
# PyPlot.ylabel("connection")

##################################
### Blade Sectional Properties ###
##################################


bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers = OWENS.readBladeData("$(path)/data/_15mTower_transient_dvawt_c_2_lcdt.bld") #reads overall blade data file
# @test isapprox(mymesh.structuralSpanLocNorm,structuralSpanLocNorm,atol=tol) #TODO: figure out how to resolve this since the old method used span length while interpolating on height (has since been updated to be consistent by using only height for both)
@test isapprox(mymesh.structuralNodeNumbers,structuralNodeNumbers,atol=tol)
@test isapprox(mymesh.structuralElNumbers[1:end-1],structuralElNumbers[1:end-1],atol=tol)

el = OWENS.readElementData(mymesh.numEl,"$(path)/data/_15mTower_transient_dvawt_c_2_lcdt.el","$(path)/data/_15mTower_transient_dvawt_c_2_lcdt.ort",bladeData) #read element data file (also reads orientation and blade data file associated with elements)

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

#Tower
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_TaperedTower.csv"
numadIn = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

#Add the full path
for (i,airfoil) in enumerate(numadIn.airfoil)
    numadIn.airfoil[i] = "$path/data/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_TaperedTower.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = OWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(mymesh.z[1:24],numadIn,precompoutput)

#Blades

# Geometry
NuMad_geom_xlscsv_file = "$path/data/NuMAD_Geom_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

#Add the full path
for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/data/airfoils/$airfoil"
end

# Materials
NuMad_mat_xlscsv_file = "$path/data/NuMAD_Materials_SNL_5MW_D_Carbon_LCDT_ThickFoils_ThinSkin.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

# Precomp Outputs Restructured into owens format
bld_precompoutput,bld_precompinput = OWENS.getPreCompOutput(numadIn_bld;plyprops)
newspan = mymesh.z[25:46].-15.0
newspan = newspan./maximum(newspan)
sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(newspan,numadIn_bld,bld_precompoutput) #TODO: why is this not aligning?

newspan2 = (newspan[1:end-1]+newspan[2:end])/2
newspan2 = newspan2.-newspan2[1]

tw_aero_new = [sectionPropsArray_bld[ii].twist[1] for ii = 1:length(sectionPropsArray_bld)]
tw_aero_spl = FLOWMath.akima(newspan2,tw_aero_new,span_loc)

ei_flap_new = [sectionPropsArray_bld[ii].EIyy[1] for ii = 1:length(sectionPropsArray_bld)]
ei_flap_spl = FLOWMath.akima(newspan2,ei_flap_new,span_loc)

ei_lag_new = [sectionPropsArray_bld[ii].EIzz[1] for ii = 1:length(sectionPropsArray_bld)]
ei_lag_spl = FLOWMath.akima(newspan2,ei_lag_new,span_loc)

gj_new = [sectionPropsArray_bld[ii].GJ[1] for ii = 1:length(sectionPropsArray_bld)]
gj_spl = FLOWMath.akima(newspan2,gj_new,span_loc)

ea_new = [sectionPropsArray_bld[ii].EA[1] for ii = 1:length(sectionPropsArray_bld)]
ea_spl = FLOWMath.akima(newspan2,ea_new,span_loc)

mass_new = [sectionPropsArray_bld[ii].rhoA[1] for ii = 1:length(sectionPropsArray_bld)]
mass_spl = FLOWMath.akima(newspan2,mass_new,span_loc)

flap_iner_new = [sectionPropsArray_bld[ii].rhoIyy[1] for ii = 1:length(sectionPropsArray_bld)]
flap_iner_spl = FLOWMath.akima(newspan2,flap_iner_new,span_loc)

lag_iner_new = [sectionPropsArray_bld[ii].rhoIzz[1] for ii = 1:length(sectionPropsArray_bld)]
lag_iner_spl = FLOWMath.akima(newspan2,lag_iner_new,span_loc)

x_cm_new = [sectionPropsArray_bld[ii].zcm[1] for ii = 1:length(sectionPropsArray_bld)]
x_cm_spl = FLOWMath.akima(newspan2,x_cm_new,span_loc)

y_cm_new = [sectionPropsArray_bld[ii].ycm[1] for ii = 1:length(sectionPropsArray_bld)]
y_cm_spl = FLOWMath.akima(newspan2,y_cm_new,span_loc)


for ii = 1:length(sectionPropsArray_bld)
    tol = 0.03 # 5% error
    tol2 = 0.07 # %
    # println(ii)
    if !(isapprox(tw_aero[ii],tw_aero_spl[ii],atol=tol))
        println("tw_aero failed at $(tol*100)% tolerance, error is $(abs((tw_aero[ii]-tw_aero_spl[ii])/tw_aero[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test isapprox(tw_aero[ii],tw_aero_spl[ii],atol=tol2)
    end
    if !(abs((ei_flap[ii]-ei_flap_spl[ii])/ei_flap[ii])<tol)
        println("ei_flap failed at $(tol*100)% tolerance, error is $(abs((ei_flap[ii]-ei_flap_spl[ii])/ei_flap[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((ei_flap[ii]-ei_flap_spl[ii])/ei_flap[ii])<tol2
    end
    if !(abs((ei_lag[ii]-ei_lag_spl[ii])/ei_lag[ii])<tol)
        println("ei_lag failed at $(tol*100)% tolerance, error is $(abs((ei_lag[ii]-ei_lag_spl[ii])/ei_lag[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((ei_lag[ii]-ei_lag_spl[ii])/ei_lag[ii])<tol2
    end
    # if !(abs((gj[ii]-gj_spl[ii])/gj[ii])<tol)
    #     println("gj failed at $(tol*100)% tolerance, error is $(abs((gj[ii]-gj_spl[ii])/gj[ii])*100)%, increasing tolerance to $(tol2*100) %")
    #     @test abs((gj[ii]-gj_spl[ii])/gj[ii])<tol2
    # end
    if !(abs((ea[ii]-ea_spl[ii])/ea[ii])<tol)
        println("ea failed at $(tol*100)% tolerance, error is $(abs((ea[ii]-ea_spl[ii])/ea[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((ea[ii]-ea_spl[ii])/ea[ii])<tol2
    end
    if !(abs((mass[ii]-mass_spl[ii])/mass[ii])<tol)
        println("mass failed at $(tol*100)% tolerance, error is $(abs((mass[ii]-mass_spl[ii])/mass[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((mass[ii]-mass_spl[ii])/mass[ii])<tol2
    end
    if !(abs((flap_iner[ii]-flap_iner_spl[ii])/flap_iner[ii])<tol)
        println("flap_iner failed at $(tol*100)% tolerance, error is $(abs((flap_iner[ii]-flap_iner_spl[ii])/flap_iner[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((flap_iner[ii]-flap_iner_spl[ii])/flap_iner[ii])<tol2
    end
    if !(abs((lag_iner[ii]-lag_iner_spl[ii])/lag_iner[ii])<tol)
        println("lag_iner failed at $(tol*100)% tolerance, error is $(abs((lag_iner[ii]-lag_iner_spl[ii])/lag_iner[ii])*100)%, increasing tolerance to $(tol2*100) %")
        @test abs((lag_iner[ii]-lag_iner_spl[ii])/lag_iner[ii])<tol2
    end
    # if !(abs((x_cm[ii]-x_cm_spl[ii])/x_cm[ii])<tol)
    #     println("x_cm failed at $(tol*100)% tolerance, error is $(abs((x_cm[ii]-x_cm_spl[ii])/x_cm[ii])*100)%, increasing tolerance to $(tol2*100) %")
    #     @test abs((x_cm[ii]-x_cm_spl[ii])/x_cm[ii])<tol2*2
    # end
    # if !(abs((y_cm[ii]-y_cm_spl[ii])/y_cm[ii])<tol)
    #     println("y_cm failed at $(tol*100)% tolerance, error is $(abs((y_cm[ii]-y_cm_spl[ii])/y_cm[ii])*100)%, increasing tolerance to $(tol2*100) %")
    #     @test abs((y_cm[ii]-y_cm_spl[ii])/y_cm[ii])<tol2
    # end

end


# PyPlot.figure()
# PyPlot.plot(span_loc, tw_aero,".-",label="tw_aero_old")
# PyPlot.plot(newspan2, tw_aero_new,".-",label="tw_aero_new")
# PyPlot.plot(span_loc, tw_aero_spl,".-",label="tw_aero_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, ei_flap,".-",label="ei_flap_old")
# PyPlot.plot(newspan2, ei_flap_new,".-",label="ei_flap_new")
# PyPlot.plot(span_loc, ei_flap_spl,".-",label="ei_flap_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, ei_lag,".-",label="ei_lag_old")
# PyPlot.plot(newspan2, ei_lag_new,".-",label="ei_lag_new")
# PyPlot.plot(span_loc, ei_lag_spl,".-",label="ei_lag_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, gj,".-",label="gj_old")
# PyPlot.plot(newspan2, gj_new,".-",label="gj_new")
# PyPlot.plot(span_loc, gj_spl,".-",label="gj_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, ea,".-",label="ea_old")
# PyPlot.plot(newspan2, ea_new,".-",label="ea_new")
# PyPlot.plot(span_loc, ea_spl,".-",label="ea_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, mass,".-",label="mass_old")
# PyPlot.plot(newspan2, mass_new,".-",label="mass_new")
# PyPlot.plot(span_loc, mass_spl,".-",label="mass_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, flap_iner,".-",label="flap_iner_old")
# PyPlot.plot(newspan2, flap_iner_new,".-",label="flap_iner_new")
# PyPlot.plot(span_loc, flap_iner_spl,".-",label="flap_iner_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, lag_iner,".-",label="lag_iner_old")
# PyPlot.plot(newspan2, lag_iner_new,".-",label="lag_iner_new")
# PyPlot.plot(span_loc, lag_iner_spl,".-",label="lag_iner_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, x_cm,".-",label="x_cm_old")
# PyPlot.plot(newspan2, x_cm_new,".-",label="x_cm_new")
# PyPlot.plot(span_loc, x_cm_spl,".-",label="x_cm_spl")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(span_loc, y_cm,".-",label="y_cm_old")
# PyPlot.plot(newspan2, y_cm_new,".-",label="y_cm_new")
# PyPlot.plot(span_loc, y_cm_spl,".-",label="y_cm_spl")
# PyPlot.legend()


# Check File Io
sectionPropsArray_bld2 = OWENS.getSectPropsFromPreComp(mymesh.z[25:47].-15.0,numadIn_bld,bld_precompoutput)

#Struts
# They are the same as the end properties of the blades

# Combined Section Props
Nremain = 8 #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;sectionPropsArray_bld2;sectionPropsArray_bld2; fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

filename = "$(path)/data/newmesh_5MW"
OWENS.saveOWENSfiles(filename,mymesh,myort,myjoint,myel,pBC,numadIn_bld)

mesh = OWENS.readMesh("$(path)/data/newmesh_5MW.mesh")
joint = DelimitedFiles.readdlm("$(path)/data/newmesh_5MW.jnt",'\t',skipstart = 0)

el = OWENS.readElementData(mymesh.numEl,"$(path)/data/newmesh_5MW.el","$(path)/data/newmesh_5MW.ort",bladeData) #read element data file (also reads orientation and blade data file associated with elements)


tol = 1e-4
@test isapprox(myel.elLen,el.elLen,atol=tol)
@test isapprox(myel.psi,el.psi,atol=tol)
@test isapprox(myel.theta,el.theta,atol=tol)
@test isapprox(myel.roll,el.roll,atol=tol)
@test isapprox(myel.rotationalEffects,el.rotationalEffects,atol=tol)

for jj = [1:23;25:length(myel.props)]
    # println(jj)
    # @test isapprox(myel.props[jj].ac,el.props[jj].ac,atol=tol) #TODO: figure this out
    @test isapprox(myel.props[jj].twist,el.props[jj].twist,atol=tol)
    @test isapprox(myel.props[jj].rhoA,el.props[jj].rhoA,atol=tol)
    @test isapprox(myel.props[jj].EIyy,el.props[jj].EIyy,atol=tol)
    @test isapprox(myel.props[jj].EIzz,el.props[jj].EIzz,atol=tol)
    @test isapprox(myel.props[jj].GJ,el.props[jj].GJ,atol=tol)
    @test isapprox(myel.props[jj].EA,el.props[jj].EA,atol=tol)
    @test isapprox(myel.props[jj].rhoIyy,el.props[jj].rhoIyy,atol=tol)
    @test isapprox(myel.props[jj].rhoIzz,el.props[jj].rhoIzz,atol=tol)
    @test isapprox(myel.props[jj].rhoJ,el.props[jj].rhoJ,atol=tol)
    @test isapprox(myel.props[jj].zcm,el.props[jj].zcm,atol=tol)
    @test isapprox(myel.props[jj].ycm,el.props[jj].ycm,atol=tol)
    # @test isapprox(myel.props[jj].a,el.props[jj].a,atol=tol) # TODO: see below
    @test isapprox(myel.props[jj].EIyz,el.props[jj].EIyz,atol=tol)
    @test isapprox(myel.props[jj].alpha1,el.props[jj].alpha1,atol=tol)
    @test isapprox(myel.props[jj].alpha2,el.props[jj].alpha2,atol=tol)
    @test isapprox(myel.props[jj].alpha3,el.props[jj].alpha3,atol=tol)
    @test isapprox(myel.props[jj].alpha4,el.props[jj].alpha4,atol=tol)
    @test isapprox(myel.props[jj].alpha5,el.props[jj].alpha5,atol=tol)
    @test isapprox(myel.props[jj].alpha6,el.props[jj].alpha6,atol=tol)
    @test isapprox(myel.props[jj].rhoIyz,el.props[jj].rhoIyz,atol=tol)
    # @test isapprox(myel.props[jj].b,el.props[jj].b,atol=tol) #TODO: the original el reader isn't using the correct column for chord, and the bld file seems to not be using the correct chord...
    @test isapprox(myel.props[jj].a0,el.props[jj].a0,atol=tol)
    # @test isapprox(myel.props[jj].aeroCenterOffset,el.props[jj].aeroCenterOffset,atol=tol) #TODO: same here
end
