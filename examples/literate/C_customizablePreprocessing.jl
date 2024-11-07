# # [Customizable Preprocessing](@id simple3)
# 
# In this example, we show the third level of what is going on behind the precompiled binary
# This includes all of the second, but also breaks out the setupOWENS function.
# This would be a good starting point if you need to make modifications to use a 
# unique mesh generation function, change how sectional properties are input, or adapt for a unique design
# an properly map the sectional properties to each element and apply unique boundary conditions, etc.
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook todo: get link working:
#-

# This example is the same as example B except that the setupOWENS function is broken out and each step defined

import OWENS
import OWENSFEA
import OWENSAero
import QuadGK
import FLOWMath
import PyPlot
#### PyPlot.pygui(true)
import OWENSOpenFASTWrappers
 

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
## runpath = path = splitdir(@__FILE__)[1]

Inp = OWENS.MasterInput("$runpath/sampleOWENS.yml")

nothing

# Unpack inputs, or you could directly input them here and bypass the file 

verbosity = 1

analysisType = Inp.analysisType
turbineType = Inp.turbineType
eta = Inp.eta
Nbld = Inp.Nbld
towerHeight = Inp.towerHeight
rho = Inp.rho
Vinf = Inp.Vinf
controlStrategy = Inp.controlStrategy
RPM = Inp.RPM
Nslices = Inp.Nslices
ntheta = Inp.ntheta
structuralModel = Inp.structuralModel
ntelem = Inp.ntelem
nbelem = Inp.nbelem
ncelem = Inp.ncelem
nselem = Inp.nselem
ifw = Inp.ifw
WindType = Inp.WindType
AeroModel = Inp.AeroModel
windINPfilename = "$(path)$(Inp.windINPfilename)"
ifw_libfile = Inp.ifw_libfile
if ifw_libfile == "nothing"
    ifw_libfile = nothing
end
Blade_Height = Inp.Blade_Height
Blade_Radius = Inp.Blade_Radius
numTS = Inp.numTS
delta_t = Inp.delta_t
NuMad_geom_xlscsv_file_twr = "$(path)$(Inp.NuMad_geom_xlscsv_file_twr)"
NuMad_mat_xlscsv_file_twr = "$(path)$(Inp.NuMad_mat_xlscsv_file_twr)"
NuMad_geom_xlscsv_file_bld = "$(path)$(Inp.NuMad_geom_xlscsv_file_bld)"
NuMad_mat_xlscsv_file_bld = "$(path)$(Inp.NuMad_mat_xlscsv_file_bld)"
NuMad_geom_xlscsv_file_strut = "$(path)$(Inp.NuMad_geom_xlscsv_file_strut)"
NuMad_mat_xlscsv_file_strut = "$(path)$(Inp.NuMad_mat_xlscsv_file_strut)"
adi_lib = Inp.adi_lib
if adi_lib == "nothing"
    adi_lib = nothing
end
adi_rootname = "$(path)$(Inp.adi_rootname)"

println("Set up Turbine")

AD15On = true
B = Nbld
R = Blade_Radius#177.2022*0.3048 #m
H = Blade_Height#1.02*R*2 #m

shapeZ = collect(LinRange(0,H,Nslices+1))
shapeX = R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)
shapeY = zeros(Nslices+1)


stack_layers_bld = nothing
stack_layers_scale = [1.0,1.0]
chord_scale = [1.0,1.0]
thickness_scale = [1.0,1.0]
Htwr_base=towerHeight
strut_mountpointbot = 0.11
strut_mountpointtop = 0.89
joint_type = 0
c_mount_ratio = 0.05
angularOffset = -pi/2
custommesh = nothing
if AeroModel=="AD" #TODO: unify flag
    AD15On=true #AD for AeroDyn, DMS for double multiple streamtube, AC for actuator cylinder
else
    AD15On=false
end
DynamicStallModel="BV"
RPI=true
cables_connected_to_blade_base = true
meshtype = "Darrieus"

nothing

# Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
# We do some intermediate calculations on the blade shape and angles

Nstrutperbld = 2 #TODO: generalize and propogate
strut_twr_mountpoint = [0.25,0.75]
strut_bld_mountpoint = [0.25,0.75]
Nbld = B
H = maximum(shapeZ) #m,
Htwr_blds = H
AD15hubR = 0.1
R = maximum(shapeX) #m,
omega = RPM / 60 * 2 * pi
tsr = omega*R/Vinf

nothing

# Here we set up the mesh using one of the pre-made meshing functions.  For this case, there is a function for the ARCUS, as 
# well as for towered VAWTs where you can have an arbitrary blade shape with connected struts, and if the blade tips touch
# the tower, then you can tell it to constrain them to the tower thus allowing for both H-VAWT and Darrieus designs.

#########################################
### Set up mesh
#########################################
if meshtype == "ARCUS" && custommesh == nothing #TODO, for all of these propogate the AeroDyn additional output requirements
    mymesh,myort,myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.create_arcus_mesh(;Htwr_base,
        Hbld = H, #blade height
        R, # m bade radius
        nblade = Nbld,
        ntelem, #tower elements
        nbelem, #blade elements
        ncelem,
        c_mount_ratio,
        bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        bshapez = shapeZ,
        AD15_ccw = true,
        joint_type, #hinged about y axis
        cables_connected_to_blade_base,
        angularOffset) #Blade shape, magnitude is irrelevant, scaled based on height and radius above
elseif (meshtype == "Darrieus" || meshtype == "H-VAWT") && custommesh == nothing
    
    if meshtype == "Darrieus"
        connectBldTips2Twr = true
    else
        connectBldTips2Twr = false
    end

    mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.create_mesh_struts(;Htwr_base,
        Htwr_blds,
        Hbld = H, #blade height
        R, # m bade radius
        AD15hubR, #TODO: hook up with AD15 file generation
        nblade = Nbld,
        ntelem, #tower elements
        nbelem, #blade elements
        nselem,
        strut_twr_mountpoint,
        strut_bld_mountpoint,
        bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        bshapez = shapeZ,
        bshapey = shapeY, # but magnitude for this is relevant
        angularOffset, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        AD15_ccw = true,
        verbosity=0, # 0 nothing, 1 basic, 2 lots: amount of printed information
        connectBldTips2Twr)
elseif custommesh != nothing
    mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = custommesh(;Htwr_base,
    Htwr_blds,
    Hbld = H, #blade height
    R, # m bade radius
    AD15hubR, #TODO: hook up with AD15 file generation
    nblade = Nbld,
    ntelem, #tower elements
    nbelem, #blade elements
    nselem,
    strut_twr_mountpoint,
    strut_bld_mountpoint,
    bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = shapeZ,
    bshapey = shapeY, # but magnitude for this is relevant
    angularOffset, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = true,
    verbosity=0, # 0 nothing, 1 basic, 2 lots: amount of printed information)
    )
else #TODO unify with HAWT
    error("please choose a valid mesh type (Darrieus, H-VAWT, ARCUS)")
end

nTwrElem = Int(mymesh.meshSeg[1])
try
    if contains(NuMad_mat_xlscsv_file_bld,"34m") || meshtype == "ARCUS" #TODO: this is really odd, 
        nTwrElem = Int(mymesh.meshSeg[1])+1
    end
catch
    nTwrElem = Int(mymesh.meshSeg[1])
end

nothing

# Here is a way that you can visualize the nodal numbers of the mesh

PyPlot.figure()
for icon = 1:length(mymesh.conn[:,1])
    idx1 = mymesh.conn[icon,1]
    idx2 = mymesh.conn[icon,2]
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    #### sleep(0.1)
end

for ijoint = 1:length(myjoint[:,1])
    idx2 = Int(myjoint[ijoint,2])
    idx1 = Int(myjoint[ijoint,3])
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",color="r",ha="center",va="center")
    PyPlot.text3D(mymesh.x[idx2].+rand()/30,mymesh.y[idx2].+rand()/30,mymesh.z[idx2].+rand()/30,"$idx2",color="r",ha="center",va="center")
    #### sleep(0.1)
end
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")


# This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

#########################################
### Set up Sectional Properties
#########################################

if !isnothing(NuMad_geom_xlscsv_file_twr)
    numadIn_twr = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_twr;section=:tower)
else
    n_web = 0
    n_stack = 2
    n_segments = 2
    span = [0.0, 6.607421057, 13.21484211, 19.82226317, 26.42968423, 33.03710529, 39.64452634, 46.2519474, 52.85936846, 59.46678951, 66.07421057, 72.68163163, 79.28905268, 85.89647374, 92.5038948, 99.11131586, 105.7187369, 112.326158, 118.933579, 125.5410001, 132.1484211]
    airfoil = ["circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular"]
    te_type = ["round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round"]
    twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    chord = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.25, 0.25, 0.25]
    xoffset = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    aerocenter = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    stack_mat_types = [8, 2]
    stack_layers = [70 3; 70 3; 70 3; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03]
    segments = [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0]
    DPtypes = ["" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""]
    skin_seq = [Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2])]
    web_seq = Array{Seq, 2}(undef, length(twist_d),0) #can be any number of stack nums, so we have to make non-square containers
    web_dp = Array{Seq, 2}(undef, length(twist_d),0) #this is fixed size square, but it's easier to do it this way

    numadIn_twr = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
end

#### Add the full path
for (i,airfoil) in enumerate(numadIn_twr.airfoil)
    numadIn_twr.airfoil[i] = "$path/airfoils/$airfoil"
end

nothing

# Here is where the material properties for the tower are either read in from the file, or directly input

if !isnothing(NuMad_mat_xlscsv_file_twr)
    plyprops_twr = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
else
    names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
    plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e9, 1.0e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
    plyprops_twr = OWENS.plyproperties(names,plies)
end

# Then this is where precomp.jl is called to get first the precomp outputs, then formatting those into the OWENS format, and then in the GXBeam.jl format for if GXBeam is used as the structural solver.

twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = OWENS.getOWENSPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
sectionPropsArray_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
stiff_twr, mass_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

nothing

# For the blades, we repeat what was done for the tower, but also include some simple design options for scaling thicknesses,
if !isnothing(NuMad_geom_xlscsv_file_bld)
    numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_bld;section=:blade)
else
    n_web = 1
    n_stack = 7
    n_segments = 12
    span = [0.0, 6.607, 13.215, 19.822, 26.43, 33.037, 39.645, 46.252, 52.859, 59.467, 66.074, 72.682, 79.289, 85.896, 92.504, 99.111, 105.719, 112.326, 118.934, 125.541, 132.148]
    airfoil = ["circular", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021"]
    te_type = ["round", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp"]
    twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    chord = [10.0, 10.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    xoffset = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    aerocenter = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
    stack_mat_types = [8, 2, 4, 8, 8, 8, 4]
    if isnothing(stack_layers_bld)
        stack_layers = [30.0 2.0 15.0 25.0 25.0 2.0 13.0; 15.0 2.0 10.0 13.0 11.0 2.0 11.0; 10.0 1.0 8.0 10.0 10.0 2.0 10.0; 8.0 1.0 6.0 9.0 10.0 1.0 9.0; 7.0 1.0 5.0 8.0 9.0 1.0 7.0; 6.0 1.0 4.0 8.0 9.0 1.0 6.0; 6.0 1.0 4.0 8.0 8.0 1.0 5.0; 6.0 1.0 4.0 7.0 7.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 7.0 8.0 3.0 5.0; 7.0 2.0 3.0 9.0 12.0 3.0 6.0; 10.0 3.0 4.0 11.0 15.0 3.0 6.0; 12.0 3.0 4.0 13.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 10.0 1.0 4.0 10.0 12.0 1.0 5.0]
    else
        stack_layers = stack_layers_bld
    end
    segments = [-1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0]
    DPtypes = ["" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""]
    skin_seq = [Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2])]
    web_seq = [Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]);;]
    web_dp = [Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]);;]

    numadIn_bld = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
end
for icol = 1:length(numadIn_bld.stack_layers[1,:])
    numadIn_bld.stack_layers[:,icol] .*= LinRange(stack_layers_scale[1],stack_layers_scale[2],length(numadIn_bld.chord))
end
numadIn_bld.chord .*= LinRange(chord_scale[1],chord_scale[2],length(numadIn_bld.chord))

for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
end

if !isnothing(NuMad_mat_xlscsv_file_bld)
    plyprops_bld = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_bld)
else
    names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
    plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
    plyprops_bld = OWENS.plyproperties(names,plies)
end

bld1start = Int(mymesh.structuralNodeNumbers[1,1]) #Get blade spanwise position
bld1end = Int(mymesh.structuralNodeNumbers[1,end])
spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

if length(thickness_scale)==2
    yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_bld.span)))
elseif length(thickness_scale)==length(numadIn_bld.span)
    yscale = thickness_scale
end

bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getOWENSPreCompOutput(numadIn_bld;yscale,plyprops = plyprops_bld)
sectionPropsArray_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
stiff_bld, mass_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

nothing

# Similarly for the struts, we do what was done for the blades
numadIn_strut = Array{OWENS.NuMad, 1}(undef, Nstrutperbld)
strut_precompoutput = Array{Array{OWENS.OWENSPreComp.Output, 1},1}(undef, Nstrutperbld)
strut_precompinput = Array{Array{OWENS.OWENSPreComp.Input, 1},1}(undef, Nstrutperbld)
plyprops_strut = Array{OWENS.plyproperties, 1}(undef, Nstrutperbld)
lam_U_strut = Array{Array{OWENS.Composites.Laminate, 2}}(undef, Nstrutperbld)
lam_L_strut = Array{Array{OWENS.Composites.Laminate, 2}}(undef, Nstrutperbld)
lam_W_strut = Array{Array{OWENS.Composites.Laminate, 2}}(undef, Nstrutperbld)
sectionPropsArray_strut = Array{Array{OWENSFEA.SectionPropsArray, 1}}(undef, Nstrutperbld)
stiff_strut = Array{Array{Array{Float64,2}, 1},1}(undef, Nstrutperbld)
mass_strut = Array{Array{Array{Float64,2}, 1},1}(undef, Nstrutperbld)
for istrut = 1:Nstrutperbld
    if !isnothing(NuMad_geom_xlscsv_file_strut)
        if typeof(NuMad_geom_xlscsv_file_strut)==String
            numadIn_strut[istrut] = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_strut)
        elseif typeof(NuMad_geom_xlscsv_file_strut) == OrderedCollections.OrderedDict{Symbol, Any}
            if length(NuMad_geom_xlscsv_file_strut[:components][:struts])==1
                numadIn_strut[istrut] = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_strut;section=:struts,subsection=1)
            else
                numadIn_strut[istrut] = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_strut;section=:struts,subsection=istrut)
            end
        else
            numadIn_strut[istrut] = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_strut[istrut])
        end
    else
        n_web = 1
        n_stack = 7
        n_segments = 12
        span = [0.0, 6.607, 13.215, 19.822, 26.43, 33.037, 39.645, 46.252, 52.859, 59.467, 66.074, 72.682, 79.289, 85.896, 92.504, 99.111, 105.719, 112.326, 118.934, 125.541, 132.148]
        airfoil = ["circular", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021"]
        te_type = ["round", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp"]
        twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        chord = [10.0, 10.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        xoffset = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        aerocenter = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
        stack_mat_types = [8, 2, 4, 8, 8, 8, 4]
        if isnothing(stack_layers_strut)
            stack_layers = [30.0 2.0 15.0 25.0 25.0 2.0 13.0; 15.0 2.0 10.0 13.0 11.0 2.0 11.0; 10.0 1.0 8.0 10.0 10.0 2.0 10.0; 8.0 1.0 6.0 9.0 10.0 1.0 9.0; 7.0 1.0 5.0 8.0 9.0 1.0 7.0; 6.0 1.0 4.0 8.0 9.0 1.0 6.0; 6.0 1.0 4.0 8.0 8.0 1.0 5.0; 6.0 1.0 4.0 7.0 7.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 7.0 8.0 3.0 5.0; 7.0 2.0 3.0 9.0 12.0 3.0 6.0; 10.0 3.0 4.0 11.0 15.0 3.0 6.0; 12.0 3.0 4.0 13.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 10.0 1.0 4.0 10.0 12.0 1.0 5.0]
        else
            stack_layers = stack_layers_strut
        end
        segments = [-1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0]
        DPtypes = ["" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""]
        skin_seq = [Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2])]
        web_seq = [Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]);;]
        web_dp = [Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]);;]

        numadIn_strut[istrut] = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end
    for icol = 1:length(numadIn_strut[istrut].stack_layers[1,:])
        numadIn_strut[istrut].stack_layers[:,icol] .*= LinRange(stack_layers_scale[1],stack_layers_scale[2],length(numadIn_strut[istrut].chord))
    end
    #### numadIn_strut[istrut].chord .*= LinRange(chord_scale[1],chord_scale[2],length(numadIn_strut[istrut].chord))

    for (i,airfoil) in enumerate(numadIn_strut[istrut].airfoil)
        numadIn_strut[istrut].airfoil[i] = "$path/airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_strut)
        plyprops_strut[istrut] = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_strut)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_strut[istrut] = OWENS.plyproperties(names,plies)
    end

    #### TODO: not straight struts
    spanpos = LinRange(0,1,nselem+1)#[0.0;cumsum(sqrt.(diff(mymesh.x[strut1start:strut1end]).^2 .+ diff(mymesh.z[strut1start:strut1end]).^2))]

    if length(thickness_scale)==2
        yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_strut[istrut].span)))
    elseif length(thickness_scale)==length(numadIn_strut[istrut].span)
        yscale = thickness_scale
    end

    strut_precompoutput[istrut],strut_precompinput[istrut],lam_U_strut[istrut],lam_L_strut[istrut],lam_W_strut[istrut] = OWENS.getOWENSPreCompOutput(numadIn_strut[istrut];yscale,plyprops = plyprops_strut[istrut])
    sectionPropsArray_strut[istrut] = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_strut[istrut],strut_precompoutput[istrut];precompinputs=strut_precompinput[istrut])
    stiff_strut[istrut], mass_strut[istrut] = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_strut[istrut],strut_precompoutput[istrut];GX=true)
end

nothing

# Here we combine the section properties into an array matching the mesh elements
bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
#### strutssecprops = collect(Iterators.flatten(fill(sectionPropsArray_strut, Nstrutperbld*Nbld)))

if meshtype == "ARCUS"
    sectionPropsArray = [sectionPropsArray_twr; bldssecprops]#; strutssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

    stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
    stiff_array = [stiff_twr; stiff_blds]#; stiff_struts]

    mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
    mass_array = [mass_twr; mass_blds]#; mass_struts]

    for icable = 1:Nbld
        sectionPropsArray = [sectionPropsArray; sectionPropsArray_strut[1]]
        stiff_array = [stiff_array;stiff_strut[1]]
        mass_array = [mass_array;mass_strut[1]]
    end
else
    sectionPropsArray = [sectionPropsArray_twr; bldssecprops]#; strutssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

    stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
    stiff_array = [stiff_twr; stiff_blds]#; stiff_struts]

    mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
    mass_array = [mass_twr; mass_blds]#; mass_struts]

    for istrut = 1:Nstrutperbld
        for ibld = 1:Nbld
            global sectionPropsArray = [sectionPropsArray; sectionPropsArray_strut[istrut]]
            global stiff_array = [stiff_array;stiff_strut[istrut]]
            global mass_array = [mass_array;mass_strut[istrut]]
        end
    end
end
rotationalEffects = ones(mymesh.numEl) #TODO: non rotating tower, or rotating blades

if length(sectionPropsArray)<mymesh.numEl
    @warn "There are more mesh elements than sectional properties, applying the last strut's sectional properties to the remaining"
    n_diff = mymesh.numEl - length(sectionPropsArray)
    sectionPropsArray = [sectionPropsArray; fill(sectionPropsArray_strut[end][2],n_diff)]
    stiff_array = [stiff_array;fill(stiff_strut[end][2],n_diff)]
    mass_array = [mass_array;fill(mass_strut[end][2],n_diff)]
end

#### store data in element object
myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)
system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_array,mass_array)

nothing

# Set up the OWENSAero aerodynamics if used
if !AD15On
    #########################################
    ### Set up aero forces
    #########################################
    #### translate from blade span to blade height between the numad definition and the vertical slice positions
    #### First get the angles from the overall geometry npoints and go to the numad npoints
    delta_xs = shapeX[2:end] - shapeX[1:end-1]
    delta_zs = shapeZ[2:end] - shapeZ[1:end-1]
    delta3D = atan.(delta_xs./delta_zs)
    delta3D_spl = safeakima(shapeZ[1:end-1]./maximum(shapeZ[1:end-1]), delta3D,LinRange(0,1,length(numadIn_bld.span)-1))
    #### now convert the numad span to a height
    bld_height_numad = cumsum(diff(numadIn_bld.span).*(1.0.-abs.(sin.(delta3D_spl))))
    bld_height_numad_unit = bld_height_numad./maximum(bld_height_numad)
    #### now we can use it to access the numad data 
    chord = safeakima(bld_height_numad_unit, numadIn_bld.chord,LinRange(bld_height_numad_unit[1],1,Nslices))
    airfoils = fill("nothing",Nslices)

    twist = safeakima(bld_height_numad_unit, numadIn_bld.twist_d.*pi/180,LinRange(bld_height_numad_unit[1],1,Nslices))

    #### Discretely assign the airfoils
    for (iheight_numad,height_numad) in enumerate(bld_height_numad_unit)
        for (iheight,height_slices) in enumerate(collect(LinRange(0,1,Nslices)))
            if airfoils[iheight]=="nothing" && height_slices<=height_numad
                airfoils[iheight] = "$(numadIn_bld.airfoil[iheight_numad]).dat"
            end
        end
    end

    OWENSAero.setupTurb(shapeX,shapeZ,B,chord,tsr,Vinf;AeroModel,DynamicStallModel,
    afname = airfoils,
    bld_y = shapeY,
    rho,
    twist, #TODO: verify twist is in same direction
    mu,
    eta,
    ifw, #TODO: propogate WindType
    turbsim_filename = windINPfilename,
    ifw_libfile,
    tau = [1e-5,1e-5],
    Aero_AddedMass_Active,
    Aero_RotAccel_Active,
    Aero_Buoyancy_Active,
    ntheta,
    Nslices,
    RPI)

    aeroForcesACDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
    deformAeroACDMS = OWENSAero.deformTurb
end
nothing

# Set up AeroDyn if used
# Here we create AeroDyn the files, first by specifying the names, then by creating the files, TODO: hook up the direct sectionPropsArray_str
# Then by initializing AeroDyn and grabbing the backend functionality with a function handle
if AD15On
    ad_input_file = "$path/ADInputFile_SingleTurbine2.dat"
    ifw_input_file = "$path/IW2.dat"
    OLAF_filename = "$path/OLAF2.dat"

    NumADBldNds = NumADStrutNds = 10 

    bldchord_spl = OWENS.safeakima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,NumADBldNds))

    #### Discretely assign the blade airfoils based on the next closest neighbor
    bld_airfoil_filenames = fill("nothing",NumADBldNds) #TODO: cable drag?
    for (ispan_numad,span_numad) in enumerate(numadIn_bld.span./maximum(numadIn_bld.span))
        for (ispan,span_slices) in enumerate(collect(LinRange(0,1,NumADBldNds)))
            if bld_airfoil_filenames[ispan]=="nothing" && span_slices<=span_numad
                bld_airfoil_filenames[ispan] = "$(numadIn_bld.airfoil[ispan_numad]).dat"
            end
        end
    end
    
    if meshtype == "ARCUS" 
        blade_filenames = ["$path/blade$i.dat" for i=1:Nbld]
        blade_chords = [bldchord_spl for i=1:Nbld]
        blade_Nnodes = [NumADBldNds for i=1:Nbld]
        airfoil_filenames = [bld_airfoil_filenames for i=1:Nbld]
        
    else
        blade_filenames = ["$path/blade$i.dat" for i=1:Nbld]
        blade_chords = [bldchord_spl for i=1:Nbld]
        blade_Nnodes = [NumADBldNds for i=1:Nbld]
        airfoil_filenames = collect(Iterators.flatten([bld_airfoil_filenames for i=1:Nbld]))
        
        for istrut = 1:Nstrutperbld
            strutchord_spl = OWENS.safeakima(numadIn_strut[istrut].span./maximum(numadIn_strut[istrut].span), numadIn_strut[istrut].chord,LinRange(0,1,NumADStrutNds))
            for ibld = 1:Nbld
                global blade_filenames = [blade_filenames;"$path/strut$(istrut)_bld$ibld.dat"]
                global blade_chords = [blade_chords;[strutchord_spl]]
                global blade_Nnodes = [blade_Nnodes;NumADStrutNds]

                #### Discretely assign the strut airfoils based on the next closest neighbor
                strut_airfoil_filenames = fill("nothing",NumADStrutNds)
                for (ispan_numad,span_numad) in enumerate(numadIn_strut[istrut].span./maximum(numadIn_strut[istrut].span))
                    for (ispan,span_slices) in enumerate(collect(LinRange(0,1,NumADBldNds)))
                        if strut_airfoil_filenames[ispan]=="nothing" && span_slices<=span_numad
                            strut_airfoil_filenames[ispan] = "$(numadIn_strut[istrut].airfoil[ispan_numad]).dat"
                        end
                    end
                end

                global airfoil_filenames = [airfoil_filenames; strut_airfoil_filenames]

            end
        end
    end

    OWENSOpenFASTWrappers.writeADinputFile(ad_input_file,blade_filenames,airfoil_filenames,OLAF_filename)

    NumADBody = length(AD15bldNdIdxRng[:,1])
    bld_len = zeros(NumADBody)
    #### bladefileissaved = false
    for (iADBody,filename) in enumerate(blade_filenames)
        strt_idx = AD15bldNdIdxRng[iADBody,1]
        end_idx = AD15bldNdIdxRng[iADBody,2]
        if end_idx<strt_idx
            tmp_end = end_idx
            end_idx = strt_idx
            strt_idx = tmp_end
        end

        #### Get the blade length
        x1 = mymesh.x[strt_idx]
        x2 = mymesh.x[end_idx]
        y1 = mymesh.y[strt_idx]
        y2 = mymesh.y[end_idx]
        z1 = mymesh.z[strt_idx]
        z2 = mymesh.z[end_idx]
        bld_len[iADBody] = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)

        #### Get the blade shape
        ADshapeZ = collect(LinRange(0,H,NumADBldNds))
        xmesh = mymesh.x[strt_idx:end_idx]
        ymesh = mymesh.y[strt_idx:end_idx]
        ADshapeX = sqrt.(xmesh.^2 .+ ymesh.^2)
        ADshapeX .-= ADshapeX[1] #get it starting at zero #TODO: make robust for blades that don't start at 0
        ADshapeXspl = OWENS.safeakima(LinRange(0,H,length(ADshapeX)),ADshapeX,ADshapeZ)
        
        if iADBody<=Nbld #&& !bladefileissaved#Note that the blades can be curved and are assumed to be oriented vertically
            ####  bladefileissaved = true
            BlSpn0=ADshapeZ
            BlCrvAC0=ADshapeXspl

            bladeangle = (iADBody-1)*2.0*pi/Nbld + angularOffset #TODO: pitch offset and twist offset that isn't from the helical

            BlSpn = ADshapeZ
            blade_twist = atan.(xmesh,ymesh).-bladeangle

            #### TODO: reevalueate these equations and make sure they are robust against varying designs
            BlCrvACinput = -ymesh.*sin(bladeangle).+xmesh.*cos(bladeangle)
            BlCrvACinput = BlCrvACinput .- BlCrvACinput[1]
            BlSwpAC = -OWENS.safeakima(LinRange(0,H,length(BlCrvACinput)),BlCrvACinput,ADshapeZ)

            BlSwpACinput = xmesh.*sin(bladeangle).+ymesh.*cos(bladeangle)
            BlSwpACinput = BlSwpACinput .- BlSwpACinput[1]
            BlCrvAC = OWENS.safeakima(LinRange(0,H,length(BlSwpACinput)),BlSwpACinput,ADshapeZ)

            BlCrvAng = zeros(blade_Nnodes[iADBody])

            BlTwistinput =(blade_twist.-blade_twist[1])*180/pi
            BlTwist = OWENS.safeakima(LinRange(0,H,length(BlTwistinput)),BlTwistinput,ADshapeZ)

            BlChord=blade_chords[iADBody]

            BlAFID=collect((iADBody-1)*NumADBldNds+1:iADBody*NumADBldNds)

        elseif iADBody>Nbld # while the arms/struts are assumed to be straight and are oriented by the mesh angle
            BlSpn=collect(LinRange(0,bld_len[iADBody],blade_Nnodes[iADBody]))
            BlCrvAC=zeros(blade_Nnodes[iADBody])
            BlSwpAC=zeros(blade_Nnodes[iADBody])
            BlCrvAng=zeros(blade_Nnodes[iADBody])
            BlTwist=zeros(blade_Nnodes[iADBody])
            BlChord=blade_chords[iADBody]
            BlAFID=collect((iADBody-1)*NumADStrutNds+1:iADBody*NumADStrutNds)
        end      
        OWENSOpenFASTWrappers.writeADbladeFile(filename;NumBlNds=blade_Nnodes[iADBody],BlSpn,BlCrvAC,BlSwpAC,BlCrvAng,BlTwist,BlChord,BlAFID)     
    end

    OWENSOpenFASTWrappers.writeOLAFfile(OLAF_filename;nNWPanel=200,nFWPanels=10)

    OWENSOpenFASTWrappers.writeIWfile(Vinf,ifw_input_file;windINPfilename)

    OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,[shapeX],[shapeZ],[B],[Htwr_base],[mymesh],[myort],[AD15bldNdIdxRng],[AD15bldElIdxRng];
            rho     = rho,
            adi_dt  = delta_t,
            adi_tmax= numTS*delta_t,
            omega   = [omega],
            adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
            adi_DT_Outs = delta_t,   # output frequency
            numTurbines = 1,
            refPos=[[0,0,0]],
            hubPos=[[0,0,0.0]],
            hubAngle=[[0,0,0]],
            nacPos=[[0,0,0]],
            adi_nstrut=[Nstrutperbld],
            adi_debug=0,
            isHAWT = false     # true for HAWT, false for crossflow or VAWT
            )

    aeroForcesAD(t,azi) = OWENS.mapAD15(t,azi,[mymesh],OWENSOpenFASTWrappers.advanceAD15;alwaysrecalc=true,verbosity=1)
    deformAeroAD=OWENSOpenFASTWrappers.deformAD15
end

nothing

# Calculate mass breakout of each material
mass_breakout_bld = OWENS.get_material_mass(plyprops_bld,numadIn_bld)
mass_breakout_blds = mass_breakout_bld.*length(mymesh.structuralNodeNumbers[:,1])
mass_breakout_twr = OWENS.get_material_mass(plyprops_twr,numadIn_twr;int_start=0.0,int_stop=Htwr_base)

################################
##### End setupOWENS
################################

nothing 

# Then the rest of this example is the same as example B
if verbosity>0
    
    println("\nBlades' Mass Breakout")
    for (i,name) in enumerate(plyprops_bld.names)
        println("$name $(mass_breakout_blds[i]) kg, $(plyprops_bld.costs[i]) \$/kg: \$$(mass_breakout_blds[i]*plyprops_bld.costs[i])")
    end
    
    println("\nTower Mass Breakout")
    for (i,name) in enumerate(plyprops_twr.names)
        println("$name $(mass_breakout_twr[i]) kg, $(plyprops_twr.costs[i]) \$/kg: \$$(mass_breakout_twr[i]*plyprops_twr.costs[i])")
    end
    
    println("Total Material Cost Blades: \$$(sum(mass_breakout_blds.*plyprops_bld.costs))")
    println("Total Material Cost Tower: \$$(sum(mass_breakout_twr.*plyprops_twr.costs))")
    println("Total Material Cost: \$$(sum(mass_breakout_blds.*plyprops_bld.costs)+ sum(mass_breakout_twr.*plyprops_twr.costs))")
    
end

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

inputs = OWENS.Inputs(;verbosity,analysisType = structuralModel,
tocp = [0.0,100000.1],
Omegaocp = [RPM,RPM] ./ 60,
tocp_Vinf = [0.0,100000.1],
Vinfocp = [Vinf,Vinf],
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2)

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
dataOutputFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = true,
numNodes = mymesh.numNodes,
RayleighAlpha = 0.05,
RayleighBeta = 0.05,
iterationType = "DI")

println("Running Unsteady")
if AD15On
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesAD,deformAero=deformAeroAD)
else
    t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
    FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
    epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
    topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesACDMS,deformAero=deformAeroACDMS)
end

if AD15On #TODO: move this into the run functions
    OWENSOpenFASTWrappers.endTurb()
end

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)


##########################################
#### Ultimate Failure #####
##########################################

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
topstrainout_tower_U,topstrainout_tower_LtopDamage_blade_U,
topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
mymesh,myel,myort,Nbld,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
Twr_LE_U_idx=1,Twr_LE_L_idx=1,
AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing,strut_precompinput,plyprops_strut,numadIn_strut,lam_U_strut,lam_L_strut) #TODO: add in ability to have material safety factors and load safety factors

nothing