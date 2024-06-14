using StaticArrays
using GXBeam
import ForwardDiff
import FiniteDiff
import OWENS
import OWENSFEA
import OWENSAero
import FLOWMath
import DelimitedFiles
using Statistics:mean
using Test
import HDF5
import YAML
using StructTypes
import OrderedCollections

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# function runprofilefunction()
path = runpath = splitdir(@__FILE__)[1]

Inp = OWENS.MasterInput("$path/RM2_OWENS_modeling_options.yaml")

windio = YAML.load_file("$path/RM1.yaml"; dicttype=OrderedCollections.OrderedDict{Symbol,Any})


# mutable struct myTest
#     can_log
#     sys_log
#     version
#     info
#     myfloat
#     myTest() = new()
# end

# yamlInputtest = YAML.load_file("$path/testdata.yaml"; dicttype=OrderedCollections.OrderedDict{Symbol,Any})
# test_dict = yamlInputtest[:tests][1]
# println(test_dict)

# StructTypes.StructType(::Type{myTest}) = StructTypes.Mutable()
# test = StructTypes.constructfrom(myTest, test_dict)



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
AModel = Inp.AModel
windINPfilename = Inp.windINPfilename
ifw_libfile = Inp.ifw_libfile
Blade_Height = Inp.Blade_Height
Blade_Radius = Inp.Blade_Radius
numTS = Inp.numTS
delta_t = Inp.delta_t
NuMad_geom_xlscsv_file_twr = Inp.NuMad_geom_xlscsv_file_twr
NuMad_mat_xlscsv_file_twr = Inp.NuMad_mat_xlscsv_file_twr
NuMad_geom_xlscsv_file_bld = Inp.NuMad_geom_xlscsv_file_bld
NuMad_mat_xlscsv_file_bld = Inp.NuMad_mat_xlscsv_file_bld
NuMad_geom_xlscsv_file_strut = Inp.NuMad_geom_xlscsv_file_strut
NuMad_mat_xlscsv_file_strut = Inp.NuMad_mat_xlscsv_file_strut
adi_lib = Inp.adi_lib
adi_rootname = Inp.adi_rootname

##############################################
# Setup
#############################################

tocp = [0.0;10.0; 1e6]
Omegaocp = [34.0; 34.0; 34.0]./60 #control inputs

tocp_Vinf = [0.0;10.0; 1e6]
Vinfocp = [10.0;10.0;10.0]

# windio
controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = FLOWMath.akima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
SNL34_unit_xz = [x_shape;;z_shape]
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])
SNL34Z = SNL34z.*Blade_Height #windio
SNL34X = SNL34x.*Blade_Radius #windio

shapeZ = SNL34Z#collect(LinRange(0,H,Nslices+1))
shapeX = SNL34X#R.*(1.0.-4.0.*(shapeZ/H.-.5).^2)#shapeX_spline(shapeZ)

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections = OWENS.setupOWENS(OWENSAero,path;
    rho, #windio
    Nslices, #modeling options
    ntheta, #modeling options
    RPM, #remove
    Vinf, #remove
    eta, #windio
    B = Nbld, #windio
    H = Blade_Height, #windio
    R = Blade_Radius, #windio
    shapeZ, #windio
    shapeX, #windio
    ifw, #modeling options
    delta_t, #modeling options
    numTS, #modeling options
    adi_lib, #remove - make smart enough to find
    adi_rootname, #modeling options
    AD15hubR = 0.0, #modeling options
    windINPfilename, #modeling options
    ifw_libfile, #remove - make smart enough to find
    NuMad_geom_xlscsv_file_twr,#windio
    NuMad_mat_xlscsv_file_twr,#windio
    NuMad_geom_xlscsv_file_bld,#windio
    NuMad_mat_xlscsv_file_bld,#windio
    NuMad_geom_xlscsv_file_strut,#windio
    NuMad_mat_xlscsv_file_strut,#windio
    Ht=towerHeight,
    ntelem, 
    nbelem, 
    ncelem,
    nselem,
    joint_type = 0,
    strut_mountpointbot = 0.03,
    strut_mountpointtop = 0.03,
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = pi/2,
    meshtype = turbineType)

top_idx = 23#Int(myjoint[7,2])
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
top_idx 1 0
top_idx 2 0
top_idx 3 0
top_idx 4 0
top_idx 5 0]

if AModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;analysisType = structuralModel,
    tocp,
    Omegaocp,
    tocp_Vinf,
    Vinfocp,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = 2,
    turbineStartup = 1,
    generatorOn = true,
    useGeneratorFunction = true,
    driveTrainOn = true,
    JgearBox = 250.0,#(2.15e3+25.7)/12*1.35582*100,
    gearRatio = 1.0,
    gearBoxEfficiency = 1.0,
    driveShaftProps = OWENS.DriveShaftProps(10000,1.5e2), #8.636e5*1.35582*0.6
    OmegaInit = Omegaocp[1])

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn=false,
    numNodes = mymesh.numNodes,
    numModes = 200,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    iterationType = "DI")

function myoptfun(xin;feamodel,mymesh,myel,inputs)
    feamodel.analysisType = "S"

    Fdof = [15,16,17]
    Fexternal = xin[1:3]
    Omega = xin[4]

    displ=zeros(mymesh.numNodes*6)
    elStorage = OWENS.OWENSFEA.initialElementCalculations(feamodel,myel,mymesh)
    displ,elStrain,staticAnalysisSuccessful,FReaction = OWENS.OWENSFEA.staticAnalysis(feamodel,mymesh,myel,displ,Omega,0.0,elStorage;Fdof,Fexternal)
    objective = FReaction[6]*Omega
    constraints = 0.1-displ[end]
    return [objective;constraints]
end

"""
element_strain(element, F, M)

Calculate the strain of a beam element given the resultant forces and moments applied on
the element expressed in the deformed beam element frame
"""
@inline function element_strain(element, F, M)
    C = element.compliance
    S11 = C[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = C[SVector{3}(1:3), SVector{3}(4:6)]
    return S11*F + S12*M
end

"""
    element_curvature(element, F, M)

Calculate the curvature of a beam element given the resultant force and moments applied on
the element expressed in the deformed beam element frame
"""
@inline function element_curvature(element, F, M)
    C = element.compliance
    S21 = C[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = C[SVector{3}(4:6), SVector{3}(4:6)]
    return S21*F + S22*M
end

function myoptfunGX(xin;assembly)
    deformedxyz = zeros(Real,length(assembly.points),3)
    system = GXBeam.StaticSystem(assembly)
    strainGX = zeros(Real,3,length(assembly.start))
    curvGX = zeros(Real,3,length(assembly.start))

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
    # fixed left side
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    # shear force on right tip
    length(assembly.start)+1 => PrescribedConditions(Fz_follower = xin[2])
    )

    # perform a static analysis
    static_analysis!(system, assembly;
    prescribed_conditions=prescribed_conditions,linear = false)

    # process state and eigenstates
    states = AssemblyState(system, assembly;
    prescribed_conditions = prescribed_conditions)

    for (ipt,point) in enumerate(states.points)
        deformedxyz[ipt,:] = point.u
    end

    for iel = 1:length(states.elements)
        strainGX[:,iel] = element_strain(assembly.elements[iel],states.elements[iel].Fi,states.elements[iel].Mi)
        curvGX[:,iel] = element_curvature(assembly.elements[iel],states.elements[iel].Fi,states.elements[iel].Mi)
    end

    objective = curvGX[3,2]*xin[1]
    constraints = 0.1 .- curvGX[2,:]
    return [objective;constraints]

end

xin = [1000.0,1000.0]

myoptfunGX0(x) = myoptfunGX(x;assembly)
J = ForwardDiff.jacobian(myoptfunGX0, xin)

function myoptfunGX2(x)
    return Float64.(myoptfunGX0(x))
end

J2 = FiniteDiff.finite_difference_jacobian(myoptfunGX2,xin)

for i1 = 1:length(J[:,1])
    for i2 = 1:length(J[1,:])
        @test isapprox(J[i1,i2],J2[i1,i2];atol=1e-5)
    end
end