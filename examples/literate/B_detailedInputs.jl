# # [Detailed Inputs](@id simple2)
# 
# In this example, we show the second level of what is going on behind the precompiled binary
# This would be appropriate if you need more customization in the run and design parameters than the
# input file currently allows, but your design still fits within the setupOWENS helper function etc.
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`B_detailedInputs.ipynb`](@__NBVIEWER_ROOT_URL__/examples/B_detailedInputs.ipynb).
#-

# First we import the packages.  If "using" was employed, then all of the functions of the packages
# specified would be made available, but "import" requires PackageName.FunctionName to be used unless
# the function was explicitely exported in the package.  If "include("filepath/filename.jl)" is used, 
# that is the same as copying and pasting.  Please see the respective page on YAML input (TODO) for a
# description of the YAML inputs

import OWENS
import OWENSAero
#### import PyPlot
# runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
runpath = path =  splitdir(@__FILE__)[1]

# Load the unified options using ModelingOptions with the complete YAML file
modelingOptions = OWENS.ModelingOptions("$runpath/sampleOWENS_complete.yml", path=path)

# Call the helper function that builds the mesh, calculates the sectional properties,
# and aligns the sectional properties to the mesh elements
# We can now use the ModelingOptions struct directly with setupOWENS

# Method 1: Using the new preprocess_unified_options_setup function
setup_options = OWENS.preprocess_unified_options_setup(modelingOptions, path)

mymesh, myel, myort, myjoint, sectionPropsArray, mass_twr, mass_bld,
stiff_twr, stiff_bld, bld_precompinput,
bld_precompoutput, plyprops_bld, numadIn_bld, lam_U_bld, lam_L_bld,
twr_precompinput, twr_precompoutput, plyprops_twr, numadIn_twr, lam_U_twr, lam_L_twr, aeroForces, deformAero,
mass_breakout_blds, mass_breakout_twr, system, assembly, sections, AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(
    path;
    setup_options = setup_options,
    verbosity = 1
)


nothing

# Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
# this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
# includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
# to work with OWENS inputs as well.

nothing

# If the sectional properties material files includes cost information, that is combined with the density 
# to estimate the overall material cost of of materials in the blades

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

nothing

# Here we apply the boundary conditions.  For this case, with a regular cantelever tower, the tower base node which is 
# 1 is constrained in all 6 degrees of freedom to have a displacement of 0.  You can change this displacement to allow for things
# like pretension, and you can apply boundary conditions to any node.

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

nothing

# There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options

if AeroModel=="AD"
    AD15On = true
else
    AD15On = false
end

inputs = OWENS.Inputs(;
    verbosity,
    analysisType = structuralModel,
    tocp = owens_options.Prescribed_RPM_time_controlpoints,
    Omegaocp = owens_options.Prescribed_RPM_RPM_controlpoints ./ 60,
    tocp_Vinf = owens_options.Prescribed_Vinf_time_controlpoints,
    Vinfocp = owens_options.Prescribed_Vinf_Vinf_controlpoints,
    numTS,
    delta_t,
    AD15On,
    aeroLoadsOn = owens_options.aeroLoadsOn
)

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

# Use OWENSFEA_Options for FEA model parameters
fea_options = modelingOptions.OWENSFEA_Options

feamodel = OWENS.FEAModel(;
    analysisType = structuralModel,
    dataOutputFilename = owens_options.dataOutputFilename,
    joint = myjoint,
    platformTurbineConnectionNodeNumber = fea_options.platformTurbineConnectionNodeNumber,
    pBC,
    nlOn = fea_options.nlOn,
    numNodes = mymesh.numNodes,
    RayleighAlpha = fea_options.RayleighAlpha,
    RayleighBeta = fea_options.RayleighBeta,
    iterationType = fea_options.iterationType
)

nothing

# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist = OWENS.Unsteady_Land(inputs;system,assembly,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

if AD15On #TODO: move this into the run functions
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

nothing

# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

nothing

# This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
# it also calculates analytical simply supported buckling safety factors

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
composite_station_idx_U_strut = [1,6,3,2,5],
composite_station_name_U_strut = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_L_strut = [1,6,3,2,5],
composite_station_name_L_strut = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_U_bld = [1,6,3,2,5],
composite_station_name_U_bld = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
composite_station_idx_L_bld = [1,6,3,2,5],
composite_station_name_L_bld = ["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
Twr_LE_U_idx=1,Twr_LE_L_idx=1,
AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors

nothing
