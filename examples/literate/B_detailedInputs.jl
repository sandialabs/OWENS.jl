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
runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
# runpath = path =  splitdir(@__FILE__)[1]

# Load the unified options using ModelingOptions with the complete YAML file
modelingOptions = OWENS.ModelingOptions("$runpath/sampleOWENS_complete.yml", path=path)

# Call the helper function that builds the mesh, calculates the sectional properties,
# and aligns the sectional properties to the mesh elements
# We can now use the ModelingOptions struct directly with setupOWENS

# Run preprocessing function to generate setup_options from modelingOptions
setup_options = OWENS.preprocess_modeling_options_setup(modelingOptions, path)

# Use the new struct output format
setup_outputs = OWENS.setupOWENS(
    path;
    setup_options = setup_options,
    verbosity = 1,
    return_componentized = true
)

nothing

# Optionally, we can run the finite element solver with gemetrically exact beam theory via GXBeam.jl
# this requires that the OWENS style inputs are converted to the GXBeam inputs.  This interface also
# includes the ability to output VTK files, which can be viewed in paraview.  We have adapted this interface
# to work with OWENS inputs as well.

nothing

# If the sectional properties material files includes cost information, that is combined with the density 
# to estimate the overall material cost of of materials in the blades
components  = setup_outputs.components

plyprops_twr = components[1].plyProps
plyprops_bld = components[2].plyProps

mass_breakout_blds = components[1].mass
mass_breakout_twr = components[2].mass

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


nothing

# Here we apply the boundary conditions.  For this case, with a regular cantelever tower, the tower base node which is 
# 1 is constrained in all 6 degrees of freedom to have a displacement of 0.  You can change this displacement to allow for things
# like pretension, and you can apply boundary conditions to any node.

# pBC = [1 1 0
# 1 2 0
# 1 3 0
# 1 4 0
# 1 5 0
# 1 6 0]

nothing

# There are inputs for the overall coupled simulation, please see the api reference for specifics on all the options


# inputs = OWENS.Inputs(;
#     verbosity = modelingOptions.OWENS_Options.verbosity,
#     analysisType = modelingOptions.OWENS_Options.structuralModel,
#     tocp = [0.0,100000.1],
#     Omegaocp = [modelingOptions.OWENSAero_Options.RPM,modelingOptions.OWENSAero_Options.RPM] ./ 60,
#     tocp_Vinf = [0.0,100000.1],
#     Vinfocp = [modelingOptions.OWENSAero_Options.Vinf,modelingOptions.OWENSAero_Options.Vinf],
#     numTS = modelingOptions.OWENS_Options.numTS,
#     delta_t = modelingOptions.OWENS_Options.delta_t,
#     AD15On = modelingOptions.OWENS_Options.AeroModel=="AD",
#     aeroLoadsOn = modelingOptions.OWENS_Options.aeroLoadsOn
# )

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)
# Use OWENSFEA_Options for FEA model parameters
# fea_options = modelingOptions.OWENSFEA_Options

# feamodel = OWENS.FEAModel(;
#     analysisType = modelingOptions.OWENS_Options.structuralModel,
#     dataOutputFilename = "none",
#     joint = myjoint,
#     platformTurbineConnectionNodeNumber = 1,
#     pBC,
#     nlOn = fea_options.nlOn,
#     numNodes = mymesh.numNodes,
#     RayleighAlpha = 0.05,
#     RayleighBeta = 0.05,
#     iterationType = "DI")

nothing

# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.

println("Running Unsteady")

unsteady_outputs = OWENS.Unsteady_Land(setup_outputs, modelingOptions, returnold=false)

# Clean up AD15 if it was used
if modelingOptions.OWENS_Options.AeroModel == "AD"
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

nothing

# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue", unsteady_outputs, setup_outputs; scaling=1)

nothing


#### Populate Components with Strain Data

# Populate the components with strain data from the unsteady simulation
for icomp in 1:length(setup_outputs.components)
    component = setup_outputs.components[icomp]
    
    # Get the element numbers for this component
    startE = component.elNumbers[1]
    stopE = component.elNumbers[end]
    
    println("Populating strain data for component: $(component.name)")
    println("  Element range: $startE to $stopE")
    println("  Strain history dimensions: $(size(unsteady_outputs.epsilon_x_hist))")
    
    # Extract strain data for this component's elements
    # The strain history arrays are organized as [4, element, time_step] where 14 is quadrature points
    # We extract the first quadrature point for each element: [1, startE:stopE, :] gives [element, time_step]
    # safetyfactor_fatigue expects [element, time_step] format
    component.e_x = unsteady_outputs.epsilon_x_hist[1, startE:stopE, :]  # [element, time_step]
    component.e_y = unsteady_outputs.epsilon_y_hist[1, startE:stopE, :]  # [element, time_step] 
    component.e_z = unsteady_outputs.epsilon_z_hist[1, startE:stopE, :]  # [element, time_step]
    component.k_x = unsteady_outputs.kappa_x_hist[1, startE:stopE, :]    # [element, time_step]
    component.k_y = unsteady_outputs.kappa_y_hist[1, startE:stopE, :]    # [element, time_step]
    component.k_z = unsteady_outputs.kappa_z_hist[1, startE:stopE, :]    # [element, time_step]
    
    println("  Populated strain data dimensions: $(size(component.e_x))")
end

nothing

#### Ultimate Failure #####

# Create PostProcessingOptions
postprocess_options = OWENS.PostProcessOptions(
    verbosity=modelingOptions.OWENS_Options.verbosity,
    usestationBld=0,
    usestationStrut=0,
    throwawayTimeSteps=1,
    calculate_fatigue=true)

# This helper function looks through all the loads and picks out the worst case safety factor in each of the stacks of composite lamina
# it also calculates analytical simply supported buckling safety factors
components = OWENS.safetyfactor_fatigue(
    setup_outputs.mymesh,
    setup_outputs.components,
    modelingOptions.OWENS_Options.delta_t;
    options=postprocess_options
)

nothing
