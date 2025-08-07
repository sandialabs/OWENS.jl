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

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
# runpath = path =  splitdir(@__FILE__)[1]

# Load options using ModelingOptions with the complete YAML file
modelingOptions = OWENS.ModelingOptions("$runpath/sampleOWENS_complete.yml", path=path)

# Preprocess the modeling options to generate setup_options
setup_options = OWENS.preprocess_modeling_options_setup(modelingOptions, path)

# Call the helper function that builds the mesh, calculates the sectional properties,
# and aligns the sectional properties to the mesh elements
# We can now use the ModelingOptions struct directly with setupOWENS
setup_outputs = OWENS.setupOWENS(
    path;
    setup_options = setup_options,
    verbosity = 1,
    return_componentized = true
)

# If the sectional properties material files includes cost information, that is combined with the density 
# to estimate the overall material cost of of materials in the blades
components  = setup_outputs.components

plyprops_twr = components[1].plyProps
plyprops_bld = components[2].plyProps

mass_breakout_blds = components[1].mass
mass_breakout_twr = components[2].mass

println("\nBlades' Mass Breakout")
for (i,name) in enumerate(plyprops_bld.names)
    cost_per_kg = plyprops_bld.costs[i]
    total_cost = mass_breakout_blds[i] * cost_per_kg
    println("$name $(mass_breakout_blds[i]) kg, $(cost_per_kg) \$/kg: \$$(total_cost)")
end

println("\nTower Mass Breakout")
for (i,name) in enumerate(plyprops_twr.names)
    cost_per_kg = plyprops_twr.costs[i]
    total_cost = mass_breakout_twr[i] * cost_per_kg
    println("$name $(mass_breakout_twr[i]) kg, $(cost_per_kg) \$/kg: \$$(total_cost)")
end

total_blade_cost = sum(mass_breakout_blds .* plyprops_bld.costs)
total_tower_cost = sum(mass_breakout_twr .* plyprops_twr.costs)
total_material_cost = total_blade_cost + total_tower_cost

println("Total Material Cost Blades: \$$(total_blade_cost)")
println("Total Material Cost Tower: \$$(total_tower_cost)")
println("Total Material Cost: \$$(total_material_cost)")

# Run the unsteady simulation

println("Running Unsteady Simulation")

unsteady_outputs = OWENS.Unsteady_Land(setup_outputs, modelingOptions, returnold=false)

# Clean up AD15 if it was used
if modelingOptions.OWENS_Options.AeroModel == "AD"
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

# Like described above, we can output vtk files viewable in paraview.  Here it is done for each time step and shows the 
# deformations.  Additionaly, there is a method to input custom values and have them show up on the vtk surface mesh
# for example, strain, or reaction force, etc.  This is described in more detail in the api reference for the function and: TODO

# Save VTK time domain files for visualization
println("Saving VTK time domain files")
# Note: This generates large output files for visualization in ParaView
OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue", unsteady_outputs, setup_outputs; scaling=1)

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
