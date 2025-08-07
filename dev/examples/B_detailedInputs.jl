import OWENS
import OWENSAero

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]

modelingOptions = OWENS.ModelingOptions("$runpath/sampleOWENS_complete.yml", path=path)

setup_options = OWENS.preprocess_modeling_options_setup(modelingOptions, path)

setup_outputs = OWENS.setupOWENS(
    path;
    setup_options = setup_options,
    verbosity = 1,
    return_componentized = true
)

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

println("Running Unsteady Simulation")

unsteady_outputs = OWENS.Unsteady_Land(setup_outputs, modelingOptions, returnold=false)

if modelingOptions.OWENS_Options.AeroModel == "AD"
    OWENS.OWENSOpenFASTWrappers.endTurb()
end

println("Saving VTK time domain files")

OWENS.OWENSFEA_VTK("$path/vtk/SNLARCUS5MW_timedomain_TNBnltrue", unsteady_outputs, setup_outputs; scaling=1)

postprocess_options = OWENS.PostProcessOptions(
    verbosity=modelingOptions.OWENS_Options.verbosity,
    usestationBld=0,
    usestationStrut=0,
    throwawayTimeSteps=1,
    calculate_fatigue=true)

components = OWENS.safetyfactor_fatigue(
    setup_outputs.mymesh,
    setup_outputs.components,
    modelingOptions.OWENS_Options.delta_t;
    options=postprocess_options
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
