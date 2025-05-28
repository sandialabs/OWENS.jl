
using OWENS
using Random
using Test

# Random stress timeseries.
amplitude, mean_stress = 150e6, 30e6
stress_timeseries = (rand(Xoshiro(1), 1000) .- 0.5) * amplitude .+ mean_stress

# Ultimate strength
ultimate_strength = 500e6

#S-N curve
failure_strength_coefficient = 250e6
fatigue_exponent = -0.12
ncycles = 1:1000:(1e5+1)
sn_log_cycles = log10.(ncycles)
sn_stress = failure_strength_coefficient * ncycles .^ fatigue_exponent
ncycles = vcat(ncycles, 1e7)
sn_log_cycles = vcat(sn_log_cycles, log10(1e7))
sn_stress = vcat(sn_stress, 1e5)


# Fatigue damage
damage = OWENS.fatigue_damage(stress_timeseries, sn_stress, sn_log_cycles, ultimate_strength; nbins_amplitude=21, nbins_mean=3, mean_correction=true)
@test isapprox(damage, 0.002632525050465263)
