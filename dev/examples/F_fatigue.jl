using OWENS
using Statistics
using Random
using Plots

# Experimental Data
ncycle_log_exp = [4, 5, 6, 7, 8.7]
ncycles_exp = 10 .^ ncycle_log_exp
mean_0 = [365.7, 260.3, 190.1, 157.2, 136.3]
mean_138 = [316.2, 219.4, 155.8, 127.8, 117.2]
mean_276 = [249.7, 173.1, 123.9, 103.1, 95.6]
mean_414 = [157.3, 114.3, 84.2, 69.7, 67.7]
mean_n138 = [missing, missing, 218.0, 176.8, 157.1];
nothing

ultimate_strength = 572e6
sn_stress = mean_0 * 1e6
sn_log_cycles = log10.(ncycles_exp)
nothing

ncyc, npercyc = 10, 10_000
time = range(0, stop=2π * ncyc, length=ncyc * npercyc)
function stress(stress_amplitude, stress_mean)
    return stress_amplitude * cos.(time) .+ stress_mean + randn(length(time)) * stress_amplitude * 0.01
end
nothing

stress_amplitude, stress_mean = sn_stress[2], 138e6
println("Stress amplitude: $(stress_amplitude*1e-6) MPa")
println("Mean stress: $(stress_mean*1e-6) MPa")

# plot
plot(time, stress(stress_amplitude, stress_mean) * 1e-6, label=:none, xlabel="Time (s)", ylabel="Stress (MPa)", title="Stress Cycle", grid=true, xlims=[time[begin], time[end]])
hline!([0,], c=:black, label=:none)
hline!([stress_mean] * 1e-6, c=:black, label=:none, linestyle=:dashdot)
hline!([stress_mean + stress_amplitude, stress_mean - stress_amplitude] * 1e-6, c=:black, label=:none, linestyle=:dash)
nothing

nbins_amplitude, nbins_mean = 21, 3
mean_levels, amplitude_levels, ncycles, amplitude_levels_effective = OWENS.rainflow_mean_corrected(stress(stress_amplitude, stress_mean), ultimate_strength; nbins_amplitude, nbins_mean)
nothing

mean_bins_mid = ceil(Int, nbins_mean / 2)
display(mean_levels * 1e-6)
println("Real mean stress: $(stress_mean*1e-6) MPa")
println("Center bin mean stress: $(mean_levels[mean_bins_mid]*1e-6) MPa")
nothing

display(amplitude_levels * 1e-6)
println("Real stress amplitude: $(stress_amplitude*1e-6) MPa")
println("Center bin stress amplitude: $(amplitude_levels[end]*1e-6) MPa")
nothing

display(ncycles)
println("Cycles at desired bin: $(ncycles[end, mean_bins_mid])")
nothing

println("Real stress amplitude: $(stress_amplitude*1e-6) MPa")
println("Effective stress amplitude: $(amplitude_levels_effective[end, mean_bins_mid]*1e-6) MPa")
nothing

log_ncycles_fail = OWENS.sn_curve_mean_corrected(sn_stress, sn_log_cycles, amplitude_levels_effective)
nothing

println("Point on effective S-N curve for mean=$(stress_mean*1e-6) MPa")
println("    N-cycles to failure: 10^$(log_ncycles_fail[end, mean_bins_mid])")
println("    Stress amplitude: $(stress_amplitude*1e-6) MPa")
nothing

# 138
stress_mean = 138e6
stress_amplitudes_138 = sn_stress[1:4]
logN_138 = zeros(0)
for stress_amplitude in stress_amplitudes_138
    mean_levels, amplitude_levels, ncycles, amplitude_levels_effective = OWENS.rainflow_mean_corrected(stress(stress_amplitude, stress_mean), ultimate_strength; nbins_amplitude, nbins_mean)
    log_ncycles_fail = OWENS.sn_curve_mean_corrected(sn_stress, sn_log_cycles, amplitude_levels_effective)
    append!(logN_138, log_ncycles_fail[end, mean_bins_mid])
end

# -138
stress_mean = -138e6
stress_amplitudes_n138 = sn_stress[3:5]
logN_n138 = zeros(0)
for stress_amplitude in stress_amplitudes_n138
    mean_levels, amplitude_levels, ncycles, amplitude_levels_effective = OWENS.rainflow_mean_corrected(stress(stress_amplitude, stress_mean), ultimate_strength; nbins_amplitude, nbins_mean)
    log_ncycles_fail = OWENS.sn_curve_mean_corrected(sn_stress, sn_log_cycles, amplitude_levels_effective)
    append!(logN_n138, log_ncycles_fail[end, mean_bins_mid])
end

# 276
stress_mean = 276e6
stress_amplitudes_276 = sn_stress[1:2]
logN_276 = zeros(0)
for stress_amplitude in stress_amplitudes_276
    mean_levels, amplitude_levels, ncycles, amplitude_levels_effective = OWENS.rainflow_mean_corrected(stress(stress_amplitude, stress_mean), ultimate_strength; nbins_amplitude, nbins_mean)
    log_ncycles_fail = OWENS.sn_curve_mean_corrected(sn_stress, sn_log_cycles, amplitude_levels_effective)
    append!(logN_276, log_ncycles_fail[end, mean_bins_mid])
end

# plot
color_0 = :black
color_138, color_276, color_414 = :blue, :red, :green
markersize, markercolor = 5, :white
# mean
plot(ncycles_exp, mean_0, label=:none, color=color_0, lw=3, xaxis=:log, title="Experimental S-N Curves for 7075-T6 Aluminum\nkt=1, axial", xlabel="Cycles to Faliures", ylabel="Stress Amplitude (MPa)", legend=:topright, grid=true)
scatter!(ncycles_exp, mean_0, label="σₘ=0", color=markercolor, markerstrokecolor=color_0, markersize=markersize)
# positive
for (mean, color, label) in zip([mean_138, mean_276, mean_414], [color_138, color_276, color_414], [138, 276, 414])
    plot!(ncycles_exp, mean, label=:none, color=color)
    scatter!(ncycles_exp, mean, label="σₘ=$label MPa", color=markercolor, markerstrokecolor=color, markersize=markersize)
end
# negative
for (mean, color, label) in zip([mean_n138,], [color_138,], [-138,])
    plot!(ncycles_exp, mean, label=:none, color=color, linestyle=:dash)
    scatter!(ncycles_exp, mean, label="σₘ=$label MPa", color=markercolor, markerstrokecolor=color, markersize=markersize)
end
# Miner's rule + Goodman's correction
plot!(10 .^ logN_138, stress_amplitudes_138 * 1e-6, color=color_138, markershape=:star5, markerstrokecolor=color_138, markersize=markersize, label="OWENS: σₘ=138 MPa")
plot!(10 .^ logN_276, stress_amplitudes_276 * 1e-6, color=color_276, markershape=:star5, markerstrokecolor=color_276, markersize=markersize, label="OWENS: σₘ=276 MPa")
plot!(10 .^ logN_n138, stress_amplitudes_n138 * 1e-6, color=color_138, linestyle=:dashdot, markershape=:star5, markerstrokecolor=color_138, markersize=markersize, label="OWENS: σₘ=-138 MPa")
plot!()
nothing

amplitude, mean = 150e6, 30e6
time = 1:1000
rng = Xoshiro(1) # seed for reproducibility
stress_timeseries = (rand(rng, length(time)) .- 0.5) * amplitude .+ mean

plot(time, stress_timeseries; xlabel="time [s]", ylabel="stress [Pa]", label=nothing, xlim=(0, Inf))
hline!([mean], label=nothing, color=:black, linestyle=:dash)
hline!([0], label=nothing, color=:black, linestyle=:solid)
nothing

# Ultimate strength
ultimate_strength = 500e6
println("Ultimate strength: ", ultimate_strength, " Pa")
nothing

# S-N curve
failure_strength_coefficient = 250e6
fatigue_exponent = -0.12
ncycles = 1:1000:(1e5+1)
sn_log_cycles = log10.(ncycles)
sn_stress = failure_strength_coefficient * ncycles .^ fatigue_exponent
ncycles = vcat(ncycles, 1e7)
sn_log_cycles = vcat(sn_log_cycles, log10(1e7))
sn_stress = vcat(sn_stress, 1e5)
plot(ncycles, sn_stress; xscale=:log10, xlabel="log(cycles)", ylabel="stress [Pa]", label=nothing)
nothing

damage = OWENS.fatigue_damage(stress_timeseries, sn_stress, sn_log_cycles, ultimate_strength; nbins_amplitude=21, nbins_mean=3, mean_correction=true)
println("Damage: ", damage)
nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
