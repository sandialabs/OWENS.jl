using Test
import OWENS

@testset "Strut strain histories use actual strut ranges" begin
    Nbld = 2
    nstrut_per_blade = 3
    n_strut_bodies = Nbld * nstrut_per_blade
    n_bodies = Nbld + n_strut_bodies
    nodes_per_body = 4
    elems_per_body = nodes_per_body - 1
    N_ts = 2

    z = Float64[]
    AD15bldNdIdxRng = zeros(Int, n_bodies, 2)
    AD15bldElIdxRng = zeros(Int, n_bodies, 2)
    for ibody = 1:n_bodies
        start_node = length(z) + 1
        append!(z, (0.0:(nodes_per_body-1)) .+ 10.0 * ibody)
        stop_node = length(z)
        start_el = (ibody - 1) * elems_per_body + 1
        stop_el = ibody * elems_per_body
        AD15bldNdIdxRng[ibody, :] = [start_node, stop_node]
        AD15bldElIdxRng[ibody, :] = [start_el, stop_el]
    end

    total_elems = n_bodies * elems_per_body
    epsilon_x_hist = zeros(4, total_elems, N_ts)
    for iel = 1:total_elems, its = 1:N_ts
        epsilon_x_hist[1, iel, its] = 100iel + its
    end

    meanepsilon_z_hist = epsilon_x_hist[1:1, :, :] .+ 1.0
    meanepsilon_y_hist = epsilon_x_hist[1:1, :, :] .+ 2.0
    kappa_x_hist = epsilon_x_hist .+ 3.0
    kappa_y_hist = epsilon_x_hist .+ 4.0
    kappa_z_hist = epsilon_x_hist .+ 5.0

    mymesh = (; z)
    strut_precompinput = [nothing, nothing, nothing]
    numadIn_strut = (; span = [0.0, 0.5, 1.0])

    histories = OWENS._interpolate_strut_strain_histories(
        mymesh,
        AD15bldNdIdxRng,
        AD15bldElIdxRng,
        Nbld,
        N_ts,
        strut_precompinput,
        numadIn_strut,
        epsilon_x_hist,
        meanepsilon_z_hist,
        meanepsilon_y_hist,
        kappa_x_hist,
        kappa_y_hist,
        kappa_z_hist,
    )

    @test size(histories.eps_x) == (n_strut_bodies, N_ts, length(strut_precompinput))

    first_strut_row = Nbld + 1
    first_strut_start_el = AD15bldElIdxRng[first_strut_row, 1]
    first_strut_stop_el = AD15bldElIdxRng[first_strut_row, 2]
    @test histories.eps_x[1, 1, :] ≈ epsilon_x_hist[
        1,
        first_strut_start_el:first_strut_stop_el,
        1,
    ]

    last_strut_row = size(AD15bldElIdxRng, 1)
    last_strut_start_el = AD15bldElIdxRng[last_strut_row, 1]
    last_strut_stop_el = AD15bldElIdxRng[last_strut_row, 2]
    @test histories.kappa_z[end, 2, :] ≈ kappa_z_hist[
        1,
        last_strut_start_el:last_strut_stop_el,
        2,
    ]
end

@testset "Non-finite fatigue input returns invalid damage" begin
    @test isnan(OWENS.fatigue_damage([NaN, 1.0], [1.0, 2.0], [3.0, 4.0], 100.0))
end
