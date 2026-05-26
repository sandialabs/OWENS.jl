using Test
using LinearAlgebra
import JLD2
import OWENS

function make_restart_topdata(; numTS = 6, ndof = 3, nelem = 2, delta_t = 0.25)
    values = Dict{Symbol,Any}(field => nothing for field in fieldnames(OWENS.TopData))

    values[:delta_t] = delta_t
    values[:numTS] = numTS
    values[:numDOFPerNode] = 6
    values[:CN2H] = Matrix{Float64}(I, 3, 3)
    values[:t] = collect(range(0.0, length = numTS, step = delta_t))
    values[:integrator] = 0.0
    values[:integrator_j] = 0.0
    values[:topDispOut] = Float64[]

    scalar_state_histories = Set([
        :aziHist,
        :OmegaHist,
        :OmegaDotHist,
        :gbHist,
        :gbDotHist,
        :gbDotDotHist,
        :genTorque,
        :genPower,
        :torqueDriveShaft,
    ])
    for field in OWENS._RESTART_STATE_HISTORY_FIELDS
        values[field] =
            field in scalar_state_histories ? zeros(Float64, numTS) : zeros(Float64, numTS, ndof)
    end
    for field in OWENS._RESTART_STEP_HISTORY_FIELDS
        values[field] = zeros(Float64, 4, nelem, numTS)
    end

    for field in (
        :u_s,
        :udot_s,
        :uddot_s,
        :u_sm1,
        :topFexternal,
        :u_sRed,
        :udot_sRed,
        :uddot_sRed,
        :u_s2,
        :udot_s2,
        :uddot_s2,
        :eta_s,
        :etadot_s,
        :etaddot_s,
        :u_j,
        :udot_j,
        :uddot_j,
        :FReactionsm1,
        :topFReaction_j,
    )
        values[field] = zeros(Float64, ndof)
    end
    values[:rbData] = zeros(Float64, ndof)

    for field in (
        :gb_s,
        :gbDot_s,
        :gbDotDot_s,
        :azi_s,
        :Omega_s,
        :OmegaDot_s,
        :genTorque_s,
        :torqueDriveShaft_s,
        :rotorSpeedForGenStart,
        :topsideMass,
        :azi_j,
        :Omega_j,
        :OmegaDot_j,
        :gb_j,
        :gbDot_j,
        :gbDotDot_j,
        :genTorque_j,
        :torqueDriveShaft_j,
    )
        values[field] = 0.0
    end
    values[:topsideMOI] = zeros(Float64, 3, 3)
    values[:topsideCG] = zeros(Float64, 3)

    return OWENS.TopData((values[field] for field in fieldnames(OWENS.TopData))...)
end

function fill_restart_state!(topdata, completed_step)
    history_index = completed_step + 1
    ndof = size(topdata.uHist, 2)

    for index = 1:history_index
        topdata.uHist[index, :] .= 10 .* index .+ (1:ndof)
        topdata.udotHist[index, :] .= 20 .* index .+ (1:ndof)
        topdata.uddotHist[index, :] .= 30 .* index .+ (1:ndof)
        topdata.FReactionHist[index, :] .= 40 .* index .+ (1:ndof)
        topdata.topFexternal_hist[index, :] .= 50 .* index .+ (1:ndof)
        topdata.aziHist[index] = index / 10
        topdata.OmegaHist[index] = 0.0
        topdata.OmegaDotHist[index] = 0.0
        topdata.genTorque[index] = 100 + index
    end

    for index = 1:completed_step
        topdata.epsilon_x_hist[:, :, index] .= index
        topdata.epsilon_y_hist[:, :, index] .= index + 10
        topdata.epsilon_z_hist[:, :, index] .= index + 20
        topdata.kappa_x_hist[:, :, index] .= index + 30
        topdata.kappa_y_hist[:, :, index] .= index + 40
        topdata.kappa_z_hist[:, :, index] .= index + 50
    end

    topdata.integrator = 1.25
    topdata.integrator_j = 2.5
    topdata.u_s .= 1000 .+ (1:ndof)
    topdata.udot_s .= 1100 .+ (1:ndof)
    topdata.uddot_s .= 1200 .+ (1:ndof)
    topdata.u_sm1 .= 900 .+ (1:ndof)
    topdata.topFReaction_j .= 1300 .+ (1:ndof)
    topdata.FReactionsm1 .= 1400 .+ (1:ndof)
    topdata.azi_s = 3.5
    topdata.Omega_s = 0.0
    topdata.genTorque_s = 77.0
    return topdata
end

@testset "compact restart state restore" begin
    completed_step = 3
    source = fill_restart_state!(make_restart_topdata(), completed_step)
    inputs = OWENS.Inputs(; generatorOn = true, omegaControl = true)

    mktempdir() do directory
        restart_file = joinpath(directory, "restart.jld2")
        OWENS._write_restart_state(restart_file, source, completed_step; inputs)

        file_data = JLD2.load(restart_file)
        @test haskey(file_data, "restart_state")
        @test !haskey(file_data, "topdata")

        restart_state = OWENS._load_restart_state(restart_file)
        target = make_restart_topdata(numTS = 9)
        target_inputs = OWENS.Inputs(; generatorOn = false, omegaControl = false)
        restore = OWENS._restore_restart!(target, restart_state; inputs = target_inputs)

        @test restore.completed_step == completed_step
        @test restore.history_index == completed_step + 1
        @test target.numTS == 9
        @test target.uHist[1:4, :] == source.uHist[1:4, :]
        @test all(iszero, target.uHist[5:end, :])
        @test target.epsilon_x_hist[:, :, 1:3] == source.epsilon_x_hist[:, :, 1:3]
        @test all(iszero, target.epsilon_x_hist[:, :, 4:end])
        @test target.u_s == source.u_s
        @test target.FReactionsm1 == source.FReactionsm1
        @test target_inputs.generatorOn
        @test target_inputs.omegaControl

        @test_throws DimensionMismatch OWENS._restore_restart!(
            make_restart_topdata(numTS = 9, ndof = 4),
            restart_state,
        )
        @test_throws ArgumentError OWENS._restore_restart!(
            make_restart_topdata(numTS = 3),
            restart_state,
        )
        @test_throws ArgumentError OWENS._restore_restart!(
            make_restart_topdata(numTS = 9, delta_t = 0.5),
            restart_state,
        )
    end
end

@testset "legacy restart inference does not depend on rotor speed" begin
    completed_step = 3
    source = fill_restart_state!(make_restart_topdata(), completed_step)
    source.OmegaHist .= 0.0

    mktempdir() do directory
        restart_file = joinpath(directory, "legacy_restart.jld2")
        JLD2.jldsave(restart_file; topdata = source)

        restart_state = OWENS._load_restart_state(restart_file)
        @test restart_state["completed_step"] == completed_step
        @test restart_state["history_index"] == completed_step + 1
    end
end

@testset "OWENSAero backend restart state restores mutable globals" begin
    OWENS.OWENSAero.setupTurb(
        [0.0, 1.0],
        [0.0, 1.0],
        2,
        [0.2],
        1.0,
        5.0;
        ntheta = 4,
        Nslices = 1,
        DynamicStallModel = "analytic",
    )
    saved_omega = [1.0, 2.0, 3.0, 4.0]
    saved_vx = [5.0, 6.0, 7.0, 8.0]
    OWENS.OWENSAero.turbslices[1].omega .= saved_omega
    OWENS.OWENSAero.envslices[1].V_x .= saved_vx
    OWENS._set_module_global!(OWENS.OWENSAero, :last_step1, 4)

    backend_state = OWENS._capture_restart_backend_states(capture_owensaero = true)

    OWENS.OWENSAero.turbslices[1].omega .= 99.0
    OWENS.OWENSAero.envslices[1].V_x .= 88.0
    OWENS._set_module_global!(OWENS.OWENSAero, :last_step1, 12)

    OWENS._restore_restart_backend_states!(backend_state)
    @test OWENS.OWENSAero.turbslices[1].omega == saved_omega
    @test OWENS.OWENSAero.envslices[1].V_x == saved_vx
    @test OWENS.OWENSAero.last_step1 == 4
end
