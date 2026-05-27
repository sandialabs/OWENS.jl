using LinearAlgebra
using Test
import OWENS

@testset "Generator and drivetrain utilities" begin
    inputs = (
        ratedTorque = 100.0,
        ratedGenSlipPerc = 10.0,
        zeroTorqueGenSpeed = 1.0,
        pulloutRatio = 2.0,
    )

    @test OWENS.simpleGenerator(inputs, -0.1) == 0.0
    @test OWENS.simpleGenerator(inputs, inputs.zeroTorqueGenSpeed) == 0.0
    @test OWENS.simpleGenerator(inputs, 1.1) ≈ inputs.ratedTorque
    @test OWENS.simpleGenerator(inputs, 0.95) ≈ -50.0
    @test OWENS.simpleGenerator(inputs, 1.3) ≈ 200.0

    drive_shaft = OWENS.DriveShaftProps(1200.0, 80.0)
    @test OWENS.calculateDriveShaftReactionTorque(drive_shaft, 0.30, 0.10, 2.0, 1.5) ≈
          280.0
end

@testset "Newmark rotor integration" begin
    unp1, udotnp1, uddotnp1, fhat =
        OWENS.timeIntegrateSubSystem(1.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.0, 0.0)
    @test unp1 ≈ 0.0025
    @test udotnp1 ≈ 0.05
    @test uddotnp1 ≈ 1.0
    @test fhat ≈ 1.0

    azi0 = 0.1
    omega0_hz = 0.25
    dt = 0.2
    torque = 1.5
    azi1, omega1_hz, omegadot1_hz, applied_torque =
        OWENS.updateRotorRotation(1.0, 0.0, 0.0, 2.0, -0.5, azi0, omega0_hz, 0.0, dt)

    omega0 = omega0_hz * 2pi
    expected_alpha = torque
    @test applied_torque ≈ torque
    @test omegadot1_hz ≈ expected_alpha / 2pi
    @test omega1_hz ≈ (omega0 + 0.5 * dt * expected_alpha) / 2pi
    @test azi1 ≈ azi0 + dt * omega0 + 0.25 * dt^2 * expected_alpha
end

@testset "Time interpolation and residual helpers" begin
    omega, omegadot, terminate =
        OWENS.omegaSpecCheck(0.75, [0.0, 0.5, 1.0, 1.5], fill(0.2, 4), 0.1)
    @test omega ≈ 0.2
    @test abs(omegadot) < 1e-12
    @test terminate === false

    omega, omegadot, terminate = OWENS.omegaSpecCheck(2.0, [0.0, 1.0], [0.2, 0.2], 0.1)
    @test omega == 0.0
    @test omegadot == 0.0
    @test terminate === true

    force_history = [
        0.0 10.0 20.0
        5.0 15.0 25.0
    ]
    forces, dofs, hydro, mooring, platform, reaction =
        OWENS.externalForcing(0.5, [0.0, 1.0, 2.0], force_history, [2, 5])
    @test forces ≈ [5.0, 10.0]
    @test dofs == [2, 5]
    @test isnothing(hydro)
    @test isnothing(mooring)
    @test isnothing(platform)
    @test isnothing(reaction)

    histories = [
        1.0 2.0 3.0
        1.0 4.0 9.0
    ]
    ts = [0.0, 1.0, 2.0]
    @test OWENS.extrap_pred_vals(histories, ts, 3.0, 0) ≈ [3.0, 9.0]
    @test OWENS.extrap_pred_vals(histories, ts, 3.0, 1) ≈ [4.0, 14.0]
    @test OWENS.extrap_pred_vals(histories, ts, 3.0, 2) ≈ [4.0, 16.0]
    @test_throws ErrorException OWENS.extrap_pred_vals(histories, ts, 3.0, 3)

    residual = OWENS.calcHydroResidual(
        [0.2, -0.1],
        [100.0, 50.0],
        [20.0, -10.0],
        [1.5, 0.25, 0.3, -0.2],
        100.0,
    )
    @test residual ≈ Float32[0.3, -0.15, 0.1, -0.1]
end

@testset "Rotation and frame utilities" begin
    @test OWENS.transMat(0.0, 0.0, 0.0) == I(3)

    transform = OWENS.transMat(0.1, -0.2, 0.3)
    @test transform' * transform ≈ I(3) atol = 1e-6
    @test det(transform) ≈ 1.0 atol = 1e-6

    z90 = OWENS.createSingleRotationDCM(90.0, 3)
    @test z90 * [1.0, 0.0, 0.0] ≈ [0.0, -1.0, 0.0] atol = 1e-12
    @test OWENS.createGeneralTransformationMatrix([90.0], [3]) ≈ z90
    @test_throws ErrorException OWENS.createSingleRotationDCM(10.0, 4)

    values = [1.0, 0.0, 0.0, 0.0, 2.0, 0.0]
    converted = OWENS.frame_convert(values, z90)
    @test converted[1:3] ≈ [0.0, -1.0, 0.0] atol = 1e-12
    @test converted[4:6] ≈ [2.0, 0.0, 0.0] atol = 1e-12

    h1, h2, h3 = OWENS.rigidBodyRotation([1.0], [0.0], [0.0], [90.0], [3])
    @test only(h1) ≈ 0.0 atol = 1e-12
    @test only(h2) ≈ 1.0 atol = 1e-12
    @test only(h3) ≈ 0.0 atol = 1e-12

    hub = OWENS.calcHubRotMat([0.0, 0.0, 0.0], pi / 2)
    @test hub ≈ [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0] atol = 1e-12

    psi, theta = OWENS.calculatePsiTheta([1.0, 0.0, 0.0])
    @test psi ≈ 0.0
    @test abs(theta) < 1e-10

    psi, theta = OWENS.calculatePsiTheta([0.0, 1.0, 0.0])
    @test psi ≈ 90.0
    @test abs(theta) < 1e-10

    psi, theta = OWENS.calculatePsiTheta([0.0, 0.0, 1.0])
    @test psi ≈ 0.0
    @test theta ≈ -90.0 atol = 1e-5
end

@testset "Interpolation and boundary-condition guardrails" begin
    x = collect(0.0:1.0:4.0)
    y = x .^ 2
    @test all(isfinite, OWENS.safeakima(x, y, [1.5, 2.5]))
    @test_throws OverflowError OWENS.safeakima(x, y, [-0.2])
    @test_logs (:warn, r"Extrapolating on akima spline") OWENS.safeakima(
        x,
        y,
        [-0.2];
        extrapolate = true,
    )

    @test OWENS.setBCs([], [], 4, 6) == []
    bcs = OWENS.setBCs([2, 5], [1, 3], 4, 6)
    @test size(bcs) == (16, 3)
    @test bcs[1:6, :] == hcat(fill(1, 6), collect(1:6), zeros(Int, 6))
    @test bcs[7:12, :] == hcat(fill(3, 6), collect(1:6), zeros(Int, 6))
    @test [2 2 0; 4 2 0; 2 5 0; 4 5 0] == bcs[13:16, :]
    @test length(unique(eachrow(bcs[:, 1:2]))) == size(bcs, 1)
end

@testset "Fatigue damage monotonicity" begin
    t = range(0.0, 4pi; length = 121)
    stress = 10.0e6 .* sin.(t)
    larger_stress = 15.0e6 .* sin.(t)
    sn_stress = [0.0, 30.0e6, 60.0e6, 90.0e6]
    sn_log_cycles = [9.0, 7.0, 5.0, 4.0]
    ultimate_strength = 200.0e6

    damage = OWENS.fatigue_damage(stress, sn_stress, sn_log_cycles, ultimate_strength)
    larger_damage =
        OWENS.fatigue_damage(larger_stress, sn_stress, sn_log_cycles, ultimate_strength)
    no_mean_correction = OWENS.fatigue_damage(
        stress,
        sn_stress,
        sn_log_cycles,
        ultimate_strength;
        mean_correction = false,
    )

    @test isfinite(damage)
    @test damage > 0.0
    @test larger_damage > damage
    @test isfinite(no_mean_correction)
    @test no_mean_correction > 0.0
    @test isnan(OWENS.fatigue_damage([NaN], sn_stress, sn_log_cycles, ultimate_strength))
end
