using Test
import OWENS
import OWENSFEA

@testset "AC/DMS aerodynamic pitching moment mapping" begin
    mesh = OWENSFEA.Mesh(
        [1, 2],
        1,
        2,
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 1.0],
        [1],
        [1 2],
        [0],
        [2],
        [0.5 0.5],
        [1 2],
        [1 1],
        0,
        1,
        zeros(3),
        zeros(3),
    )
    el = OWENSFEA.El(
        [(twist = [0.0, 0.0],)],
        [1.0],
        [0.0],
        [0.0],
        [0.0],
        [1.0],
    )

    function advance_stub(t; azi, alwaysrecalc)
        nblade = 1
        nslices = 3
        nsteps = 1
        blade_slice_step = zeros(nblade, nslices, nsteps)
        slice_step = zeros(nslices, nsteps)
        step = zeros(nsteps)
        return (
            0.0,
            blade_slice_step,
            blade_slice_step,
            blade_slice_step,
            blade_slice_step,
            slice_step,
            slice_step,
            blade_slice_step,
            slice_step,
            step,
            nsteps,
            step,
            step,
            step,
            step,
            step,
            step,
            step,
            step,
            step,
            [0.25, 0.5, 0.75],
            zeros(nblade, nslices),
            blade_slice_step,
            blade_slice_step,
            step,
            slice_step,
            slice_step,
            slice_step,
            slice_step,
            zeros(nslices, nsteps, 3),
            slice_step,
            fill(12.0, nblade, nslices, nsteps),
        )
    end

    forces, dofs = OWENS.mapACDMS(0.0, 0.0, mesh, el, advance_stub)

    @test dofs == collect(1:12)
    @test forces[1:3, 1] == zeros(3)
    @test forces[10, 1] ≈ 6.0 atol = 1e-12

    buoyancy_only_result = (ntuple(_ -> nothing, 28)..., zeros(1, 3, 1, 3))
    @test OWENS._optional_acdms_m25(buoyancy_only_result, (1, 3, 1)) === nothing

    m25 = reshape(collect(1.0:3.0), 1, 3, 1)
    with_m25_result = (ntuple(_ -> nothing, 30)..., m25)
    @test OWENS._optional_acdms_m25(with_m25_result, size(m25)) === m25
end
