using LinearAlgebra
using OWENS
using OWENSPreComp
using Test

struct GXSectionNumadStub
    span
    aerocenter
    twist_d
    chord
end

struct GXSectionAirfoilStub
    xnode
    ynode
end

function synthetic_precomp_input()
    x_upper = collect(range(0.0, 1.0; length = 21))
    y_upper = 0.06 .* sin.(pi .* x_upper)
    x_lower = reverse(x_upper[1:(end-1)])
    y_lower = -0.06 .* sin.(pi .* x_lower)

    xnode = [x_upper; x_lower]
    ynode = [y_upper; y_lower]

    return OWENSPreComp.Input(
        1.0,
        0.0,
        0.0,
        0.25,
        xnode,
        ynode,
        [70.0e9, 12.0e9],
        [70.0e9, 12.0e9],
        [26.0e9, 4.5e9],
        [0.33, 0.30],
        [2700.0, 1600.0],
        [0.0, 0.45, 1.0],
        [1, 2],
        [1.0, 1.0, 1.0],
        [0.002, 0.0015, 0.001],
        [0.0, 45.0, -45.0],
        [1, 2, 2],
        [0.0, 1.0],
        [1],
        [1.0],
        [0.002],
        [0.0],
        [1],
        [0.35],
        [1],
        [1.0],
        [0.0015],
        [0.0],
        [1],
    )
end

@testset "GXBeam sectional property backend" begin
    pc_input = synthetic_precomp_input()

    section = OWENS.gxbeam_sectional_properties_from_precomp(pc_input; wns = 2)
    @test length(section.nodes) > 0
    @test length(section.elements) > 0
    @test size(section.stiffness) == (6, 6)
    @test size(section.mass) == (6, 6)
    @test isapprox(section.stiffness, section.stiffness'; rtol = 1e-10, atol = 1e-6)
    @test isapprox(section.mass, section.mass'; rtol = 1e-10, atol = 1e-10)
    @test minimum(eigvals(Symmetric(section.stiffness))) > 0.0
    @test section.mass[1, 1] > 0.0

    generalized_strain = [1.0e-5, 2.0e-6, -1.0e-6, 3.0e-5, 1.0e-4, -2.0e-4]
    strain_b, stress_b, strain_p, stress_p =
        OWENS.recover_gxbeam_sectional_strain(section, generalized_strain)
    @test size(strain_b) == (6, length(section.elements))
    @test size(stress_b) == (6, length(section.elements))
    @test size(strain_p) == (6, length(section.elements))
    @test size(stress_p) == (6, length(section.elements))
    @test all(isfinite, strain_b)
    @test all(isfinite, stress_b)
    @test all(isfinite, strain_p)
    @test all(isfinite, stress_p)

    numad = GXSectionNumadStub([0.0, 1.0], [0.25, 0.25], [0.0, 0.0], [1.0, 1.0])
    precomp_output = OWENSPreComp.Output(
        1.0e6,
        2.0e6,
        3.0e5,
        4.0e7,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        8.0,
        0.02,
        0.03,
        0.0,
        0.0,
        0.0,
    )
    stiff, mass, element_sections = OWENS.getSectPropsFromOWENSPreComp(
        [0.0, 0.5, 1.0],
        numad,
        [precomp_output, precomp_output];
        GX = true,
        precompinputs = [pc_input, pc_input],
        AddedMass_Coeff_Ca = 0.0,
        sectional_property_source = :gxbeam,
        return_sectional_mesh = true,
        gxbeam_section_kwargs = (; wns = 2),
    )

    @test length(stiff) == 2
    @test length(mass) == 2
    @test length(element_sections) == 2
    @test isapprox(stiff[1], section.stiffness)
    @test isapprox(mass[1], section.mass)
    @test element_sections[1] isa OWENS.GXBeamSectionalProperties

    doubled_history(value) = [value 2 * value; value 2 * value]
    component = OWENS.Component(;
        name = "gxbeam-test-component",
        elNumbers = [1, 2],
        gxBeamSectionalProperties = element_sections,
        e_x = doubled_history(generalized_strain[1]),
        e_y = doubled_history(generalized_strain[2]),
        e_z = doubled_history(generalized_strain[3]),
        k_x = doubled_history(generalized_strain[4]),
        k_y = doubled_history(generalized_strain[5]),
        k_z = doubled_history(generalized_strain[6]),
    )

    OWENS.populate_gxbeam_sectional_recovery!(component)
    @test size(component.gxBeamSectionalRecovery) == (2, 2)
    @test component.gxBeamSectionalRecovery[1, 1] isa OWENS.GXBeamSectionalRecovery
    @test isapprox(component.gxBeamSectionalRecovery[1, 1].beam_strain, strain_b)
    @test isapprox(component.gxBeamSectionalRecovery[1, 1].beam_stress, stress_b)
    @test isapprox(component.gxBeamSectionalRecovery[1, 1].ply_strain, strain_p)
    @test isapprox(component.gxBeamSectionalRecovery[1, 1].ply_stress, stress_p)
    @test !isapprox(
        component.gxBeamSectionalRecovery[1, 2].beam_stress,
        component.gxBeamSectionalRecovery[1, 1].beam_stress,
    )
end

@testset "GXBeam sectional property input validation" begin
    lower_surface_first = GXSectionAirfoilStub(
        [1.0, 0.7, 0.2, 0.4, 1.0],
        [0.0, -0.04, -0.02, 0.05, 0.0],
    )
    xaf, yaf = OWENS._precomp_airfoil_to_gxbeam(lower_surface_first)
    @test xaf == reverse(lower_surface_first.xnode)
    @test yaf == reverse(lower_surface_first.ynode)
    @test xaf[1] == 1.0
    @test xaf[end] == 1.0
    @test yaf[2] > 0.0

    upper_surface_first = GXSectionAirfoilStub(
        [1.0, 0.6, 0.2, 0.5, 1.0],
        [0.0, 0.05, 0.02, -0.04, 0.0],
    )
    xaf, yaf = OWENS._precomp_airfoil_to_gxbeam(upper_surface_first)
    @test xaf == upper_surface_first.xnode
    @test yaf == upper_surface_first.ynode

    missing_upper_surface = GXSectionAirfoilStub([1.0, 0.6, 0.0], [0.0, -0.04, 0.0])
    @test_throws ArgumentError OWENS._precomp_airfoil_to_gxbeam(missing_upper_surface)

    pc_input = synthetic_precomp_input()
    section = OWENS.gxbeam_sectional_properties_from_precomp(pc_input; wns = 2)
    strain_history = zeros(1, 2)

    missing_histories = OWENS.Component(;
        name = "missing-gxbeam-histories",
        gxBeamSectionalProperties = [section],
        e_x = strain_history,
    )
    @test_throws ArgumentError OWENS.recover_gxbeam_component_sectional_strain(
        missing_histories,
    )

    mismatched_histories = OWENS.Component(;
        name = "mismatched-gxbeam-histories",
        gxBeamSectionalProperties = [section],
        e_x = strain_history,
        e_y = zeros(2, 2),
        e_z = strain_history,
        k_x = strain_history,
        k_y = strain_history,
        k_z = strain_history,
    )
    @test_throws DimensionMismatch OWENS.recover_gxbeam_component_sectional_strain(
        mismatched_histories,
    )

    mismatched_sections = OWENS.Component(;
        name = "mismatched-gxbeam-sections",
        gxBeamSectionalProperties = [section, section],
        e_x = strain_history,
        e_y = strain_history,
        e_z = strain_history,
        k_x = strain_history,
        k_y = strain_history,
        k_z = strain_history,
    )
    @test_throws DimensionMismatch OWENS.recover_gxbeam_component_sectional_strain(
        mismatched_sections,
    )
end
