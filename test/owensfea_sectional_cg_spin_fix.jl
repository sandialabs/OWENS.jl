if !isdefined(@__MODULE__, :_owensfea_sectional_cg_steady_spin_fixed)
    function _owensfea_sectional_cg_steady_spin_fixed()
        L = 2.0
        nelem = 2
        section_ycm = 0.3
        spin_hz = 5.0
        b = 0.05
        h = 0.02
        area = b*h
        E = 2.1e17
        nu = 0.28
        G = E/(2*(1 + nu))
        rho = 7800.0
        Iyy = b*h^3/12
        Izz = b^3*h/12
        J = Iyy + Izz
        rhoA = rho*area
        pair(value) = fill(value, 2)

        section_props = Array{OWENSFEA.SectionPropsArray, 1}(undef, nelem)
        for i = 1:nelem
            section_props[i] = OWENSFEA.SectionPropsArray(
                pair(0.0), pair(0.0),
                fill(rhoA, 2),
                fill(E*Iyy, 2),
                fill(E*Izz, 2),
                fill(G*J, 2),
                fill(E*area, 2),
                fill(rho*Iyy, 2),
                fill(rho*Izz, 2),
                fill(rho*J, 2),
                pair(0.0), pair(section_ycm), pair(0.0), pair(0.0),
                pair(0.0), pair(0.0), pair(0.0), pair(0.0), pair(0.0), pair(0.0),
                pair(0.0), pair(0.0), pair(0.0), pair(0.0),
                nothing, nothing,
                pair(0.0), pair(0.0),
                fill(G*area, 2),
                fill(G*area, 2),
            )
        end

        x = collect(range(0.0, L, length=nelem + 1))
        conn = hcat(collect(1:nelem), collect(2:nelem + 1))
        mesh = OWENSFEA.Mesh(
            collect(1:nelem + 1),
            nelem,
            nelem + 1,
            x,
            zeros(nelem + 1),
            zeros(nelem + 1),
            collect(1:nelem),
            conn,
            zeros(Int, nelem),
            [nelem],
            zeros(1, 1),
            zeros(Int, 1, 1),
            zeros(Int, 1, 1),
        )
        el = OWENSFEA.El(
            section_props,
            fill(L/nelem, nelem),
            zeros(nelem),
            zeros(nelem),
            zeros(nelem),
            ones(nelem),
        )
        root_fixed = [
            1 1 0
            1 2 0
            1 3 0
            1 4 0
            1 5 0
            1 6 0
        ]
        feamodel = OWENSFEA.FEAModel(;
            analysisType="S",
            dataOutputFilename="none",
            joint=zeros(0, 8),
            pBC=root_fixed,
            numNodes=mesh.numNodes,
            nlOn=false,
            gravityOn=false,
            iterationType="LINEAR",
            maxNumLoadSteps=3,
        )
        el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
        displ0 = zeros(mesh.numNodes*6)
        _, _, success, reaction = redirect_stdout(devnull) do
            OWENSFEA.staticAnalysis(feamodel, mesh, el, displ0, spin_hz, spin_hz, el_storage)
        end
        moment_scale = abs(rhoA*(2*pi*spin_hz)^2*section_ycm*L^2)
        return success && isapprox(reaction[6], 0.0; atol=1e-6*moment_scale)
    end
end
