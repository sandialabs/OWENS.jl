using Test
import DelimitedFiles
import OWENS
import OWENSOpenFASTWrappers

function _ordered_body_nodes(indices)
    start_node, stop_node = indices
    return start_node <= stop_node ? (start_node:stop_node) : (stop_node:start_node)
end

function _read_openfast_output_table(filename)
    lines = split(replace(read(filename, String), Char(0) => ' '), '\n')
    header_index = findfirst(line -> startswith(strip(line), "Time"), lines)
    isnothing(header_index) &&
        throw(ArgumentError("OpenFAST output header was not found in $filename"))

    channels = split(strip(lines[header_index]))
    units = split(strip(lines[header_index+1]))
    rows = Vector{Vector{Float64}}()
    for line in lines[(header_index+2):end]
        tokens = split(strip(line))
        length(tokens) == length(channels) || continue
        values = tryparse.(Float64, tokens)
        any(isnothing, values) && continue
        push!(rows, Float64.(values))
    end
    isempty(rows) && throw(ArgumentError("OpenFAST output table has no numeric rows"))

    data = Matrix{Float64}(undef, length(rows), length(channels))
    for (irow, row) in enumerate(rows)
        data[irow, :] .= row
    end
    return (channels = channels, units = units, data = data)
end

function _write_helical_aerodyn_blade_files!(
    blade_files,
    mesh,
    body_nodes,
    height,
    n_blades,
    angular_offset,
    num_ad_nodes,
)
    for (ibody, filename) in enumerate(blade_files)
        nodes = _ordered_body_nodes(body_nodes[ibody, :])
        if ibody <= n_blades
            blade_angle = (ibody - 1) * 2pi / n_blades + angular_offset
            geometry = OWENS._aerodyn_blade_geometry_from_mesh(
                mesh.x[nodes],
                mesh.y[nodes],
                height,
                num_ad_nodes,
                blade_angle,
            )
            OWENSOpenFASTWrappers.writeADbladeFile(
                filename;
                NumBlNds = num_ad_nodes,
                BlSpn = geometry.BlSpn,
                BlCrvAC = geometry.BlCrvAC,
                BlSwpAC = geometry.BlSwpAC,
                BlCrvAng = geometry.BlCrvAng,
                BlTwist = geometry.BlTwist,
                BlChord = fill(0.25, num_ad_nodes),
                BlAFID = fill(1, num_ad_nodes),
            )
        else
            node_start = first(nodes)
            node_stop = last(nodes)
            body_length = sqrt(
                (mesh.x[node_stop] - mesh.x[node_start])^2 +
                (mesh.y[node_stop] - mesh.y[node_start])^2 +
                (mesh.z[node_stop] - mesh.z[node_start])^2,
            )
            OWENSOpenFASTWrappers.writeADbladeFile(
                filename;
                NumBlNds = num_ad_nodes,
                BlSpn = collect(range(0.0, body_length; length = num_ad_nodes)),
                BlCrvAC = zeros(num_ad_nodes),
                BlSwpAC = zeros(num_ad_nodes),
                BlCrvAng = zeros(num_ad_nodes),
                BlTwist = zeros(num_ad_nodes),
                BlChord = fill(0.15, num_ad_nodes),
                BlAFID = fill(1, num_ad_nodes),
            )
        end
    end
end

@testset "native AeroDyn initializes generated helical VAWT files" begin
    test_dir = @__DIR__
    mktempdir() do workdir
        height = 2.0
        radius = 1.5
        n_blades = 3
        n_struts_per_blade = 2
        omega = 4.0
        dt = 0.01
        tmax = 0.03
        hub_height = 2.0
        angular_offset = 0.0
        num_ad_nodes = 7
        shape_z = collect(range(0.0, height; length = 5))
        shape_x = fill(radius, length(shape_z))
        helical_y = [0.0, 0.3, 0.6, 0.3, 0.0]

        mesh, ort, _, body_nodes, body_elements = OWENS.create_mesh_struts(;
            Htwr_base = hub_height - height / 2,
            Htwr_blds = height / 2,
            Hbld = height,
            R = radius,
            AD15hubR = 0.2,
            nblade = n_blades,
            ntelem = 4,
            nbelem = 4,
            nselem = 2,
            strut_twr_mountpoint = [0.25, 0.75],
            strut_bld_mountpoint = [0.25, 0.75],
            bshapex = shape_x,
            bshapez = shape_z,
            bshapey = helical_y,
            angularOffset = angular_offset,
            AD15_ccw = true,
            verbosity = 0,
            connectBldTips2Twr = false,
        )
        @test size(body_nodes) == (9, 2)
        @test size(body_elements) == (9, 2)

        blade_files = [joinpath(workdir, "body_$ibody.dat") for ibody = 1:9]
        _write_helical_aerodyn_blade_files!(
            blade_files,
            mesh,
            body_nodes,
            height,
            n_blades,
            angular_offset,
            num_ad_nodes,
        )

        blade_1_table = DelimitedFiles.readdlm(blade_files[1]; skipstart = 6)
        @test size(blade_1_table) == (num_ad_nodes, 7)
        @test Float64.(blade_1_table[:, 1]) == [
            0.0,
            0.3333333333333333,
            0.6666666666666666,
            1.0,
            1.3333333333333333,
            1.6666666666666667,
            2.0,
        ]
        @test Float64.(blade_1_table[:, 2]) ≈ zeros(num_ad_nodes) atol = eps()
        @test Float64.(blade_1_table[:, 3]) ≈ [
            -0.0,
            -0.29999879999999995,
            -0.3,
            -0.6,
            -0.30000120000480013,
            -0.3,
            5.293955920339377e-23,
        ] atol = 1e-15
        @test Float64.(blade_1_table[:, 5]) ≈ [
            -11.309888400336881,
            0.0,
            4.407368333049016e-5,
            10.491521086014927,
            8.814752939883142e-5,
            4.407368333049016e-5,
            -11.309888400336883,
        ] atol = 1e-12

        ad_file = joinpath(workdir, "AD_helical.dat")
        ifw_file = joinpath(workdir, "IW_helical.dat")
        olaf_file = joinpath(workdir, "OLAF_helical.dat")
        rootname = joinpath(workdir, "helical_native")
        airfoil_file = joinpath(test_dir, "airfoils", "NACA_0018_AllRe.dat")

        OWENSOpenFASTWrappers.writeADinputFile(
            ad_file,
            blade_files,
            airfoil_file,
            olaf_file;
            rho = 1.225,
        )
        OWENSOpenFASTWrappers.writeOLAFfile(olaf_file; nNWPanel = 30, nFWPanels = 0)
        OWENSOpenFASTWrappers.writeIWfile(
            8.0,
            ifw_file;
            RefHt = hub_height,
            RefLength = radius,
            WindType = 1,
        )

        OWENSOpenFASTWrappers.setupTurb(
            nothing,
            ad_file,
            ifw_file,
            rootname,
            [shape_x],
            [shape_z],
            [n_blades],
            [hub_height],
            [mesh],
            [ort],
            [body_nodes],
            [body_elements];
            rho = 1.225,
            adi_dt = dt,
            adi_tmax = tmax,
            omega = [omega],
            adi_wrOuts = 1,
            adi_DT_Outs = dt,
            WrVTK = 0,
            numTurbines = 1,
            hubPos = [[0.0, 0.0, hub_height]],
            hubAngle = [[0.0, 0.0, 0.0]],
            adi_nstrut = [n_struts_per_blade],
            isHAWT = false,
        )
        try
            num_channels = OWENSOpenFASTWrappers.turbenv.num_channels
            initial_output = vec(OWENSOpenFASTWrappers.adiCalcOutput(0.0, num_channels))
            @test num_channels == 1974
            @test length(initial_output) == num_channels
            @test eltype(initial_output) === Float32
            @test all(isfinite, initial_output)

            for t = 0.0:dt:(tmax-dt)
                azimuth = omega * (t + dt)
                zero_state = zeros(mesh.numNodes * 6)
                OWENSOpenFASTWrappers.deformAD15(
                    [zero_state],
                    [zero_state],
                    [zero_state],
                    [azimuth],
                    [omega],
                    [0.0],
                    [[0.0, 0.0, hub_height]],
                    [[0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 0.0, 0.0, 0.0, omega]],
                    [zeros(6)],
                )
                OWENSOpenFASTWrappers.advanceAD15(t, [mesh], [azimuth])
            end
        finally
            OWENSOpenFASTWrappers.endTurb()
        end

        table = _read_openfast_output_table("$rootname.out")
        @test length(table.channels) == 1975
        @test size(table.data, 1) == 4
        @test "RtAeroPwr" in table.channels
        @test "RtSpeed" in table.channels
        @test "RtTSR" in table.channels
        @test "AB1N003Alpha" in table.channels
        @test "AB1N003Fx" in table.channels
        @test "AB2N003Alpha" in table.channels
        @test "AB3N003Alpha" in table.channels

        channel_value(name) = table.data[end, findfirst(==(name), table.channels)]
        @test table.data[end, 1] == 0.0
        response_channels = [
            "RtAeroPwr",
            "RtSpeed",
            "RtTSR",
            "AB1N003Alpha",
            "AB1N003Fx",
            "AB2N003Alpha",
            "AB3N003Alpha",
        ]
        response_values = channel_value.(response_channels)
        @test all(isfinite, response_values)
        @test channel_value("RtSpeed") > 0.0
        @test abs(channel_value("RtAeroPwr")) > 1.0
        @test abs(channel_value("AB1N003Fx")) > 1.0
    end
end
