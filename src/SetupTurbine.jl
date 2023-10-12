function setupOWENS(VAWTAero,path;
    rho = 1.225,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    shapeY = collect(LinRange(0,H,Nslices+1)),
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2),#shapeX_spline(shapeY)
    ifw=false,
    turbsim_filename="$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0,1.0],
    chord_scale = [1.0,1.0],
    thickness_scale = [1.0,1.0],
    Ht=2.0,
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    nselem = 5,
    strut_mountpointbot = 0.01,
    strut_mountpointtop = 0.01,
    joint_type = 2,
    c_mount_ratio = 0.05,
    angularOffset = -pi/2,
    AModel="DMS",
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = "Darrieus") #Darrieus, H-VAWT, ARCUS

    Nbld = B
    H = maximum(shapeY) #m,
    R = maximum(shapeX) #m,
    omega = RPM / 60 * 2 * pi
    tsr = omega*R/Vinf

    shapeX_spline = FLOWMath.Akima(shapeY, shapeX)
    bladelen = sum(sqrt.((shapeX[2:end].-shapeX[1:end-1]).^2 .+ (shapeY[2:end].-shapeY[1:end-1]).^2 ))
    # println("bladelen 161.06: $bladelen")
    h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
    h_elem = (shapeY[2:end] - shapeY[1:end-1])
    h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

    RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, H, atol=1e-10)
    RefArea = RefArea_half*2
    delta_xs = shapeX[2:end] - shapeX[1:end-1]
    delta_zs = shapeY[2:end] - shapeY[1:end-1]

    delta3D = atan.(delta_xs./delta_zs)

    #########################################
    ### Set up mesh
    #########################################
    if meshtype == "ARCUS" #TODO, for all of these propogate the AeroDyn additional output requirements
        mymesh,myort,myjoint = OWENS.create_arcus_mesh(;Ht,
            Hb = H, #blade height
            R, # m bade radius
            nblade = Nbld,
            ntelem, #tower elements
            nbelem, #blade elements
            ncelem,
            c_mount_ratio,
            bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
            bshapez = shapeY,
            joint_type, #hinged about y axis
            cables_connected_to_blade_base,
            angularOffset) #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    elseif meshtype == "Darrieus"
        mymesh,myort,myjoint = OWENS.create_mesh_struts(;Ht,
            Hb = H, #blade height
            R, # m bade radius
            nblade = Nbld,
            ntelem, #tower elements
            nbelem, #blade elements
            nselem,
            strut_mountpointbot, # This puts struts at top and bottom
            strut_mountpointtop, # This puts struts at top and bottom
            bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
            bshapez = shapeY,
            bshapey = zeros(nbelem+1), # but magnitude for this is relevant
            angularOffset) #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        
    elseif meshtype == "H-VAWT"
        # mymesh,myort,myjoint = OWENS.create_mesh_HVAWT(;Ht = 15.0,
        #     Hb = 147.148-15.0, #blade height
        #     R = 54.014, # m bade radius
        #     nblade = 3,
        #     ntelem = 30, #tower elements
        #     nbelem = 30, #blade elements
        #     nselem = 4,
        #     strut_mountpoint = 0.01, # This puts struts at top and bottom
        #     bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        #     bshapez = zeros(nbelem+1),
        #     angularOffset = 0.0) #Blade shape, magnitude is irrelevant, scaled based on height and radius above 
    else #TODO unify with HAWT
        error("please choose a valid mesh type (Darrieus, H-VAWT, ARCUS)")
    end
    
    nTwrElem = Int(mymesh.meshSeg[1])+1

    # PyPlot.figure()
    # PyPlot.plot(mymesh.x,mymesh.z,"b-")
    #  for myi = 1:length(mymesh.x)
    #      PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
    #      PyPlot.draw()
    #      sleep(0.1)
    #  end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # # PyPlot.axis("equal")

    #########################################
    ### Set up Sectional Properties
    #########################################
    # println("Calculate/Set up sectional properties")
    #Tower
    if !isnothing(NuMad_geom_xlscsv_file_twr)
        numadIn_twr = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_twr)
    else
        n_web = 0
        n_stack = 2
        n_segments = 2
        span = [0.0, 6.607421057, 13.21484211, 19.82226317, 26.42968423, 33.03710529, 39.64452634, 46.2519474, 52.85936846, 59.46678951, 66.07421057, 72.68163163, 79.28905268, 85.89647374, 92.5038948, 99.11131586, 105.7187369, 112.326158, 118.933579, 125.5410001, 132.1484211]
        airfoil = ["circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular"]
        te_type = ["round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round"]
        twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        chord = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.25, 0.25, 0.25]
        xoffset = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        aerocenter = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
        stack_mat_types = [8, 2]
        stack_layers = [70 3; 70 3; 70 3; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03]
        segments = [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0]
        DPtypes = ["" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""]
        skin_seq = [Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2]); Seq([2, 1, 2]) Seq([2, 1, 2])]
        web_seq = Array{Seq, 2}(undef, length(twist_d),0) #can be any number of stack nums, so we have to make non-square containers
        web_dp = Array{Seq, 2}(undef, length(twist_d),0) #this is fixed size square, but it's easier to do it this way

        numadIn_twr = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end

    #Add the full path
    for (i,airfoil) in enumerate(numadIn_twr.airfoil)
        numadIn_twr.airfoil[i] = "$path/airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_twr)
        plyprops_twr = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e9, 1.0e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_twr = OWENS.plyproperties(names,plies)
    end

    twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = OWENS.getPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
    sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
    stiff_twr, mass_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

    #Blades
    if !isnothing(NuMad_geom_xlscsv_file_bld)
        numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_bld)
    else
        n_web = 1
        n_stack = 7
        n_segments = 12
        span = [0.0, 6.607, 13.215, 19.822, 26.43, 33.037, 39.645, 46.252, 52.859, 59.467, 66.074, 72.682, 79.289, 85.896, 92.504, 99.111, 105.719, 112.326, 118.934, 125.541, 132.148]
        airfoil = ["circular", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021", "NACA_0021"]
        te_type = ["round", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp"]
        twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        chord = [10.0, 10.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        xoffset = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        aerocenter = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
        stack_mat_types = [8, 2, 4, 8, 8, 8, 4]
        if isnothing(stack_layers_bld)
            stack_layers = [30.0 2.0 15.0 25.0 25.0 2.0 13.0; 15.0 2.0 10.0 13.0 11.0 2.0 11.0; 10.0 1.0 8.0 10.0 10.0 2.0 10.0; 8.0 1.0 6.0 9.0 10.0 1.0 9.0; 7.0 1.0 5.0 8.0 9.0 1.0 7.0; 6.0 1.0 4.0 8.0 9.0 1.0 6.0; 6.0 1.0 4.0 8.0 8.0 1.0 5.0; 6.0 1.0 4.0 7.0 7.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 7.0 8.0 3.0 5.0; 7.0 2.0 3.0 9.0 12.0 3.0 6.0; 10.0 3.0 4.0 11.0 15.0 3.0 6.0; 12.0 3.0 4.0 13.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 10.0 1.0 4.0 10.0 12.0 1.0 5.0]
        else
            stack_layers = stack_layers_bld
        end
        segments = [-1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0]
        DPtypes = ["" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""]
        skin_seq = [Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2])]
        web_seq = [Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]);;]
        web_dp = [Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]);;]

        numadIn_bld = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end
    for icol = 1:length(numadIn_bld.stack_layers[1,:])
        numadIn_bld.stack_layers[:,icol] .*= LinRange(stack_layers_scale[1],stack_layers_scale[2],length(numadIn_bld.chord))
    end
    numadIn_bld.chord .*= LinRange(chord_scale[1],chord_scale[2],length(numadIn_bld.chord))

    for (i,airfoil) in enumerate(numadIn_bld.airfoil)
        numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_bld)
        plyprops_bld = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_bld)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_bld = OWENS.plyproperties(names,plies)
    end

    # Get blade spanwise position
    bld1start = Int(mymesh.structuralNodeNumbers[1,1])
    bld1end = Int(mymesh.structuralNodeNumbers[1,end])
    spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

    if length(thickness_scale)==2
        yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_bld.span)))
    elseif length(thickness_scale)==length(numadIn_bld.span)
        yscale = thickness_scale
    end

    bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getPreCompOutput(numadIn_bld;yscale,plyprops = plyprops_bld)
    sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
    stiff_bld, mass_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

    #Struts
    # They are the same as the end properties of the blades

    # Combined Section Props
    Nremain = sum(Int,mymesh.meshSeg[Nbld+1+1:end]) #strut elements remain
    bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
    cable_secprop = sectionPropsArray_twr[end]
    sectionPropsArray = [fill(sectionPropsArray_twr[1],length(sectionPropsArray_twr));bldssecprops; fill(cable_secprop,Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
    # sectionPropsArray = fill(sectionPropsArray[end],length(sectionPropsArray))
    rotationalEffects = ones(mymesh.numEl)

    #store data in element object
    myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

    #########################################
    ### Set up aero forces
    #########################################
    # println("Initialize Aerodynamics")
    chord_spl = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,Nslices))
    VAWTAero.setupTurb(shapeX,shapeY,B,chord_spl,tsr,Vinf;AModel,DSModel,
    afname = "$path/airfoils/NACA_0021.dat", #TODO: map to the numad input
    ifw,
    turbsim_filename,
    ifw_libfile,
    ntheta,Nslices,rho,eta,RPI)

    aeroForces(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
    
    
    # Calculate mass breakout of each material
    mass_breakout_bld = get_material_mass(plyprops_bld,numadIn_bld)
    mass_breakout_blds = mass_breakout_bld.*length(mymesh.structuralNodeNumbers[:,1])
    mass_breakout_twr = get_material_mass(plyprops_twr,numadIn_twr;int_start=0.0,int_stop=Ht)

    return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,RefArea,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,RefArea,
    mass_breakout_blds,mass_breakout_twr
end


function get_material_mass(plyprops_in,numadIn;int_start=numadIn.span[1],int_stop=numadIn.span[end])
    # Get Relative contribution to mass by setting all but one of the materials to zero.
    mass_component_material = zeros(length(plyprops_in.names))
    for imat = 1:length(plyprops_in.names)
        # Initialize array since it is immutable
        plies = Array{Composites.Material}(undef,length(plyprops_in.names))

        # Fill it in, setting all the rho's to zero except the one matching imat
        for imat2 = 1:length(plyprops_in.names)
            if imat != imat2
                plies[imat2] = Composites.Material(plyprops_in.plies[imat2].e1,
                plyprops_in.plies[imat2].e2,
                plyprops_in.plies[imat2].g12,
                plyprops_in.plies[imat2].nu12,
                0.0, #rho
                plyprops_in.plies[imat2].xt,
                plyprops_in.plies[imat2].xc,
                plyprops_in.plies[imat2].yt,
                plyprops_in.plies[imat2].yc,
                plyprops_in.plies[imat2].s,
                plyprops_in.plies[imat2].t)
            else
                plies[imat2] = plyprops_in.plies[imat2]
            end
        end

        plyprops = OWENS.plyproperties(plyprops_in.names,plies)
        # Get the precomp output
        precompoutput,_,_,_,_ = OWENS.getPreCompOutput(numadIn;plyprops)
        mass_array = [precompoutput[iter].mass for iter=1:length(precompoutput)]
        # Spline and integrate that output across the span
        mass_spl = FLOWMath.Akima(numadIn.span,mass_array)
        mass_component_material[imat], error = QuadGK.quadgk(mass_spl, int_start, int_stop, atol=1e-10)

    end
    return mass_component_material
end

function setupOWENShawt(VAWTAero,path;
    rho = 1.225,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    hubR = 2.0,
    shapeY = collect(LinRange(0,H,Nslices+1)),
    shapeX = R.*(1.0.-4.0.*(shapeY/H.-.5).^2),#shapeX_spline(shapeY)
    ifw=false,
    turbsim_filename="$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0,1.0],
    chord_scale = [1.0,1.0],
    thickness_scale = [1.0,1.0],
    Ht=2.0,
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    joint_type = 2,
    c_mount_ratio = 0.05,
    AModel="DMS",
    DSModel="BV",
    RPI=true,
    biwing=false,
    hub_depth = 15.0, #Hub Beam Depth
    R_root = 10.0, # m biwing radius
    R_biwing = 30.0, # outer radius
    R_tip = 54.014, # outer radius
    nbelem_root = 30, #biwing elements, for each 
    nbelem_biwing = 30, #tip elements
    nbelem_tip = 30, #tip elements
    bshapex_root = LinRange(0.0,R_root,nbelem_root+1), #Blade shape, magnitude is relevant
    bshapez_root = zeros(nbelem_root+1), #Blade shape, magnitude is relevant
    bshapex_biwing_U = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_U = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_biwing_L = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_L = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_tip = LinRange(R_biwing,R_tip,nbelem_tip+1), #Blade shape, magnitude is relevant
    bshapez_tip = zeros(nbelem_tip+1), #Blade shape, magnitude is relevant
    AD15_ccw = false,
    angularOffset = 0.0
    )

    Nbld = B
    H = maximum(shapeY) #m,
    R = maximum(shapeX) #m,
    omega = RPM / 60 * 2 * pi
    tsr = omega*R/Vinf

    shapeX_spline = FLOWMath.Akima(shapeY, shapeX)
    bladelen = sum(sqrt.((shapeX[2:end].-shapeX[1:end-1]).^2 .+ (shapeY[2:end].-shapeY[1:end-1]).^2 ))
    # println("bladelen 161.06: $bladelen")
    h_frac = (shapeY[2:end] - shapeY[1:end-1])./shapeY[end];
    h_elem = (shapeY[2:end] - shapeY[1:end-1])
    h = (shapeY[2:end] + shapeY[1:end-1])/2.0;

    RefArea = pi*R^2
    #########################################
    ### Set up mesh
    #########################################
    # println("Create Mesh")
    if !biwing
        mymesh,myort,myjoint,bladeIdx,bladeElem = create_hawt_mesh(;hub_depth=Ht,
        tip_precone = H, #blade height
        R, # m bade radius
        hubR,
        nblade = Nbld,
        AD15_ccw,
        ntelem, #tower elements
        nbelem, #blade elements
        bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
        bshapez = shapeY,
        joint_type, #hinged about y axis
        angularOffset) #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    else
        mymesh,myort,myjoint,bladeIdx,bladeElem = create_hawt_biwing_mesh(;
            nblade = Nbld,
            hub_depth,
            R_root,
            R_biwing,
            R_tip,
            ntelem,
            nbelem_root,
            nbelem_biwing,
            nbelem_tip,
            bshapex_root,
            bshapez_root,
            bshapex_biwing_U,
            bshapez_biwing_U,
            bshapex_biwing_L,
            bshapez_biwing_L,
            bshapex_tip,
            bshapez_tip,
            joint_type,
            angularOffset = -pi/2)
    end

    nTwrElem = Int(mymesh.meshSeg[1])+1

    # PyPlot.figure()
    # PyPlot.plot(mymesh.x,mymesh.z,"b-")
    #  for myi = 1:length(mymesh.x)
    #      PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
    #      PyPlot.draw()
    #      sleep(0.1)
    #  end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.axis("equal")

    #########################################
    ### Set up Sectional Properties
    #########################################
    # println("Calculate/Set up sectional properties")
    #Tower
    if !isnothing(NuMad_geom_xlscsv_file_twr)
        numadIn_twr = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_twr)
    else
        n_web = 0
        n_stack = 2
        n_segments = 2
        span = [0.0, 6.607421057, 13.21484211, 19.82226317, 26.42968423, 33.03710529, 39.64452634, 46.2519474, 52.85936846, 59.46678951, 66.07421057, 72.68163163, 79.28905268, 85.89647374, 92.5038948, 99.11131586, 105.7187369, 112.326158, 118.933579, 125.5410001, 132.1484211]
        airfoil = ["circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular", "circular"]
        te_type = ["round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round", "round"]
        twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        chord = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.25, 0.25, 0.25]
        xoffset = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        aerocenter = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
        stack_mat_types = [8, 2]
        stack_layers = [70 3; 70 3; 70 3; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03; 30.303 3.03]
        segments = [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0]
        DPtypes = ["" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""; "" "" ""]
        skin_seq = [OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]); OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2])]
        web_seq = Array{OWENS.Seq, 2}(undef, length(twist_d),0) #can be any number of stack nums, so we have to make non-square containers
        web_dp = Array{OWENS.Seq, 2}(undef, length(twist_d),0) #this is fixed size square, but it's easier to do it this way

        numadIn_twr = OWENS.NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end

    #Add the full path
    for (i,airfoil) in enumerate(numadIn_twr.airfoil)
        numadIn_twr.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_twr)
        plyprops_twr = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e9, 1.0e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_twr = OWENS.plyproperties(names,plies)
    end

    twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = OWENS.getPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
    sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
    stiff_twr, mass_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

    #Blades
    if !isnothing(NuMad_geom_xlscsv_file_bld)
        numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_bld)
    else
        n_web = 1
        n_stack = 7
        n_segments = 12
        span = [0.0, 6.607, 13.215, 19.822, 26.43, 33.037, 39.645, 46.252, 52.859, 59.467, 66.074, 72.682, 79.289, 85.896, 92.504, 99.111, 105.719, 112.326, 118.934, 125.541, 132.148]
        airfoil = ["Cylinder1", "Cylinder1", "Cylinder1", "Cylinder2", "Cylinder2", "Cylinder2", "DU21_A17", "DU21_A17", "DU21_A17", "DU25_A17", "DU25_A17", "DU25_A17", "DU30_A17", "DU30_A17", "DU30_A17", "DU35_A17", "DU35_A17", "DU35_A17", "DU40_A17", "DU40_A17", "NACA64_A17"]
        te_type = ["round", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp", "sharp"]
        twist_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        chord = [10.0, 10.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        xoffset = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        aerocenter = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
        stack_mat_types = [8, 2, 4, 8, 8, 8, 4]
        if isnothing(stack_layers_bld)
            stack_layers = [30.0 2.0 15.0 25.0 25.0 2.0 13.0; 15.0 2.0 10.0 13.0 11.0 2.0 11.0; 10.0 1.0 8.0 10.0 10.0 2.0 10.0; 8.0 1.0 6.0 9.0 10.0 1.0 9.0; 7.0 1.0 5.0 8.0 9.0 1.0 7.0; 6.0 1.0 4.0 8.0 9.0 1.0 6.0; 6.0 1.0 4.0 8.0 8.0 1.0 5.0; 6.0 1.0 4.0 7.0 7.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 7.0 8.0 3.0 5.0; 7.0 2.0 3.0 9.0 12.0 3.0 6.0; 10.0 3.0 4.0 11.0 15.0 3.0 6.0; 12.0 3.0 4.0 13.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 10.0 1.0 4.0 10.0 12.0 1.0 5.0]
        else
            stack_layers = stack_layers_bld
        end
        segments = [-1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0]
        DPtypes = ["" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""]
        skin_seq = [OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2]); OWENS.Seq([2, 5, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 1, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 4, 2]) OWENS.Seq([2, 3, 2]) OWENS.Seq([2, 5, 2])]
        web_seq = [OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]); OWENS.Seq([6, 7, 6]);;]
        web_dp = [OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]); OWENS.Seq([9, 3, 3, 9]);;]

        numadIn_bld = OWENS.NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end
    for icol = 1:length(numadIn_bld.stack_layers[1,:])
        numadIn_bld.stack_layers[:,icol] .*= LinRange(stack_layers_scale[1],stack_layers_scale[2],length(numadIn_bld.chord))
    end
    numadIn_bld.chord .*= LinRange(chord_scale[1],chord_scale[2],length(numadIn_bld.chord))

    for (i,airfoil) in enumerate(numadIn_bld.airfoil)
        numadIn_bld.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_bld)
        plyprops_bld = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_bld)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_bld = OWENS.plyproperties(names,plies)
    end

    # Get blade spanwise position
    bld1start = Int(mymesh.structuralNodeNumbers[1,1])
    bld1end = Int(mymesh.structuralNodeNumbers[1,end])
    spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

    if length(thickness_scale)==2
        yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_bld.span)))
    elseif length(thickness_scale)==length(numadIn_bld.span)
        yscale = thickness_scale
    end

    bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getPreCompOutput(numadIn_bld;yscale,plyprops = plyprops_bld)
    sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
    stiff_bld, mass_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

    #Struts
    # They are the same as the end properties of the blades

    # Combined Section Props
    bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
    sectionPropsArray = [fill(sectionPropsArray_twr[1],length(sectionPropsArray_twr));bldssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
    # sectionPropsArray = fill(sectionPropsArray[end],length(sectionPropsArray))
    rotationalEffects = ones(mymesh.numEl)

    #store data in element object
    myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

    #########################################
    ### Set up aero forces
    #########################################
    # println("Initialize Aerodynamics")
    # chord_spl = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,Nslices))
    # VAWTAero.setupTurb(shapeX,shapeY,B,chord_spl,tsr,Vinf;AModel,DSModel,
    # afname = "$path/Airfoils/NACA_0021.dat", #TODO: map to the numad input
    # ifw,
    # turbsim_filename,
    # ifw_libfile,
    # ntheta,Nslices,rho,eta,RPI)

    aeroForces(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,VAWTAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
    
    
    # Calculate mass breakout of each material

    function get_material_mass(plyprops_in,numadIn;int_start=numadIn.span[1],int_stop=numadIn.span[end])
        # Get Relative contribution to mass by setting all but one of the materials to zero.
        mass_component_material = zeros(length(plyprops_in.names))
        for imat = 1:length(plyprops_in.names)
            # Initialize array since it is immutable
            plies = Array{Composites.Material}(undef,length(plyprops_in.names))

            # Fill it in, setting all the rho's to zero except the one matching imat
            for imat2 = 1:length(plyprops_in.names)
                if imat != imat2
                    plies[imat2] = Composites.Material(plyprops_in.plies[imat2].e1,
                    plyprops_in.plies[imat2].e2,
                    plyprops_in.plies[imat2].g12,
                    plyprops_in.plies[imat2].nu12,
                    0.0, #rho
                    plyprops_in.plies[imat2].xt,
                    plyprops_in.plies[imat2].xc,
                    plyprops_in.plies[imat2].yt,
                    plyprops_in.plies[imat2].yc,
                    plyprops_in.plies[imat2].s,
                    plyprops_in.plies[imat2].t)
                else
                    plies[imat2] = plyprops_in.plies[imat2]
                end
            end

            plyprops = OWENS.plyproperties(plyprops_in.names,plies)
            # Get the precomp output
            precompoutput,_,_,_,_ = OWENS.getPreCompOutput(numadIn;plyprops)
            mass_array = [precompoutput[iter].mass for iter=1:length(precompoutput)]
            # Spline and integrate that output across the span
            mass_spl = FLOWMath.Akima(numadIn.span,mass_array)
            mass_component_material[imat], error = QuadGK.quadgk(mass_spl, int_start, int_stop, atol=1e-10)

        end
        return mass_component_material
    end

    mass_breakout_bld = get_material_mass(plyprops_bld,numadIn_bld)
    mass_breakout_blds = mass_breakout_bld.*length(mymesh.structuralNodeNumbers[:,1])
    mass_breakout_twr = get_material_mass(plyprops_twr,numadIn_twr;int_start=0.0,int_stop=Ht)

    system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld;nblade=Nbld,struts=false)#,VTKmeshfilename="vtk_fvw/HAWT1")


    return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,RefArea,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,RefArea,
    mass_breakout_blds,mass_breakout_twr,bladeIdx,bladeElem,system,assembly,sections 
end