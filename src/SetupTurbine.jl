function setupOWENS(OWENSAero,path;
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
    AD15hubR = 0.1,
    WindType=1,
    delta_t = 0.01,
    numTS = 100,
    adi_lib = "$(path)../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    windINPfilename="$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    NuMad_geom_xlscsv_file_strut = nothing,
    NuMad_mat_xlscsv_file_strut = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0,1.0],
    chord_scale = [1.0,1.0],
    thickness_scale = [1.0,1.0],
    Ht=2.0,
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    nselem = 5,
    strut_mountpointbot = 0.11,
    strut_mountpointtop = 0.11,
    joint_type = 2,
    c_mount_ratio = 0.05,
    angularOffset = -pi/2,
    AModel="DMS",
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    meshtype = "Darrieus") #Darrieus, H-VAWT, ARCUS

    if AModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    # Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
    # We do some intermediate calculations on the blade shape and angles

    Nstrutperbld = 2 #TODO: generalize and propogate

    Nbld = B
    H = maximum(shapeY) #m,
    R = maximum(shapeX) #m,
    omega = RPM / 60 * 2 * pi
    tsr = omega*R/Vinf
    
    nothing
    
    # Here we set up the mesh using one of the pre-made meshing functions.  For this case, there is a function for the ARCUS, as 
    # well as for towered VAWTs where you can have an arbitrary blade shape with connected struts, and if the blade tips touch
    # the tower, then you can tell it to constrain them to the tower thus allowing for both H-VAWT and Darrieus designs.
    
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
    elseif meshtype == "Darrieus" || meshtype == "H-VAWT"
        
        if meshtype == "Darrieus"
            connectBldTips2Twr = true
        else
            connectBldTips2Twr = false
        end
    
        mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.create_mesh_struts(;Ht,
            Hb = H, #blade height
            R, # m bade radius
            AD15hubR, #TODO: hook up with AD15 file generation
            nblade = Nbld,
            ntelem, #tower elements
            nbelem, #blade elements
            nselem,
            strut_twr_mountpointbot = strut_mountpointbot, # This puts struts at top and bottom
            strut_twr_mountpointtop = strut_mountpointtop, # This puts struts at top and bottom
            strut_bld_mountpointbot = strut_mountpointbot, # This puts struts at top and bottom
            strut_bld_mountpointtop = strut_mountpointtop, # This puts struts at top and bottom
            bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
            bshapez = shapeY,
            bshapey = zeros(nbelem+1), # but magnitude for this is relevant
            angularOffset, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
            AD15_ccw = true,
            verbosity=0, # 0 nothing, 1 basic, 2 lots: amount of printed information
            connectBldTips2Twr)
    else #TODO unify with HAWT
        error("please choose a valid mesh type (Darrieus, H-VAWT, ARCUS)")
    end
    
    nTwrElem = Int(mymesh.meshSeg[1])
    if contains(NuMad_mat_xlscsv_file_bld,"34m") #TODO: this is really odd, 
        nTwrElem = Int(mymesh.meshSeg[1])+1
    end
    
    nothing

    # Here is a way that you can visualize the nodal numbers of the mesh

    # PyPlot.figure()
    # PyPlot.plot(mymesh.x,mymesh.z,"b-")
    # for myi = 1:length(mymesh.x)
    #     PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
    #     PyPlot.draw()
    #     #sleep(0.1)
    # end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")

    # This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

    #########################################
    ### Set up Sectional Properties
    #########################################

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

    nothing

    # Here is where the material properties for the tower are either read in from the file, or directly input

    if !isnothing(NuMad_mat_xlscsv_file_twr)
        plyprops_twr = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e9, 1.0e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_twr = OWENS.plyproperties(names,plies)
    end

    # Then this is where precomp.jl is called to get first the precomp outputs, then formatting those into the OWENS format, and then in the GXBeam.jl format for if GXBeam is used as the structural solver.

    twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = OWENS.getOWENSPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
    sectionPropsArray_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
    stiff_twr, mass_twr = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

    nothing

    # For the blades, we repeat what was done for the tower, but also include some simple design options for scaling thicknesses,
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

    bld1start = Int(mymesh.structuralNodeNumbers[1,1]) #Get blade spanwise position
    bld1end = Int(mymesh.structuralNodeNumbers[1,end])
    spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

    if length(thickness_scale)==2
        yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_bld.span)))
    elseif length(thickness_scale)==length(numadIn_bld.span)
        yscale = thickness_scale
    end

    bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getOWENSPreCompOutput(numadIn_bld;yscale,plyprops = plyprops_bld)
    sectionPropsArray_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
    stiff_bld, mass_bld = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

    nothing

    # Similarly for the struts, we do what was done for the blades
    if !isnothing(NuMad_geom_xlscsv_file_strut)
        numadIn_strut = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file_strut)
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
        if isnothing(stack_layers_strut)
            stack_layers = [30.0 2.0 15.0 25.0 25.0 2.0 13.0; 15.0 2.0 10.0 13.0 11.0 2.0 11.0; 10.0 1.0 8.0 10.0 10.0 2.0 10.0; 8.0 1.0 6.0 9.0 10.0 1.0 9.0; 7.0 1.0 5.0 8.0 9.0 1.0 7.0; 6.0 1.0 4.0 8.0 9.0 1.0 6.0; 6.0 1.0 4.0 8.0 8.0 1.0 5.0; 6.0 1.0 4.0 7.0 7.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 8.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 1.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 6.0 6.0 2.0 5.0; 7.0 1.0 3.0 7.0 8.0 3.0 5.0; 7.0 2.0 3.0 9.0 12.0 3.0 6.0; 10.0 3.0 4.0 11.0 15.0 3.0 6.0; 12.0 3.0 4.0 13.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 12.0 3.0 4.0 15.0 15.0 3.0 6.0; 10.0 1.0 4.0 10.0 12.0 1.0 5.0]
        else
            stack_layers = stack_layers_strut
        end
        segments = [-1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0; -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0]
        DPtypes = ["" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""; "" "" "" "" "" "" "" "" "" "" "" "" ""]
        skin_seq = [Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]); Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2])]
        web_seq = [Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]);;]
        web_dp = [Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]);;]

        numadIn_strut = NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
    end
    for icol = 1:length(numadIn_strut.stack_layers[1,:])
        numadIn_strut.stack_layers[:,icol] .*= LinRange(stack_layers_scale[1],stack_layers_scale[2],length(numadIn_strut.chord))
    end
    numadIn_strut.chord .*= LinRange(chord_scale[1],chord_scale[2],length(numadIn_strut.chord))

    for (i,airfoil) in enumerate(numadIn_strut.airfoil)
        numadIn_strut.airfoil[i] = "$path/airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_strut)
        plyprops_strut = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_strut)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_strut = OWENS.plyproperties(names,plies)
    end

    #TODO: not straight struts
    spanpos = LinRange(0,1,nselem+1)#[0.0;cumsum(sqrt.(diff(mymesh.x[strut1start:strut1end]).^2 .+ diff(mymesh.z[strut1start:strut1end]).^2))]

    if length(thickness_scale)==2
        yscale = collect(LinRange(thickness_scale[1],thickness_scale[2],length(numadIn_strut.span)))
    elseif length(thickness_scale)==length(numadIn_strut.span)
        yscale = thickness_scale
    end

    strut_precompoutput,strut_precompinput,lam_U_strut,lam_L_strut,lam_W_strut = OWENS.getOWENSPreCompOutput(numadIn_strut;yscale,plyprops = plyprops_strut)
    sectionPropsArray_strut = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_strut,strut_precompoutput;precompinputs=strut_precompinput)
    stiff_strut, mass_strut = OWENS.getSectPropsFromOWENSPreComp(spanpos,numadIn_strut,strut_precompoutput;GX=true)

    nothing

    # Here we combine the section properties into an array matching the mesh elements
    bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
    strutssecprops = collect(Iterators.flatten(fill(sectionPropsArray_strut, Nstrutperbld*Nbld)))

    if meshtype == "ARCUS"
        cable_secprop = sectionPropsArray_twr[end]
        Nremain = sum(Int,mymesh.meshSeg[Nbld+1+1:end]) #strut elements remain
        sectionPropsArray = [fill(sectionPropsArray_twr[1],length(sectionPropsArray_twr));bldssecprops; fill(cable_secprop,Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

        # GXBeam sectional properties
        stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
        stiff_cables = fill(stiff_twr[end],Nremain)
        stiff_array = [stiff_twr; stiff_blds; stiff_cables]

        mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
        mass_cables = fill(mass_twr[end],Nremain)
        mass_array = [mass_twr; mass_blds; mass_cables]
    else
        sectionPropsArray = [sectionPropsArray_twr; bldssecprops; strutssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
        
        # GXBeam sectional properties
        stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
        stiff_struts = collect(Iterators.flatten(fill(stiff_strut, Nstrutperbld*Nbld)))
        stiff_array = [stiff_twr; stiff_blds; stiff_struts]

        mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
        mass_struts = collect(Iterators.flatten(fill(mass_strut, Nstrutperbld*Nbld)))
        mass_array = [mass_twr; mass_blds; mass_struts]
    end
    rotationalEffects = ones(mymesh.numEl) #TODO: non rotating tower, or rotating blades

    #store data in element object
    myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)
    system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_array,mass_array)

    nothing

    # Set up the OWENSAero aerodynamics if used
    if !AD15On
        #########################################
        ### Set up aero forces
        #########################################
        chord_spl = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,Nslices))

        T1 = round(Int,(5.8/H)*Nslices)
        T2 = round(Int,(11.1/H)*Nslices)
        T3 = round(Int,(29.0/H)*Nslices)
        T4 = round(Int,(34.7/H)*Nslices)
        airfoils = fill("$(path)/airfoils/NACA_0021.dat",Nslices)
        airfoils[T1:T4] .= "$(path)/airfoils/Sandia_001850.dat"
        
        chord = fill(1.22,Nslices) #TODO: link chord to numad and height as opposed to span
        chord[T1:T4] .= 1.07
        chord[T2:T3] .= 0.9191

        OWENSAero.setupTurb(shapeX,shapeY,B,chord,tsr,Vinf;AModel,DSModel,
        afname = airfoils, #TODO: map to the numad input
        rho,
        eta,
        ifw, #TODO: propogate WindType
        turbsim_filename = windINPfilename,
        ifw_libfile,
        tau = [1e-5,1e-5],
        ntheta,
        Nslices,
        RPI)

        aeroForcesACDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
        deformAeroACDMS = OWENSAero.deformTurb
    end
    nothing

    # Set up AeroDyn if used
    # Here we create AeroDyn the files, first by specifying the names, then by creating the files, TODO: hook up the direct sectionPropsArray_str
    # Then by initializing AeroDyn and grabbing the backend functionality with a function handle
    if AD15On
        ad_input_file="$path/ADInputFile_SingleTurbine2.dat"
        ifw_input_file="$path/IW2.dat"
        blade_filename="$path/blade2.dat"
        lower_strut_filename="$path/lower_arm2.dat"
        upper_strut_filename="$path/upper_arm2.dat"
        airfoil_filenames = "$path/airfoils/NACA_0018_AllRe.dat"
        OLAF_filename = "$path/OLAF2.dat"

        NumADBldNds = NumADStrutNds = 10 

        bldchord_spl = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,NumADBldNds))
        
        if meshtype == "ARCUS" 
            blade_filenames = [blade_filename for i=1:Nbld]
            blade_chords = [bldchord_spl for i=1:Nbld]
            blade_Nnodes = [NumADBldNds for i=1:Nbld]
        else
            blade_filenames = [[blade_filename for i=1:Nbld];[lower_strut_filename for i=1:Nbld];[upper_strut_filename for i=1:Nbld]]
            strutchord_spl = FLOWMath.akima(numadIn_strut.span./maximum(numadIn_strut.span), numadIn_strut.chord,LinRange(0,1,NumADStrutNds))
            blade_chords = [[bldchord_spl for i=1:Nbld];[strutchord_spl for i=1:Nbld];[strutchord_spl for i=1:Nbld]]
            blade_Nnodes = [[NumADBldNds for i=1:Nbld];[NumADStrutNds for i=1:Nbld];[NumADStrutNds for i=1:Nbld]]
        end

        OWENSOpenFASTWrappers.writeADinputFile(ad_input_file,blade_filenames,airfoil_filenames,OLAF_filename)

        NumADBody = length(AD15bldNdIdxRng[:,1])
        bld_len = zeros(NumADBody)
        for (iADBody,filename) in enumerate(blade_filenames)
            strt_idx = AD15bldNdIdxRng[iADBody,1]
            end_idx = AD15bldNdIdxRng[iADBody,2]
            if end_idx<strt_idx
                tmp_end = end_idx
                end_idx = strt_idx
                strt_idx = tmp_end
            end

            #Get the blade length
            x1 = mymesh.x[strt_idx]
            x2 = mymesh.x[end_idx]
            y1 = mymesh.y[strt_idx]
            y2 = mymesh.y[end_idx]
            z1 = mymesh.z[strt_idx]
            z2 = mymesh.z[end_idx]
            bld_len[iADBody] = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)

            #Get the blade shape
            ADshapeY = collect(LinRange(0,H,NumADBldNds))
            xmesh = mymesh.x[strt_idx:end_idx]
            ymesh = mymesh.y[strt_idx:end_idx]
            ADshapeX = sqrt.(xmesh.^2 .+ ymesh.^2)
            ADshapeX .-= ADshapeX[1] #get it starting at zero #TODO: make robust for blades that don't start at 0
            ADshapeXspl = FLOWMath.akima(LinRange(0,H,length(ADshapeX)),ADshapeX,ADshapeY)
            
            if iADBody<=Nbld #Note that the blades can be curved and are assumed to be oriented vertically
                BlSpn=ADshapeY
                BlCrvAC=ADshapeXspl
            else # while the arms/struts are assumed to be straight and are oriented by the mesh angle
                BlSpn=collect(LinRange(0,bld_len[iADBody],blade_Nnodes[iADBody]))
                BlCrvAC=zeros(blade_Nnodes[iADBody])
            end
            BlSwpAC=zeros(blade_Nnodes[iADBody])
            BlCrvAng=zeros(blade_Nnodes[iADBody])
            BlTwist=zeros(blade_Nnodes[iADBody])
            BlChord=blade_chords[iADBody]
            BlAFID=ones(Int,blade_Nnodes[iADBody])
            OWENSOpenFASTWrappers.writeADbladeFile(filename;NumBlNds=blade_Nnodes[iADBody],BlSpn,BlCrvAC,BlSwpAC,BlCrvAng,BlTwist,BlChord,BlAFID)
        end

        OWENSOpenFASTWrappers.writeOLAFfile(OLAF_filename;nNWPanel=200,nFWPanels=10)

        OWENSOpenFASTWrappers.writeIWfile(Vinf,ifw_input_file;windINPfilename)

        OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,[shapeX],[shapeY],[B],[Ht],[mymesh],[myort],[AD15bldNdIdxRng],[AD15bldElIdxRng];
                rho     = rho,
                adi_dt  = delta_t,
                adi_tmax= numTS*delta_t,
                omega   = [omega],
                adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
                adi_DT_Outs = delta_t,   # output frequency
                numTurbines = 1,
                refPos=[[0,0,0]],
                hubPos=[[0,0,0.0]],
                hubAngle=[[0,0,0]],
                nacPos=[[0,0,0]],
                adi_nstrut=[2],
                adi_debug=0,
                isHAWT = false     # true for HAWT, false for crossflow or VAWT
                )

        aeroForcesAD(t,azi) = OWENS.mapAD15(t,azi,[mymesh],OWENSOpenFASTWrappers.advanceAD15;alwaysrecalc=true,verbosity=1)
        deformAeroAD=OWENSOpenFASTWrappers.deformAD15
    end

    nothing

    # Calculate mass breakout of each material
    mass_breakout_bld = OWENS.get_material_mass(plyprops_bld,numadIn_bld)
    mass_breakout_blds = mass_breakout_bld.*length(mymesh.structuralNodeNumbers[:,1])
    mass_breakout_twr = OWENS.get_material_mass(plyprops_twr,numadIn_twr;int_start=0.0,int_stop=Ht)

    if AD15On
        return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
        stiff_twr, stiff_bld,bld_precompinput,
        bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
        twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForcesAD,deformAeroAD,
        mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng 
    else
        return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
        stiff_twr, stiff_bld,bld_precompinput,
        bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
        twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForcesACDMS,deformAeroACDMS,
        mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng 
    end
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

        plyprops = plyproperties(plyprops_in.names,plies)
        # Get the precomp output
        precompoutput,_,_,_,_ = getOWENSPreCompOutput(numadIn;plyprops)
        mass_array = [precompoutput[iter].mass for iter=1:length(precompoutput)]
        # Spline and integrate that output across the span
        mass_spl = FLOWMath.Akima(numadIn.span,mass_array)
        mass_component_material[imat], error = QuadGK.quadgk(mass_spl, int_start, int_stop, atol=1e-10)

    end
    return mass_component_material
end

function setupOWENShawt(OWENSAero,path;
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
    windINPfilename="$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
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

    nTwrElem = Int(mymesh.meshSeg[1])+2

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
        numadIn_twr = readNuMadGeomCSV(NuMad_geom_xlscsv_file_twr)
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
        numadIn_twr.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_twr)
        plyprops_twr = readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e9, 1.0e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_twr = plyproperties(names,plies)
    end

    twr_precompoutput,twr_precompinput,lam_U_twr,lam_L_twr,lam_W_twr = getOWENSPreCompOutput(numadIn_twr;plyprops = plyprops_twr)
    sectionPropsArray_twr = getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;precompinputs=twr_precompinput)
    stiff_twr, mass_twr = getSectPropsFromOWENSPreComp(LinRange(0,1,nTwrElem),numadIn_twr,twr_precompoutput;GX=true)

    #Blades
    if !isnothing(NuMad_geom_xlscsv_file_bld)
        numadIn_bld = readNuMadGeomCSV(NuMad_geom_xlscsv_file_bld)
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
        numadIn_bld.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_bld)
        plyprops_bld = readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_bld)
    else
        names = ["CLA_5500", "CBX_2400", "ETLX_2400", "Airex_C70_55", "EBX_2400_x10", "ETLX_2400_x10", "Airex_C70_55_x10", "CFP-baseline"]
        plies = [Composites.Material{Float64}(9.824e10, 5.102e9, 4.274e9, 0.3, 1540.0, 8.75634139e8, 5.92949102e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(1.4931e10, 1.4931e10, 2.389e10, 0.3, 1530.0, 4.55053962e8, 4.55053962e8, 1.0e8, 1.0e8, 1.0e8, 0.0008100000000000001), Composites.Material{Float64}(2.0333e10, 9.305e9, 4.756e9, 0.3, 1900.0, 5.30896289e8, 5.30896289e8, 1.0e8, 1.0e8, 1.0e8, 0.00066), Composites.Material{Float64}(4.5e7, 4.5e7, 2.2e7, 0.2, 59.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 0.001), Composites.Material{Float64}(9.824e11, 5.102e10, 4.274e10, 0.3, 15300.0, 4.55053962e9, 4.55053962e9, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.4931e11, 1.4931e11, 2.389e11, 0.3, 19000.0, 5.30896289e9, 5.30896289e9, 1.0e8, 1.0e8, 1.0e8, 8.0e-5), Composites.Material{Float64}(2.03335e11, 9.3051e10, 4.756e10, 0.2, 590.0, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 1.0e8, 7.000000000000001e-5), Composites.Material{Float64}(1.576e11, 9.1e9, 3.3e9, 0.263, 1600.0, 2.236e9, 1.528e9, 1.0e8, 1.0e8, 1.0e8, 0.00066)]
        plyprops_bld = plyproperties(names,plies)
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

    bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = getOWENSPreCompOutput(numadIn_bld;yscale,plyprops = plyprops_bld)
    sectionPropsArray_bld = getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
    stiff_bld, mass_bld = getSectPropsFromOWENSPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

    #Struts
    # They are the same as the end properties of the blades

    # Combined Section Props
    bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
    sectionPropsArray = [fill(sectionPropsArray_twr[1],length(sectionPropsArray_twr));bldssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
    # sectionPropsArray = fill(sectionPropsArray[end],length(sectionPropsArray))
    rotationalEffects = ones(mymesh.numEl)

    # GXBeam sectional properties
    stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
    stiff_array = [stiff_twr; stiff_blds]

    mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
    mass_array = [mass_twr; mass_blds]

    system, assembly, sections = owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_array,mass_array)

    #store data in element object
    myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

    #########################################
    ### Set up aero forces
    #########################################
    # println("Initialize Aerodynamics")
    # chord_spl = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,Nslices))
    # OWENSAero.setupTurb(shapeX,shapeY,B,chord_spl,tsr,Vinf;AModel,DSModel,
    # afname = "$path/Airfoils/NACA_0021.dat", #TODO: map to the numad input
    # ifw,
    # windINPfilename,
    # ifw_libfile,
    # ntheta,Nslices,rho,eta,RPI)

    aeroForces(t,azi) = mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
    
    
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

            plyprops = plyproperties(plyprops_in.names,plies)
            # Get the precomp output
            precompoutput,_,_,_,_ = getOWENSPreCompOutput(numadIn;plyprops)
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

    return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,
    mass_breakout_blds,mass_breakout_twr,bladeIdx,bladeElem,system,assembly,sections 
end