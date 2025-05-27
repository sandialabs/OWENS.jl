
function addSectionalPropertiesComponent(sectionPropsArray,stiff_array,mass_array,component,path,rho,AddedMass_Coeff_Ca;name=nothing)

    input_layup = component.input_layup
    input_materials = component.input_materials
    start_el_idx = component.elNumbers[1]
    end_el_idx = component.elNumbers[end]

    nElem = end_el_idx-start_el_idx+1+1
    subsection = nothing
    if !isnothing(input_layup)
        if contains(component.name,"tower")
            section=:tower
        elseif contains(component.name,"blade")
            section=:blade
        elseif contains(component.name,"strut")
            section=:struts
            if typeof(input_layup) == OrderedCollections.OrderedDict{Symbol, Any}
                if length(input_layup[:components][:struts])==1
                    subsection = 1
                else
                    subsection = parse(Int,component.name[end]) #This assumes that you have a numad file for each strut, and that you have 9 or fewer struts
                end
            else
                subsection = nothing
            end
        elseif contains(component.name,"intra_cable")
            section=:intra_cable
        elseif contains(component.name,"guy")
            section=:guy
        end
        numadIn = OWENS.readNuMadGeomCSV(input_layup;section,subsection) #note that this function will run either the numad input or the windio input depending on the input type
    else
        @error "Input data or numad file must be defined for this method"
    end

    #### Add the full path
    for (i,airfoil) in enumerate(numadIn.airfoil)
        numadIn.airfoil[i] = "$path/airfoils/$airfoil"
    end

    # Here is where the material properties for the tower are either read in from the file, or directly input
    if !isnothing(input_materials)
        plyprops = OWENS.readNuMadMaterialsCSV(input_materials)
    else
        @error "Input data or numad file must be defined for this method"
    end

    # Then this is where precomp.jl is called to get first the precomp outputs, then formatting those into the OWENS format, and then in the GXBeam.jl format for if GXBeam is used as the structural solver.
    precompoutput,precompinput,lam_U,lam_L,lam_W = OWENS.getOWENSPreCompOutput(numadIn;plyprops = plyprops)
    sectionPropsArray_component = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nElem),numadIn,precompoutput;precompinputs=precompinput,fluid_density=rho,AddedMass_Coeff_Ca)
    stiff, mass = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,nElem),numadIn,precompoutput;GX=true,precompinputs=precompinput,fluid_density=rho,AddedMass_Coeff_Ca)

    sectionPropsArray = [sectionPropsArray; sectionPropsArray_component]
    stiff_array = [stiff_array;stiff]
    mass_array = [mass_array;mass]


    component.lam_U = lam_U
    component.lam_L = lam_L
    component.lam_W = lam_W
    component.preCompInput = precompinput
    component.preCompOutput = precompoutput
    component.plyProps = plyprops
    component.nuMadIn = numadIn
    component.sectionProps = sectionPropsArray_component
    component.stiff_matrix = stiff
    component.mass_matrix = mass

    return sectionPropsArray,stiff_array,mass_array,component
end

# * `strut_twr_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the tower # This puts struts at top and bottom, as a fraction of the blade position
# * `strut_bld_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the blade # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
# NuMad_geom_xlscsv_file_strut = nothing, or sting for the filename, or array of strings matching length of Nstruts (length(strut_twr_mountpoint))
function setupOWENS(OWENSAero,path;
    rho = 1.225,
    mu = 1.7894e-5,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    shapeZ = collect(LinRange(0,H,Nslices+1)),
    shapeX = R.*(1.0.-4.0.*(shapeZ/H.-.5).^2),#shapeX_spline(shapeZ)
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    nselem = 5,
    shapeY = zeros(Nslices+1),
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
    NuMad_geom_xlscsv_file_intra_blade_cable = nothing,
    NuMad_mat_xlscsv_file_intra_blade_cable = nothing,
    NuMad_geom_xlscsv_file_guys = nothing,
    NuMad_mat_xlscsv_file_guys = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0,1.0],
    chord_scale = [1.0,1.0],
    thickness_scale = [1.0,1.0],
    Htwr_base = 2.0,
    Htwr_blds = H,
    strut_twr_mountpoint = [0.25,0.75],
    strut_bld_mountpoint = [0.25,0.75],
    joint_type = 2,
    c_mount_ratio = 0.05,
    angularOffset = -pi/2,
    AeroModel="DMS",
    DynamicStallModel="BV",
    RPI=true,
    Aero_AddedMass_Active = false,
    Aero_RotAccel_Active = false,
    Aero_Buoyancy_Active = false,
    cables_connected_to_blade_base = true,
    meshtype = "Darrieus",
    custommesh = nothing,
    AddedMass_Coeff_Ca = 0.0,
    verbosity=1,
    VTKmeshfilename = nothing,
    return_componentized = false)

    custom_mesh_outputs = []

    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    if meshtype == "Darrieus"
        connectBldTips2Twr = true
    else
        connectBldTips2Twr = false
    end

    if minimum(shapeZ)!=0
        @error "blade shapeZ must start at 0.0"
    end

    if AddedMass_Coeff_Ca>0.0
        centrifugal_force_flag = true
    else
        centrifugal_force_flag = false
    end

    # Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
    # We do some intermediate calculations on the blade shape and angles

    Nstrutperbld = length(strut_twr_mountpoint)

    if typeof(NuMad_geom_xlscsv_file_strut)==String || typeof(NuMad_geom_xlscsv_file_strut)==OrderedCollections.OrderedDict{Symbol, Any}
        NuMad_geom_xlscsv_file_strut = fill(NuMad_geom_xlscsv_file_strut,Nstrutperbld)
    end
    

    Nbld = B
    H = maximum(shapeZ) #m,
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
   
    # if isnothing(custommesh)
    #     error("this function requires a custom mesh be input")
    # end

    mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs = create_mesh_struts(;Htwr_base,
    Htwr_blds,
    Hbld = H, #blade height
    R, # m bade radius
    AD15hubR, #TODO: hook up with AD15 file generation
    nblade = Nbld,
    ntelem, #tower elements
    nbelem, #blade elements
    nselem,
    strut_twr_mountpoint,
    strut_bld_mountpoint,
    bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = shapeZ,
    bshapey = shapeY, # but magnitude for this is relevant
    angularOffset, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = true,
    verbosity, # 0 nothing, 1 basic, 2 lots: amount of printed information)
    connectBldTips2Twr,
    
    )

    nTwrElem = Int(mymesh.meshSeg[1])+1
    # try
    #     if contains(NuMad_mat_xlscsv_file_bld,"34m") || meshtype == "ARCUS" #TODO: this is really odd, 
    #         nTwrElem = Int(mymesh.meshSeg[1])+1
    #     end
    # catch
    #     nTwrElem = Int(mymesh.meshSeg[1])
    # end
    
    nothing

    # Here is a way that you can visualize the nodal numbers of the mesh

    # PyPlot.figure()
    # for icon = 1:length(mymesh.conn[:,1])
    #     idx1 = mymesh.conn[icon,1]
    #     idx2 = mymesh.conn[icon,2]
    #     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    #     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    #     # sleep(0.1)
    # end
    
    # for ijoint = 1:length(myjoint[:,1])
    #     idx2 = Int(myjoint[ijoint,2])
    #     idx1 = Int(myjoint[ijoint,3])
    #     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    #     # sleep(0.1)
    # end

    # Unpack to get the element numbers for the extra components
    components,s2t_idx,intra_blade_cable_startidx_el,intra_blade_cable_endidx_el,topel_idx,s2b_idx = custom_mesh_outputs

    # This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

    #########################################
    ### Set up Sectional Properties
    #########################################
    sectionPropsArray = []
    stiff_array = []
    mass_array = []
    rotationalEffects = ones(mymesh.numEl)

    for icomponent = 1:size(components)[1]
        if contains(components[icomponent].name,"tower")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_twr
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_twr
        elseif contains(components[icomponent].name,"blade")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_bld
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_bld
        elseif contains(components[icomponent].name,"strut")
            istrut = parse(Int,components[icomponent].name[end]) #This assumes that you have a numad file for each strut, and that you have 9 or fewer struts
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_strut[istrut]
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_strut
        elseif contains(components[icomponent].name,"intra_cable")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_intra_blade_cable
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_intra_blade_cable
        elseif contains(components[icomponent].name,"guy")
            components[icomponent].input_layup = NuMad_geom_xlscsv_file_guys
            components[icomponent].input_materials = NuMad_mat_xlscsv_file_guys
            rotationalEffects[components[icomponent].elNumbers] .= 0.0 #turn rotational effects off for guy wires
        end

        if verbosity>1
            println("Calculating Sectional Properties for Component $icomponent $(components[icomponent].name)")
        elseif verbosity>2
            println("Using for the layup, $(components[icomponent].input_layup) and for the materials, $(components[icomponent].input_materials)")
        end

        sectionPropsArray,stiff_array,mass_array,components[icomponent] = addSectionalPropertiesComponent(sectionPropsArray,stiff_array,mass_array,components[icomponent],path,rho,AddedMass_Coeff_Ca;name=nothing)

        components[icomponent].mass = OWENS.get_material_mass(components[icomponent].plyProps,components[icomponent].nuMadIn)
        components[icomponent].cost = ["$name $(components[icomponent].mass[i]) kg, $(components[icomponent].plyProps.costs[i]) \$/kg: \$$(components[icomponent].mass[i]*components[icomponent].plyProps.costs[i])" for (i,name) in enumerate(components[icomponent].plyProps.names)]
    end

    if verbosity>1
        println("\nTotal Mass: $(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]])) kg")
        println("Total Material Cost: \$$(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]]))")
    end

    # # Plot the rotational effects' status
    # for iel = 1:length(mymesh.conn[:,1])
    #     if rotationalEffects[iel] == 1
    #         idx1 = mymesh.conn[iel,1]
    #         idx2 = mymesh.conn[iel,2]
    #         PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],linewidth=10.0,"b.-")
    #         # PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    #         # sleep(0.1)
    #     end
    # end

    if length(sectionPropsArray)<mymesh.numEl
        @warn "There are more mesh elements $(mymesh.numEl) than sectional properties $(length(sectionPropsArray)), applying the last element's sectional properties to the remaining"
        n_diff = mymesh.numEl - length(sectionPropsArray)
        sectionPropsArray = [sectionPropsArray; fill(sectionPropsArray[end],n_diff)]
        stiff_array = [stiff_array;fill(stiff_array[end],n_diff)]
        mass_array = [mass_array;fill(mass_array[end],n_diff)]
    end

    #### store data in element object
    myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

    system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_array,mass_array;VTKmeshfilename)


    function calcMass(sectionPropsArray,myort)
        mass = 0.0
        lens = zeros(length(sectionPropsArray))
        rhoAs = zeros(length(sectionPropsArray))
        for (i,sectionProp) in enumerate(sectionPropsArray)
            lenEl = 0
            try
                lenEl = myort.Length[i]
            catch
                lenEl = myort.Length[i-1]
            end
            lens[i] = lenEl
            rhoA = sectionProp.rhoA[1]
            rhoAs[i] = rhoA
            mass += lenEl*rhoA
        end
        span_array = cumsum(lens)
        spl_rhoa = FLOWMath.Akima(span_array,rhoAs)
        mass2, error = QuadGK.quadgk(spl_rhoa, span_array[1], span_array[end], atol=1e-10)
        return mass, mass2
    end

    turb_masskg1,turb_masskg = calcMass(myel.props,myort)
    if verbosity>1
        println("\nMass of Turbine: $turb_masskg kg, $turb_masskg1 kg")
    end

    # Set up the OWENSAero aerodynamics if used
    blade_component_index = findfirst(s -> s.name[1:end-1] == "blade", components) #This assumes 9 or fewer blades, not including struts
    numadIn_bld = components[blade_component_index].nuMadIn

    strut_component_index = findall(s -> contains(s.name, "strut"), components) #This assumes 9 or fewer blades, not including struts
    numadIn_strut = [components[strut_component_index1].nuMadIn for strut_component_index1 = strut_component_index[1:Nbld:end]]
    if !AD15On
        #########################################
        ### Set up aero forces
        #########################################
        #### translate from blade span to blade height between the numad definition and the vertical slice positions
        #### First get the angles from the overall geometry npoints and go to the numad npoints
        delta_xs = shapeX[2:end] - shapeX[1:end-1]
        delta_zs = shapeZ[2:end] - shapeZ[1:end-1]
        delta3D = atan.(delta_xs./delta_zs)
        delta3D_spl = OWENS.safeakima(shapeZ[1:end-1]./maximum(shapeZ[1:end-1]), delta3D,LinRange(0,1,length(numadIn_bld.span)-1))
        #### now convert the numad span to a height
        bld_height_numad = cumsum(diff(numadIn_bld.span).*(1.0.-abs.(sin.(delta3D_spl))))
        bld_height_numad_unit = bld_height_numad./maximum(bld_height_numad)
        #### now we can use it to access the numad data 
        chord = OWENS.safeakima(bld_height_numad_unit, numadIn_bld.chord,LinRange(bld_height_numad_unit[1],1,Nslices))
        airfoils = fill("nothing",Nslices)

        twist = OWENS.safeakima(bld_height_numad_unit, numadIn_bld.twist_d.*pi/180,LinRange(bld_height_numad_unit[1],1,Nslices))

        # Discretely assign the airfoils
        for (iheight_numad,height_numad) in enumerate(bld_height_numad_unit)
            for (iheight,height_slices) in enumerate(collect(LinRange(0,1,Nslices)))
                if airfoils[iheight]=="nothing" && height_slices<=height_numad
                    airfoils[iheight] = "$(numadIn_bld.airfoil[iheight_numad]).dat"
                end
            end
        end
        
        # Map the element wise mass to the input aero shape
        mass_bld = []
        for component in components
            if contains(component.name,"blade1")
                mass_bld = component.mass_matrix
                break
            end
        end

        rhoA_el = [mass_bld[i][1,1] for i = 1:length(mass_bld)]
        if AD15bldNdIdxRng[1,2]<AD15bldNdIdxRng[1,1]
            bld_z_node = mymesh.z[AD15bldNdIdxRng[1,2]:AD15bldNdIdxRng[1,1]]
        else
            bld_z_node = mymesh.z[AD15bldNdIdxRng[1,1]:AD15bldNdIdxRng[1,2]]
        end
        bld_z_el = bld_z_node[1:end-1] .+ (bld_z_node[2:end].-bld_z_node[1:end-1])./2
        rhoA_in = FLOWMath.akima(bld_z_el,rhoA_el,shapeZ)

        OWENSAero.setupTurb(shapeX,shapeZ,B,chord,tsr,Vinf;AeroModel,DynamicStallModel,
        afname = airfoils,
        bld_y = shapeY,
        rho,
        twist, #TODO: verify twist is in same direction
        mu,
        eta,
        ifw, #TODO: propogate WindType
        turbsim_filename = windINPfilename,
        ifw_libfile,
        tau = [1e-5,1e-5],
        Aero_AddedMass_Active,
        Aero_RotAccel_Active,
        Aero_Buoyancy_Active,
        centrifugal_force_flag,
        ntheta,
        Nslices,
        RPI,
        rhoA_in)

        aeroForcesACDMS(t,azi) = OWENS.mapACDMS(t,azi,mymesh,myel,OWENSAero.AdvanceTurbineInterpolate;alwaysrecalc=true)
        deformAeroACDMS = OWENSAero.deformTurb
    end
    nothing

    # Set up AeroDyn if used
    # Here we create AeroDyn the files, first by specifying the names, then by creating the files, TODO: hook up the direct sectionPropsArray_str
    # Then by initializing AeroDyn and grabbing the backend functionality with a function handle
    if AD15On
        ad_input_file = "$path/ADInputFile_SingleTurbine2.dat"
        ifw_input_file = "$path/IW2.dat"
        OLAF_filename = "$path/OLAF2.dat"

        NumADBldNds = NumADStrutNds = 10 

        bldchord_spl = OWENS.safeakima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,NumADBldNds))

        # Discretely assign the blade airfoils based on the next closest neighbor
        bld_airfoil_filenames = fill("nothing",NumADBldNds) #TODO: cable drag?
        for (ispan_numad,span_numad) in enumerate(numadIn_bld.span./maximum(numadIn_bld.span))
            for (ispan,span_slices) in enumerate(collect(LinRange(0,1,NumADBldNds)))
                if bld_airfoil_filenames[ispan]=="nothing" && span_slices<=span_numad
                    bld_airfoil_filenames[ispan] = "$(numadIn_bld.airfoil[ispan_numad]).dat"
                end
            end
        end
        
        if meshtype == "ARCUS" 
            blade_filenames = ["$path/blade$i.dat" for i=1:Nbld]
            blade_chords = [bldchord_spl for i=1:Nbld]
            blade_Nnodes = [NumADBldNds for i=1:Nbld]
            airfoil_filenames = [bld_airfoil_filenames for i=1:Nbld]
            
        else
            blade_filenames = ["$path/blade$i.dat" for i=1:Nbld]
            blade_chords = [bldchord_spl for i=1:Nbld]
            blade_Nnodes = [NumADBldNds for i=1:Nbld]
            airfoil_filenames = collect(Iterators.flatten([bld_airfoil_filenames for i=1:Nbld]))
            
            for istrut = 1:Nstrutperbld
                strutchord_spl = OWENS.safeakima(numadIn_strut[istrut].span./maximum(numadIn_strut[istrut].span), numadIn_strut[istrut].chord,LinRange(0,1,NumADStrutNds))
                for ibld = 1:Nbld
                    blade_filenames = [blade_filenames;"$path/strut$(istrut)_bld$ibld.dat"]
                    blade_chords = [blade_chords;[strutchord_spl]]
                    blade_Nnodes = [blade_Nnodes;NumADStrutNds]

                    # Discretely assign the strut airfoils based on the next closest neighbor
                    strut_airfoil_filenames = fill("nothing",NumADStrutNds)
                    for (ispan_numad,span_numad) in enumerate(numadIn_strut[istrut].span./maximum(numadIn_strut[istrut].span))
                        for (ispan,span_slices) in enumerate(collect(LinRange(0,1,NumADBldNds)))
                            if strut_airfoil_filenames[ispan]=="nothing" && span_slices<=span_numad
                                strut_airfoil_filenames[ispan] = "$(numadIn_strut[istrut].airfoil[ispan_numad]).dat"
                            end
                        end
                    end

                    airfoil_filenames = [airfoil_filenames; strut_airfoil_filenames]

                end
            end
        end

        OWENSOpenFASTWrappers.writeADinputFile(ad_input_file,blade_filenames,airfoil_filenames,OLAF_filename;rho)

        NumADBody = length(AD15bldNdIdxRng[:,1])
        bld_len = zeros(NumADBody)
        # bladefileissaved = false
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
            ADshapeZ = collect(LinRange(0,H,NumADBldNds))
            xmesh = mymesh.x[strt_idx:end_idx]
            ymesh = mymesh.y[strt_idx:end_idx]
            ADshapeX = sqrt.(xmesh.^2 .+ ymesh.^2)
            ADshapeX .-= ADshapeX[1] #get it starting at zero #TODO: make robust for blades that don't start at 0
            ADshapeXspl = OWENS.safeakima(LinRange(0,H,length(ADshapeX)),ADshapeX,ADshapeZ)
            
            if iADBody<=Nbld #&& !bladefileissaved#Note that the blades can be curved and are assumed to be oriented vertically
                # bladefileissaved = true
                BlSpn0=ADshapeZ
                BlCrvAC0=ADshapeXspl

                bladeangle = (iADBody-1)*2.0*pi/Nbld + angularOffset #TODO: pitch offset and twist offset that isn't from the helical

                BlSpn = ADshapeZ
                blade_twist = atan.(xmesh,ymesh).-bladeangle

                #TODO: reevalueate these equations and make sure they are robust against varying designs
                # BlCrvACinput = xmesh.*cos(bladeangle) .+ ymesh.*sin(bladeangle)
                # BlCrvACinput = BlCrvACinput .- BlCrvACinput[1]
                # BlCrvAC = OWENS.safeakima(LinRange(0,H,length(BlCrvACinput)),BlCrvACinput,ADshapeZ)

                # BlSwpACinput = -ymesh.*cos(bladeangle) .+ xmesh.*sin(bladeangle)
                # BlSwpACinput = BlSwpACinput .- BlSwpACinput[1]
                # BlSwpAC = OWENS.safeakima(LinRange(0,H,length(BlSwpACinput)),BlSwpACinput,ADshapeZ)

                # BlCrvAng = zeros(blade_Nnodes[iADBody])

                # BlTwistinput =(blade_twist.-blade_twist[1])*180/pi
                # BlTwist = OWENS.safeakima(LinRange(0,H,length(BlTwistinput)),BlTwistinput,ADshapeZ)


                BlCrvACinput = -ymesh.*sin(bladeangle).+xmesh.*cos(bladeangle)
                BlCrvACinput = BlCrvACinput .- BlCrvACinput[1]
                BlSwpAC = -safeakima(LinRange(0,H,length(BlCrvACinput)),BlCrvACinput,ADshapeZ)

                BlSwpACinput = xmesh.*sin(bladeangle).+ymesh.*cos(bladeangle)
                BlSwpACinput = BlSwpACinput .- BlSwpACinput[1]
                BlCrvAC = safeakima(LinRange(0,H,length(BlSwpACinput)),BlSwpACinput,ADshapeZ)

                BlCrvAng = zeros(blade_Nnodes[iADBody])

                BlTwistinput =(blade_twist.-blade_twist[2])*180/pi
                BlTwist = safeakima(LinRange(0,H,length(BlTwistinput)),BlTwistinput,ADshapeZ)

                BlChord=blade_chords[iADBody]

                BlAFID=collect((iADBody-1)*NumADBldNds+1:iADBody*NumADBldNds)

            elseif iADBody>Nbld # while the arms/struts are assumed to be straight and are oriented by the mesh angle
                BlSpn=collect(LinRange(0,bld_len[iADBody],blade_Nnodes[iADBody]))
                BlCrvAC=zeros(blade_Nnodes[iADBody])
                BlSwpAC=zeros(blade_Nnodes[iADBody])
                BlCrvAng=zeros(blade_Nnodes[iADBody])
                BlTwist=zeros(blade_Nnodes[iADBody])
                BlChord=blade_chords[iADBody]
                BlAFID=collect((iADBody-1)*NumADStrutNds+1:iADBody*NumADStrutNds)
            end      
            OWENSOpenFASTWrappers.writeADbladeFile(filename;NumBlNds=blade_Nnodes[iADBody],BlSpn,BlCrvAC,BlSwpAC,BlCrvAng,BlTwist,BlChord,BlAFID)     
        end

        OWENSOpenFASTWrappers.writeOLAFfile(OLAF_filename;nNWPanel=200,nFWPanels=10)

        OWENSOpenFASTWrappers.writeIWfile(Vinf,ifw_input_file;WindType,windINPfilename)

        OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,[shapeX],[shapeZ],[B],[Htwr_base],[mymesh],[myort],[AD15bldNdIdxRng],[AD15bldElIdxRng];
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
                adi_nstrut=[Nstrutperbld],
                adi_debug=0,
                isHAWT = false     # true for HAWT, false for crossflow or VAWT
                )

        aeroForcesAD(t,azi) = OWENS.mapAD15(t,azi,[mymesh],OWENSOpenFASTWrappers.advanceAD15;alwaysrecalc=true,verbosity)
        deformAeroAD=OWENSOpenFASTWrappers.deformAD15
    end

    if return_componentized
        if AD15On
            return mymesh,myel,myort,myjoint,components,aeroForcesAD,deformAeroAD,system, assembly, sections
        else
            return mymesh,myel,myort,myjoint,components,aeroForcesACDMS,deformAeroACDMS,system, assembly, sections
        end
    else
        @warn "Not using the componetized model is being depreciated, please consider updating to return_componentized=true"

        mass_twr = []
        stiff_twr = []
        twr_precompinput = []
        twr_precompoutput = []
        plyprops_twr = []
        numadIn_twr = []
        lam_U_twr = []
        lam_L_twr = []
        mass_breakout_twr = []
        mass_bld = []
        stiff_bld = []
        bld_precompinput = []
        bld_precompoutput = []
        plyprops_bld = []
        numadIn_bld = []
        lam_U_bld = []
        lam_L_bld = []
        mass_breakout_blds = []

        for component in components
            if contains(component.name,"tower")
                mass_twr = component.mass_matrix
                stiff_twr = component.stiff_matrix
                twr_precompinput = component.preCompInput
                twr_precompoutput = component.preCompOutput
                plyprops_twr = component.plyProps
                numadIn_twr = component.nuMadIn
                lam_U_twr = component.lam_U
                lam_L_twr = component.lam_L
                mass_breakout_twr = component.mass
            end

            if contains(component.name,"blade1")
                mass_bld = component.mass_matrix
                stiff_bld = component.stiff_matrix
                bld_precompinput = component.preCompInput
                bld_precompoutput = component.preCompOutput
                plyprops_bld = component.plyProps
                numadIn_bld = component.nuMadIn
                lam_U_bld = component.lam_U
                lam_L_bld = component.lam_L
                mass_breakout_blds = component.mass
            end
        end
        
        if AD15On
            return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
            stiff_twr, stiff_bld,bld_precompinput,
            bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
            twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForcesAD,deformAeroAD,
            mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs,stiff_array,mass_array
        else
            return mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
            stiff_twr, stiff_bld,bld_precompinput,
            bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
            twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForcesACDMS,deformAeroACDMS,
            mass_breakout_blds,mass_breakout_twr,system,assembly,sections,AD15bldNdIdxRng, AD15bldElIdxRng, custom_mesh_outputs,stiff_array,mass_array
        end
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
