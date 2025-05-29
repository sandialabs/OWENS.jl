"""
    setupOWENS_config(path::String; mesh_config = default_mesh_config(), tower_config = default_tower_config(), blade_config = default_blade_config(), material_config = default_material_config(), aero_config = default_aero_config(), VTKmeshfilename = nothing, verbosity = 1, return_componentized = false)

Set up and configure an OWENS turbine mesh.

This function handles the setup process for an OWENS turbine model, including mesh generation,
sectional properties calculation, and aerodynamic model initialization. It processes various configuration
parameters and returns the assembled system components for analysis.

# Arguments
- `path::String`: Base path for file operations and data access

# Keyword Arguments
- `mesh_config`: Configuration for mesh generation (default: default_mesh_config())
- `tower_config`: Configuration for tower properties (default: default_tower_config())
- `blade_config`: Configuration for blade properties (default: default_blade_config())
- `material_config`: Configuration for material properties (default: default_material_config())
- `aero_config`: Configuration for aerodynamic properties (default: default_aero_config())
- `VTKmeshfilename`: Optional filename for VTK mesh output
- `verbosity`: Level of output verbosity (default: 1)
- `return_componentized`: Whether to return componentized model (default: false)

# Returns
When `return_componentized=true`:
- `mymesh`: Mesh properties
- `myel`: Element properties
- `myort`: Orientation properties
- `myjoint`: Joint properties
- `components`: Component definitions
- `aeroForces`: Aerodynamic forces
- `deformAero`: Aerodynamic deformation
- `system`: System properties
- `assembly`: Assembly properties
- `sections`: Section properties

When `return_componentized=false`:
Returns an extended set of properties including mass matrices, stiffness matrices, and various component-specific properties.

# Notes
- The function performs mesh setup, sectional properties calculation, and aerodynamic model initialization
- Supports both AD15 and ACDMS aerodynamic models
- Handles both componentized and non-componentized return formats
- Calculates total mass and material costs when verbosity > 1
"""
function setupOWENS_config(
    path::String;
    mesh_config = default_mesh_config(),
    tower_config = default_tower_config(),
    blade_config = default_blade_config(),
    material_config = default_material_config(),
    aero_config = default_aero_config(),
    VTKmeshfilename = nothing,
    verbosity = 1,
    return_componentized = false,
)
    custom_mesh_outputs = []

    # Unpack the mesh config
    meshtype = mesh_config.meshtype
    mesh_config.connectBldTips2Twr = meshtype == "Darrieus"


    # Unpack the blade config
    if minimum(blade_config.shapeZ)!=0
        @error "blade shapeZ must start at 0.0"
    end

    # Unpack the aero config
    aero_config.AD15On = aero_config.AeroModel == "AD"

    aero_config.centrifugal_force_flag =  material_config.AddedMass_Coeff_Ca>0.0

    # Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
    # We do some intermediate calculations on the blade shape and angles

    Nstrutperbld = length(tower_config.strut_twr_mountpoint)

    if typeof(tower_config.NuMad_geom_xlscsv_file_strut)==String || typeof(tower_config.NuMad_geom_xlscsv_file_strut)==OrderedCollections.OrderedDict{Symbol, Any}
        tower_config.NuMad_geom_xlscsv_file_strut = fill(tower_config.NuMad_geom_xlscsv_file_strut,Nstrutperbld)
    end

    nothing

    # Here we set up the mesh using one of the pre-made meshing functions.  For this case, there is a function for the ARCUS, as 
    # well as for towered VAWTs where you can have an arbitrary blade shape with connected struts, and if the blade tips touch
    # the tower, then you can tell it to constrain them to the tower thus allowing for both H-VAWT and Darrieus designs.

    #########################################
    ### Set up mesh
    #########################################
    mesh_props = setup_mesh(mesh_config, blade_config, tower_config, aero_config, verbosity)

    # Unpack the mesh properties
    mymesh = mesh_props.mymesh
    myort = mesh_props.myort
    myjoint = mesh_props.myjoint
    AD15bldNdIdxRng = mesh_props.AD15bldNdIdxRng
    AD15bldElIdxRng = mesh_props.AD15bldElIdxRng

    components,
    s2t_idx,
    intra_blade_cable_startidx_el,
    intra_blade_cable_endidx_el,
    topel_idx,
    s2b_idx = mesh_props.custom_mesh_outputs

    # nTwrElem = Int(mymesh.meshSeg[1])+1
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

    # This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

    #########################################
    ### Set up Sectional Properties
    #########################################

    sectional_props = setup_sectional_props(
        mesh_props,
        material_config,
        aero_config,
        tower_config,
        blade_config,
        path,
        # TODO: Add these
        nothing, # NuMad_geom_xlscsv_file_intra_blade_cable
        nothing, # NuMad_mat_xlscsv_file_intra_blade_cable
        nothing, # NuMad_geom_xlscsv_file_guys
        nothing, # NuMad_mat_xlscsv_file_guys
        verbosity
    )
    sectionPropsArray = sectional_props.sectionPropsArray
    stiff_array = sectional_props.stiff_array
    mass_array = sectional_props.mass_array
    rotationalEffects = sectional_props.rotationalEffects

    if verbosity>1
        println(
            "\nTotal Mass: $(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]])) kg",
        )
        println(
            "Total Material Cost: \$$(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]]))",
        )
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
        sectionPropsArray = [sectionPropsArray; fill(sectionPropsArray[end], n_diff)]
        stiff_array = [stiff_array; fill(stiff_array[end], n_diff)]
        mass_array = [mass_array; fill(mass_array[end], n_diff)]
    end

    #### store data in element object
    myel = OWENSFEA.El(
        sectionPropsArray,
        myort.Length,
        myort.Psi_d,
        myort.Theta_d,
        myort.Twist_d,
        rotationalEffects,
    )

    system, assembly, sections = OWENS.owens_to_gx(
        mymesh,
        myort,
        myjoint,
        sectionPropsArray,
        stiff_array,
        mass_array;
        VTKmeshfilename,
    )

    function calcMass(sectionPropsArray, myort)
        mass = 0.0
        lens = zeros(length(sectionPropsArray))
        rhoAs = zeros(length(sectionPropsArray))
        for (i, sectionProp) in enumerate(sectionPropsArray)
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
        spl_rhoa = FLOWMath.Akima(span_array, rhoAs)
        mass2, error = QuadGK.quadgk(spl_rhoa, span_array[1], span_array[end], atol = 1e-10)
        return mass, mass2
    end

    turb_masskg1, turb_masskg = calcMass(myel.props, myort)
    if verbosity>1
        println("\nMass of Turbine: $turb_masskg kg, $turb_masskg1 kg")
    end

    # Set up the OWENSAero aerodynamics if used
    
    aeroForcesAD, deformAeroAD, aeroForcesACDMS, deformAeroACDMS = setup_aerodynamic_model(
        blade_config,
        aero_config,
        tower_config,
        mesh_config,
        mesh_props,
        components,
        path,
        myel,
        verbosity,
    )

    # Initialize mass and stiffness variables
    # mass_twr = nothing
    # mass_bld = nothing
    # stiff_twr = nothing
    # stiff_bld = nothing
    # bld_precompinput = nothing
    # bld_precompoutput = nothing
    # plyprops_bld = nothing
    # lam_U_bld = nothing
    # lam_L_bld = nothing
    # twr_precompinput = nothing
    # twr_precompoutput = nothing
    # plyprops_twr = nothing
    # lam_U_twr = nothing
    # lam_L_twr = nothing
    # mass_breakout_blds = nothing
    # mass_breakout_twr = nothing

    # # Calculate mass breakouts
    # if !isnothing(mass_bld)
    #     mass_breakout_blds = sum(mass_bld)
    # end
    # if !isnothing(mass_twr)
    #     mass_breakout_twr = sum(mass_twr)
    # end

    # Return values based on componentized flag
    if return_componentized
        if aero_config.AD15On
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
        
        if aero_config.AD15On
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
