
"""
readNuMadGeomCSV(NuMad_geom_file)

Parameters defining the rotor (apply to all sections).

**Arguments**
- `NuMad_geom_file::String`: name of the numad excel CSV file being read (!!! THE NUMAD TAB MUST BE SAVED AS A CSV FOR THIS TO WORK !!!)
- `NuMad_geom_file::OrderedCollections.OrderedDict{Symbol, Any}`: Alternatively, the already loaded in dictionary of windio inputs


**Returns**
- `Output::NuMad`: numad structure as defined in the NuMad structure docstrings.
"""
function readNuMadGeomCSV(NuMad_geom_file::OrderedCollections.OrderedDict{Symbol, Any};section=:blade,span=nothing)

    # Reuse the input file as the dictionary input
    sec_Dict = NuMad_geom_file[:components][section]
    
    # internal_structure_2d_fem:
    # reference_axis:
    #     x:
    #         grid: [0.0, 0.03333333333333333, 0.06666666666666667, 0.1, 0.13333333333333333, 0.16666666666666666, 0.2, 0.23333333333333334, 0.26666666666666666, 0.3, 0.3333333333333333, 0.36666666666666664, 0.4, 0.43333333333333335, 0.4666666666666667, 0.5, 0.5333333333333333, 0.5666666666666667, 0.6, 0.6333333333333333, 0.6666666666666666, 0.7, 0.7333333333333333, 0.7666666666666667, 0.8, 0.8333333333333334, 0.8666666666666667, 0.9, 0.9333333333333333, 0.9666666666666667, 1.0]
    #         values: [0.0, 6.961447494399997, 13.442795161599998, 19.444043001599994, 24.965191014399995, 30.006239199999996, 34.5671875584, 38.6480360896, 42.248784793599995, 45.3694336704, 48.00998272, 50.1704319424, 51.8507813376, 53.0510309056, 53.7711806464, 54.01123056, 53.7711806464, 53.0510309056, 51.8507813376, 50.1704319424, 48.009982720000004, 45.36943367040001, 42.24878479360001, 38.6480360896, 34.56718755839999, 30.006239199999996, 24.965191014399995, 19.444043001599994, 13.442795161599998, 6.961447494399997, 0.0]
    #     y:
    #         grid: [0.0, 0.03333333333333333, 0.06666666666666667, 0.1, 0.13333333333333333, 0.16666666666666666, 0.2, 0.23333333333333334, 0.26666666666666666, 0.3, 0.3333333333333333, 0.36666666666666664, 0.4, 0.43333333333333335, 0.4666666666666667, 0.5, 0.5333333333333333, 0.5666666666666667, 0.6, 0.6333333333333333, 0.6666666666666666, 0.7, 0.7333333333333333, 0.7666666666666667, 0.8, 0.8333333333333334, 0.8666666666666667, 0.9, 0.9333333333333333, 0.9666666666666667, 1.0]
    #         values: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #     z:
    #         grid: [0.0, 0.03333333333333333, 0.06666666666666667, 0.1, 0.13333333333333333, 0.16666666666666666, 0.2, 0.23333333333333334, 0.26666666666666666, 0.3, 0.3333333333333333, 0.36666666666666664, 0.4, 0.43333333333333335, 0.4666666666666667, 0.5, 0.5333333333333333, 0.5666666666666667, 0.6, 0.6333333333333333, 0.6666666666666666, 0.7, 0.7333333333333333, 0.7666666666666667, 0.8, 0.8333333333333334, 0.8666666666666667, 0.9, 0.9333333333333333, 0.9666666666666667, 1.0]
    #         values: [0.0, 3.67276364, 7.34552728, 11.01829092, 14.69105456, 18.363818199999997, 22.03658184, 25.70934548, 29.38210912, 33.054872759999995, 36.727636399999994, 40.400400039999994, 44.07316368, 47.74592732, 51.41869096, 55.0914546, 58.76421824, 62.43698188, 66.10974551999999, 69.78250915999999, 73.45527279999999, 77.12803643999999, 80.80080007999999, 84.47356372, 88.14632736, 91.819091, 95.49185464, 99.16461828, 102.83738192, 106.51014556, 110.1829092]
    
    # Note that the input for the composite inputs is like a square blanket with shortened sections to make it any "shape" of composite inputs.  While it is possible to not define things along the blade/strut/tower etc depending on height, that makes non-square matrices which is more complex to code for and has not been propogated throughout.  It is also not condusive to continuous changes to the composite input.  Rather, just define the inputs as square (same number of chordwise stations along each spanwise position) and just set the thicknesses or distances between control points to be numerically 0.

    # So, the windio does not direcly have stack sequences, so if you want to do a sandwitch with foam inside, you have to define the layer material twice.  The ordered dictionary keeps everything in sequence, so we just assume the layers are defined from outer to inner.
    # Pull together all the keypoints, and then number them, and then read everything in and apply to square arrays, with zero thickness/layers where it isn't defined.
    # note that the keypoints are 0 to 1, trailing edge to trailing edge, with a leading edge position note.

    #TODO: unit span as much as is possible?

    airfoil_grid = sec_Dict[:outer_shape_bem][:airfoil_position][:grid]
    if isnothing(span)
        span = airfoil_grid
        println("Custom span is not specified in OWENS input, using WindIO airfoil grid as common span that all the other values are splined to")
    end

    if span[1]!=0.0 && span[end]!=1.0
        @error "Span definition must encompass the entire blade from root 0, to tip 1"
    end

    airfoil = Array{String,1}(undef,length(span))
    airfoil_names = sec_Dict[:outer_shape_bem][:airfoil_position][:labels]
    # spline the airfoils used to the current overall grid, then round to enable mapping to the discrete airfoil inputs
    airfoil_station_numbers = round.(Int,FLOWMath.akima(airfoil_grid,collect(1:length(airfoil_names)),span))
    # Now map the resulting airfoils to the new grid
    for istation = 1:length(airfoil_station_numbers)
        for iaf = 1:length(airfoil_names)
            if iaf == airfoil_station_numbers[istation]
                airfoil[istation] = airfoil_names[iaf]
                break
            end
        end
    end

    n_web = size(sec_Dict[:internal_structure_2d_fem][:webs])[1] # number of shear webs
    n_stack = size(sec_Dict[:internal_structure_2d_fem][:layers])[1] # number of stacks, note the to comply with the current windio, all layers are considered a stack

    te_type = nothing # this is unused and is not in the windio file.

    twist_grid = sec_Dict[:outer_shape_bem][:twist][:grid]
    twist_vals = sec_Dict[:outer_shape_bem][:twist][:values]
    twist_d = FLOWMath.akima(twist_grid,twist_vals,span) .* 180/pi

    chord_grid = sec_Dict[:outer_shape_bem][:chord][:grid]
    chord_vals = sec_Dict[:outer_shape_bem][:chord][:values]
    chord = FLOWMath.akima(chord_grid,chord_vals,span)

    pitch_axis_grid = sec_Dict[:outer_shape_bem][:pitch_axis][:grid]
    pitch_axis_vals = sec_Dict[:outer_shape_bem][:pitch_axis][:values]
    pitch_axis = FLOWMath.akima(pitch_axis_grid,pitch_axis_vals,span)

    xoffset = pitch_axis
    aerocenter = pitch_axis #TODO: this was originally used for the automated flutter analysis within the original OWENS code, which is not currently implemented

    DPtypes = nothing #currently unused

    # To get the number of segments, we need to figure out the chordwise segment control points, but this is tricky since each layer can have a different set...
    # So, let's get all of them and combine
    # Note that we assume that the leading edge is always an arc value of 0.5
    # We need to enforce common grids, otherwise the problem is impossible.  Pick out the first one, and if the others don't match, hard error
    # Then, with a common grid, we gather the starting and stopping values and pick out the unique for each layer
    # Then, we select the unique control points across all layers and apply to all span locations so it is a square matrix
    # Then, we use the starting and stopping locations to apply the zero thickness locations, and fill in the intermediate thickness
    # Then, the stack sequence for each position will be 1:N_layers

    # This is insanity, revert to just requiring each grid within a layer to being the same and only arc start and stop positions are accepted

    segments_bld = Array{Any,2}(undef,n_stack,2)
    notweb = zeros(Int,n_stack)
    for istack = 1:n_stack
        layer_Dict = sec_Dict[:internal_structure_2d_fem][:layers][istack]
        println(istack)
        println(layer_Dict[:name])
        if !(contains(layer_Dict[:name],"web"))
            notweb[istack] = 1
           
            if haskey(layer_Dict,:n_plies) && haskey(layer_Dict,:start_nd_arc) && haskey(layer_Dict,:end_nd_arc) && haskey(layer_Dict[:start_nd_arc],:grid) && haskey(layer_Dict[:end_nd_arc],:grid)
                material = layer_Dict[:material]

                n_plies_grid = layer_Dict[:n_plies][:grid]
                # n_plies_vals = layer_Dict[:n_plies][:values]

                fiber_orientation_grid = layer_Dict[:fiber_orientation][:grid]
                # fiber_orientation_vals = layer_Dict[:fiber_orientation][:values]

                start_nd_arc_grid = layer_Dict[:start_nd_arc][:grid]
                start_nd_arc_vals = layer_Dict[:start_nd_arc][:values]

                end_nd_arc_grid = layer_Dict[:end_nd_arc][:grid]
                end_nd_arc_vals = layer_Dict[:end_nd_arc][:values]
                    
                if (n_plies_grid != fiber_orientation_grid) && (n_plies_grid != start_nd_arc_grid) && (n_plies_grid != end_nd_arc_grid)
                    @error "specified grids within a layer must be the same"
                end
            else
                @error "For OWENS, to reduce the combinatorial number of input options and corner cases, please change the WindIO composite layer definitions to:
1) Only use start_nd_arc and end_nd_arc with grid and values, all other input options are not supported (mid and width, start and width, end and width, pitch axis and width, and the additional combinations of mid fixed, start fixed, end fixed (that is 12 different input combinations instead of 1).)
2) Use the same grid within a layer for each data type of the layer (otherwise, this requires error checking that the same starting and ending positions are used for each value as it is splined)
3) For webs, do likewise, and pitch axis with angle is not supported
4) If multiple methods are defined for a layer, only the start and end arc positions are used
5) Material n_plies instead of thickness is used"
            end

            segments_bld[istack,1] = start_nd_arc_vals
            segments_bld[istack,2] = end_nd_arc_vals
        end # if not web
    end # each layer

    segments_web = Array{Any,2}(undef,n_web,2)
    println("Webs")
    for iweb = 1:n_web # it is a web
        layer_Dict = sec_Dict[:internal_structure_2d_fem][:webs][iweb]
        # either we have a rotation and offset
        println(iweb)
        println(layer_Dict[:name])
 
        if haskey(layer_Dict,:start_nd_arc) && haskey(layer_Dict,:end_nd_arc) && haskey(layer_Dict[:start_nd_arc],:grid) && haskey(layer_Dict[:end_nd_arc],:grid)
            println("web start and end arc points")

            start_nd_arc_grid = layer_Dict[:start_nd_arc][:grid]
            start_nd_arc_vals = layer_Dict[:start_nd_arc][:values]
            start_nd_arc = FLOWMath.akima(start_nd_arc_grid,start_nd_arc_vals,span)

            end_nd_arc_grid = layer_Dict[:end_nd_arc][:grid]
            end_nd_arc_vals = layer_Dict[:end_nd_arc][:values]
            end_nd_arc = FLOWMath.akima(end_nd_arc_grid,end_nd_arc_vals,span)
        else 
            println("web rotation and offset")

            @error "Please specify web locations using the start and end arc points as defined in the windio and not the offset and rotation (this is not yet implemented)"

        end

        segments_web[iweb,1] = start_nd_arc
        segments_web[iweb,2] = end_nd_arc
    end


    # Combine to get the common keypoints since we are doing a square matrix for each span location and each chordwise location
    segments_temp2 = [segments_bld[notweb.==1,:];segments_web] 
    # determine unique keypoints
    start_nd_arc = segments_temp2[1,1]
    end_nd_arc = segments_temp2[1,2]
    common_segments_unsrt = unique([start_nd_arc;end_nd_arc])
    for istack = 2:length(segments_temp2[:,1])
        start_nd_arc = segments_temp2[istack,1]
        end_nd_arc = segments_temp2[istack,2]

        common_segments_unsrt = unique([common_segments_unsrt;start_nd_arc;end_nd_arc;[0.0,0.5,1.0]]) #add on ending points and leading edge to ensure they are there
    end
    common_segments = sort(common_segments_unsrt)
    n_segments = length(common_segments)

    segments_windio = zeros(length(span),n_segments)
    segments = zeros(length(span),n_segments) #note the reverse in definition between windio and numad!!!
    for ispan = 1:length(span)
        segments_windio[ispan,:] = common_segments
        segments[ispan,:] = reverse(-((common_segments .* 2.0) .- 1.0)) # convert from 0-0.5-1 trailing_lp-leading-trailing_hp, to Numad's 1-0-(-1), and then reverse, so that the high pressure side trailing edge is first
    end

    # Now that we have the common segments, we need an array of bit logic for each layer, matching the size of the segments, and we need to say if the stack is active for the given segment
    # Note that the High pressure trailing edge point is not included, but rather the low pressure is used, so all the matrices are one column shorter, but the segments are read in and the correction is applied internally

    stacks_active_windio_bld = zeros(Int,length(span),n_segments,n_stack)
    stacks_active_windio_web = zeros(Int,length(span),n_web,n_stack)
    
    # Also get the material types and layers used for each stack
    layer_mat_names = Array{String,1}(undef,n_stack)
    stack_mat_types = zeros(Int,n_stack) 
    N_materials = length(NuMad_geom_file[:materials])
    input_material_names = [NuMad_geom_file[:materials][imat][:name] for imat = 1:N_materials]
    stack_layers = zeros(length(span),n_stack) 

    for istack = 1:n_stack 
        layer_Dict = sec_Dict[:internal_structure_2d_fem][:layers][istack]
        
        println(istack)
        println(layer_Dict[:name])

        # Fill in the material types, which are in order based on the material names inputs
        stack_mat_types[istack] = findfirst(x->x==layer_Dict[:material],input_material_names)

        # get the number of plies
        n_plies_grid = layer_Dict[:n_plies][:grid]
        n_plies_vals = layer_Dict[:n_plies][:values]
        stack_layers[:,istack] = FLOWMath.akima(n_plies_grid,n_plies_vals,span) #note that since we cut off layers in the stack sequences below, extrapolated values are not used here

        if !(contains(layer_Dict[:name],"web")) 
            start_nd_arc_grid = layer_Dict[:start_nd_arc][:grid]
            start_nd_arc_vals = layer_Dict[:start_nd_arc][:values]
            start_nd_arc = FLOWMath.akima(start_nd_arc_grid,start_nd_arc_vals,span)

            end_nd_arc_grid = layer_Dict[:end_nd_arc][:grid]
            end_nd_arc_vals = layer_Dict[:end_nd_arc][:values]
            end_nd_arc = FLOWMath.akima(end_nd_arc_grid,end_nd_arc_vals,span)

            for ispan = 1:length(span)
                for iseg=1:length(common_segments)
                    # check that the layer is active for the span position, and that it is active for the chordwise position for the given span
                    if (span[ispan]>=start_nd_arc_grid[1] && span[ispan]<=end_nd_arc_grid[end]) && (common_segments[iseg]>=start_nd_arc[ispan] && common_segments[iseg]<=end_nd_arc[ispan])
                        stacks_active_windio_bld[ispan,iseg,istack] = 1
                    end
                end
            end
        else   
            for ispan = 1:length(span)
                for iweb=0:n_web-1 #this is per the windio standard...

                    # check that the layer is active for the span position, and that it is active for the chordwise position for the given span
                    if contains(layer_Dict[:web],"$iweb") && (span[ispan]>=n_plies_grid[1] && span[ispan]<=n_plies_grid[end])
                        stacks_active_windio_web[ispan,iweb+1,istack] = 1
                    end
                end
            end
        end

    end

    # flip the array to align with numad
    stacks_active_bld = reverse(stacks_active_windio_bld,dims=2)
    stacks_active_web = reverse(stacks_active_windio_web,dims=2)

    # Now take that logic and create the stack sequences
    skin_seq = Array{OWENS.Seq, 2}(undef, length(span),length(common_segments)) #can be any number of stack nums, so we have to make non-square containers

    for sta_idx = 1:length(span)
        # sta_idx = 1
        for seg_idx = 1:length(common_segments)
            # seg_idx = 1
            stack_array = []
            for istack = 1:n_stack       
                if stacks_active_bld[sta_idx,seg_idx,istack] == 1
                    stack_array = Int.([stack_array;istack])
                end
            end
            skin_seq[sta_idx,seg_idx] = OWENS.Seq(stack_array)
        end
    end
    
    # Now do the same but for the web
    web_seq = Array{OWENS.Seq, 2}(undef, length(span),n_web) #can be any number of stack nums, so we have to make non-square containers
    web_dp = Array{OWENS.Seq, 2}(undef, length(span),n_web) #this is fixed size square, but it's easier to do it this way

    for iweb = 1:n_web
        # iweb = 1
        start_nd_arc = segments_web[iweb,1]
        end_nd_arc = segments_web[iweb,2]
        for ispan = 1:length(span)
            # ispan = 1
            stack_array = []
            for istack = 1:n_stack       
                if stacks_active_web[ispan,iweb,istack] == 1
                    stack_array = Int.([stack_array;istack])
                end
            end
            if isempty(stack_array)
                @error "Please define shear webs for the entire span of the blade and set thickness to 0 where they are not included"
            end
            web_seq[ispan,iweb] = OWENS.Seq(stack_array)

            

            web_lp_idx_windio = findfirst(x->isapprox(start_nd_arc[ispan],x),segments_windio[ispan,:])
            web_hp_idx_windio = findfirst(x->isapprox(end_nd_arc[ispan],x),segments_windio[ispan,:])
            # now to from windio chordwise station to numad chordwise station by reversing the numbering
            # 1 2 3 4 5 6 7 8
            # 8 7 6 5 4 3 2 1
            web_lp_idx_temp = length(segments_windio[ispan,:])-web_lp_idx_windio + 1
            web_hp_idx_temp = length(segments_windio[ispan,:])-web_hp_idx_windio + 1

            # then since numad reuses the high pressure trailing edge position as the low pressure trailing edge (that is windio's last position), it is the first in numad, so offset everything by 1 to align
            web_lp_idx = web_lp_idx_temp
            web_hp_idx = web_hp_idx_temp

            web_dp[ispan,iweb] = OWENS.Seq([web_lp_idx;web_hp_idx;web_hp_idx;web_lp_idx])
        end
    end

    return NuMad(n_web,n_stack,n_segments-1,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
end

    # segments_bld = Array{Any,2}(undef,n_stack,2)
    # notweb = zeros(Int,n_stack)
    # for istack = 1:n_stack
    #     layer_Dict = sec_Dict[:internal_structure_2d_fem][:layers][istack]
    #     println(istack)
    #     println(layer_Dict[:name])
    #     if !(contains(layer_Dict[:name],"web"))
    #         notweb[istack] = 1
    #         tryfailed = false
    #         # midpoint and width
    #         try # midpoint is defined
    #             println("# midpoint is defined")
    #             try # midpoint has grid and vals
    #                 println("# midpoint has grid and vals")
    #                 midpoint_nd_arc_grid = layer_Dict[:midpoint_nd_arc][:grid]
    #                 midpoint_nd_arc_vals = layer_Dict[:midpoint_nd_arc][:values]
    #                 # midpoint_nd_arc = FLOWMath.akima(midpoint_nd_arc_grid,midpoint_nd_arc_vals,span)
    #             catch # midpoint is LE or TE
    #                 println("# midpoint is LE or TE")
    #                 if layer_Dict[:midpoint_nd_arc][:fixed] == "LE"
    #                     midpoint_nd_arc_vals = ones(length(span)).*0.5
    #                 elseif layer_Dict[:midpoint_nd_arc][:fixed] == "TE"
    #                     midpoint_nd_arc_vals = ones(length(span)).*0.0
    #                 else
    #                     @error "midpoint_nd_arc fixed must be LE or TE"
    #                 end
    #                 midpoint_nd_arc_grid = span
    #             end

    #             width_grid = layer_Dict[:width][:grid]
    #             width_vals = layer_Dict[:width][:values]
    #             width = FLOWMath.akima(width_grid,width_vals,midpoint_nd_arc_grid)./(chord*2) # unifying around start and end arc points requires width in meters (as defined by windio to be converted to the arc definition where)
    #             start_nd_arc = midpoint_nd_arc_vals - width./2
    #             end_nd_arc = midpoint_nd_arc_vals + width./2
    #             tryfailed = false
    #         catch
    #             tryfailed = true
    #         end

    #         if tryfailed
    #             # start and width, or start and end
    #             try # start is defined
    #                 println("# start is defined")

    #                 try #start had grid and vals
    #                     println("#start had grid and vals")
    #                     start_nd_arc_grid = layer_Dict[:start_nd_arc][:grid]
    #                     start_nd_arc_vals = layer_Dict[:start_nd_arc][:values]
    #                     # start_nd_arc = FLOWMath.akima(start_nd_arc_grid,start_nd_arc_vals,span)
    #                 catch # start is LE or TE
    #                     println("# start is LE or TE")
    #                     if layer_Dict[:start_nd_arc][:fixed] == "LE"
    #                         start_nd_arc_vals = ones(length(span)).*0.5
    #                     elseif layer_Dict[:start_nd_arc][:fixed] == "TE"
    #                         start_nd_arc_vals = ones(length(span)).*0.0
    #                     else
    #                         @error "start_nd_arc fixed must be LE or TE"
    #                     end
    #                     start_nd_arc_grid = span
    #                 end

    #                 try # if width is defined
    #                     println("# if width is defined")
    #                     width_grid = layer_Dict[:width][:grid]
    #                     if minimum(width_grid)>minimum(start_nd_arc_grid) || maximum(width_grid)<maximum(start_nd_arc_grid)
    #                         @error "width grid must be the same span as th start grid, and if the start_nd_arc is fixed, that is 0 to 1"
    #                     end
    #                     width_vals = layer_Dict[:width][:values]
    #                     width = FLOWMath.akima(width_grid,width_vals,start_nd_arc_grid)./(chord*2) # unifying around start and end arc points requires width in meters (as defined by windio) to be converted to the arc definition where 0.5 is leading edge, so span is approx 2x chord.
    #                     end_nd_arc_vals = start_nd_arc_vals+width

    #                 catch #otherwise, end must be specified
    #                     println("#otherwise, end must be specified")
    #                     try #end has grid and vals
    #                         # end_nd_arc_grid = layer_Dict[:end_nd_arc][:grid]
    #                         end_nd_arc_vals = layer_Dict[:end_nd_arc][:values]
    #                         # end_nd_arc = FLOWMath.akima(end_nd_arc_grid,end_nd_arc_vals,span)
    #                     catch # end is LE or TE
    #                         println("# end is LE or TE")
    #                         if layer_Dict[:end_nd_arc][:fixed] == "LE"
    #                             end_nd_arc_vals = ones(length(span)).*0.5
    #                         elseif layer_Dict[:end_nd_arc][:fixed] == "TE"
    #                             end_nd_arc_vals = ones(length(span)).*0.0
    #                         else
    #                             @error "end_nd_arc fixed must be LE or TE"
    #                         end
    #                     end
                        
    #                 end
    #                 tryfailed = false
    #             catch
    #                 tryfailed = true
    #             end
    #         end

    #         if tryfailed
    #             # end and width
    #             try
    #                 println("#end must be specified")
    #                 try #end has grid and vals
    #                     end_nd_arc_grid = layer_Dict[:end_nd_arc][:grid]
    #                     end_nd_arc_vals = layer_Dict[:end_nd_arc][:values]
    #                     # end_nd_arc = FLOWMath.akima(end_nd_arc_grid,end_nd_arc_vals,span)
    #                 catch # end is LE or TE
    #                     println("# end is LE or TE")
    #                     if layer_Dict[:end_nd_arc][:fixed] == "LE"
    #                         end_nd_arc = ones(length(span)).*0.5
    #                     elseif layer_Dict[:end_nd_arc][:fixed] == "TE"
    #                         end_nd_arc = ones(length(span)).*0.0
    #                     else
    #                         @error "end_nd_arc fixed must be LE or TE"
    #                     end
    #                     end_nd_arc_grid = span
    #                 end
    #                 println("# if width is defined")
    #                 width_grid = layer_Dict[:width][:grid]
    #                 if minimum(width_grid)>minimum(end_nd_arc_grid) || maximum(width_grid)<maximum(end_nd_arc_grid)
    #                     @error "width grid must be the same span as th end grid, and if the end_nd_arc is fixed, that is 0 to 1"
    #                 end
    #                 width_vals = layer_Dict[:width][:values]
    #                 width = FLOWMath.akima(width_grid,width_vals,end_nd_arc_grid)./(chord*2) # unifying around start and end arc points requires width in meters (as defined by windio to be converted to the arc definition where)
    #                 start_nd_arc_vals = end_nd_arc_vals-width

    #                 tryfailed = false
    #             catch
    #                 tryfailed = true
    #             end
    #         end

    #         if tryfailed
    #             # offset_y_pa and width
    #             println("offset_y_pa and width")
    #             offset_y_pa_grid = layer_Dict[:offset_y_pa][:grid]
    #             if minimum(offset_y_pa_grid)>minimum(span) || maximum(offset_y_pa_grid)<maximum(span)
    #                 @error "offset_y_pa_grid grid must span 0 to 1, set thickness to 0 for the layer for the span positions it is not present"
    #             end
    #             offset_y_pa_vals = layer_Dict[:offset_y_pa][:values]
    #             offset_y_pa = FLOWMath.akima(offset_y_pa_grid,offset_y_pa_vals,span) # in meters, need to convert to side

    #             side = layer_Dict[:side]

    #             if side == "suction"
    #                 midpoint_nd_arc = 0.5 .- (pitch_axis+offset_y_pa)./(chord*2) # 0.5 is the arc length where the reference axis is, minus the distance to the defined position normalized by the double arc length of the chord for the upper and lower halves
    #             else
    #                 midpoint_nd_arc = 0.5 .+ (pitch_axis+offset_y_pa)./(chord*2) # 0.5 is the arc length where the reference axis is, plus the distance to the defined position normalized by the double arc length of the chord for the upper and lower halves
    #             end

    #             width_grid = layer_Dict[:width][:grid]
    #             width_vals = layer_Dict[:width][:values]
    #             width = FLOWMath.akima(width_grid,width_vals,span)./(chord*2) # unifying around start and end arc points requires width in meters (as defined by windio to be converted to the arc definition where)
    #             start_nd_arc_vals = midpoint_nd_arc .- width./2
    #             end_nd_arc_vals = midpoint_nd_arc .+ width./2

    #             @warn "Please consider directly specifying the start and ending node arc positions for the highest accuracy"
    #             tryfailed = false
    #         end

    #         if tryfailed
    #             @error "Composite layer definitions do not comply with windio standard, please review"
    #         end

    #         if maximum(start_nd_arc_vals .< 0.0) || maximum(end_nd_arc_vals .> 1.0)
    #             @error "composite length definition exceeds arc length of 0 to 1 where 0 is the trailing suction side, 0.5 is leading edge, and 1.0 is trailing pressure side"
    #         end

    #         segments_bld[istack,1] = start_nd_arc_vals
    #         segments_bld[istack,2] = end_nd_arc_vals
    #     end # if not web
    # end # each layer

    # rotation_grid = layer_Dict[:rotation][:grid]
    # rotation_vals = layer_Dict[:rotation][:values]
    # rotation = FLOWMath.akima(rotation_grid,rotation_vals,span) # in meters, need to convert to side

    # offset_y_pa_grid = layer_Dict[:offset_y_pa][:grid]
    # offset_y_pa_vals = layer_Dict[:offset_y_pa][:values]
    # offset_y_pa = FLOWMath.akima(offset_y_pa_grid,offset_y_pa_vals,span) # in meters, need to convert to side

    # width_grid = layer_Dict[:width][:grid]
    # width_vals = layer_Dict[:width][:values]
    # width = FLOWMath.akima(width_grid,width_vals,span)./(chord*2) # unifying around start and end arc points requires width in meters (as defined by windio to be converted to the arc definition where)
    # start_nd_arc = midpoint_nd_arc - width./2
    # end_nd_arc = midpoint_nd_arc + width./2

    # this requires having the airfoil coordinates, which has not been implemented yet at this pitch_system_mass_cost_coeff


   

function readNuMadGeomCSV(NuMad_geom_file::String)
    #TODO: add composite orientation
    csvdata = DelimitedFiles.readdlm(NuMad_geom_file,',',skipstart = 0)

    n_station = length(csvdata[4:end,1])- sum(isempty.(csvdata[4:end,1]))
    n_web = Int(csvdata[1,6])
    n_stack = Int(csvdata[1,8])
    n_segments = Int(csvdata[2,8])
    span = Float64.(csvdata[4:n_station+3,2])
    #TODO: interpolations
    airfoil = csvdata[4:n_station+3,3]
    te_type = csvdata[4:n_station+3,4]
    twist_d = Float64.(csvdata[4:n_station+3,5])
    chord = Float64.(csvdata[4:n_station+3,6])
    xoffset = Float64.(csvdata[4:n_station+3,7])
    aerocenter = Float64.(csvdata[4:n_station+3,8])

    # Read stack info
    stack_idx_end = 10+n_stack-1
    stack_mat_types = Int.(csvdata[2,10:stack_idx_end])
    stack_layers = Float64.(csvdata[4:n_station+3,10:stack_idx_end])

    seg_idx_end = stack_idx_end+n_segments+1
    segments = Float64.(csvdata[4:n_station+3,stack_idx_end+1:seg_idx_end])

    DP_idx_end = seg_idx_end+n_segments+1
    DPtypes = csvdata[4:n_station+3,seg_idx_end+1:DP_idx_end]

    skin_seq = Array{Seq, 2}(undef, n_station,n_segments) #can be any number of stack nums, so we have to make non-square containers

    skin_idx_end = DP_idx_end+n_segments

    for sta_idx = 1:n_station
        # sta_idx = 1
        for seg_idx = 1:n_segments
            # seg_idx = 1
            str = split(csvdata[3+sta_idx,DP_idx_end+seg_idx],",")
            if length(str)>1 && str[2]=="" #allow for single number
                str = str[1]
            end
            skin_seq[sta_idx,seg_idx] = Seq([parse(Int,x) for x in str])
        end
    end

    web_seq = Array{Seq, 2}(undef, n_station,n_web) #can be any number of stack nums, so we have to make non-square containers
    web_dp = Array{Seq, 2}(undef, n_station,n_web) #this is fixed size square, but it's easier to do it this way

    for web_idx = 1:n_web
        # web_idx = 1
        for sta_idx = 1:n_station
            # sta_idx = 1
            str = split(csvdata[3+sta_idx,skin_idx_end+web_idx*2-1],",")
            if !isempty(str[1])
                if str[2]=="" #allow for single number
                    str = str[1]
                end
                web_seq[sta_idx,web_idx] = Seq([parse(Int,x) for x in str])

                str = split(csvdata[3+sta_idx,skin_idx_end+web_idx*2],",")
                web_dp[sta_idx,web_idx] = Seq([parse(Int,x) for x in str])
            end
        end
    end

    return NuMad(n_web,n_stack,n_segments,span,airfoil,te_type,twist_d,chord,xoffset,aerocenter,stack_mat_types,stack_layers,segments,DPtypes,skin_seq,web_seq,web_dp)
end

function readNuMadMaterialsCSV(NuMad_materials_xlscsv_file::OrderedCollections.OrderedDict{Symbol, Any})
    mat_Dict = NuMad_materials_xlscsv_file[:materials]

    N_materials = size(mat_Dict)[1]
    # Error checking on input file
    for  imat = 1:N_materials
        if !haskey(mat_Dict[imat],:name) || !haskey(mat_Dict[imat],:ply_t) || !haskey(mat_Dict[imat],:E) || !haskey(mat_Dict[imat],:G) || !haskey(mat_Dict[imat],:nu) || !haskey(mat_Dict[imat],:rho) || !haskey(mat_Dict[imat],:Xt) || !haskey(mat_Dict[imat],:Xc) || !haskey(mat_Dict[imat],:S)
            @error "Materials must specify name, ply_t, E, G, nu, rho, Xt, Xc, S,
Additionally, these are very useful: unit_cost, S, and m.  
Note that S is currently implemented as an array of floats corresponding to spline control points for stress in MPa, it should start near the ultimate strength and end near 0.  Then m is a corresponding array of floats that gives the Log10 of the cycles to failure, with a very large number corresponding to 0 stress (i.e. 20).  The number of control points for each material should be the same.
Also Note: the material names defined must match those used in the layers
And all materials are currently expecting orthotropic inputs, so please input multiples for steel etc."
        end
    end
    
    # TODO: propogate the unused orthogonal terms?
    names = [mat_Dict[imat][:name] for imat = 1:N_materials]
    plythickness = [mat_Dict[imat][:ply_t] for imat = 1:N_materials]
    e1 = [mat_Dict[imat][:E][1] for imat = 1:N_materials]
    e2 = [mat_Dict[imat][:E][2] for imat = 1:N_materials]
    g12 = [mat_Dict[imat][:G][1] for imat = 1:N_materials]
    anu = [mat_Dict[imat][:nu][1] for imat = 1:N_materials]
    rho = [mat_Dict[imat][:rho] for imat = 1:N_materials]
    xt = [mat_Dict[imat][:Xt][1] for imat = 1:N_materials]
    xc = [mat_Dict[imat][:Xc][1] for imat = 1:N_materials]
    yt = [mat_Dict[imat][:Xt][2] for imat = 1:N_materials]
    yc = [mat_Dict[imat][:Xc][2] for imat = 1:N_materials]
    s = [mat_Dict[imat][:S][1] for imat = 1:N_materials]
    costs = [mat_Dict[imat][:unit_cost] for imat = 1:N_materials]
    SN_stressMpa = collect(cat([mat_Dict[imat][:A] for imat = 1:N_materials][:,:]...,dims=2)')
    Log_SN_cycles2Fail = collect(cat([mat_Dict[imat][:m] for imat = 1:N_materials][:,:]...,dims=2)')

    return plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness),costs,SN_stressMpa,Log_SN_cycles2Fail)
end
"""
readNuMadMaterialsCSV(NuMad_materials_xlscsv_file)

Parameters defining the rotor materials.

**Arguments**
- `NuMad_materials_xlscsv_file::String`: name of the numad excel CSV file being read (!!! THE NUMAD TAB MUST BE SAVED AS A CSV FOR THIS TO WORK !!!)


**Returns**
- `Output::plyproperties`: plyproperties structure as defined in the plyproperties structure docstrings.
"""
function readNuMadMaterialsCSV(NuMad_materials_xlscsv_file)


    csvdata = DelimitedFiles.readdlm(NuMad_materials_xlscsv_file,',',skipstart = 0)

    data_start = 0
    data_end = length(csvdata[:,1])
    for ii = 1:length(csvdata[:,1])
        if csvdata[ii,1]=="Material ID"
            data_start = ii+1
        end

        if data_start!=0 && ii>data_start && !isa(csvdata[ii,1],Number) # in case they have a file with excess rows
            data_end = ii-1
            break
        end
    end

    names = csvdata[data_start:data_end,2]
    plythickness = Float64.(csvdata[data_start:data_end,4]) .* 1e-3 #meters
    e1 = Float64.(csvdata[data_start:data_end,5]) .* 1e6
    e2 = Float64.(csvdata[data_start:data_end,6]) .* 1e6
    g12 = Float64.(csvdata[data_start:data_end,8]) .* 1e6
    anu = Float64.(csvdata[data_start:data_end,11])  #ratio
    rho = Float64.(csvdata[data_start:data_end,14])  #g/cc * 1000 #kg/m3
    xt = abs.(Float64.(csvdata[data_start:data_end,15]) .* 1e6)  #pa
    xc = abs.(Float64.(csvdata[data_start:data_end,16]) .* 1e6)  #pa, abs since composites is looking for positive failure values and handles the negative.
    if length(csvdata[4,:])>17
        yt = abs.(Float64.(csvdata[data_start:data_end,17]) .* 1e6)  #pa
        yc = abs.(Float64.(csvdata[data_start:data_end,18]) .* 1e6)  #pa, abs since composites is looking for positive failure values and handles the negative.
        costs = Float64.(csvdata[data_start:data_end,21]) #$/kg
        try
            SN_stressMpa = Float64.(csvdata[data_start:data_end,23:28])
            Log_SN_cycles2Fail = Float64.(csvdata[data_start:data_end,29:34])
        catch
            SN_stressMpa = collect(cat(fill(collect(LinRange(1e12,0,6)),length(names))[:,:]...,dims=2)')
            Log_SN_cycles2Fail = collect(cat(fill(collect(LinRange(0,7,6)),length(names))[:,:]...,dims=2)')
            @warn "Data for SN curve control points not found in material file columns 23:28 for stress in Mpa, 29:33 for cycles in log10"
        end
    else
        yt = ones(length(e1)) .* 100.0e6 #made up
        yc = ones(length(e1)) .* 100.0e6  #made up, abs since composites is looking for positive failure values and handles the negative.
        costs = zeros(length(e1)) #$/kg
        SN_stressMpa = collect(cat(fill(collect(LinRange(1e12,0,6)),length(names))[:,:]...,dims=2)')
        Log_SN_cycles2Fail = collect(cat(fill(collect(LinRange(0,7,6)),length(names))[:,:]...,dims=2)')
    end

    s = abs.(ones(length(e1)) .* 100.0e6)  #made up, abs since composites.jl is expecting positive failure values

    return plyproperties(names,Composites.Material.(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness),costs,SN_stressMpa,Log_SN_cycles2Fail)

end


function saveOWENSfiles(filename,mymesh,myort,myjoint,myEl,pBC,numadIn_bld)
    # Create Filenames

    # Save mesh
    DelimitedFiles.open("$filename.mesh", "w") do io
           DelimitedFiles.writedlm(io, [mymesh.numNodes mymesh.numEl], ' ')
           DelimitedFiles.writedlm(io, [round.(Int,mymesh.nodeNum) mymesh.x mymesh.y mymesh.z], ' ')
           DelimitedFiles.writedlm(io, [round.(Int,mymesh.elNum) zeros(Int,length(mymesh.elNum)).+2 round.(Int,mymesh.conn)], ' ')
           DelimitedFiles.writedlm(io, " ", ' ')
           DelimitedFiles.writedlm(io, [length(mymesh.meshSeg) [meshSeg for meshSeg in mymesh.meshSeg]'], ' ')
       end

    # Save El

    DelimitedFiles.open("$filename.el", "w") do io
        for ii = 1:mymesh.numEl
            DelimitedFiles.writedlm(io, [mymesh.z[ii]/maximum(mymesh.z) -(myEl.props[ii].ac[1].+0.5)./2 myEl.props[ii].twist[1] myEl.props[ii].rhoA[1] myEl.props[ii].EIyy[1] myEl.props[ii].EIzz[1] myEl.props[ii].GJ[1] myEl.props[ii].EA[1] myEl.props[ii].alpha1[1] myEl.props[ii].rhoIyy[1] myEl.props[ii].rhoIzz[1] 0.0 0.0 myEl.props[ii].zcm[1] myEl.props[ii].ycm[1] 0.0 myEl.props[ii].a[1]], ' ')
            DelimitedFiles.writedlm(io, [mymesh.z[ii]/maximum(mymesh.z) -(myEl.props[ii].ac[2].+0.5)./2 myEl.props[ii].twist[2] myEl.props[ii].rhoA[2] myEl.props[ii].EIyy[2] myEl.props[ii].EIzz[2] myEl.props[ii].GJ[2] myEl.props[ii].EA[2] myEl.props[ii].alpha1[2] myEl.props[ii].rhoIyy[2] myEl.props[ii].rhoIzz[2] 0.0 0.0 myEl.props[ii].zcm[2] myEl.props[ii].ycm[2] 0.0 myEl.props[ii].a[2]], ' ')
        end
       end

    # Save Joint
    DelimitedFiles.open("$filename.jnt", "w") do io
           DelimitedFiles.writedlm(io, myjoint, '\t')
       end

    # Save Ort
    DelimitedFiles.open("$filename.ort", "w") do io
           DelimitedFiles.writedlm(io, [myort.elNum[:,1] myort.Psi_d myort.Theta_d myort.Twist_d myort.Length myort.Offset'], ' ')
       end

    # Save pBC
    DelimitedFiles.open("$filename.bc", "w") do io
           DelimitedFiles.writedlm(io, pBC, '\t')
       end

    # Save Blade
    # Used
    chord = FLOWMath.akima(numadIn_bld.span./maximum(numadIn_bld.span),numadIn_bld.chord,mymesh.structuralSpanLocNorm[1,:])

    bldArray = zeros(length(mymesh.structuralSpanLocNorm),16)
    row = 1
    for ibld = 1:length(mymesh.structuralSpanLocNorm[:,1])
        for ispan = 1:length(mymesh.structuralSpanLocNorm[1,:])
            bldArray[row,1] = ibld
            bldArray[row,2] = mymesh.structuralSpanLocNorm[ibld,ispan]
            bldArray[row,3] = mymesh.structuralNodeNumbers[ibld,ispan]
            bldArray[row,4] = mymesh.structuralElNumbers[ibld,ispan]
            bldArray[row,14] = chord[ispan]
            bldArray[row,16] = 2*pi
            row += 1
        end
    end
    DelimitedFiles.open("$filename.bld", "w") do io
           DelimitedFiles.writedlm(io, bldArray, '\t')
       end

    return nothing
end


"""

    readMesh(filename)

Reads the mesh file and stores data in the mesh object.

input:
* `filename::string` string containing mesh path/to/filename.mesh

output:
* `mesh::OWENSFEA.Mesh` see OWENSFEA.Mesh
"""
function readMesh(filename)

    fid = open(filename,"r")   #open mesh file

    # temp = fscanf(fid,'#i',2)   #read in number of nodes and number of elements
    line = readline(fid)
    temp = split(line)

    numNodes = parse(Int,temp[1])
    numEl = parse(Int,temp[2])

    nodeNum = zeros(numNodes,1)
    x = zeros(numNodes,1)
    y = zeros(numNodes,1)
    z = zeros(numNodes,1)

    conn = zeros(numEl,2)
    elNum = zeros(numEl,1)

    for i=1:numNodes            # read in node number and node coordinates
        line = readline(fid)
        temp = split(line)
        nodeNum[i] = parse(Float64,temp[1])
        x[i] = parse(Float64,temp[2])
        y[i] = parse(Float64,temp[3])
        z[i] = parse(Float64,temp[4])
    end

    for i=1:numEl               # read in element number and connectivity list
        line = readline(fid)
        temp = split(line)
        elNum[i] = parse(Float64,temp[1])

        conn[i,1] = parse(Float64,temp[3])
        conn[i,2] = parse(Float64,temp[4])
    end

    line = readline(fid) #get blank line
    line = readline(fid)
    temp = split(line)
    numComponents = parse(Int,temp[1])
    meshSeg = zeros(Int,numComponents)
    for i=1:numComponents
        meshSeg[i] = parse(Int,temp[i+1])
    end

    close(fid)  #close mesh file

    mesh = OWENSFEA.Mesh(nodeNum,
    numEl,
    numNodes,
    x,
    y,
    z,
    elNum,
    conn,
    zeros(Int,numEl),
    meshSeg,
    0,
    0,
    0)

    return mesh

end

"""

    readBCdata(bcfilename,numNodes,numDofPerNode)

This function reads the boundray condition file and stores data in the
boundary condition object.

#Input
* `bcfilename::string`:    string containing boundary condition filename
* `numNodes::int`:      number of nodes in structural model
* `numDofPerNode::int`: number of degrees of freedom per node

#Output
* `BC::OWENSFEA.BC_struct`:   see OWENSFEA.BC_struct, object containing boundary condition data
"""
function readBCdata(bcfilename,numNodes,numDofPerNode)

    fid = open(bcfilename)       #open boundary condition file
    numpBC = parse(Int,readline(fid)) #read in number of boundary conditions (displacement boundary conditions)
    pBC = zeros(Int,numpBC,3)         #initialize boundary conditions
    for i=1:numpBC

        line = readline(fid)

        # Find where all of the delimiters are
        #first two are boundary condition node number and local DOF number
        #third is boundary condition value (typically zero)
        delimiter_idx = [0;collect.(Int,findall(" ",line));length(line)+1]
        # Extract the data from the beginning to the last delimiter
        for k = 2:length(delimiter_idx)
            pBC[i,k-1] = Int(parse(Float64,line[delimiter_idx[k-1][1]+1:delimiter_idx[k][1]-1]))
        end

    end

    totalNumDof = numNodes*numDofPerNode

    numsBC = 0
    nummBC = 0

    close(fid)

    #create a vector denoting constrained DOFs in the model (0 unconstrained, 1
    #constrained)


    #calculate constrained dof vector
    isConstrained = zeros(totalNumDof,1)
    constDof = (pBC[:,1].-1)*numDofPerNode + pBC[:,2]
    index = 1
    for i=1:numNodes
        for j=1:numDofPerNode
            if ((i-1)*numDofPerNode + j in constDof)
                isConstrained[index] = 1
            end
            index = index + 1
        end
    end

    BC = OWENSFEA.BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    [],
    [])

    return BC

end

"""

    readBladeData(filename)

This function reads blade data from file

#Input
* `filename::string`:   string containing /path/to/bladedata.bld

#Output
* `bladeData::BladeData`:  see ?BladeData object containing blade data
"""
function readBladeData(filename)

    a = DelimitedFiles.readdlm(filename,'\t',skipstart = 0)

    bladeNum = a[:,1]

    numBlades = Int(maximum(bladeNum))
    #     numStruts = min(bladeNum)
    #     if (numStruts>0)
    #         numStruts = 0
    #     else
    #         numStruts = abs(numStruts)
    #     end

    strutStartIndex = 0
    for i=1:length(bladeNum)
        if (isempty(a[i,end]))
            strutStartIndex = i
            break
        end
    end



    if (strutStartIndex!=0)
        #         strutDataBlock = a(strutStartIndex:end,:)
        #         [strutEntries, _] = size(strutDataBlock)
        #         numNodesPerStrut = strutEntries/numStruts
        #         numElPerStrut = numNodesPerStrut - 1
    else
        temp=size(a)

        strutStartIndex = temp[1] + 1
    end

    bladeDataBlock = a[1:strutStartIndex-1,:]
    bladeEntries, _ = size(bladeDataBlock)
    numNodesPerBlade = round(Int,bladeEntries/numBlades)

    structuralSpanLocNorm = zeros(numBlades,numNodesPerBlade)
    structuralNodeNumbers = zeros(numBlades,numNodesPerBlade)
    structuralElNumbers = zeros(numBlades,numNodesPerBlade)
    for i=1:numBlades
        structuralSpanLocNorm[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,2]./bladeDataBlock[i*numNodesPerBlade,2]
        structuralNodeNumbers[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,3]
        structuralElNumbers[i,:] = bladeDataBlock[(i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,4]
    end

    bladeData = BladeData(numBlades,  #assign data to bladeData object
    bladeDataBlock[:,1],
    bladeDataBlock[:,2],
    bladeDataBlock[:,3],
    bladeDataBlock[:,4],
    bladeDataBlock[:,5:end])

    return bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers

end

"""

    readElementData(numElements,elfile,ortfile,bldfile

Reads element data and stores data in the element data
object.

#Input
* `numElements::int`:  number of elements in structural mesh
* `elfile::string`:       element data path/to/filename
* `ortfile::string`:      element orientation path/to/filename
* `bldfile::string`:      blade data path/to/filename

#Output
* `el::OWENSFEA.El`:       see OWENSFEA.El    element data object
"""
function readElementData(numElements,elfile,ortfile,bladeData_struct)

    fid = open(elfile,"r") #open element data file

    ac = zeros(2)
    twist = zeros(2)
    rhoA = zeros(2)
    EIyy = zeros(2)
    EIzz = zeros(2)
    GJ = zeros(2)
    EA = zeros(2)
    rhoIyy = zeros(2)
    rhoIzz = zeros(2)
    rhoJ = zeros(2)
    zcm = zeros(2)
    ycm = zeros(2)
    a = zeros(2)
    EIyz = zeros(2)
    alpha1 = zeros(2)
    alpha2 = zeros(2)
    alpha3 = zeros(2)
    alpha4 = zeros(2)
    alpha5 = zeros(2)
    alpha6 = zeros(2)
    rhoIyz = zeros(2)
    b = zeros(2)
    a0 = zeros(2)
    aeroCenterOffset = zeros(2)

    sectionPropsArray = Array{OWENSFEA.SectionPropsArray, 1}(undef, numElements)

    data1 = zeros(1,17)
    data2 = zeros(1,17)
    for i=1:numElements
        data1=parse.(Float64,split(readline(fid))) #read element data
        data2=parse.(Float64,split(readline(fid)))

        #structural properties
        ac = -([data1[2], data2[2]].-0.5) #TODO: why are we doing it this way???
        twist=[data1[3], data2[3]]
        rhoA = [data1[4], data2[4]]
        EIyy = [data1[5], data2[5]]
        EIzz = [data1[6], data2[6]]
        if (minimum(abs.(EIyy - EIzz)) < 1.0e-3)
            EIzz = EIzz.*1.0001
        end
        GJ = [data1[7], data2[7]]
        EA = [data1[8], data2[8]]
        alpha1 = [data1[9], data2[9]]

        rhoIyy = [data1[10], data2[10]]
        rhoIzz = [data1[11], data2[11]]
        rhoJ = [(data1[10]+data1[11]), (data2[10]+data2[11])]
        zcm = [data1[14], data2[14]]
        ycm = [data1[15], data2[15]]
        a = [data1[17], data2[17]]

        #coupling factors
        EIyz = [0.0, 0.0]
        alpha1 = [0.0, 0.0]
        alpha2 = [0.0, 0.0]
        alpha3 = [0.0, 0.0]
        alpha4 = [0.0, 0.0]
        alpha5 = [0.0, 0.0]
        alpha6 = [0.0, 0.0]
        rhoIyz = [0.0, 0.0]
        b = [0.0, 0.0]
        a0 = [2*pi, 2*pi]

        sectionPropsArray[i] = OWENSFEA.SectionPropsArray(ac,twist,rhoA,EIyy,EIzz,GJ,EA,rhoIyy,rhoIzz,rhoJ,zcm,ycm,a,EIyz,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,rhoIyz,b,a0,aeroCenterOffset)

    end
    close(fid) #close element file

    nodeNum = Int.(bladeData_struct.nodeNum)  #node number associated with blade section
    elNum = Int.(bladeData_struct.elementNum)    #element number associated with blade section
    bladeData = bladeData_struct.remaining  #blade data

    chord = zeros(maximum(Int,nodeNum),1)
    for i=1:length(elNum)
        chord[nodeNum[i]] = bladeData[i,10]  #store chord of blade sections
    end

    for i=1:length(elNum)
        if (elNum[i]!=-1)

            sectionPropsArray[elNum[i]].b = 0.5.*[chord[nodeNum[i]], chord[nodeNum[i+1]]] #element semi chord
            sectionPropsArray[elNum[i]].a0 = [bladeData[i,12], bladeData[i+1,12]]         #element lift curve slope (needed for flutter analysis)

            #convert "a" to semichord fraction aft of halfchord
            sectionPropsArray[elNum[i]].a = (sectionPropsArray[elNum[i]].a + 0.25*(2*sectionPropsArray[elNum[i]].b) - sectionPropsArray[elNum[i]].b)./sectionPropsArray[elNum[i]].b

            #convert "ac" to semichord fraction foreward of halfchord TODO: why are we doing it this way???
            sectionPropsArray[elNum[i]].ac = (sectionPropsArray[elNum[i]].ac).*2

            #physical aero center offset from elastic axis
            sectionPropsArray[elNum[i]].aeroCenterOffset = (sectionPropsArray[elNum[i]].ac).*sectionPropsArray[elNum[i]].b - sectionPropsArray[elNum[i]].a
        end
    end


    # println("EIyz, rhoIyz deactivated")

    #read element orientation data
    elLen = zeros(numElements)
    psi = zeros(numElements)
    theta = zeros(numElements)
    roll = zeros(numElements)
    fid = open(ortfile,"r")
    for i=1:numElements
        temp = parse.(Float64,split(readline(fid)))
        elLen[i]=temp[5]
        psi[i]=temp[2]
        theta[i]=temp[3]
        roll[i]=temp[4]
    end
    close(fid) #close ort file

    rotationalEffects = ones(numElements)

    #store data in element object
    el = OWENSFEA.El(sectionPropsArray,elLen,psi,theta,roll,rotationalEffects)

    return el

end

"""

    readGeneratorProps(generatorfilename)

This function reads generator properties from file.

#Input
* `generatorfilenanme::string`:  = string containing path/to/generatorfile

#Output
* `genprops`:    = model object containing generator properties
"""
function readGeneratorProps(generatorfilename)

    # fid = fopen(generatorfilename) #open generator property file
    # if (fid!=-1) #if file can be opened
    #         genprops.ratedTorque = fscanf(fid,'#f',1) #store rated torque
    #         dum = fgetl(fid)
    #         genprops.zeroTorqueGenSpeed = fscanf(fid,'#f',1) #store zero torque generator zpeed
    #         dum = fgetl(fid)
    #         genprops.pulloutRatio = fscanf(fid,'#f',1) #store pullout ratio
    #         dum = fgetl(fid)
    #         genprops.ratedGenSlipPerc= fscanf(fid,'#f',1) #store rated generator slip percentage
    #         dum = fgetl(fid)
    #
    #         fclose(fid) #close generator propery file
    genprops = 0.0
    println("Generator File to Generator properties not yet implemented")

    genprops = 0.0 #if generator property file does not exist, set object to null

    return genprops
end

"""

    writeOwensNDL(fileRoot, nodes, cmkType, cmkValues)

writes a nodal input file

#Intput
* `fileRoot::string`: string path to desired location with name but no extension
* `nodes::int`: node numbers for C/M/K
* `cmkType::string`: "C" "M" or "K"
* `cmkValues::float`: C/M/K value

#Output
* `none`:

"""
function writeOwensNDL(fileRoot, nodes, cmkType, cmkValues)

    # open the BC file to save the boundary conditions to
    BCfile = string(fileRoot, ".ndl")    #construct file name

    open(BCfile, "w") do file
        # write out the boundary conditions into the file
        for nn = 1:length(nodes)
            # [row, col, val] = find(cmkValues[nn])
            indices = findall(x->x!=0,cmkValues[:,:,nn])
            # println(indices)
            for ii = 1:length(indices)
                row = indices[ii][1]
                col = indices[ii][2]
                write(file, "$(nodes[nn]) $(cmkType[nn]) $(row) $(col) $(cmkValues[row,col,nn])\n")
            end
        end
    end
end


"""
Internal, reads modal file and returns freq, damp, and modeshapes
"""
function readResultsModalOut(resultsFile,numNodes)
    data = DelimitedFiles.readdlm(resultsFile,'\t',skipstart = 0)

    nmodes = round(Int,min(length(data[:,1])/(numNodes*2+6),30))

    freq = zeros(nmodes)
    damp = zeros(nmodes)
    U_x_0 = zeros(numNodes,nmodes)
    U_y_0 = zeros(numNodes,nmodes)
    U_z_0 = zeros(numNodes,nmodes)
    theta_x_0 = zeros(numNodes,nmodes)
    theta_y_0 = zeros(numNodes,nmodes)
    theta_z_0 = zeros(numNodes,nmodes)
    U_x_90 = zeros(numNodes,nmodes)
    U_y_90 = zeros(numNodes,nmodes)
    U_z_90 = zeros(numNodes,nmodes)
    theta_x_90 = zeros(numNodes,nmodes)
    theta_y_90 = zeros(numNodes,nmodes)
    theta_z_90 = zeros(numNodes,nmodes)

    for i_mode = 1:nmodes
        i_line = (i_mode-1)*(numNodes*2+7)+1

        freq[i_mode] = parse(Float64,(split(data[i_line+1,1])[2])[1:end-1])
        damp[i_mode] = parse(Float64,(split(data[i_line+2,1])[2])[1:end-1])

        # 0 degree shapes, with the max value scaled to 1

        temp = Float64.(data[i_line+5:i_line+4+numNodes,1])
        U_x_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,2])
        U_y_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,3])
        U_z_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,4])
        theta_x_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,5])
        theta_y_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,6])
        theta_z_0[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())

        i_line = i_line+numNodes+2 #90 degree shapes, with the max value scaled to 1

        temp = Float64.(data[i_line+5:i_line+4+numNodes,1])
        U_x_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,2])
        U_y_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,3])
        U_z_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,4])
        theta_x_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,5])
        theta_y_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
        temp = Float64.(data[i_line+5:i_line+4+numNodes,6])
        theta_z_90[:,i_mode] = temp#./max(maximum(abs.(temp)),eps())
    end
    return freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
end
