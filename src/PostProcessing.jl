"""
`my_getABD(matid::AbstractArray{<:Integer,1},
nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
theta::AbstractArray{<:Real,1},
q::AbstractArray{<:AbstractArray{<:Real,2},1}), offset::Real`

`my_getABD(lam::Laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1}), offset::Real`

Returns A, B, and D matrices

# Arguments:
* `matid::AbstractArray{<:Integer,1}`: material id of each lamina
* `nply::AbstractArray{<:Integer,1}`: number of plies in each lamina
* `tply::AbstractArray{<:Real,1}`: thickness of a ply (m) in each lamina
* `theta::AbstractArray{<:Real,1}`: orientation (deg) of each lamina
* `q::AbstractArray{<:AbstractArray{<:Real,2}}`: Stiffness matrix of each lamina
* `offset::Real`: Optional, used if neutral axis is not centered on the ply. If used, specify distance from neutral axis to upper surface (m).
"""
function my_getABD!(A::AbstractArray{<:Real,2}, B::AbstractArray{<:Real,2}, D::AbstractArray{<:Real,2},
    matid::AbstractArray{<:Integer,1}, nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, q::AbstractArray{<:AbstractArray{<:Real,2},1}, offset = 0.0)

    typeof(A)

    R = eltype(tply)

    nlam = length(nply)

    z = my_getz(tply, nply, offset)

    qbar = zeros(R, 3, 3)

    # Loop through layers filling in ABD matrix
    for k = 1:nlam
        # Rotate material stiffness properties by specifed angle
        imat = matid[k]
        qbar .= q[imat]
        Composites.rotQ!(qbar, theta[k])

        # Find lumped stiffness properties for the laminate
        for i = 1:3
            for j = 1:3
                A[i,j] = A[i,j]+qbar[i,j]*(z[k+1]-z[k])
                B[i,j] = B[i,j]+qbar[i,j]/2.0*(z[k+1]^2.0-z[k]^2.0)
                D[i,j] = D[i,j]+qbar[i,j]/3.0*(z[k+1]^3.0-z[k]^3.0)
            end
        end
    end

    return A, B, D
end

function my_getABD(matid::AbstractArray{<:Integer,1}, nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, q::AbstractArray{<:AbstractArray{<:Real,2},1}, offset = 0.0)

    R = eltype(tply)

    A = zeros(R, 3, 3)
    B = zeros(R, 3, 3)
    D = zeros(R, 3, 3)

    return my_getABD!(A, B, D, matid, nply, tply, theta, q, offset)
end

my_getABD(lam, q::AbstractArray{<:AbstractArray{<:Real,2},1}, offset = 0.0) =
my_getABD(lam.matid, lam.nply, lam.tply, lam.theta, q, offset)

"""
`getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1}, offset::Real)`

`getz(lam::Laminate, offset::Real)`

Returns a laminate's z-coordinates (coordinates of top and bottom of laminas)
given the thickness of plies in each lamina and the number of plies in each
lamina. Offset is optional, defaults to 0.0, is the distance from the neutal axis
to the upper ply if the laminate is not centered about the neutral axis of a beam.
"""
function my_getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1}, offset = 0.0)
    ttotal = 0.0
    tlam = zeros(eltype(tply), length(tply)+1)
    for i = 1:length(tply)
        ttotal += tply[i]*nply[i]
        tlam[i+1] = ttotal
    end
    tlam .-= ttotal/2.0
    if offset == 0.0
        return tlam
    else
        return tlam .+ (offset-ttotal/2.0)
    end
end
my_getz(lam, offset = 0.0) = my_getz(lam.tply, lam.nply, offset)


"""
`my_getplystrain(nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
theta::AbstractArray{<:Real,1}, resultantstrain::AbstractArray{<:Real,1}, offset::Real)`

`my_getplystrain(lam::Laminate, resultantstrain::AbstractArray{<:Real,1}, offset::Real)`

# Optional Argument
* `offset::Real`: Optional, used if neutral axis is not centered on the ply. If used, specify distance from neutral axis to upper surface (m).

Calculates strains in each ply aligned with principal material direction.
"""
function my_getplystrain(nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, resultantstrain::AbstractArray{<:Real,1}, offset = zeros(3))

    # number of laminas
    nlam = length(nply)

    # coordinates of top and bottom of plies
    z = my_getz(tply, nply, offset[3])

    # local strains at top and bottom of each lamina
    # localstrain = [resultantstrain[1:3]+z[k]*resultantstrain[4:6] for k=1:nlam+1]
    localstrain = zeros(3,nlam+1)
    # 4 is curvature about x axis, 5 is about y axis
    localstrain[1,:] = [resultantstrain[1]+z[k]*resultantstrain[5]+offset[2]*resultantstrain[6] for k=1:nlam+1]
    localstrain[2,:] = [resultantstrain[2]+z[k]*resultantstrain[4] for k=1:nlam+1]
    # localstrain[3,:] = [resultantstrain[3]+z[k]*resultantstrain[6] for k=1:nlam+1] #TODO: gamma_xy


    # rotate to align with material principal axes
    lowerplystrain = [Composites.rotstrain(localstrain[:,k], theta[k]) for k = 1:nlam]
    upperplystrain = [Composites.rotstrain(localstrain[:,k+1],theta[k]) for k = 1:nlam]

    return lowerplystrain, upperplystrain
end
my_getplystrain(lam, resultantstrain, offset = zeros(3)) =
my_getplystrain(lam.nply, lam.tply, lam.theta, resultantstrain, offset)

# Rainflow code originally from https://github.com/dhoegh/Rainflow.jl under MIT, but that code has been orphaned and the pull requests are not being evaluated.


""" This function sorts out points where the slope changes sign"""
function sort_peaks(signal::AbstractArray{Float64,1}, dt=collect(1.:length(signal)))
    slope = diff(signal)
    # Determines if the point is local extremum
    is_extremum = vcat(true, (slope[1:end-1].*slope[2:end]).<=0., true)
    return signal[is_extremum] , dt[is_extremum]
end

struct Cycle  # This is the information stored for each cycle found
    count::Float64
    range::Float64
    mean::Float64
    Rvalue::Float64   # value
    v_s::Float64   # value start
    t_s::Float64   # time start
    v_e::Float64   # value end
    t_e::Float64   # time end
end

show(io::IO,x::Cycle) = print(io, "Cycle: count=", x.count, ", range=",x.range, ", mean=",x.mean, ", R=", x.Rvalue)

function cycle(count::Float64, v_s::Float64, t_s::Float64, v_e::Float64, t_e::Float64)
    Cycle(count, abs(v_s-v_e), (v_s+v_e)/2, min(v_s,v_e)/max(v_s,v_e), v_s, t_s, v_e, t_e)
end

""" Count the cycles from the data """
function count_cycles(peaks::Array{Float64,1},t::Array{Float64,1})
    list = copy(peaks) # Makes a copy because they will be sorted in the vectors
    time = copy(t)
    currentindex = 1
    nextindex = 2
    cycles = Cycle[]
    sizehint!(cycles, length(list)) # This reduces the memory consumption a bit
    @inbounds begin
    while length(list) > (currentindex+1)
            currentvalue = abs(list[currentindex+1]-list[currentindex])
            nextvalue = abs(list[nextindex+1]-list[nextindex])
        if nextvalue > currentvalue
            if currentindex == 1 # This case counts a half and cycle deletes the poit that is counted
                push!(cycles,cycle(0.5 ,list[currentindex], time[currentindex], list[currentindex+1],time[currentindex+1]))
                popfirst!(list) # Removes the first entrance in ext and time
                popfirst!(time)
            else # This case counts one cycle and deletes the point that is counted
                push!(cycles,cycle(1. ,list[currentindex], time[currentindex], list[currentindex+1], time[currentindex+1]))
                deleteat!(list, currentindex:(currentindex+1)) # Removes the i and i+1 entrance in ext and time
                deleteat!(time, currentindex:(currentindex+1))
            end
            currentindex = 1
            nextindex = 2
        else
            currentindex += 1
            nextindex += 1
        end
    end
    for currentindex=1:length(list)-1 # This counts the rest of the points that have not been counted as a half cycle

        push!(cycles,cycle(0.5 ,list[currentindex], time[currentindex], list[currentindex+1],time[currentindex+1]))
    end
    end
    return cycles
end

mutable struct Cycles_bounds #
    min_mean::Float64
    max_mean::Float64
    max_range::Float64
    min_R::Float64
    max_R::Float64
end

show(io::IO,x::Cycles_bounds) = print(io, "Cycles_bounds : min mean value=", x.min_mean, ", max mean value=", x.max_mean, ", max range=",x.max_range, ", min R=", x.min_R, ", max R=",x.max_R)

""" Find the minimum and maximum mean value and maximum range from a vector of cycles"""
function find_boundary_vals(cycles::Array{Cycle,1})
    bounds = Cycles_bounds(Inf, -Inf, -Inf, Inf, -Inf)
    for cycle in cycles
        cycle.mean > bounds.max_mean && setfield!(bounds, :max_mean, cycle.mean)
        cycle.mean < bounds.min_mean && setfield!(bounds, :min_mean, cycle.mean)
        cycle.range > bounds.max_range && setfield!(bounds, :max_range, cycle.range)
        cycle.Rvalue < bounds.min_R && setfield!(bounds, :min_R, cycle.Rvalue)
        cycle.Rvalue > bounds.max_R && setfield!(bounds, :max_R, cycle.Rvalue)
    end
    return bounds
end

""" Returns the range index where the value is found """
function find_range(interval::Array{T,1},value) where {T <: Real}
    for i=1:length(interval)-1
        if interval[i] <= value < interval[i+1]
            return i
        elseif value == interval[i+1]
            return i
        elseif isapprox(interval[i],interval[i+1];atol=1e-7)
            return i
        end
    end
    println(interval[i])
    println(interval[i+1])
    println(value)
    error("The value is not in range")
end

Interval{T} = Union{Array{T,1}, StepRangeLen{T}}

""" Sums the cycle count given intervals of range_intervals and mean_intervals. The range_intervals and mean_intervals are given in fraction of range size"""
function sum_cycles(cycles::Array{Cycle,1}, range_intervals::Interval{T}, mean_intervals::Interval{T}) where {T <: Real}
    bounds = find_boundary_vals(cycles)
    bins = zeros(length(range_intervals)-1, length(mean_intervals)-1)
    range_in = (range_intervals*bounds.max_range)/100
    mean_in = (mean_intervals*(bounds.max_mean-bounds.min_mean))/100
    mean_in = mean_in .+ bounds.min_mean
    issorted(mean_intervals) || error("The array needs to be sorted in raising order")
    issorted(range_intervals) || error("The array needs to be sorted in raising order")
    nr_digits = 14  # The rounding is performed due to numerical noise in the floats when comparing
    mean_i = collect(mean_in)
    range_i = collect(range_in)
    #ensure the cycles are in the intervals by adding a small error to the end values of the interal.
    error_m = (bounds.max_mean-bounds.min_mean)*1e-14
    mean_i[end]+=error_m
    mean_i[1]-=error_m
    error_r = bounds.max_range*1e-14
    range_i[end]+=error_r
    range_i[1]-=error_r
    #show(mean_intervals)
    for cycle in cycles
        i = find_range(range_i, cycle.range)
        j = find_range(mean_i, cycle.mean)
        bins[i,j] += cycle.count
    end
    return bins, mean_i, range_i
end

function sum_cycles(cycles::Array{Cycle,1}, nr_ranges::Int=10, nr_means::Int=1)
    range_intervals = range(0,stop=100,length=nr_ranges+1)
    mean_intervals = range(0,stop=100,length=nr_means+1)
    sum_cycles(cycles, range_intervals, mean_intervals)
end

"""
    rainflow(signal;nbins_range=10,nbins_mean=10)

Convenience function that returns the binned cycles with the corresponding ranges and means

# Inputs
* `signal::Array{<:Real,1}`: data input
* `nbins_range::Array{<:Int,1}`: Number of bins for range
* `nbins_mean::Array{<:Int,1}`: Number of bins for mean
* `m`` :    Wohler exponent (default is 3)
* `Teq`` : The equivalent number of load cycles (default is 1, but normally the time duration in seconds is used)

# Outputs:
* `Ncycles::Array{<:Real,2}`: Summed/binned cycles with columns corresponding to mean levels and rows corresponding to range levels
* `meanIntervals::Array{<:Real,1}`: Mean levels corresponging with bins columns
* `rangeIntervals::Array{<:Real,1}`: Range levels corresponging with bins rows
* `equivalentLoad::Array{<:Real,1}`: Design equivalent load for each mean level

"""
function rainflow(signal;nbins_range=10,nbins_mean=10,m=3,Teq=1)
    extremum, t = sort_peaks(signal) # Sorts the signal to local extrema's, could optionally take a time vector
    cycles = count_cycles(extremum, t) # find all the cycles in the data
    range_intervals_percent = collect(LinRange(0,100,nbins_range+1)) # User defined intervals can also be specified
    mean_intervals_perc = collect(LinRange(0,100,nbins_mean+1))
    Ncycles,meanIntervals,rangeIntervals = sum_cycles(cycles, range_intervals_percent, mean_intervals_perc) # Sums the cycles in the given intervals

    # get DEL
    rangeIntervalLevels = (rangeIntervals[1:end-1] + rangeIntervals[2:end])./2 # as opposed to the edges
    equivalentLoad = zeros(nbins_mean)
    for imean = 1:nbins_mean
        DELs = rangeIntervalLevels.^m .* Ncycles[:,imean] ./ Teq
        equivalentLoad[imean] = sum(DELs) ^ (1/m)
    end

    return Ncycles,meanIntervals,rangeIntervals,equivalentLoad
end

##########################################
#### Composite Failure & Buckling ########
##########################################

function calcSF(total_t,stress,SF_ult,SF_buck,lencomposites_span,plyprops,
    precompinput,precompoutput,lam_in,eps_x,eps_z,eps_y,kappa_x,
    kappa_y,kappa_z,numadIn;failmethod = "maxstress",CLT=false,upper=true,layer=-1,calculate_fatigue=true)
    topstrainout = zeros(length(eps_x[:,1]),lencomposites_span,length(lam_in[1,:]),9) # time, span, lam, x,y, Assumes you use zero plies for sections that aren't used

    damage = zeros(lencomposites_span,length(lam_in[1,:])) #span length, with number of laminates, with number of plies NOTE: this assumes the number of plies is constant across all span and laminate locations
    for i_station = 1:lencomposites_span
        # i_station = 1
        
        if upper # Get laminate location starting and ending points
            lam_y_loc_le = precompinput[i_station].xsec_nodeU.*precompinput[i_station].chord
        else
            lam_y_loc_le = precompinput[i_station].xsec_nodeL.*precompinput[i_station].chord
        end
        lam_y_loc_le = (lam_y_loc_le[2:end]+lam_y_loc_le[1:end-1])./2 # Get the center locations

        # Get thicknesses
        refY = precompinput[i_station].le_loc*precompinput[i_station].chord # This is the distance from the leading edge to the reference axis.  Everything else is relative to the reference axis
        thickness_precomp_lag_te = (precompinput[i_station].chord-(refY+precompoutput[i_station].y_sc))
        thickness_precomp_lag_le = precompinput[i_station].chord - thickness_precomp_lag_te

        for j_lam = 1:length(lam_in[i_station,:])
            
            idx_y = round(Int,lam_y_loc_le[j_lam]/precompinput[i_station].chord*length(precompinput[i_station].ynode)/2)

            if upper
                thickness_precomp_flap = precompinput[i_station].ynode[idx_y]*precompinput[i_station].chord - precompoutput[i_station].x_sc
            else # Since the airfoil data wraps around (and is splined to always be the same number of points top and bottom), we index backwards 
                thickness_precomp_flap = precompinput[i_station].ynode[end-idx_y+1]*precompinput[i_station].chord - precompoutput[i_station].x_sc
            end
            offsetz = -thickness_precomp_flap #Negative due to sign convention of strain; i.e. if the blade starts at the bottom of the turbine, and goes out and bends up, the top should be in compression and the bottom in tension, which is opposite to the strain convention

            # Get y-distance from shear center to laminate location
            # Get location of j_lam relative to the shear center
            if j_lam == 1 # And address leading and trailing edge cases where we want the worst case and not center of lamina
                offsety = thickness_precomp_lag_le
            elseif j_lam == length(lam_in[i_station,:])
                offsety = -thickness_precomp_lag_te
            else
                offsety = -(lam_y_loc_le[j_lam]-(refY+precompoutput[i_station].y_sc))
            end

            offset = [0.0,offsety,offsetz]

            lam = lam_in[i_station,j_lam]

            stress_eachlayer = zeros(length(eps_x[:,1]),length(lam.nply),3) #all timesteps for each stack layer

            mat_idx = [numadIn.stack_mat_types[imat] for imat in lam.matid] # Map from the stack number to the material number
            materials = [plyprops.plies[imat] for imat in mat_idx] # Then get the actual material used
            SN_stressMpa = [plyprops.SN_stressMpa[imat,:] for imat in mat_idx]
            Log_SN_cycles2Fail = [plyprops.Log_SN_cycles2Fail[imat,:] for imat in mat_idx]
            # Map the stack number to the material number

            q = Composites.getQ.(materials,lam.theta)

            for its = 1:length(eps_x[:,1])
                # j_lam = 5
                # its = 1
                # println("its $its j_lam $j_lam i_station $i_station")
                resultantstrain = [
                eps_x[its,i_station],
                eps_z[its,i_station],
                eps_y[its,i_station],
                kappa_x[its,i_station],
                kappa_y[its,i_station],
                kappa_z[its,i_station]]


                lowerplystrain, upperplystrain = my_getplystrain(lam, resultantstrain, offset)
                topstrainout[its,i_station,j_lam,1:3] = upperplystrain[1]
                topstrainout[its,i_station,j_lam,4:end] = resultantstrain[:]
                upperplystress = q.*upperplystrain  #TODO: bottom side of plies needed for inter-laminate failure?
                out = Composites.getmatfail.(upperplystress,materials,failmethod)
                fail = [out[iii][1] for iii = 1:length(out)]
                sf = [out[iii][2] for iii = 1:length(out)]

                for ilayer = 1:length(lam.nply)
                    stress_eachlayer[its,ilayer,:] = upperplystress[ilayer]
                end

                if failmethod == "maxstress"
                    try
                        SF_ult[its,i_station,j_lam] = minimum(minimum([sf[isf][sf[isf].>0.0] for isf = 1:length(sf)])) # Pick out the worst case safety factor that is positive - negative means it is in the wrong direction for the failure criteria
                    catch
                    end
                else
                    @warn "Use maxstress for now, or inspect safety factors for each layer manually, need to identify which layer is failing"
                    if layer == -1
                        SF_ult[its,i_station,j_lam] = minimum(sf) #TODO; use findmin to identify which layer is failing
                    else
                        SF_ult[its,i_station,j_lam] = minimum(sf[layer])
                    end
                end
                if layer == -1
                    if length(upperplystress)>1
                        buck_layer = 1 #TODO: propogate throughout if a gel coat is being included
                    else
                        buck_layer = 1
                    end
                    stress[its,i_station,j_lam,:] = upperplystress[buck_layer]
                else
                    stress[its,i_station,j_lam,:] = upperplystress[layer]
                end

                # b is panel width
                # R is cylinder Radius
                # h is laminate total height

                # We are supplying the materials directly aligned with Q (i.e. repeated Q values aligned with each lamina) as opposed to by material
                A, B, D = my_getABD(collect(1:length(q)), lam.nply, lam.tply,lam.theta, q, offset[3])

                b = (numadIn.segments[i_station,j_lam+2]-numadIn.segments[i_station,j_lam+1])*numadIn.chord[i_station]
                R = numadIn.chord[i_station]/2
                laminate_h = sum(precompinput[i_station].t_lamU.*precompinput[i_station].n_pliesU)
                plate_uniaxial_local_buckling_stress = Composites.plate_uniaxial_local_buckling(A,D,b)
                plate_bending_local_buckling_stress = Composites.plate_bending_local_buckling(D,b)
                plate_shear_local_buckling_stress = Composites.plate_shear_local_buckling(A,D,b)

                weakest_buck = minimum([plate_uniaxial_local_buckling_stress,plate_bending_local_buckling_stress,plate_shear_local_buckling_stress])
                compression_load = min(stress[its,i_station,j_lam,1], stress[its,i_station,j_lam,2]) # if they are both negative, then it will give the worse compression load
                SF_buck[its,i_station,j_lam] = weakest_buck/(-compression_load)

                # cylinder_uniaxial_local_buckling_stress = Composites.cylinder_uniaxial_local_buckling(A, D, R, laminate_h)
                # cylinder_bending_local_buckling_stress = Composites.cylinder_bending_local_buckling(A, D, R, laminate_h)
            end

            # Miner's Damage
            if calculate_fatigue
                damage_layers = zeros(length(lam.nply))
                for ilayer = 1:length(lam.nply)
                    #TODO: non-principal stress
                    stressForFatigue = stress_eachlayer[:,ilayer,1]
                    ultimate_strength = materials[ilayer].xt
                    SN_stress = SN_stressMpa[ilayer] .* 1e6
                    Log_SN_cycles2Fail1 = Log_SN_cycles2Fail[ilayer]

                    # reverse order of cycles
                    if issorted(Log_SN_cycles2Fail1)
                        reverse!(SN_stress)
                        reverse!(Log_SN_cycles2Fail1)
                    end

                    damage_layers[ilayer] = fatigue_damage(stressForFatigue, SN_stress, Log_SN_cycles2Fail1, ultimate_strength)/total_t*60*60 #Damage to damage rate per hour
                end
                damage[i_station,j_lam] = maximum(damage_layers)
            end
        end
    end
    return topstrainout,damage
end

function fatigue_damage(stress, sn_stress, sn_log_cycles, ultimate_strength; nbins_amplitude=20, nbins_mean=nothing, mean_correction=true, wohler_exp=3, equiv_cycles=1)
    # default values
    isnothing(nbins_mean) && (nbins_mean = (mean_correction ? 10 : 1))

    ncycles, mean_bins, amplitude_bins, _ = rainflow(stress; nbins_range=nbins_amplitude, nbins_mean, m=wohler_exp, Teq=equiv_cycles)
    amplitude_levels = (amplitude_bins[1:end-1] .+ amplitude_bins[2:end]) ./ 2 # bin centers
    mean_levels = (mean_bins[1:end-1] .+ mean_bins[2:end]) ./ 2 # bin centers
    if mean_correction
        ncycles = [(ncycles...)] # flatten matrix
        effective_amplitude_levels = [iamplitude / (1 - imean / ultimate_strength) for iamplitude in amplitude_levels, imean in mean_levels]
        effective_amplitude_levels = [(effective_amplitude_levels...)] # flatten matrix
    else
        ncycles = sum(ncycles, dims=2) # sum over mean bins
        effective_amplitude_levels = amplitude_levels
    end
    log_ncycles_fail = safeakima(sn_stress, sn_log_cycles, effective_amplitude_levels) # interpolation
    ncycles_fail = 10.0 .^ log_ncycles_fail
    return sum(ncycles ./ ncycles_fail)
end

function printSF(verbosity,SF_ult,SF_buck,composite_station_idx, composite_station_name,lencomposites_span,lam_used,damage,delta_t,total_t;useStation=0)
    # verbosity: 0 Nothing, 1 summary, 2 summary and spar, 3 everything except bucking, 4 everything
    #Ultimate
    SF_ultmin,idx = findmin(SF_ult)
    damagemax, idxdamage = findmax(damage)
    idx1 = idx[1]
    idx2 = idx[2]
    idx3 = idx[3]
    idxdamage1 = idxdamage[1]
    idxdamage2 = idxdamage[2]

    if useStation != 0
        println("Using Specified composite station $useStation")
        SF_ultmin,idx = findmin(SF_ult[:,useStation,:])
        idx1 = idx[1]
        idx2 = useStation
        idx3 = idx[2]
        damagemax, idxdamage = findmax(damage[useStation,:])
        idxdamage1 = useStation
        idxdamage2 = idxdamage[1]
    end

    damageperhour = damage[idxdamage1,idxdamage2]/total_t*60*60

    if verbosity>0
        println("Minimum Safety Factor on Surface: $SF_ultmin")
        println("At time $(idx1*delta_t)s at composite station $(idx2) of $(lencomposites_span) at lam $(idx3) of $(length(lam_used[idx2,:]))")

        println("Maximum Damage per hr: $damageperhour")
        println("At composite station $(idxdamage1) of $(lencomposites_span) at lam $(idxdamage2) of $(length(lam_used[idxdamage1,:]))")
    end

    #Buckling
    if !isempty(SF_buck[SF_buck.>0.0])
        SF_buck[SF_buck.<0.0] .= 1e6
        # SF_buck[:,:,LE_idx] .= 1e6 #ignore leading edge
        # SF_buck[:,:,TE_idx] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck)
        if verbosity>2
            println("\nWorst buckling safety factor $(minbuck_sf)")
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(lencomposites_span) at lam $(minbuck_sfidx[3]) of $(length(lam_used[minbuck_sfidx[2],:]))")
        end
        if verbosity>3
            println("Buckling")
            for istation = 1:lencomposites_span
                println(minimum(SF_buck[minbuck_sfidx[1],istation,:]))
            end
        end
    else
        if verbosity>2
            println("Buckling not a factor, no sections in compression")
        end
    end

    function printlamInfo(SF_ultin,damagein,lam_idx,name,lencomposites_span,verbosity,delta_t,total_t)

        mymin,idx = findmin(SF_ultin[:,:,lam_idx])
        damagemax, idxdamage = findmax(damagein[:,lam_idx])
        damageperhour = damagemax/total_t*60*60

        idx1 = idx[1]
        idx2 = idx[2]
        idxdamage1 = idxdamage[1]

        if useStation != 0
            println("Using Specified composite station $useStation")
            idx2 = useStation
            idxdamage1 = useStation
        end

        if verbosity>1
            println("\n$name SF min: $mymin")
            println("At time $(idx1*delta_t)s at composite station $(idx2) of $(lencomposites_span)")

            println("$name Damage max per hour: $damagemax")
            println("At composite station $(idxdamage[1]) of $(lencomposites_span)")
        end

        if verbosity>2
            println("\n$name SF")
            for SF in SF_ultin[idx1,:,lam_idx]
                println(SF)
            end

            println("\n$name Damage")
            for DMGE in damagein[:,lam_idx]
                println(DMGE)
            end
        end
    end

    if verbosity>1
        for (i_iter,station_idx) in enumerate(composite_station_idx)
            printlamInfo(SF_ult,damage,station_idx,composite_station_name[i_iter],lencomposites_span,verbosity,delta_t,total_t)
        end
    end
end

function printsf_twr(verbosity,lam_twr,SF_ult_T,SF_buck_T,lencomposites_span,Twr_LE_idx,damage,delta_t,total_t)

    #Ultimate
    mymin,idx = findmin(SF_ult_T)
    damagemax, idxdamage = findmax(damage)
    damageperhour = damagemax/total_t*60*60
    if verbosity>0
        println("\nMinimum Safety Factor on tower Surface: $(minimum(SF_ult_T))")
        println("At time $(idx[1]*delta_t)s at composite station $(idx[2]) of $(lencomposites_span) at lam $(idx[3]) of $(length(lam_twr[idx[2],:]))")
        println("Maximum Damage per hr: $damageperhour")
        println("At composite station $(idxdamage[1]) of $(lencomposites_span) at lam $(idxdamage[2]) of $(length(lam_twr[idxdamage[1],:]))")


    end
    #Buckling
    if !isempty(SF_buck_T[SF_buck_T.>0.0])
        SF_buck_T[SF_buck_T.<0.0] .= 1e6
        # SF_buck_T[:,:,1] .= 1e6 #ignore leading edge
        # SF_buck_T[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_T)
        if verbosity>3
            println("\nWorst buckling safety factor $(minbuck_sf)")
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(lencomposites_span) at lam $(minbuck_sfidx[3]) of $(length(lam_twr[minbuck_sfidx[2],:]))")
        elseif verbosity>3
            println("Buckling")
            for istation = 1:lencomposites_span
                println(minimum(SF_buck_T[minbuck_sfidx[1],istation,:]))
            end
        end

    else
        if verbosity>3
            println("Buckling not a factor, no sections in compression")
        end
    end

    if verbosity>3
        println("\nLeading Edge SF")
        for SF in SF_ult_T[idx[1],:,Twr_LE_idx]
            println(SF)
        end

        println("\nLeading Edge Damage")
        for DMGE in damage[:,Twr_LE_idx]
            println(DMGE)
        end
    end
end

function extractSF(bld_precompinput,bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,Nbld,epsilon_x_hist_ps,kappa_y_hist_ps,kappa_z_hist_ps,epsilon_z_hist_ps,kappa_x_hist_ps,epsilon_y_hist_ps;
    epsilon_x_hist_1=nothing,kappa_y_hist_1=nothing,kappa_z_hist_1=nothing,
    epsilon_z_hist_1=nothing,kappa_x_hist_1=nothing,epsilon_y_hist_1=nothing,verbosity=2,usestationBld=0,usestationStrut=0,
    composite_station_idx_U_strut = [],#[1,6,3,2,5],
    composite_station_name_U_strut = [],#["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
    composite_station_idx_L_strut = [],#[1,6,3,2,5],
    composite_station_name_L_strut = [],#["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
    composite_station_idx_U_bld = [],#[1,6,3,2,5],
    composite_station_name_U_bld = [],#["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
    composite_station_idx_L_bld = [],#[1,6,3,2,5],
    composite_station_name_L_bld = [],#["Leading Edge","Trailing Edge","Spar Cap","Front Panel","Rear Panel"],
    Twr_LE_U_idx=1,Twr_LE_L_idx=1,throwawayTimeSteps=1,delta_t=0.001,AD15bldNdIdxRng=nothing,AD15bldElIdxRng=nothing,
    strut_precompoutput=nothing,strut_precompinput=nothing,plyprops_strut=nothing,numadIn_strut=nothing,lam_U_strut=nothing,lam_L_strut=nothing,calculate_fatigue=true)

    @warn "extractSF is being depreciated in favor of safetyfactor_fatigue, which uses the componetized workflow"

    # Linearly Superimpose the Strains
    epsilon_x_hist = epsilon_x_hist_ps#copy(epsilon_x_hist_ps).*0.0
    kappa_y_hist = kappa_y_hist_ps#copy(kappa_y_hist_ps).*0.0
    kappa_z_hist = kappa_z_hist_ps#copy(kappa_z_hist_ps).*0.0
    epsilon_z_hist = epsilon_z_hist_ps#copy(epsilon_z_hist_ps).*0.0
    kappa_x_hist = kappa_x_hist_ps#copy(kappa_x_hist_ps).*0.0
    epsilon_y_hist = epsilon_y_hist_ps#copy(epsilon_y_hist_ps).*0.0

    if !isnothing(epsilon_x_hist_1)
        for ipt = 1:4
            for jts = throwawayTimeSteps:length(epsilon_x_hist_ps[1,1,:])
                epsilon_x_hist[ipt,:,jts] = epsilon_x_hist_ps[ipt,:,jts] + epsilon_x_hist_1[ipt,:,end]
                kappa_y_hist[ipt,:,jts] = kappa_y_hist_ps[ipt,:,jts] + kappa_y_hist_1[ipt,:,end]
                kappa_z_hist[ipt,:,jts] = kappa_z_hist_ps[ipt,:,jts] + kappa_z_hist_1[ipt,:,end]
                epsilon_z_hist[ipt,:,jts] = epsilon_z_hist_ps[ipt,:,jts] + epsilon_z_hist_1[ipt,:,end]
                kappa_x_hist[ipt,:,jts] = kappa_x_hist_ps[ipt,:,jts] + kappa_x_hist_1[ipt,:,end]
                epsilon_y_hist[ipt,:,jts] = epsilon_y_hist_ps[ipt,:,jts] + epsilon_y_hist_1[ipt,:,end]
            end
        end
    else
        epsilon_x_hist[:,:,throwawayTimeSteps:end] = epsilon_x_hist_ps[:,:,throwawayTimeSteps:end]
        kappa_y_hist[:,:,throwawayTimeSteps:end] = kappa_y_hist_ps[:,:,throwawayTimeSteps:end]
        kappa_z_hist[:,:,throwawayTimeSteps:end] = kappa_z_hist_ps[:,:,throwawayTimeSteps:end]
        epsilon_z_hist[:,:,throwawayTimeSteps:end] = epsilon_z_hist_ps[:,:,throwawayTimeSteps:end]
        kappa_x_hist[:,:,throwawayTimeSteps:end] = kappa_x_hist_ps[:,:,throwawayTimeSteps:end]
        epsilon_y_hist[:,:,throwawayTimeSteps:end] = epsilon_y_hist_ps[:,:,throwawayTimeSteps:end]

        epsilon_x_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
        kappa_y_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
        kappa_z_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
        epsilon_z_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
        kappa_x_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
        epsilon_y_hist[:,:,1:throwawayTimeSteps-1] .= 0.0
    end

    total_t = length(epsilon_x_hist[1,1,throwawayTimeSteps:end])*delta_t # for damage rate
    #############################################
    #### Get strain values at the blades ########
    #############################################

    meanepsilon_z_hist = Statistics.mean(epsilon_z_hist,dims=1)
    meanepsilon_y_hist = Statistics.mean(epsilon_y_hist,dims=1)

    N_ts = length(epsilon_x_hist[1,1,:])
    eps_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    eps_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    eps_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    for ibld = 1:Nbld
        startN = Int(AD15bldNdIdxRng[ibld,1])
        stopN = Int(AD15bldNdIdxRng[ibld,2])
        if stopN < startN
            temp = stopN
            stopN = startN
            startN = temp
        end

        startE = Int(AD15bldElIdxRng[ibld,1])
        stopE = Int(AD15bldElIdxRng[ibld,2])
        if stopE < startE
            temp = stopE
            stopE = startE
            startE = temp
        end

        mesh_span_bld = mymesh.z[startN:stopN].-mymesh.z[startN]
        mesh_span_bld_el = safeakima(LinRange(0,1,length(mesh_span_bld)),mesh_span_bld,LinRange(0,1,stopE-startE+1))
        composites_span_bld = numadIn_bld.span/maximum(numadIn_bld.span)*maximum(mesh_span_bld_el) #TODO: this assumes struts and blades are straight, should be mapped for all potential cases, also need to input hub offset

        for its = 1:N_ts
            # Interpolate to the composite inputs #TODO: verify node vs el in strain
            eps_x_bld[ibld,its,:] = safeakima(mesh_span_bld_el,epsilon_x_hist[1,startE:stopE,its],composites_span_bld)
            eps_z_bld[ibld,its,:] = safeakima(mesh_span_bld_el,meanepsilon_z_hist[1,startE:stopE,its],composites_span_bld)
            eps_y_bld[ibld,its,:] = safeakima(mesh_span_bld_el,meanepsilon_y_hist[1,startE:stopE,its],composites_span_bld)
            kappa_x_bld[ibld,its,:] = safeakima(mesh_span_bld_el,kappa_x_hist[1,startE:stopE,its],composites_span_bld)
            kappa_y_bld[ibld,its,:] = safeakima(mesh_span_bld_el,kappa_y_hist[1,startE:stopE,its],composites_span_bld)
            kappa_z_bld[ibld,its,:] = safeakima(mesh_span_bld_el,kappa_z_hist[1,startE:stopE,its],composites_span_bld)
        end
    end

    if !isnothing(strut_precompoutput)
        ####################################################
        #### Get strain values at the struts ########
        ####################################################

        eps_x_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))
        eps_z_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))
        eps_y_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))
        kappa_x_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))
        kappa_y_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))
        kappa_z_strut = zeros(Nbld*2,N_ts,length(strut_precompinput))

        istrut = 0
        for ibld = Nbld+1:length(AD15bldNdIdxRng[:,1])
            istrut += 1

            startN = Int(AD15bldNdIdxRng[ibld,1])
            stopN = Int(AD15bldNdIdxRng[ibld,2])
            if stopN < startN
                temp = stopN
                stopN = startN
                startN = temp
            end

            startE = Int(AD15bldElIdxRng[ibld,1])
            stopE = Int(AD15bldElIdxRng[ibld,2])
            if stopE < startE
                temp = stopE
                stopE = startE
                startE = temp
            end

            mesh_span_strut = abs.(mymesh.z[startN:stopN].-mymesh.z[startN])
            mesh_span_strut_el = safeakima(LinRange(0,1,length(mesh_span_strut)),mesh_span_strut,LinRange(0,1,stopE-startE+1))
            composites_span_strut = numadIn_strut.span/maximum(numadIn_strut.span)*maximum(mesh_span_strut_el) #TODO: this assumes struts and blades are straight, should be mapped for all potential cases, also need to input hub offset

            for its = 1:N_ts
                # Interpolate to the composite inputs #TODO: verify node vs el in strain
                eps_x_strut[istrut,its,:] = safeakima(mesh_span_strut_el,epsilon_x_hist[1,startE:stopE,its],composites_span_strut)
                eps_z_strut[istrut,its,:] = safeakima(mesh_span_strut_el,meanepsilon_z_hist[1,startE:stopE,its],composites_span_strut)
                eps_y_strut[istrut,its,:] = safeakima(mesh_span_strut_el,meanepsilon_y_hist[1,startE:stopE,its],composites_span_strut)
                kappa_x_strut[istrut,its,:] = safeakima(mesh_span_strut_el,kappa_x_hist[1,startE:stopE,its],composites_span_strut)
                kappa_y_strut[istrut,its,:] = safeakima(mesh_span_strut_el,kappa_y_hist[1,startE:stopE,its],composites_span_strut)
                kappa_z_strut[istrut,its,:] = safeakima(mesh_span_strut_el,kappa_z_hist[1,startE:stopE,its],composites_span_strut)
            end
        end
    end

    ##########################################
    #### Get strain values at the tower #####
    ##########################################

    N_ts = length(epsilon_x_hist[1,1,:])
    eps_x_twr = zeros(1,N_ts,length(twr_precompinput))
    eps_z_twr = zeros(1,N_ts,length(twr_precompinput))
    eps_y_twr = zeros(1,N_ts,length(twr_precompinput))
    kappa_x_twr = zeros(1,N_ts,length(twr_precompinput))
    kappa_y_twr = zeros(1,N_ts,length(twr_precompinput))
    kappa_z_twr = zeros(1,N_ts,length(twr_precompinput))

    start = 1
    stop = length(mymesh.type[mymesh.type.==1])
    x = mymesh.z[start:stop]
    x = x.-x[1] #zero
    x = x./x[end] #normalize
    mesh_span_twr = mymesh.z[start:stop].-mymesh.z[start]
    composites_span_twr = safeakima(LinRange(0,1,length(mesh_span_twr)),mesh_span_twr,LinRange(0,1,length(twr_precompinput)))
    for its = 1:N_ts
        # Interpolate to the composite inputs
        eps_x_twr[1,its,:] = safeakima(mesh_span_twr,epsilon_x_hist[1,start:stop,its],composites_span_twr)
        eps_z_twr[1,its,:] = safeakima(mesh_span_twr,meanepsilon_z_hist[1,start:stop,its],composites_span_twr)
        eps_y_twr[1,its,:] = safeakima(mesh_span_twr,meanepsilon_y_hist[1,start:stop,its],composites_span_twr)
        kappa_x_twr[1,its,:] = safeakima(mesh_span_twr,kappa_x_hist[1,start:stop,its],composites_span_twr)
        kappa_y_twr[1,its,:] = safeakima(mesh_span_twr,kappa_y_hist[1,start:stop,its],composites_span_twr)
        kappa_z_twr[1,its,:] = safeakima(mesh_span_twr,kappa_z_hist[1,start:stop,its],composites_span_twr)

    end

    ##########################################
    #### Calculate Stress At the Blades
    ##########################################

    stress_U = zeros(N_ts,length(bld_precompinput),length(lam_U_bld[1,:]),3)
    SF_ult_U = zeros(N_ts,length(bld_precompinput),length(lam_U_bld[1,:]))
    SF_buck_U = zeros(N_ts,length(bld_precompinput),length(lam_U_bld[1,:]))

    topstrainout_blade_U,topDamage_blade_U = calcSF(total_t,stress_U,SF_ult_U,SF_buck_U,length(bld_precompinput),plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_U_bld,eps_x_bld[1,:,:],eps_z_bld[1,:,:],eps_y_bld[1,:,:],kappa_x_bld[1,:,:],
    kappa_y_bld[1,:,:],kappa_z_bld[1,:,:],numadIn_bld;failmethod = "maxstress",upper=true,calculate_fatigue)

    if verbosity>0
        println("Composite Ultimate and Buckling Safety Factors")
        println("\n\nUPPER BLADE SURFACE")
    end
    if usestationBld !=0
        println("maximum blade stress: $(maximum(stress_U[:,usestationBld,8,1]))")
        println("minimum blade stress: $(minimum(stress_U[:,usestationBld,8,1]))")
    end

    printSF(verbosity,SF_ult_U,SF_buck_U,composite_station_idx_U_bld, composite_station_name_U_bld,length(bld_precompinput),lam_U_bld,topDamage_blade_U,delta_t,total_t;useStation=usestationBld)

    stress_L = zeros(N_ts,length(bld_precompinput),length(lam_U_bld[1,:]),3)
    SF_ult_L = zeros(N_ts,length(bld_precompinput),length(lam_L_bld[1,:]))
    SF_buck_L = zeros(N_ts,length(bld_precompinput),length(lam_L_bld[1,:]))

    topstrainout_blade_L,topDamage_blade_L = calcSF(total_t,stress_L,SF_ult_L,SF_buck_L,length(bld_precompinput),plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_L_bld,eps_x_bld[1,:,:],eps_z_bld[1,:,:],eps_y_bld[1,:,:],kappa_x_bld[1,:,:],
    kappa_y_bld[1,:,:],kappa_z_bld[1,:,:],numadIn_bld;failmethod = "maxstress",upper=false,calculate_fatigue)

    if verbosity>0
        println("\n\nLOWER BLADE SURFACE")
    end
    if usestationBld !=0
        println("maximum blade stress: $(maximum(stress_L[:,usestationBld,8,1]))")
        println("minimum blade stress: $(minimum(stress_L[:,usestationBld,8,1]))")
    end
    printSF(verbosity,SF_ult_L,SF_buck_L,composite_station_idx_L_bld, composite_station_name_L_bld,length(bld_precompinput),lam_L_bld,topDamage_blade_L,delta_t,total_t;useStation=usestationBld)

    if !isnothing(strut_precompoutput)
        ##########################################
        #### Calculate Stress At the Struts
        ##########################################

        stress_U_strut = zeros(N_ts,length(strut_precompinput),length(lam_U_strut[1,:]),3)
        SF_ult_U_strut = zeros(N_ts,length(strut_precompinput),length(lam_U_strut[1,:]))
        SF_buck_U_strut = zeros(N_ts,length(strut_precompinput),length(lam_U_strut[1,:]))

        topstrainout_strut_U,topDamage_strut_U = calcSF(total_t,stress_U_strut,SF_ult_U_strut,SF_buck_U_strut,length(strut_precompinput),plyprops_strut,
        strut_precompinput,strut_precompoutput,lam_U_strut,eps_x_strut[1,:,:],eps_z_strut[1,:,:],eps_y_strut[1,:,:],kappa_x_strut[1,:,:],
        kappa_y_strut[1,:,:],kappa_z_strut[1,:,:],numadIn_strut;failmethod = "maxstress",upper=true,calculate_fatigue)

        if verbosity>0
            println("Composite Ultimate and Buckling Safety Factors")
            println("\n\nUPPER STRUT SURFACE")
        end

        if usestationStrut !=0
            println("maximum strut stress: $(maximum(stress_U_strut[:,usestationStrut,8,1]))")
            println("minimum strut stress: $(minimum(stress_U_strut[:,usestationStrut,8,1]))")
        end
        printSF(verbosity,SF_ult_U_strut,SF_buck_U_strut,composite_station_idx_U_strut, composite_station_name_U_strut,length(strut_precompinput),lam_U_strut,topDamage_strut_U,delta_t,total_t;useStation=usestationStrut)

        stress_L_strut = zeros(N_ts,length(strut_precompinput),length(lam_U_strut[1,:]),3)
        SF_ult_L_strut = zeros(N_ts,length(strut_precompinput),length(lam_L_strut[1,:]))
        SF_buck_L_strut = zeros(N_ts,length(strut_precompinput),length(lam_L_strut[1,:]))

        topstrainout_strut_L,topDamage_strut_L = calcSF(total_t,stress_L_strut,SF_ult_L_strut,SF_buck_L_strut,length(strut_precompinput),plyprops_strut,
        strut_precompinput,strut_precompoutput,lam_L_strut,eps_x_strut[1,:,:],eps_z_strut[1,:,:],eps_y_strut[1,:,:],kappa_x_strut[1,:,:],
        kappa_y_strut[1,:,:],kappa_z_strut[1,:,:],numadIn_strut;failmethod = "maxstress",upper=false,calculate_fatigue)

        if verbosity>0
            println("\n\nLOWER STRUT SURFACE")
        end
        if usestationStrut !=0
            println("maximum strut stress: $(maximum(stress_L_strut[:,usestationStrut,8,1]))")
            println("minimum strut stress: $(minimum(stress_L_strut[:,usestationStrut,8,1]))")
        end
        printSF(verbosity,SF_ult_L_strut,SF_buck_L_strut,composite_station_idx_L_strut, composite_station_name_L_strut,length(strut_precompinput),lam_L_strut,topDamage_strut_L,delta_t,total_t;useStation=usestationStrut)
    else
        topDamage_strut_U = nothing
        topDamage_strut_L  = nothing
        stress_U_strut = nothing
        SF_ult_U_strut = nothing
        SF_buck_U_strut = nothing
        stress_L_strut = nothing
        SF_ult_L_strut = nothing
        SF_buck_L_strut = nothing
    end

    ##########################################
    #### Calculate Stress At the Tower
    ##########################################

    stress_TU = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]),3)
    SF_ult_TU = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]))
    SF_buck_TU = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]))

    topstrainout_tower_U,topDamage_tower_U = calcSF(total_t,stress_TU,SF_ult_TU,SF_buck_TU,length(twr_precompinput),plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr[1,:,:],eps_z_twr[1,:,:],eps_y_twr[1,:,:],kappa_x_twr[1,:,:],
    kappa_y_twr[1,:,:],kappa_z_twr[1,:,:],numadIn_twr;failmethod = "maxstress",upper=true,calculate_fatigue)

    println("\n\nUPPER TOWER")
    printsf_twr(verbosity,lam_U_twr,SF_ult_TU,SF_buck_TU,length(twr_precompinput),Twr_LE_U_idx,topDamage_tower_U,delta_t,total_t)

    stress_TL = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]),3)
    SF_ult_TL = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]))
    SF_buck_TL = zeros(N_ts,length(twr_precompinput),length(lam_U_twr[1,:]))

    topstrainout_tower_L,topDamage_tower_L = calcSF(total_t,stress_TL,SF_ult_TL,SF_buck_TL,length(twr_precompinput),plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr[1,:,:],eps_z_twr[1,:,:],eps_y_twr[1,:,:],kappa_x_twr[1,:,:],
    kappa_y_twr[1,:,:],kappa_z_twr[1,:,:],numadIn_twr;failmethod = "maxstress",upper=false,calculate_fatigue)

    println("\n\nLower TOWER")
    printsf_twr(verbosity,lam_L_twr,SF_ult_TL,SF_buck_TL,length(twr_precompinput),Twr_LE_L_idx,topDamage_tower_L,delta_t,total_t)

    ##########################################
    #### Calculate Mass
    ##########################################

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

    _,turb_masskg = calcMass(myel.props,myort)
    if verbosity>1
        println("\nMass of Turbine: $turb_masskg kg")
    end

    return turb_masskg,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,
    stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,
    topstrainout_blade_L,topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
    topDamage_blade_L,topDamage_tower_U,topDamage_tower_L,topDamage_strut_U,topDamage_strut_L,
    stress_U_strut,SF_ult_U_strut,SF_buck_U_strut,stress_L_strut,SF_ult_L_strut,SF_buck_L_strut
end

"""
    PostProcessOptions

Options structure for controlling post-processing behavior in the OWENS code.

# Fields
* `verbosity::Int`: Controls the level of output verbosity (0=none, 1=warnings, 2=summary, 3=detailed)
* `usestationBld::Int`: Index of the blade station to use for detailed output (0=all stations)
* `usestationStrut::Int`: Index of the strut station to use for detailed output (0=all stations)
* `throwawayTimeSteps::Int`: Number of initial time steps to discard from analysis
* `calculate_fatigue::Bool`: Whether to calculate fatigue damage

# Constructor
```julia
PostProcessOptions(;verbosity=2,
    usestationBld=0,
    usestationStrut=0,
    throwawayTimeSteps=1,
    calculate_fatigue=true)
```
"""
struct PostProcessOptions
    verbosity
    usestationBld
    usestationStrut
    throwawayTimeSteps
    calculate_fatigue
end
function PostProcessOptions(;verbosity=2,
    usestationBld=0,
    usestationStrut=0,
    throwawayTimeSteps=1,
    calculate_fatigue=true)
    return PostProcessOptions(verbosity,usestationBld,usestationStrut,throwawayTimeSteps,calculate_fatigue)
end

function safetyfactor_fatigue(mymesh,components,delta_t; options=PostProcessOptions())
    
    
    verbosity = options.verbosity
    usestationBld = options.usestationBld
    usestationStrut = options.usestationStrut
    throwawayTimeSteps = options.throwawayTimeSteps
    calculate_fatigue = options.calculate_fatigue
    
    for icomp = 1:size(components)[1]

        composite_station_idx_U = []
        composite_station_name_U = []
        composite_station_idx_L = []  
        composite_station_name_L = []
        
        numadIn = components[icomp].nuMadIn
        preCompInput = components[icomp].preCompInput
        preCompOutput = components[icomp].preCompOutput
        plyprops = components[icomp].plyProps
        lam_U = components[icomp].lam_U
        lam_L = components[icomp].lam_L
        e_x = components[icomp].e_x
        e_y = components[icomp].e_y
        e_z = components[icomp].e_z
        k_x = components[icomp].k_x
        k_y = components[icomp].k_y
        k_z = components[icomp].k_z
        
        total_t = length(e_x[1,throwawayTimeSteps:end])*delta_t # for damage rate

        ###############################
        #### Get strain values ########
        ###############################

        N_ts = length(e_x[1,:])
        e_x_numad = zeros(N_ts,length(preCompInput))
        e_z_numad = zeros(N_ts,length(preCompInput))
        e_y_numad = zeros(N_ts,length(preCompInput))
        k_x_numad = zeros(N_ts,length(preCompInput))
        k_y_numad = zeros(N_ts,length(preCompInput))
        k_z_numad = zeros(N_ts,length(preCompInput))
    
        startN = components[icomp].nodeNumbers[1]
        stopN = components[icomp].nodeNumbers[end]

        startE = components[icomp].elNumbers[1]
        stopE = components[icomp].elNumbers[end]

        mesh_spanx = mymesh.x[startN:stopN].-mymesh.x[startN]
        mesh_spany = mymesh.y[startN:stopN].-mymesh.y[startN]
        mesh_spanz = mymesh.z[startN:stopN].-mymesh.z[startN]
        mesh_span = sqrt.(mesh_spanx.^2 .+ mesh_spany.^2 .+ mesh_spanz.^2)
        mesh_span_el = OWENS.safeakima(LinRange(0,1,length(mesh_span)),mesh_span,LinRange(0,1,stopE-startE+1))
        span_numad = numadIn.span/maximum(numadIn.span)*maximum(mesh_span_el)

        for its = 1:N_ts
            # Interpolate to the composite inputs #TODO: verify node vs el in strain
            e_x_numad[its,:] = OWENS.safeakima(mesh_span_el,e_x[:,its],span_numad)
            e_z_numad[its,:] = OWENS.safeakima(mesh_span_el,e_y[:,its],span_numad)
            e_y_numad[its,:] = OWENS.safeakima(mesh_span_el,e_z[:,its],span_numad)
            k_x_numad[its,:] = OWENS.safeakima(mesh_span_el,k_x[:,its],span_numad)
            k_y_numad[its,:] = OWENS.safeakima(mesh_span_el,k_y[:,its],span_numad)
            k_z_numad[its,:] = OWENS.safeakima(mesh_span_el,k_z[:,its],span_numad)
        end

        ##########################################
        #### Calculate Stress
        ##########################################
        stress_U = zeros(N_ts,length(preCompInput),length(lam_U[1,:]),3)
        SF_ult_U = zeros(N_ts,length(preCompInput),length(lam_U[1,:]))
        SF_buck_U = zeros(N_ts,length(preCompInput),length(lam_U[1,:]))

        topstrainout_U,topDamage_U = calcSF(total_t,stress_U,SF_ult_U,SF_buck_U,length(preCompInput),plyprops,
        preCompInput,preCompOutput,lam_U,e_x_numad,e_z_numad,e_y_numad,k_x_numad,
        k_y_numad,k_z_numad,numadIn;failmethod = "maxstress",upper=true,calculate_fatigue)

        if verbosity>0
            println("\n\n$(components[icomp].name) Composite Ultimate and Buckling Safety Factors")
            println("UPPER SURFACE")
        end
        if usestationBld !=0
            println("maximum $(components[icomp].name) stress: $(maximum(stress_U[:,usestationBld,8,1]))")
            println("minimum $(components[icomp].name) stress: $(minimum(stress_U[:,usestationBld,8,1]))")
        end

        OWENS.printSF(verbosity,SF_ult_U,SF_buck_U,composite_station_idx_U, composite_station_name_U,length(preCompInput),lam_U,topDamage_U,delta_t,total_t;useStation=usestationBld)

        stress_L = zeros(N_ts,length(preCompInput),length(lam_L[1,:]),3)
        SF_ult_L = zeros(N_ts,length(preCompInput),length(lam_L[1,:]))
        SF_buck_L = zeros(N_ts,length(preCompInput),length(lam_L[1,:]))

        topstrainout_L,topDamage_L = calcSF(total_t,stress_L,SF_ult_L,SF_buck_L,length(preCompInput),plyprops,
        preCompInput,preCompOutput,lam_L,e_x_numad,e_z_numad,e_y_numad,k_x_numad,
        k_y_numad,k_z_numad,numadIn;failmethod = "maxstress",upper=false,calculate_fatigue)

        if verbosity>0
            println("\n\nLOWER SURFACE")
        end
        if usestationBld !=0
            println("maximum $(components[icomp].name) stress: $(maximum(stress_L[:,usestationBld,8,1]))")
            println("minimum $(components[icomp].name) stress: $(minimum(stress_L[:,usestationBld,8,1]))")
        end
        OWENS.printSF(verbosity,SF_ult_L,SF_buck_L,composite_station_idx_L, composite_station_name_L,length(preCompInput),lam_L,topDamage_L,delta_t,total_t;useStation=usestationBld)

        components[icomp].stress_U = stress_U
        components[icomp].stress_L = stress_L
        components[icomp].strain_U = topstrainout_U
        components[icomp].strain_L = topstrainout_L
        components[icomp].ultsafetyfactor_U = SF_ult_U
        components[icomp].ultsafetyfactor_L = SF_ult_L
        components[icomp].bucksafetyfactor_U = SF_buck_U
        components[icomp].bucksafetyfactor_L = SF_buck_L
        components[icomp].damage_U = topDamage_U
        components[icomp].damage_L = topDamage_L

    end

    return components
end