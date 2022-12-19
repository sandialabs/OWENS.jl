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

##########################################
#### Composite Failure & Buckling ########
##########################################

function calcSF(stress,SF_ult,SF_buck,composites_span,plyprops,
    precompinput,precompoutput,lam_in,eps_x,eps_z,eps_y,kappa_x,
    kappa_y,kappa_z,numadIn;failmethod = "maxstress",CLT=false,upper=true,layer=-1)
    topstrainout = zeros(length(eps_x[1,:,1]),length(composites_span),length(lam_in[1,:]),9) # time, span, lam, x,y, Assumes you use zero plies for sections that aren't used
    for i_station = 1:length(composites_span)
        # i_station = 1
        ibld = 1

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

            for its = 1:length(eps_x[1,:,1])
                # j_lam = 5
                # its = 1
                # println("its $its j_lam $j_lam i_station $i_station")
                resultantstrain = [
                eps_x[ibld,its,i_station],
                eps_z[ibld,its,i_station],
                eps_y[ibld,its,i_station],
                kappa_x[ibld,its,i_station],
                kappa_y[ibld,its,i_station],
                kappa_z[ibld,its,i_station]]

                lam = lam_in[i_station,j_lam] #TODO: make this as an input

                mat_idx = [numadIn.stack_mat_types[imat] for imat in lam.matid] # Map from the stack number to the material number
                materials = [plyprops.plies[imat] for imat in mat_idx] # Then get the actual material used
                
                # Map the stack number to the material number

                q = Composites.getQ.(materials,lam.theta)

                #TODO: if we are doing the lower surface, should we use the lower plystrain?
                lowerplystrain, upperplystrain = my_getplystrain(lam, resultantstrain, offset)
                topstrainout[its,i_station,j_lam,1:3] = upperplystrain[1]
                topstrainout[its,i_station,j_lam,4:end] = resultantstrain[:]
                upperplystress = q.*upperplystrain #TODO: bootom side of plies needed for inter-laminate failure?
                out = Composites.getmatfail.(upperplystress,materials,failmethod)
                fail = [out[iii][1] for iii = 1:length(out)]
                sf = [out[iii][2] for iii = 1:length(out)]

                if failmethod == "maxstress"
                    SF_ult[its,i_station,j_lam] = minimum(minimum([sf[isf][sf[isf].>0.0] for isf = 1:length(sf)])) # Pick out the worst case safety factor that is positive - negative means it is in the wrong direction for the failure criteria
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
                        buck_layer = 2
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
        end
    end
    return topstrainout
end

function printSF(verbosity,SF_ult,SF_buck,LE_idx,TE_idx,SparCap_idx,ForePanel_idx,AftPanel_idx,composites_span_bld,lam_used)
    #Ultimate
    mymin,idx = findmin(SF_ult)
    if verbosity>0
        println("\nMinimum Safety Factor on Surface: $(minimum(SF_ult))")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld)) at lam $(idx[3]) of $(length(lam_used[idx[2],:]))")
    end
    #Buckling
    if !isempty(SF_buck[SF_buck.>0.0])
        SF_buck[SF_buck.<0.0] .= 1e6
        SF_buck[:,:,LE_idx] .= 1e6 #ignore leading edge
        SF_buck[:,:,TE_idx] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck)
        if verbosity>0
            println("\nWorst buckling safety factor $(minbuck_sf)")
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_bld)) at lam $(minbuck_sfidx[3]) of $(length(lam_used[minbuck_sfidx[2],:]))")
        end
        if verbosity>1
            println("Buckling")
            for istation = 1:length(composites_span_bld)
                println(minimum(SF_buck[minbuck_sfidx[1],istation,:]))
            end
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end

    mymin,idx = findmin(SF_ult[:,:,SparCap_idx])
    if verbosity>0
        println("\nSpar Cap SF min: $mymin")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("\nSpar")
        for SF in SF_ult[idx[1],:,SparCap_idx]
            println(SF)
        end
    end

    mymin,idx = findmin(SF_ult[:,:,LE_idx])
    if verbosity>0
        println("\nLeading Edge SF min: $mymin")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Leading Edge")
        for SF in SF_ult[idx[1],:,LE_idx]
            println(SF)
        end
    end

    mymin,idx = findmin(SF_ult[:,:,TE_idx])
    if verbosity>0
        println("\nTrailling Edge SF min: $mymin")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Trailing Edge")
        for SF in SF_ult[idx[1],:,TE_idx]
            println(SF)
        end
    end

    mymin,idx = findmin(SF_ult[:,:,ForePanel_idx])
    if verbosity>0
        println("\nFore Panel SF min: $mymin")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Fore Panel")
        for SF in SF_ult[idx[1],:,ForePanel_idx]
            println(SF)
        end
    end

    mymin,idx = findmin(SF_ult[:,:,AftPanel_idx])
    if verbosity>0
        println("\nAft Panel SF min: $mymin")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Aft Panel")
        for SF in SF_ult[idx[1],:,AftPanel_idx]
            println(SF)
        end
    end
end


function extractSF(bld_precompinput,bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,Nbld,epsilon_x_hist_ps,kappa_y_hist_ps,kappa_z_hist_ps,epsilon_z_hist_ps,kappa_x_hist_ps,epsilon_y_hist_ps;
    epsilon_x_hist_1=nothing,kappa_y_hist_1=nothing,kappa_z_hist_1=nothing,
    epsilon_z_hist_1=nothing,kappa_x_hist_1=nothing,epsilon_y_hist_1=nothing,verbosity=2,
    LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
    LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
    Twr_LE_U_idx=1,Twr_LE_L_idx=1)

    # Linearly Superimpose the Strains
    epsilon_x_hist = copy(epsilon_x_hist_ps).*0.0
    kappa_y_hist = copy(kappa_y_hist_ps).*0.0
    kappa_z_hist = copy(kappa_z_hist_ps).*0.0
    epsilon_z_hist = copy(epsilon_z_hist_ps).*0.0
    kappa_x_hist = copy(kappa_x_hist_ps).*0.0
    epsilon_y_hist = copy(epsilon_y_hist_ps).*0.0

    if epsilon_x_hist_1!=nothing
        for ipt = 1:4
            for jts = 1:length(epsilon_x_hist_ps[1,1,:])
                epsilon_x_hist[ipt,:,jts] = epsilon_x_hist_ps[ipt,:,jts] + epsilon_x_hist_1[ipt,:,end]
                kappa_y_hist[ipt,:,jts] = kappa_y_hist_ps[ipt,:,jts] + kappa_y_hist_1[ipt,:,end]
                kappa_z_hist[ipt,:,jts] = kappa_z_hist_ps[ipt,:,jts] + kappa_z_hist_1[ipt,:,end]
                epsilon_z_hist[ipt,:,jts] = epsilon_z_hist_ps[ipt,:,jts] + epsilon_z_hist_1[ipt,:,end]
                kappa_x_hist[ipt,:,jts] = kappa_x_hist_ps[ipt,:,jts] + kappa_x_hist_1[ipt,:,end]
                epsilon_y_hist[ipt,:,jts] = epsilon_y_hist_ps[ipt,:,jts] + epsilon_y_hist_1[ipt,:,end]
            end
        end
    else
        epsilon_x_hist = epsilon_x_hist_ps
        kappa_y_hist = kappa_y_hist_ps
        kappa_z_hist = kappa_z_hist_ps
        epsilon_z_hist = epsilon_z_hist_ps
        kappa_x_hist = kappa_x_hist_ps
        epsilon_y_hist = epsilon_y_hist_ps
    end

    ##########################################
    #### Get strain values at the blades #####
    ##########################################

    meanepsilon_z_hist = Statistics.mean(epsilon_z_hist,dims=1)
    meanepsilon_y_hist = Statistics.mean(epsilon_y_hist,dims=1)

    N_ts = length(epsilon_x_hist[1,1,:])
    eps_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    eps_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    eps_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_x_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_y_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    kappa_z_bld = zeros(Nbld,N_ts,length(bld_precompinput))
    mesh_span_bld = zeros(length(mymesh.structuralNodeNumbers[1,:]))
    composites_span_bld = zeros(length(bld_precompinput))
    for ibld = 1:Nbld
        start = Int(mymesh.structuralNodeNumbers[ibld,1])
        stop = Int(mymesh.structuralNodeNumbers[ibld,end])
        x = mymesh.z[start:stop]
        x = x.-x[1] #zero
        x = x./x[end] #normalize
        mesh_span_bld[:] = mymesh.z[start:stop].-mymesh.z[start]
        global composites_span_bld = FLOWMath.akima(LinRange(0,1,length(mesh_span_bld)),mesh_span_bld,LinRange(0,1,length(bld_precompinput)))
        for its = 1:N_ts
            # Interpolate to the composite inputs
            eps_x_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,epsilon_x_hist[1,start:stop,its],composites_span_bld)
            eps_z_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,meanepsilon_z_hist[1,start:stop,its],composites_span_bld)
            eps_y_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,meanepsilon_y_hist[1,start:stop,its],composites_span_bld)
            kappa_x_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_x_hist[1,start:stop,its],composites_span_bld)
            kappa_y_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_y_hist[1,start:stop,its],composites_span_bld)
            kappa_z_bld[ibld,its,:] = FLOWMath.akima(mesh_span_bld,kappa_z_hist[1,start:stop,its],composites_span_bld)
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
    composites_span_twr = FLOWMath.akima(LinRange(0,1,length(mesh_span_twr)),mesh_span_twr,LinRange(0,1,length(twr_precompinput)))
    for its = 1:N_ts
        # Interpolate to the composite inputs
        eps_x_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,epsilon_x_hist[1,start:stop,its],composites_span_twr)
        eps_z_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,meanepsilon_z_hist[1,start:stop,its],composites_span_twr)
        eps_y_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,meanepsilon_y_hist[1,start:stop,its],composites_span_twr)
        kappa_x_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_x_hist[1,start:stop,its],composites_span_twr)
        kappa_y_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_y_hist[1,start:stop,its],composites_span_twr)
        kappa_z_twr[1,its,:] = FLOWMath.akima(mesh_span_twr,kappa_z_hist[1,start:stop,its],composites_span_twr)

    end

    ##########################################
    #### Calculate Stress At the Blades
    ##########################################

    stress_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]),3)
    SF_ult_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]))
    SF_buck_U = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]))

    topstrainout_blade_U = calcSF(stress_U,SF_ult_U,SF_buck_U,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_U_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "maxstress",upper=true)

    if verbosity>0
        println("Composite Ultimate and Buckling Safety Factors")
        println("\n\nUPPER BLADE SURFACE")
    end
    printSF(verbosity,SF_ult_U,SF_buck_U,LE_U_idx,TE_U_idx,SparCapU_idx,ForePanelU_idx,AftPanelU_idx,composites_span_bld,lam_U_bld)

    stress_L = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]),3)
    SF_ult_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))
    SF_buck_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))

    topstrainout_blade_L = calcSF(stress_L,SF_ult_L,SF_buck_L,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_L_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "maxstress",upper=false)

    if verbosity>0
        println("\n\nLOWER BLADE SURFACE")
    end
    printSF(verbosity,SF_ult_L,SF_buck_L,LE_L_idx,TE_L_idx,SparCapU_idx,ForePanelU_idx,AftPanelU_idx,composites_span_bld,lam_L_bld)

    ##########################################
    #### Calculate Stress At the Tower
    ##########################################

    stress_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
    SF_ult_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
    SF_buck_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

    calcSF(stress_TU,SF_ult_TU,SF_buck_TU,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "maxstress",upper=true)


    if verbosity>0
        println("\n\nUPPER TOWER")
    end
    #Ultimate
    mymin,idx = findmin(SF_ult_TU)
    if verbosity>0
        println("\nMinimum Safety Factor on tower Surface: $(minimum(SF_ult_TU))")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_U_twr[idx[2],:]))")
    end
    #Buckling
    if !isempty(SF_buck_TU[SF_buck_TU.>0.0])
        SF_buck_TU[SF_buck_TU.<0.0] .= 1e6
        # SF_buck_TU[:,:,1] .= 1e6 #ignore leading edge
        # SF_buck_TU[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_TU)
        if verbosity>0
            println("\nWorst buckling safety factor $(minbuck_sf)")
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_U_twr[minbuck_sfidx[2],:]))")
        elseif verbosity>1
            println("Buckling")
            for istation = 1:length(composites_span_twr)
                println(minimum(SF_buck_TU[minbuck_sfidx[1],istation,:]))
            end
        end
    
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end

    if verbosity>1
        println("\nLeading Edge")
        for SF in SF_ult_TU[idx[1],:,Twr_LE_U_idx]
            println(SF)
        end
    end



    stress_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
    SF_ult_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
    SF_buck_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

    calcSF(stress_TL,SF_ult_TL,SF_buck_TL,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "maxstress",upper=false)


    if verbosity>0
        println("\n\nLOWER TOWER")
    end
    #Ultimate
    mymin,idx = findmin(SF_ult_TL)
    if verbosity>0
        println("\nMinimum Safety Factor on tower Surface: $(minimum(SF_ult_TL))")
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_L_twr[idx[2],:]))")
    end
    #Buckling
    if !isempty(SF_buck_TL[SF_buck_TL.>0.0])
        SF_buck_TL[SF_buck_TL.<0.0] .= 1e6
        # SF_buck_TL[:,:,1] .= 1e6 #ignore leading edge
        # SF_buck_TL[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_TL)
        if verbosity>0
            println("\nWorst buckling safety factor $(minbuck_sf)")
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_L_twr[minbuck_sfidx[2],:]))")
        end
        if verbosity>1
            println("Buckling")
            for istation = 1:length(composites_span_twr)
                println(minimum(SF_buck_TL[minbuck_sfidx[1],istation,:]))
            end
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end

    if verbosity>1
        println("\nLeading Edge")
        for SF in SF_ult_TL[idx[1],:,Twr_LE_L_idx]
            println(SF)
        end
    end

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
    if verbosity>0
        println("\nMass of Turbine: $turb_masskg kg")
    end

    return turb_masskg,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,
    stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L
end
