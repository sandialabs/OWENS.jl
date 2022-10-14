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
    kappa_y,kappa_z,numadIn;failmethod = "maxstress",upper=true,layer=-1)
    for i_station = 1:length(composites_span)
        # i_station = 1
        ibld = 1

        for j_lam = 1:length(lam_in[i_station,:])

            # Get thickness
            refY = precompinput[i_station].le_loc*precompinput[i_station].chord
            thickness_precomp_lag_te = (precompinput[i_station].chord-(refY+precompoutput[i_station].y_sc))
            thickness_precomp_lag_le = precompinput[i_station].chord - thickness_precomp_lag_te
            crossover_fraction = thickness_precomp_lag_le/precompinput[i_station].chord
            idx_y = round(Int,j_lam/length(lam_in[i_station,:])*length(precompinput[i_station].ynode)/2) #TODO: lower
            thickness_precomp_flap = precompinput[i_station].ynode[idx_y]*precompinput[i_station].chord - precompoutput[i_station].x_sc
            if upper
                offsetz = thickness_precomp_flap
            else
                offsetz = -thickness_precomp_flap
            end

            if (j_lam-1)/length(lam_in[i_station,:]) < crossover_fraction
                offsety = thickness_precomp_lag_le * 1-((j_lam-1)/length(lam_in[i_station,:])/crossover_fraction)
            else
                offsety = -thickness_precomp_lag_te * ((j_lam)/length(lam_in[i_station,:])-(1-crossover_fraction))/crossover_fraction
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

                materials = [plyprops.plies[imat] for imat in lam.matid]
                q = Composites.getQ.(materials,lam.theta)

                lowerplystrain, upperplystrain = my_getplystrain(lam, resultantstrain, offset)
                upperplystress = q.*upperplystrain #TODO: bootom side of plies needed for inter-laminate failure?
                out = Composites.getmatfail.(upperplystress,materials,failmethod)
                fail = [out[iii][1] for iii = 1:length(out)]
                sf = [out[iii][2] for iii = 1:length(out)]

                if failmethod == "maxstress"
                    SF_ult[its,i_station,j_lam] = minimum([minimum(abs.(sf[i])) for i = 1:length(sf)])
                else
                    if layer == -1
                        SF_ult[its,i_station,j_lam] = minimum(sf) #TODO; use findmin to identify which layer is failing
                    else
                        SF_ult[its,i_station,j_lam] = minimum(sf[layer])
                    end
                end
                stress[its,i_station,j_lam,:] = upperplystress[2]

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
end


function extractSF(bld_precompinput,bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
    mymesh,myel,myort,Nbld,epsilon_x_hist_ps,kappa_y_hist_ps,kappa_z_hist_ps,epsilon_z_hist_ps,kappa_x_hist_ps,epsilon_y_hist_ps,
    epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1;verbosity=2)

    # Linearly Superimpose the Strains
    epsilon_x_hist = zero(epsilon_x_hist_ps)
    kappa_y_hist = zero(kappa_y_hist_ps)
    kappa_z_hist = zero(kappa_z_hist_ps)
    epsilon_z_hist = zero(epsilon_z_hist_ps)
    kappa_x_hist = zero(kappa_x_hist_ps)
    epsilon_y_hist = zero(epsilon_y_hist_ps)

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

    calcSF(stress_U,SF_ult_U,SF_buck_U,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_U_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "tsaiwu",upper=true)

    if verbosity>0
        println("\nUPPER SURFACE")
    end
    if !isempty(SF_buck_U[SF_buck_U.>0.0])
        SF_buck_U[SF_buck_U.<0.0] .= 1e6
        SF_buck_U[:,:,1] .= 1e6 #ignore leading edge
        SF_buck_U[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_U)
        if verbosity>0
            println("Worst buckling safety factor $(minbuck_sf)")
        end
        if verbosity>0
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_bld)) at lam $(minbuck_sfidx[3]) of $(length(lam_U_bld[minbuck_sfidx[2],:]))")
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end
    if verbosity>0
        println("Minimum Safety Factor on Blade Surface: $(minimum(SF_ult_U))")
    end
    mymin,idx = findmin(SF_ult_U)
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld)) at lam $(idx[3]) of $(length(lam_U_bld[idx[2],:]))")
    end
    mymin,idx = findmin(SF_ult_U[:,:,3])
    if verbosity>0
        println("Spar Cap SF min: $mymin")
    end
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end
    mymin,idx = findmin(SF_ult_U[:,:,1])
    if verbosity>0
        println("Leading Edge SF min: $mymin")
    end
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Upper Spar")
    end
    for SF in SF_ult_U[idx[1],:,3]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Upper Leading Edge")
    end
    for SF in SF_ult_U[idx[1],:,1]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Upper Trailing Edge")
    end
    for SF in SF_ult_U[idx[1],:,6]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Upper Fore Panel")
    end
    for SF in SF_ult_U[idx[1],:,2]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Upper Aft Panel")
    end
    for SF in SF_ult_U[idx[1],:,5]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Upper Buckling")
    end
    for istation = 1:length(composites_span_bld)
        if verbosity>1
            println(minimum(SF_buck_U[idx[1],istation,:]))
        end
    end


    stress_L = zeros(N_ts,length(composites_span_bld),length(lam_U_bld[1,:]),3)
    SF_ult_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))
    SF_buck_L = zeros(N_ts,length(composites_span_bld),length(lam_L_bld[1,:]))

    calcSF(stress_L,SF_ult_L,SF_buck_L,composites_span_bld,plyprops_bld,
    bld_precompinput,bld_precompoutput,lam_L_bld,eps_x_bld,eps_z_bld,eps_y_bld,kappa_x_bld,
    kappa_y_bld,kappa_z_bld,numadIn_bld;failmethod = "tsaiwu",upper=false)

    if verbosity>0
        println("\nLOWER SURFACE")
    end
    if !isempty(SF_buck_L[SF_buck_L.>0.0])
        SF_buck_L[SF_buck_L.<0.0] .= 1e6
        SF_buck_L[:,:,1] .= 1e6 #ignore leading edge
        SF_buck_L[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_L)
        if verbosity>0
            println("Worst buckling safety factor $(minbuck_sf)")
        end
        if verbosity>0
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_bld)) at lam $(minbuck_sfidx[3]) of $(length(lam_L_bld[minbuck_sfidx[2],:]))")
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end
    if verbosity>0
        println("Minimum Safety Factor on Blade Surface: $(minimum(SF_ult_L))")
    end
    mymin,idx = findmin(SF_ult_L)
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld)) at lam $(idx[3]) of $(length(lam_L_bld[idx[2],:]))")
    end
    mymin,idx = findmin(SF_ult_L[:,:,3])
    if verbosity>0
        println("Spar Cap SF min: $mymin")
    end
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end
    mymin,idx = findmin(SF_ult_L[:,:,1])
    if verbosity>0
        println("Leading Edge SF min: $mymin")
    end
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_bld))")
    end

    if verbosity>1
        println("Lower Spar")
    end
    for SF in SF_ult_L[idx[1],:,3]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Lower Leading Edge")
    end
    for SF in SF_ult_L[idx[1],:,1]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Lower Trailing Edge")
    end
    for SF in SF_ult_L[idx[1],:,6]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Lower Fore Panel")
    end
    for SF in SF_ult_L[idx[1],:,2]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Lower Aft Panel")
    end
    for SF in SF_ult_L[idx[1],:,5]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Lower Buckling")
    end
    for istation = 1:length(composites_span_bld)
        if verbosity>1
            println(minimum(SF_buck_L[idx[1],istation,:]))
        end
    end

    ##########################################
    #### Calculate Stress At the Tower
    ##########################################

    stress_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
    SF_ult_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
    SF_buck_TU = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

    calcSF(stress_TU,SF_ult_TU,SF_buck_TU,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "tsaiwu",upper=true)


    if verbosity>0
        println("\nUPPER TOWER")
    end
    if !isempty(SF_buck_TU[SF_buck_TU.>0.0])
        SF_buck_TU[SF_buck_TU.<0.0] .= 1e6
        # SF_buck_TU[:,:,1] .= 1e6 #ignore leading edge
        # SF_buck_TU[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_TU)
        if verbosity>0
            println("Worst buckling safety factor $(minbuck_sf)")
        end
        if verbosity>0
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_U_twr[minbuck_sfidx[2],:]))")
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end
    if verbosity>0
        println("Minimum Safety Factor on tower Surface: $(minimum(SF_ult_TU))")
    end
    mymin,idx = findmin(SF_ult_TU)
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_U_twr[idx[2],:]))")
    end

    if verbosity>1
        println("Leading Edge")
    end
    for SF in SF_ult_TU[idx[1],:,1]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Buckling")
    end
    for istation = 1:length(composites_span_twr)
        if verbosity>1
            println(minimum(SF_buck_TU[minbuck_sfidx[1],istation,:]))
        end
    end


    stress_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]),3)
    SF_ult_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))
    SF_buck_TL = zeros(N_ts,length(composites_span_twr),length(lam_U_twr[1,:]))

    calcSF(stress_TL,SF_ult_TL,SF_buck_TL,composites_span_twr,plyprops_twr,
    twr_precompinput,twr_precompoutput,lam_U_twr,eps_x_twr,eps_z_twr,eps_y_twr,kappa_x_twr,
    kappa_y_twr,kappa_z_twr,numadIn_twr;failmethod = "tsaiwu",upper=false)


    if verbosity>0
        println("\nLOWER TOWER")
    end
    if !isempty(SF_buck_TL[SF_buck_TL.>0.0])
        SF_buck_TL[SF_buck_TL.<0.0] .= 1e6
        # SF_buck_TL[:,:,1] .= 1e6 #ignore leading edge
        # SF_buck_TL[:,:,6] .= 1e6 #ignore trailing edge
        minbuck_sf,minbuck_sfidx = findmin(SF_buck_TL)
        if verbosity>0
            println("Worst buckling safety factor $(minbuck_sf)")
        end
        if verbosity>0
            println("At time $(minbuck_sfidx[1]*0.05)s at composite station $(minbuck_sfidx[2]) of $(length(composites_span_twr)) at lam $(minbuck_sfidx[3]) of $(length(lam_L_twr[minbuck_sfidx[2],:]))")
        end
    else
        if verbosity>0
            println("Buckling not a factor, no sections in compression")
        end
    end
    if verbosity>0
        println("Minimum Safety Factor on tower Surface: $(minimum(SF_ult_TL))")
    end
    mymin,idx = findmin(SF_ult_TL)
    if verbosity>0
        println("At time $(idx[1]*0.05)s at composite station $(idx[2]) of $(length(composites_span_twr)) at lam $(idx[3]) of $(length(lam_L_twr[idx[2],:]))")
    end

    if verbosity>1
        println("Leading Edge")
    end
    for SF in SF_ult_TL[idx[1],:,1]
        if verbosity>1
            println(SF)
        end
    end

    if verbosity>1
        println("Buckling")
    end
    for istation = 1:length(composites_span_twr)
        if verbosity>1
            println(minimum(SF_buck_TL[minbuck_sfidx[1],istation,:]))
        end
    end

    ##########################################
    #### Calculate Mass
    ##########################################

    function massOWENS(sectionPropsArray,myort)
        mass = 0.0
        for (i,sectionProp) in enumerate(sectionPropsArray)
            lenEl = 0
            try
                lenEl = myort.Length[i]
            catch
                lenEl = myort.Length[i-1]
            end
            rhoA = sectionProp.rhoA[1]
            mass += lenEl*rhoA
        end
        return mass
    end

    massOwens = massOWENS(myel.props,myort)
    if verbosity>0
        println("Mass ARCUS 5MW Turbine: $massOwens")
    end
    if verbosity>1
        println("Mass Baseline 5MW: 181888.58246407832")
        println("That's a $((181888.58246407832-massOwens)/181888.58246407832*100)% reduction in mass")
    end

    return massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL
end
