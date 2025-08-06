########################################
############ GXBeam Setup ##############
########################################

const e1 = SVector(1, 0, 0)
const e2 = SVector(0, 1, 0)
const e3 = SVector(0, 0, 1)
function owens_to_gx(
    mymesh,
    myort,
    myjoint,
    sectionPropsArray,
    stiff_array,
    mass_array;
    VTKmeshfilename = nothing,
    damp_coef = 0.05,
    inf_stiff = false,
)

    # Space out the mesh so no two points are exactly on top of one another.
    Noverlap = 1e8
    points = [[mymesh.x[i], mymesh.y[i], mymesh.z[i]] for i = 1:length(mymesh.x)]
    # while Noverlap > 1
    #     Noverlap = 1
    #     for ipt = 2:length(points)
    #         indices = findall(x->isapprox(points[ipt],x;atol=1e-5),points)
    #         if !isempty(indices) && length(indices) > 1
    #             # println("$Noverlap")
    #             points[ipt].+=1e-5
    #             Noverlap += 1
    #         end
    #     end
    # end
    start = Int.(vcat(mymesh.conn[:, 1], myjoint[:, 2]))
    stop = Int.(vcat(mymesh.conn[:, 2], myjoint[:, 3]))

    NT = promote_type(eltype(myort.Psi_d), eltype(myort.Theta_d), eltype(myort.Twist_d))

    frames = Vector{Matrix{NT}}(undef, length(start))
    xaf = zeros(NT, length(start), length(sectionPropsArray[1].xaf))
    yaf = zeros(NT, length(start), length(sectionPropsArray[1].xaf))

    for ielem = 1:length(start)

        nodeNum = start[ielem] # Get node number
        elNum = findfirst(x->x==nodeNum, mymesh.conn[:, 1]) # Get element number
        if isnothing(elNum) || elNum>length(sectionPropsArray)
            elNum = findfirst(x->x==nodeNum, mymesh.conn[:, 2]) # Get element number
        end
        yaw = (myort.Psi_d[elNum]) * pi/180 #deg to rad
        pitch = (myort.Theta_d[elNum]) * pi/180 #deg to rad
        roll = myort.Twist_d[elNum] * pi/180 #deg to rad
        cy = cos(yaw)
        cp = cos(pitch)
        cr = cos(roll)
        sy = sin(yaw)
        sp = sin(pitch)
        sr = sin(roll)
        frame_y = [
            cy -sy 0
            sy cy 0
            0 0 1
        ]
        frame_p = [
            cp 0 sp
            0 1 0
            -sp 0 cp
        ]
        frame_r = [
            1 0 0
            0 cr -sr
            0 sr cr
        ]

        frames[ielem] = frame_y * frame_p * frame_r
        if ielem > length(mymesh.conn[:, 1]) # joint frame of reference
            frames[ielem] = [1 0 0; 0 1 0; 0 0 1]
        end

        xaf[ielem, :] .= sectionPropsArray[elNum].xaf
        yaf[ielem, :] .= sectionPropsArray[elNum].yaf
    end

    #Combine everything
    Nremain_joints = length(myjoint[:, 2])
    stiff_joints = fill(zeros(6, 6), Nremain_joints) #infinite mass and stiffness for joints
    mass_joints = fill(zeros(6, 6), Nremain_joints)

    stiff = vcat(stiff_array, stiff_joints)
    mass = vcat(mass_array, mass_joints)

    if inf_stiff
        stiff = zero.(stiff)
        mass = zero.(mass)
    end

    damping = fill(
        [damp_coef, damp_coef, damp_coef, damp_coef, damp_coef, damp_coef],
        length(mass),
    )

    # create assembly of interconnected nonlinear beams
    assembly = GXBeam.Assembly(
        points,
        start,
        stop;
        stiffness = stiff,
        frames = frames,
        damping,
        mass,
    )

    system = GXBeam.DynamicSystem(assembly)#;prescribed_points=keys(prescribed_conditionsidx))

    sections = zeros(NT, 3, length(sectionPropsArray[1].xaf), length(assembly.points))

    for ipt = 1:length(assembly.points)

        elNum = findfirst(x->x==ipt, mymesh.conn[:, 1]) # Get element number
        if isnothing(elNum)
            elNum = findfirst(x->x==ipt, mymesh.conn[:, 2]) # Get element number
        end

        yaw = (myort.Psi_d[elNum]) * pi/180 #deg to rad
        pitch = (myort.Theta_d[elNum]) * pi/180 #deg to rad
        roll = myort.Twist_d[elNum] * pi/180 #deg to rad
        cy = cos(yaw)
        cp = cos(pitch)
        cr = cos(roll)
        sy = sin(yaw)
        sp = sin(pitch)
        sr = sin(roll)
        frame_y = [
            cy -sy 0
            sy cy 0
            0 0 1
        ]
        frame_p = [
            cp 0 sp
            0 1 0
            -sp 0 cp
        ]
        frame_r = [
            1 0 0
            0 cr -sr
            0 sr cr
        ]

        frames2 = frame_y * frame_p * frame_r

        myaf =
            frames2*[zero(xaf[elNum, :]) (xaf[elNum, :] .- maximum(xaf[elNum, :])/2) yaf[
                elNum,
                :,
            ]]'
        sections[1, :, ipt] = myaf[1, :]
        sections[2, :, ipt] = myaf[2, :]
        sections[3, :, ipt] = myaf[3, :]

    end

    if !isnothing(VTKmeshfilename)
        try #this should error if someone on windows uses backslash '\'
            lastforwardslash = findlast(x->x=='/', VTKmeshfilename)
            filepath = VTKmeshfilename[1:(lastforwardslash-1)]
            if !isdir(filepath)
                mkdir(filepath)
            end
        catch
            @info "Please manually create the directory to house $VTKmeshfilename"
        end
        mywrite_vtk(VTKmeshfilename, assembly; sections)
    end

    return system, assembly, sections, frames, points, start, stop, stiff, mass
end


################################################################
################ SAVE VTK TIME DOMAIN OUTPUT ###################
################################################################
"""

    OWENSVTK(VTKsaveName,rundata,system,assembly,sections,mymesh,myel;tsave_idx=1:length(rundata.t,))

Formats and outputs OWENS data into VTK format

#Intput

#Output
* `none`:

"""
function OWENSVTK(
    VTKsaveName,
    rundata,
    system,
    assembly,
    sections,
    mymesh,
    myel;
    stress = nothing,
    tsave_idx = 1:length(rundata.t),
)

    return OWENSVTK(
        VTKsaveName,
        rundata.t,
        rundata.uHist,
        system,
        assembly,
        sections,
        rundata.aziHist,
        mymesh,
        myel,
        rundata.epsilon_x_hist,
        rundata.epsilon_y_hist,
        rundata.epsilon_z_hist,
        rundata.kappa_x_hist,
        rundata.kappa_y_hist,
        rundata.kappa_z_hist,
        rundata.FReactionHist,
        rundata.topFexternal_hist;
        stress,
        tsave_idx,
    )
end

"""

    OWENSVTK(savename,t,uHist,system,assembly,sections,aziHist,mymesh,myel,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FReactionHist)

Formats and outputs OWENS data into VTK format

#Intput

#Output
* `none`:

"""
function OWENSVTK(
    VTKsaveName,
    t,
    uHist,
    system,
    assembly,
    sections,
    aziHist,
    mymesh,
    myel,
    epsilon_x_hist,
    epsilon_y_hist,
    epsilon_z_hist,
    kappa_x_hist,
    kappa_y_hist,
    kappa_z_hist,
    FReactionHist,
    topFexternal_hist;
    stress = nothing,
    tsave_idx = 1:length(t),
)

    println("Saving VTK time domain files")
    userPointNames=[
        "EA",
        "rhoA",
        "EIyy",
        "EIzz",
        "e_x",
        "e_y",
        "e_z",
        "k_x",
        "k_y",
        "k_z",
        "Fx_Reaction",
        "Fy_Reaction",
        "Fz_Reaction",
        "Mx_Reaction",
        "My_Reaction",
        "Mz_Reaction",
        "Fx_Applied",
        "Fy_Applied",
        "Fz_Applied",
        "Mx_Applied",
        "My_Applied",
        "Mz_Applied",
        "Principal_Surface_Layer_Stress",
    ]#,"Fx","Fy","Fz","Mx","My","Mz"]
    # userPointData[iname,it,ipt] = Float64

    # map el props to points using con
    userPointData = zeros(length(userPointNames), length(tsave_idx), mymesh.numNodes)
    EA_points = zeros(mymesh.numNodes)
    rhoA_points = zeros(mymesh.numNodes)
    EIyy_points = zeros(mymesh.numNodes)
    EIzz_points = zeros(mymesh.numNodes)

    # Time-invariant data
    for iel = 1:length(myel.props)
        # iel = 1
        nodes = mymesh.conn[iel, :]
        EA_points[Int.(nodes)] = myel.props[iel].EA
        rhoA_points[Int.(nodes)] = myel.props[iel].rhoA
        EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
        EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
    end


    epsilon_x_histused = mean(epsilon_x_hist; dims = 1)
    epsilon_y_histused = mean(epsilon_y_hist; dims = 1)
    epsilon_z_histused = mean(epsilon_z_hist; dims = 1)
    kappa_x_histused = mean(kappa_x_hist; dims = 1)
    kappa_y_histused = mean(kappa_y_hist; dims = 1)
    kappa_z_histused = mean(kappa_z_hist; dims = 1)

    # fill in the big matrix
    for it_save = 1:length(tsave_idx)
        it = tsave_idx[it_save]
        userPointData[1, it_save, :] = EA_points
        userPointData[2, it_save, :] = rhoA_points
        userPointData[3, it_save, :] = EIyy_points
        userPointData[4, it_save, :] = EIzz_points
        for iel = 1:length(myel.props)
            nodes = mymesh.conn[iel, :]
            userPointData[5, it_save, Int.(nodes)] .= epsilon_x_histused[1, iel, it]
            userPointData[6, it_save, Int.(nodes)] .= epsilon_y_histused[1, iel, it]
            userPointData[7, it_save, Int.(nodes)] .= epsilon_z_histused[1, iel, it]
            userPointData[8, it_save, Int.(nodes)] .= kappa_x_histused[1, iel, it]
            userPointData[9, it_save, Int.(nodes)] .= kappa_y_histused[1, iel, it]
            userPointData[10, it_save, Int.(nodes)] .= kappa_z_histused[1, iel, it]
        end
        userPointData[11, it_save, :] .= FReactionHist[it, 1:6:end]
        userPointData[12, it_save, :] .= FReactionHist[it, 2:6:end]
        userPointData[13, it_save, :] .= FReactionHist[it, 3:6:end]
        userPointData[14, it_save, :] .= FReactionHist[it, 4:6:end]
        userPointData[15, it_save, :] .= FReactionHist[it, 5:6:end]
        userPointData[16, it_save, :] .= FReactionHist[it, 6:6:end]

        userPointData[17, it_save, :] .= topFexternal_hist[it, 1:6:end]
        userPointData[18, it_save, :] .= topFexternal_hist[it, 2:6:end]
        userPointData[19, it_save, :] .= topFexternal_hist[it, 3:6:end]
        userPointData[20, it_save, :] .= topFexternal_hist[it, 4:6:end]
        userPointData[21, it_save, :] .= topFexternal_hist[it, 5:6:end]
        userPointData[22, it_save, :] .= topFexternal_hist[it, 6:6:end]
    end

    azi=aziHist#./aziHist*1e-6
    # VTKsaveName = "$path/vtk/$(windINPfilename[1:end-4])"
    OWENS.OWENSFEA_VTK(
        VTKsaveName,
        t[tsave_idx],
        uHist[tsave_idx, :],
        system,
        assembly,
        sections;
        scaling = 1,
        azi = azi[tsave_idx],
        userPointNames,
        userPointData,
        stress,
    )

end

function OWENSFEA_VTK(
    filename,
    tvec,
    uHist,
    system,
    assembly,
    sections;
    scaling = 1,
    azi = zero(tvec),
    delta_x = zero(tvec),
    delta_y = zero(tvec),
    delta_z = zero(tvec),
    userPointNames = nothing,
    userPointData = nothing,
    stress = nothing,
)

    history = Vector{GXBeam.AssemblyState{eltype(system)}}(undef, length(tvec))
    uHist = uHist'
    for isave = 1:length(tvec)
        disp_x = [uHist[i, isave] for i = 1:6:length(uHist[:, isave])]
        disp_y = [uHist[i, isave] for i = 2:6:length(uHist[:, isave])]
        disp_z = [uHist[i, isave] for i = 3:6:length(uHist[:, isave])]
        curv_x = [uHist[i, isave] for i = 4:6:length(uHist[:, isave])]
        curv_y = [uHist[i, isave] for i = 5:6:length(uHist[:, isave])]
        curv_z = [uHist[i, isave] for i = 6:6:length(uHist[:, isave])]

        disp_matrix = [disp_x;; disp_y;; disp_z]
        curv_matrix = [curv_x;; curv_y;; curv_z]

        TF = promote_type(eltype(system), eltype(system.x))
        points = Vector{GXBeam.PointState{TF}}(undef, length(assembly.points))

        for ipoint = 1:length(points)

            u = SVector(
                Float64(disp_matrix[ipoint, 1]),
                Float64(disp_matrix[ipoint, 2]),
                Float64(disp_matrix[ipoint, 3]),
            )
            udot = SVector(0.0, 0.0, 0.0)
            theta = SVector(
                Float64(curv_matrix[ipoint, 1]),
                Float64(curv_matrix[ipoint, 2]),
                Float64(curv_matrix[ipoint, 3]),
            )
            thetadot = SVector(0.0, 0.0, 0.0)
            V = SVector(0.0, 0.0, 0.0)
            Vdot = SVector(0.0, 0.0, 0.0)
            Omega = SVector(0.0, 0.0, 0.0)
            Omegadot = SVector(0.0, 0.0, 0.0)
            F = SVector(0.0, 0.0, 0.0)
            M = SVector(0.0, 0.0, 0.0)

            points[ipoint] = GXBeam.PointState(
                u,#::SVector{3, TF}
                udot,
                theta,#::SVector{3, TF}
                thetadot,
                V,
                Vdot,
                Omega,
                Omegadot,
                F,#::SVector{3, TF}
                M,
            )#::SVector{3, TF}
        end

        TF = promote_type(eltype(system), eltype(system.x))
        elements = Vector{GXBeam.ElementState{TF}}(undef, length(assembly.elements))

        for ielem = 1:length(elements)

            u = SVector(0.0, 0.0, 0.0) # Linear displacement
            udot = SVector(0.0, 0.0, 0.0) # Linear displacement rate
            theta = SVector(0.0, 0.0, 0.0) # Angular displacement (Wiener-Milenkovic parameters)
            thetadot = SVector(0.0, 0.0, 0.0) # Angular displacement rate
            V = SVector(0.0, 0.0, 0.0) # Linear velocity
            Vdot = SVector(0.0, 0.0, 0.0) # Linear velocity rate
            Omega = SVector(0.0, 0.0, 0.0) # Angular velocity
            Omegadot = SVector(0.0, 0.0, 0.0) # Angular velocity rate
            Fi = SVector(0.0, 0.0, 0.0) # Internal forces
            Mi = SVector(0.0, 0.0, 0.0) # Internal moments

            elements[ielem] = GXBeam.ElementState(
                u,
                udot,
                theta,
                thetadot,
                V,
                Vdot,
                Omega,
                Omegadot,
                Fi,
                Mi,
            )#::SVector{3, TF}
        end

        history[isave] = GXBeam.AssemblyState(points, elements)
    end

    try #this should error if someone on windows uses backslash '\'
        lastforwardslash = findlast(x->x=='/', filename)
        filepath = filename[1:(lastforwardslash-1)]
        if !isdir(filepath)
            mkdir(filepath)
        end
    catch
        @info "Please manually create the directory to house $filename"
    end

    mywrite_vtk(
        filename,
        assembly,
        history,
        tvec;
        scaling,
        sections,
        theta_z = azi,
        delta_x,
        delta_y,
        delta_z,
        userPointNames,
        userPointData,
        stress,
    )
end

"""
    OWENSFEA_VTK(filename, topData::OWENS.TopData, setupOutputs::OWENS.SetupOutputs; kwargs...)

Wrapper function for OWENSFEA_VTK that takes TopData and SetupOutputs structures directly.
This function unpacks the required data from the structures and calls the original OWENSFEA_VTK function.

# Arguments
- `filename`: Output filename for the VTK file
- `topData::OWENS.TopData`: Structure containing time history data and simulation results
- `setupOutputs::OWENS.SetupOutputs`: Structure containing system setup and assembly data

# Keyword Arguments
- `scaling`: Scaling factor for deflections (default: 1)
- `delta_x`, `delta_y`, `delta_z`: Additional displacement offsets (default: zero vectors)
- `userPointNames`: Names for user-defined points (default: nothing)
- `userPointData`: Data for user-defined points (default: nothing)
- `stress`: Stress data (default: nothing)

# Returns
Calls the original OWENSFEA_VTK function with unpacked data from the input structures.
"""
function OWENSFEA_VTK(
    filename,
    topData,
    setupOutputs;
    scaling = 1,
    delta_x = nothing,
    delta_y = nothing,
    delta_z = nothing,
    userPointNames = nothing,
    userPointData = nothing,
    stress = nothing,
)
    # Extract required data from TopData
    tvec = topData.t
    uHist = topData.uHist
    azi = topData.aziHist
    
    # Extract required data from SetupOutputs
    system = setupOutputs.system
    assembly = setupOutputs.assembly
    sections = setupOutputs.sections
    
    # Set default values for delta vectors if not provided
    if delta_x === nothing
        delta_x = zero(tvec)
    end
    if delta_y === nothing
        delta_y = zero(tvec)
    end
    if delta_z === nothing
        delta_z = zero(tvec)
    end
    
    # Call the original function with unpacked data
    return OWENSFEA_VTK(
        filename,
        tvec,
        uHist,
        system,
        assembly,
        sections;
        scaling = scaling,
        azi = azi,
        delta_x = delta_x,
        delta_y = delta_y,
        delta_z = delta_z,
        userPointNames = userPointNames,
        userPointData = userPointData,
        stress = stress,
    )
end

"""
    mywrite_vtk(name, assembly::Assembly; kwargs...)
    mywrite_vtk(name, assembly::Assembly, state::AssemblyState; kwargs...)
    mywrite_vtk(name, assembly::Assembly, history::Vector{<:AssemblyState}], dt;
        kwargs...)

Write the deformed geometry (and associated data) to a VTK file for visualization
using ParaView.

The `state` argument may be omitted to write the original geometry to a VTK file
without any associated data.

If the solution time `history` is provided, the time step must also be provided

# Keyword Arguments
 - `sections = nothing`: Cross section geometry corresponding to each point,
    defined in a frame aligned with the body frame but centered around the
    corresponding point. Defined as an array with shape `(3, ncross, np)` where `ncross`
    is the number of points in each cross section and `np` is the number of points.
 - `scaling=1.0`: Parameter to scale the deflections (only valid if state is provided)
 - `metadata=Dict()`: Dictionary of metadata for the file(s)
"""
function mywrite_vtk(name, assembly; sections = nothing, metadata = Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i, ip] = assembly.points[ip][i]
            end
        end

        # create cells
        cells = [
            MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for
            i = 1:nelem
        ]

    else

        li = LinearIndices((ncross, npoint))

        # extract cell point locations
        points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
        if ndims(sections) > 2
            for ip = 1:npoint
                for ic = 1:ncross
                    points[:, li[ic, ip]] = assembly.points[ip] + sections[:, ic, ip]
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    points[:, li[ic, ip]] = assembly.points[ip] + sections[:, ic]
                end
            end
        end

        # construct triangle strip for each beam element
        cells = Vector{MeshCell{VTKCellType,Vector{Int64}}}(undef, nelem)
        for ielem = 1:nelem
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ielem]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ielem]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
        end

    end

    # write vtk file
    vtk_grid(name, points, cells) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add local axis data
        axis_name = ["x-axis", "y-axis", "z-axis"]
        axis_vector = [e1, e2, e3]
        for i = 1:3
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                data[:, ielem] .= assembly.elements[ielem].Cab*axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
        end
    end

    return nothing
end

function mywrite_vtk(
    name,
    assembly,
    state;
    sections = nothing,
    scaling = 1.0,
    metadata = Dict(),
)

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i, ip] = assembly.points[ip][i] + scaling*state.points[ip].u[i]
            end
        end

        # create cells
        cells = [
            MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for
            i = 1:nelem
        ]
    else
        li = LinearIndices((ncross, npoint))

        # extract cell point locations
        points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
        if ndims(sections) > 2
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct = wiener_milenkovic(c)'
                    points[:, li[ic, ip]] = assembly.points[ip] + u + Ct*sections[:, ic, ip]
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct = wiener_milenkovic(c)'
                    points[:, li[ic, ip]] = assembly.points[ip] + u + Ct*sections[:, ic]
                end
            end
        end

        # construct triangle strip for each beam element
        cells = Vector{MeshCell{VTKCellType,Vector{Int64}}}(undef, nelem)
        for ielem = 1:nelem
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ielem]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ielem]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
        end
    end

    # write vtk file
    vtk_grid(name, points, cells) do vtkfile

        vtkfile["scaling"] = scaling

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add local axis data
        axis_name = ["x-axis", "y-axis", "z-axis"]
        axis_vector = [e1, e2, e3]
        for i = 1:3
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                CtCab =
                    wiener_milenkovic(
                        state.elements[ielem].theta,
                    )'*assembly.elements[ielem].Cab
                data[:, ielem] .= CtCab * axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
        end

        # add point data
        for field in fieldnames(GXBeam.PointState)
            data = Matrix{eltype(assembly)}(undef, 3, npoint*ncross)
            li = LinearIndices((ncross, npoint))
            for ip = 1:npoint
                data[:, li[:, ip]] .= getproperty(state.points[ip], field)
            end
            vtkfile[string(field), VTKPointData()] = data
        end

        # add cell data
        for field in fieldnames(ElementState)
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                data[:, ielem] .= getproperty(state.elements[ielem], field)
            end
            vtkfile[string(field), VTKCellData()] = data
        end

    end

    return nothing
end

function mywrite_vtk(
    name,
    assembly,
    history,
    t;
    sections = nothing,
    scaling = 1.0,
    metadata = Dict(),
    theta_x = zero(t),
    theta_y = zero(t),
    theta_z = zero(t),
    delta_x = zero(t),
    delta_y = zero(t),
    delta_z = zero(t),
    userPointData = nothing,
    userPointNames = nothing,
    stress = nothing,
) #userPointData[iname,it,ipt] = Float64, userPointNames=["data1","data2","data3"]

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    paraview_collection(name) do pvd

        for (current_step, state) in enumerate(history)

            if isnothing(sections)
                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i, ip] =
                            assembly.points[ip][i] + scaling*state.points[ip].u[i]
                    end
                end

                # create cells
                cells = [
                    MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]])
                    for i = 1:nelem
                ]
            else

                li = LinearIndices((ncross, npoint))

                # extract cell point locations
                points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                if ndims(sections) > 2
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct = GXBeam.wiener_milenkovic(c)'
                            points[:, li[ic, ip]] =
                                assembly.points[ip] + u + Ct*sections[:, ic, ip]
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct = GXBeam.wiener_milenkovic(c)'
                            points[:, li[ic, ip]] =
                                assembly.points[ip] + u + Ct*sections[:, ic]
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType,Vector{Int64}}}(undef, nelem)
                for ielem = 1:nelem
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ielem]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ielem]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
                end
            end

            # Rotate
            cz = cos(theta_z[current_step])
            cy = cos(theta_y[current_step])
            cx = cos(theta_x[current_step])
            sz = sin(theta_z[current_step])
            sy = sin(theta_y[current_step])
            sx = sin(theta_x[current_step])

            Rz = [
                cz -sz 0
                sz cz 0
                0 0 1
            ]

            Ry = [
                cy 0 sy
                0 1 0
                -sy 0 cy
            ]

            Rx = [
                1 0 0
                0 cx -sx
                0 sx cx
            ]

            pointstemp = Rz * Ry * Rx * points
            pointstemp[:, 1] .+= delta_x[current_step]
            pointstemp[:, 2] .+= delta_y[current_step]
            pointstemp[:, 3] .+= delta_z[current_step]

            # write vtk file
            vtkfile = vtk_grid(name*"-step$current_step", pointstemp, cells)

            # add metadata
            vtkfile["scaling"] = scaling
            vtkfile["time"] = t[current_step]
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # add local axis data
            axis_name = ["x-axis", "y-axis", "z-axis"]
            axis_vector = [e1, e2, e3]
            for i = 1:3
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    CtCab =
                        GXBeam.wiener_milenkovic(state.elements[ielem].theta)' *
                        assembly.elements[ielem].Cab
                    data[:, ielem] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
            end

            # add point data
            for field in fieldnames(GXBeam.PointState)
                data = Matrix{eltype(assembly)}(undef, 3, npoint*ncross)
                li = LinearIndices((ncross, npoint))
                for ip = 1:npoint
                    data[:, li[:, ip]] .= getproperty(state.points[ip], field)
                end
                vtkfile[string(field), VTKPointData()] = data
            end


            # add userPointData
            # userPointData[iname,it,ipt] = Float64
            if !isnothing(userPointData)
                for (iname, dataname) in enumerate(userPointNames)
                    data = Matrix{eltype(assembly)}(undef, 1, npoint*ncross)
                    li = LinearIndices((ncross, npoint))
                    if contains(dataname, "Stress")
                        if !isnothing(stress)
                            for ip = 1:npoint
                                data[:, li[:, ip]] = stress[current_step, ip, :]
                            end
                        end
                    else
                        for ip = 1:npoint
                            data[:, li[:, ip]] .= userPointData[iname, current_step, ip]
                        end
                    end
                    vtkfile[dataname, VTKPointData()] = data
                end
            end

            # add cell data
            for field in fieldnames(GXBeam.ElementState)
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    data[:, ielem] .= getproperty(state.elements[ielem], field)
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t[current_step]] = vtkfile
        end
    end

    return nothing
end

"""
    mywrite_vtk(name, assembly::Assembly, [state::AssemblyState, ]λ::Number,
        eigenstate::AssemblyState; scaling=1.0, mode_scaling=1.0, cycles=1,
        steps=100)

Write a series of files corresponding to the elastic motion of the `assembly`
about the deformed state encoded in `state` defined by the eigenvalue `λ` and
the eigenvector encoded in `eigenstate` over the time period specified by `time`.

The steady-state deflections can be scaled with `scaling` and the eigenmode
deflections can be scaled using `mode_scaling`.

The current time is encoded in the metadata tag "time"
"""
function mywrite_vtk(
    name,
    assembly,
    state,
    λ,
    eigenstate;
    sections = nothing,
    scaling = 1.0,
    mode_scaling = 1.0,
    cycles = 1,
    steps = 100,
)

    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    damping = real(λ)
    period = 2*pi/imag(λ)

    if damping <= 0.0
        start = -period*cycles
        stop = 0.0
    else
        start = 0.0
        stop = period*cycles
    end

    time=range(start, stop, length = steps)

    paraview_collection(name) do pvd

        for (current_step, t) in enumerate(time)

            if isnothing(sections)

                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i, ip] =
                            assembly.points[ip][i] +
                            scaling*state.points[ip].u[i] +
                            mode_scaling*real(eigenstate.points[ip].u[i]*exp(λ*t))
                    end
                end

                # create cells
                cells = [
                    MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]])
                    for i = 1:nelem
                ]

            else

                li = LinearIndices((ncross, npoint))

                # extract cell point locations
                points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                if ndims(sections) > 2
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = mode_scaling*real.(eigenstate.points[ip].u .* exp(λ*t))
                            c = mode_scaling*real.(eigenstate.points[ip].theta .* exp(λ*t))
                            Ct = GXBeam.wiener_milenkovic(c)'
                            points[:, li[ic, ip]] =
                                assembly.points[ip] + u + Ct*sections[:, ic, ip]
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = mode_scaling*real.(eigenstate.points[ip].u .* exp(λ*t))
                            c = mode_scaling*real.(eigenstate.points[ip].theta .* exp(λ*t))
                            Ct = GXBeam.wiener_milenkovic(c)'
                            points[:, li[ic, ip]] =
                                assembly.points[ip] + u + Ct*sections[:, ic]
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType,Vector{Int64}}}(undef, nelem)
                for ielem = 1:nelem
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ielem]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ielem]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
                end

            end

            # write vtk file
            vtkfile = vtk_grid(name*"-step$(current_step)", points, cells)

            # add metadata
            vtkfile["scaling"] = scaling
            vtkfile["time"] = t
            vtkfile["phase"] = 2*pi*t/period

            # add local axis data
            axis_name = ["x-axis", "y-axis", "z-axis"]
            axis_vector = [e1, e2, e3]
            for i = 1:3
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    CtCab =
                        GXBeam.wiener_milenkovic(
                            real(state.elements[ielem].theta .* exp(λ*t)),
                        )' * assembly.elements[ielem].Cab
                    data[:, ielem] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
            end

            # add point data
            for field in fieldnames(GXBeam.PointState)
                data = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                li = LinearIndices((ncross, npoint))
                for ip = 1:npoint
                    data[:, li[:, ip]] .=
                        getproperty(state.points[ip], field) +
                        mode_scaling*real(
                            getproperty(eigenstate.points[ip], field) .* exp(λ*t),
                        )
                end
                vtkfile[string(field), VTKPointData()] = data
            end

            # add cell data
            for field in fieldnames(GXBeam.ElementState)
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    data[:, ielem] .=
                        getproperty(state.elements[ielem], field) +
                        mode_scaling*real(
                            getproperty(eigenstate.elements[ielem], field) .* exp(λ*t),
                        )
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t] = vtkfile
        end
    end

    return nothing
end
