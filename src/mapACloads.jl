function mapACloads(u_jLast,udot_j,Omega_j,t,PEy,QCy,NElem,NBlade,RefR,mesh,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,el,turbine3D,env,step_AC,us_param)

N_blade_nodes = length(structuralSpanLocNorm[1,:])+1
# Initialize bladeForce
bladeForce_N = zeros(NBlade,NElem)
bladeForce_T = zeros(NBlade,NElem)
bladeForce_M25 = zeros(NBlade,NElem)

for k = 1:length(turbine3D)
    #TODO: Incorporate deflections and changes in omega ->
    #r,twist,delta,omega all need to be vectors aligning with the ntheta
    #discretizations of the cylinder

    #TODO: ensure that the deflections aren't compounding after a
    #revolution.  They may be wrong

    #TODO: Verify units everywhere

    circ_step_num = floor((step_AC-1)/turbine3D[k].ntheta*turbine3D[k].B)
    circular_step = step_AC-circ_step_num*turbine3D[k].ntheta/turbine3D[k].B
    idx_sub = Int.(collect(circular_step:turbine3D[k].ntheta/turbine3D[k].B:turbine3D[k].ntheta-turbine3D[k].ntheta/turbine3D[k].B+1+circular_step))

    #TODO: this is hard coded for 2 blades, need to simplify
    # Interpolate the deformations onto the aero model for the current step
    norm_disp_h = LinRange(0,1,N_blade_nodes)
    # 1 = Z deformation - not modeled in 2D AC method
    # 2 = Tangential deformation - no real effect on AC model
    offset = 3
    turbine3D[k].r[idx_sub[1]] = turbine3D[k].r[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    turbine3D[k].r[idx_sub[2]] = turbine3D[k].r[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    offset = 4
    turbine3D[k].twist[idx_sub[1]] = turbine3D[k].twist[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    turbine3D[k].twist[idx_sub[2]] = turbine3D[k].twist[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    offset = 5
    turbine3D[k].delta[idx_sub[1]] = turbine3D[k].delta[idx_sub[1]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    turbine3D[k].delta[idx_sub[2]] = turbine3D[k].delta[idx_sub[2]] + FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    # 6 = Sweep deformation, not modeled in AC method - assuming it is small so that it doesn't spill over into the next step/theta discretization

    turbine3D[k].omega[:] .= Omega_j

    # Interpolate deformation induced velocities onto the aero model for the most current step
    offset = 1
    env[k].V_vert[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    env[k].V_vert[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    offset = 2
    env[k].V_tang[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    env[k].V_tang[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    offset = 3
    env[k].V_rad[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    env[k].V_rad[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    offset = 4
    env[k].V_twist[idx_sub[1]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes+offset:6:N_blade_nodes+N_blade_nodes*6-6+offset],k/length(turbine3D))
    env[k].V_twist[idx_sub[2]] = FLOWMath.akima(norm_disp_h,u_jLast[N_blade_nodes*2+offset:6:N_blade_nodes*2+N_blade_nodes*6-6+offset],k/length(turbine3D))
    # 5 = Change in delta angle, not modeled in 2D AC method
    # 6 = Change in sweep angle, not modeled in 2D AC method

    # envin = env[k]
    # turbine_in = turbine3D[k]
    # mat"[$Rp, $Tp, $Zp, $Mp, $envout] = actuatorcylinder_substep($turbine_in, $envin,$us_param, $step_AC,$alpha_rad,$cl_af,$cd_af)"
    # env[k].rho = envout["rho"]
    # env[k].mu = envout["mu"]
    # env[k].G_amp = envout["G_amp"]
    # env[k].gusttime = envout["gusttime"]
    # env[k].gustX0 = envout["gustX0"]
    # env[k].N_Rev = envout["N_Rev"]
    # env[k].idx_sub = envout["idx_sub"]
    # env[k].wsave = envout["wsave"]
    # env[k].Vinf_nominal = envout["Vinf_nominal"]
    # env[k].V_vert = envout["V_vert"]
    # env[k].V_tang = envout["V_tang"]
    # env[k].V_rad = envout["V_rad"]
    # env[k].V_twist = envout["V_twist"]
    # env[k].Vinf = envout["Vinf"]
    # env[k].V_wake_old = envout["V_wake_old"]
    # env[k].steplast = envout["steplast"]
    # turbine3D[k].ntheta = envout["ntheta"]
    # env[k].tau = envout["tau"]

    Q, Rp, Tp, Zp, Vinf_used, alpha, cl, cd, Vloc, Re = VAWTAero.unsteady_step(turbine3D[k],env[k],us_param,step_AC)
    Mp = zeros(length(Zp)) #TODO: fix moment CALCS in VAWTAero



    for j=1:NBlade
        delta = turbine3D[k].delta[idx_sub[j]]
        bladeForce_N[j,k] = -Rp[j]*cos(delta) + -Zp[j]*sin(delta)
        bladeForce_T[j,k] = -Tp[j] #TODO: fix RPI's difficulty to converge when running in reverse (may have to change RPI indexing)
        bladeForce_M25[j,k] = Mp[j]
    end
end



# scatter(t,bladeForce[1].T[floor(blade[j].NElem/2)])
# hold on
# pause[0.001)
#define these from params file
ft2m = 1 / 3.281

#     RefAR = cactusGeom.RefAR*ft2m*ft2m
RefR = RefR*ft2m

spanLocNorm = zeros(NBlade,NElem)
for i=1:NBlade
    spanLocNorm[i,:] = PEy[1:NElem[1,1],1].*RefR[1,1]/(QCy[NElem[1,1]+1,1]*RefR[1,1])
end

#Initialize structuralLoad

    structuralLoad_N = zeros(NBlade,length(structuralElNumbers[1,:]))
    structuralLoad_T = zeros(NBlade,length(structuralElNumbers[1,:]))
    structuralLoad_M25 = zeros(NBlade,length(structuralElNumbers[1,:]))

for i=1:NBlade
    structuralLoad_N[i,:] = FLOWMath.linear(spanLocNorm[i,:],bladeForce_N[i,:],structuralSpanLocNorm[i,:])
    structuralLoad_T[i,:] = FLOWMath.linear(spanLocNorm[i,:],bladeForce_T[i,:],structuralSpanLocNorm[i,:])
    structuralLoad_M25[i,:]= FLOWMath.linear(spanLocNorm[i,:],bladeForce_M25[i,:],structuralSpanLocNorm[i,:])
end

_,numNodesPerBlade = size(structuralNodeNumbers)

#integrate over elements

#read element data in

numDofPerNode = 6
#     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6))
for j = 1:NBlade
    for k = 1:numNodesPerBlade-1
        #get element data
        # orientation angle,xloc,sectionProps,element order]
        elNum = Int(structuralElNumbers[j,k])
        #get dof map
        node1 = Int(structuralNodeNumbers[j,k])
        node2 = Int(structuralNodeNumbers[j,k+1])
        dofList = [(node1-1)*numDofPerNode.+(1:6), (node2-1)*numDofPerNode.+(1:6)]

        elementOrder = 1
        x = [mesh.x[node1], mesh.x[node2]]
        elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
        xloc = [0 elLength]
        twist = el.props[elNum].twist
        sweepAngle = el.psi[elNum]
        coneAngle = el.theta[elNum]
        rollAngle = el.roll[elNum]

        extDistF2Node =  [structuralLoad_T[j,k],   structuralLoad_T[j,k+1]]
        extDistF3Node = -[structuralLoad_N[j,k],   structuralLoad_N[j,k+1]]
        extDistF4Node = -[structuralLoad_M25[j,k], structuralLoad_M25[j,k+1]]

        mat"[$Fe] = calculateLoadVecFromDistForce($elementOrder,$x,$xloc,$twist,$sweepAngle,$coneAngle,$rollAngle,$extDistF2Node,$extDistF3Node,$extDistF4Node)"

        #asssembly
        for m = 1:length(dofList)
            Fg[dofList[m]] =  Fg[dofList[m]].+Fe[m]
        end

    end
end

ForceDof = Float64.(1:length(Fg))

return Fg,ForceDof,env

end
