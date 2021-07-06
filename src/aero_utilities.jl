function mapACloads(u_jLast,udot_j,Omega_j,t,PEy,QCy,NElem,NBlade,RefR,mesh,el,turbine3D,env,step_AC,us_param)

    structuralSpanLocNorm = mesh.structuralSpanLocNorm
    structuralNodeNumbers = mesh.structuralNodeNumbers
    structuralElNumbers = mesh.structuralElNumbers

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

        Q, Rp, Tp, Zp, Vinf_used, alpha, cl, cd, Vloc, Re = VAWTAero.Unsteady_Step(turbine3D[k],env[k],us_param,step_AC)
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
            twist_d = el.props[elNum].twist
            sweepAngle_d = el.psi[elNum]
            coneAngle_d = el.theta[elNum]
            rollAngle_d = el.roll[elNum]

            extDistF2Node =  [structuralLoad_T[j,k],   structuralLoad_T[j,k+1]]
            extDistF3Node = -[structuralLoad_N[j,k],   structuralLoad_N[j,k+1]]
            extDistF4Node = -[structuralLoad_M25[j,k], structuralLoad_M25[j,k+1]]

            Fe = calculateLoadVecFromDistForce(elementOrder,x,xloc,twist_d,sweepAngle_d,coneAngle_d,rollAngle_d,extDistF2Node,extDistF3Node,extDistF4Node)

            #asssembly
            for m = 1:length(dofList)
                Fg[dofList[m]] =  Fg[dofList[m]].+Fe[m]
            end

        end
    end

    ForceDof = Float64.(1:length(Fg))

    return Fg,ForceDof,env

end

function mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)

    cactusGeom = readCactusGeom(geomFn)
    blade = cactusGeom.blade

    # aero_data = importCactusFile(loadsFn,1,2002,22,',')
    aero_data = DelimitedFiles.readdlm(loadsFn,',',skipstart = 1)
    #define these from params file
    ft2m = 1 / 3.281
    rho = 1.225
    #     RefAR = cactusGeom.RefAR*ft2m*ft2m
    RefR = cactusGeom.RefR*ft2m
    V = 25 #m/s

    normTime = aero_data[:,1]

    numAeroEl = 0
    for i=1:cactusGeom.NBlade
        numAeroEl = numAeroEl + cactusGeom.blade[i].NElem
    end

    len,_ = size(aero_data)

    numAeroTS = Int(len/numAeroEl)
    # time = zeros(len/numAeroEl,1)
    time = normTime[1:Int(numAeroEl):end,1].*RefR[1]./V[1]

    urel = aero_data[:,12]
    uloc = urel.*V

    cn = aero_data[:,17]
    ct = aero_data[:,18]
    cm25 = aero_data[:,15]

    #     cl = aero_data[:,13]
    #     cd = aero_data[:,14]
    #
    #     cx = aero_data[:,19]
    #     cy = aero_data[:,20]
    #     cz = aero_data[:,21]

    #calculate element areas
    #     Fx = zeros(len)
    #     Fy = Fx
    #     Fz = Fx

    NperSpan = zeros(len)
    TperSpan = zeros(len)
    M25perSpan = zeros(len)
    Mecc = zeros(len)

    for i=1:len

        NperSpan[i] =  cn[i]  * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
        TperSpan[i] =  ct[i]  * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
        M25perSpan[i] = cm25[i] * 0.5*rho*uloc[i]^2*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)*blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR
        momentArm = 0.0
        Mecc[i] = NperSpan[i] * momentArm
    end

    # Initialize bladeForces
    N = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)
    T = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)
    M25 = zeros(cactusGeom.NBlade,numAeroTS,blade[1].NElem)

    index = 1
    for i=1:numAeroTS
        for j=1:cactusGeom.NBlade
            for k=1:blade[j].NElem
                N[j,i,k] = NperSpan[index]
                T[j,i,k] = TperSpan[index]
                M25[j,i,k] = M25perSpan[index]
                index = index + 1
            end
        end
    end

    spanLocNorm = zeros(cactusGeom.NBlade,blade[1].NElem)
    for i=1:cactusGeom.NBlade
        spanLocNorm[i,:] = blade[i].PEy[1:blade[1].NElem[1,1],1].*RefR[1,1]/(blade[i].QCy[blade[1].NElem[1,1]+1,1]*RefR[1,1])
    end

    bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers = readBladeData(bldFn)

    #Initialize structuralLoad
    struct_N = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_T = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_M25 = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))

    for i=1:cactusGeom.NBlade
        for j=1:numAeroTS
            struct_N[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],N[i,j,:],structuralSpanLocNorm[i,:])
            struct_T[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],T[i,j,:],structuralSpanLocNorm[i,:])
            struct_M25[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],M25[i,j,:],structuralSpanLocNorm[i,:])
        end
    end

    _,numNodesPerBlade = size(structuralNodeNumbers)

    #integrate over elements

    #read element aero_data in
    mesh = readMesh(meshFn)
    el = readElementData(mesh.numEl,elFn,ortFn,bladeData)
    numDofPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6),numAeroTS)
    for i=1:numAeroTS
        for j = 1:cactusGeom.NBlade
            for k = 1:numNodesPerBlade-1
                #get element aero_data
                # orientation angle,xloc,sectionProps,element order]
                elNum = Int(structuralElNumbers[j,k])
                #get dof map
                node1 = Int(structuralNodeNumbers[j,k])
                node2 = Int(structuralNodeNumbers[j,k+1])
                dofList = [(node1-1)*numDofPerNode.+(1:6) (node2-1)*numDofPerNode.+(1:6)]

                elementOrder = 1
                x = [mesh.x[node1], mesh.x[node2]]
                elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
                xloc = [0 elLength]
                twist = el.props[elNum].twist
                sweepAngle = el.psi[elNum]
                coneAngle = el.theta[elNum]
                rollAngle = el.roll[elNum]

                extDistF2Node =  [struct_T[j,i,k]    struct_T[j,i,k+1]]
                extDistF3Node = -[struct_N[j,i,k]    struct_N[j,i,k+1]]
                extDistF4Node = -[struct_M25[j,i,k]  struct_M25[j,i,k+1]]

                Fe = calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)

                #asssembly
                for m = 1:length(dofList)
                    Fg[dofList[m],i] =  Fg[dofList[m],i]+Fe[m]
                end

            end
        end
    end

    #reduce Fg to nonzero components
    #assumes any loaded DOF will never be identically zero throughout time
    #history
    ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
    ForceDof = zeros(sum(Fg[:,1].!=0),1)
    index = 1
    for i=1:Int(maximum(maximum(structuralNodeNumbers))*6)
        if !isempty(findall(x->x!=0,Fg[i,:]))

            ForceValHist[index,:] = Fg[i,:]
            ForceDof[index] = i
            index = index + 1
        end
    end

    return time,ForceValHist,ForceDof,cactusGeom
end

function calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)
    #calculateTimoshenkoElementNL performs nonlinear element calculations
    #   [output] = calculateTimoshenkoElementNL(input,elStorage)
    #
    #   This function performs nonlinear element calculations.
    #
    #      input:
    #      input      = object containing element input
    #      elStorage  = obect containing precalculated element aero_data
    #
    #      output:
    #      output     = object containing element aero_data

    #--------------------------------------------
    numGP = 4   #number of gauss points for full integration
    #calculate quad points
    xi,weight = GyricFEA.getGP(numGP)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    F1 = zeros(numNodesPerEl,1)
    F3 = zero(F1)
    F2 = zero(F1)
    F4 = zero(F1)
    F5 = zero(F1)
    F6 = zero(F1)

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg = rollAngle + 0.5*(twist[1] + twist[2])
    lambda = GyricFEA.calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,_,Jac = GyricFEA.calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = N
        N2 = N
        N3 = N
        N4 = N
        N5 = N
        N6 = N
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        extDistF1 = 0
        extDistF2 = GyricFEA.interpolateVal(extDistF2Node,N2)
        extDistF3 = GyricFEA.interpolateVal(extDistF3Node,N3)
        extDistF4 = GyricFEA.interpolateVal(extDistF4Node,N4)
        extDistF5 = 0
        extDistF6 = 0

        #distributed/body force load calculations
        F1 = calculateVec1(extDistF1,integrationFactor,N1,F1)
        F2 = calculateVec1(extDistF2,integrationFactor,N2,F2)
        F3 = calculateVec1(extDistF3,integrationFactor,N3,F3)
        F4 = calculateVec1(extDistF4,integrationFactor,N4,F4)
        F5 = calculateVec1(extDistF5,integrationFactor,N5,F5)
        F6 = calculateVec1(extDistF6,integrationFactor,N6,F6)


    end #END OF INTEGRATION LOOP

    #compile element force vector
    Fe = GyricFEA.mapVector([F1;F2;F3;F4;F5;F6])

    # transform matrices for sweep
    # Note,a negative sweep angle, will sweep away from the direction of
    # positive rotation
    lambdaTran = lambda'
    # lambdaTran = sparse(lambdaTran)
    Fe = lambdaTran*Fe

    return Fe

end

#Element calculation functions
function calculateVec1(f,integrationFactor,N,F)
    #This function is a general routine to calculate an element vector
    len=length(N)
    for i=1:len
        F[i] = F[i] + f*N[i]*integrationFactor
    end
    return F
end


function makeBCdata(pBC,numNodes,numDofPerNode,reducedDOFList,jointTransform)
    #readBDdata  reads boundary condition file
    #   [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
    #
    #   This function reads the boundray condition file and stores data in the
    #   boundary condition object.
    #
    #      input:
    #      bcfilename    = string containing boundary condition filename
    #      numNodes      = number of nodes in structural model
    #      numDofPerNode = number of degrees of freedom per node

    #      output:
    #      BC            = object containing boundary condition data

    totalNumDof = numNodes*numDofPerNode

    numsBC = 0
    nummBC = 0

    #create a vector denoting constrained DOFs in the model (0 unconstrained, 1
    #constrained)

    #calculate constrained dof vector
    isConstrained = zeros(totalNumDof)
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
    numpBC = length(pBC[:,1])

    map = GyricFEA.calculateBCMap(numpBC,pBC,numDofPerNode,reducedDOFList)
    numReducedDof = length(jointTransform[1,:])
    redVectorMap = GyricFEA.constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,numpBC,pBC,isConstrained) #create a map between reduced and full DOF lists

    BC = GyricFEA.BC_struct(numpBC,
    pBC,
    numsBC,
    nummBC,
    isConstrained,
    map,
    redVectorMap)

    return BC

end
