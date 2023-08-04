"""
    mapAD15(t,mesh)

map AD15 forces to OWENS mesh dofs

# Inputs
* `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
* `mesh::GyricFEA.Mesh`: see ?GyricFEA.Mesh
* `mesh::GyricFEA.El`: see ?GyricFEA.El

# Outputs:
* `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
* `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)

"""
function mapAD15(t,azi_j,mesh,advanceAD15;numAeroTS = 1,alwaysrecalc=true,verbosity=0)
    Nturb = length(mesh)
    n_steps,Fx,Fy,Fz,Mx,My,Mz = advanceAD15(t,mesh,azi_j)
    #TODO: multiple turbine Fx, etc.

    # NOTE on AD15 advanceTurb values (Fx,Fy,Fz,Mx,My,Mz)
    #       - forces/moments are in hub coordinates (converted in advanceAD15)
    #       - array length is the number of OWENS mesh points
    #       - This includes the struts (and could include tower when we add that to the AD15 interface)

    
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    ForceValHist = [zeros(Int(mesh[iturb].numNodes*6),numAeroTS) for iturb = 1:Nturb]
    # DOFs are sequential through all nodes
    ForceDof=[collect(1:1:mesh[iturb].numNodes*6) for iturb = 1:Nturb]

    for iturb = 1:Nturb
        # Map loads over from advanceTurb
        for i=1:mesh[iturb].numNodes
            ForceValHist[iturb][(i-1)*6+1,:] = Fx[iturb][i,1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+2,:] = Fy[iturb][i,1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+3,:] = Fz[iturb][i,1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+4,:] = Mx[iturb][i,1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+5,:] = My[iturb][i,1:numAeroTS]
            ForceValHist[iturb][(i-1)*6+6,:] = Mz[iturb][i,1:numAeroTS]
        end
    end
    
    return ForceValHist,ForceDof
end


"""
    mapACDMS(t,mesh,el)

map VAWTAero forces to OWENS mesh dofs

# Inputs
* `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
* `mesh::GyricFEA.Mesh`: see ?GyricFEA.Mesh
* `mesh::GyricFEA.El`: see ?GyricFEA.El

# Outputs:
* `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
* `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)

"""
function mapACDMS(t,azi_j,mesh,el,advanceTurb;numAeroTS = 1,alwaysrecalc=true,outputfile=nothing,offsetmomentarm=0.0)
    CP,Rp,Tp,Zp,alpha,cl,cd_af,Vloc,Re,thetavec,n_steps,Fx_base,Fy_base,Fz_base,
    Mx_base,My_base,Mz_base,power,power2,rev_step,z3Dnorm,delta,Xp,Yp = advanceTurb(t;azi=azi_j+3*pi/2,alwaysrecalc) #add 3pi/2 to align aero with structural azimuth

    NBlade = length(Rp[:,1,1])
    Nslices = length(Rp[1,:,1])

    # Initialize bladeForces
    N = zeros(NBlade,numAeroTS,Nslices)
    T = zeros(NBlade,numAeroTS,Nslices)
    X = zeros(NBlade,numAeroTS,Nslices)
    Y = zeros(NBlade,numAeroTS,Nslices)
    Z = zeros(NBlade,numAeroTS,Nslices)
    M25 = zeros(NBlade,numAeroTS,Nslices)

    for iTS=1:numAeroTS
        if numAeroTS == 1
            t_idx = length(Rp[1,1,:])
        else
            t_idx = iTS
        end
        for jbld=1:NBlade
            for islice=1:Nslices
                N[jbld,iTS,islice] = Rp[jbld,islice,t_idx] #Normal force on the structure is inward, VAWTAero normal is inward positive #we multiply by cos(delta) to go from force per height to force per span, and then divide by cos(delta) to go from radial to normal, so they cancel
                T[jbld,iTS,islice] = -Tp[jbld,islice,t_idx]*cos(delta[jbld,islice]) ##Tangential force on the structure is against turbine rotation, VAWTAero tangential is with rotation positive # multiply by delta to convert from force per height to force per span
                M25[jbld,iTS,islice] = Rp[jbld,islice,t_idx]*offsetmomentarm #0.0#M25perSpan[index]
                X[jbld,iTS,islice] = Xp[jbld,islice,t_idx]#*cos(-azi_j) + Yp[jbld,islice,t_idx]*sin(-azi_j) #*cos(delta[jbld,islice]) 
                Y[jbld,iTS,islice] = Yp[jbld,islice,t_idx]#*sin(-azi_j) + Yp[jbld,islice,t_idx]*cos(-azi_j)#*cos(delta[jbld,islice]) 
                Z[jbld,iTS,islice] = Zp[jbld,islice,t_idx]#*cos(delta[jbld,islice])
            end
        end
    end

    spanLocNorm = zeros(NBlade,Nslices)

    for i=1:NBlade
        spanLocNorm[i,:] = z3Dnorm #Note that the lookup for this is not the span position, but the vertical position
    end

    structuralSpanLocNorm = mesh.structuralSpanLocNorm # this is also just the blade z position
    structuralNodeNumbers = mesh.structuralNodeNumbers
    structuralElNumbers = mesh.structuralElNumbers

    #Initialize structuralLoad
    struct_N = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_T = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_M25 = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_X = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_Y = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_Z = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    if maximum(structuralSpanLocNorm)>1.0 || minimum(structuralSpanLocNorm)<0.0
        @warn "extrapolating on akima spline, unexpected behavior may occur (very large numbers)."
    end
    for i=1:NBlade
        for j=1:numAeroTS
            # AC and DMS calculate inbetween aero slices, so we add the 0 and 1 normed values here to ensure we don't extrapolate
            struct_N[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;N[i,j,:];0.0],structuralSpanLocNorm[i,:])
            struct_T[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;T[i,j,:];0.0],structuralSpanLocNorm[i,:])
            struct_M25[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;M25[i,j,:];0.0],structuralSpanLocNorm[i,:])
            struct_X[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;X[i,j,:];0.0],structuralSpanLocNorm[i,:])
            struct_Y[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;Y[i,j,:];0.0],structuralSpanLocNorm[i,:])
            struct_Z[i,j,:] = FLOWMath.akima([0.0;spanLocNorm[i,:];1.0],[0.0;Z[i,j,:];0.0],structuralSpanLocNorm[i,:])
        end
    end

    _,numNodesPerBlade = size(structuralNodeNumbers)

    #integrate over elements

    #read element aero_data in
    numDofPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(Int(mesh.numNodes*6),numAeroTS)
    Fg_global = zeros(Int(mesh.numNodes*6),numAeroTS)
    for i=1:numAeroTS
        for j = 1:NBlade
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
                elHeight = abs(mesh.z[node2]-mesh.z[node1])
                xloc = [0 elLength]
                twist = el.props[elNum].twist
                sweepAngle = el.psi[elNum]
                coneAngle = el.theta[elNum]
                rollAngle = el.roll[elNum]

                extDistF2Node = [struct_T[j,i,k]    struct_T[j,i,k+1]]
                extDistF3Node = [struct_N[j,i,k]    struct_N[j,i,k+1]]
                extDistF4Node = [struct_M25[j,i,k]  struct_M25[j,i,k+1]]

                Fe = GyricFEA.calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)
                Fe_global = [struct_X[j,i,k]*elHeight,struct_Y[j,i,k]*elHeight,struct_Z[j,i,k]*elHeight,0.0,0.0,0.0,
                            struct_X[j,i,k+1]*elHeight,struct_Y[j,i,k+1]*elHeight,struct_Z[j,i,k+1]*elHeight,0.0,0.0,0.0]

                #assembly
                for m = 1:length(dofList)
                    Fg[dofList[m],i] =  Fg[dofList[m],i]+Fe[m]
                    Fg_global[dofList[m],i] =  Fg_global[dofList[m],i]+Fe_global[m]
                end

            end
        end
    end

    #reduce Fg to nonzero components
    #assumes any loaded DOF will never be identically zero throughout time
    #history
    # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
    # ForceDof = zeros(sum(Fg[:,1].!=0),1)
    ForceValHist = zeros(length(Fg[:,1]),length(Fg[1,:]))
    ForceDof = zeros(Int,length(Fg[:,1]),1)
    index = 1
    for i=1:Int(mesh.numNodes*6)
        # if !isempty(findall(x->x!=0,Fg[i,:]))

            ForceValHist[index,:] = Fg[i,:]
            ForceDof[index] = i
            index = index + 1
        # end
    end

    if outputfile!=nothing
        DelimitedFiles.open(string("$(outputfile)_fullmesh.txt"), "a") do io
            if t==0
                header1 = ["t" "azi" "nodenum" "Fx" "Fy" "Fz" "Mx" "My" "Mz"]
                header2 = ["(s)" "(rad)" "(#)" "(m/s)" "(N)" "(N)" "(N)" "(N-m)" "(N-m)" "(N-m)"]
            
                DelimitedFiles.writedlm(io, header1, '\t')
                DelimitedFiles.writedlm(io, header2, '\t')
            end
            Fx = ForceValHist[1:6:end,end]
            Fy = ForceValHist[2:6:end,end]
            Fz = ForceValHist[3:6:end,end]
            Mx = ForceValHist[4:6:end,end]
            My = ForceValHist[5:6:end,end]
            Mz = ForceValHist[6:6:end,end]
            for inode = 1:Int(mesh.numNodes)

                data = [t azi_j inode Fx[inode] Fy[inode] Fz[inode] Mx[inode] My[inode] Mz[inode]]
        
                DelimitedFiles.writedlm(io, data, '\t')
        
            end
        
        end

        for j = 1:NBlade
            DelimitedFiles.open(string("$(outputfile)_blade$j.txt"), "a") do io
                if t==0
                    header1 = ["t" "azi" "nodenum" "Fx" "Fy" "Fz" "Mx" "My" "Mz"]
                    header2 = ["(s)" "(rad)" "(#)" "(m/s)" "(N)" "(N)" "(N)" "(N-m)" "(N-m)" "(N-m)"]
                
                    DelimitedFiles.writedlm(io, header1, '\t')
                    DelimitedFiles.writedlm(io, header2, '\t')
                end
                Fx = ForceValHist[1:6:end,end]
                Fy = ForceValHist[2:6:end,end]
                Fz = ForceValHist[3:6:end,end]
                Mx = ForceValHist[4:6:end,end]
                My = ForceValHist[5:6:end,end]
                Mz = ForceValHist[6:6:end,end]
                for k = 1:numNodesPerBlade
                    #get element aero_data
                    # orientation angle,xloc,sectionProps,element order]
                    elNum = Int(structuralElNumbers[j,k])
                    #get dof map
                    inode = Int(structuralNodeNumbers[j,k])
    
                    data = [t azi_j inode Fx[inode] Fy[inode] Fz[inode] Mx[inode] My[inode] Mz[inode]]
            
                    DelimitedFiles.writedlm(io, data, '\t')
            
                end
            
            end
        end
    end

    # return Fexternal, Fdof
    return ForceValHist[:,1:numAeroTS],ForceDof,Fg_global,ForceDof,ForceValHist[:,1:numAeroTS],z3Dnorm
end

# """
#     mapACDMS(t,mesh,el,loadsFn)
#
# map VAWTAero forces to OWENS mesh dofs using a file of loads TODO: merge the two functions together
#
# # Inputs
# * `t::float`: time at which to get the loads (can be called repeatedly at the same time or for large time gaps, will infill run as needed)
# * `mesh::GyricFEA.Mesh`: see ?GyricFEA.Mesh
# * `mesh::GyricFEA.El`: see ?GyricFEA.El
# * `loadsFn::string`: path/name to loads filename, in cactus Element_Data format
#
#
# # Outputs:
# * `ForceValHist::Array(<:float)`: Force or moment (N, N-m) at the time corresponding to the time specified
# * `ForceDof::Array(<:int)`: DOF numbers cooresponding to forces (i.e. mesh element 1 has dofs 1-6, 2 has dofs 7-12, etc)
#
# """
# function mapCACTUSFILE_minimalio(t,mesh,el,loadsFn)
#     # TODO: if this function is used, either export the required data from VAWTAero to get rid of these globals, or put the function back into VAWTAero scope
#     global turbslices
#     global envslices
#     NBlade = turbslices[1].B
#
#     aero_data = DelimitedFiles.readdlm(loadsFn,',',skipstart = 1)
#
#     #define these from params file
#     rho = envslices[1].rho
#
#     RefR = turbslices[1].R
#     NElem = length(turbslices)
#     #TODO: get from VAWTAero structs
#     # chord = RefR.*[1.11093e-01,1.11093e-01,8.20390e-02,6.35940e-02,5.68145e-02,5.55467e-02,5.88008e-02,6.82860e-02,9.11773e-02,1.11093e-01,1.11093e-01]
#     # chord = (chord[1:end-1]+chord[2:end])./2
#     chord = [turb.chord for turb in turbslices]
#     V = envslices[1].V_x[1] #m/s #TODO: get nominal vinf
#     global z3Dnorm
#     # z3Dnorm = 1/2.44680.*[1.22340e-01,3.67020e-01,6.11700e-01,8.56380e-01,1.10106e+00,1.34574e+00,1.59042e+00,1.83510e+00,2.07978e+00,2.32446e+00]
#
#     normTime = aero_data[:,1]
#
#     numAeroEl = 0
#     for i=1:NBlade
#         numAeroEl = numAeroEl + NElem
#     end
#
#     len,_ = size(aero_data)
#
#     numAeroTS = Int(len/numAeroEl)
#
#     time = normTime[1:Int(numAeroEl):end,1].*RefR[1]./V[1]
#
#     urel = aero_data[:,15]
#     uloc = urel.*V
#
#     cn = aero_data[:,24]
#     ct = aero_data[:,25]
#     cm25 = aero_data[:,22]
#
#     NperSpan = zeros(len)
#     TperSpan = zeros(len)
#     M25perSpan = zeros(len)
#
#     for i=1:len
#         NperSpan[i] =  cn[i]  * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
#         TperSpan[i] =  ct[i]  * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)
#         M25perSpan[i] = cm25[i] * 0.5*rho*uloc[i]^2#*(blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR)*blade[Int(aero_data[i,3])].ECtoR[Int(aero_data[i,4])]*RefR
#     end
#
#     # Initialize bladeForces
#     N = zeros(NBlade,numAeroTS,NElem)
#     T = zeros(NBlade,numAeroTS,NElem)
#     M25 = zeros(NBlade,numAeroTS,NElem)
#
#     index = 1
#     for i=1:numAeroTS
#         for j=1:NBlade
#             for k=1:NElem
#                 N[j,i,k] = NperSpan[index]
#                 T[j,i,k] = TperSpan[index]
#                 M25[j,i,k] = M25perSpan[index]
#                 index = index + 1
#             end
#         end
#     end
#
#     #Apply chord since it was pulled out above
#     for i=1:numAeroTS
#         for j=1:NBlade
#             N[j,i,:] = N[j,i,:].*chord
#             T[j,i,:] = T[j,i,:].*chord
#             M25[j,i,:] = M25[j,i,:].*chord
#         end
#     end
#
#     spanLocNorm = zeros(NBlade,NElem)
#
#     for i=1:NBlade
#         spanLocNorm[i,:] = z3Dnorm
#     end
#
#     structuralSpanLocNorm = mesh.structuralSpanLocNorm
#     structuralNodeNumbers = mesh.structuralNodeNumbers
#     structuralElNumbers = mesh.structuralElNumbers
#
#     #Initialize structuralLoad
#     struct_N = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#     struct_T = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#     struct_M25 = zeros(NBlade,numAeroTS,length(structuralElNumbers[1,:]))
#
#     for i=1:NBlade
#         for j=1:numAeroTS
#             struct_N[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],N[i,j,:],structuralSpanLocNorm[i,:])
#             struct_T[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],T[i,j,:],structuralSpanLocNorm[i,:])
#             struct_M25[i,j,:] = FLOWMath.linear(spanLocNorm[i,:],M25[i,j,:],structuralSpanLocNorm[i,:])
#         end
#     end
#
#     _,numNodesPerBlade = size(structuralNodeNumbers)
#
#     #integrate over elements
#
#     #read element aero_data in
#     numDofPerNode = 6
#     #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
#     Fg = zeros(Int(max(maximum(structuralNodeNumbers))*6),numAeroTS)
#     for i=1:numAeroTS
#         for j = 1:NBlade
#             for k = 1:numNodesPerBlade-1
#                 #get element aero_data
#                 # orientation angle,xloc,sectionProps,element order]
#                 elNum = Int(structuralElNumbers[j,k])
#                 #get dof map
#                 node1 = Int(structuralNodeNumbers[j,k])
#                 node2 = Int(structuralNodeNumbers[j,k+1])
#                 dofList = [(node1-1)*numDofPerNode.+(1:6) (node2-1)*numDofPerNode.+(1:6)]
#
#                 elementOrder = 1
#                 x = [mesh.x[node1], mesh.x[node2]]
#                 elLength = sqrt((mesh.x[node2]-mesh.x[node1])^2 + (mesh.y[node2]-mesh.y[node1])^2 + (mesh.z[node2]-mesh.z[node1])^2)
#                 xloc = [0 elLength]
#                 twist = el.props[elNum].twist
#                 sweepAngle = el.psi[elNum]
#                 coneAngle = el.theta[elNum]
#                 rollAngle = el.roll[elNum]
#
#                 extDistF2Node =  [struct_T[j,i,k]    struct_T[j,i,k+1]]
#                 extDistF3Node = -[struct_N[j,i,k]    struct_N[j,i,k+1]]
#                 extDistF4Node = -[struct_M25[j,i,k]  struct_M25[j,i,k+1]]
#
#                 Fe = GyricFEA.calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)
#
#                 #asssembly
#                 for m = 1:length(dofList)
#                     Fg[dofList[m],i] =  Fg[dofList[m],i]+Fe[m]
#                 end
#
#             end
#         end
#     end
#
#     #reduce Fg to nonzero components
#     #assumes any loaded DOF will never be identically zero throughout time
#     #history
#     # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
#     # ForceDof = zeros(sum(Fg[:,1].!=0),1)
#     ForceValHist = zeros(length(Fg[:,1]),length(Fg[1,:]))
#     ForceDof = zeros(length(Fg[:,1]),1)
#     index = 1
#     for i=1:Int(maximum(maximum(structuralNodeNumbers))*6)
#         # if !isempty(findall(x->x!=0,Fg[i,:]))
#
#             ForceValHist[index,:] = Fg[i,:]
#             ForceDof[index] = i
#             index = index + 1
#         # end
#     end
#
#     #TODO: wrap the function at this level for time so you don't read in the file each time
#     Fexternal = zeros(length(ForceDof))
#     for i = 1:length(ForceDof)
#         Fexternal[i] = FLOWMath.linear(time,ForceValHist[i,:],t)
#     end
#
#     return Fexternal, ForceDof
# end
