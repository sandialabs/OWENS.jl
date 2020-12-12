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
    println("Here")
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
    ForceDof = zeros(Int,sum(Fg[:,1].!=0),1)
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
    xi,weight = getGP(numGP)

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
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,_,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = N
        N2 = N
        N3 = N
        N4 = N
        N5 = N
        N6 = N
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        extDistF1 = 0
        extDistF2 = interpolateVal(extDistF2Node,N2)
        extDistF3 = interpolateVal(extDistF3Node,N3)
        extDistF4 = interpolateVal(extDistF4Node,N4)
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
    Fe = mapVector([F1;F2;F3;F4;F5;F6])

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

function mapVector(Ftemp)
    # function to form total force vector and transform to desired
    # DOF mapping
    a=length(Ftemp)
    Fel=zeros(a)

    # #declare map
    map = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]

    for i=1:a
        I=map[i]
        Fel[I] = Ftemp[i]
    end
    return Fel
end
