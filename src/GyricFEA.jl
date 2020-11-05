function initialElementCalculations(model,el,mesh)
    #initialElementCalculations  performs intitial element calculations
    #   [elStorage] = initialElementCalculations(model,el,mesh)
    #
    #   This function performs initial element calculation for use later in
    #   analysis for efficiency gains.
    #
    #      input:
    #      model               = object containing model information
    #      el                  = object containing element information
    #      mesh                = object containing mesh information
    #
    #      output:
    #      elStorage           = object containing stored element data

    #initial element calculation
    numNodesPerEl = 2

    # elStorage = repmat(single_elStorage,mesh.numEl,1)
    elStorage = Array{ElStorage, 1}(undef, mesh.numEl)
    # elStorage = {}

    for i=1:mesh.numEl
        #Calculate Ke and Fe for element i
        elementOrder = model.elementOrder #assign for element i
        modalFlag = true
        xloc = [0.0 el.elLen[i]]
        sectionProps = el.props[i]
        sweepAngle = el.psi[i]
        coneAngle = el.theta[i]
        rollAngle = el.roll[i]
        aeroSweepAngle = 0.0

        elx = zeros(numNodesPerEl)
        ely = zeros(numNodesPerEl)
        elz = zeros(numNodesPerEl)
        for j=1:numNodesPerEl
            #get element cooridnates
            elx[j] = mesh.x[Int(mesh.conn[i,j])]
            ely[j] = mesh.y[Int(mesh.conn[i,j])]
            elz[j] = mesh.z[Int(mesh.conn[i,j])]
        end

        #get concentrated terms associated with elemetn
        massConc,_,_,_,_,_ = ConcMassAssociatedWithElement(mesh.conn[i,:],model.joint,model.nodalTerms.concMass,model.nodalTerms.concStiff,model.nodalTerms.concLoad)

        concMassFlag = !isempty(findall(x->x!=0,massConc))

        Omega = 0.0

        elStorage[i] = calculateTimoshenkoElementInitialRun(elementOrder,modalFlag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,elx,ely,elz,concMassFlag,massConc,Omega) #initial element calculations for storage
    
    end
    return elStorage
end

function calculateTimoshenkoElementInitialRun(elementOrder,modalFlag,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,x,y,z,concMassFlag,concMass,Omega)
    #calculateTimoshenkoElementInitialRun performs initial element calculations
    #   [elStorage] = calculateTimoshenkoElementInitialRun(input)
    #
    #   This function performs initial element calculations and stores them for
    #   later use and efficiency gains.
    #
    #      input:
    #      input      = object containing element input
    #
    #      output:
    #      elStorage  = object containing stored element data
    #
    #-------- assign input block ----------------
    # elementOrder   = input.elementOrder
    # x              = input.x
    # y              = input.y
    # z              = input.z
    # xloc           = input.xloc
    #
    # sectionProps   = input.sectionProps
    # sweepAngle     = input.sweepAngle
    # coneAngle      = input.coneAngle
    # rollAngle      = input.rollAngle
    #
    # concMass = input.concMass
    # concMassFlag = input.concMassFlag

    # CN2H           = eye(3,3)

    #--------------------------------------------

    numGP = 4

    #calculate quad points
    xi,weight = getGP(numGP)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    K11 = zeros(numNodesPerEl,2)
    K12 = zero(K11)
    K13 = zero(K11)
    K14 = zero(K11)
    K15 = zero(K11)
    K16 = zero(K11)
    K22 = zero(K11)
    K23 = zero(K11)
    K24 = zero(K11)
    K25 = zero(K11)
    K26 = zero(K11)
    K33 = zero(K11)
    K34 = zero(K11)
    K35 = zero(K11)
    K36 = zero(K11)
    K44 = zero(K11)
    K45 = zero(K11)
    K46 = zero(K11)
    K55 = zero(K11)
    K56 = zero(K11)
    K66 = zero(K11)

    S11 = zero(K11)
    S12 = zero(K11)
    S13 = zero(K11)
    S14_1 = zero(K11)
    S14_2 = zero(K11)
    S15 = zero(K11)
    S16 = zero(K11)
    S22 = zero(K11)
    S23 = zero(K11)
    S24_1 = zero(K11)
    S24_2 = zero(K11)
    S25 = zero(K11)
    S26 = zero(K11)
    S33 = zero(K11)
    S34_1 = zero(K11)
    S34_2 = zero(K11)
    S35 = zero(K11)
    S36 = zero(K11)
    S44_1 = zero(K11)
    S44_2 = zero(K11)
    S44_3 = zero(K11)
    S45_1 = zero(K11)
    S45_2 = zero(K11)
    S46_1 = zero(K11)
    S46_2 = zero(K11)
    S55 = zero(K11)
    S56 = zero(K11)
    S66 = zero(K11)

    #     F1 = zeros(numNodesPerEl,1)
    #     F3 = F1
    #     F2 = F1
    #     F4 = F1
    #     F5 = F1
    #     F6 = F1


    M11 = zero(K11)
    M15 = zero(K11)
    M16 = zero(K11)
    M22 = zero(K11)
    M24 = zero(K11)
    M33 = zero(K11)
    M34 = zero(K11)
    M44 = zero(K11)
    M55 = zero(K11)
    M56 = zero(K11)
    M66 = zero(K11)

    C12 = zero(K12)
    C13 = zero(K13)
    C14_1 = zero(K14)
    C14_2 = zero(K14)
    C23 = zero(K23)
    C24 = zero(K24)
    C34 = zero(K34)
    C25 = zero(K11)
    C26 = zero(K11)
    C35 = zero(K11)
    C36 = zero(K11)
    C45_1 = zero(K11)
    C45_2 = zero(K11)
    C46_1 = zero(K11)
    C46_2 = zero(K11)

    elementMass = 0.0
    elementItens = zeros(3,3)
    elxm = zeros(3)

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = N
        p_N1_x = p_N_x
        N2 = N
        p_N2_x = p_N_x
        N3 = N
        p_N3_x = p_N_x
        N4 = N
        p_N4_x = p_N_x
        N5 = N
        p_N5_x = p_N_x
        N6 = N
        p_N6_x = p_N_x
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        EA   = interpolateVal(sectionProps.EA,N) #struct stiffness terms
        EIyy = interpolateVal(sectionProps.EIyy,N)
        EIzz = interpolateVal(sectionProps.EIzz,N)
        GJ   = interpolateVal(sectionProps.GJ,N)
        EIyz = interpolateVal(sectionProps.EIyz,N)

        couple16 = 0.0   #int(Ey dA)    v bend - extension
        couple15 = 0.0   #int(Ez dA)    w bend - extension
        couple45 = 0.0   #int(Gz dA)    w bend - twist
        couple46 = 0.0   #int(Gz dA)    v bend - twist
        couple14 = 0.0   # extension twist
        couple34 = couple46
        couple24 = couple45

        rhoA   = interpolateVal(sectionProps.rhoA,N) #struct mass terms
        rhoIyy = interpolateVal(sectionProps.rhoIyy,N)
        rhoIzz = interpolateVal(sectionProps.rhoIzz,N)
        rhoJ   = interpolateVal(sectionProps.rhoJ,N)
        rhoIyz = interpolateVal(sectionProps.rhoIyz,N)

        vprime = 0.0 #set to zero to deactivate nonlinearites from initial element calculations
        wprime = 0.0

        ycm = interpolateVal(sectionProps.ycm,N)
        zcm = interpolateVal(sectionProps.zcm,N)

        xgp      = interpolateVal(x,N1)
        ygp      = interpolateVal(y,N1)
        zgp      = interpolateVal(z,N1)

        #.... end interpolate value at quad points ........

        #adjust moments of inertia for offsets
        rhoIyy = rhoIyy + rhoA*zcm^2
        rhoIzz = rhoIzz + rhoA*ycm^2
        rhoIyz = rhoIyz + rhoA*ycm*zcm
        rhoJ   = rhoJ + rhoA*(ycm^2 + zcm^2)

        #Calculate strutural stiffness sub matrices
        K11 = calculateElement1(EA,integrationFactor,p_N1_x,p_N1_x,K11)
        K12 = calculateElement1(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12)
        K13 = calculateElement1(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13)
        K14 = calculateElement1(couple14,integrationFactor,p_N1_x,p_N4_x,K14)
        K15 = calculateElement1(couple15,integrationFactor,p_N1_x,p_N5_x,K15)
        K16 = calculateElement1(-couple16,integrationFactor,p_N1_x,p_N6_x,K16)
        K22 = calculateElement1(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22)
        K24 = calculateElement1(-couple24,integrationFactor,p_N2_x,p_N4_x,K24)
        K33 = calculateElement1(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33)
        K34 = calculateElement1(couple34,integrationFactor,p_N3_x,p_N4_x,K34)
        K44 = calculateElement1(GJ,integrationFactor,p_N4_x,p_N4_x,K44)
        K45 = calculateElement1(couple45,integrationFactor,p_N4_x,N5,K45)
        K46 = calculateElement1(couple46,integrationFactor,p_N4_x,N6,K46)
        K55 = calculateElement1(EIyy,integrationFactor,p_N5_x,p_N5_x,K55)
        K56 = calculateElement1(-EIyz,integrationFactor,p_N5_x,p_N6_x,K56)
        K66 = calculateElement1(EIzz,integrationFactor,p_N6_x,p_N6_x,K66)

        #Calculate structural mass sub matrices
        M11 = calculateElement1(rhoA,integrationFactor,N1,N1,M11)
        M15 = calculateElement1(rhoA*zcm,integrationFactor,N1,N5,M15)
        M16 = calculateElement1(-rhoA*ycm,integrationFactor,N1,N6,M16)
        M22 = calculateElement1(rhoA,integrationFactor,N2,N2,M22)
        M24 = calculateElement1(-rhoA*zcm,integrationFactor,N2,N4,M24)
        M33 = calculateElement1(rhoA,integrationFactor,N3,N3,M33)
        M34 = calculateElement1(rhoA*ycm,integrationFactor,N3,N4,M34)
        M44 = calculateElement1(rhoJ,integrationFactor,N4,N4,M44)
        M55 = calculateElement1(rhoIyy,integrationFactor,N5,N5,M55)
        M56 = calculateElement1(-rhoIyz,integrationFactor,N5,N6,M56)
        M66 = calculateElement1(rhoIzz,integrationFactor,N6,N6,M66)

        #Calculate Centrifugal load vector and gravity load vector
        #eventually incorporate lambda into gp level to account for variable
        #twist

        O1 = 1 #these are set to unity to get coefficients for omega components
        O2 = 1
        O3 = 1

        posLocal = lambda[1:3,1:3]*[xgp, ygp, zgp]
        xbarlocal = posLocal[1]
        ybarlocal = posLocal[2]
        zbarlocal = posLocal[3]

        #        g=9.81 #gravitational acceleration [m/s^2]
        #        a_x = 0 #acceleration of body in x and y (hardwired to zero for now)
        #        a_y = 0
        #        a_z = -g
        #        fx = rhoA*a_x #let these loads be defined in the inertial frame
        #        fy = rhoA*a_y
        #        fz = rhoA*a_z
        #        rvec = [ 0 ycm zcm]
        #
        #        fi_hub = CN2H*[fxfyfz]
        #
        #        disLoadgpLocal = lambda(1:3,1:3)*fi_hub
        #        disMomentgp = cross(rvec,disLoadgpLocal)

        #        f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal) - disLoadgpLocal(1)    #omega dot loading not
        #        [F1] = calculateVec1(f1,integrationFactor,N1,F1)
        #        f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2) - disLoadgpLocal(2)
        #        [F2] = calculateVec1(f2,integrationFactor,N2,F2)
        #        f3 = rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal) - disLoadgpLocal(3)
        #        [F3] = calculateVec1(f3,integrationFactor,N3,F3)
        #        f4 = rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2))...
        #                   + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3)) - disMomentgp(1)
        #        [F4] = calculateVec1(f4,integrationFactor,N4,F4)
        #        f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3) - disMomentgp(2)
        #        [F5] = calculateVec1(f5,integrationFactor,N5,F5)
        #        f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2))) - disMomentgp(3)
        #        [F6] = calculateVec1(f6,integrationFactor,N6,F6)
        #
        #Gyric matrix (Coriolis)
        C12 = calculateElement1(-2*rhoA*O3,integrationFactor,N1,N2,C12)
        C13 = calculateElement1(2*rhoA*O2,integrationFactor,N1,N3,C13)

        C14_1 = calculateElement1(2*rhoA*(ycm*O2),integrationFactor,N1,N4,C14_1)
        C14_2 = calculateElement1(2*rhoA*(zcm*O3),integrationFactor,N1,N4,C14_2)

        C23 = calculateElement1(-2*rhoA*O1,integrationFactor,N2,N3,C23)
        C24 = calculateElement1(-2*rhoA*ycm*O1,integrationFactor,N2,N4,C24)
        C25 = calculateElement1(2*rhoA*zcm*O3,integrationFactor,N2,N5,C25)
        C26 = calculateElement1(-2*rhoA*ycm*O3,integrationFactor,N2,N6,C26)
        C34 = calculateElement1(-2*rhoA*zcm*O1,integrationFactor,N3,N4,C34)
        C35 = calculateElement1(-2*rhoA*zcm*O2,integrationFactor,N3,N5,C35)
        C36 = calculateElement1(2*rhoA*ycm*O2,integrationFactor,N3,N6,C36)

        C45_1 = calculateElement1(-2*(rhoIyy*O3),integrationFactor,N4,N5,C45_1)
        C45_2 = calculateElement1(-2*(rhoIyz*O2),integrationFactor,N4,N5,C45_2)

        C46_1 = calculateElement1(2*(rhoIzz*O2),integrationFactor,N4,N6,C46_1)
        C46_2 = calculateElement1(2*(rhoIyz*O3),integrationFactor,N4,N6,C46_2)

        #Spin softening matrix
        S11 = calculateElement1(-rhoA*(O2^2+O3^2),integrationFactor,N1,N1,S11)
        S12 = calculateElement1(rhoA*O1*O2,integrationFactor,N1,N2,S12)
        S13 = calculateElement1(rhoA*O1*O3,integrationFactor,N1,N3,S13)

        S14_1 = calculateElement1(rhoA*(ycm*O1*O3),integrationFactor,N1,N4,S14_1)
        S14_2 = calculateElement1(rhoA*(-zcm*O1*O2),integrationFactor,N1,N4,S14_2)

        S15 = calculateElement1(-rhoA*zcm*(O2^2+O3^2),integrationFactor,N1,N5,S15)
        S16 = calculateElement1(rhoA*ycm*(O2^2+O3^2),integrationFactor,N1,N6,S16)
        S22 = calculateElement1(-rhoA*(O1^2+O3^2),integrationFactor,N2,N2,S22)
        S23 = calculateElement1(rhoA*O2*O3,integrationFactor,N2,N3,S23)

        S24_1 = calculateElement1(rhoA*zcm*(O1^2+O3^2),integrationFactor,N2,N4,S24_1)
        S24_2 = calculateElement1(rhoA*ycm*O2*O3,integrationFactor,N2,N4,S24_2)

        S25 = calculateElement1(rhoA*zcm*O1*O2,integrationFactor,N2,N5,S25)
        S26 = calculateElement1(-rhoA*ycm*O1*O2,integrationFactor,N2,N6,S26)
        S33 = calculateElement1(-rhoA*(O1^2+O2^2),integrationFactor,N3,N3,S33)

        S34_1 = calculateElement1(-rhoA*(ycm*(O1^2+O2^2)),integrationFactor,N3,N4,S34_1)
        S34_2 = calculateElement1(-rhoA*(zcm*O2*O3),integrationFactor,N3,N4,S34_2)


        S35 = calculateElement1(rhoA*zcm*O1*O3,integrationFactor,N3,N5,S35)
        S36 = calculateElement1(-rhoA*ycm*O1*O3,integrationFactor,N3,N6,S36)

        S44_1 = calculateElement1(-(rhoIyy*(O1^2+O3^2)),integrationFactor,N4,N4,S44_1)
        S44_2 = calculateElement1(-(rhoIzz*(O1^2+O2^2)),integrationFactor,N4,N4,S44_2)
        S44_3 = calculateElement1(-(2*rhoIyz*O2*O3),integrationFactor,N4,N4,S44_3)

        S45_1 = calculateElement1(rhoIyz*O1*O3,integrationFactor,N4,N5,S45_1)
        S45_2 = calculateElement1(-rhoIyy*O1*O2,integrationFactor,N4,N5,S45_2)

        S46_1 = calculateElement1(rhoIyz*O1*O2,integrationFactor,N4,N6,S46_1)
        S46_2 = calculateElement1(-rhoIzz*O1*O3,integrationFactor,N4,N6,S46_2)


        S55 = calculateElement1(-rhoIyy*(O2^2+O3^2),integrationFactor,N5,N5,S55)
        S56 = calculateElement1(rhoIyz*(O2^2+O3^2),integrationFactor,N5,N6,S56)
        S66 = calculateElement1(-rhoIzz*(O2^2+O3^2),integrationFactor,N6,N6,S66)

        elementMass,elementItens,elxm = calculateElementMass(rhoA,rhoIyy,rhoIzz,rhoIyz,rhoJ,ycm,zcm,xbarlocal,ybarlocal,zbarlocal,integrationFactor,elementMass,elementItens,elxm)

    end

    #######################################
    #Reduced integration loop
    numGP = 1
    xi,weight = getGP(numGP)

    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N5 = N
        N6 = N
        p_N2_x = p_N_x
        p_N3_x = p_N_x
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        EA   = interpolateVal(sectionProps.EA,N) #struct stiffness terms
        GA = EA/2.6*5/6
        #.... end interpolate value at quad points ........

        #Calculate strutural stiffness sub matrices
        K22 = calculateElement1(GA,integrationFactor,p_N2_x,p_N2_x,K22)
        K26 = calculateElement1(-GA,integrationFactor,p_N2_x,N6,K26)
        K33 = calculateElement1(GA,integrationFactor,p_N3_x,p_N3_x,K33)
        K35 = calculateElement1(GA,integrationFactor,p_N3_x,N5,K35)
        K55 = calculateElement1(GA,integrationFactor,N5,N5,K55)
        K66 = calculateElement1(GA,integrationFactor,N6,N6,K66)

    end

    lamSlim = lambda[1:3,1:3]
    lamSlimTran = lamSlim'
    elementMOI = lamSlimTran*elementItens*lamSlim[1:3,1:3]
    elxm = lamSlimTran*elxm
    #
    #
    if (concMassFlag)
        #modify element mass, moi, and xm to account for concentrated terms
        elementMass = elementMass + sum(concMass[1,:])

        elementMOI[1,1] = elementMOI[1,1] + concMass[1,1]*(y[1]^2 + z[1]^2)+ concMass[1,2]*(y[2]^2 + z[2]^2) + concMass[2,1] + concMass[2,2]
        elementMOI[2,2] = elementMOI[2,2] + concMass[1,1]*(x[1]^2 + z[1]^2)+ concMass[1,2]*(x[2]^2 + z[2]^2) + concMass[3,1] + concMass[3,2]
        elementMOI[3,3] = elementMOI[3,3] + concMass[1,1]*(x[1]^2 + y[1]^2)+ concMass[1,2]*(x[2]^2 + y[2]^2) + concMass[4,1] + concMass[4,2]
        elementMOI[1,2] = elementMOI[1,2] - concMass[1,1]*x[1]*y[1] - concMass[1,2]*x[2]*y[2]
        elementMOI[1,3] = elementMOI[1,3] - concMass[1,1]*x[1]*z[1] - concMass[1,2]*x[2]*z[2]
        elementMOI[2,1] = elementMOI[2,1] - concMass[1,1]*x[1]*y[1] - concMass[1,2]*x[2]*y[2]
        elementMOI[2,3] = elementMOI[2,3] - concMass[1,1]*y[1]*z[1] - concMass[1,2]*y[2]*z[2]
        elementMOI[3,1] = elementMOI[3,1] - concMass[1,1]*x[1]*z[1] - concMass[1,2]*x[2]*z[2]
        elementMOI[3,2] = elementMOI[3,2] - concMass[1,1]*y[1]*z[1] - concMass[1,2]*y[2]*z[2]

        elxm[1] = elxm[1] + concMass[1,1]*x[1] + concMass[1,2]*x[2]
        elxm[2] = elxm[2] + concMass[1,1]*y[1] + concMass[1,2]*y[2]
        elxm[3] = elxm[3] + concMass[1,1]*z[1] + concMass[1,2]*z[2]
    end

    #store element mass properties
    return ElStorage(K11,     #Store structural stiffness "K" into elementStorage
          K12,
          K13,
          K14,
          K15,
          K16,
          K22,
          K23,
          K24,
          K25,
          K26,
          K33,
          K34,
          K35,
          K36,
          K44,
          K45,
          K46,
          K55,
          K56,
          K66,
          M11, #Store structural stiffness "M" into elementStorage
          M15,
          M16,
          M22,
          M24,
          M33,
          M34,
          M44,
          M55,
          M56,
          M66,
          0.5*S11, #Store spin softening coefficient "S" into element storage
          S12,
          S13,
          0.5*S15,
          0.5*S16,
          0.5*S22,
          S23,
          S25,
          S26,
          0.5*S33,
          S35,
          S36,
          0.5*S55,
          0.5*S56,
          0.5*S66,
          S14_1,
          S14_2,
          S24_1,
          S24_2,
          S34_1,
          S34_2,
          S45_1,
          S45_2,
          S46_1,
          S46_2,
          S44_1,
          S44_2,
          S44_3,
          C12, #Store coriolis coefficient "C" into element sotrage
          C13,
          C23,
          C24,
          C25,
          C26,
          C34,
          C35,
          C36,
          C14_1,
          C14_2,
          C45_1,
          C45_2,
          C46_1,
          C46_2,
          elementMass,
          elementMOI,
          elxm)
end
