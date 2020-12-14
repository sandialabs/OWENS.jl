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

        #        f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal) - disLoadgpLocal[1]    #omega dot loading not
        #        [F1] = calculateVec1(f1,integrationFactor,N1,F1)
        #        f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2) - disLoadgpLocal[2]
        #        [F2] = calculateVec1(f2,integrationFactor,N2,F2)
        #        f3 = rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal) - disLoadgpLocal[3]
        #        [F3] = calculateVec1(f3,integrationFactor,N3,F3)
        #        f4 = rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2))...
        #                   + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3)) - disMomentgp[1]
        #        [F4] = calculateVec1(f4,integrationFactor,N4,F4)
        #        f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3) - disMomentgp[2]
        #        [F5] = calculateVec1(f5,integrationFactor,N5,F5)
        #        f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2))) - disMomentgp[3]
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

function  structuralDynamicsTransient(model,mesh,el,dispData,Omega,OmegaDot,time,delta_t,elStorage,Fexternal,Fdof,CN2H,rbData)

    #   [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,...
    #                             dispData,Omega,OmegaDot,time,delta_t,...
    #                             elStorage,Fexternal,Fdof,CN2H,rbData)
    #
    #   This function performs transient structural dynamics analysis.
    #
    #   input:
    #   model      = object containing model data
    #   mesh       = object containing mesh data
    #   el         = object containing element data
    #   dispData   = object containing displacement data
    #   Omega      = rotor speed (Hz)
    #   OmegaDot   = rotor acceleratin (Hz)
    #   time       = current simulation time
    #   delta_t    = time step size
    #   elStorage  = object containing stored element data
    #   Fexternal  = vector containing external force values
    #   Fdof       = vector containing global DOF numbering associated with
    #                external force values
    #   CN2H       = transformation matrix from inertial frame to hub frame
    #   rbData     = vector containing rigid body displacement, velocity, and
    #                acceleration
    #
    #   output:
    #   dispOut       = object containing displacement data at end of time step
    #   FReaction_sp1 = vector containing reaction force at turbine base at
    #                   end of time step

    #-------- get model information -----------
    numEl = mesh.numEl
    x = mesh.x
    y = mesh.y
    z = mesh.z
    conn = Int.(mesh.conn) #TODO: make the generator output ints
    numNodes = length(x)
    elementOrder = model.elementOrder
    BC = model.BC

    numNodesPerEl = elementOrder + 1
    numDOFPerNode = 6
    totalNumDOF = numNodes * numDOFPerNode
    # [~,numReducedDOF]=size(model.jointTransform)
    nodalTerms = model.nodalTerms
    nodalTermsCopy = deepcopy(nodalTerms)
    #-----------------------------------------

    #initialize displacements, tolerance, uNorm, iteration count for nonlinear
    #iteration
    unorm = 1e6
    tol = model.nlParams.tolerance
    maxIterations = model.nlParams.maxIterations
    iterationCount = 0

    elx=zeros(numNodesPerEl)
    ely=zeros(numNodesPerEl)
    elz=zeros(numNodesPerEl,2)
    eldisp = zeros(numNodesPerEl*numDOFPerNode)
    eldisp_sm1 = zeros(numNodesPerEl*numDOFPerNode)
    eldispdot = zero(eldisp)
    eldispddot = zero(eldisp)
    eldispiter = zero(eldisp)
    # if (model.nlOn)
    #      iterationType = model.nlParams.iterationType
    # else
    #      iterationType = 'LINEAR'
    # end

    iterationType = "DI"
    analysisType = model.analysisType

    if occursin("TNB",analysisType)
        #------ newmark integration parameters ---------
        alpha = 0.5
        gamma = 0.5
        beta = 0.5*gamma

        delta_t = delta_t
        a1 = alpha*delta_t
        a2 = (1.0-alpha)*delta_t
        a3 = 1.0/(beta*delta_t*delta_t)
        a3 = a3
        a4 = a3*delta_t
        a5 = 1.0/gamma-1.0
        a6 = alpha/(beta*delta_t)
        a7 = alpha/beta - 1.0
        a8 = delta_t*(alpha/gamma-1.0)

        timeInt = TimeInt(delta_t,a1,a2,a3,a4,a5,a6,a7,a8)

        disp_s = copy(dispData.displ_s)
        dispdot_s = copy(dispData.displdot_s)
        dispddot_s = copy(dispData.displddot_s)

        displddot_im1 = copy(dispddot_s)
        displdot_im1 = copy(dispdot_s)
        displ_im1 = copy(disp_s)

    elseif occursin("TD",analysisType)
        #------ dean integration parameters -------------
        alpha = 0.25

        delta_t = delta_t
        a1 = alpha*delta_t^2
        a2 = (1-2*alpha)*delta_t^2
        a3 = delta_t/2.0
        a4 = delta_t*delta_t

        timeInt = TimeInt(delta_t,a1,a2,a3,a4,0.0,0.0,0.0,0.0)

        disp_s = dispData.displ_s
        #     disp_sm1 = dispData.displ_sm1
        #-------------------------------------------
    else
        error("analysis type not supported, choose another")
    end

    #-----------------------------------------------
    # Initialize elInput, and DO NOT redundantly re-assign the memory in the
    # while and for loops below.

    elInput = ElInput(elementOrder,
    true, #modalFlag,
    timeInt,
    zeros(2), #xloc,
    el.props[1], #sectionProps,
    0.0, #sweepAngle,
    0.0, #coneAngle,
    0.0, #rollAngle,
    0.0, #aeroSweepAngle,
    iterationType,
    model.nlOn, #useDisp,
    false, #preStress,
    false, #aeroElasticOn,
    false, #aeroForceOn,
    0.0, #loadStepPrev,
    1.0, #loadStep,
    model.nlParams.maxNumLoadSteps,
    model.nlParams.maxIterations, #MAXIT
    model.nlParams.tolerance, #tolerance
    analysisType,
    zeros(12), #disp,
    eldispdot, #dispdot,
    eldispddot, #dispddot,
    eldispiter, #displ_iter,
    zeros(4,2), #concMass,
    zeros(6,2), #concStiff,
    zeros(6,2), #concLoad,
    eldisp_sm1, #dispm1,
    elx, #x,
    ely, #y,
    elz, #z,
    model.gravityOn,
    model.RayleighAlpha,
    model.RayleighBeta,
    rbData[1:3], #accelVec,
    rbData[4:6], #omegaVec,
    rbData[7:9], #omegaDotVec,
    Omega,
    OmegaDot,
    CN2H,
    model.airDensity,
    0.0, #freq,
    true) #firstIteration

    while (unorm>tol && iterationCount < maxIterations) #iteration loop
        #------- intitialization -----------------
        Kg = zeros(totalNumDOF,totalNumDOF) #initialize global stiffness and force vector
        Fg = zeros(totalNumDOF)

        nodalTerms = deepcopy(nodalTermsCopy)
        #-------------------------------------------

        #---- element  calculation and assembly ----------------------------------
        for i=1:numEl
            #Calculate Ke and Fe for element i
            index = 1                           #initialize element data
            elInput.xloc = [0.0 el.elLen[i]]
            elInput.sectionProps = el.props[i]
            elInput.sweepAngle = el.psi[i]
            elInput.coneAngle = el.theta[i]
            elInput.rollAngle = el.roll[i]
            if (iterationCount == 0)
                elInput.firstIteration = true
            else
                elInput.firstIteration = false
            end

            for j=1:numNodesPerEl

                #get element cooridnates
                elx[j] = x[conn[i,j]]
                ely[j] = y[conn[i,j]]
                elz[j] = z[conn[i,j]]

                #get element nodal displacements at s and s-1 time step
                for k=1:numDOFPerNode
                    #                 if occursin("TD",analysisType)
                    #                     eldisp[index] = disp_s((conn[i,j]-1)*numDOFPerNode + k)
                    #                     eldisp_sm1[index] = disp_sm1((conn[i,j]-1)*numDOFPerNode + k)
                    #                     eldispiter[index] = displ_iter((conn[i,j]-1)*numDOFPerNode + k)
                    #                 end
                    if occursin("TNB",analysisType)
                        eldispiter[index] = displ_im1[(conn[i,j]-1)*numDOFPerNode + k]
                        if (occursin("NR",iterationType))
                            eldisp[index] = displ_im1[(conn[i,j]-1)*numDOFPerNode + k]
                            eldispdot[index] = displdot_im1[(conn[i,j]-1)*numDOFPerNode + k]
                            eldispddot[index] = displddot_im1[(conn[i,j]-1)*numDOFPerNode + k]
                        elseif (occursin("DI",iterationType)||occursin("LINEAR",iterationType))
                            eldisp[index] = disp_s[(conn[i,j]-1)*numDOFPerNode + k]
                            eldispdot[index] = dispdot_s[(conn[i,j]-1)*numDOFPerNode + k]
                            eldispddot[index] = dispddot_s[(conn[i,j]-1)*numDOFPerNode + k]
                        end
                    end
                    index = index + 1
                end
            end

            #get concentrated terms associated with elemetn
            massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff = ConcMassAssociatedWithElement(conn[i,:],model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad)

            elInput.concMass = massConc
            elInput.concStiff = stiffConc
            elInput.concLoad = loadConc
            elInput.disp = eldisp

            # specific to 'TD', but must be declared
            elInput.dispm1= eldisp_sm1

            # specific to 'TNB' , but must be declared
            elInput.dispdot = eldispdot
            elInput.dispddot = eldispddot

            elInput.x = elx
            elInput.y = ely
            elInput.z = elz

            if (el.rotationalEffects)
                elInput.Omega = Omega
                elInput.OmegaDot = OmegaDot
            else
                elInput.Omega = 0.0
                elInput.OmegaDot = 0.0
            end

            elInput.displ_iter = eldispiter

            # Juno.@enter calculateTimoshenkoElementNL(elInput,elStorage[i]) #calculate timoshenko element
            in2 = elStorage[i]
            mat"$elOutput = calculateTimoshenkoElementNL($elInput,$in2)" #calculate timoshenko element

            conin = conn[i,:]


            Kg,Fg = assembly(elOutput["Ke"],elOutput["Fe"],conn[i,:],numNodesPerEl,numDOFPerNode,Kg,Fg) #assemble element stiffness matrix and force vector

            #         Erestotal = Erestotal + elOutput["Eres"]
            #................................................
        end #for
        #------- end element calculation and assembly ------------------

        ##
        #----------------------------------------------------------------------

        ##
        #Apply external loads to structure
        for i=1:length(Fexternal)
            if occursin("TD",analysisType)
                Fg[Fdof[i]] = Fg[Fdof[i]] + Fexternal[i]*delta_t^2
            end
            if occursin("TNB",analysisType)
                Fg[Fdof[i]] = Fg[Fdof[i]] + Fexternal[i]
            end
        end

        #------ apply constraints on system -----------------------------------
        Kg = applyConstraints(Kg,model.jointTransform)
        Fg = applyConstraintsVec(Fg,model.jointTransform)

        #----------------------------------------------------------------------
        ##

        #Apply BCs to global system
        KgTotal,FgTotal = applyBC(Kg,Fg,BC,numDOFPerNode)

        solution = KgTotal\FgTotal  #solve for displacements

        solution = model.jointTransform*solution #transform to full dof listing

        if model.nlOn  #calculate norm between current iteration and last iteration
            if occursin("NR",iterationType)
                unorm = calcUnorm(displ_im1+solution,displ_im1)
            else
                unorm = calcUnorm(solution,displ_im1)
            end
        else
            unorm = 0.0
        end

        if occursin("NR",iterationType)
            #if newton raphson update u, udot, uddot at each iteration
            displ_im1 = displ_im1 + solution
            cap_delta_displ = displ_im1 - dispData.displ_s
            displddot_im1 = timeInt.a3*(cap_delta_displ) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s
            displdot_im1  = -timeInt.a7*dispData.displdot_s -timeInt.a8*dispData.displddot_s + timeInt.a6*(cap_delta_displ)
        elseif (occursin("DI",iterationType)||occursin("LINEAR",iterationType))
            displ_im1 = solution
        else
            error("iteration type not supported, choose another")
        end

        iterationCount = iterationCount + 1
    end #While

    #Calculate reaction at turbine base (hardwired to node number 1)
    reactionNodeNumber = model.platformTurbineConnectionNodeNumber


    mat"$FReaction = calculateReactionForceAtNode($reactionNodeNumber,$model,$mesh,$el,$elStorage,$timeInt,$dispData,$displ_im1,$rbData,$Omega,$OmegaDot,$CN2H)"

    #Calculate strain
    elStrain = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ_im1,model.nlOn)
    if (iterationCount>=maxIterations)
        error("Maximum iterations exceeded.")
    end

    FReaction_sp1 = FReaction
    displ_sp1 = displ_im1

    # Specific to TNB, but must be declared
    displddot_sp1 = timeInt.a3*(displ_sp1-dispData.displ_s) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s #store velocity vector in dispOut
    displdot_sp1 = dispData.displdot_s + timeInt.a2*dispData.displddot_s + timeInt.a1*displddot_sp1                    #store acceleration vector in dispOut

    dispOut = DispOut(elStrain,displ_sp1,displddot_sp1,displdot_sp1)

    return elStrain,dispOut,FReaction_sp1
end


function applyConstraints(Kg,transMatrix)
    #This function transforms a matrix by the transformation matrix to
    #enforce joint constraints
    # Kg = transMatrix'*(Kg*transMatrix)
    return transMatrix'*(Kg*transMatrix)
end

function applyConstraintsVec(Fg,transMatrix)
    #This function transforms a vector by the transformation matrix to
    #enforce joint constraints
    return transMatrix'*Fg
end

function calcUnorm(unew,uold)
    #This function calculates a relative norm between two vectors: unew and
    #uold
    return LinearAlgebra.norm(unew-uold)/LinearAlgebra.norm(unew)
end

function calculateTimoshenkoElementNL(input,elStorage)
    #calculateTimoshenkoElementNL performs nonlinear element calculations
    #   [output] = calculateTimoshenkoElementNL(input,elStorage)
    #
    #   This function performs nonlinear element calculations.
    #
    #      input:
    #      input      = object containing element input
    #      elStorage  = obect containing precalculated element data
    #
    #      output:
    #      output     = object containing element data

    ###-------- assign input block ----------------
    elementOrder   = input.elementOrder
    x              = input.x
    y              = input.y
    z              = input.z
    xloc           = input.xloc
    disp           = input.disp
    sectionProps   = input.sectionProps
    sweepAngle     = input.sweepAngle
    coneAngle      = input.coneAngle
    rollAngle      = input.rollAngle
    Omega          = input.Omega
    OmegaDot       = input.OmegaDot
    concStiff      = input.concStiff
    concMass       = input.concMass
    concLoad       = input.concLoad
    # modalFlag      = input.modalFlag
    omegaVec       = zeros(3)
    omegaDotVec    = zeros(3)
    accelVec       = zeros(3)
    CN2H           = 1.0*LinearAlgebra.I(3) #same as eye(3) #initialize CN2H to identity for static or modal analysis

    useDisp        = input.useDisp
    preStress      = input.preStress
    aeroElasticOn  = input.aeroElasticOn
    aeroForceOn    = input.aeroForceOn
    iterationType  = input.iterationType
    disp_iter      = input.displ_iter
    omegaVec       = input.omegaVec
    CN2H           = input.CN2H
    omegaDotVec    = input.omegaDotVec
    timeInt        = input.timeInt
    analysisType = input.analysisType
    dispm1 = zeros(12) #declare type
    dispdot = zeros(12)#declare type
    dispddot = zeros(12)#declare type

    #options for Dean integrator
    if (occursin("TD",analysisType))
        dispm1         = input.dispm1
    elseif (occursin("TNB",analysisType))#options for newmark beta integrator
        accelVec       = input.accelVec
        dispdot      = input.dispdot
        dispddot     = input.dispddot
    end

    #--------------------------------------------
    #setting for modal analysis flag
    if (occursin("M",analysisType))
        disp_iter=copy(disp)
    end

    #setting for initial reduced order model calculations
    if (occursin("RM0",analysisType))
        disp_iter=copy(disp)
        omegaVec = input.omegaVec
        omegaDotVec = input.omegaDotVec
        accelVec = input.accelVec
    end

    #settings if aeroelastic analysis is active
    if (aeroElasticOn)
        freq = input.freq
    else
        freq = 0.0 #Not used, but must be declared
    end


    numGP = 4   #number of gauss points for full integration
    numGPRI = 1 #number of gauss points for reduced integration
    #calculate quad points
    xi,weight = getGP(numGP)
    xiRI,weightRI = getGP(numGPRI)

    #Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x)

    F1 = zeros(numNodesPerEl,1)
    F3 = zero(F1)
    F2 = zero(F1)
    F4 = zero(F1)
    F5 = zero(F1)
    F6 = zero(F1)

    SS22 = zeros(numNodesPerEl) #initialize pre-stress (stress stiffening matrices)
    SS33 = zero(SS22)

    #initialize nonlinear element matrices, only used if (useDisp)
    K12NL = zeros(numNodesPerEl)
    K13NL = zero(K12NL)
    K22NL = zero(K12NL)
    K23NL = zero(K12NL)
    K33NL = zero(K12NL)
    K22NLhat = zero(K22NL)
    K33NLhat = zero(K33NL)
    K23NLhat = zero(K23NL)


    #initialize aeroelastic matrices only used if aeroElasticOn, but must declare type
    K33Aero = zeros(numNodesPerEl)
    K34Aero = zero(K33Aero)
    C33Aero = zero(K33Aero)
    C34Aero = zero(K33Aero)
    K43Aero = zero(K33Aero)
    K44Aero = zero(K33Aero)
    C43Aero = zero(K33Aero)
    C44Aero = zero(K33Aero)

    C33 = zeros(numNodesPerEl,numNodesPerEl)
    C44 = zero(C33)

    #Convert frequencies from Hz to radians
    Omega = 2*pi*Omega
    OmegaDot = 2*pi*OmegaDot

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node

    twistAvg = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)
    lambdaSlim = lambda[1:3,1:3]

    dispLocal = lambda*disp_iter

    uNode = [dispLocal[1] dispLocal[7]]
    vNode = [dispLocal[2] dispLocal[8]]
    wNode = [dispLocal[3] dispLocal[9]]
    #     theta_xNode = [dispLocal(4)  dispLocal(10)]
    #     theta_yNode = [dispLocal(5)  dispLocal(11)]
    #     theta_zNode = [dispLocal(6)  dispLocal(12)]


    omega_x=omegaVec[1]
    omega_y=omegaVec[2]
    omega_z = omegaVec[3] + Omega
    omegaDot_x=omegaDotVec[1]
    omegaDot_y=omegaDotVec[2]
    omegaDot_z = omegaDotVec[3] + OmegaDot
    Ohub = [omega_x;omega_y;omega_z]
    ODotHub = [omegaDot_x;omegaDot_y;omegaDot_z]
    Oel = lambdaSlim*Ohub
    ODotel = lambdaSlim*ODotHub
    O1 = Oel[1]
    O2 = Oel[2]
    O3 = Oel[3]

    O1dot = ODotel[1]
    O2dot = ODotel[2]
    O3dot = ODotel[3]

    if (input.gravityOn)
        g = 9.81 		  #gravitational acceleration [m/s^2]
    else
        g = 0.0
    end

    a_x = accelVec[1] #acceleration of body in hub frame (from platform rigid body motion)
    a_y = accelVec[2]
    a_z = accelVec[3]

    a_x_n = 0.0 #accelerations in inertial frame
    a_y_n = 0.0
    a_z_n = g
    a_temp = CN2H*[a_x_n; a_y_n; a_z_n]

    a_x = a_x + a_temp[1]
    a_y = a_y + a_temp[2]
    a_z = a_z + a_temp[3]


    #Integration loop
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xi[i],xloc)
        N1 = copy(N);  p_N1_x = copy(p_N_x)
        N2 = copy(N);  p_N2_x = copy(p_N_x)
        N3 = copy(N);  p_N3_x = copy(p_N_x)
        N4 = copy(N);
        N5 = copy(N);
        N6 = copy(N);
        integrationFactor = Jac * weight[i]

        #..... interpolate for value at quad point .....
        rhoA   = interpolateVal(sectionProps.rhoA,N) #struct mass terms

        # Only used if (useDisp || preStress)
        EA   = interpolateVal(sectionProps.EA,N)
        uprime = interpolateVal(uNode,p_N1_x)
        vprime = interpolateVal(vNode,p_N2_x)
        wprime = interpolateVal(wNode,p_N3_x)
        Faxial = EA*(uprime+0.5*vprime^2+0.5*wprime^2)


        #mass center offsets
        ycm = interpolateVal(sectionProps.ycm,N)
        zcm = interpolateVal(sectionProps.zcm,N)

        xgp      = interpolateVal(x,N1)
        ygp      = interpolateVal(y,N1)
        zgp      = interpolateVal(z,N1)



        if (aeroElasticOn || aeroForceOn)
            #aerodynamic props/data
            airDensity = input.airDensity
            acgp     = interpolateVal(sectionProps.ac,N1)
            a0gp     = interpolateVal([sectionProps.a0],N1)
            bgp      = interpolateVal(sectionProps.b,N1)
            agp      = interpolateVal(sectionProps.a,N1)
            radiusgp = sqrt(xgp^2+ygp^2) # radiusgp = 0 for tower
            Uinfgp   = Omega*radiusgp
            twistgp = interpolateVal(sectionProps.twist,N1)
            #.... end interpolate value at quad points ........
        else
            airDensity = 0.0 #Not used but must declare type
            acgp     = 0.0 #Not used but must declare type
            a0gp     = 0.0 #Not used but must declare type
            bgp      = 0.0 #Not used but must declare type
            agp      = 0.0 #Not used but must declare type
            Uinfgp   = 0.0 #Not used but must declare type
            twistgp = 0.0 #Not used but must declare type
        end

        if (aeroElasticOn)
            kgp      = freq*bgp/Uinfgp
            Theogp   = calculateTheo(kgp)
        else
            kgp = 0.0 #Not used but must declare type
            Theogp = 0.0 #Not used but must declare type
        end


        #Calculate Centrifugal load vector and gravity load vector
        #eventually incorporate lambda into gp level to account for variable
        #twist
        posLocal = lambdaSlim*[xgp;ygp;zgp]
        xbarlocal = posLocal[1]
        ybarlocal = posLocal[2]
        zbarlocal = posLocal[3]

        fx = rhoA*a_x #let these loads be defined in the inertial frame
        fy = rhoA*a_y
        fz = rhoA*a_z
        rvec = [ 0; ycm; zcm]

        fi_hub = [fx;fy;fz]

        disLoadgpLocal = lambdaSlim*fi_hub
        cpskew= [0 -rvec[3] rvec[2]
                rvec[3] 0 -rvec[1]
                -rvec[2] rvec[1] 0]
        disMomentgp = cpskew*disLoadgpLocal

        if (preStress) #stress-stiffening/pre-stress calculations
            SS22 = calculateElement1(Faxial,integrationFactor,p_N2_x,p_N2_x,SS22)
            SS33 = calculateElement1(Faxial,integrationFactor,p_N3_x,p_N3_x,SS33)
        end


        #calculate static aerodynamic load
        sectionAeroLift = 0.0
        sectionAeroMoment = 0.0
        if (aeroForceOn)
            cl = a0gp*twistgp*pi/180.0
            qinf = 0.5*airDensity*Uinfgp^2*(2.0*bgp)
            sectionAeroLift = qinf*cl
            sectionAeroMoment = sectionAeroLift*(acgp+agp)
        end

        #distributed/body force load calculations
        f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal + O3dot*ybarlocal - O2dot*zbarlocal) - disLoadgpLocal[1]
        F1 = calculateVec1(f1,integrationFactor,N1,F1)
        f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2 + O1dot*zbarlocal - O3dot*xbarlocal) - disLoadgpLocal[2]
        F2 = calculateVec1(f2,integrationFactor,N2,F2)
        f3 = sectionAeroLift + rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal + O2dot*xbarlocal - O1dot*ybarlocal) - disLoadgpLocal[3]
        F3 = calculateVec1(f3,integrationFactor,N3,F3)
        f4 = sectionAeroMoment + rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2)) + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3) + ycm*(O2dot*xbarlocal - O1dot*ybarlocal) - zcm*(O1dot*zbarlocal - O3dot*xbarlocal)) - disMomentgp[1]
        F4 = calculateVec1(f4,integrationFactor,N4,F4)
        f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3 - O2dot*zbarlocal + O3dot*ybarlocal) - disMomentgp[2]
        F5 = calculateVec1(f5,integrationFactor,N5,F5)
        f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2)) - O3dot*ybarlocal + O2dot*zbarlocal) - disMomentgp[3]
        F6 = calculateVec1(f6,integrationFactor,N6,F6)

        if (aeroElasticOn && (bgp != 0)) #aeroelastic calculations
            #This is a real valued aeroelastic representation from
            #Wright and Cooper
            #This version assumes aerodynamic center at quarter chord of
            #airfoil
            Fgp = real(Theogp)
            Ggp = imag(Theogp)
            lcsrat = a0gp/(2*pi)
            Fgp = Fgp*lcsrat
            Ggp = Ggp*lcsrat
            agp = agp/bgp

            if (Uinfgp==0)
                kgp = 1
            end
            lz = -2*pi*(-0.5*kgp^2-Ggp*kgp) #leading negative for difference in unsteady z and w-flap direction
            lzdot = -2*pi*Fgp #leading negative for difference in unsteady z and w-flap direction
            ltheta = -2*pi*(0.5*kgp^2*agp + Fgp - Ggp*kgp*(0.5-agp)) #leading negative for difference in pitch and torsion angles
            if (kgp==0)
                lthetadot = -2*pi*(.5 + Fgp*(.5-agp))
            else
                lthetadot = -2*pi*(.5 + Fgp*(.5-agp) + Ggp/kgp)
            end

            mz = -2*pi*(-0.5*kgp^2*agp-kgp*(agp+acgp)*Ggp) #same as above for lz,lzdot,leheta,lthetadot
            mzdot = -2*pi*(agp+acgp)*Fgp
            mtheta = -2*pi*(0.5*kgp^2*(1.0/8.0+agp^2)+Fgp*(agp+acgp)-kgp*Ggp*(agp+0.5)*(0.5-agp))
            if (kgp==0)
                mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp))
            else
                mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp)+Ggp/kgp*(agp+acgp))
            end

            k33fac = airDensity*Uinfgp^2*lz
            k34fac = airDensity*Uinfgp^2*bgp*ltheta
            k43fac = -airDensity*Uinfgp^2*bgp*mz #leading negative
            k44fac = -airDensity*Uinfgp^2*bgp^2*mtheta #leading negative

            c33fac = airDensity*Uinfgp*bgp*lzdot
            c34fac = airDensity*Uinfgp*bgp^2*lthetadot
            c43fac = -airDensity*Uinfgp*bgp^2*mzdot #leading negative
            c44fac = -airDensity*Uinfgp*bgp^3*mthetadot #leading negative

            K33Aero = calculateElement1(-k33fac,integrationFactor,N3,N3,K33Aero)
            K34Aero = calculateElement1(-k34fac,integrationFactor,N3,N4,K34Aero)
            C34Aero = calculateElement1(-c34fac,integrationFactor,N3,N4,C34Aero)
            C33Aero = calculateElement1(-c33fac,integrationFactor,N3,N3,C33Aero)

            K43Aero = calculateElement1(-k43fac,integrationFactor,N4,N3,K43Aero)
            K44Aero = calculateElement1(-k44fac,integrationFactor,N4,N4,K44Aero)
            C43Aero = calculateElement1(-c43fac,integrationFactor,N4,N3,C43Aero)
            C44Aero = calculateElement1(-c44fac,integrationFactor,N4,N4,C44Aero)

        end

    end #END OF INTEGRATION LOOP

    #Integration loop
    for i=1:numGPRI
        #Calculate shape functions at quad point i
        N,p_N_x,Jac = calculateShapeFunctions(elementOrder,xiRI[i],xloc)
        p_N1_x = copy(p_N_x)
        p_N2_x = copy(p_N_x)
        p_N3_x = copy(p_N_x)

        integrationFactor = Jac * weightRI[i]

        #..... interpolate for value at quad point .....

        if (useDisp || preStress)
            EA   = interpolateVal(sectionProps.EA,N)
            uprime = interpolateVal(uNode,p_N1_x)
            vprime = interpolateVal(vNode,p_N2_x)
            wprime = interpolateVal(wNode,p_N3_x)
        end


        if (useDisp)
            #nonlinear element calculations
            K12NL = calculateElement1(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12NL)
            K13NL = calculateElement1(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13NL)
            K22NL = calculateElement1(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22NL)
            K23NL = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NL)
            K33NL = calculateElement1(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33NL)

            #K12NLhat = K12
            #K13NLhat = K13
            #nonlinear element tangent matrix component calculations
            # T_ij = K_ij + Khat_ij
            if (occursin("NR",iterationType))
                K22NLhat = calculateElement1(EA*(uprime + vprime^2 + 0.5*wprime^2),integrationFactor,p_N2_x,p_N2_x,K22NLhat)
                K33NLhat = calculateElement1(EA*(uprime + wprime^2 + 0.5*vprime^2),integrationFactor,p_N3_x,p_N3_x,K33NLhat)
                K23NLhat = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NLhat)
            end

        end

    end #END OF REDUCED INTEGRATION LOOP

    #unpack stored element stiffness data
    K11 = elStorage.K11
    K12 = elStorage.K12
    K13 = elStorage.K13
    K14 = elStorage.K14
    K15 = elStorage.K15
    K16 = elStorage.K16
    K22 = elStorage.K22
    K23 = elStorage.K23
    K24 = elStorage.K24
    K25 = elStorage.K25
    K26 = elStorage.K26
    K33 = elStorage.K33
    K34 = elStorage.K34
    K35 = elStorage.K35
    K36 = elStorage.K36
    K44 = elStorage.K44
    K45 = elStorage.K45
    K46 = elStorage.K46
    K55 = elStorage.K55
    K56 = elStorage.K56
    K66 = elStorage.K66

    if (useDisp) #modify stiffness matrices to account for nonlinear effects
        K21 = K12' + 2*K12NL'
        K12 = K12 + K12NL
        K31 = K13' + 2*K13NL'
        K13 = K13 + K13NL
        K22 = K22 + K22NL
        K23 = K23 + K23NL
        K33 = K33 + K33NL
    else
        K21 = K12'
        K31 = K13'
    end

    # Only used if (useDisp)
    K12hat  = K12
    K13hat  = K13

    #unpack stored element mass data
    M11 = elStorage.M11
    M15 = elStorage.M15
    M16 = elStorage.M16
    M22 = elStorage.M22
    M24 = elStorage.M24
    M33 = elStorage.M33
    M34 = elStorage.M34
    M44 = elStorage.M44
    M55 = elStorage.M55
    M56 = elStorage.M56
    M66 = elStorage.M66

    #unpack and scale stored element spin softening data
    S11 = elStorage.S11.*(O2^2 + O3^2)
    S12 = elStorage.S12.*O1*O2
    S13 = elStorage.S13.*O1*O3
    S15 = elStorage.S15.*(O2^2+O3^2)
    S16 = elStorage.S16.*(O2^2+O3^2)
    S22 = elStorage.S22.*(O1^2 + O3^2)
    S23 = elStorage.S23.*O2*O3
    S25 = elStorage.S25.*(O1*O2)
    S26 = elStorage.S26.*(O1*O2)
    S33 = elStorage.S33.*(O1^2+O2^2)
    S35 = elStorage.S35.*O1*O3
    S36 = elStorage.S36.*O1*O3
    S55 = elStorage.S55.*(O2^2+O3^2)
    S56 = elStorage.S56.*(O2^2+O3^2)
    S66 = elStorage.S66.*(O2^2+O3^2)
    S14 = elStorage.S14_1.*O1*O3 + elStorage.S14_2.*O1*O2
    S24 = elStorage.S24_1.*(O1^2+O3^2) + elStorage.S24_2.*O2*O3
    S34 = elStorage.S34_1.*(O1^2+O2^2) + elStorage.S34_2.*O2*O3
    S45 = elStorage.S45_1.*O1*O3 + elStorage.S45_2.*O1*O2
    S46 = elStorage.S46_1.*O1*O2 + elStorage.S46_2.*O1*O3
    S44 = elStorage.S44_1.*(O1^2+O3^2) + elStorage.S44_2.*(O1^2+O2^2) + elStorage.S44_3.*O2*O3

    #unpack and scale stored element Corilois data
    C12 = elStorage.C12.*O3
    C13 = elStorage.C13.*O2
    C23 = elStorage.C23.*O1
    C24 = elStorage.C24.*O1
    C25 = elStorage.C25.*O3
    C26 = elStorage.C26.*O3
    C34 = elStorage.C34.*O1
    C35 = elStorage.C35.*O2
    C36 = elStorage.C36.*O2
    C14 = elStorage.C14_1.*O2 + elStorage.C14_2.*O3
    C45 = elStorage.C45_1.*O3 + elStorage.C45_2.*O2
    C46 = elStorage.C46_1.*O2 + elStorage.C46_2.*O3

    #unpack and scale stored element Circulatory data
    H12 = 0.5*elStorage.C12.*O3dot
    H13 = 0.5*elStorage.C13.*O2dot
    H23 = 0.5*elStorage.C23.*O1dot
    H24 = 0.5*elStorage.C24.*O1dot
    H25 = 0.5*elStorage.C25.*O3dot
    H26 = 0.5*elStorage.C26.*O3dot
    H34 = 0.5*elStorage.C34.*O1dot
    H35 = 0.5*elStorage.C35.*O2dot
    H36 = 0.5*elStorage.C36.*O2dot
    H14 = 0.5*(elStorage.C14_1.*O2dot + elStorage.C14_2.*O3dot)
    H45 = 0.5*(elStorage.C45_1.*O3dot + elStorage.C45_2.*O2dot)
    H46 = 0.5*(elStorage.C46_1.*O2dot + elStorage.C46_2.*O3dot)


    #compile stiffness matrix without rotational effects
    Kenr = mapMatrixNonSym([K11 K12 K13 K14 K15 K16
    K21 K22 K23 K24 K25 K26
    K31 K23' K33 K34 K35 K36
    K13' K24' K34' K44 K45 K46
    K15' K25' K35' K45' K55 K56
    K16' K26' K36' K46' K56' K66])


    #add spin softening and circulatory effects to stiffness marix
    K11 = K11 .+ S11
    K21 = K21 .+ S12' .- H12'
    K12 = K12 .+ S12 .+ H12
    K31 = K31 .+ S13' .- H13'
    K13 = K13 .+ S13 .+ H13
    K41 = K14' .+ S14' .- H14'
    K14 = K14 .+ S14 .+ H14
    K15 = K15 .+ S15
    K16 = K16 .+ S16
    K22 = K22 .+ S22 .+ SS22
    K32 = K23' .+ S23' .- H23'
    K23 = K23 .+ S23 .+ H23
    K42 = K24'.+ S24' .- H24'
    K24 = K24 .+ S24 .+ H24
    K52 = K25'.+ S25' .- H25'
    K25 = K25 .+ S25 .+ H25
    K62 = K26' .+ S26' .- H26'
    K26 = K26 .+ S26 .+ H26
    K33 = K33 .+ S33 .+ SS33
    K43 = K34' .+ S34' .- H34'
    K34 = K34 .+ S34 .+ H34
    K53 = K35' .+ S35' .- H35'
    K35 = K35 .+ S35 .+ H35
    K63 = K36' .+ S36' .- H36'
    K36 = K36 .+ S36 .+ H36
    K44 = K44 .+ S44
    K54 = K45' .+ S45' .- H45'
    K45 = K45 .+ S45 .+ H45
    K64 = K46' .+ S46' .- H46'
    K46 = K46 .+ S46 .+ H46
    K55 = K55 .+ S55
    K56 = K56 .+ S56
    K66 = K66 .+ S66


    C43 = -C34'
    if (aeroElasticOn) #modify element matrices for aeroelastic effects
        K33 = K33 + K33Aero
        K34 = K34 + K34Aero
        C33 = C33 + C33Aero
        C34 = C34 + C34Aero
        K43 = K43 + K43Aero
        K44 = K44 + K44Aero
        C43 = C43 + C43Aero
        C44 = C44 + C44Aero
    end

    #---------------------------------------------
    zm=zeros(2,2)

    #compile stiffness matrix with rotational effects
    Ke = mapMatrixNonSym([K11 K12 K13 K14 K15 K16
    K21 K22 K23 K24 K25 K26
    K31 K32 K33 K34 K35 K36
    K41 K42 K43 K44 K45 K46
    K15' K52 K53 K54 K55 K56
    K16' K62 K63 K64 K56' K66])

    Kehat = 0.0 # Declare type
    if (useDisp && occursin("NR",iterationType))
        #compile component of tangent matrix
        Kehat = mapMatrixNonSym([zm K12hat K13hat zm zm zm
        zm K22NLhat K23NLhat zm zm zm
        zm K23NLhat' K33NLhat zm zm zm
        zm zm zm zm zm zm
        zm zm zm zm zm zm
        zm zm zm zm zm zm])
    end

    #compile Coriolis/damping matrix

    Ce = mapMatrixNonSym([zm C12 C13 C14 zm zm
    -C12' zm C23 C24 C25 C26
    -C13' -C23' C33 C34 C35 C36
    -C14' -C24' C43 C44 C45 C46
    zm -C25' -C35' -C45' zm zm
    zm -C26' -C36' -C46' zm zm])

    #compile mass matrix
    Me = mapMatrixNonSym([M11 zm zm zm M15 M16
    zm M22 zm M24 zm zm
    zm zm M33 M34 zm zm
    zm M24' M34' M44 zm zm
    M15' zm zm zm M55 M56
    M16' zm zm zm M56' M66])

    #account for rayleigh damping
    alpha = input.RayleighAlpha
    beta = input.RayleighBeta

    CeRayleigh = alpha.*Kenr + beta.*Me
    Ce = Ce + CeRayleigh

    #compile element force vector
    Fe = mapVector([F1;F2;F3;F4;F5;F6])

    # transform matrices for sweep
    # Note,a negative sweep angle, will sweep away from the direction of
    # positive rotation
    lambdaTran = lambda'

    # lambda = sparse(lambda) #TODO: see if using sparse arrays actually speeds things up
    # lambdaTran = sparse(lambdaTran)

    Me = lambdaTran*Me*lambda
    Ce = lambdaTran*Ce*lambda

    Ke = lambdaTran*Ke*lambda
    if (useDisp && occursin("NR",iterationType))
        Kehat =  lambdaTran*Kehat*lambda
    end

    Fe = lambdaTran*Fe

    ##

    ##concentrated mass
    #NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
    # if some ycm or zcm offset from the node was accounted for in concentrated mass terms

    concMassFlag = !isempty(findall(x->x!=0,concMass))
    concStiffFlag = !isempty(findall(x->x!=0,concStiff))
    concLoadFlag = !isempty(findall(x->x!=0,concLoad))
    if (concMassFlag)
        #modify Me for concentrated mass
        Me[1,1] = Me[1,1] + concMass[1,1]
        Me[2,2] = Me[2,2] + concMass[1,1]
        Me[3,3] = Me[3,3] + concMass[1,1]
        Me[4,4] = Me[4,4] + concMass[2,1]
        Me[5,5] = Me[5,5] + concMass[3,1]
        Me[6,6] = Me[6,6] + concMass[4,1]

        Me[7,7] = Me[7,7] + concMass[1,2]
        Me[8,8] = Me[8,8] + concMass[1,2]
        Me[9,9] = Me[9,9] + concMass[1,2]
        Me[10,10] = Me[10,10] + concMass[2,2]
        Me[11,11] = Me[11,11] + concMass[3,2]
        Me[12,12] = Me[12,12] + concMass[4,2]

        #modify Ce for concentrated mass
        Ce[1,2] = Ce[1,2] - 2*concMass[1,1]*omega_z
        Ce[2,1] = Ce[2,1] + 2*concMass[1,1]*omega_z
        Ce[1,3] = Ce[1,3] + 2*concMass[1,1]*omega_y
        Ce[3,1] = Ce[3,1] - 2*concMass[1,1]*omega_y
        Ce[2,3] = Ce[2,3] - 2*concMass[1,1]*omega_x
        Ce[3,2] = Ce[3,2] + 2*concMass[1,1]*omega_x
        Ce[7,8] = Ce[7,8] - 2*concMass[1,2]*omega_z
        Ce[8,7] = Ce[8,7] + 2*concMass[1,2]*omega_z
        Ce[7,9] = Ce[7,9] + 2*concMass[1,2]*omega_y
        Ce[9,7] = Ce[9,7] - 2*concMass[1,2]*omega_y
        Ce[8,9] = Ce[8,9] - 2*concMass[1,2]*omega_x
        Ce[9,8] = Ce[9,8] + 2*concMass[1,2]*omega_x
    end

    if (concMassFlag || concStiffFlag)
        #modify Ke for concentrated mass
        Ke[1,1] = Ke[1,1] + concStiff[1,1] - concMass[1,1]*(omega_y^2 + omega_z^2)
        Ke[1,2] = Ke[1,2] + concMass[1,1]*omega_x*omega_y - concMass[1,1]*omegaDot_z
        Ke[2,1] = Ke[2,1] + concMass[1,1]*omega_x*omega_y + concMass[1,1]*omegaDot_z
        Ke[1,3] = Ke[1,3] + concMass[1,1]*omega_x*omega_z + concMass[1,1]*omegaDot_y
        Ke[3,1] = Ke[3,1] + concMass[1,1]*omega_x*omega_z - concMass[1,1]*omegaDot_y
        Ke[2,3] = Ke[2,3] + concMass[1,1]*omega_y*omega_z - concMass[1,1]*omegaDot_x
        Ke[3,2] = Ke[3,2] + concMass[1,1]*omega_y*omega_z + concMass[1,1]*omegaDot_x
        Ke[2,2] = Ke[2,2] + concStiff[2,1] - concMass[1,1]*(omega_x^2 + omega_z^2)
        Ke[3,3] = Ke[3,3] + concStiff[3,1] - concMass[1,1]*(omega_x^2 + omega_y^2)
        Ke[4,4] = Ke[4,4] + concStiff[4,1]
        Ke[5,5] = Ke[5,5] + concStiff[5,1]
        Ke[6,6] = Ke[6,6] + concStiff[6,1]
        Ke[7,7] = Ke[7,7] + concStiff[1,2] - concMass[1,2]*(omega_y^2 + omega_z^2)
        Ke[7,8] = Ke[7,8] + concMass[1,2]*omega_x*omega_y - concMass[1,2]*omegaDot_z
        Ke[8,7] = Ke[8,7] + concMass[1,2]*omega_x*omega_y + concMass[1,2]*omegaDot_z
        Ke[7,9] = Ke[7,9] + concMass[1,2]*omega_x*omega_z + concMass[1,2]*omegaDot_y
        Ke[9,7] = Ke[9,7] + concMass[1,2]*omega_x*omega_z - concMass[1,2]*omegaDot_y
        Ke[8,9] = Ke[8,9] + concMass[1,2]*omega_y*omega_z - concMass[1,2]*omegaDot_x
        Ke[9,8] = Ke[9,8] + concMass[1,2]*omega_y*omega_z + concMass[1,2]*omegaDot_x
        Ke[8,8] = Ke[8,8] + concStiff[2,2] - concMass[1,2]*(omega_x^2 + omega_z^2)
        Ke[9,9] = Ke[9,9] + concStiff[3,2] - concMass[1,2]*(omega_x^2 + omega_y^2)
        Ke[10,10] = Ke[10,10] + concStiff[4,2]
        Ke[11,11] = Ke[11,11] + concStiff[5,2]
        Ke[12,12] = Ke[12,12] + concStiff[6,2]
    end

    #modify Fe for  concentrated load
    if (concMassFlag || concLoadFlag)
        Fe[1] = Fe[1] + concLoad[1,1] + concMass[1,1]*(x[1]*(omega_y^2 + omega_z^2)-omega_x*omega_y*y[1] - omega_x*omega_z*z[1]) + concMass[1,1]*(y[1]*omegaDot_z-z[1]*omegaDot_y)  -  concMass[1,1]*a_x
        Fe[2] = Fe[2] + concLoad[2,1] + concMass[1,1]*(y[1]*(omega_x^2 + omega_z^2)-omega_y*omega_z*z[1] - omega_y*omega_x*x[1]) + concMass[1,1]*(z[1]*omegaDot_x-x[1]*omegaDot_z)  -  concMass[1,1]*a_y
        Fe[3] = Fe[3] + concLoad[3,1] + concMass[1,1]*(z[1]*(omega_x^2 + omega_y^2)-omega_z*omega_x*x[1] - omega_z*omega_y*y[1]) + concMass[1,1]*(x[1]*omegaDot_y-y[1]*omegaDot_x)  -  concMass[1,1]*a_z
        Fe[4] = Fe[4] + concLoad[4,1]
        Fe[5] = Fe[5] + concLoad[5,1]
        Fe[6] = Fe[6] + concLoad[6,1]
        Fe[7] = Fe[7] + concLoad[1,2] + concMass[1,2]*(x[2]*(omega_y^2 + omega_z^2)- omega_x*omega_y*y[2] - omega_x*omega_z*z[2]) + concMass[1,2]*(y[2]*omegaDot_z-z[2]*omegaDot_y) -  concMass[1,2]*a_x
        Fe[8] = Fe[8] + concLoad[2,2] + concMass[1,2]*(y[2]*(omega_x^2 + omega_z^2)-omega_y*omega_z*z[2] - omega_y*omega_x*x[2]) + concMass[1,2]*(z[2]*omegaDot_x-x[2]*omegaDot_z)  -  concMass[1,2]*a_y
        Fe[9] = Fe[9] + concLoad[3,2] + concMass[1,2]*(z[2]*(omega_x^2 + omega_y^2)-omega_z*omega_x*x[2] - omega_z*omega_y*y[2]) + concMass[1,2]*(x[2]*omegaDot_y-y[2]*omegaDot_x)  -  concMass[1,2]*a_z
        Fe[10] = Fe[10] + concLoad[4,2]
        Fe[11] = Fe[11] + concLoad[5,2]
        Fe[12] = Fe[12] + concLoad[6,2]
    end


    ##
    # Declare Types
    Fhate = zeros(12)
    FhatLessConc = zeros(12)
    if (occursin("TD",analysisType)) #calculate effective stiffness matrix and force vector for Dean integrator

        a1 = timeInt.a1
        a2 = timeInt.a2
        a3 = timeInt.a3
        a4 = timeInt.a4

        xn=disp[1:12]
        xnm1=dispm1[1:12]
        A = 2.0.*xn - xnm1
        B = -a1.*xnm1 - a2.*xn
        D = a3.*xnm1

        Khate = Ke*a1 + a3.*Ce + Me
        Fhate = Fe*a4 + Me*(A') + Ke*(B') + Ce*(D')

        FhatLessConc =   Fhate - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]].*a4

        #........................................................

        #..........................................................

        Ke = copy(Khate)
        Fe = copy(Fhate)
    end

    if (occursin("TNB",analysisType)) #calculate effective stiffness matrix and load vector for Newmark-Beta integrator
        #     a1 = timeInt.a1
        #     a2 = timeInt.a2
        a3 = timeInt.a3
        a4 = timeInt.a4
        a5 = timeInt.a5
        a6 = timeInt.a6
        a7 = timeInt.a7
        a8 = timeInt.a8

        u=copy(disp)
        udot=copy(dispdot)
        uddot=copy(dispddot)
        if (occursin("NR",iterationType))    #considerations if newton raphson iteration is used
            if (input.firstIteration)
                A = a3*u + a4*udot + a5*uddot
                B = a6*u + a7*udot + a8*uddot
                Fhate = Fe + Me*(A') + Ce*(B') - Ke*u'
            else
                Fhate = Fe  - Me*uddot' - Ce*udot' - (Ke)*u'
            end
        elseif (occursin("DI",iterationType)||occursin("LINEAR",iterationType))   #considerations if direct iteration is used or linear analysis
            A = a3*u + a4*udot + a5*uddot
            B = a6*u + a7*udot + a8*uddot
            Fhate = Fe + Me*(A) + Ce*(B)
        end

        Khate = Ke + a3.*Me + a6.*Ce
        if (occursin("NR",iterationType)) #considerations if newton raphson iteration is used
            Khate = Kehat + Khate
        end

        FhatLessConc =   Fhate - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]]

        #........................................................

        Ke = copy(Khate)
        Fe = copy(Fhate)

    end

    FhatLessConc  = zeros(size(Fe))

    if (occursin("M",analysisType))
        FhatLessConc =   Fe - [concLoad[1,1]
        concLoad[2,1]
        concLoad[3,1]
        concLoad[4,1]
        concLoad[5,1]
        concLoad[6,1]
        concLoad[1,2]
        concLoad[2,2]
        concLoad[3,2]
        concLoad[4,2]
        concLoad[5,2]
        concLoad[6,2]]

        if (occursin("DI",iterationType))
            Fe = Fe*input.loadStep
        end
    end

    if ((occursin("M",analysisType) || occursin("S",analysisType)) && occursin("NR",iterationType)) #considerations for newton-raphson iteration
        Fe = Fe*input.loadStep - Ke*disp_iter'
        Ke = Ke + Kehat
    end

    ###----- assign output block ----------------

    if !(occursin("M",analysisType)||occursin("RM0",analysisType))
        Me = zeros(size(Me))
        Ce = zeros(size(Ce))
    end

    if !(occursin("TD",analysisType) || occursin("TNB",analysisType))
        FhatLessConc  = zeros(size(Fe))
    end
    ###------------------------------------------
    return ElOutput(FhatLessConc,Ke,Fe,Me,Ce)
end

function mapMatrixNonSym(Ktemp)
    ###----- function to form total stifness matrix and transform to desired
    # DOF mapping

    T = [1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

    #map to FEA numbering
    Kel = T'*Ktemp*T

    #declare map
    # map = [1, 7, 2, 8, 3, 9,...
    #       4, 10, 5, 11, 6, 12];
    #
    # #map to FEA numbering
    # for i=1:a
    #     I=map[i];
    #     for j=1:a
    #         J=map(j);
    #         Kel(I,J) = Ktemp(i,j);
    #     end
    # end

    return Kel

end

function calculateTimoshenkoElementStrain(elementOrder,nlOn,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,disp)
    #calculateTimoshenkoElementNL performs nonlinear element calculations
    #   [output] = calculateTimoshenkoElementNL(input,elStorage)
    #
    #   This function performs nonlinear element calculations.
    #
    #      input:
    #      input      = object containing element input
    #      elStorage  = obect containing precalculated element data
    #
    #      output:
    #      output     = object containing element data

    ###------- assign input block ----------------
    # elementOrder   = input.elementOrder
    # xloc           = input.xloc
    # disp           = input.disp
    # sectionProps   = input.sectionProps
    # sweepAngle     = input.sweepAngle
    # coneAngle      = input.coneAngle
    # rollAngle      = input.rollAngle
    # nlOn = input.nlOn
    ###--------------------------------------------

    numGP = 4   #number of gauss points for full integration
    #calculate quad points
    xi,_ = getGP(numGP)

    p_disp_x = zeros(numGP,6)

    #Initialize element sub matrices and sub vectors

    #Sort displacement vector
    #Written for 2 node element with 6 dof per node
    twistAvg = rollAngle + 0.5*(sectionProps.twist[1] + sectionProps.twist[2])
    lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0)

    dispLocal = lambda*disp'

    uNode = [dispLocal[1] dispLocal[7]]
    vNode = [dispLocal[2] dispLocal[8]]
    wNode = [dispLocal[3] dispLocal[9]]
    theta_xNode = [dispLocal[4]  dispLocal[10]]
    theta_yNode = [dispLocal[5]  dispLocal[11]]
    theta_zNode = [dispLocal[6]  dispLocal[12]]

    #Integration loop
    eps_xx_0 = zeros(numGP)
    eps_xx_z = zeros(numGP)
    eps_xx_y = zeros(numGP)
    gam_xz_0 = zeros(numGP)
    gam_xz_y = zeros(numGP)
    gam_xy_0 = zeros(numGP)
    gam_xy_z = zeros(numGP)
    for i=1:numGP
        #Calculate shape functions at quad point i
        N,p_N_x,_ = calculateShapeFunctions(elementOrder,xi[i],xloc)
        #N1 = N
        #N2 = N
        #N3 = N
        #N4 = N
        N5 = copy(N)
        N6 = copy(N)
        p_N1_x = copy(p_N_x)
        p_N2_x = copy(p_N_x)
        p_N3_x = copy(p_N_x)
        p_N4_x = copy(p_N_x)
        p_N5_x = copy(p_N_x)
        p_N6_x = copy(p_N_x)

        #calculate displacement derivatives at quad point i
        uprime = interpolateVal(uNode,p_N1_x)
        vprime = interpolateVal(vNode,p_N2_x)
        wprime = interpolateVal(wNode,p_N3_x)
        theta_x_prime = interpolateVal(theta_xNode,p_N4_x)
        theta_y_prime = interpolateVal(theta_yNode,p_N5_x)
        theta_y_gp = interpolateVal(theta_yNode,N5)
        theta_z_prime = interpolateVal(theta_zNode,p_N6_x)
        theta_z_gp = interpolateVal(theta_zNode,N6)
        p_disp_x[i,:] = [uprime, vprime, wprime, theta_x_prime, theta_y_prime, theta_z_prime]

        if nlOn
            eps_xx_0[i] = uprime + 0.5*(wprime^2 + vprime^2)
        else
            eps_xx_0[i] = uprime
        end
        eps_xx_z[i] = theta_y_prime
        eps_xx_y[i] = -theta_z_prime
        gam_xz_0[i] =  theta_y_gp + wprime
        gam_xz_y[i] =  theta_x_prime
        gam_xy_0[i] = -theta_z_gp + vprime
        gam_xy_z[i] =  -theta_x_prime
    end #END OF INTEGRATION LOOP

    return ElStrain(eps_xx_0,eps_xx_z,eps_xx_y,gam_xz_0,gam_xz_y,gam_xy_0,gam_xy_z)
end

function elementPostProcess(elementNumber,model,mesh,el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
    #elementPostProcess post processes element for reaction force
    #   [Fpp] = elementPostProcess(elementNumber,model,mesh,el,elStorage,....
    #           timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
    #
    #   This function calculates the reaction force associated with an element.
    #
    #   input:
    #   elementNumber  = node number joint constraints are desired at
    #   model          = object containing model data
    #   mesh           = object containing mesh data
    #   elStorage      = object containing stored element data
    #   el             = object containing element data
    #   timeInt        = object containing time integration parameters
    #   dispData       = object containing displacement data
    #   displ_iter     = converged displacement solution
    #   rbData         = vector containing rigid body displacement, velocity,
    #                     and acceleration
    #   Omega          = rotor speed (Hz)
    #   OmegaDot       = rotor acceleratin (Hz)
    #   CN2H           = transformation matrix from inertial frame to hub frame
    #
    #   output:
    #   Fpp            = vector containing reaction force vector associated
    #                    with element


    #some initializations
    elementOrder = model.elementOrder
    numNodesPerEl = elementOrder + 1
    numDOFPerNode = 6
    elx=zeros(numNodesPerEl)
    ely=zeros(numNodesPerEl)
    elz=zeros(numNodesPerEl,2)

    eldisp = zeros(numNodesPerEl*numDOFPerNode)
    eldisp_sm1 = zeros(numNodesPerEl*numDOFPerNode)
    eldispdot = zero(eldisp)
    eldispddot = zero(eldisp)
    eldispiter = zero(eldisp)
    disp_s = dispData.displ_s
    # not always used, but must be declared, specific to TNB and ROM
    dispdot_s = dispData.displdot_s
    dispddot_s = dispData.displddot_s
    #unpack displacement information
    analysisType = model.analysisType

    #unpack mesh information
    x = mesh.x
    y = mesh.y
    z = mesh.z
    conn = Int.(mesh.conn) #TODO: change this at the source

    # construct elInput for elementNumber
    index = 1
    elementOrder = elementOrder
    modalFlag = true
    timeInt = timeInt
    xloc = [0.0 el.elLen[elementNumber]]
    sectionProps = el.props[elementNumber]
    sweepAngle = el.psi[elementNumber]
    coneAngle = el.theta[elementNumber]
    rollAngle = el.roll[elementNumber]
    aeroSweepAngle = 0.0
    iterationType = "DI"
    useDisp = model.nlOn
    preStress = false
    aeroElasticOn = false
    aeroForceOn = false
    nlParams = model.nlParams

    # if nlParams.adaptiveLoadSteppingFlag
    loadStepPrev = 0.0
    loadStep = 1.0
    # else
    #     loadStepPrev = 0.0
    #     loadStep = nlParams.prescribedLoadStep[1]
    #     fprintf("Prescribed load step: #f\n",nlParams.prescribedLoadStep[1])
    # end

    maxNumLoadSteps = nlParams.maxNumLoadSteps
    MAXIT = nlParams.maxIterations
    tolerance = nlParams.tolerance

    if (occursin("TNB",analysisType) || occursin("ROM",analysisType))
        analysisType = "TNB"
    else#if occursin("S",analysisType)
        analysisType = "M"
        # else
        #     error("analysisType not supported, choose another")
    end

    #unpack connectivity list, nodal terms, etc.
    nodalTerms = model.nodalTerms

    for j=1:numNodesPerEl
        for k=1:numDOFPerNode
            #get element cooridnates
            elx[j] = x[conn[elementNumber,j]]
            ely[j] = y[conn[elementNumber,j]]
            elz[j] = z[conn[elementNumber,j]]
            if (occursin("TNB",analysisType) || occursin("ROM",analysisType))
                eldisp[index] = disp_s[(conn[elementNumber,j]-1)*numDOFPerNode + k]
                eldispdot[index] = dispdot_s[(conn[elementNumber,j]-1)*numDOFPerNode + k]
                eldispddot[index] = dispddot_s[(conn[elementNumber,j]-1)*numDOFPerNode + k]
                eldispiter[index] = displ_iter[(conn[elementNumber,j]-1)*numDOFPerNode + k]
            elseif occursin("S",analysisType)
                eldisp[index] = displ_iter[(conn[elementNumber,j]-1)*numDOFPerNode + k]
            end
            index = index + 1
        end
    end

    disp = copy(eldisp)
    dispdot = copy(eldispdot) #specific to TNB and ROM
    dispddot = copy(eldispddot) #specific to TNB and ROM

    if (occursin("TNB",analysisType) || occursin("ROM",analysisType))
        displ_iter = copy(eldispiter)
    else#if occursin("S",analysisType)
        displ_iter = copy(eldisp)
        eldispiter = copy(eldisp)
        # else
        #     error("analysisType not supported, choose another")
    end

    #get concentrated terms associated with element
    massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff = ConcMassAssociatedWithElement(conn[elementNumber,:],model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad)

    concMass = massConc
    concStiff = stiffConc
    concLoad = loadConc
    disp = eldisp
    dispm1 = eldisp_sm1
    x = elx
    y = ely
    z = elz
    gravityOn = model.gravityOn

    RayleighAlpha = model.RayleighAlpha
    RayleighBeta = model.RayleighBeta

    if (occursin("TNB",analysisType) || occursin("ROM",analysisType))
        accelVec = rbData[1:3]
        omegaVec = rbData[4:6]
        omegaDotVec = rbData[7:9]
    else #if occursin("S",analysisType)
        accelVec = zeros(3)
        omegaVec = zeros(3)
        omegaDotVec = zeros(3)
        # else
        #     error("analysisType not supported, choose another")
    end

    if el.rotationalEffects
        Omega = Omega
        OmegaDot = OmegaDot
    else
        Omega = 0.0
        OmegaDot = 0.0
    end

    # Specific to TNB and ROM, but must be declared
    CN2H = CN2H


    airDensity = model.airDensity
    freq = 0.0 #Is not used for this model type, but must be declared.
    firstIteration = false

    elInput = ElInput(elementOrder,modalFlag,timeInt,xloc,sectionProps,sweepAngle,coneAngle,rollAngle,aeroSweepAngle,iterationType,useDisp,preStress,aeroElasticOn,aeroForceOn,loadStepPrev,loadStep,maxNumLoadSteps,MAXIT,tolerance,analysisType,disp,dispdot,dispddot,displ_iter,concMass,concStiff,concLoad,dispm1,x,y,z,gravityOn,RayleighAlpha,RayleighBeta,accelVec,omegaVec,omegaDotVec,Omega,OmegaDot,CN2H,airDensity,freq,firstIteration)

    #calculate element stiffness matrix and force vector
    #(or effective stiffness matrix and force vector from time integration)
    in2 = elStorage[elementNumber]
    mat"$elOutput = calculateTimoshenkoElementNL($elInput,$in2)"

    #post process for reaction force
    FhatEl1PP = elOutput["Ke"]*eldispiter
    if occursin("TD",analysisType)
        denom = timeInt.a4
    else #if (occursin("TNB",analysisType) || occursin("ROM",analysisType))||occursin("S",analysisType)
        denom = 1.0
        # else
        #     error("analysisType not supported, choose another")
    end
    Fpp = (FhatEl1PP - elOutput["FhatLessConc"])./denom
    ###----------------------------------------------------------------------
    return Fpp
end
