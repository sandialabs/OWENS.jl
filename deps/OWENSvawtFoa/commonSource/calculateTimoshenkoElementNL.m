function [output] = calculateTimoshenkoElementNL(input,elStorage)
%calculateTimoshenkoElementNL performs nonlinear element calculations
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [output] = calculateTimoshenkoElementNL(input,elStorage)
%                    
%   This function performs nonlinear element calculations.
%
%      input:
%      input      = object containing element input   
%      elStorage  = obect containing precalculated element data
% 
%      output:
%      output     = object containing element data

%-------- assign input block ----------------
elementOrder   = input.elementOrder;
x              = input.x;
y              = input.y;
z              = input.z;
xloc           = input.xloc;
disp           = input.disp;
sectionProps   = input.sectionProps;
sweepAngle     = input.sweepAngle;
coneAngle      = input.coneAngle;
rollAngle      = input.rollAngle;
Omega          = input.Omega;
OmegaDot       = input.OmegaDot;
concStiff      = input.concStiff;
concMass       = input.concMass;
concLoad       = input.concLoad;
modalFlag      = input.modalFlag;
omegaVec       = zeros(3,1);
omegaDotVec    = zeros(3,1);
accelVec       = zeros(3,1);
CN2H           = eye(3); %initialize CN2H to identity for static or modal analysis

useDisp        = input.useDisp;
preStress      = input.preStress;
aeroElasticOn  = input.aeroElasticOn;
aeroForceOn    = input.aeroForceOn;
iterationType  = input.iterationType;

%options for Dean integrator
if(strcmp(input.analysisType,'TD'))
    timeInt        = input.timeInt;
    CN2H           = input.CN2H;
    omegaVec       = input.omegaVec;
    omegaDotVec    = input.omegaDotVec;
    dispm1         = input.dispm1;
    disp_iter      = input.displ_iter;
end

%options for newmark beta integrator
if(strcmp(input.analysisType,'TNB'))
    timeInt        = input.timeInt;
    CN2H           = input.CN2H;
    omegaVec       = input.omegaVec;
    omegaDotVec    = input.omegaDotVec;
	accelVec       = input.accelVec;
    dispdot      = input.dispdot;
    dispddot     = input.dispddot;
    disp_iter    = input.displ_iter;
end

%--------------------------------------------
%setting for modal analysis flag
if(strcmp(input.analysisType,'M'))
    disp_iter=disp;
end

%setting for initial reduced order model calculations
if(strcmp(input.analysisType,'RM0'))
    disp_iter=disp;
    omegaVec = input.omegaVec;
    omegaDotVec = input.omegaDotVec;
    accelVec = input.accelVec;
end

%settings if aeroelastic analysis is active
if(aeroElasticOn)
   freq = input.freq;
end


numGP = 4;   %number of gauss points for full integration
numGPRI = 1; %number of gauss points for reduced integration 
%calculate quad points
    [xi,weight] = getGP(numGP);
    [xiRI,weightRI] = getGP(numGPRI);
    
%Initialize element sub matrices and sub vectors
    numNodesPerEl = length(x);
 
    F1 = zeros(numNodesPerEl,1);
    F3 = F1;
    F2 = F1;
    F4 = F1;
    F5 = F1;
    F6 = F1;
    
    SS22 = zeros(numNodesPerEl); %initialize pre-stress (stress stiffening matrices)
    SS33 = SS22;

    if(useDisp) %initialize nonlinear element matrices
        K12NL = zeros(numNodesPerEl);
        K13NL = K12NL;
        K22NL = K12NL;
        K23NL = K12NL;
        K33NL = K12NL;
        K22NLhat = K22NL;
        K33NLhat = K33NL;
        K23NLhat = K23NL;
    end
    
    if(aeroElasticOn) %initialize aeroelastic matrices
        K33Aero = zeros(numNodesPerEl);
        K34Aero = K33Aero;
        C33Aero = K33Aero;
        C34Aero = K33Aero;
        K43Aero = K33Aero;
        K44Aero = K33Aero;
        C43Aero = K33Aero;
        C44Aero = K33Aero;
    end
            C33 = zeros(numNodesPerEl);
            C44 = C33;
    
%Convert frequencies from Hz to radians
    Omega = 2*pi*Omega;
    OmegaDot = 2*pi*OmegaDot;

%Sort displacement vector
%Written for 2 node element with 6 dof per node
twistAvg = rollAngle + 0.5*(sectionProps.twist(1) + sectionProps.twist(2));
lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0);
lambdaSlim = lambda(1:3,1:3);

dispLocal = lambda*disp_iter';       

    uNode = [dispLocal(1) dispLocal(7)];
    vNode = [dispLocal(2) dispLocal(8)];
    wNode = [dispLocal(3) dispLocal(9)];
    theta_xNode = [dispLocal(4)  dispLocal(10)];
    theta_yNode = [dispLocal(5)  dispLocal(11)];
    theta_zNode = [dispLocal(6)  dispLocal(12)];


       omega_x=omegaVec(1); omega_y=omegaVec(2); omega_z = omegaVec(3) + Omega;
       omegaDot_x=omegaDotVec(1); omegaDot_y=omegaDotVec(2); omegaDot_z = omegaDotVec(3) + OmegaDot;
       Ohub = [omega_x;omega_y;omega_z];
       ODotHub = [omegaDot_x;omegaDot_y;omegaDot_z];
       Oel = lambdaSlim*Ohub;
       ODotel = lambdaSlim*ODotHub;
       O1 = Oel(1);
       O2 = Oel(2);
       O3 = Oel(3);
       
       O1dot = ODotel(1);
       O2dot = ODotel(2);
       O3dot = ODotel(3);
       
       if(input.gravityOn)
           g = 9.81; 		  %gravitational acceleration [m/s^2]
       else
           g = 0.0;
       end
       
       a_x = accelVec(1); %acceleration of body in hub frame (from platform rigid body motion)
       a_y = accelVec(2);
       a_z = accelVec(3);
	   
	   a_x_n = 0.0; %accelerations in inertial frame
	   a_y_n = 0.0;
	   a_z_n = g; 
	   a_temp = CN2H*[a_x_n; a_y_n; a_z_n];
	   
	   a_x = a_x + a_temp(1);
	   a_y = a_y + a_temp(2);
	   a_z = a_z + a_temp(3);
	   

%Integration loop
    for i=1:numGP
        %Calculate shape functions at quad point i
       [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi(i),xloc);
       N1 = N;  p_N1_x = p_N_x;
       N2 = N;  p_N2_x = p_N_x;
       N3 = N;  p_N3_x = p_N_x;
       N4 = N; 
       N5 = N;
       N6 = N;
       integrationFactor = Jac * weight(i);

       %..... interpolate for value at quad point .....
       rhoA   = interpolateVal(sectionProps.rhoA,N); %struct mass terms
       
       if(useDisp || preStress)
            EA   = interpolateVal(sectionProps.EA,N);
            uprime = interpolateVal(uNode,p_N1_x);
            vprime = interpolateVal(vNode,p_N2_x);
            wprime = interpolateVal(wNode,p_N3_x);
            Faxial = EA*(uprime+0.5*vprime^2+0.5*wprime^2);
       end
       
       %mass center offsets
       ycm = interpolateVal(sectionProps.ycm,N);
       zcm = interpolateVal(sectionProps.zcm,N);
       
       xgp      = interpolateVal(x,N1);
       ygp      = interpolateVal(y,N1);
       zgp      = interpolateVal(z,N1);
       
       if(aeroElasticOn || aeroForceOn)
       %aerodynamic props/data
       airDensity = input.airDensity;
       acgp     = interpolateVal(sectionProps.ac,N1);
       a0gp     = interpolateVal([sectionProps.a0],N1);
       bgp      = interpolateVal(sectionProps.b,N1);
       agp      = interpolateVal(sectionProps.a,N1);
       radiusgp = sqrt(xgp^2+ygp^2); % radiusgp = 0 for tower
       Uinfgp   = Omega*radiusgp;
       twistgp = interpolateVal(sectionProps.twist,N1);
       %.... end interpolate value at quad points ........
       end
       
       if(aeroElasticOn)
          kgp      = freq*bgp/Uinfgp;
          Theogp   = calculateTheo(kgp);
       end
           
     
       %Calculate Centrifugal load vector and gravity load vector
       %eventually incorporate lambda into gp level to account for variable
       %twist
       posLocal = lambdaSlim*[xgp;ygp;zgp];
       xbarlocal = posLocal(1);
       ybarlocal = posLocal(2);
       zbarlocal = posLocal(3);
       
       fx = rhoA*a_x; %let these loads be defined in the inertial frame
       fy = rhoA*a_y;
       fz = rhoA*a_z;
       rvec = [ 0; ycm; zcm];
       
      fi_hub = [fx;fy;fz];
       
       disLoadgpLocal = lambdaSlim*fi_hub;
       cpskew= [0,-rvec(3),rvec(2);rvec(3),0,-rvec(1);-rvec(2),rvec(1),0];
       disMomentgp = cpskew*disLoadgpLocal;
       
       if(preStress) %stress-stiffening/pre-stress calculations
           [SS22] = calculateElement1(Faxial,integrationFactor,p_N2_x,p_N2_x,SS22);
           [SS33] = calculateElement1(Faxial,integrationFactor,p_N3_x,p_N3_x,SS33);
       end
       
       
       %calculate static aerodynamic load
       sectionAeroLift = 0.0;
       sectionAeroMoment = 0.0;
       if(aeroForceOn)
          cl = a0gp*twistgp*pi/180.0; 
          qinf = 0.5*airDensity*Uinfgp^2*(2.0*bgp);
          sectionAeroLift = qinf*cl;
          sectionAeroMoment = sectionAeroLift*(acgp+agp);
       end
       
       %distributed/body force load calculations
       f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal + O3dot*ybarlocal - O2dot*zbarlocal) - disLoadgpLocal(1);
       [F1] = calculateVec1(f1,integrationFactor,N1,F1);
       f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2 + O1dot*zbarlocal - O3dot*xbarlocal) - disLoadgpLocal(2);
       [F2] = calculateVec1(f2,integrationFactor,N2,F2);
       f3 = sectionAeroLift + rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal + O2dot*xbarlocal - O1dot*ybarlocal) - disLoadgpLocal(3);
       [F3] = calculateVec1(f3,integrationFactor,N3,F3);
       f4 = sectionAeroMoment + rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2))...
                  + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3) + ycm*(O2dot*xbarlocal - O1dot*ybarlocal) - zcm*(O1dot*zbarlocal - O3dot*xbarlocal)) - disMomentgp(1);
       [F4] = calculateVec1(f4,integrationFactor,N4,F4);
       f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3 - O2dot*zbarlocal + O3dot*ybarlocal) - disMomentgp(2);
       [F5] = calculateVec1(f5,integrationFactor,N5,F5); 
       f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2)) - O3dot*ybarlocal + O2dot*zbarlocal) - disMomentgp(3);
       [F6] = calculateVec1(f6,integrationFactor,N6,F6); 
    
   if(aeroElasticOn && (bgp ~= 0)) %aeroelastic calculations
          %This is a real valued aeroelastic representation from
          %Wright and Cooper
          %This version assumes aerodynamic center at quarter chord of
          %airfoil
        Fgp = real(Theogp);
        Ggp = imag(Theogp);
        lcsrat = a0gp/(2*pi);
        Fgp = Fgp*lcsrat;
        Ggp = Ggp*lcsrat;
        agp = agp/bgp;

    if(Uinfgp==0)
        kgp = 1;
    end
    lz = -2*pi*(-0.5*kgp^2-Ggp*kgp); %leading negative for difference in unsteady z and w-flap direction
    lzdot = -2*pi*Fgp; %leading negative for difference in unsteady z and w-flap direction
    ltheta = -2*pi*(0.5*kgp^2*agp + Fgp - Ggp*kgp*(0.5-agp)); %leading negative for difference in pitch and torsion angles
    if(kgp==0)
        lthetadot = -2*pi*(.5 + Fgp*(.5-agp));
    else
        lthetadot = -2*pi*(.5 + Fgp*(.5-agp) + Ggp/kgp);
    end

    mz = -2*pi*(-0.5*kgp^2*agp-kgp*(agp+acgp)*Ggp); %same as above for lz,lzdot,leheta,lthetadot
    mzdot = -2*pi*(agp+acgp)*Fgp;
    mtheta = -2*pi*(0.5*kgp^2*(1.0/8.0+agp^2)+Fgp*(agp+acgp)-kgp*Ggp*(agp+0.5)*(0.5-agp));
    if(kgp==0)
        mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp));
    else
        mthetadot = -2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+acgp)*(.5-agp)+Ggp/kgp*(agp+acgp));
    end

     k33fac = airDensity*Uinfgp^2*lz;
     k34fac = airDensity*Uinfgp^2*bgp*ltheta;
     k43fac = -airDensity*Uinfgp^2*bgp*mz; %leading negative
     k44fac = -airDensity*Uinfgp^2*bgp^2*mtheta; %leading negative

     c33fac = airDensity*Uinfgp*bgp*lzdot;
     c34fac = airDensity*Uinfgp*bgp^2*lthetadot;
     c43fac = -airDensity*Uinfgp*bgp^2*mzdot; %leading negative
     c44fac = -airDensity*Uinfgp*bgp^3*mthetadot; %leading negative

     [K33Aero] = calculateElement1(-k33fac,integrationFactor,N3,N3,K33Aero);
     [K34Aero] = calculateElement1(-k34fac,integrationFactor,N3,N4,K34Aero);
     [C34Aero] = calculateElement1(-c34fac,integrationFactor,N3,N4,C34Aero);
     [C33Aero] = calculateElement1(-c33fac,integrationFactor,N3,N3,C33Aero);

     [K43Aero] = calculateElement1(-k43fac,integrationFactor,N4,N3,K43Aero);
     [K44Aero] = calculateElement1(-k44fac,integrationFactor,N4,N4,K44Aero);
     [C43Aero] = calculateElement1(-c43fac,integrationFactor,N4,N3,C43Aero);
     [C44Aero] = calculateElement1(-c44fac,integrationFactor,N4,N4,C44Aero);

    end   
       
    end %END OF INTEGRATION LOOP
    
%Integration loop
    for i=1:numGPRI
        %Calculate shape functions at quad point i
       [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xiRI(i),xloc);
       p_N1_x = p_N_x;
       p_N2_x = p_N_x;
       p_N3_x = p_N_x;

       integrationFactor = Jac * weightRI(i);

       %..... interpolate for value at quad point .....
   
       if(useDisp || preStress)
            EA   = interpolateVal(sectionProps.EA,N);
            uprime = interpolateVal(uNode,p_N1_x);
            vprime = interpolateVal(vNode,p_N2_x);
            wprime = interpolateVal(wNode,p_N3_x);
       end
       
           
       if(useDisp)
            %nonlinear element calculations
            [K12NL] = calculateElement1(0.5*EA*vprime,integrationFactor,p_N1_x,p_N2_x,K12NL);
            [K13NL] = calculateElement1(0.5*EA*wprime,integrationFactor,p_N1_x,p_N3_x,K13NL);
            [K22NL] = calculateElement1(0.5*EA*vprime^2,integrationFactor,p_N2_x,p_N2_x,K22NL);
            [K23NL] = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NL);
            [K33NL] = calculateElement1(0.5*EA*wprime^2,integrationFactor,p_N3_x,p_N3_x,K33NL);
            
            %K12NLhat = K12;
            %K13NLhat = K13;
            %nonlinear element tangent matrix component calculations
            % T_ij = K_ij + Khat_ij
            if(strcmp(iterationType,'NR'))
                [K22NLhat] = calculateElement1(EA*(uprime + vprime^2 + 0.5*wprime^2),integrationFactor,p_N2_x,p_N2_x,K22NLhat);
                [K33NLhat] = calculateElement1(EA*(uprime + wprime^2 + 0.5*vprime^2),integrationFactor,p_N3_x,p_N3_x,K33NLhat);
                [K23NLhat] = calculateElement1(0.5*EA*vprime*wprime,integrationFactor,p_N2_x,p_N3_x,K23NLhat);
            end
            
       end
       
   end %END OF REDUCED INTEGRATION LOOP
    
    %unpack stored element stiffness data
    K11 = elStorage.K11;
    K12 = elStorage.K12;
    K13 = elStorage.K13;
    K14 = elStorage.K14;
    K15 = elStorage.K15;
    K16 = elStorage.K16;
    K22 = elStorage.K22;
    K23 = elStorage.K23;
    K24 = elStorage.K24;
    K25 = elStorage.K25;
    K26 = elStorage.K26;
    K33 = elStorage.K33;
    K34 = elStorage.K34;
    K35 = elStorage.K35;
    K36 = elStorage.K36;
    K44 = elStorage.K44;
    K45 = elStorage.K45;
    K46 = elStorage.K46;
    K55 = elStorage.K55;
    K56 = elStorage.K56;
    K66 = elStorage.K66;
    
    if(useDisp) %modify stiffness matrices to account for nonlinear effects
       K21 = K12' + 2*K12NL';
       K12 = K12 + K12NL; 
       K31 = K13' + 2*K13NL';
       K13 = K13 + K13NL;
       K22 = K22 + K22NL;
       K23 = K23 + K23NL;
       K33 = K33 + K33NL;
	   K12hat  = K12;
	   K13hat  = K13;
    else
        K21 = K12';
        K31 = K13';
    end
    
    %unpack stored element mass data
    M11 = elStorage.M11;
    M15 = elStorage.M15;
    M16 = elStorage.M16;
    M22 = elStorage.M22;
    M24 = elStorage.M24;
    M33 = elStorage.M33;
    M34 = elStorage.M34;
    M44 = elStorage.M44;
    M55 = elStorage.M55;
    M56 = elStorage.M56;
    M66 = elStorage.M66;
    
    %unpack and scale stored element spin softening data
    S11 = elStorage.S11.*(O2^2 + O3^2);
    S12 = elStorage.S12.*O1*O2;
    S13 = elStorage.S13.*O1*O3;
    S15 = elStorage.S15.*(O2^2+O3^2);
    S16 = elStorage.S16.*(O2^2+O3^2);
    S22 = elStorage.S22.*(O1^2 + O3^2);
    S23 = elStorage.S23.*O2*O3;
    S25 = elStorage.S25.*(O1*O2);
    S26 = elStorage.S26.*(O1*O2);
    S33 = elStorage.S33.*(O1^2+O2^2);
    S35 = elStorage.S35.*O1*O3;
    S36 = elStorage.S36.*O1*O3;
    S55 = elStorage.S55.*(O2^2+O3^2);
    S56 = elStorage.S56.*(O2^2+O3^2);
    S66 = elStorage.S66.*(O2^2+O3^2);
    S14 = elStorage.S14_1.*O1*O3 + elStorage.S14_2.*O1*O2;
    S24 = elStorage.S24_1.*(O1^2+O3^2) + elStorage.S24_2.*O2*O3;
    S34 = elStorage.S34_1.*(O1^2+O2^2) + elStorage.S34_2.*O2*O3;
    S45 = elStorage.S45_1.*O1*O3 + elStorage.S45_2.*O1*O2;
    S46 = elStorage.S46_1.*O1*O2 + elStorage.S46_2.*O1*O3;
    S44 = elStorage.S44_1.*(O1^2+O3^2) + elStorage.S44_2.*(O1^2+O2^2) + elStorage.S44_3.*O2*O3;
       
    %unpack and scale stored element Corilois data
    C12 = elStorage.C12.*O3;
    C13 = elStorage.C13.*O2;
    C23 = elStorage.C23.*O1;
    C24 = elStorage.C24.*O1;
    C25 = elStorage.C25.*O3;
    C26 = elStorage.C26.*O3;
    C34 = elStorage.C34.*O1;
    C35 = elStorage.C35.*O2;
    C36 = elStorage.C36.*O2;
    C14 = elStorage.C14_1.*O2 + elStorage.C14_2.*O3;
    C45 = elStorage.C45_1.*O3 + elStorage.C45_2.*O2;
    C46 = elStorage.C46_1.*O2 + elStorage.C46_2.*O3;
    
    %unpack and scale stored element Circulatory data
    H12 = 0.5*elStorage.C12.*O3dot;
    H13 = 0.5*elStorage.C13.*O2dot;
    H23 = 0.5*elStorage.C23.*O1dot;
    H24 = 0.5*elStorage.C24.*O1dot;
    H25 = 0.5*elStorage.C25.*O3dot;
    H26 = 0.5*elStorage.C26.*O3dot;
    H34 = 0.5*elStorage.C34.*O1dot;
    H35 = 0.5*elStorage.C35.*O2dot;
    H36 = 0.5*elStorage.C36.*O2dot;
    H14 = 0.5*(elStorage.C14_1.*O2dot + elStorage.C14_2.*O3dot);
    H45 = 0.5*(elStorage.C45_1.*O3dot + elStorage.C45_2.*O2dot);
    H46 = 0.5*(elStorage.C46_1.*O2dot + elStorage.C46_2.*O3dot);
    
    
    %compile stiffness matrix without rotational effects
    [Kenr] = mapMatrixNonSym([K11,K12,K13,K14,K15,K16;
                           K21,K22,K23,K24,K25,K26;
                           K31,K23',K33,K34,K35,K36;
                           K13',K24',K34',K44,K45,K46;
                           K15',K25',K35',K45',K55,K56;
                           K16',K26',K36',K46',K56',K66]);   

                       
    %add spin softening and circulatory effects to stiffness marix   
    K11 = K11 + S11;
    K21 = K21 + S12' - H12';
    K12 = K12 + S12 + H12;
    K31 = K31 + S13' - H13';
    K13 = K13 + S13 + H13;
    K41 = K14' + S14' - H14';
    K14 = K14 + S14 + H14;
    K15 = K15 + S15;
    K16 = K16 + S16;
    K22 = K22 + S22 + SS22;
    K32 = K23' +S23' - H23';
    K23 = K23 + S23 + H23;
    K42 = K24'+ S24' - H24';
    K24 = K24 + S24 + H24;
    K52 = K25'+ S25' - H25';
    K25 = K25 + S25 + H25;
    K62 = K26' + S26' - H26';
    K26 = K26 + S26 + H26;
    K33 = K33 + S33 + SS33;
    K43 = K34' + S34' - H34';
    K34 = K34 + S34 + H34;
    K53 = K35' + S35' - H35';
    K35 = K35 + S35 + H35;
    K63 = K36' + S36' - H36';
    K36 = K36 + S36 + H36;
    K44 = K44 + S44;
    K54 = K45' + S45' - H45';
    K45 = K45 + S45 + H45;
    K64 = K46' + S46' - H46';
    K46 = K46 + S46 + H46;
    K55 = K55 + S55;
    K56 = K56 + S56;
    K66 = K66 + S66;
    
    
    C43 = -C34';
   if(aeroElasticOn) %modify element matrices for aeroelastic effects
        K33 = K33 + K33Aero;
        K34 = K34 + K34Aero;
        C33 = C33 + C33Aero;
        C34 = C34 + C34Aero;
        K43 = K43 + K43Aero;
        K44 = K44 + K44Aero;
        C43 = C43 + C43Aero;
        C44 = C44 + C44Aero;
   end
    
    %---------------------------------------------
    zm=zeros(2,2);
    
    %compile stiffness matrix with rotational effects
    [Ke] = mapMatrixNonSym([K11,K12,K13,K14,K15,K16;
                           K21,K22,K23,K24,K25,K26;
                           K31,K32,K33,K34,K35,K36;
                           K41,K42,K43,K44,K45,K46;
                           K15',K52,K53,K54,K55,K56;
                           K16',K62,K63,K64,K56',K66]);
    
                       
    if(useDisp && strcmp(iterationType,'NR'))
    %compile component of tangent matrix
    [Kehat] = mapMatrixNonSym([zm,K12hat,K13hat,zm,zm,zm;
                               zm,K22NLhat,K23NLhat,zm,zm,zm;
                               zm,K23NLhat',K33NLhat,zm,zm,zm;
                               zm,zm,zm,zm,zm,zm;
                               zm,zm,zm,zm,zm,zm;
                               zm,zm,zm,zm,zm,zm;]);
    end
                       
    %compile Coriolis/damping matrix
    [Ce] = mapMatrixNonSym([zm,C12,C13,C14,zm,zm;
                           -C12',zm,C23,C24,C25,C26;
                           -C13',-C23',C33,C34,C35,C36;
                           -C14',-C24',C43,C44,C45,C46;
                           zm,-C25',-C35',-C45',zm,zm;
                           zm,-C26',-C36',-C46',zm,zm]);
    
    %compile mass matrix
    [Me] = mapMatrixNonSym([M11,zm,zm,zm,M15,M16;
                           zm,M22,zm,M24,zm,zm;
                           zm,zm,M33,M34,zm,zm;
                           zm,M24',M34',M44,zm,zm;
                           M15',zm,zm,zm,M55,M56;
                           M16',zm,zm,zm,M56',M66]);

    %account for rayleigh damping
    alpha = input.RayleighAlpha;
    beta = input.RayleighBeta;
    
    CeRayleigh = alpha.*Kenr + beta.*Me;       
    Ce = Ce + CeRayleigh;
    
    %compile element force vector
    [Fe] = mapVector([F1;F2;F3;F4;F5;F6]);
    
   % transform matrices for sweep
   % Note,a negative sweep angle, will sweep away from the direction of
   % positive rotation
    lambdaTran = lambda';
    
    lambda = sparse(lambda);
    lambdaTran = sparse(lambdaTran);

    Me = lambdaTran*Me*lambda;
    Ce = lambdaTran*Ce*lambda;

    Ke = lambdaTran*Ke*lambda;
    if(useDisp && strcmp(iterationType,'NR'))
        Kehat =  lambdaTran*Kehat*lambda;
    end

    Fe = lambdaTran*Fe;

  %%
  
  %%concentrated mass 
    %NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
    % if some ycm or zcm offset from the node was accounted for in concentrated mass terms  
    concMassFlag = ~(isempty(find(concMass,1)));
    concStiffFlag = ~(isempty(find(concStiff,1)));
    concLoadFlag = ~(isempty(find(concLoad,1)));
    if(concMassFlag)
    %modify Me for concentrated mass
    Me(1,1) = Me(1,1) + concMass(1,1);
    Me(2,2) = Me(2,2) + concMass(1,1);
    Me(3,3) = Me(3,3) + concMass(1,1);
    Me(4,4) = Me(4,4) + concMass(2,1);
    Me(5,5) = Me(5,5) + concMass(3,1);
    Me(6,6) = Me(6,6) + concMass(4,1);
    
    Me(7,7) = Me(7,7) + concMass(1,2);
    Me(8,8) = Me(8,8) + concMass(1,2);
    Me(9,9) = Me(9,9) + concMass(1,2);
    Me(10,10) = Me(10,10) + concMass(2,2);
    Me(11,11) = Me(11,11) + concMass(3,2);
    Me(12,12) = Me(12,12) + concMass(4,2);
    
    %modify Ce for concentrated mass
    Ce(1,2) = Ce(1,2) - 2*concMass(1,1)*omega_z;
    Ce(2,1) = Ce(2,1) + 2*concMass(1,1)*omega_z;
    Ce(1,3) = Ce(1,3) + 2*concMass(1,1)*omega_y;
    Ce(3,1) = Ce(3,1) - 2*concMass(1,1)*omega_y;
    Ce(2,3) = Ce(2,3) - 2*concMass(1,1)*omega_x;
    Ce(3,2) = Ce(3,2) + 2*concMass(1,1)*omega_x;
    Ce(7,8) = Ce(7,8) - 2*concMass(1,2)*omega_z;
    Ce(8,7) = Ce(8,7) + 2*concMass(1,2)*omega_z;
    Ce(7,9) = Ce(7,9) + 2*concMass(1,2)*omega_y;
    Ce(9,7) = Ce(9,7) - 2*concMass(1,2)*omega_y;
    Ce(8,9) = Ce(8,9) - 2*concMass(1,2)*omega_x;
    Ce(9,8) = Ce(9,8) + 2*concMass(1,2)*omega_x;
    end
    
    if(concMassFlag || concStiffFlag)
    %modify Ke for concentrated mass
    Ke(1,1) = Ke(1,1) + concStiff(1,1) - concMass(1,1)*(omega_y^2 + omega_z^2);
    Ke(1,2) = Ke(1,2) + concMass(1,1)*omega_x*omega_y - concMass(1,1)*omegaDot_z;
    Ke(2,1) = Ke(2,1) + concMass(1,1)*omega_x*omega_y + concMass(1,1)*omegaDot_z;
    Ke(1,3) = Ke(1,3) + concMass(1,1)*omega_x*omega_z + concMass(1,1)*omegaDot_y;
    Ke(3,1) = Ke(3,1) + concMass(1,1)*omega_x*omega_z - concMass(1,1)*omegaDot_y;
    Ke(2,3) = Ke(2,3) + concMass(1,1)*omega_y*omega_z - concMass(1,1)*omegaDot_x;
    Ke(3,2) = Ke(3,2) + concMass(1,1)*omega_y*omega_z + concMass(1,1)*omegaDot_x;
    Ke(2,2) = Ke(2,2) + concStiff(2,1) - concMass(1,1)*(omega_x^2 + omega_z^2);
    Ke(3,3) = Ke(3,3) + concStiff(3,1) - concMass(1,1)*(omega_x^2 + omega_y^2);
    Ke(4,4) = Ke(4,4) + concStiff(4,1);
    Ke(5,5) = Ke(5,5) + concStiff(5,1);
    Ke(6,6) = Ke(6,6) + concStiff(6,1);
    Ke(7,7) = Ke(7,7) + concStiff(1,2) - concMass(1,2)*(omega_y^2 + omega_z^2);
    Ke(7,8) = Ke(7,8) + concMass(1,2)*omega_x*omega_y - concMass(1,2)*omegaDot_z;
    Ke(8,7) = Ke(8,7) + concMass(1,2)*omega_x*omega_y + concMass(1,2)*omegaDot_z;
    Ke(7,9) = Ke(7,9) + concMass(1,2)*omega_x*omega_z + concMass(1,2)*omegaDot_y;
    Ke(9,7) = Ke(9,7) + concMass(1,2)*omega_x*omega_z - concMass(1,2)*omegaDot_y;
    Ke(8,9) = Ke(8,9) + concMass(1,2)*omega_y*omega_z - concMass(1,2)*omegaDot_x;
    Ke(9,8) = Ke(9,8) + concMass(1,2)*omega_y*omega_z + concMass(1,2)*omegaDot_x;
    Ke(8,8) = Ke(8,8) + concStiff(2,2) - concMass(1,2)*(omega_x^2 + omega_z^2);
    Ke(9,9) = Ke(9,9) + concStiff(3,2) - concMass(1,2)*(omega_x^2 + omega_y^2);
    Ke(10,10) = Ke(10,10) + concStiff(4,2);
    Ke(11,11) = Ke(11,11) + concStiff(5,2);
    Ke(12,12) = Ke(12,12) + concStiff(6,2);
    end

    %modify Fe for  concentrated load
    if(concMassFlag || concLoadFlag)
    Fe(1) = Fe(1) + concLoad(1,1) + concMass(1,1)*(x(1)*(omega_y^2 + omega_z^2)-omega_x*omega_y*y(1) - omega_x*omega_z*z(1)) + concMass(1,1)*(y(1)*omegaDot_z-z(1)*omegaDot_y)  -  concMass(1,1)*a_x;
    Fe(2) = Fe(2) + concLoad(2,1) + concMass(1,1)*(y(1)*(omega_x^2 + omega_z^2)-omega_y*omega_z*z(1) - omega_y*omega_x*x(1)) + concMass(1,1)*(z(1)*omegaDot_x-x(1)*omegaDot_z)  -  concMass(1,1)*a_y;
    Fe(3) = Fe(3) + concLoad(3,1) + concMass(1,1)*(z(1)*(omega_x^2 + omega_y^2)-omega_z*omega_x*x(1) - omega_z*omega_y*y(1)) + concMass(1,1)*(x(1)*omegaDot_y-y(1)*omegaDot_x)  -  concMass(1,1)*a_z;
    Fe(4) = Fe(4) + concLoad(4,1);
    Fe(5) = Fe(5) + concLoad(5,1);
    Fe(6) = Fe(6) + concLoad(6,1); 
    Fe(7) = Fe(7) + concLoad(1,2) + concMass(1,2)*(x(2)*(omega_y^2 + omega_z^2)- omega_x*omega_y*y(2) - omega_x*omega_z*z(2)) + concMass(1,2)*(y(2)*omegaDot_z-z(2)*omegaDot_y) -  concMass(1,2)*a_x;
    Fe(8) = Fe(8) + concLoad(2,2) + concMass(1,2)*(y(2)*(omega_x^2 + omega_z^2)-omega_y*omega_z*z(2) - omega_y*omega_x*x(2)) + concMass(1,2)*(z(2)*omegaDot_x-x(2)*omegaDot_z)  -  concMass(1,2)*a_y;
    Fe(9) = Fe(9) + concLoad(3,2) + concMass(1,2)*(z(2)*(omega_x^2 + omega_y^2)-omega_z*omega_x*x(2) - omega_z*omega_y*y(2)) + concMass(1,2)*(x(2)*omegaDot_y-y(2)*omegaDot_x)  -  concMass(1,2)*a_z;
    Fe(10) = Fe(10) + concLoad(4,2);
    Fe(11) = Fe(11) + concLoad(5,2);
    Fe(12) = Fe(12) + concLoad(6,2);
    end
    
  
    %%

    if(strcmp(input.analysisType,'TD')) %calculate effective stiffness matrix and force vector for Dean integrator
  
    a1 = timeInt.a1;
    a2 = timeInt.a2;
    a3 = timeInt.a3;
    a4 = timeInt.a4;
    
    xn=disp(1:12); xnm1=dispm1(1:12);
    A = 2.0.*xn - xnm1;
    B = -a1.*xnm1 - a2.*xn;
    D = a3.*xnm1;
    
    Khate = Ke*a1 + a3.*Ce + Me;
    Fhate = Fe*a4 + Me*(A') + Ke*(B') + Ce*(D');
    
	FhatLessConc =   Fhate - [concLoad(1,1);
							  concLoad(2,1);
							  concLoad(3,1);
							  concLoad(4,1);
							  concLoad(5,1);
							  concLoad(6,1); 
							  concLoad(1,2);
							  concLoad(2,2);
							  concLoad(3,2);
                              concLoad(4,2);
							  concLoad(5,2);
							  concLoad(6,2)].*a4;
	
    %........................................................
    
    %..........................................................
    
    Ke = Khate;
    Fe = Fhate;
  
  
    end
    
    if(strcmp(input.analysisType,'TNB')) %calculate effective stiffness matrix and load vector for Newmark-Beta integrator
    a1 = timeInt.a1;
    a2 = timeInt.a2;
    a3 = timeInt.a3;
    a4 = timeInt.a4;
    a5 = timeInt.a5;
    a6 = timeInt.a6;
    a7 = timeInt.a7;
    a8 = timeInt.a8;

    u=disp; udot=dispdot; uddot=dispddot;
    if(strcmp(iterationType,'NR'))    %considerations if newton raphson iteration is used
        if(input.firstIteration)
            A = a3*u + a4*udot + a5*uddot;
            B = a6*u + a7*udot + a8*uddot;
            Fhate = Fe + Me*(A') + Ce*(B') - Ke*u';
        else
            Fhate = Fe  - Me*uddot' - Ce*udot' - (Ke)*u';
        end
    elseif(strcmp(iterationType,'DI')||strcmp(iterationType,'LINEAR'))   %considerations if direct iteration is used or linear analysis
           A = a3*u + a4*udot + a5*uddot;
           B = a6*u + a7*udot + a8*uddot;
           Fhate = Fe + Me*(A') + Ce*(B');
    end
     
        Khate = Ke + a3.*Me + a6.*Ce;
    if(strcmp(iterationType,'NR')) %considerations if newton raphson iteration is used
       Khate = Kehat + Khate; 
    end
     
	FhatLessConc =   Fhate - [concLoad(1,1);
							  concLoad(2,1);
							  concLoad(3,1);
							  concLoad(4,1);
							  concLoad(5,1);
							  concLoad(6,1); 
							  concLoad(1,2);
							  concLoad(2,2);
							  concLoad(3,2);
                              concLoad(4,2);
							  concLoad(5,2);
							  concLoad(6,2)];

    %........................................................
   
    Ke = Khate;
    Fe = Fhate;
    
    end
    
    if(strcmp(input.analysisType,'M'))
        FhatLessConc =   Fe - [concLoad(1,1);
							  concLoad(2,1);
							  concLoad(3,1);
							  concLoad(4,1);
							  concLoad(5,1);
							  concLoad(6,1); 
							  concLoad(1,2);
							  concLoad(2,2);
							  concLoad(3,2);
                              concLoad(4,2);
							  concLoad(5,2);
							  concLoad(6,2)];
                          
        output.FhatLessConc = FhatLessConc;
        
        if(strcmp(iterationType,'DI'))
           Fe = Fe*input.loadStep; 
        end
    end
    
    if((strcmp(input.analysisType,'M') || strcmp(input.analysisType,'S')) && strcmp(iterationType,'NR')) %considerations for newton-raphson iteration
       Fe = Fe*input.loadStep - Ke*disp_iter';
           Ke = Ke + Kehat;
    end
    
    
    
    
    %----- assign output block ----------------
    output.Ke = Ke;
    output.Fe = Fe;
        
    if(strcmp(input.analysisType,'M')||strcmp(input.analysisType,'RM0'))
        output.Me = Me;
        output.Ce = Ce;
    end
    
    if(strcmp(input.analysisType,'TD') || strcmp(input.analysisType,'TNB'))
        output.FhatLessConc = FhatLessConc;
    end
    %------------------------------------------

end

function [valGP] = interpolateVal(valNode,N)
%This function interpolates a value using distinct values at valNode
%and the corresponding shape function N.
    valGP = 0.0;
    for i=1:length(N)
        valGP = valGP + N(i)*valNode(i);
    end
end

%Element calculation functions---------------------------------

function [K] = calculateElement1(EA,integrationFactor,N1,N2,K)
%This function is a general routine to calculate an element matrix
len1 = length(N1);
len2 = length(N2);
    for i=1:len1
        for j=1:len2
        K(i,j) = K(i,j) + EA*N1(i)*N2(j)*integrationFactor;
        end
    end

end

function [F] = calculateVec1(f,integrationFactor,N,F)
%This function is a general routine to calculate an element vector
    len=length(N);
    for i=1:len
        F(i) = F(i) + f*N(i)*integrationFactor;
    end

end

function [Kel] = mapMatrixNonSym(Ktemp)
%----- function to form total stifness matrix and transform to desired
% DOF mapping

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

%map to FEA numbering
Kel = T'*Ktemp*T;

%declare map
% map = [1, 7, 2, 8, 3, 9,...
%       4, 10, 5, 11, 6, 12];
% 
% %map to FEA numbering
% for i=1:a
%     I=map(i);
%     for j=1:a
%         J=map(j);
%         Kel(I,J) = Ktemp(i,j);
%     end
% end

end

function [Fel] = mapVector(Ftemp)
%----- function to form total force vector and transform to desired
% DOF mapping
a=length(Ftemp);
Fel=zeros(a,1);
% 
% %declare map
map = [1, 7, 2, 8, 3, 9,...
      4, 10, 5, 11, 6, 12];

for i=1:a
    I=map(i);
    Fel(I) = Ftemp(i);
end

end
% %-------------------------------------------------------------------------


