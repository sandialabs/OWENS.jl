function [FgAero] = getAeroLoads(bladeData,aeroloadfile,time,sectionPropsArray,totalNumDof)
#getAeroLoads  gets aero loads from a CACTUS run to apply to structure
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [FgAero] = getAeroLoads(bladeData,aeroloadfile,time,sectionPropsArray,totalNumDof)
#
#   This function gets aero loads from a previous, uncoupled CACTUS run for
#   application to a structure
#
#   input:
#   bladeData          = blade data object
#   aeroloadfile      = .mat file containing aerodynamic loads
#   time              = current simulation time loads are requested at
#   sectionPropsArray = object containing arrays of section properties
#   totalNumDof       = total number of degrees of freedom in structural
#                       model
#
#   output:
#   FgAero      = force vector of aerodynamic loads

load(aeroloadfile);  #loads aerodynic forces
numDofPerNode = 6;
numBlades = bladeData.numBlades;
bladeNodes = bladeData.nodeNum;
FgAero = zeros(totalNumDof,1); #initializes aerodynamic laod vector
aeroStat = length(h);
index = 1;
gamma = zeros(aeroStat-1,1);
hAero = zeros(aeroStat-1,1);
for j=1:aeroStat-1
    gamma(index) = -atan2((h(j+1)-h(j)),(r(j+1)-r(j)))*180/pi; #h and r are loaded in with aeroload file
    hAero(index) = 0.5*(h(j)+h(j+1));
    index = index+1;
end

[numBladeEl,~] = size(Fn);
# numBlades = numBladeEl/(aeroStat-1);

#...... populate gamma array ........
gammaArray = [];
for i=1:numBlades
    gammaArray = cat(2,gammaArray,gamma);
end
gamma = gammaArray;
#.....................................
clear h;


#get xAC (aerodynamic center) at mid points of structural elements
numNodesPerBlade = length(bladeData.bladeNum)/bladeData.numBlades;
hStructural = bladeData.h;
xAC = zeros(length(hStructural),1);
hStructuralMid = zeros(length(hStructural),1);
index = 1;
for i=1:length(hStructural)
    if(bladeData.elementNum(i) != -1)
        xAC(index) = 0.5*sum(sectionPropsArray{bladeData.elementNum(i)}.aeroCenterOffset);
        hStructuralMid(index) = 0.5*(hStructural(i) + hStructural(i+1));
        index = index+1;
    end
end


#get normal and tangential forces on aero blade element
FtCurrent = zeros(numBladeEl,1);
FnCurrent = zeros(numBladeEl,1);
xACaero = zeros(numBladeEl,1);
Fh1 = zeros(numBladeEl,1);
Fh2 = zeros(numBladeEl,1);
Fh3 = zeros(numBladeEl,1);
MPitch = zeros(numBladeEl,1);
Mh1 = zeros(numBladeEl,1);
Mh2 = zeros(numBladeEl,1);
Mh3 = zeros(numBladeEl,1);
for i=1:numBladeEl
    FtCurrent(i) = interp1(t,Ft(i,:),time,'linear','extrap'); #normal and gangential loads at current time
    FnCurrent(i) = interp1(t,Fn(i,:),time,'linear','extrap');

    if(mod(i,aeroStat-1)!=0 && i!=1)
        aeroIndex = i - floor(i/(aeroStat-1))*(aeroStat-1);
    elseif(mod(i,aeroStat-1)==0)
        aeroIndex = aeroStat-1;
    else
        aeroIndex = i;
    end

    #map xAC to aero element domain
    xACaero(i) = interp1(hStructuralMid(1:length(hStructuralMid)/2),xAC(1:length(xAC)/2),hAero(aeroIndex),'linear','extrap');
    FlocalTemp = [0;-FtCurrent(i);FnCurrent(i)];

    #mapping value for 2 and 3 blades
    if(numBlades == 2)
        if(i/(aeroStat-1) <= 1)
            azi = 180.0;
        elseif(i/(aeroStat-1) <=2)
            azi = 0.0;
        end
    end

    if(numBlades == 3)
        if(i/(aeroStat-1) <= 1)
            azi = 180.0;
        elseif(i/(aeroStat-1) <=2)
            azi = -60.0;
        elseif(i/(aeroStat-1) <=3)
            azi = 60.0;
        end
    end

    #transform from blade element frame to the hub frame (still on aero element)
    #     lambda = calculateLambda(0.0,gamma(i)*pi/180, azi*pi/180.0); #this maps to the "hub" frame
    lambda = calculateLambda(azi*pi/180.0,gamma(i)*pi/180,0.0);
    FhubTemp = lambda(1:3,1:3)'*FlocalTemp; #'

    Fh1(i) =  FhubTemp(1);
    Fh2(i) =  FhubTemp(2);
    Fh3(i) =  FhubTemp(3);

    MPitch(i)=xACaero(i)*FnCurrent(i); #positive xAC is forward of elastic axis(towards leading edge)
    MlocalTemp = [MPitch(i);0;0];
    MhubTemp = lambda(1:3,1:3)'*MlocalTemp; #'

    Mh1(i) =  MhubTemp(1);
    Mh2(i) =  MhubTemp(2);
    Mh3(i) =  MhubTemp(3);

end

#map from aero element to structural element
for i=1:numBlades
    Fh1structural =  interp1(hAero,Fh1((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');
    Fh2structural =  interp1(hAero,Fh2((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');
    Fh3structural =  interp1(hAero,Fh3((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');
    Mh1structural =  interp1(hAero,Mh1((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');
    Mh2structural =  interp1(hAero,Mh2((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');
    Mh3structural =  interp1(hAero,Mh3((i-1)*(aeroStat-1)+1:i*(aeroStat-1)),hStructural((i-1)*(numNodesPerBlade-1)+1*i:i*(numNodesPerBlade)),'linear','extrap');

    #place aerodynamic loads in correct DOF of force vector
    singleBladeNodes = bladeNodes((i-1)*(numNodesPerBlade)+1:i*(numNodesPerBlade));
    for j=1:length(singleBladeNodes)
        Fxdof = (singleBladeNodes(j)-1)*numDofPerNode+1;
        Fydof = (singleBladeNodes(j)-1)*numDofPerNode+2;
        Fzdof = (singleBladeNodes(j)-1)*numDofPerNode+3;
        Mxdof = (singleBladeNodes(j)-1)*numDofPerNode+4;
        Mydof = (singleBladeNodes(j)-1)*numDofPerNode+5;
        Mzdof = (singleBladeNodes(j)-1)*numDofPerNode+6;

        FgAero(Fxdof) = Fh1structural(j);
        FgAero(Fydof) = Fh2structural(j);
        FgAero(Fzdof) = Fh3structural(j);
        FgAero(Mxdof) = Mh1structural(j);
        FgAero(Mydof) = Mh2structural(j);
        FgAero(Mzdof) = Mh3structural(j);
    end
end

end
