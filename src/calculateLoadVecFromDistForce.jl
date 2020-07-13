function [output] = calculateLoadVecFromDistForce(input)
#calculateTimoshenkoElementNL performs nonlinear element calculations
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
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

#-------- assign input block ----------------
elementOrder   = input.elementOrder;
x              = input.x;
xloc           = input.xloc;

sectionProps   = input.sectionProps;
sweepAngle     = input.sweepAngle;
coneAngle      = input.coneAngle;
rollAngle      = input.rollAngle;

extDistF2Node = input.extDistF2Node;
extDistF3Node = input.extDistF3Node;
extDistF4Node = input.extDistF4Node;

#--------------------------------------------
numGP = 4;   #number of gauss points for full integration
#calculate quad points
[xi,weight] = getGP(numGP);

#Initialize element sub matrices and sub vectors
numNodesPerEl = length(x);

F1 = zeros(numNodesPerEl,1);
F3 = F1;
F2 = F1;
F4 = F1;
F5 = F1;
F6 = F1;

#Sort displacement vector
#Written for 2 node element with 6 dof per node
twistAvg = rollAngle + 0.5*(sectionProps.twist(1) + sectionProps.twist(2));
lambda = calculateLambda(sweepAngle*pi/180.0,coneAngle*pi/180.0,twistAvg.*pi/180.0);

#Integration loop
for i=1:numGP
    #Calculate shape functions at quad point i
    [N,~,Jac] = calculateShapeFunctions(elementOrder,xi(i),xloc);
    N1 = N;
    N2 = N;
    N3 = N;
    N4 = N;
    N5 = N;
    N6 = N;
    integrationFactor = Jac * weight(i);
    
    #..... interpolate for value at quad point .....
    extDistF1 = 0;
    extDistF2 = interpolateVal(extDistF2Node,N2);
    extDistF3 = interpolateVal(extDistF3Node,N3);
    extDistF4 = interpolateVal(extDistF4Node,N4);
    extDistF5 = 0;
    extDistF6 = 0;
    
    #distributed/body force load calculations
    [F1] = calculateVec1(extDistF1,integrationFactor,N1,F1);
    [F2] = calculateVec1(extDistF2,integrationFactor,N2,F2);
    [F3] = calculateVec1(extDistF3,integrationFactor,N3,F3);
    [F4] = calculateVec1(extDistF4,integrationFactor,N4,F4);
    [F5] = calculateVec1(extDistF5,integrationFactor,N5,F5);
    [F6] = calculateVec1(extDistF6,integrationFactor,N6,F6);
    
    
end #END OF INTEGRATION LOOP

#---------------------------------------------
#compile element force vector
[Fe] = mapVector([F1;F2;F3;F4;F5;F6]);

# transform matrices for sweep
# Note,a negative sweep angle, will sweep away from the direction of
# positive rotation
lambdaTran = lambda';
lambdaTran = sparse(lambdaTran);
Fe = lambdaTran*Fe;

#----- assign output block ----------------
output.Fe = Fe;
#------------------------------------------

end

function [valGP] = interpolateVal(valNode,N)
#This function interpolates a value using distinct values at valNode
#and the corresponding shape function N.
valGP = 0.0;
for i=1:length(N)
    valGP = valGP + N(i)*valNode(i);
end
end

#Element calculation functions---------------------------------
function [F] = calculateVec1(f,integrationFactor,N,F)
#This function is a general routine to calculate an element vector
len=length(N);
for i=1:len
    F(i) = F(i) + f*N(i)*integrationFactor;
end

end

function [Fel] = mapVector(Ftemp)
#----- function to form total force vector and transform to desired
# DOF mapping
a=length(Ftemp);
Fel=zeros(a,1);
#
# #declare map
map = [1, 7, 2, 8, 3, 9,...
    4, 10, 5, 11, 6, 12];

for i=1:a
    I=map(i);
    Fel(I) = Ftemp(i);
end

end
# #-------------------------------------------------------------------------



