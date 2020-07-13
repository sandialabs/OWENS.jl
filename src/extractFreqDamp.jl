function [freq,damp,phase1,phase2,sortedModes] = extractFreqDamp(val,vec,numDOFPerNode,jointTransform,reducedDOFList,BC,analysisType)
#extractFreqDamp   extract frequency, damping, mode shapes from eigsolution
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [freq,damp,phase1,phase2,sortedModes] = extractFreqDamp(val,vec,...
#                                         numDOFPerNode,jointTransform,...
#                                         reducedDOFList,BC,analysisType)
#
#   This function calculates the eigenvalues and vectors of a structural
#   dynamic system.
#
#   input:
#   val            = eigenvalue
#   vec            = eigenvector
#   numDOFPerNode  = number of degrees of freedom per node
#   jointTransform = joint transformation matrix from reduced to full DOF
#                    list
#   reducedDOFList = listing of reduced DOFs
#   BC             = boundary condition object containing boundary
#                    condition info
#   analysisType   = analysis type
#
#   output:
#   freq        = modal frequency
#   damp        = modal damping
#   phase1      = in phase mode shape (real part of mode shape)
#   phase2      = out of phase mode shape (imaginary part of mode shape)
#   sortedModes = total, complex mode shape

freq = abs(imag(val))/(2*pi);   #damped frequency for state space system
damp = -real(val)/abs(imag(val)); #modal damping

if(abs(imag(val)) < 1.0e-4)     #if imaginary part of eigenvalue is numeric zero treat as spring-mass system
    freq = sqrt(abs(real(val)))/(2*pi);
    damp = 0.0;
end

# if(~strcmp(analysisType,'FA'))  #for all but automated flutter analysis
#          [len,numModeShapes] = size(vec);
[dispReduc] = constructReducedDispVecFromEigVec(vec,reducedDOFList,BC); #construct mode shape vector with boundary conditinos
dispOrig = jointTransform*dispReduc; #transform from reduced DOF list to full DOF list
lenOrig=length(dispOrig);

sortedModes0=zeros(lenOrig/numDOFPerNode,numDOFPerNode);
sortedModes = complex(sortedModes0,0);
for i=1:lenOrig/numDOFPerNode     #construct matrix of nodal DOF values from full DOF eigenvector
    for j=1:numDOFPerNode
        sortedModes(i,j) = dispOrig((i-1)*6+j);
    end
end

phase1 = real(sortedModes);  #phase 1 is real part of modeshape (0 deg in phase)
phase2 = imag(sortedModes);  #phase 2 is imag part of modeshape (90 deg out of phase)

max1=max(max(abs(phase1))); #find maximum values for modeshape normalization
max2=max(max(abs(phase2)));


if(max1>max2)
    maxOverall = max1;
else
    maxOverall = max2;
end

if(maxOverall == 0)
    maxOverall = 1;
end

phase1 = phase1./maxOverall;  #normalize modeshapes
phase2 = phase2./maxOverall;

if(abs(min(min(phase1))+1)<1.0e-4 || abs(min(min(phase2))+1)<1.0e-4)
    phase1 = -1*phase1;
    phase2 = -1*phase2;
    
end

# else  #return null mode shapes if mode shapes not requested
#     phase1 = [];
#     phase2 = [];
#     sortedModes = [];
# end

end


function [vec1Red] = constructReducedDispVecFromEigVec(vec1,reducedDOFList,BC)
#This function takes the original mode shape and modifies it to
#account for boundary conditions
bcdoflist=zeros(1,BC.numpBC);
#form pBC DOF list
for i=1:BC.numpBC
    bcnodenum = BC.pBC(i,1);
    bcdofnum = BC.pBC(i,2);
    bcdoflist(i) = (bcnodenum-1)*6 + bcdofnum;
end

index = 1;
len = length(vec1);
vec1Red0 = zeros(length(reducedDOFList),1);
vec1Red = complex(vec1Red0,0);
for i=1:length(reducedDOFList)
    if(~ismember(reducedDOFList(i),bcdoflist))
        vec1Red(i) = vec1(len/2+index,1);
        index = index + 1;
    else
        vec1Red(i) = 0;
    end
end

end

