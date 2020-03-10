function executeDynamic(rotorSpeedHz,timestep,nlFlag)

d=dir('*.owens');
inputFileName = d.name;
fnstring = inputFileName(1:end-6);
[numModesForROM] = getNumModesForROM(6);

oad('aeroLoads.mat', 'timeArray');
load('aeroLoads.mat', 'ForceValHist')
timeDuration = timeArray(end);
numTimesteps = timeDuration/timestep;
numTimesteps = floor(numTimesteps);

spinUpOn = nlFlag;
owens(inputFileName,'ROM',timestep,numTimesteps,numModesForROM,spinUpOn,0,[0,1e6],[rotorSpeedHz rotorSpeedHz]);
copyfile([fnstring,'.mat'],[fnstring,'_aeroLinNewWay.mat']);

% spinUpOn = true;
% owens(inputFileName,'ROM',timestep,numTimesteps,numModesForROM,spinUpOn,0,[0,1e6],[rotorSpeedHz rotorSpeedHz]);
% copyfile(dtemp.name,[dtemp.name(1:end-4),'_aeroNL.mat']);
% % 
% useAeroLoads = false;
% save aeroLoads timeArray ForceVal ForceDof useAeroLoads
% % 
% spinUpOn = false;
% owens(inputFileName,'TNB',timestep,numTimesteps,numModesForROM,spinUpOn,0,[0,1e6],[rotorSpeedHz rotorSpeedHz]);
% dtemp = dir('*vawt*t.mat');
% copyfile(dtemp.name,[dtemp.name(1:end-4),'_noAeroLin.mat']);
% 
% spinUpOn = true;
% owens(inputFileName,'ROM',timestep,numTimesteps,numModesForROM,spinUpOn,0,[0,1e6],[rotorSpeedHz rotorSpeedHz]);
% copyfile(dtemp.name,[dtemp.name(1:end-4),'_noAeroNL.mat']);

end

function [numModesForRom] = getNumModesForROM(numDofPerNode)
	
	d1 = dir('*.mesh');
	d2 = dir('*.jnt');
	d3 = dir('*.bc');
	
	mesh = readMesh(d1.name);
    jnt = readJointData(d2.name);
    numJoint = size(jnt,1);
    bc = readBCdata(d3.name,mesh.numNodes,numDofPerNode);
    
    %for now assume all joints are weld constraints
    numModesForRom = mesh.numNodes*numDofPerNode - numJoint*numDofPerNode - bc.numpBC;
    
end