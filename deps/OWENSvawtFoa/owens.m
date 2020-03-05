function [freq,damp]=owens(varargin)

%owens Startup function for the OWENS toolkit
% **********************************************************************
% *                   Part of the SNL OWENS toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freq,damp]=owens(varargin)
%
%   This function is a start up function for launching various analysis
%   modes of the OWENS toolkit.
%
%      input:
%      varargin      = input parameter list
%         varargin{1} is the .owens file associated with analysis
%         varargin{2} is a string describing analysis type
%                     'S' = static analysis
%                     'M' = modal analysis
%                     'TNB' = transient analysis with Newmark-Beta time
%                     integration
%                     'ROM' = reduced order model transient analysis
%      output:
%      freq         = array of modal frequencies (when applicable to
%                     analysis type)
%      damp         = array of modal damping (when applicable to analysis
%                     type)
%      displ        = array containing converged solution for static
%                     displacement


inputfile = varargin{1};            %input file initialization
analysisType = varargin{2};         %anaysis type intialization
model.analysisType = analysisType;
model.turbineStartup = 0;           %initialization of turbine startup,
model.aeroElasticOn = false;        % aeroElastic flags, and air density
model.airDensity = 0;
model.gravityOn = true;             %flag to activate gravity loading in structural dynamics/static simulations

if(strcmp(analysisType,'S')) %STATIC ANALYSIS
   Omega = varargin{3};            %initialization of rotor speed (Hz)
   model.nlOn= varargin{4};        %flag for nonlinear elastic calculation
   if(length(varargin)>4)                %sets initial guess for nonlinear calculations
       displInitGuess = varargin{5};
   end
%    if(length(varargin)>5)                %sets air density if simple thin
%        model.airDensity = varargin{6};   % airfoil theory loading desired
%    else
       model.airDensity = 1.2041;
%    end


elseif(strcmp(analysisType,'M')) %MODAL ANALYSIS
   Omega = varargin{3};              %initialization of rotor speed (Hz)
   model.spinUpOn = varargin{4};     %flag for pre-stressed modal analysis
   model.nlOn = true;
   if(length(varargin)>4)
       model.numModesToExtract = varargin{5}; %number of modes to extract
   else                                       %lowest freq modes extracted first
       model.numModesToExtract = 20;
   end
   if(length(varargin)>5)
       displInitGuess = varargin{6};    %sets initial guess for nonlinear calculations
   end
%    if(length(varargin)>5)
%        model.airDensity = varargin{6};
%    else
       model.airDensity = 1.2041;
%    end

elseif(strcmp(analysisType,'TNB')||strcmp(analysisType,'TD')) %TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration)
   model.delta_t = varargin{3};  % time step size
   model.numTS = varargin{4};    % number of time steps
   model.nlOn = varargin{5};     % flag for nonlinear elastic calculation
   turbineOpFlag = varargin{6};           %turbine operation flag
   if(turbineOpFlag == 1) %generator start up operation mode
       model.turbineStartup = turbineOpFlag;
       model.OmegaInit = varargin{7};   %initial rotor speed (Hz)
   elseif(turbineOpFlag == 2) %self starting operation mode
       model.turbineStartup = turbineOpFlag;
       model.OmegaInit = varargin{7};   %initial rotor speed (Hz)
       model.OmegaGenStart = varargin{8}; %rotor speed at which generator activates (Hz)
   else                       %specified rotor speed profile
       model.turbineStartup = 0;
       if(length(varargin) == 6)
           model.usingRotorSpeedFunction = true; %set flag to use user specified rotor speed function
           [~,model.OmegaInit,~] = getRotorPosSpeedAccelAtTime(-1.0,0.0,0);
       else
           %this option uses a discretely specified rotor speed profile
           model.usingRotorSpeedFunction = false; %set flag to not use user specified rotor speed function
           model.tocp = varargin{7}; %time points for rotor speed provfile
           model.Omegaocp = varargin{8}; %rotor speed value at time points (Hz)
           model.OmegaInit = model.Omegaocp(1);
       end
   end

elseif(strcmp(analysisType,'ROM')) %REDUCED ORDER MODEL FOR TRANSIENT ANALYSIS
   model.delta_t = varargin{3}; %time step size
   model.numTS = varargin{4};   %number of time steps
   model.numModesForROM = varargin{5}; %number of lower system modes to include in ROM
   model.nlOn = varargin{6};    %flag for nonlinear elastic calculation
   turbineOpFlag = varargin{7};
   if(turbineOpFlag == 1) %generator start up operation mode
       model.turbineStartup = turbineOpFlag;
       model.OmegaInit = varargin{8}; %initial rotor speed
   elseif(turbineOpFlag == 2) % self starting operation mode
       model.turbineStartup = turbineOpFlag;
       model.OmegaInit = varargin{8}; %initial rotor speed (Hz)
       model.OmegaGenStart = varargin{9}; %rotor speed at which generator activates (Hz)
   else                         %specified rotor speed profile
       model.turbineStartup = 0;
       if(length(varargin) == 7)
           model.usingRotorSpeedFunction = true; %set flag to use user specified rotor speed function
           [~,model.OmegaInit,~] = getRotorPosSpeedAccelAtTime(-1.0,0.0,0);
       else
           %this option uses a discretely specified rotor speed profile
           model.usingRotorSpeedFunction = false; %set flag to not use user specified rotor speed function
           model.tocp = varargin{8}; %time points for rotor speed provfile
           model.Omegaocp = varargin{9}; %rotor speed value at time points (Hz)
           model.OmegaInit = model.Omegaocp(1);
       end
   end

elseif(strcmp(analysisType,'F'))  %MANUAL FLUTTER ANALYSIS
    Omega = varargin{3};   %rotor speed (Hz)
    model.spinUpOn = varargin{4}; %flag for pre-stressed modal analysis
    model.guessFreq = varargin{5}; %``guess'' modal frequency
    model.aeroElasticOn = true;

    if(length(varargin)>5)   %air density initialization
        model.airDensity = varargin{6};
    else
        model.airDensity = 1.2041;
    end
    if(length(varargin)>6)   %number of lower system modes to extract
        model.numModesToExtract = varargin{7};
    else
        model.numModesToExtract = 20;
    end

elseif(strcmp(analysisType,'FA')) %AUTOMATED FLUTTER ANALYSIS
   omegaArray = varargin{3};    %array of rotor speed values(Hz)
   model.spinUpOn = varargin{4}; %flag for pre-stressed modal analysis
   model.aeroElasticOn = true;

   if(length(varargin)>4)    %air density initializatio
       model.airDensity = varargin{5};
   else
       model.airDensity = 1.2041;
   end
   if(length(varargin)>5)    %number of lower system modes to extract
       model.numModesToExtract = varargin{6};
   else
       model.numModesToExtract = 20;
   end

else
    error('Analysis type not recognized.');
end


fid = fopen(inputfile,'r'); %reads in model file names from .owens file
last_delimiter = find(or(inputfile == '/', inputfile == '\')); %'
fdirectory = inputfile(1:last_delimiter(end));
meshfilename    = [fdirectory fscanf(fid,'%s',1)]; %mesh file name
eldatafilename  = [fdirectory fscanf(fid,'%s',1)]; %element data file name
ortdatafilename = [fdirectory fscanf(fid,'%s',1)]; %element orientation file name
jntdatafilename = [fdirectory fscanf(fid,'%s',1)]; %joint data file name
ndldatafilename = [fdirectory fscanf(fid,'%s',1)]; %concentrated nodal data file name
bcdatafilename  = [fdirectory fscanf(fid,'%s',1)]; %boundary condition file name
platformFlag    = fscanf(fid,'%i',1);
platfilename    = [fdirectory fscanf(fid,'%s',1)];
initcondfilename = [fdirectory fscanf(fid,'%s',1)]; %initial condition filename
aeroFlag        = fscanf(fid,'%i',1); %flag for activating aerodynamic analysis
blddatafilename = [fdirectory fscanf(fid,'%s',1)]; %blade data file name
model.aeroloadfile = [fdirectory fscanf(fid,'%s',1)]; %.mat file containing CACTUS aerodynamic loads

driveShaftFlag = fscanf(fid,'%i',1); %flag to include drive shaft effects
driveshaftfilename = [fdirectory fscanf(fid,'%s',1)]; %drive shaft file name

generatorfilename = [fdirectory fscanf(fid,'%s',1)]; %generator file name
rayleighDamping   = fscanf(fid,'%f',2); %read in alpha/beta for rayleigh damping
if(isempty(rayleighDamping))
    model.RayleighAlpha = 0.0;
    model.RayleighBeta = 0.0;
else
    model.RayleighAlpha = rayleighDamping(1);
    model.RayleighBeta = rayleighDamping(2);
end

fclose(fid); % close .owens file

%model definitions
model.elementOrder = 1; %linear element order
%--------------------------------------------
[mesh] = readMesh(meshfilename); %read mesh file
numDofPerNode = 6;
[model.BC] = readBCdata(bcdatafilename,mesh.numNodes,numDofPerNode); %read boundary condition file
[el] = readElementData(mesh.numEl,eldatafilename,ortdatafilename,blddatafilename); %read element data file (also reads orientation and blade data file associated with elements)
[model.joint] = readJointData(jntdatafilename); %read joint data file
rbarFileName = [inputfile(1:end-6),'.rbar']; %setrbarfile
[model.joint] = readRBarFile(rbarFileName,model.joint,mesh); %read rbar file name
[model.nodalTerms] = readNodalTerms(ndldatafilename); %read concentrated nodal terms file
[model] = readPlatformFile(model,platformFlag,platfilename);

if(strcmp(analysisType,'TNB')||strcmp(analysisType,'TD')||strcmp(analysisType,'ROM')) %for transient analysis...

    [model.initCond] = readInitCond(initcondfilename); %read initial conditions

    if(aeroFlag)
        model.aeroLoadsOn = true;
    else
        model.aeroLoadsOn = false;
    end

    [model.bladeData] = readBladeData(blddatafilename); %reads overall blade data file

    [model] = readDriveShaftProps(model,driveShaftFlag,driveshaftfilename); %reads drive shaft properties

    if(str2double(generatorfilename)==1.0)
        model.useGeneratorFunction = true;
    else
        model.useGeneratorFunction = false;
        [model.generatorProps] = readGeneratorProps(generatorfilename); %reads generator properties
    end

end

[model.outFilename] = generateOutputFilename(inputfile,analysisType); %generates an output filename for analysis results %TODO: map to the output location instead of input

[model.jointTransform,model.reducedDOFList] = createJointTransform(model.joint,mesh.numNodes,6); %creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints
[model.BC.map] = calculateBCMap(model.BC.numpBC,model.BC.pBC,numDofPerNode,model.reducedDOFList); %create boundary condition map from original DOF numbering to reduced/constrained DOF numbering
[~,numReducedDof] = size(model.jointTransform);
[model.BC.redVectorMap] = constructReducedDispVectorMap(mesh.numNodes,numDofPerNode,numReducedDof,model.BC); %create a map between reduced and full DOF lists

if(strcmp(analysisType,'S')) %EXECUTE STATIC ANALYSIS
   [model.nlParams] = readNLParamsFile(inputfile);
    if(length(varargin)<=4 || ~model.nlOn)                %sets initial guess for nonlinear calculations
       displInitGuess = zeros(mesh.numNodes*6,1);
   end

    OmegaStart = 0.0;
    staticExec(model,mesh,el,displInitGuess,Omega,OmegaStart);
end

if(strcmp(analysisType,'M') || strcmp(analysisType,'F')) %EXECUTE MODAL OR MANUAL FLUTTER ANALYSIS
   [model.nlParams] = readNLParamsFile(inputfile);
   if(length(varargin)<=5 || ~model.nlOn)
       displInitGuess = zeros(mesh.numNodes*6,1);
   end
    OmegaStart = 0.0;
    [freq,damp]=modalExec(model,mesh,el,displInitGuess,Omega,OmegaStart);
end

if(strcmp(analysisType,'FA')) %EXECUTE AUTOMATED FLUTTER ANALYSIS
    displ = zeros(mesh.numNodes*6,1);
    OmegaStart = 0.0;
    [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart);
end

if(strcmp(analysisType,'TNB')||strcmp(analysisType,'TD')||strcmp(analysisType,'ROM')) %EXECUTE TRANSIENT ANALYSIS
    [model.nlParams] = readNLParamsFile(inputfile);
    model.analysisType = analysisType;
    transientExec(model,mesh,el);
    freq=0;
end

end

function [outputfilename] = generateOutputFilename(inputfilename,analysisType)
    %This function generates an output file name depending on the analysis type

    %find the last '.' in inputfilename - helps to extract the prefix in the .owens
    index_all = find(inputfilename == '.');
    index = index_all(end);

    if(strcmp(analysisType,'M')||strcmp(analysisType,'F')||strcmp(analysisType,'FA')) %output filename (*.out) for modal/flutter analysis
        outputfilename = [inputfilename(1:index-1),'.out'];
    elseif(strcmp(analysisType,'S')) %output file name (*_static.mat) for static analysis
        outputfilename = [inputfilename(1:index-1),'_static.mat'];
    elseif(strcmp(analysisType,'TNB')||strcmp(analysisType,'TD')||strcmp(analysisType,'ROM')) %output filename (*.mat) for transient analysis
        outputfilename = [inputfilename(1:index-1),'.mat'];
    end

end

function [redVectorMap] = constructReducedDispVectorMap(numNodes,numDofPerNode,numReducedDof,BC)
    %This function creates a map of unconstrained DOFs between a full
    %listing and reduced listing (aftger constraints have been applied)

    bcdoflist=[];

    %create a listing of constrained DOFs from boundary condition file
    for i=1:BC.numpBC
        bcnodenum = BC.pBC(i,1);
        bcdofnum = BC.pBC(i,2);
        bcdoflist(i) = (bcnodenum-1)*numDofPerNode + bcdofnum;
    end

    dofList = calculateReducedDOFVector(numNodes,numDofPerNode,BC.isConstrained); %calculate a reduced (unconstrained) DOF vector

    redVectorMap = zeros(numReducedDof,1);

    for i=1:numReducedDof

        if(ismember(i,bcdoflist))              %creates a map of unconstrained reduced DOFs
             redVectorMap(i) = -1.0;
        else
            index = find(ismembc(dofList,i));
            redVectorMap(i) = index;
        end

    end

end

function [dofVector] = calculateReducedDOFVector(numNodes,numDofPerNode,isConstrained)
    %This function searches over all DOFs in a structural model and
    %determines and returns "dofVector" containing only unconstrained DOFs

    index = 1;
    dofVector=[];

    %loop over all DOFs in the model checking if constrained by BC or not
    for i=1:numNodes
       for j=1:numDofPerNode
          if(isConstrained((i-1)*numDofPerNode + j))
              constrained = true;
          else
              constrained = false;
          end

           if(constrained == false)
                dofVector(index) = (i-1)*numDofPerNode + j; %DOF vector only contains unconstrained DOFs
                index = index + 1;
           end
       end
    end

end
