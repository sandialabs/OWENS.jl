path = splitdir(@__FILE__)[1]

import OWENS
import OWENSFEA

using Test
import HDF5
import DelimitedFiles
import FLOWMath


## DEFINE helper functions/structs

# Cactus Related Structs
"""
Internal, struct containing the CACTUS geometry file data
"""
mutable struct CactusGeom
    NBlade
    NStrut
    RotN
    RotP
    RefAR
    RefR
    blade
    strut
end

"""
    owens(owensfile,analysisType;
        delta_t=2e-3,
        numTS=100,
        tocp=[0.0,1.1],
        Omegaocp=[0.0,1.0],
        OmegaInit=0.0,
        OmegaGenStart=0.0,
        usingRotorSpeedFunction=false,
        nlOn=true,
        Omega=0.0,
        turbineStartup=0,
        spinUpOn=false,
        numModesToExtract=20,
        displInitGuess=0.0,
        airDensity=1.2041,
        aeroElasticOn = false,
        guessFreq = 0,
        gravityOn = true,
        generatorOn = false,
        omegaControl = false,
        iterationType = "NR", # nlParams
        adaptiveLoadSteppingFlag = true,
        tolerance = 1.0000e-06,
        maxIterations = 50,
        maxNumLoadSteps = 20,
        minLoadStepDelta = 0.0500,
        minLoadStep = 0.0500,
        prescribedLoadStep = 0.0,
        elementOrder = 1,
        numDOFPerNode = 6,
        platformActive = false,
        platformTurbineConnectionNodeNumber = 1,
        JgearBox =0.0,
        gearRatio = 1.0,
        gearBoxEfficiency = 1.0,
        useGeneratorFunction = false,
        generatorProps = 0.0,
        driveTrainOn = false)

Original serial and file reading method of running an analysis.

#Inputs
See ?OWENS.Model and ?OWENSFEA.FEAModel

#Outputs
See ?OWENS.Unsteady, ?OWENSFEA.Modal

"""
function owens(owensfile,analysisType;
    delta_t=2e-3,
    numTS=100,
    tocp=[0.0,1.1],
    Omegaocp=[0.0,1.0],
    OmegaInit=0.0,
    OmegaGenStart=0.0,
    usingRotorSpeedFunction=false,
    nlOn=true,
    Omega=0.0,
    turbineStartup=0,
    spinUpOn=false,
    numModesToExtract=20,
    displInitGuess=0.0,
    airDensity=1.2041,
    aeroElasticOn = false,        # aeroElastic flags, and air density,
    guessFreq = 0,          #``guess"" modal frequency
    gravityOn = true,             #flag to activate gravity loading in structural dynamics/static simulations,
    generatorOn = false,
    omegaControl = false,
    iterationType = "NR", # nlParams
    adaptiveLoadSteppingFlag = true,
    tolerance = 1.0000e-06,
    maxIterations = 50,
    maxNumLoadSteps = 20,
    minLoadStepDelta = 0.0500,
    minLoadStep = 0.0500,
    prescribedLoadStep = 0.0,
    elementOrder = 1, #linear element order
    numDOFPerNode = 6,
    platformActive = false,
    interpOrder = 2,
    platformTurbineConnectionNodeNumber = 1,
    JgearBox =0.0,
    gearRatio = 1.0,             #set gear ratio and efficiency to 1
    gearBoxEfficiency = 1.0,
    useGeneratorFunction = false,
    generatorProps = 0.0,
    driveTrainOn = false,          #set drive shaft unactive
    hydrodynLib = "none",
    moordynLib = "none")

    # if(analysisType=="S") #STATIC ANALYSIS
    #     Omega = varargin{3}            #initialization of rotor speed (Hz)
    #     inputs.nlOn= varargin{4}        #flag for nonlinear elastic calculation
    #     if(length(varargin)>4)                #sets initial guess for nonlinear calculations
    #         displInitGuess = varargin{5}
    #     end
    #     #    if(length(varargin)>5)                #sets air density if simple thin
    #     #        inputs.airDensity = varargin{6}   # airfoil theory loading desired
    #     #    else
    #     inputs.airDensity = 1.2041
    #     #    end

    if (analysisType=="TNB" || analysisType=="TD") #TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration)
        if usingRotorSpeedFunction
            _,OmegaInit,_ = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
        else
            #this option uses a discretely specified rotor speed profile
            OmegaInit = Omegaocp[1] #TODO: simplify this downstream since the data doesn't need to be repeated
        end
    end

    # elseif(analysisType=="ROM") #REDUCED ORDER inputs FOR TRANSIENT ANALYSIS
    #     inputs.delta_t = varargin{3} #time step size
    #     inputs.numTS = varargin{4}   #number of time steps
    #     inputs.numModesForROM = varargin{5} #number of lower system modes to include in ROM
    #     inputs.nlOn = varargin{6}    #flag for nonlinear elastic calculation
    #     turbineOpFlag = varargin{7}
    #     if(turbineOpFlag == 1) #generator start up operation mode
    #         inputs.OmegaInit = varargin{8} #initial rotor speed
    #     elseif(turbineOpFlag == 2) # self starting operation mode
    #         inputs.OmegaInit = varargin{8} #initial rotor speed (Hz)
    #         inputs.OmegaGenStart = varargin{9} #rotor speed at which generator activates (Hz)
    #     else                         #specified rotor speed profile
    #         if(length(varargin) == 7)
    #             inputs.usingRotorSpeedFunction = true #set flag to use user specified rotor speed function
    #             [~,inputs.OmegaInit,~] = getRotorPosSpeedAccelAtTime(-1.0,0.0,0)
    #         else
    #             #this option uses a discretely specified rotor speed profile
    #             inputs.usingRotorSpeedFunction = false #set flag to not use user specified rotor speed function
    #             inputs.tocp = varargin{8} #time points for rotor speed provfile
    #             Omegaocp = varargin{9} #rotor speed value at time points (Hz)
    #             inputs.Omegaocp = Omegaocp
    #             inputs.OmegaInit = Omegaocp(1)
    #         end
    #     end
    # elseif(analysisType=="F")  #MANUAL FLUTTER ANALYSIS
    #     Omega = varargin{3}   #rotor speed (Hz)
    #     inputs.spinUpOn = varargin{4} #flag for pre-stressed modal analysis
    #     inputs.guessFreq = varargin{5} #``guess"" modal frequency
    #     inputs.aeroElasticOn = true
    #     inputs.nlOn = true
    #
    #     if(length(varargin)>5)   #air density initialization
    #         inputs.airDensity = varargin{6}
    #     else
    #         inputs.airDensity = 1.2041
    #     end
    #     if(length(varargin)>6)   #number of lower system modes to extract
    #         inputs.numModesToExtract = varargin{7}
    #     else
    #         inputs.numModesToExtract = 20
    #     end
    #
    # elseif(analysisType=="FA") #AUTOMATED FLUTTER ANALYSIS
    #     omegaArray = varargin{3}    #array of rotor speed values(Hz)
    #     inputs.spinUpOn = varargin{4} #flag for pre-stressed modal analysis
    #     inputs.aeroElasticOn = true
    #     inputs.nlOn = true
    #
    #     if(length(varargin)>4)    #air density initializatio
    #         inputs.airDensity = varargin{5}
    #     else
    #         inputs.airDensity = 1.2041
    #     end
    #     if(length(varargin)>5)    #number of lower system modes to extract
    #         inputs.numModesToExtract = varargin{6}
    #     else
    #         inputs.numModesToExtract = 20
    #     end
    #
    # else
    #     error("Analysis type not recognized.")
    # end

    fid = open(owensfile,"r") #reads in inputs file names from .owens file

    last_delimiter   = findall(r"\W", owensfile) #Find the file directory
    fdirectory       = owensfile[1:last_delimiter[end-1][1]]

    meshfilename     = string(fdirectory, readline(fid)) #mesh file name
    eldatafilename   = string(fdirectory, readline(fid)) #element data file name
    ortdatafilename  = string(fdirectory, readline(fid)) #element orientation file name
    jntdatafilename  = string(fdirectory, readline(fid)) #joint data file name
    ndldatafilename  = string(fdirectory, readline(fid)) #concentrated nodal data file name
    bcdatafilename   = string(fdirectory, readline(fid)) #boundary condition file name
    line             = readline(fid)
    platformFlag     = real(parse(Int,line[1:2]))
    platfilename     = string(fdirectory, line[3:end])

    initcondfilename = string(fdirectory, readline(fid)) #initial condition filename

    line             = readline(fid)
    delimiter_idx    = findall(" ",line)

    aeroLoadsOn      = Bool(real(parse(Int,line[1]))) #flag for activating aerodynamic analysis

    blddatafilename  = string(fdirectory, line[delimiter_idx[1][1]+1:delimiter_idx[2][1]-1]) #blade data file name
    aeroloadfile = string(fdirectory, line[delimiter_idx[2][1]+1:end]) #.csv file containing CACTUS aerodynamic loads
    if analysisType!="M"
        owensfile = string(blddatafilename[1:end-4], ".owens") #TODO: this is redundant and confusing since it is specified at the beginning, clean up
    end
    line             = readline(fid) #flag to include drive shaft effects
    driveShaftFlag   = real(parse(Int,line[1:2]))
    driveshaftfilename = string(fdirectory, line[3:end]) #drive shaft file name

    generatorfilename = string(fdirectory, readline(fid)) #generator file name

    line = readline(fid)
    delimiter_idx = findall(" ",line)
    platformActive = Bool(real(parse(Int,line[1]))) #flag for activating hydrodynamic analysis
    potflowfile  = string(fdirectory, line[delimiter_idx[1][1]+1:delimiter_idx[2][1]-1]) # potential flow file prefix
    interpOrder = real(parse(Int,line[delimiter_idx[2][1]+1])) # interpolation order for HD/MD libraries
    line = readline(fid)
    rayleighDamping = split(line)

    if (isempty(rayleighDamping))
        RayleighAlpha = 0.0
        RayleighBeta = 0.0
    else
        RayleighAlpha = parse(Float64,rayleighDamping[1])
        RayleighBeta = parse(Float64,rayleighDamping[2])
    end

    close(fid) # close .owens file

    #inputs definitions
    #--------------------------------------------
    mesh = OWENS.readMesh(meshfilename) #read mesh file
    bladeData,_,_,_ = OWENS.readBladeData(blddatafilename) #reads overall blade data file
    BC = OWENS.readBCdata(bcdatafilename,mesh.numNodes,numDOFPerNode) #read boundary condition file
    el = OWENS.readElementData(mesh.numEl,eldatafilename,ortdatafilename,bladeData) #read element data file (also reads orientation and blade data file associated with elements)
    joint = DelimitedFiles.readdlm(jntdatafilename,'\t',skipstart = 0) #readJointData(jntdatafilename) #read joint data file
    nodalTerms = OWENSFEA.applyConcentratedTerms(mesh.numNodes, 6;filename=ndldatafilename) #read concentrated nodal terms file
    initCond = []

    #     [inputs] = readDriveShaftProps(inputs,driveShaftFlag,driveshaftfilename) #reads drive shaft properties

    driveShaftProps = OWENS.DriveShaftProps(0.0,0.0)       #set drive shat properties to 0

    if (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #for transient analysis...

        if !(occursin("[",generatorfilename)) #If there isn't a file
            useGeneratorFunction = true
            generatorProps = 0.0
        else
            useGeneratorFunction = false
            generatorProps = OWENS.readGeneratorProps(generatorfilename) #reads generator properties
        end

    end

    dataOutputFilename = generateOutputFilename(owensfile,analysisType) #generates an output filename for analysis results #TODO: map to the output location instead of input
    jointTransform, reducedDOFList = OWENSFEA.createJointTransform(joint,mesh.numNodes,6) #creates a joint transform to constrain inputs degrees of freedom (DOF) consistent with joint constraints
    numReducedDof = length(jointTransform[1,:])

    nlParams = OWENSFEA.NlParams(iterationType,adaptiveLoadSteppingFlag,tolerance,
    maxIterations,maxNumLoadSteps,minLoadStepDelta,minLoadStep,prescribedLoadStep)

    inputs = OWENS.Inputs(;analysisType,turbineStartup,usingRotorSpeedFunction,tocp,numTS,delta_t,Omegaocp,
    aeroLoadsOn,driveTrainOn,generatorOn,platformActive=false,topsideOn=true,interpOrder,hd_input_file=[],md_input_file=[],JgearBox,gearRatio,gearBoxEfficiency,
    useGeneratorFunction,generatorProps,OmegaGenStart,omegaControl,OmegaInit,
    aeroloadfile,owensfile,potflowfile=[],dataOutputFilename,bladeData,driveShaftProps)

    feamodel = OWENSFEA.FEAModel(;analysisType,
    initCond,
    aeroElasticOn,
    guessFreq,
    airDensity,
    gravityOn,
    nlOn,
    spinUpOn,
    dataOutputFilename,
    RayleighAlpha,
    RayleighBeta,
    elementOrder,
    joint,
    platformTurbineConnectionNodeNumber,
    jointTransform,
    reducedDOFList,
    nlParams,
    numNodes = mesh.numNodes,
    numModes = numModesToExtract,
    pBC=BC.pBC,
    nodalTerms)

    #     if(analysisType=="S") #EXECUTE STATIC ANALYSIS
    #         if(length(varargin)<=4 || ~inputs.nlOn)                #sets initial guess for nonlinear calculations
    #             displInitGuess = zeros(mesh.numNodes*6,1)
    #         end
    #
    #         OmegaStart = 0.0
    #         staticExec(inputs,mesh,el,displInitGuess,Omega,OmegaStart)
    #     end
    #
    if (analysisType == "M" || analysisType == "F") #EXECUTE MODAL OR MANUAL FLUTTER ANALYSIS
        if (displInitGuess==0 || !nlOn)
            displInitGuess = zeros(mesh.numNodes*6)
        end
        OmegaStart = 0.0
        freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90=OWENSFEA.modal(feamodel,mesh,el;displ=displInitGuess,Omega,OmegaStart)
        return freq,damp,imagCompSign,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90
    end
    #
    #     if(analysisType=="FA") #EXECUTE AUTOMATED FLUTTER ANALYSIS
    #         displ = zeros(mesh.numNodes*6,1)
    #         OmegaStart = 0.0
    #         [freq,damp]=ModalAuto(inputs,mesh,el,displ,omegaArray,OmegaStart)
    #     end
    #

    if ((hydrodynLib!="none")||(moordynLib!="none")) # Map shared library paths for external dependencies
        bin = OWENS.Bin(hydrodynLib, moordynLib)
    else
        bin = nothing
    end

    if (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #EXECUTE TRANSIENT ANALYSIS
        aeroLoadsFile_root = inputs.aeroloadfile[1:end-16] #cut off the _ElementData.csv
        OWENSfile_root = inputs.owensfile[1:end-6] #cut off the .owens

        geomFn = string(aeroLoadsFile_root, ".geom")
        loadsFn = string(aeroLoadsFile_root, "_ElementData.csv")
        bldFn = string(OWENSfile_root, ".bld")
        elFn = string(OWENSfile_root, ".el")
        ortFn = string(OWENSfile_root, ".ort")
        meshFn = string(OWENSfile_root, ".mesh")

        aerotimeArray,aeroForceValHist,aeroForceDof,cactusGeom = mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)

        aeroForces(t,azi) = OWENS.externalForcing(t+delta_t,aerotimeArray,aeroForceValHist,aeroForceDof)

        deformAero(azi;newOmega=-1,newVinf=-1,bld_x=-1,bld_z=-1,bld_twist=-1,accel_flap_in=-1,accel_edge_in=-1,gravity=-1) = 0.0 #placeholder function
        
        t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,
        FReactionHist,FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,
        uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,
        kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
        topFexternal_hist,rbDataHist = OWENS.Unsteady(inputs;topModel=feamodel,topMesh=mesh,topEl=el,bin,aero=aeroForces,deformAero)
            
        OWENS.outputData(;mymesh=mesh,inputs,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist)

        return inputs
    end

    #
    # end
    #

end

"""
    readCactusGeom(geom_fn)

Reads a CACTUS geometry file and creates a CactusGeom structure containing blade and strut geometry data.

The CACTUS geometry file format is as follows:
- First 6 lines contain global parameters:
  1. Number of blades
  2. Number of struts
  3. Rotor normal vector (3 components)
  4. Rotor position vector (3 components)
  5. Reference area
  6. Reference radius
- Following sections contain blade data (repeated for each blade):
  - Number of elements
  - Flip normal flag
  - Quarter chord points (x,y,z)
  - Tangent vectors (x,y,z)
  - Chord to radius ratios
  - Panel edge points (x,y,z)
  - Tangent vectors at edges (x,y,z)
  - Normal vectors at edges (x,y,z)
  - Span vectors at edges (x,y,z)
  - Element chord to radius ratios
  - Element area ratios
  - Section indices
- Final sections contain strut data (repeated for each strut):
  - Number of elements
  - Thickness to chord ratio
  - Mid-chord points (x,y,z)
  - Chord to radius ratios
  - Panel edge points (x,y,z)
  - Span vectors at edges (x,y,z)
  - Element chord to radius ratios
  - Element area ratios
  - Blade indices (start/end)

# Arguments
* `geom_fn::String`: Path to the CACTUS geometry file

# Returns
* `CactusGeom`: Structure containing:
  - `NBlade`: Number of blades
  - `NStrut`: Number of struts
  - `RotN`: Rotor normal vector
  - `RotP`: Rotor position vector
  - `RefAR`: Reference area
  - `RefR`: Reference radius
  - `blade`: Array of Blade structures
  - `strut`: Array of Strut structures
"""
function readCactusGeom(geom_fn)

    data = DelimitedFiles.readdlm(geom_fn,skipstart = 0)

    NBlade = Int(data[1,2])
    NStrut = Int(data[2,2])
    RotN = Float64.(data[3,2:4])
    RotP = Float64.(data[4,2:4])
    RefAR = Float64(data[5,2])
    RefR = Float64(data[6,2])

    blade = Array{OWENS.Blade, 1}(undef, NBlade)

    idx = 9
    for bld_idx = 1:NBlade

        NElem = Int(data[idx+0,2])
        FlipN = Int(data[idx+1,2])
        QCx = Float64.(data[idx+2,2:NElem+2])
        QCy = Float64.(data[idx+3,2:NElem+2])
        QCz = Float64.(data[idx+4,2:NElem+2])
        tx = Float64.(data[idx+5,2:NElem+2])
        ty = Float64.(data[idx+6,2:NElem+2])
        tz = Float64.(data[idx+7,2:NElem+2])
        CtoR = Float64.(data[idx+8,2:NElem+2])
        PEx = Float64.(data[idx+9,2:NElem+1])
        PEy = Float64.(data[idx+10,2:NElem+1])
        PEz = Float64.(data[idx+11,2:NElem+1])
        tEx = Float64.(data[idx+12,2:NElem+1])
        tEy = Float64.(data[idx+13,2:NElem+1])
        tEz = Float64.(data[idx+14,2:NElem+1])
        nEx = Float64.(data[idx+15,2:NElem+1])
        nEy = Float64.(data[idx+16,2:NElem+1])
        nEz = Float64.(data[idx+17,2:NElem+1])
        sEx = Float64.(data[idx+18,2:NElem+1])
        sEy = Float64.(data[idx+19,2:NElem+1])
        sEz = Float64.(data[idx+20,2:NElem+1])
        ECtoR = Float64.(data[idx+21,2:NElem+1])
        EAreaR = Float64.(data[idx+22,2:NElem+1])
        iSect = Float64.(data[idx+23,2:NElem+1])
        idx += 25
        blade[bld_idx] = OWENS.Blade(NElem,FlipN,QCx,QCy,QCz,tx,ty,tz,CtoR,PEx,PEy,PEz,tEx,tEy,tEz,nEx,nEy,nEz,sEx,sEy,sEz,ECtoR,EAreaR,iSect)
    end

    strut = Array{OWENS.Strut, 1}(undef, NStrut)

    for strut_idx = 1:NStrut

        NElem = Int(data[idx+0,2])
        TtoC = Float64(data[idx+1,2])
        MCx = Float64.(data[idx+2,2:NElem+2])
        MCy = Float64.(data[idx+3,2:NElem+2])
        MCz = Float64.(data[idx+4,2:NElem+2])
        CtoR = Float64.(data[idx+5,2:NElem+2])
        PEx = Float64.(data[idx+6,2:NElem+1])
        PEy = Float64.(data[idx+7,2:NElem+1])
        PEz = Float64.(data[idx+8,2:NElem+1])
        sEx = Float64.(data[idx+9,2:NElem+1])
        sEy = Float64.(data[idx+10,2:NElem+1])
        sEz = Float64.(data[idx+11,2:NElem+1])
        ECtoR = Float64.(data[idx+12,2:NElem+1])
        EAreaR = Float64.(data[idx+13,2:NElem+1])
        BIndS = Int(data[idx+14,2])
        EIndS = Int(data[idx+15,2])
        BIndE = Int(data[idx+16,2])
        EIndE = Int(data[idx+17,2])

        idx += 19
        strut[strut_idx] = OWENS.Strut(NElem,TtoC,MCx,MCy,MCz,CtoR,PEx,PEy,PEz,sEx,sEy,sEz,ECtoR,EAreaR,BIndS,EIndS,BIndE,EIndE)
    end

    return CactusGeom(NBlade,NStrut,RotN,RotP,RefAR,RefR,blade,strut)
end

"""
Internal, takes cactus loads and geometry and OWENS mesh and maps the loads to the blades' FEA dofs
"""
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
    V = 25.0#6.787243728 #m/s #27.148993200000003
    println("Velocity for the cacus test case is hardcoded here at $V")

    normTime = aero_data[:,1]

    numAeroEl = 0
    for i=1:cactusGeom.NBlade
        numAeroEl = numAeroEl + cactusGeom.blade[i].NElem
    end

    len,_ = size(aero_data)

    numAeroTS = Int(len/numAeroEl)
    # time = zeros(len/numAeroEl,1)
    time = normTime[1:Int(numAeroEl):end,1].*RefR[1]./V[1]

    urel = aero_data[:,15]
    uloc = urel.*V

    cn = aero_data[:,24]
    ct = aero_data[:,25]
    cm25 = aero_data[:,22]

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

    bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers = OWENS.readBladeData(bldFn)

    #Initialize structuralLoad
    struct_N = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_T = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))
    struct_M25 = zeros(cactusGeom.NBlade,numAeroTS,length(structuralElNumbers[1,:]))

    # TODO: the interpolation input here is blade height vs blade span length, should be updated to be consistent
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
    mesh = OWENS.readMesh(meshFn)
    el = OWENS.readElementData(mesh.numEl,elFn,ortFn,bladeData)
    numDOFPerNode = 6
    #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
    Fg = zeros(mesh.numNodes*6,numAeroTS)
    for i=1:numAeroTS
        for j = 1:cactusGeom.NBlade
            for k = 1:numNodesPerBlade-1
                #get element aero_data
                # orientation angle,xloc,sectionProps,element order]
                elNum = Int(structuralElNumbers[j,k])
                #get dof map
                node1 = Int(structuralNodeNumbers[j,k])
                node2 = Int(structuralNodeNumbers[j,k+1])
                dofList = [(node1-1)*numDOFPerNode.+(1:6) (node2-1)*numDOFPerNode.+(1:6)]

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

                Fe = OWENSFEA.calculateLoadVecFromDistForce(elementOrder,x,xloc,twist,sweepAngle,coneAngle,rollAngle,extDistF2Node,extDistF3Node,extDistF4Node)

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
    # ForceValHist = zeros(sum(Fg[:,1].!=0),length(Fg[1,:]))
    # ForceDof = zeros(sum(Fg[:,1].!=0),1)
    ForceValHist = zeros(length(Fg[:,1]),length(Fg[1,:]))
    ForceDof = zeros(length(Fg[:,1]),1)
    index = 1
    for i=1:Int(mesh.numNodes*6)
        # if !isempty(findall(x->x!=0,Fg[i,:]))

            ForceValHist[index,:] = Fg[i,:]
            ForceDof[index] = i
            index = index + 1
        # end
    end
    return time,ForceValHist,ForceDof,cactusGeom
end

"""
Internal, generates an output file name depending on the analysis type
"""
function generateOutputFilename(owensfilename,analysisType)

    #find the last "." in owensfilename - helps to extract the prefix in the .owens
    index = findlast(".",owensfilename)[1]

    if (analysisType=="M"||analysisType=="F"||analysisType=="FA") #output filename (*.out) for modal/flutter analysis
        outputfilename = string(owensfilename[1:index-1],".out")
    elseif (analysisType=="S") #output file name (*_static.mat) for static analysis
        outputfilename = string(owensfilename[1:index-1],"_static.mat")
    elseif (analysisType=="TNB"||analysisType=="TD"||analysisType=="ROM") #output filename (*.mat) for transient analysis
        outputfilename = string(owensfilename[1:index-1],".mat")
    end

    return outputfilename

end

# import PyPlot
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(4.5, 3))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=3.0)
# PyPlot.rc("legend", frameon=false)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]


test_transient = true
test_modal = true
test_flutter = false
println("Starting")

## ****************** COORDINATE SYSTEM DEFINITION *********************** #
# 1 - x -  surge (pointed downstream)
# 2 - y -  sway  (right hand rule)
# 3 - z -  heave (pointed upwards)
# 4 - Ox - roll
# 5 - Oy - pitch
# 6 - Oz - yaw
# *********************************************************************** #

# use this benchmark file
bmOwens = "$path/data/_15mTower_transient_dvawt_c_2_lcdt"
# append this name to the end of the saved files
outFileExt = "_15mTowerExt_NOcentStiff"

# convert rotational stiffness to N-m/rad
StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]
convRotStiff = [1 1 1 180/pi 180/pi 180/pi]
platStiffDiag = StiffDiag_Nm_deg .* convRotStiff

# filename root to save the created nodal file
platfileRoot = "$path/data/1_FourColumnSemi_2ndPass"
platMassDiag = [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9]
platStiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6]

# external dependencies
hydrodynLib = []#"$path/../bin/HydroDyn_c_lib_x64"
moordynLib = []#"$path/../bin/MoorDyn_c_lib_x64"

# define the filename saving convention
fname = string(platfileRoot,outFileExt)

# *********************************************************************
# perform operations for the nodal file generation
# *********************************************************************
MassVal = platMassDiag
StiffVal = platStiffDiag

# nodes = [1 1]
# cmkType = ["M6" "K6"]
# cmkValues = zeros(6,6,2)
# for dd = 1:6
#     # set up mass matrix
#     cmkValues[dd,dd,1] = MassVal[dd]
#     # set up stiffness matrix
#     cmkValues[dd,dd,2] = StiffVal[dd]
# end
# OWENS.writeOwensNDL(fname, nodes, cmkType, cmkValues)

# *********************************************************************
# perform operations for the aerodynamic forces file generation
# *********************************************************************
CACTUSfileRoot = "$path/data/DVAWT_2B_LCDT"
OWENSfileRoot = bmOwens

operatingRPM = 7.2 # rpm
Nrpm = 10    # number of rpm stations
maxRPM = 10
Nmodes = 4  # number of modes to calculate/extract:
timeStep = 2e-3
timeSim = 0.1       # [sec]
n_t = timeSim/timeStep # length of time vector
timeArray = [0, timeSim+1]
rpmArray  = [operatingRPM, operatingRPM]
omegaArrayHz = rpmArray ./ 60
omegaArrayHz2 = [7.1]/60

## ************************************************************************
# perform the transient simulations using OWENS
# *************************************************************************
if test_transient
    println("Running Transient")

    owens(string(fname, ".owens"),"TNB";iterationType="DI", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz, hydrodynLib=hydrodynLib, moordynLib=moordynLib)

    # Perform Tests
    tol = 1e-5
    n_t = 50
    old_filename = "$path/data/UNIT_TEST_15mTower_transient_dvawt_c_2_lcdt.h5"
    new_filename = "$path/data/_15mTower_transient_dvawt_c_2_lcdt.h5"

    if (time() - mtime(new_filename)) >0.70
        println("Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.")
    end

    old_t = HDF5.h5read(old_filename,"t")
    old_aziHist = HDF5.h5read(old_filename,"aziHist")
    old_OmegaHist = HDF5.h5read(old_filename,"OmegaHist")
    old_OmegaDotHist = HDF5.h5read(old_filename,"OmegaDotHist")
    old_gbHist = HDF5.h5read(old_filename,"gbHist")
    old_gbDotHist = HDF5.h5read(old_filename,"gbDotHist")
    old_gbDotDotHist = HDF5.h5read(old_filename,"gbDotDotHist")
    old_FReactionHist = HDF5.h5read(old_filename,"FReactionHist")
    old_genTorque = HDF5.h5read(old_filename,"genTorque")
    old_genPower = HDF5.h5read(old_filename,"genPower")
    old_torqueDriveShaft = HDF5.h5read(old_filename,"torqueDriveShaft")
    old_uHist = HDF5.h5read(old_filename,"uHist")
    old_eps_xx_0_hist = HDF5.h5read(old_filename,"eps_xx_0_hist")
    old_eps_xx_z_hist = HDF5.h5read(old_filename,"eps_xx_z_hist")
    old_eps_xx_y_hist = HDF5.h5read(old_filename,"eps_xx_y_hist")
    old_gam_xz_0_hist = HDF5.h5read(old_filename,"gam_xz_0_hist")
    old_gam_xz_y_hist = HDF5.h5read(old_filename,"gam_xz_y_hist")
    old_gam_xy_0_hist = HDF5.h5read(old_filename,"gam_xy_0_hist")

    t = HDF5.h5read(new_filename,"t")
    aziHist = HDF5.h5read(new_filename,"aziHist")
    OmegaHist = HDF5.h5read(new_filename,"OmegaHist")
    OmegaDotHist = HDF5.h5read(new_filename,"OmegaDotHist")
    gbHist = HDF5.h5read(new_filename,"gbHist")
    gbDotHist = HDF5.h5read(new_filename,"gbDotHist")
    gbDotDotHist = HDF5.h5read(new_filename,"gbDotDotHist")
    FReactionHist = HDF5.h5read(new_filename,"FReactionHist")
    genTorque = HDF5.h5read(new_filename,"genTorque")
    genPower = HDF5.h5read(new_filename,"genPower")
    torqueDriveShaft = HDF5.h5read(new_filename,"torqueDriveShaft")
    uHist = HDF5.h5read(new_filename,"uHist")
    eps_xx_0_hist = HDF5.h5read(new_filename,"epsilon_x_hist")
    eps_xx_z_hist = HDF5.h5read(new_filename,"kappa_y_hist")
    eps_xx_y_hist = HDF5.h5read(new_filename,"kappa_z_hist")
    gam_xz_0_hist = HDF5.h5read(new_filename,"epsilon_z_hist")
    gam_xz_y_hist = HDF5.h5read(new_filename,"kappa_x_hist")
    gam_xy_0_hist = HDF5.h5read(new_filename,"epsilon_y_hist")
    maxT_idx = length(t)
    @test isapprox(old_t[1:maxT_idx],t,atol = tol)
    @test isapprox(old_aziHist[1:maxT_idx],aziHist,atol = tol)
    @test isapprox(old_OmegaHist[1:maxT_idx],OmegaHist,atol = tol)
    @test isapprox(old_OmegaDotHist[1:maxT_idx],OmegaDotHist,atol = tol)
    @test isapprox(old_gbHist[1:maxT_idx],gbHist,atol = tol)
    @test isapprox(old_gbDotHist[1:maxT_idx],gbDotHist,atol = tol)
    @test isapprox(old_gbDotDotHist[1:maxT_idx],gbDotDotHist,atol = tol)
    # for ii = 1:6
    #     PyPlot.figure()
    #     PyPlot.plot(LinRange(0,1,length(old_FReactionHist[:,1])),old_FReactionHist[:,ii],label="old")
    #     PyPlot.plot(LinRange(0,1,length(FReactionHist[:,1])),FReactionHist[:,ii],label="new")
    #     PyPlot.ylabel("Freaction $ii")
    #     PyPlot.legend()
    # end
    for ii = 1:length(FReactionHist[:,1])
        for jj = 1:6#length(FReactionHist[1,:])
            local digits = floor(log10(abs(old_FReactionHist[ii,jj]))) #this way if the tol is 1e-5, then we are actually looking at significant digits, much better than comparing 1e-5 on a 1e6 large number, that's 11 significant digits!
            @test isapprox(old_FReactionHist[ii,jj],FReactionHist[ii,jj],atol=tol*10^(digits+3))
        end
    end
    @test isapprox(old_genTorque[1:maxT_idx],genTorque,atol = tol)
    @test isapprox(old_genPower[1:maxT_idx],genPower,atol = tol)
    @test isapprox(old_torqueDriveShaft[1:maxT_idx],torqueDriveShaft,atol = tol)
    @test isapprox(old_uHist[:,1:maxT_idx],collect(uHist'),atol = 1e-4)
    @test isapprox(old_eps_xx_0_hist[:,:,1:maxT_idx],eps_xx_0_hist,atol = tol)
    @test isapprox(old_eps_xx_z_hist[:,:,1:maxT_idx],eps_xx_z_hist,atol = tol)
    @test isapprox(old_eps_xx_y_hist[:,:,1:maxT_idx],eps_xx_y_hist,atol = tol)
    @test isapprox(old_gam_xz_0_hist[:,:,1:maxT_idx],gam_xz_0_hist,atol = tol)
    @test isapprox(old_gam_xz_y_hist[:,:,1:maxT_idx],gam_xz_y_hist,atol = tol)
    @test isapprox(old_gam_xy_0_hist[:,:,1:maxT_idx],gam_xy_0_hist,atol = tol)
end

# *********************************************************************
# run a modal analysis of the platform design
# *********************************************************************
if test_modal

    owens(string(fname, ".owens"),"M", Omega=0.5*maxRPM*2*pi/60, spinUpOn=false, numModesToExtract=Nmodes)

    # if verify_modal
    old_filename = "$path/data/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out"
    new_filename = "$path/data/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out"

    #Reading function

    numNodes = 82#mesh.numNodes

    freqOLD,dampOLD,U_x_0OLD,U_y_0OLD,U_z_0OLD,theta_x_0OLD,theta_y_0OLD,theta_z_0OLD,U_x_90OLD,U_y_90OLD,U_z_90OLD,theta_x_90OLD,theta_y_90OLD,theta_z_90OLD = OWENS.readResultsModalOut(old_filename,numNodes)
    freq,damp,U_x_0,U_y_0,U_z_0,theta_x_0,theta_y_0,theta_z_0,U_x_90,U_y_90,U_z_90,theta_x_90,theta_y_90,theta_z_90 = OWENS.readResultsModalOut(new_filename,numNodes)

    # tol = 1e-6 # This is covered by the campbell diagram test
    # for imode = 1:length(freq)
    #     used_tol = max(tol*freq[imode],tol) #don't enforce 1e-6 precision on a 1e6 number when we want 6 digits and not 12 digits of precision, also limit it for small number errors
    #     @test isapprox(freqOLD[imode],freq[imode],atol = used_tol)
    #     used_tol = max(tol*damp[imode],tol)
    #     @test isapprox(dampOLD[imode],damp[imode],atol = used_tol)
    # end

    tol = 1e-1
    U_x_0pass = 0
    U_y_0pass = 0
    U_z_0pass = 0
    theta_x_0pass = 0
    theta_y_0pass = 0
    theta_z_0pass = 0
    U_x_90pass = 0
    U_y_90pass = 0
    U_z_90pass = 0
    theta_x_90pass = 0
    theta_y_90pass = 0
    theta_z_90pass = 0

    for imode = 1:length(U_x_0OLD[1,:])
        # println("mode: $imode")
        # for inode = 1:numNodes
        # println("node $inode")
        for inode = 1:length(U_x_0OLD[:,1])
            if isapprox(abs.(U_x_0OLD[inode,imode]),abs.(U_x_0[inode,imode]),atol = tol)
                global U_x_0pass += 1
            end
            if isapprox(abs.(U_y_0OLD[inode,imode]),abs.(U_y_0[inode,imode]),atol = tol)
                global U_y_0pass += 1
            end
            if isapprox(abs.(U_z_0OLD[inode,imode]),abs.(U_z_0[inode,imode]),atol = tol)
                global U_z_0pass += 1
            end
            if isapprox(abs.(theta_x_0OLD[inode,imode]),abs.(theta_x_0[inode,imode]),atol = tol)
                global theta_x_0pass += 1
            end
            if isapprox(abs.(theta_y_0OLD[inode,imode]),abs.(theta_y_0[inode,imode]),atol = tol)
                global theta_y_0pass += 1
            end
            if isapprox(abs.(theta_z_0OLD[inode,imode]),abs.(theta_z_0[inode,imode]),atol = tol)
                global theta_z_0pass += 1
            end
            if isapprox(abs.(U_x_90OLD[inode,imode]),abs.(U_x_90[inode,imode]),atol = tol)
                global U_x_90pass += 1
            end
            if isapprox(abs.(U_y_90OLD[inode,imode]),abs.(U_y_90[inode,imode]),atol = tol)
                global U_y_90pass += 1
            end
            if isapprox(abs.(U_z_90OLD[inode,imode]),abs.(U_z_90[inode,imode]),atol = tol)
                global U_z_90pass += 1
            end
            if isapprox(abs.(theta_x_90OLD[inode,imode]),abs.(theta_x_90[inode,imode]),atol = tol)
                global theta_x_90pass += 1
            end
            if isapprox(abs.(theta_y_90OLD[inode,imode]),abs.(theta_y_90[inode,imode]),atol = tol)
                global theta_y_90pass += 1
            end
            if isapprox(abs.(theta_z_90OLD[inode,imode]),abs.(theta_z_90[inode,imode]),atol = tol)
                global theta_z_90pass += 1
            end
        end
    end

    # at least 80 percent of the modeshapes are identical indicates (despite the recripocity of the solutions) that the analysis is adequate
    tol2 = 0.8
    @test U_x_0pass/length(U_x_0OLD)>tol2
    @test U_y_0pass/length(U_x_0OLD)>tol2
    @test U_z_0pass/length(U_x_0OLD)>tol2
    @test theta_x_0pass/length(U_x_0OLD)>tol2
    @test theta_y_0pass/length(U_x_0OLD)>tol2
    @test theta_z_0pass/length(U_x_0OLD)>tol2
    @test U_x_90pass/length(U_x_0OLD)>tol2
    @test U_y_90pass/length(U_x_0OLD)>tol2
    @test U_z_90pass/length(U_x_0OLD)>tol2
    @test theta_x_90pass/length(U_x_0OLD)>tol2
    @test theta_y_90pass/length(U_x_0OLD)>tol2
    @test theta_z_90pass/length(U_x_0OLD)>tol2
    # end

end

if test_flutter
    # freq,damp=owens(string(fname ".owens"),"FA", omegaArrayHz2, true, 1.2041, Nmodes)
end

println("Function Finished")
