global hdlib
global sym_calcoutput
global sym_updatestates
global sym_end
global backup_Vx

path,_ = splitdir(@__FILE__)

mutable struct HD_Error
error_status
error_message
end

function HD_Init(hdlib_filename, output_root_name; hd_input_file="none", WtrDens=1025, WtrDpth=200, MSL2SWL=0,
    WaveMod=2, WaveStMod=0, WaveHs=2.0, WaveTp=6.0, WavePkShp=1.0, WaveDir=0.0, WaveSeed=123456789, gravity = 9.81,
    PotFile="$path/../test/tlpmit",
    CurrMod=0, CurrSSV0=0, CurrSSDir="DEFAULT", CurrNSRef=20, CurrNSV0=0, CurrNSDir=0, CurrDIV=0, CurrDIDir=0,
    ptfm_ref_pos_x=0.0, ptfm_ref_pos_y=0.0, num_node_pts=1,
    init_node_pos=zeros(6), interp_order=1, t_initial=0.0, dt=0.01, t_max=60.0)

    global hd_abort_error_level = 4


    if hd_input_file == "none"
        # Where the input is manipulated
        WtrDens_str = "      $WtrDens    WtrDens        - Water density (kg/m^3)"
        WtrDpth_str = "      $WtrDpth    WtrDpth        - Water depth (meters)"
        MSL2SWL_str = "      $MSL2SWL    MSL2SWL        - Offset between still-water level and mean sea level (meters) [positive upward; unused when WaveMod = 6; must be zero if PotMod=1 or 2]"
        WaveMod_str = "      $WaveMod    WaveMod        - Incident wave kinematics model {0: none=still water, 1: regular (periodic), 1P#: regular with user-specified phase, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: White noise spectrum (irregular), 4: user-defined spectrum from routine UserWaveSpctrm (irregular), 5: Externally generated wave-elevation time series, 6: Externally generated full wave-kinematics time series [option 6 is invalid for PotMod/=0]} (switch)"
        WaveStMod_str = "  $WaveStMod    WaveStMod      - Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching} (switch) [unused when WaveMod=0 or when PotMod/=0]"
        WaveHs_str = "       $WaveHs     WaveHs         - Significant wave height of incident waves (meters) [used only when WaveMod=1, 2, or 3]"
        WaveTp_str = "       $WaveTp     WaveTp         - Peak-spectral period of incident waves       (sec) [used only when WaveMod=1 or 2]"
        WavePkShp_str = "  $WavePkShp    WavePkShp      - Peak-shape parameter of incident wave spectrum (-) or DEFAULT (string) [used only when WaveMod=2; use 1.0 for Pierson-Moskowitz]"
        WaveDir_str = "      $WaveDir    WaveDir        - Incident wave propagation heading direction                         (degrees) [unused when WaveMod=0 or 6]"
        WaveSeed1_str = "   $WaveSeed    WaveSeed(1)    - First  random seed of incident waves [-2147483648 to 2147483647]    (-)       [unused when WaveMod=0, 5, or 6]"
        PotFile_str = "  \"$PotFile\"    PotFile       - Root name of potential-flow model data; WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3) (quoted string) [1 to NBody if NBodyMod>1] [MAKE SURE THE FREQUENCIES INHERENT IN THESE WAMIT FILES SPAN THE PHYSICALLY-SIGNIFICANT RANGE OF FREQUENCIES FOR THE GIVEN PLATFORM; THEY MUST CONTAIN THE ZERO- AND INFINITE-FREQUENCY LIMITS!]"

        input_string_array = [
            "------- HydroDyn v2.03.* Input File --------------------------------------------",
            "NREL 5.0 MW offshore baseline floating platform HydroDyn input properties for the TLP.",
            "   False         Echo           - Echo the input file data (flag)",
        "---------------------- ENVIRONMENTAL CONDITIONS --------------------------------",
                WtrDens_str,
                WtrDpth_str,
                MSL2SWL_str,
        "---------------------- WAVES ---------------------------------------------------",
                WaveMod_str,
                WaveStMod_str,
                "     4600   WaveTMax       - Analysis time for incident wave calculations (sec) [unused when WaveMod=0; determines WaveDOmega=2Pi/WaveTMax in the IFFT]",
                "      0.2   WaveDT         - Time step for incident wave calculations     (sec) [unused when WaveMod=0; 0.1<=WaveDT<=1.0 recommended; determines WaveOmegaMax=Pi/WaveDT in the IFFT]",
                WaveHs_str,
                WaveTp_str,
                WavePkShp_str,
                " 0.314159   WvLowCOff      - Low  cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]",
                " 1.570796   WvHiCOff       - High cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s) [unused when WaveMod=0, 1, or 6]",
                WaveDir_str,
                "        0   WaveDirMod     - Directional spreading function {0: none, 1: COS2S}                  (-)       [only used when WaveMod=2,3, or 4]",
                "        1   WaveDirSpread  - Wave direction spreading coefficient ( > 0 )                        (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1]",
                "        1   WaveNDir       - Number of wave directions                                           (-)       [only used when WaveMod=2,3, or 4 and WaveDirMod=1; odd number only]",
                "        0   WaveDirRange   - Range of wave directions (full range: WaveDir +/- 1/2*WaveDirRange) (degrees) [only used when WaveMod=2,3,or 4 and WaveDirMod=1]",
                WaveSeed1_str,
            "   \"RANLUX\"   WaveSeed(2)    - Second random seed of incident waves [-2147483648 to 2147483647] for intrinsic pRNG, or an alternative pRNG: \"RanLux\"    (-)       [unused when WaveMod=0, 5, or 6]",
            "    FALSE        WaveNDAmp      - Flag for normally distributed amplitudes                            (flag)    [only used when WaveMod=2, 3, or 4]",
            "   \"\"         WvKinFile      - Root name of externally generated wave data file(s)        (quoted string)    [used only when WaveMod=5 or 6]",
                "        1   NWaveElev      - Number of points where the incident wave elevations can be computed (-)       [maximum of 9 output locations]",
                "        0   WaveElevxi     - List of xi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]",
                "        0   WaveElevyi     - List of yi-coordinates for points where the incident wave elevations can be output (meters) [NWaveElev points, separated by commas or white space; usused if NWaveElev = 0]",
        "---------------------- 2ND-ORDER WAVES ----------------------------------------- [unused with WaveMod=0 or 6]",
            "   False        WvDiffQTF      - Full difference-frequency 2nd-order wave kinematics (flag)",
            "   False        WvSumQTF       - Full summation-frequency  2nd-order wave kinematics (flag)",
                "        0   WvLowCOffD     - Low  frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]",
                " 1.256637   WvHiCOffD      - High frequency cutoff used in the difference-frequencies (rad/s) [Only used with a difference-frequency method]",
                " 0.618319   WvLowCOffS     - Low  frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]",
                " 3.141593   WvHiCOffS      - High frequency cutoff used in the summation-frequencies  (rad/s) [Only used with a summation-frequency  method]",
        "---------------------- CURRENT ------------------------------------------------- [unused with WaveMod=6]",
                "       0    CurrMod        - Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent} (switch)",
                "       0    CurrSSV0       - Sub-surface current velocity at still water level  (m/s) [used only when CurrMod=1]",
                "\"DEFAULT\"  CurrSSDir      - Sub-surface current heading direction (degrees) or DEFAULT (string) [used only when CurrMod=1]",
                "       20   CurrNSRef      - Near-surface current reference depth            (meters) [used only when CurrMod=1]",
                "       0    CurrNSV0       - Near-surface current velocity at still water level (m/s) [used only when CurrMod=1]",
                "       0    CurrNSDir      - Near-surface current heading direction         (degrees) [used only when CurrMod=1]",
                "       0    CurrDIV        - Depth-independent current velocity                 (m/s) [used only when CurrMod=1]",
        "---------------------- FLOATING PLATFORM --------------------------------------- [unused with WaveMod=6]",
                "        1   PotMod         - Potential-flow model {0: none=no potential flow, 1: frequency-to-time-domain transforms based on WAMIT output, 2: fluid-impulse theory (FIT)} (switch)",
                "        1   ExctnMod       - Wave-excitation model {0: no wave-excitation calculation, 1: DFT, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ssexctn INPUT FILE]",
                "        1   RdtnMod        - Radiation memory-effect model {0: no memory-effect calculation, 1: convolution, 2: state-space} (switch) [only used when PotMod=1; STATE-SPACE REQUIRES *.ss INPUT FILE]",
                "       60   RdtnTMax       - Analysis time for wave radiation kernel calculations (sec) [only used when PotMod=1 and RdtnMod>0; determines RdtnDOmega=Pi/RdtnTMax in the cosine transform; MAKE SURE THIS IS LONG ENOUGH FOR THE RADIATION IMPULSE RESPONSE FUNCTIONS TO DECAY TO NEAR-ZERO FOR THE GIVEN PLATFORM!]",
        " \"default\"        RdtnDT         - Time step for wave radiation kernel calculations (sec) [only used when PotMod=1 and ExctnMod>0 or RdtnMod>0; DT<=RdtnDT<=0.1 recommended; determines RdtnOmegaMax=Pi/RdtnDT in the cosine transform]",
                "        1   NBody          - Number of WAMIT bodies to be used (-) [>=1; only used when PotMod=1. If NBodyMod=1, the WAMIT data contains a vector of size 6*NBody x 1 and matrices of size 6*NBody x 6*NBody; if NBodyMod>1, there are NBody sets of WAMIT data each with a vector of size 6 x 1 and matrices of size 6 x 6]",
                "        1   NBodyMod       - Body coupling model {1: include coupling terms between each body and NBody in HydroDyn equals NBODY in WAMIT, 2: neglect coupling terms between each body and NBODY=1 with XBODY=0 in WAMIT, 3: Neglect coupling terms between each body and NBODY=1 with XBODY=/0 in WAMIT} (switch) [only used when PotMod=1]",
                PotFile_str,
                "        1   WAMITULEN      - Characteristic body length scale used to redimensionalize WAMIT output (meters) [1 to NBody if NBodyMod>1] [only used when PotMod=1]",
                "      0.0   PtfmRefxt      - The xt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]",
                "      0.0   PtfmRefyt      - The yt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1]",
                "      0.0   PtfmRefzt      - The zt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] [only used when PotMod=1. If NBodyMod=2,PtfmRefzt=0.0]",
                "      0.0   PtfmRefztRot   - The rotation about zt of the body reference frame(s) from xt/yt (degrees) [1 to NBody] [only used when PotMod=1]",
                "    13917   PtfmVol0       - Displaced volume of water when the body is in its undisplaced position (m^3) [1 to NBody] [only used when PotMod=1; USE THE SAME VALUE COMPUTED BY WAMIT AS OUTPUT IN THE .OUT FILE!]",
                "      0.0   PtfmCOBxt      - The xt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]",
                "      0.0   PtfmCOByt      - The yt offset of the center of buoyancy (COB) from (0,0) (meters) [1 to NBody] [only used when PotMod=1]",
    "---------------------- 2ND-ORDER FLOATING PLATFORM FORCES ---------------------- [unused with WaveMod=0 or 6, or PotMod=0 or 2]",
                "        0   MnDrift        - Mean-drift 2nd-order forces computed                                       {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. If NBody>1, MnDrift  /=8]",
                "        0   NewmanApp      - Mean- and slow-drift 2nd-order forces computed with Newman's approximation {0: None; [7, 8, 9, 10, 11, or 12]: WAMIT file to use} [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero. If NBody>1, NewmanApp/=8. Used only when WaveDirMod=0]",
                "       12   DiffQTF        - Full difference-frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}          [Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero]",
                "       12   SumQTF         - Full summation -frequency 2nd-order forces computed with full QTF          {0: None; [10, 11, or 12]: WAMIT file to use}",
    "---------------------- PLATFORM ADDITIONAL STIFFNESS AND DAMPING  -------------- [unused with PotMod=0 or 2]",
                "        0   AddF0    - Additional preload (N, N-m) [If NBodyMod=1, one size 6*NBody x 1 vector; if NBodyMod>1, NBody size 6 x 1 vectors]",
                "        0",
                "        0",
                "        0",
                "        0",
                "        0",
                "        0             0             0             0             0             0   AddCLin  - Additional linear stiffness (N/m, N/rad, N-m/m, N-m/rad)                     [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0   AddBLin  - Additional linear damping(N/(m/s), N/(rad/s), N-m/(m/s), N-m/(rad/s))        [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0   AddBQuad - Additional quadratic drag(N/(m/s)^2, N/(rad/s)^2, N-m(m/s)^2, N-m/(rad/s)^2) [If NBodyMod=1, one size 6*NBody x 6*NBody matrix; if NBodyMod>1, NBody size 6 x 6 matrices]",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
                "        0             0             0             0             0             0",
    "---------------------- AXIAL COEFFICIENTS --------------------------------------",
                "        2   NAxCoef        - Number of axial coefficients (-)",
    "AxCoefID  AxCd     AxCa     AxCp",
        "(-)    (-)      (-)      (-)",
    "      1     0.00     0.00     1.00",
    "      2     9.60     0.00     1.00",
    "---------------------- MEMBER JOINTS -------------------------------------------",
                "       44   NJoints        - Number of joints (-)   [must be exactly 0 or at least 2]",
    "JointID   Jointxi     Jointyi     Jointzi  JointAxID   JointOvrlp   [JointOvrlp= 0: do nothing at joint, 1: eliminate overlaps by calculating super member]",
    "    (-)     (m)         (m)         (m)        (-)       (switch)",
    " 1     0.00000     0.00000   -20.00000      1            0",
    " 2     0.00000     0.00000    10.00000      1            0",
    " 3    14.43376    25.00000   -14.00000      1            0",
    " 4    14.43376    25.00000    12.00000      1            0",
    " 5   -28.86751     0.00000   -14.00000      1            0",
    " 6   -28.86751     0.00000    12.00000      1            0",
    " 7    14.43376   -25.00000   -14.00000      1            0",
    " 8    14.43376   -25.00000    12.00000      1            0",
    " 9    14.43375    25.00000   -20.00000      2            0",
    "10   -28.86750     0.00000   -20.00000      2            0",
    "11    14.43375   -25.00000   -20.00000      2            0",
    "12     9.23760    22.00000    10.00000      1            0",
    "13   -23.67130     3.00000    10.00000      1            0",
    "14   -23.67130    -3.00000    10.00000      1            0",
    "15     9.23760   -22.00000    10.00000      1            0",
    "16    14.43375   -19.00000    10.00000      1            0",
    "17    14.43375    19.00000    10.00000      1            0",
    "18     4.04145    19.00000   -17.00000      1            0",
    "19   -18.47520     6.00000   -17.00000      1            0",
    "20   -18.47520    -6.00000   -17.00000      1            0",
    "21     4.04145   -19.00000   -17.00000      1            0",
    "22    14.43375   -13.00000   -17.00000      1            0",
    "23    14.43375    13.00000   -17.00000      1            0",
    "24     1.62500     2.81500    10.00000      1            0",
    "25    11.43376    19.80385    10.00000      1            0",
    "26    -3.25000     0.00000    10.00000      1            0",
    "27   -22.87000     0.00000    10.00000      1            0",
    "28     1.62500    -2.81500    10.00000      1            0",
    "29    11.43376   -19.80385    10.00000      1            0",
    "30     1.62500     2.81500   -17.00000      1            0",
    "31     8.43376    14.60770   -17.00000      1            0",
    "32    -3.25000     0.00000   -17.00000      1            0",
    "33   -16.87000     0.00000   -17.00000      1            0",
    "34     1.62500    -2.81500   -17.00000      1            0",
    "35     8.43376   -14.60770   -17.00000      1            0",
    "36     1.62500     2.81500   -16.20000      1            0",
    "37    11.43376    19.80385     9.13000      1            0",
    "38    -3.25000     0.00000   -16.20000      1            0",
    "39   -22.87000     0.00000     9.13000      1            0",
    "40     1.62500    -2.81500   -16.20000      1            0",
    "41    11.43376   -19.80385     9.13000      1            0",
    "42    14.43376    25.00000   -19.94000      1            0",
    "43   -28.86751     0.00000   -19.94000      1            0",
    "44    14.43376   -25.00000   -19.94000      1            0",
    "---------------------- MEMBER CROSS-SECTION PROPERTIES -------------------------",
                "        4   NPropSets      - Number of member property sets (-)",
    "PropSetID    PropD         PropThck",
    "     (-)        (m)            (m)",
    "1        6.50000        0.03000          ! Main Column",
    "2       12.00000        0.06000          ! Upper Columns",
    "3       24.00000        0.06000          ! Base Columns",
    "4        1.60000        0.01750          ! Pontoons",
    "---------------------- SIMPLE HYDRODYNAMIC COEFFICIENTS (model 1) --------------",
    "SimplCd    SimplCdMG    SimplCa    SimplCaMG    SimplCp    SimplCpMG   SimplAxCd  SimplAxCdMG  SimplAxCa  SimplAxCaMG  SimplAxCp   SimplAxCpMG",
    "    (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)         (-)",
    "   0.00        0.00        0.00        0.00        1.00        1.00        0.00        0.00        0.00        0.00        1.00        1.00",
    "---------------------- DEPTH-BASED HYDRODYNAMIC COEFFICIENTS (model 2) ---------",
                "       0   NCoefDpth       - Number of depth-dependent coefficients (-)",
    "Dpth      DpthCd   DpthCdMG   DpthCa   DpthCaMG       DpthCp   DpthCpMG   DpthAxCd   DpthAxCdMG   DpthAxCa   DpthAxCaMG   DpthAxCp   DpthAxCpMG",
    "   (m)       (-)      (-)        (-)      (-)            (-)      (-)        (-)        (-)          (-)        (-)          (-)        (-)",
    "---------------------- MEMBER-BASED HYDRODYNAMIC COEFFICIENTS (model 3) --------",
                "      25   NCoefMembers       - Number of member-based coefficients (-)",
    "MemberID    MemberCd1     MemberCd2    MemberCdMG1   MemberCdMG2    MemberCa1     MemberCa2    MemberCaMG1   MemberCaMG2    MemberCp1     MemberCp2    MemberCpMG1   MemberCpMG2   MemberAxCd1   MemberAxCd2  MemberAxCdMG1 MemberAxCdMG2  MemberAxCa1   MemberAxCa2  MemberAxCaMG1 MemberAxCaMG2  MemberAxCp1  MemberAxCp2   MemberAxCpMG1   MemberAxCpMG2",
    "     (-)         (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)           (-)",
    "    1          0.56          0.56          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Main Column",
    "    2          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 1",
    "    3          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 2",
    "    4          0.61          0.61          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Upper Column 3",
    "    5          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 1",
    "    6          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 2",
    "    7          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base Column 3",
    "   23          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 1",
    "   24          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 2",
    "   25          0.68          0.68          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Base column cap 3",
    "    8          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 1",
    "    9          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 2",
    "   10          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Upper 3",
    "   11          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 1",
    "   12          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 2",
    "   13          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Delta Pontoon, Lower 3",
    "   14          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 1",
    "   15          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 2",
    "   16          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Upper 3",
    "   17          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 1",
    "   18          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 2",
    "   19          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Y Pontoon, Lower 3",
    "   20          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 1",
    "   21          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 2",
    "   22          0.63          0.63          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00          0.00         ! Cross Brace 3",
    "-------------------- MEMBERS -------------------------------------------------",
                "      25   NMembers       - Number of members (-)",
    "MemberID  MJointID1  MJointID2  MPropSetID1  MPropSetID2  MDivSize   MCoefMod  PropPot   [MCoefMod=1: use simple coeff table, 2: use depth-based coeff table, 3: use member-based coeff table] [ PropPot/=0 if member is modeled with potential-flow theory]",
    "    (-)        (-)        (-)         (-)          (-)        (m)      (switch)   (flag)",
    "     1         1          2           1            1         1.0000      3        TRUE           ! Main Column",
    "     2         3          4           2            2         1.0000      3        TRUE           ! Upper Column 1",
    "     3         5          6           2            2         1.0000      3        TRUE           ! Upper Column 2",
    "     4         7          8           2            2         1.0000      3        TRUE           ! Upper Column 3",
    "     5        42          3           3            3         1.0000      3        TRUE           ! Base Column 1",
    "     6        43          5           3            3         1.0000      3        TRUE           ! Base Column 2",
    "     7        44          7           3            3         1.0000      3        TRUE           ! Base Column 3",
    "    23         9         42           3            3         1.0000      3        TRUE           ! Base column cap 1",
    "    24        10         43           3            3         1.0000      3        TRUE           ! Base column cap 2",
    "    25        11         44           3            3         1.0000      3        TRUE           ! Base column cap 3",
    "     8        12         13           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 1",
    "     9        14         15           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 2",
    "    10        16         17           4            4         1.0000      3        TRUE           ! Delta Pontoon, Upper 3",
    "    11        18         19           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 1",
    "    12        20         21           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 2",
    "    13        22         23           4            4         1.0000      3        TRUE           ! Delta Pontoon, Lower 3",
    "    14        24         25           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 1",
    "    15        26         27           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 2",
    "    16        28         29           4            4         1.0000      3        TRUE           ! Y Pontoon, Upper 3",
    "    17        30         31           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 1",
    "    18        32         33           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 2",
    "    19        34         35           4            4         1.0000      3        TRUE           ! Y Pontoon, Lower 3",
    "    20        36         37           4            4         1.0000      3        TRUE           ! Cross Brace 1",
    "    21        38         39           4            4         1.0000      3        TRUE           ! Cross Brace 2",
    "    22        40         41           4            4         1.0000      3        TRUE           ! Cross Brace 3",
    "---------------------- FILLED MEMBERS ------------------------------------------",
                "       2   NFillGroups     - Number of filled member groups (-) [If FillDens = DEFAULT, then FillDens = WtrDens; FillFSLoc is related to MSL2SWL]",
    "FillNumM FillMList             FillFSLoc     FillDens",
    "    (-)      (-)                   (m)           (kg/m^3)",
    "     3   2   3   4    -6.17           1025",
    "     3   5   6   7   -14.89           1025",
    "---------------------- MARINE GROWTH -------------------------------------------",
                "       0   NMGDepths      - Number of marine-growth depths specified (-)",
    "MGDpth     MGThck       MGDens",
    "   (m)        (m)         (kg/m^3)",
    "---------------------- MEMBER OUTPUT LIST --------------------------------------",
                "       0   NMOutputs      - Number of member outputs (-) [must be < 10]",
    "MemberID   NOutLoc    NodeLocs [NOutLoc < 10; node locations are normalized distance from the start of the member, and must be >=0 and <= 1] [unused if NMOutputs=0]",
    "     (-)        (-)        (-)",
    "---------------------- JOINT OUTPUT LIST ---------------------------------------",
                "       0   NJOutputs      - Number of joint outputs [Must be < 10]",
                "       0   JOutLst        - List of JointIDs which are to be output (-)[unused if NJOutputs=0]",
    "---------------------- OUTPUT --------------------------------------------------",
        "  True             HDSum          - Output a summary file [flag]",
        "  False            OutAll         - Output all user-specified member and joint loads (only at each member end, not interior locations) [flag]",
        "               3   OutSwtch       - Output requested channels to: [1=Hydrodyn.out, 2=GlueCode.out, 3=both files]",
        " \"E15.7e2\"       OutFmt         - Output format for numerical results (quoted string) [not checked for validity]",
        "\"A11\"            OutSFmt        - Output format for header strings (quoted string) [not checked for validity]",
    "---------------------- OUTPUT CHANNELS -----------------------------------------",
        "\"Wave1Elev\"               - Wave elevation at the platform reference point (  0,  0)",
        "\"HydroFxi\"",
        "\"HydroFyi\"",
        "\"HydroFzi\"",
        "\"HydroMxi\"",
        "\"HydroMyi\"",
        "\"HydroMzi\"",
        "\"B1Surge\"",
        "\"B1Sway\"",
        "\"B1Heave\"",
        "\"B1Roll\"",
        "\"B1Pitch\"",
        "\"B1Yaw\"",
        "\"B1TVxi\"",
        "\"B1TVyi\"",
        "\"B1TVzi\"",
        "\"B1RVxi\"",
        "\"B1RVyi\"",
        "\"B1RVzi\"",
        "\"B1TAxi\"",
        "\"B1TAyi\"",
        "\"B1TAzi\"",
        "\"B1RAxi\"",
        "\"B1RAyi\"",
        "\"B1RAzi\"",
        "\"B1WvsFxi\"",
        "\"B1WvsFyi\"",
        "\"B1WvsFzi\"",
        "\"B1WvsMxi\"",
        "\"B1WvsMyi\"",
        "\"B1WvsMzi\"",
        "\"B1HDSFxi\"",
        "\"B1HDSFyi\"",
        "\"B1HDSFzi\"",
        "\"B1HDSMxi\"",
        "\"B1HDSMyi\"",
        "\"B1HDSMzi\"",
        "\"B1RdtFxi\"",
        "\"B1RdtFyi\"",
        "\"B1RdtFzi\"",
        "\"B1RdtMxi\"",
        "\"B1RdtMyi\"",
        "\"B1RdtMzi\"",
        "END of output channels and end of file. (the word 'END' must appear in the first 3 columns of this line)",
    ]

    else
        println("Reading HydroDyn data from $hd_input_file.")
        fid = open(hd_input_file, "r") 
        input_string_array = readlines(fid)
        close(fid)
    end

    input_string        = join(input_string_array, "\0")
    input_string_length = length(input_string)

    # Allocate Outputs
    num_channels = [0]
    channel_names = string(repeat(" ", 20 * 4000))
    channel_units = string(repeat(" ", 20 * 4000))

    global hdlib = Libdl.dlopen(hdlib_filename) # Open the library explicitly.
    global hd_active = true

    hd_sym_init = Libdl.dlsym(hdlib, :HydroDyn_C_Init)   # Get a symbol for the function to call.
    global hd_sym_calcoutput = Libdl.dlsym(hdlib, :HydroDyn_C_CalcOutput)   # Get a symbol for the function to call.
    global hd_sym_updatestates = Libdl.dlsym(hdlib, :HydroDyn_C_UpdateStates)
    global hd_sym_end = Libdl.dlsym(hdlib, :HydroDyn_C_End) # !!! "c" is capitalized in library, change if errors
    global hd_err = HD_Error([0], string(repeat(" ", 1025)))
    
    ccall(hd_sym_init,Cint,
        (Cstring,           # IN: output_root_name
        Ptr{Ptr{Cchar}},    # IN: input_string_array
        Ref{Cint},          # IN: input_string_length
        Ref{Cfloat},        # IN: gravity
        Ref{Cfloat},        # IN: WtrDens
        Ref{Cfloat},        # IN: WtrDpth
        Ref{Cfloat},        # IN: MSL2SWL
        Ref{Cfloat},        # IN: ptfm_ref_pos_x
        Ref{Cfloat},        # IN: ptfm_ref_pos_y
        Ref{Cint},          # IN: num_node_pts
        Ref{Cfloat},        # IN: init_node_pos
        Ref{Cint},          # IN: interp_order
        Ref{Cdouble},       # IN: t_initial
        Ref{Cdouble},       # IN: dt
        Ref{Cdouble},       # IN: t_max
        Ptr{Cint},          # OUT: num_channels
        Cstring,            # OUT: channel_names
        Cstring,            # OUT: channel_units
        Ptr{Cint},          # OUT: error_status
        Cstring),           # OUT: error_message
        output_root_name,
        [input_string],
        input_string_length,
        gravity,
        WtrDens,
        WtrDpth,
        MSL2SWL,
        ptfm_ref_pos_x,
        ptfm_ref_pos_y,
        num_node_pts,
        Cfloat.(init_node_pos),
        interp_order,
        t_initial,
        dt,
        t_max,
        num_channels,
        channel_names,
        channel_units,
        hd_err.error_status,
        hd_err.error_message)

    hd_check_error()

end

function HD_CalcOutput(time, node_pos, node_vel, node_acc, node_force, out_channel_vals; num_node_pts=1)

    # error_message = string(repeat(" ", 1025))
    # error_status = [0]

    if hd_active
        
        ccall(hd_sym_calcoutput,Cint,
            (Ptr{Cdouble},      # IN: time
            Ref{Cint},          # IN: num_node_pts
            Ref{Cfloat},        # IN: node_pos
            Ref{Cfloat},        # IN: node_vel
            Ref{Cfloat},        # IN: node_acc
            Ptr{Cfloat},        # OUT: node_force
            Ptr{Cfloat},        # OUT: out_channel_vals
            Ptr{Cint},          # OUT: error_status
            Cstring),      # OUT: error_message 
            [time],
            num_node_pts,
            Cfloat.(node_pos),
            Cfloat.(node_vel),
            Cfloat.(node_acc),
            node_force,
            out_channel_vals,
            hd_err.error_status,
            hd_err.error_message) 

        hd_check_error()

    else
        error("HydroDyn instance has not been initialized. Use HD_Init() function.")
    end

    return node_force, out_channel_vals
end

function HD_UpdateStates(time, next_time, node_pos, node_vel, node_acc; num_node_pts=1)

    if hd_active

        ccall(hd_sym_updatestates,Cint,
            (Ptr{Cdouble},      # IN: time
            Ptr{Cdouble},       # IN: next_time
            Ref{Cint},          # IN: num_node_pts
            Ref{Cfloat},        # IN: node_pos
            Ref{Cfloat},        # IN: node_vel
            Ref{Cfloat},        # IN: node_acc
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message 
            [time],
            [next_time],
            num_node_pts,
            Cfloat.(node_pos),
            Cfloat.(node_vel),
            Cfloat.(node_acc),
            hd_err.error_status,
            hd_err.error_message) 

        hd_check_error()

    else
        error("HydroDyn instance has not been initialized. Use HD_Init() function.")
    end
end

function HD_End()

    if hd_active

        global hd_active = false
        ccall(hd_sym_end,Cint,
        (Ptr{Cint},         # OUT: ErrStat_C
        Cstring),           # OUT: ErrMsg_C
        hd_err.error_status,
        hd_err.error_message)

        Libdl.dlclose(hdlib) # Close the library explicitly.

    end
end

function hd_check_error()
    if hd_err.error_status[1] == 0
        hd_err.error_status = [0] # reset error status/message
        hd_err.error_message = string(repeat(" ", 1025))
    elseif hd_err.error_status[1] < hd_abort_error_level
        @warn("Error status " * string(hd_err.error_status[1]) * ": " * string(hd_err.error_message))
        hd_err.error_status = [0] # reset error status/message
        hd_err.error_message = string(repeat(" ", 1025))
    else
        @warn("Error status " * string(hd_err.error_status[1]) * ": " * string(hd_err.error_message))
        hd_err.error_status = [0] # reset error status/message
        hd_err.error_message = string(repeat(" ", 1025))
        HD_End()
        error("HydroDyn terminated prematurely.")
    end
end
