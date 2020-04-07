//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: transientExec.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//

// Include Files
#include "transientExec.h"
#include "externalForcing.h"
#include "eye.h"
#include "fileManager.h"
#include "fwrite.h"
#include "initialElementCalculations.h"
#include "norm.h"
#include "processAeroLoadsBLE.h"
#include "rt_nonfinite.h"
#include "sprintf.h"
#include "strcmp.h"
#include "structuralDynamicsTransient.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include "tic.h"
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <string.h>

// Function Declarations
static void omegaSpecCheck(double tCurrent, const double tocp[2], double
  *OmegaCurrent, double *OmegaDotCurrent, boolean_T *terminateSimulation);

// Function Definitions

//
// Arguments    : double tCurrent
//                const double tocp[2]
//                double *OmegaCurrent
//                double *OmegaDotCurrent
//                boolean_T *terminateSimulation
// Return Type  : void
//
static void omegaSpecCheck(double tCurrent, const double tocp[2], double
  *OmegaCurrent, double *OmegaDotCurrent, boolean_T *terminateSimulation)
{
  int k;
  int exitg1;
  double Vq;
  double b_Vq;
  if (1.1 < tCurrent) {
    *terminateSimulation = true;
    *OmegaCurrent = 0.0;
    *OmegaDotCurrent = 0.0;
  } else {
    k = 0;
    do {
      exitg1 = 0;
      if (k < 2) {
        if (rtIsNaN(tocp[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        *OmegaCurrent = rtNaN;
        if (!rtIsNaN(tCurrent)) {
          *OmegaCurrent = 0.12000000000000001;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);

    // interpolated discreteized profile for current omega
    // calculate current rotor acceleration
    k = 0;
    do {
      exitg1 = 0;
      if (k < 2) {
        if (rtIsNaN(tocp[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        Vq = rtNaN;
        if (!rtIsNaN(tCurrent - 0.001)) {
          Vq = 0.12000000000000001;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);

    k = 0;
    do {
      exitg1 = 0;
      if (k < 2) {
        if (rtIsNaN(tocp[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_Vq = rtNaN;
        if ((!rtIsNaN(tCurrent + 0.001)) && (!(tCurrent + 0.001 > 1.1))) {
          b_Vq = 0.12000000000000001;
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);

    *OmegaDotCurrent = (b_Vq - Vq) / 0.002;
    *terminateSimulation = false;
  }
}

//
// transientExec performs modular transient analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//
//    transientExec(model,mesh,el)
//
//    %This function is an executable function for transient analysis. It
//    provides the interface of various external module with transient
//    structural dynamics analysis capability.
//
//    input:
//    model       = object containing model data
//    mesh        = object containing mesh data
//    el          = object containing element data
//
//
//    output: (NONE)
// Arguments    : const char model_analysisType[3]
//                const double model_tocp[2]
//                const emxArray_char_T *model_aeroloadfile
//                const emxArray_char_T *model_owensfile
//                double model_RayleighAlpha
//                double model_RayleighBeta
//                double model_BC_numpBC
//                const emxArray_real_T *model_BC_pBC
//                const emxArray_real_T *model_joint
//                const char model_outFilename_data[]
//                const int model_outFilename_size[2]
//                const emxArray_real_T *model_jointTransform
//                const g_struct_T mesh
//                const h_struct_T el
// Return Type  : void
//
void transientExec(const char model_analysisType[3], const double model_tocp[2],
                   const emxArray_char_T *model_aeroloadfile, const
                   emxArray_char_T *model_owensfile, double model_RayleighAlpha,
                   double model_RayleighBeta, double model_BC_numpBC, const
                   emxArray_real_T *model_BC_pBC, const emxArray_real_T
                   *model_joint, const char model_outFilename_data[], const int
                   model_outFilename_size[2], const emxArray_real_T
                   *model_jointTransform, const g_struct_T mesh, const
                   h_struct_T el)
{
  emxArray_real_T *aeroLoads_ForceValHist;
  emxArray_real_T *aeroLoads_ForceDof;
  emxArray_real_T *udot_j;
  emxArray_real_T *uddot_j;
  b_emxArray_struct_T *dispOut_elStrain;
  double aeroLoads_timeArray_data[2002];
  int aeroLoads_timeArray_size[1];
  int i;
  int loop_ub;
  emxArray_real_T *u_s;
  int platNorm;
  emxArray_real_T *udot_s;
  emxArray_real_T *uddot_s;
  int gbNorm;
  emxArray_real_T *uHist;
  double t_data[51];
  double FReactionHist_data[306];
  b_emxArray_struct_T *strainHist;
  double aziHist_data[51];
  double OmegaHist_data[51];
  double OmegaDotHist_data[51];
  double gbHist_data[51];
  double gbDotHist_data[51];
  double genPower_data[51];
  double gb_s;
  double gbDot_s;
  double azi_s;
  int b_i;
  double FReactionsm1[6];
  c_emxArray_struct_T *elStorage;
  emxArray_real_T *dispOut_displ_sp1;
  emxArray_real_T *dispOut_displddot_sp1;
  emxArray_real_T *dispOut_displdot_sp1;
  emxArray_real_T *u_j;
  emxArray_real_T *u_jLast;
  emxArray_real_T *Fexternal_sub;
  i_struct_T expl_temp;
  j_struct_T dispOut;
  emxArray_real_T *b_Fexternal_sub;
  emxArray_real_T *b_udot_j;
  boolean_T exitg1;
  double omegaCurrent;
  double OmegaDotCurrent;
  boolean_T terminateSimulation;
  double azi_j;
  int numIterations;
  double uNorm;
  double aziNorm;
  int b_model_outFilename_size[2];
  char b_model_outFilename_data[76];
  emxArray_char_T c_model_outFilename_data;
  signed char input_sizes_idx_0;
  emxArray_char_T *line;
  emxArray_char_T *b_r;
  emxArray_char_T *b_r1;
  emxArray_char_T *b_r2;
  emxArray_char_T *r3;
  emxArray_char_T *r4;
  emxArray_char_T *r5;
  emxArray_char_T *r6;
  emxArray_char_T *r7;
  emxArray_char_T *r8;
  emxArray_char_T *r9;
  emxArray_char_T *r10;
  emxArray_char_T *r11;
  emxArray_char_T *r12;
  emxArray_char_T *r13;
  emxArray_char_T *r14;
  emxArray_char_T *r15;
  emxArray_char_T *r16;
  static const char b_cv[198] = { 't', ',', 'a', 'z', 'i', 'H', 'i', 's', 't',
    ',', 'O', 'm', 'e', 'g', 'a', 'H', 'i', 's', 't', ',', 'O', 'm', 'e', 'g',
    'a', 'D', 'o', 't', 'H', 'i', 's', 't', ',', 'g', 'b', 'H', 'i', 's', 't',
    ',', 'g', 'b', 'D', 'o', 't', 'H', 'i', 's', 't', ',', 'g', 'b', 'D', 'o',
    't', 'D', 'o', 't', 'H', 'i', 's', 't', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '1', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '2', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '3', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '4', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '5', ',', 'F', 'R', 'e', 'a', 'c', 't',
    'i', 'o', 'n', 'H', 'i', 's', 't', '6', ',', 'r', 'i', 'g', 'i', 'd', 'D',
    'o', 'f', ',', 'g', 'e', 'n', 'T', 'o', 'r', 'q', 'u', 'e', ',', 'g', 'e',
    'n', 'P', 'o', 'w', 'e', 'r', ',', 't', 'o', 'r', 'q', 'u', 'e', 'D', 'r',
    'i', 'v', 'e', 'S', 'h', 'a', 'f', 't', '\x0a' };

  int sizes_idx_0;
  int c_model_outFilename_size[2];
  signed char input_sizes_idx_1;
  char d_model_outFilename_data[82];
  emxArray_char_T e_model_outFilename_data;
  static const char cv1[10] = { '_', 'u', 'H', 'i', 's', 't', '.', 't', 'x', 't'
  };

  double structureMOI[9];
  int d_model_outFilename_size[2];
  int i1;
  char f_model_outFilename_data[87];
  emxArray_char_T g_model_outFilename_data;
  static const char cv2[15] = { '_', 's', 't', 'r', 'a', 'i', 'n', 'H', 'i', 's',
    't', '.', 't', 'x', 't' };

  double CP2H_tmp[9];
  double b_CP2H_tmp[9];
  static const char cv3[9] = { 'e', 'p', 's', '_', 'x', 'x', '_', '0', ' ' };

  static const char cv4[9] = { 'e', 'p', 's', '_', 'x', 'x', '_', 'z', ' ' };

  static const char cv5[9] = { 'e', 'p', 's', '_', 'x', 'x', '_', 'y', ' ' };

  static const char cv6[9] = { 'g', 'a', 'm', '_', 'x', 'z', '_', '0', ' ' };

  static const char cv7[9] = { 'g', 'a', 'm', '_', 'x', 'z', '_', 'y', ' ' };

  static const char cv8[9] = { 'g', 'a', 'm', '_', 'x', 'y', '_', '0', ' ' };

  static const char cv9[9] = { 'g', 'a', 'm', '_', 'x', 'y', '_', 'z', ' ' };

  emxInit_real_T(&aeroLoads_ForceValHist, 2);
  emxInit_real_T(&aeroLoads_ForceDof, 1);
  emxInit_real_T(&udot_j, 1);
  emxInit_real_T(&uddot_j, 1);
  emxInit_struct_T1(&dispOut_elStrain, 2);

  //  activate platform module
  // ............... flags for module activation ....................
  // modularIteration
  //  Get AeroLoads
  processAeroLoadsBLE(model_aeroloadfile, model_owensfile,
                      aeroLoads_timeArray_data, aeroLoads_timeArray_size,
                      aeroLoads_ForceValHist, aeroLoads_ForceDof);

  //  Declare Variable Type, are set later
  i = udot_j->size[0];
  udot_j->size[0] = 1;
  emxEnsureCapacity_real_T(udot_j, i);
  udot_j->data[0] = 0.0;
  i = uddot_j->size[0];
  uddot_j->size[0] = 1;
  emxEnsureCapacity_real_T(uddot_j, i);
  uddot_j->data[0] = 0.0;
  i = dispOut_elStrain->size[0] * dispOut_elStrain->size[1];
  dispOut_elStrain->size[0] = 1;
  loop_ub = static_cast<int>(mesh.numEl);
  dispOut_elStrain->size[1] = loop_ub;
  emxEnsureCapacity_struct_T1(dispOut_elStrain, i);
  for (i = 0; i < loop_ub; i++) {
    dispOut_elStrain->data[i] = r1;
  }

  emxInit_real_T(&u_s, 1);

  // ................................................................
  //  Rotor mode initialization
  // ..........................................................................
  // Initial rotor speed (Hz)
  //      Omega = OmegaInitial;
  // ensures generator always off for practical purposes
  // ..........................................................................
  //  if(model.turbineStartup ==1 || model.turbineStartup==2)
  //      plotGenSpeedVsTorque([0:.01:5],model.generatorProps);
  //  end
  //
  //  state initialization
  // ......... specify initial conditions .......................
  platNorm = static_cast<int>((mesh.numNodes * 6.0));
  i = u_s->size[0];
  u_s->size[0] = platNorm;
  emxEnsureCapacity_real_T(u_s, i);
  for (i = 0; i < platNorm; i++) {
    u_s->data[i] = 0.0;
  }

  emxInit_real_T(&udot_s, 1);

  // setInitialConditions sets initial conditions
  //  **********************************************************************
  //  *                   Part of the SNL OWENS Toolkit                    *
  //  * Developed by Sandia National Laboratories Wind Energy Technologies *
  //  *             See license.txt for disclaimer information             *
  //  **********************************************************************
  //    [u] =  setInitialConditions(initCond,u,numDOFPerNode)
  //
  //    This function reads initial conditions from file
  //
  //    input:
  //    initCond      = array containing initial conditions
  //                      initCond(i,1) = node number for init cond i
  //                      initCond(i,2) = local DOF number for init cond i
  //                      initCond(i,3) = value for init cond i
  //    u             = displacement vector
  //    numDOFPerNode = number of degrees of freedom per node
  //
  //    output:
  //     u             = displacement vector modified for initial conditions
  // get number of specified initial conditions
  // unspecified initial conditions are assumed to
  // be zero
  i = udot_s->size[0];
  udot_s->size[0] = platNorm;
  emxEnsureCapacity_real_T(udot_s, i);
  for (i = 0; i < platNorm; i++) {
    udot_s->data[i] = 0.0;
  }

  emxInit_real_T(&uddot_s, 1);
  i = uddot_s->size[0];
  uddot_s->size[0] = udot_s->size[0];
  emxEnsureCapacity_real_T(uddot_s, i);
  gbNorm = udot_s->size[0];
  for (i = 0; i < gbNorm; i++) {
    uddot_s->data[i] = udot_s->data[i];
  }

  emxInit_real_T(&uHist, 2);

  // ............................................................
  // define number of time steps
  // define time step size
  i = uHist->size[0] * uHist->size[1];
  uHist->size[0] = platNorm;
  uHist->size[1] = 51;
  emxEnsureCapacity_real_T(uHist, i);
  gbNorm = platNorm * 51;
  for (i = 0; i < gbNorm; i++) {
    uHist->data[i] = 0.0;
  }

  for (i = 0; i < platNorm; i++) {
    uHist->data[i] = 0.0;
  }

  // store initial condition
  // initialize omega_platform, omega_platform_dot, omegaPlatHist
  //  omega_platform = zeros(3,1);
  //  omega_platform_dot = zeros(3,1);
  //  omegaPlatHist(:,1) = omega_platform;
  std::memset(&t_data[0], 0, 51U * sizeof(double));
  std::memset(&FReactionHist_data[0], 0, 306U * sizeof(double));
  emxInit_struct_T1(&strainHist, 2);

  //  strainHist(numTS+1) = struct();
  i = strainHist->size[0] * strainHist->size[1];
  strainHist->size[0] = loop_ub;
  strainHist->size[1] = 50;
  emxEnsureCapacity_struct_T1(strainHist, i);
  loop_ub *= 50;
  for (i = 0; i < loop_ub; i++) {
    strainHist->data[i] = r1;
  }

  // genTorque = zeros(1,numTS+1);
  std::memset(&aziHist_data[0], 0, 51U * sizeof(double));
  std::memset(&OmegaHist_data[0], 0, 51U * sizeof(double));
  std::memset(&OmegaDotHist_data[0], 0, 51U * sizeof(double));
  std::memset(&gbHist_data[0], 0, 51U * sizeof(double));
  std::memset(&gbDotHist_data[0], 0, 51U * sizeof(double));
  std::memset(&genPower_data[0], 0, 51U * sizeof(double));
  t_data[0] = 0.0;

  // initialize various states and variables
  gb_s = 0.0;
  gbDot_s = 0.0;
  azi_s = 0.0;

  //  azi_sm1 = -Omega*delta_t*2*pi;
  aziHist_data[0] = 0.0;
  OmegaHist_data[0] = 0.12000000000000001;
  OmegaDotHist_data[0] = 0.0;
  for (b_i = 0; b_i < 6; b_i++) {
    FReactionsm1[b_i] = 0.0;
    FReactionHist_data[51 * b_i] = 0.0;
  }

  gbHist_data[0] = 0.0;
  gbDotHist_data[0] = 0.0;

  //
  tic();

  //  structural dynamics initialization
  // ..........................................................................
  emxInit_struct_T2(&elStorage, 2);
  if (!b_strcmp(model_analysisType)) {
    initialElementCalculations(model_joint, el.props, el.elLen, el.psi, el.theta,
      el.roll, mesh.numEl, mesh.x, mesh.y, mesh.z, mesh.conn, elStorage);

    // perform initial element calculations for conventional structural dynamics analysis 
  } else {
    // initialize reduced order model
  }

  // calculate structural/platform moi
  // ..........................................................................
  //  Main Loop - iterate for a solution at each time step, i
  b_i = 1;
  emxInit_real_T(&dispOut_displ_sp1, 1);
  emxInit_real_T(&dispOut_displddot_sp1, 1);
  emxInit_real_T(&dispOut_displdot_sp1, 1);
  emxInit_real_T(&u_j, 1);
  emxInit_real_T(&u_jLast, 1);
  emxInit_real_T(&Fexternal_sub, 2);
  emxInitStruct_struct_T2(&expl_temp);
  emxInitStruct_struct_T3(&dispOut);
  emxInit_real_T(&b_Fexternal_sub, 2);
  emxInit_real_T(&b_udot_j, 2);
  exitg1 = false;
  while ((!exitg1) && (b_i - 1 < 50)) {
    //      i %TODO add verbose printing
    //     %% check for specified rotor speed at t(i) + delta_t
    // use discreteized rotor speed profile function
    omegaSpecCheck(t_data[b_i - 1] + 0.002, model_tocp, &omegaCurrent,
                   &OmegaDotCurrent, &terminateSimulation);
    if (terminateSimulation) {
      exitg1 = true;
    } else {
      //     %%
      //     %% initialize "j" Gauss-Sidel iteration
      i = u_j->size[0];
      u_j->size[0] = u_s->size[0];
      emxEnsureCapacity_real_T(u_j, i);
      loop_ub = u_s->size[0];
      for (i = 0; i < loop_ub; i++) {
        u_j->data[i] = u_s->data[i];
      }

      azi_j = azi_s;

      // initialize  platform module related variables only used if(model.hydroOn) 
      //  		Accel_j = Accel;
      //  		Accel_jLast = Accel;
      // gauss-seidel iteration tolerance for various modules
      // max iteration for various modules
      numIterations = 1;
      uNorm = 1.0E+6;
      platNorm = 1000000;
      aziNorm = 1.0E+6;
      gbNorm = 1000000;

      // initialize norms for various module states
      //     %%
      while (((uNorm > 1.0E-8) || (platNorm > 1.0E-8) || (aziNorm > 1.0E-8) ||
              (gbNorm > 1.0E-8)) && (numIterations < 50)) {
        // module gauss-seidel iteration loop
        // calculate CP2H (platform frame to hub frame transformation matrix)
        aziNorm = std::sin(azi_j);
        uNorm = std::cos(azi_j);

        // .... inertial frame to platform transformation matrix ...........
        // .........................................
        //          %% evaluate platform module
        //          %%====================================================
        //          if(model.hydroOn)
        //              % 	Accel_jLast= Accel_j;
        //              Ywec_jLast = Ywec_j;
        //              if(model.platformTurbineYawInteraction == 0)
        //                  FReaction0 = [-FReactionHist(i,1:5)'; 0.0]; %'
        //                  FReaction1 =  [-FReaction_j(1:5); 0.0];
        //              elseif(model.platformTurbineYawInteraction == 1)
        //                  FReaction0 = (-FReactionHist(i,1:6)');
        //                  FReaction1 =  (-FReaction_j(1:6));
        //              elseif(model.platformTurbineYawInteraction == 2)
        //                  FReaction0 = [-FReactionHist(i,1:5)'; genTorque_s];
        //                  FReaction1 =  [-FReaction_j(1:5); genTorque_j];
        //              else
        //                  error('PlatformTurbineYawInteraction flag not recognized.'); 
        //              end
        //              [rbData,Ywec_j,~] = platformModule([t(i) t(i)+delta_t],Ywec(i,:),CP2H,FReaction0,FReaction1,d_input_streamPlatform,d_output_streamPlatform); 
        //          end
        //          %====================================================
        //         %%
        //         %% evaluate generator module
        // ===== generator module ===========================
        // ==================================================
        //         %%
        //         %% evaluate drivetrain module
        //  %===== drivetrain module ==========================
        gb_s = azi_j;
        gbDot_s = omegaCurrent * 2.0 * 3.1415926535897931;

        // ==================================================
        //         %% rotor speed update
        // ===== update rotor speed =========================
        azi_j = azi_s + omegaCurrent * 0.002 * 2.0 * 3.1415926535897931;

        // ===================================================
        //         %%
        //         %% evaluate aerodynamic module (CACTUS ONE-WAY)
        // ======= aerodynamics module ======================
        // ==================================================
        //         %% evaluate aerodynamic module (TU DELFT)
        // ======= aerodynamics module ======================
        // ==================================================
        //         %%
        //         %% compile external forcing on rotor
        // compile forces to supply to structural dynamics solver
        externalForcing(t_data[b_i - 1] + 0.002, aeroLoads_timeArray_data,
                        aeroLoads_timeArray_size, aeroLoads_ForceValHist,
                        aeroLoads_ForceDof, Fexternal_sub, udot_j);
        if (Fexternal_sub->size[1] != 0) {
          gbNorm = Fexternal_sub->size[1];
        } else {
          gbNorm = 0;
        }

        if ((gbNorm == 0) || (Fexternal_sub->size[1] != 0)) {
          input_sizes_idx_0 = 1;
        } else {
          input_sizes_idx_0 = 0;
        }

        if (udot_j->size[0] != 0) {
          sizes_idx_0 = udot_j->size[0];
        } else {
          sizes_idx_0 = 0;
        }

        if ((sizes_idx_0 == 0) || (udot_j->size[0] != 0)) {
          input_sizes_idx_1 = 1;
        } else {
          input_sizes_idx_1 = 0;
        }

        //         %% evaluate structural dynamics
        // call structural dynamics solver
        // initialization of structural dynamics displacements, velocities, accelerations, etc. 
        //              dispData.displ_sm1 = u_sm1; %Not even used
        // specific to TNB, only used if that is the case, and must be declared
        //          if(strcmp(model.analysisType,'ROM'))
        //              dispData.displ_s = u_s;
        //              dispData.displdot_s = udot_s;
        //              dispData.displddot_s = uddot_s;
        //
        //              dispData.eta_s     = eta_s;
        //              dispData.etadot_s  = etadot_s;
        //              dispData.etaddot_s = etaddot_s;
        //          end
        if (!b_strcmp(model_analysisType)) {
          //  evalulate structural dynamics using conventional representation
          eye(structureMOI);
          i = expl_temp.displddot_s->size[0];
          expl_temp.displddot_s->size[0] = uddot_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displddot_s, i);
          loop_ub = uddot_s->size[0];
          for (i = 0; i < loop_ub; i++) {
            expl_temp.displddot_s->data[i] = uddot_s->data[i];
          }

          i = expl_temp.displdot_s->size[0];
          expl_temp.displdot_s->size[0] = udot_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displdot_s, i);
          loop_ub = udot_s->size[0];
          for (i = 0; i < loop_ub; i++) {
            expl_temp.displdot_s->data[i] = udot_s->data[i];
          }

          i = expl_temp.displ_s->size[0];
          expl_temp.displ_s->size[0] = u_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displ_s, i);
          loop_ub = u_s->size[0];
          for (i = 0; i < loop_ub; i++) {
            expl_temp.displ_s->data[i] = u_s->data[i];
          }

          platNorm = input_sizes_idx_0;
          CP2H_tmp[0] = uNorm;
          CP2H_tmp[3] = aziNorm;
          CP2H_tmp[6] = 0.0;
          CP2H_tmp[1] = -aziNorm;
          CP2H_tmp[4] = uNorm;
          CP2H_tmp[7] = 0.0;
          CP2H_tmp[2] = 0.0;
          CP2H_tmp[5] = 0.0;
          CP2H_tmp[8] = 1.0;
          i = b_Fexternal_sub->size[0] * b_Fexternal_sub->size[1];
          b_Fexternal_sub->size[0] = input_sizes_idx_0;
          b_Fexternal_sub->size[1] = gbNorm;
          emxEnsureCapacity_real_T(b_Fexternal_sub, i);
          for (i = 0; i < gbNorm; i++) {
            for (i1 = 0; i1 < platNorm; i1++) {
              b_Fexternal_sub->data[b_Fexternal_sub->size[0] * i] =
                Fexternal_sub->data[input_sizes_idx_0 * i];
            }
          }

          i = b_udot_j->size[0] * b_udot_j->size[1];
          b_udot_j->size[0] = sizes_idx_0;
          b_udot_j->size[1] = input_sizes_idx_1;
          emxEnsureCapacity_real_T(b_udot_j, i);
          loop_ub = input_sizes_idx_1;
          for (i = 0; i < loop_ub; i++) {
            for (i1 = 0; i1 < sizes_idx_0; i1++) {
              b_udot_j->data[i1] = udot_j->data[i1];
            }
          }

          for (i = 0; i < 3; i++) {
            aziNorm = CP2H_tmp[i + 3];
            i1 = static_cast<int>(CP2H_tmp[i + 6]);
            for (platNorm = 0; platNorm < 3; platNorm++) {
              b_CP2H_tmp[i + 3 * platNorm] = (CP2H_tmp[i] * structureMOI[3 *
                platNorm] + aziNorm * structureMOI[3 * platNorm + 1]) +
                static_cast<double>(i1) * structureMOI[3 * platNorm + 2];
            }
          }

          structuralDynamicsTransient(model_analysisType, model_RayleighAlpha,
            model_RayleighBeta, model_BC_numpBC, model_BC_pBC, model_joint,
            model_jointTransform, mesh.numEl, mesh.x, mesh.y, mesh.z, mesh.conn,
            el, expl_temp, omegaCurrent, OmegaDotCurrent, elStorage,
            b_Fexternal_sub, b_udot_j, b_CP2H_tmp, &dispOut, FReactionsm1);
          i = dispOut_elStrain->size[0] * dispOut_elStrain->size[1];
          dispOut_elStrain->size[0] = 1;
          dispOut_elStrain->size[1] = dispOut.elStrain->size[1];
          emxEnsureCapacity_struct_T1(dispOut_elStrain, i);
          loop_ub = dispOut.elStrain->size[0] * dispOut.elStrain->size[1];
          for (i = 0; i < loop_ub; i++) {
            dispOut_elStrain->data[i] = dispOut.elStrain->data[i];
          }

          i = dispOut_displ_sp1->size[0];
          dispOut_displ_sp1->size[0] = dispOut.displ_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displ_sp1, i);
          loop_ub = dispOut.displ_sp1->size[0];
          for (i = 0; i < loop_ub; i++) {
            dispOut_displ_sp1->data[i] = dispOut.displ_sp1->data[i];
          }

          i = dispOut_displddot_sp1->size[0];
          dispOut_displddot_sp1->size[0] = dispOut.displddot_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displddot_sp1, i);
          loop_ub = dispOut.displddot_sp1->size[0];
          for (i = 0; i < loop_ub; i++) {
            dispOut_displddot_sp1->data[i] = dispOut.displddot_sp1->data[i];
          }

          i = dispOut_displdot_sp1->size[0];
          dispOut_displdot_sp1->size[0] = dispOut.displdot_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displdot_sp1, i);
          loop_ub = dispOut.displdot_sp1->size[0];
          for (i = 0; i < loop_ub; i++) {
            dispOut_displdot_sp1->data[i] = dispOut.displdot_sp1->data[i];
          }
        }

        // update last iteration displacement vector
        i = u_jLast->size[0];
        u_jLast->size[0] = u_j->size[0];
        emxEnsureCapacity_real_T(u_jLast, i);
        loop_ub = u_j->size[0];
        for (i = 0; i < loop_ub; i++) {
          u_jLast->data[i] = u_j->data[i];
        }

        i = u_j->size[0];
        u_j->size[0] = dispOut_displ_sp1->size[0];
        emxEnsureCapacity_real_T(u_j, i);
        loop_ub = dispOut_displ_sp1->size[0];
        for (i = 0; i < loop_ub; i++) {
          u_j->data[i] = dispOut_displ_sp1->data[i];
        }

        // update current estimates of velocity, acceleration
        //  Only used for TNB, but must be declared
        i = udot_j->size[0];
        udot_j->size[0] = dispOut_displdot_sp1->size[0];
        emxEnsureCapacity_real_T(udot_j, i);
        loop_ub = dispOut_displdot_sp1->size[0];
        for (i = 0; i < loop_ub; i++) {
          udot_j->data[i] = dispOut_displdot_sp1->data[i];
        }

        i = uddot_j->size[0];
        uddot_j->size[0] = dispOut_displddot_sp1->size[0];
        emxEnsureCapacity_real_T(uddot_j, i);
        loop_ub = dispOut_displddot_sp1->size[0];
        for (i = 0; i < loop_ub; i++) {
          uddot_j->data[i] = dispOut_displddot_sp1->data[i];
        }

        //         %%
        //         %% calculate norms
        i = u_jLast->size[0];
        u_jLast->size[0] = dispOut_displ_sp1->size[0];
        emxEnsureCapacity_real_T(u_jLast, i);
        loop_ub = dispOut_displ_sp1->size[0];
        for (i = 0; i < loop_ub; i++) {
          u_jLast->data[i] = dispOut_displ_sp1->data[i] - u_jLast->data[i];
        }

        uNorm = b_norm(u_jLast) / b_norm(dispOut_displ_sp1);

        // structural dynamics displacement iteration norm
        aziNorm = std::abs(azi_j - gb_s) / std::abs(azi_j);

        // rotor azimuth iteration norm
        platNorm = 0;
        gbNorm = 0;
        numIterations++;
      }

      // end iteration while loop
      //     %% calculate converged generator torque/power
      genPower_data[b_i] = 0.0 * (gbDot_s * 2.0 * 3.1415926535897931);

      //     %% update timestepping variables and other states, store in history arrays 
      if (c_strcmp(model_analysisType)) {
        i = u_s->size[0];
        u_s->size[0] = u_j->size[0];
        emxEnsureCapacity_real_T(u_s, i);
        loop_ub = u_j->size[0];
        for (i = 0; i < loop_ub; i++) {
          u_s->data[i] = u_j->data[i];
        }

        i = udot_s->size[0];
        udot_s->size[0] = udot_j->size[0];
        emxEnsureCapacity_real_T(udot_s, i);
        loop_ub = udot_j->size[0];
        for (i = 0; i < loop_ub; i++) {
          udot_s->data[i] = udot_j->data[i];
        }

        i = uddot_s->size[0];
        uddot_s->size[0] = uddot_j->size[0];
        emxEnsureCapacity_real_T(uddot_s, i);
        loop_ub = uddot_j->size[0];
        for (i = 0; i < loop_ub; i++) {
          uddot_s->data[i] = uddot_j->data[i];
        }
      }

      loop_ub = u_s->size[0];
      for (i = 0; i < loop_ub; i++) {
        uHist->data[i + uHist->size[0] * b_i] = u_s->data[i];
      }

      for (i = 0; i < 6; i++) {
        FReactionHist_data[b_i + 51 * i] = FReactionsm1[i];
      }

      loop_ub = strainHist->size[0];
      for (i = 0; i < loop_ub; i++) {
        strainHist->data[i + strainHist->size[0] * (b_i - 1)] =
          dispOut_elStrain->data[i];
      }

      t_data[b_i] = t_data[b_i - 1] + 0.002;
      azi_s = azi_j;
      aziHist_data[b_i] = azi_j;
      OmegaHist_data[b_i] = omegaCurrent;
      OmegaDotHist_data[b_i] = OmegaDotCurrent;
      gbHist_data[b_i] = gb_s;
      gbDotHist_data[b_i] = gbDot_s;

      // genTorque(i+1) = genTorque_s;
      for (i = 0; i < 6; i++) {
        FReactionHist_data[b_i + 51 * i] = FReactionsm1[i];
      }

      //     %%
      //     %% check rotor speed for generator operation
      //     %%
      b_i++;
    }
  }

  emxFree_real_T(&b_udot_j);
  emxFree_real_T(&b_Fexternal_sub);
  emxFreeStruct_struct_T3(&dispOut);
  emxFreeStruct_struct_T2(&expl_temp);
  emxFree_real_T(&Fexternal_sub);
  emxFree_real_T(&u_jLast);
  emxFree_real_T(&u_j);
  emxFree_struct_T2(&elStorage);
  emxFree_real_T(&uddot_s);
  emxFree_real_T(&udot_s);
  emxFree_real_T(&u_s);
  emxFree_real_T(&dispOut_displdot_sp1);
  emxFree_real_T(&dispOut_displddot_sp1);
  emxFree_real_T(&dispOut_displ_sp1);
  emxFree_struct_T1(&dispOut_elStrain);
  emxFree_real_T(&uddot_j);
  emxFree_real_T(&udot_j);
  emxFree_real_T(&aeroLoads_ForceDof);
  emxFree_real_T(&aeroLoads_ForceValHist);

  // end timestep loop
  //  kill platform module process
  //  kill aerodynamic module process
  //
  // toc
  //  save aeroOutputArray
  // save simulation data in .mat file
  //  save(model.outFilename,'t','uHist','aziHist','OmegaHist','OmegaDotHist','gbHist','gbDotHist','gbDotDotHist','FReactionHist','rigidDof','genTorque','genPower','torqueDriveShaft','strainHist'); 
  //  fprintf('%s\n','Output Saving Not Currently Implemented')
  // Writefile
  printf("%s\n", "Writing Verification File");
  fflush(stdout);
  loop_ub = model_outFilename_size[1];
  b_model_outFilename_size[0] = 1;
  b_model_outFilename_size[1] = model_outFilename_size[1];
  if (0 <= loop_ub - 4) {
    std::memcpy(&b_model_outFilename_data[0], &model_outFilename_data[0],
                (loop_ub + -3) * sizeof(char));
  }

  b_model_outFilename_data[model_outFilename_size[1] - 3] = 't';
  b_model_outFilename_data[model_outFilename_size[1] - 2] = 'x';
  b_model_outFilename_data[model_outFilename_size[1] - 1] = 't';
  c_model_outFilename_data.data = &b_model_outFilename_data[0];
  c_model_outFilename_data.size = &b_model_outFilename_size[0];
  c_model_outFilename_data.allocatedSize = 76;
  c_model_outFilename_data.numDimensions = 2;
  c_model_outFilename_data.canFreeData = false;
  input_sizes_idx_0 = c_cfopen(&c_model_outFilename_data, "wb+");

  // open/create new for writing and discard existing data
  emxInit_char_T(&line, 2);
  emxInit_char_T(&b_r, 2);
  emxInit_char_T(&b_r1, 2);
  emxInit_char_T(&b_r2, 2);
  emxInit_char_T(&r3, 2);
  emxInit_char_T(&r4, 2);
  emxInit_char_T(&r5, 2);
  emxInit_char_T(&r6, 2);
  emxInit_char_T(&r7, 2);
  emxInit_char_T(&r8, 2);
  emxInit_char_T(&r9, 2);
  emxInit_char_T(&r10, 2);
  emxInit_char_T(&r11, 2);
  emxInit_char_T(&r12, 2);
  emxInit_char_T(&r13, 2);
  emxInit_char_T(&r14, 2);
  emxInit_char_T(&r15, 2);
  emxInit_char_T(&r16, 2);
  for (b_i = 0; b_i < 52; b_i++) {
    if (b_i == 0) {
      i = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] = 198;
      emxEnsureCapacity_char_T(line, i);
      for (i = 0; i < 198; i++) {
        line->data[i] = b_cv[i];
      }
    } else {
      b_sprintf(t_data[b_i - 1], b_r);
      b_sprintf(aziHist_data[b_i - 1], b_r1);
      b_sprintf(OmegaHist_data[b_i - 1], b_r2);
      b_sprintf(OmegaDotHist_data[b_i - 1], r3);
      b_sprintf(gbHist_data[b_i - 1], r4);
      b_sprintf(gbDotHist_data[b_i - 1], r5);
      b_sprintf(0.0, r6);
      b_sprintf(FReactionHist_data[b_i - 1], r7);
      b_sprintf(FReactionHist_data[b_i + 50], r8);
      b_sprintf(FReactionHist_data[b_i + 101], r9);
      b_sprintf(FReactionHist_data[b_i + 152], r10);
      b_sprintf(FReactionHist_data[b_i + 203], r11);
      b_sprintf(FReactionHist_data[b_i + 254], r12);
      b_sprintf(0.0, r13);
      b_sprintf(0.0, r14);
      b_sprintf(genPower_data[b_i - 1], r15);
      b_sprintf(0.0, r16);
      i = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] = ((((((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->
        size[1]) + r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1]) +
        r7->size[1]) + r8->size[1]) + r9->size[1]) + r10->size[1]) + r11->size[1])
                           + r12->size[1]) + r13->size[1]) + r14->size[1]) +
                        r15->size[1]) + r16->size[1]) + 17;
      emxEnsureCapacity_char_T(line, i);
      loop_ub = b_r->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[i] = b_r->data[i];
      }

      line->data[b_r->size[1]] = ',';
      loop_ub = b_r1->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(i + b_r->size[1]) + 1] = b_r1->data[i];
      }

      line->data[(b_r->size[1] + b_r1->size[1]) + 1] = ',';
      loop_ub = b_r2->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[((i + b_r->size[1]) + b_r1->size[1]) + 2] = b_r2->data[i];
      }

      line->data[((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) + 2] = ',';
      loop_ub = r3->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) + 3] =
          r3->data[i];
      }

      line->data[(((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) + r3->size[1])
        + 3] = ',';
      loop_ub = r4->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) +
                    r3->size[1]) + 4] = r4->data[i];
      }

      line->data[((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) + r3->size[1])
                  + r4->size[1]) + 4] = ',';
      loop_ub = r5->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) +
                     r3->size[1]) + r4->size[1]) + 5] = r5->data[i];
      }

      line->data[(((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) + r3->size
                    [1]) + r4->size[1]) + r5->size[1]) + 5] = ',';
      loop_ub = r6->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) +
                      r3->size[1]) + r4->size[1]) + r5->size[1]) + 6] = r6->
          data[i];
      }

      line->data[((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) + r3->
                     size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1]) + 6] =
        ',';
      loop_ub = r7->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) +
                       r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
          + 7] = r7->data[i];
      }

      line->data[(((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                      r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
                  + r7->size[1]) + 7] = ',';
      loop_ub = r8->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[((((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1]) +
                        r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
                    + r7->size[1]) + 8] = r8->data[i];
      }

      line->data[((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                       r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
                   + r7->size[1]) + r8->size[1]) + 8] = ',';
      loop_ub = r9->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1])
                         + r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->
                      size[1]) + r7->size[1]) + r8->size[1]) + 9] = r9->data[i];
      }

      line->data[(((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                        r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
                    + r7->size[1]) + r8->size[1]) + r9->size[1]) + 9] = ',';
      loop_ub = r10->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[((((((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1])
                          + r3->size[1]) + r4->size[1]) + r5->size[1]) +
                       r6->size[1]) + r7->size[1]) + r8->size[1]) + r9->size[1])
          + 10] = r10->data[i];
      }

      line->data[((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                         r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1])
                     + r7->size[1]) + r8->size[1]) + r9->size[1]) + r10->size[1])
        + 10] = ',';
      loop_ub = r11->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + b_r->size[1]) + b_r1->size[1]) + b_r2->size[1])
                           + r3->size[1]) + r4->size[1]) + r5->size[1]) +
                        r6->size[1]) + r7->size[1]) + r8->size[1]) + r9->size[1])
                    + r10->size[1]) + 11] = r11->data[i];
      }

      line->data[(((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                          r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size
                       [1]) + r7->size[1]) + r8->size[1]) + r9->size[1]) +
                   r10->size[1]) + r11->size[1]) + 11] = ',';
      loop_ub = r12->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + (b_r->size[1] + b_r1->size[1])) + b_r2->size[1])
                            + r3->size[1]) + r4->size[1]) + r5->size[1]) +
                         r6->size[1]) + r7->size[1]) + r8->size[1]) + r9->size[1])
                     + r10->size[1]) + r11->size[1]) + 12] = r12->data[i];
      }

      line->data[((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                           r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->
                        size[1]) + r7->size[1]) + r8->size[1]) + r9->size[1]) +
                    r10->size[1]) + r11->size[1]) + r12->size[1]) + 12] = ',';
      loop_ub = r13->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + ((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]))
                             + r3->size[1]) + r4->size[1]) + r5->size[1]) +
                          r6->size[1]) + r7->size[1]) + r8->size[1]) + r9->size
                       [1]) + r10->size[1]) + r11->size[1]) + r12->size[1]) + 13]
          = r13->data[i];
      }

      line->data[(((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
                            r3->size[1]) + r4->size[1]) + r5->size[1]) +
                         r6->size[1]) + r7->size[1]) + r8->size[1]) + r9->size[1])
                     + r10->size[1]) + r11->size[1]) + r12->size[1]) + r13->
                  size[1]) + 13] = ',';
      loop_ub = r14->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + (((b_r->size[1] + b_r1->size[1]) + b_r2->size
          [1]) + r3->size[1])) + r4->size[1]) + r5->size[1]) + r6->size[1]) +
                          r7->size[1]) + r8->size[1]) + r9->size[1]) + r10->
                       size[1]) + r11->size[1]) + r12->size[1]) + r13->size[1])
          + 14] = r14->data[i];
      }

      line->data[((((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
        r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1]) + r7->size[1])
                        + r8->size[1]) + r9->size[1]) + r10->size[1]) +
                     r11->size[1]) + r12->size[1]) + r13->size[1]) + r14->size[1])
        + 14] = ',';
      loop_ub = r15->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + ((((b_r->size[1] + b_r1->size[1]) + b_r2->
          size[1]) + r3->size[1]) + r4->size[1])) + r5->size[1]) + r6->size[1])
                           + r7->size[1]) + r8->size[1]) + r9->size[1]) +
                        r10->size[1]) + r11->size[1]) + r12->size[1]) +
                     r13->size[1]) + r14->size[1]) + 15] = r15->data[i];
      }

      line->data[(((((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1]) +
        r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1]) + r7->size[1])
                         + r8->size[1]) + r9->size[1]) + r10->size[1]) +
                      r11->size[1]) + r12->size[1]) + r13->size[1]) + r14->size
                   [1]) + r15->size[1]) + 15] = ',';
      loop_ub = r16->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(((((((((((i + (((((b_r->size[1] + b_r1->size[1]) +
          b_r2->size[1]) + r3->size[1]) + r4->size[1]) + r5->size[1])) +
                             r6->size[1]) + r7->size[1]) + r8->size[1]) +
                          r9->size[1]) + r10->size[1]) + r11->size[1]) +
                       r12->size[1]) + r13->size[1]) + r14->size[1]) + r15->
                    size[1]) + 16] = r16->data[i];
      }

      line->data[((((((((((((((((b_r->size[1] + b_r1->size[1]) + b_r2->size[1])
        + r3->size[1]) + r4->size[1]) + r5->size[1]) + r6->size[1]) + r7->size[1])
                          + r8->size[1]) + r9->size[1]) + r10->size[1]) +
                       r11->size[1]) + r12->size[1]) + r13->size[1]) + r14->
                    size[1]) + r15->size[1]) + r16->size[1]) + 16] = '\x0a';
    }

    b_fwrite(static_cast<double>(input_sizes_idx_0), line);
  }

  emxFree_char_T(&r16);
  emxFree_char_T(&r15);
  emxFree_char_T(&r14);
  emxFree_char_T(&r13);
  emxFree_char_T(&r12);
  emxFree_char_T(&r11);
  emxFree_char_T(&r10);
  emxFree_char_T(&r9);
  emxFree_char_T(&r8);
  emxFree_char_T(&r7);
  emxFree_char_T(&r6);
  emxFree_char_T(&r5);
  emxFree_char_T(&r4);
  emxFree_char_T(&r3);
  emxFree_char_T(&b_r2);
  cfclose(static_cast<double>(input_sizes_idx_0));

  //  Save uHist
  if (1 > model_outFilename_size[1] - 4) {
    loop_ub = 0;
  } else {
    loop_ub = model_outFilename_size[1] - 4;
  }

  c_model_outFilename_size[0] = 1;
  c_model_outFilename_size[1] = loop_ub + 10;
  if (0 <= loop_ub - 1) {
    std::memcpy(&d_model_outFilename_data[0], &model_outFilename_data[0],
                loop_ub * sizeof(char));
  }

  for (i = 0; i < 10; i++) {
    d_model_outFilename_data[i + loop_ub] = cv1[i];
  }

  e_model_outFilename_data.data = &d_model_outFilename_data[0];
  e_model_outFilename_data.size = &c_model_outFilename_size[0];
  e_model_outFilename_data.allocatedSize = 82;
  e_model_outFilename_data.numDimensions = 2;
  e_model_outFilename_data.canFreeData = false;
  input_sizes_idx_0 = c_cfopen(&e_model_outFilename_data, "wb+");

  // open/create new for writing and discard existing data
  for (b_i = 0; b_i < 52; b_i++) {
    if (b_i == 0) {
      c_fwrite(static_cast<double>(input_sizes_idx_0));
    } else {
      i = uHist->size[0] - 2;
      for (platNorm = 0; platNorm <= i; platNorm++) {
        b_sprintf(t_data[b_i - 1], b_r);
        b_sprintf(uHist->data[platNorm + uHist->size[0] * (b_i - 1)], b_r1);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = (b_r->size[1] + b_r1->size[1]) + 1;
        emxEnsureCapacity_char_T(line, i1);
        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1] = b_r->data[i1];
        }

        line->data[b_r->size[1]] = ' ';
        loop_ub = b_r1->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[(i1 + b_r->size[1]) + 1] = b_r1->data[i1];
        }

        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
      }

      b_sprintf(t_data[b_i - 1], b_r);
      b_sprintf(uHist->data[(uHist->size[0] + uHist->size[0] * (b_i - 1)) - 1],
                b_r1);
      i = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] = (b_r->size[1] + b_r1->size[1]) + 2;
      emxEnsureCapacity_char_T(line, i);
      loop_ub = b_r->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[i] = b_r->data[i];
      }

      line->data[b_r->size[1]] = ' ';
      loop_ub = b_r1->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(i + b_r->size[1]) + 1] = b_r1->data[i];
      }

      line->data[(b_r->size[1] + b_r1->size[1]) + 1] = '\x0a';
      b_fwrite(static_cast<double>(input_sizes_idx_0), line);
    }
  }

  emxFree_real_T(&uHist);
  cfclose(static_cast<double>(input_sizes_idx_0));

  //  Save strainHist
  if (1 > model_outFilename_size[1] - 4) {
    loop_ub = 0;
  } else {
    loop_ub = model_outFilename_size[1] - 4;
  }

  d_model_outFilename_size[0] = 1;
  d_model_outFilename_size[1] = loop_ub + 15;
  if (0 <= loop_ub - 1) {
    std::memcpy(&f_model_outFilename_data[0], &model_outFilename_data[0],
                loop_ub * sizeof(char));
  }

  for (i = 0; i < 15; i++) {
    f_model_outFilename_data[i + loop_ub] = cv2[i];
  }

  g_model_outFilename_data.data = &f_model_outFilename_data[0];
  g_model_outFilename_data.size = &d_model_outFilename_size[0];
  g_model_outFilename_data.allocatedSize = 87;
  g_model_outFilename_data.numDimensions = 2;
  g_model_outFilename_data.canFreeData = false;
  input_sizes_idx_0 = c_cfopen(&g_model_outFilename_data, "wb+");

  // open/create new for writing and discard existing data
  for (b_i = 0; b_i < 51; b_i++) {
    if (b_i == 0) {
      d_fwrite(static_cast<double>(input_sizes_idx_0));
      b_sprintf(51.0, b_r);
      b_sprintf(static_cast<double>(strainHist->size[0]), b_r1);
      i = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] = (b_r->size[1] + b_r1->size[1]) + 2;
      emxEnsureCapacity_char_T(line, i);
      loop_ub = b_r->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[i] = b_r->data[i];
      }

      line->data[b_r->size[1]] = ' ';
      loop_ub = b_r1->size[1];
      for (i = 0; i < loop_ub; i++) {
        line->data[(i + b_r->size[1]) + 1] = b_r1->data[i];
      }

      line->data[(b_r->size[1] + b_r1->size[1]) + 1] = '\x0a';
      b_fwrite(static_cast<double>(input_sizes_idx_0), line);
    } else {
      i = strainHist->size[0] - 1;
      for (platNorm = 0; platNorm <= i; platNorm++) {
        e_fwrite(static_cast<double>(input_sizes_idx_0));
        b_sprintf(t_data[b_i - 1], b_r);
        b_sprintf(static_cast<double>(platNorm) + 1.0, b_r1);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = (b_r->size[1] + b_r1->size[1]) + 2;
        emxEnsureCapacity_char_T(line, i1);
        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1] = b_r->data[i1];
        }

        line->data[b_r->size[1]] = ' ';
        loop_ub = b_r1->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[(i1 + b_r->size[1]) + 1] = b_r1->data[i1];
        }

        line->data[(b_r->size[1] + b_r1->size[1]) + 1] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .eps_xx_0[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv3[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  eps_xx_0[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv3[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .eps_xx_z[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv4[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  eps_xx_z[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv4[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .eps_xx_y[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv5[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  eps_xx_y[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv5[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .gam_xz_0[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv6[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  gam_xz_0[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv6[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .gam_xz_y[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv7[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  gam_xz_y[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv7[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .gam_xy_0[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv8[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  gam_xy_0[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv8[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        for (gbNorm = 0; gbNorm < 3; gbNorm++) {
          b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)]
                    .gam_xy_z[gbNorm], b_r);
          i1 = line->size[0] * line->size[1];
          line->size[0] = 1;
          line->size[1] = b_r->size[1] + 9;
          emxEnsureCapacity_char_T(line, i1);
          for (i1 = 0; i1 < 9; i1++) {
            line->data[i1] = cv9[i1];
          }

          loop_ub = b_r->size[1];
          for (i1 = 0; i1 < loop_ub; i1++) {
            line->data[i1 + 9] = b_r->data[i1];
          }

          b_fwrite(static_cast<double>(input_sizes_idx_0), line);
        }

        b_sprintf(strainHist->data[platNorm + strainHist->size[0] * (b_i - 1)].
                  gam_xy_z[3], b_r);
        i1 = line->size[0] * line->size[1];
        line->size[0] = 1;
        line->size[1] = b_r->size[1] + 10;
        emxEnsureCapacity_char_T(line, i1);
        for (i1 = 0; i1 < 9; i1++) {
          line->data[i1] = cv9[i1];
        }

        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          line->data[i1 + 9] = b_r->data[i1];
        }

        line->data[b_r->size[1] + 9] = '\x0a';
        b_fwrite(static_cast<double>(input_sizes_idx_0), line);
      }
    }
  }

  emxFree_char_T(&b_r1);
  emxFree_char_T(&b_r);
  emxFree_char_T(&line);
  emxFree_struct_T1(&strainHist);
  cfclose(static_cast<double>(input_sizes_idx_0));
}

//
// File trailer for transientExec.cpp
//
// [EOF]
//
