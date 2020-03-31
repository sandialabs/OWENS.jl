//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: transientExec.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "transientExec.h"
#include "externalForcing.h"
#include "initialElementCalculations.h"
#include "mapCactusLoadsFile.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include "structuralDynamicsTransient.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include "tic.h"
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <string.h>

// Function Definitions

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
//                const emxArray_real_T *model_jointTransform
//                const e_struct_T mesh
//                const f_struct_T el
// Return Type  : void
//
void transientExec(const char model_analysisType[3], const double model_tocp[2],
                   const emxArray_char_T *model_aeroloadfile, const
                   emxArray_char_T *model_owensfile, double model_RayleighAlpha,
                   double model_RayleighBeta, double model_BC_numpBC, const
                   emxArray_real_T *model_BC_pBC, const emxArray_real_T
                   *model_joint, const emxArray_real_T *model_jointTransform,
                   const e_struct_T mesh, const f_struct_T el)
{
  emxArray_char_T *b_model_aeroloadfile;
  int gbNorm;
  int platNorm;
  int i;
  emxArray_char_T *c_model_aeroloadfile;
  static const char cv[5] = { '.', 'g', 'e', 'o', 'm' };

  emxArray_char_T *b_model_owensfile;
  static const char cv1[16] = { '_', 'E', 'l', 'e', 'm', 'e', 'n', 't', 'D', 'a',
    't', 'a', '.', 'c', 's', 'v' };

  emxArray_char_T *c_model_owensfile;
  emxArray_char_T *d_model_owensfile;
  emxArray_char_T *e_model_owensfile;
  emxArray_real_T *aeroLoads_ForceValHist;
  static const char cv2[5] = { '.', 'm', 'e', 's', 'h' };

  emxArray_real_T *aeroLoads_ForceDof;
  emxArray_real_T *udot_j;
  emxArray_real_T *uddot_j;
  emxArray_real_T *u_s;
  double aeroLoads_timeArray_data[2002];
  int aeroLoads_timeArray_size[1];
  emxArray_real_T *udot_s;
  emxArray_real_T *uddot_s;
  double t_data[251];
  double azi_s;
  b_emxArray_struct_T *elStorage;
  int b_i;
  emxArray_real_T *dispOut_displ_sp1;
  emxArray_real_T *dispOut_displddot_sp1;
  emxArray_real_T *dispOut_displdot_sp1;
  emxArray_real_T *u_j;
  emxArray_real_T *u_jLast;
  emxArray_real_T *Fexternal_sub;
  g_struct_T expl_temp;
  h_struct_T dispOut;
  emxArray_real_T *b_Fexternal_sub;
  emxArray_real_T *b_udot_j;
  boolean_T exitg1;
  int exitg2;
  double omegaCurrent;
  double Vq;
  double b_Vq;
  double OmegaDotCurrent;
  double azi_j;
  int numIterations;
  double uNorm;
  double aziNorm;
  signed char CN2P[9];
  double azi_jLast;
  int sizes_idx_1;
  signed char input_sizes_idx_0;
  int sizes_idx_0;
  signed char input_sizes_idx_1;
  double CP2H_tmp[9];
  int i1;
  double b_CP2H_tmp[9];
  double FReactionsm1[6];
  emxInit_char_T(&b_model_aeroloadfile, 2);

  //  activate platform module
  // ............... flags for module activation ....................
  // modularIteration
  //  Get AeroLoads
  if (1 > model_aeroloadfile->size[1] - 16) {
    gbNorm = 0;
  } else {
    gbNorm = model_aeroloadfile->size[1] - 16;
  }

  // cut off the _ElementData.csv
  if (1 > model_owensfile->size[1] - 6) {
    platNorm = 0;
  } else {
    platNorm = model_owensfile->size[1] - 6;
  }

  // cut off the .owens
  i = b_model_aeroloadfile->size[0] * b_model_aeroloadfile->size[1];
  b_model_aeroloadfile->size[0] = 1;
  b_model_aeroloadfile->size[1] = gbNorm + 5;
  emxEnsureCapacity_char_T(b_model_aeroloadfile, i);
  for (i = 0; i < gbNorm; i++) {
    b_model_aeroloadfile->data[i] = model_aeroloadfile->data[i];
  }

  for (i = 0; i < 5; i++) {
    b_model_aeroloadfile->data[i + gbNorm] = cv[i];
  }

  emxInit_char_T(&c_model_aeroloadfile, 2);
  i = c_model_aeroloadfile->size[0] * c_model_aeroloadfile->size[1];
  c_model_aeroloadfile->size[0] = 1;
  c_model_aeroloadfile->size[1] = gbNorm + 16;
  emxEnsureCapacity_char_T(c_model_aeroloadfile, i);
  for (i = 0; i < gbNorm; i++) {
    c_model_aeroloadfile->data[i] = model_aeroloadfile->data[i];
  }

  for (i = 0; i < 16; i++) {
    c_model_aeroloadfile->data[i + gbNorm] = cv1[i];
  }

  emxInit_char_T(&b_model_owensfile, 2);
  i = b_model_owensfile->size[0] * b_model_owensfile->size[1];
  b_model_owensfile->size[0] = 1;
  b_model_owensfile->size[1] = platNorm + 4;
  emxEnsureCapacity_char_T(b_model_owensfile, i);
  for (i = 0; i < platNorm; i++) {
    b_model_owensfile->data[i] = model_owensfile->data[i];
  }

  emxInit_char_T(&c_model_owensfile, 2);
  b_model_owensfile->data[platNorm] = '.';
  b_model_owensfile->data[platNorm + 1] = 'b';
  b_model_owensfile->data[platNorm + 2] = 'l';
  b_model_owensfile->data[platNorm + 3] = 'd';
  i = c_model_owensfile->size[0] * c_model_owensfile->size[1];
  c_model_owensfile->size[0] = 1;
  c_model_owensfile->size[1] = platNorm + 3;
  emxEnsureCapacity_char_T(c_model_owensfile, i);
  for (i = 0; i < platNorm; i++) {
    c_model_owensfile->data[i] = model_owensfile->data[i];
  }

  emxInit_char_T(&d_model_owensfile, 2);
  c_model_owensfile->data[platNorm] = '.';
  c_model_owensfile->data[platNorm + 1] = 'e';
  c_model_owensfile->data[platNorm + 2] = 'l';
  i = d_model_owensfile->size[0] * d_model_owensfile->size[1];
  d_model_owensfile->size[0] = 1;
  d_model_owensfile->size[1] = platNorm + 4;
  emxEnsureCapacity_char_T(d_model_owensfile, i);
  for (i = 0; i < platNorm; i++) {
    d_model_owensfile->data[i] = model_owensfile->data[i];
  }

  emxInit_char_T(&e_model_owensfile, 2);
  d_model_owensfile->data[platNorm] = '.';
  d_model_owensfile->data[platNorm + 1] = 'o';
  d_model_owensfile->data[platNorm + 2] = 'r';
  d_model_owensfile->data[platNorm + 3] = 't';
  i = e_model_owensfile->size[0] * e_model_owensfile->size[1];
  e_model_owensfile->size[0] = 1;
  e_model_owensfile->size[1] = platNorm + 5;
  emxEnsureCapacity_char_T(e_model_owensfile, i);
  for (i = 0; i < platNorm; i++) {
    e_model_owensfile->data[i] = model_owensfile->data[i];
  }

  for (i = 0; i < 5; i++) {
    e_model_owensfile->data[i + platNorm] = cv2[i];
  }

  emxInit_real_T(&aeroLoads_ForceValHist, 2);
  emxInit_real_T(&aeroLoads_ForceDof, 1);
  emxInit_real_T(&udot_j, 1);
  emxInit_real_T(&uddot_j, 1);
  emxInit_real_T(&u_s, 1);
  mapCactusLoadsFile(b_model_aeroloadfile, c_model_aeroloadfile,
                     b_model_owensfile, c_model_owensfile, d_model_owensfile,
                     e_model_owensfile, aeroLoads_timeArray_data,
                     aeroLoads_timeArray_size, aeroLoads_ForceValHist,
                     aeroLoads_ForceDof);

  //  method constrained by not passing the .mat filename
  //  save('aeroLoads.mat','timeArray','ForceValHist','ForceDof');
  //  disp('New aeroLoads.mat file saved.')
  //  Declare Variable Type, are set later
  i = udot_j->size[0];
  udot_j->size[0] = 1;
  emxEnsureCapacity_real_T(udot_j, i);
  udot_j->data[0] = 0.0;
  i = uddot_j->size[0];
  uddot_j->size[0] = 1;
  emxEnsureCapacity_real_T(uddot_j, i);
  uddot_j->data[0] = 0.0;

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
  gbNorm = static_cast<int>((mesh.numNodes * 6.0));
  i = u_s->size[0];
  u_s->size[0] = gbNorm;
  emxEnsureCapacity_real_T(u_s, i);
  emxFree_char_T(&e_model_owensfile);
  emxFree_char_T(&d_model_owensfile);
  emxFree_char_T(&c_model_owensfile);
  emxFree_char_T(&b_model_owensfile);
  emxFree_char_T(&c_model_aeroloadfile);
  emxFree_char_T(&b_model_aeroloadfile);
  for (i = 0; i < gbNorm; i++) {
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
  udot_s->size[0] = gbNorm;
  emxEnsureCapacity_real_T(udot_s, i);
  for (i = 0; i < gbNorm; i++) {
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

  // ............................................................
  // define number of time steps
  // define time step size
  // store initial condition
  // initialize omega_platform, omega_platform_dot, omegaPlatHist
  //  omega_platform = zeros(3,1);
  //  omega_platform_dot = zeros(3,1);
  //  omegaPlatHist(:,1) = omega_platform;
  std::memset(&t_data[0], 0, 251U * sizeof(double));

  //  strainHist(numTS+1) = struct();
  // genTorque = zeros(1,numTS+1);
  t_data[0] = 0.0;

  // initialize various states and variables
  azi_s = 0.0;

  //  azi_sm1 = -Omega*delta_t*2*pi;
  //
  tic();

  //  structural dynamics initialization
  // ..........................................................................
  emxInit_struct_T1(&elStorage, 2);
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
  b_i = 0;
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
  while ((!exitg1) && (b_i < 250)) {
    //      i %TODO add verbose printing
    if (std::fmod(static_cast<double>(b_i) + 1.0, 100.0) == 0.0) {
      // print command that displays progress of time stepping
      // TODO: add string conversion of i
      printf("%s\n", "Iteration: ");
      fflush(stdout);
    }

    //     %% check for specified rotor speed at t(i) + delta_t
    // use discreteized rotor speed profile function
    if (1.5 < t_data[b_i] + 0.002) {
      exitg1 = true;
    } else {
      platNorm = 0;
      do {
        exitg2 = 0;
        if (platNorm < 2) {
          if (rtIsNaN(model_tocp[platNorm])) {
            exitg2 = 1;
          } else {
            platNorm++;
          }
        } else {
          omegaCurrent = rtNaN;
          if (!(t_data[b_i] + 0.002 > 1.5)) {
            omegaCurrent = 0.12000000000000001;
          }

          exitg2 = 1;
        }
      } while (exitg2 == 0);

      // interpolated discreteized profile for current omega
      // calculate current rotor acceleration
      platNorm = 0;
      do {
        exitg2 = 0;
        if (platNorm < 2) {
          if (rtIsNaN(model_tocp[platNorm])) {
            exitg2 = 1;
          } else {
            platNorm++;
          }
        } else {
          Vq = rtNaN;
          if (!((t_data[b_i] + 0.002) - 0.001 > 1.5)) {
            Vq = 0.12000000000000001;
          }

          exitg2 = 1;
        }
      } while (exitg2 == 0);

      platNorm = 0;
      do {
        exitg2 = 0;
        if (platNorm < 2) {
          if (rtIsNaN(model_tocp[platNorm])) {
            exitg2 = 1;
          } else {
            platNorm++;
          }
        } else {
          b_Vq = rtNaN;
          if (!((t_data[b_i] + 0.002) + 0.001 > 1.5)) {
            b_Vq = 0.12000000000000001;
          }

          exitg2 = 1;
        }
      } while (exitg2 == 0);

      OmegaDotCurrent = (b_Vq - Vq) / 0.002;

      //     %%
      //     %% initialize "j" Gauss-Sidel iteration
      i = u_j->size[0];
      u_j->size[0] = u_s->size[0];
      emxEnsureCapacity_real_T(u_j, i);
      gbNorm = u_s->size[0];
      for (i = 0; i < gbNorm; i++) {
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
        for (i = 0; i < 9; i++) {
          CN2P[i] = 0;
        }

        CN2P[0] = 1;
        CN2P[4] = 1;
        CN2P[8] = 1;

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
        // ==================================================
        //         %% rotor speed update
        // ===== update rotor speed =========================
        azi_jLast = azi_j;
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
        externalForcing(t_data[b_i] + 0.002, aeroLoads_timeArray_data,
                        aeroLoads_timeArray_size, aeroLoads_ForceValHist,
                        aeroLoads_ForceDof, Fexternal_sub, udot_j);
        if (Fexternal_sub->size[1] != 0) {
          sizes_idx_1 = Fexternal_sub->size[1];
        } else {
          sizes_idx_1 = 0;
        }

        if ((sizes_idx_1 == 0) || (Fexternal_sub->size[1] != 0)) {
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
          i = expl_temp.displddot_s->size[0];
          expl_temp.displddot_s->size[0] = uddot_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displddot_s, i);
          gbNorm = uddot_s->size[0];
          for (i = 0; i < gbNorm; i++) {
            expl_temp.displddot_s->data[i] = uddot_s->data[i];
          }

          i = expl_temp.displdot_s->size[0];
          expl_temp.displdot_s->size[0] = udot_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displdot_s, i);
          gbNorm = udot_s->size[0];
          for (i = 0; i < gbNorm; i++) {
            expl_temp.displdot_s->data[i] = udot_s->data[i];
          }

          i = expl_temp.displ_s->size[0];
          expl_temp.displ_s->size[0] = u_s->size[0];
          emxEnsureCapacity_real_T(expl_temp.displ_s, i);
          gbNorm = u_s->size[0];
          for (i = 0; i < gbNorm; i++) {
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
          b_Fexternal_sub->size[1] = sizes_idx_1;
          emxEnsureCapacity_real_T(b_Fexternal_sub, i);
          for (i = 0; i < sizes_idx_1; i++) {
            for (i1 = 0; i1 < platNorm; i1++) {
              b_Fexternal_sub->data[b_Fexternal_sub->size[0] * i] =
                Fexternal_sub->data[input_sizes_idx_0 * i];
            }
          }

          i = b_udot_j->size[0] * b_udot_j->size[1];
          b_udot_j->size[0] = sizes_idx_0;
          b_udot_j->size[1] = input_sizes_idx_1;
          emxEnsureCapacity_real_T(b_udot_j, i);
          gbNorm = input_sizes_idx_1;
          for (i = 0; i < gbNorm; i++) {
            for (i1 = 0; i1 < sizes_idx_0; i1++) {
              b_udot_j->data[i1] = udot_j->data[i1];
            }
          }

          for (i = 0; i < 3; i++) {
            aziNorm = CP2H_tmp[i + 3];
            i1 = static_cast<int>(CP2H_tmp[i + 6]);
            for (platNorm = 0; platNorm < 3; platNorm++) {
              b_CP2H_tmp[i + 3 * platNorm] = (CP2H_tmp[i] * static_cast<double>
                (CN2P[3 * platNorm]) + aziNorm * static_cast<double>(CN2P[3 *
                platNorm + 1])) + static_cast<double>((i1 * CN2P[3 * platNorm +
                2]));
            }
          }

          structuralDynamicsTransient(model_analysisType, model_RayleighAlpha,
            model_RayleighBeta, model_BC_numpBC, model_BC_pBC, model_joint,
            model_jointTransform, mesh.numEl, mesh.x, mesh.y, mesh.z, mesh.conn,
            el, expl_temp, omegaCurrent, OmegaDotCurrent, elStorage,
            b_Fexternal_sub, b_udot_j, b_CP2H_tmp, &dispOut, FReactionsm1);
          i = dispOut_displ_sp1->size[0];
          dispOut_displ_sp1->size[0] = dispOut.displ_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displ_sp1, i);
          gbNorm = dispOut.displ_sp1->size[0];
          for (i = 0; i < gbNorm; i++) {
            dispOut_displ_sp1->data[i] = dispOut.displ_sp1->data[i];
          }

          i = dispOut_displddot_sp1->size[0];
          dispOut_displddot_sp1->size[0] = dispOut.displddot_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displddot_sp1, i);
          gbNorm = dispOut.displddot_sp1->size[0];
          for (i = 0; i < gbNorm; i++) {
            dispOut_displddot_sp1->data[i] = dispOut.displddot_sp1->data[i];
          }

          i = dispOut_displdot_sp1->size[0];
          dispOut_displdot_sp1->size[0] = dispOut.displdot_sp1->size[0];
          emxEnsureCapacity_real_T(dispOut_displdot_sp1, i);
          gbNorm = dispOut.displdot_sp1->size[0];
          for (i = 0; i < gbNorm; i++) {
            dispOut_displdot_sp1->data[i] = dispOut.displdot_sp1->data[i];
          }
        }

        // update last iteration displacement vector
        i = u_jLast->size[0];
        u_jLast->size[0] = u_j->size[0];
        emxEnsureCapacity_real_T(u_jLast, i);
        gbNorm = u_j->size[0];
        for (i = 0; i < gbNorm; i++) {
          u_jLast->data[i] = u_j->data[i];
        }

        i = u_j->size[0];
        u_j->size[0] = dispOut_displ_sp1->size[0];
        emxEnsureCapacity_real_T(u_j, i);
        gbNorm = dispOut_displ_sp1->size[0];
        for (i = 0; i < gbNorm; i++) {
          u_j->data[i] = dispOut_displ_sp1->data[i];
        }

        // update current estimates of velocity, acceleration
        //  Only used for TNB, but must be declared
        i = udot_j->size[0];
        udot_j->size[0] = dispOut_displdot_sp1->size[0];
        emxEnsureCapacity_real_T(udot_j, i);
        gbNorm = dispOut_displdot_sp1->size[0];
        for (i = 0; i < gbNorm; i++) {
          udot_j->data[i] = dispOut_displdot_sp1->data[i];
        }

        i = uddot_j->size[0];
        uddot_j->size[0] = dispOut_displddot_sp1->size[0];
        emxEnsureCapacity_real_T(uddot_j, i);
        gbNorm = dispOut_displddot_sp1->size[0];
        for (i = 0; i < gbNorm; i++) {
          uddot_j->data[i] = dispOut_displddot_sp1->data[i];
        }

        //         %%
        //         %% calculate norms
        i = u_jLast->size[0];
        u_jLast->size[0] = dispOut_displ_sp1->size[0];
        emxEnsureCapacity_real_T(u_jLast, i);
        gbNorm = dispOut_displ_sp1->size[0];
        for (i = 0; i < gbNorm; i++) {
          u_jLast->data[i] = dispOut_displ_sp1->data[i] - u_jLast->data[i];
        }

        uNorm = b_norm(u_jLast) / b_norm(dispOut_displ_sp1);

        // structural dynamics displacement iteration norm
        aziNorm = std::abs(azi_j - azi_jLast) / std::abs(azi_j);

        // rotor azimuth iteration norm
        platNorm = 0;
        gbNorm = 0;
        numIterations++;
      }

      // end iteration while loop
      //     %% calculate converged generator torque/power
      //     %% update timestepping variables and other states, store in history arrays 
      if (c_strcmp(model_analysisType)) {
        i = u_s->size[0];
        u_s->size[0] = u_j->size[0];
        emxEnsureCapacity_real_T(u_s, i);
        gbNorm = u_j->size[0];
        for (i = 0; i < gbNorm; i++) {
          u_s->data[i] = u_j->data[i];
        }

        i = udot_s->size[0];
        udot_s->size[0] = udot_j->size[0];
        emxEnsureCapacity_real_T(udot_s, i);
        gbNorm = udot_j->size[0];
        for (i = 0; i < gbNorm; i++) {
          udot_s->data[i] = udot_j->data[i];
        }

        i = uddot_s->size[0];
        uddot_s->size[0] = uddot_j->size[0];
        emxEnsureCapacity_real_T(uddot_s, i);
        gbNorm = uddot_j->size[0];
        for (i = 0; i < gbNorm; i++) {
          uddot_s->data[i] = uddot_j->data[i];
        }
      }

      t_data[b_i + 1] = t_data[b_i] + 0.002;
      azi_s = azi_j;

      // genTorque(i+1) = genTorque_s;
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
  emxFree_struct_T1(&elStorage);
  emxFree_real_T(&uddot_s);
  emxFree_real_T(&udot_s);
  emxFree_real_T(&u_s);
  emxFree_real_T(&dispOut_displdot_sp1);
  emxFree_real_T(&dispOut_displddot_sp1);
  emxFree_real_T(&dispOut_displ_sp1);
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
  printf("%s\n", "NOT WRITING Verification File");
  fflush(stdout);
}

//
// File trailer for transientExec.cpp
//
// [EOF]
//
