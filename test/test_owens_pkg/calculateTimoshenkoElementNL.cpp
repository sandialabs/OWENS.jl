//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateTimoshenkoElementNL.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "calculateTimoshenkoElementNL.h"
#include "calculateLambda.h"
#include "calculateShapeFunctions.h"
#include "find.h"
#include "mtimes.h"
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "strcmp.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Declarations
static void mapMatrixNonSym(const double Ktemp[144], double Kel[144]);

// Function Definitions

//
// ----- function to form total stifness matrix and transform to desired
//  DOF mapping
// Arguments    : const double Ktemp[144]
//                double Kel[144]
// Return Type  : void
//
static void mapMatrixNonSym(const double Ktemp[144], double Kel[144])
{
  emxArray_real_T *y_d;
  int i;
  emxArray_int32_T *y_colidx;
  emxArray_int32_T *y_rowidx;
  emxArray_real_T *b_y_d;
  int ctr;
  int col;
  emxArray_int32_T *b_y_colidx;
  int cend;
  static const signed char x[144] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };

  emxArray_int32_T *b_y_rowidx;
  emxArray_real_T *t1_d;
  emxArray_int32_T *t1_colidx;
  emxArray_int32_T *t1_rowidx;
  static const signed char b_x[144] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };

  emxArray_real_T *t2_d;
  emxArray_int32_T *t2_colidx;
  emxArray_int32_T *t2_rowidx;
  emxInit_real_T(&y_d, 1);

  // map to FEA numbering
  i = y_d->size[0];
  y_d->size[0] = 12;
  emxEnsureCapacity_real_T(y_d, i);
  for (i = 0; i < 12; i++) {
    y_d->data[i] = 0.0;
  }

  emxInit_int32_T(&y_colidx, 1);
  i = y_colidx->size[0];
  y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(y_colidx, i);
  for (i = 0; i < 13; i++) {
    y_colidx->data[i] = 0;
  }

  emxInit_int32_T(&y_rowidx, 1);
  y_colidx->data[0] = 1;
  i = y_rowidx->size[0];
  y_rowidx->size[0] = 12;
  emxEnsureCapacity_int32_T(y_rowidx, i);
  for (i = 0; i < 12; i++) {
    y_rowidx->data[i] = 0;
  }

  emxInit_real_T(&b_y_d, 1);
  y_rowidx->data[0] = 1;
  ctr = 0;
  i = b_y_d->size[0];
  b_y_d->size[0] = 12;
  emxEnsureCapacity_real_T(b_y_d, i);
  for (col = 0; col < 12; col++) {
    for (cend = 0; cend < 12; cend++) {
      if (x[cend + 12 * col] != 0) {
        y_rowidx->data[ctr] = cend + 1;
        y_d->data[ctr] = 1.0;
        ctr++;
      }
    }

    y_colidx->data[col + 1] = ctr + 1;
    b_y_d->data[col] = 0.0;
  }

  emxInit_int32_T(&b_y_colidx, 1);
  i = b_y_colidx->size[0];
  b_y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(b_y_colidx, i);
  for (i = 0; i < 13; i++) {
    b_y_colidx->data[i] = 0;
  }

  emxInit_int32_T(&b_y_rowidx, 1);
  b_y_colidx->data[0] = 1;
  i = b_y_rowidx->size[0];
  b_y_rowidx->size[0] = 12;
  emxEnsureCapacity_int32_T(b_y_rowidx, i);
  for (i = 0; i < 12; i++) {
    b_y_rowidx->data[i] = 0;
  }

  b_y_rowidx->data[0] = 1;
  ctr = 0;
  for (col = 0; col < 12; col++) {
    for (cend = 0; cend < 12; cend++) {
      if (b_x[cend + 12 * col] != 0) {
        b_y_rowidx->data[ctr] = cend + 1;
        b_y_d->data[ctr] = 1.0;
        ctr++;
      }
    }

    b_y_colidx->data[col + 1] = ctr + 1;
  }

  emxInit_real_T(&t1_d, 1);
  emxInit_int32_T(&t1_colidx, 1);
  emxInit_int32_T(&t1_rowidx, 1);
  emxInit_real_T(&t2_d, 1);
  emxInit_int32_T(&t2_colidx, 1);
  emxInit_int32_T(&t2_rowidx, 1);
  sparse(Ktemp, t1_d, t1_colidx, t1_rowidx);
  sparse_mtimes(y_d, y_colidx, y_rowidx, t1_d, t1_colidx, t1_rowidx, t2_d,
                t2_colidx, t2_rowidx);
  sparse_mtimes(t2_d, t2_colidx, t2_rowidx, b_y_d, b_y_colidx, b_y_rowidx, y_d,
                y_colidx, y_rowidx);
  emxFree_int32_T(&t2_rowidx);
  emxFree_int32_T(&t2_colidx);
  emxFree_real_T(&t2_d);
  emxFree_int32_T(&t1_rowidx);
  emxFree_int32_T(&t1_colidx);
  emxFree_real_T(&t1_d);
  emxFree_int32_T(&b_y_rowidx);
  emxFree_int32_T(&b_y_colidx);
  emxFree_real_T(&b_y_d);
  std::memset(&Kel[0], 0, 144U * sizeof(double));
  for (col = 0; col < 12; col++) {
    cend = y_colidx->data[col + 1] - 1;
    i = y_colidx->data[col];
    for (ctr = i; ctr <= cend; ctr++) {
      Kel[(y_rowidx->data[ctr - 1] + 12 * col) - 1] = y_d->data[ctr - 1];
    }
  }

  emxFree_int32_T(&y_rowidx);
  emxFree_int32_T(&y_colidx);
  emxFree_real_T(&y_d);

  // declare map
  //  map = [1, 7, 2, 8, 3, 9,...
  //        4, 10, 5, 11, 6, 12];
  //
  //  %map to FEA numbering
  //  for i=1:a
  //      I=map(i);
  //      for j=1:a
  //          J=map(j);
  //          Kel(I,J) = Ktemp(i,j);
  //      end
  //  end
}

//
// calculateTimoshenkoElementNL performs nonlinear element calculations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [output] = calculateTimoshenkoElementNL(input,elStorage)
//
//    This function performs nonlinear element calculations.
//
//       input:
//       input      = object containing element input
//       elStorage  = obect containing precalculated element data
//
//       output:
//       output     = object containing element data
// Arguments    : const p_struct_T *input
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void b_calculateTimoshenkoElementNL(const p_struct_T *input, const f_struct_T
  *elStorage, o_struct_T *output)
{
  double disp_iter_data[12];
  double dispm1_data[12];
  int dispdot_size_idx_1;
  double dispdot_data[12];
  int dispddot_size_idx_1;
  double dispddot_data[12];
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double Omega;
  double OmegaDot;
  double lambda[144];
  int i;
  double O1;
  double Oel[3];
  double K15_idx_0;
  double O2;
  double H34_idx_3;
  double O3;
  double O1dot;
  double ODotel[3];
  double a_temp[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double H34_idx_0;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double integrationFactor;
  double c_tmp;
  double b_c_tmp;
  double S12_idx_0;
  double S13_idx_0;
  double S23_idx_0;
  double S25_idx_0;
  double S26_idx_0;
  double S35_idx_0;
  double S14_idx_0;
  double S24_idx_0;
  double S34_idx_0;
  double C34_idx_0;
  double H12_idx_0;
  double zcm;
  double rhoA;
  double ycm;
  double H35_idx_0;
  double H36_idx_0;
  double H14_idx_0;
  double disMomentgp[3];
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double posLocal[3];
  double S35_idx_1;
  double disLoadgpLocal[3];
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C34_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double H23_idx_1;
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C34_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C34_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double H46_idx_3;
  double ktemp2[144];
  double Kenr[144];
  double c_c_tmp;
  double K15_idx_1;
  double K16_idx_1;
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double ktemp2_tmp;
  double b_ktemp2_tmp;
  double c_ktemp2_tmp;
  double d_ktemp2_tmp;
  double e_ktemp2_tmp;
  double f_ktemp2_tmp;
  double g_ktemp2_tmp;
  double Khate[144];
  double Ce[144];
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fhate[12];
  double Fe[12];
  int tmp_size[1];
  boolean_T concMassFlag;
  double FhatLessConc[12];
  double b_Khate[12];
  double b_data[12];

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  std::memcpy(&disp_iter_data[0], &input->displ_iter.data[0], 12U * sizeof
              (double));
  dispm1_data[0] = 0.0;

  // declare type
  dispdot_size_idx_1 = 1;
  dispdot_data[0] = 0.0;

  // declare type
  dispddot_size_idx_1 = 1;
  dispddot_data[0] = 0.0;

  // declare type
  // options for Dean integrator
  if (g_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&dispm1_data[0], &input->dispm1.data[0], 12U * sizeof(double));
  } else {
    if (h_strcmp(input->analysisType.data, input->analysisType.size)) {
      // options for newmark beta integrator
      dispdot_size_idx_1 = 12;
      dispddot_size_idx_1 = 12;
      std::memcpy(&dispdot_data[0], &input->dispdot.data[0], 12U * sizeof(double));
      std::memcpy(&dispddot_data[0], &input->dispddot.data[0], 12U * sizeof
                  (double));
    }
  }

  // --------------------------------------------
  // setting for modal analysis flag
  if (i_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&disp_iter_data[0], &input->disp.data[0], 12U * sizeof(double));
  }

  // setting for initial reduced order model calculations
  if (j_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&disp_iter_data[0], &input->disp.data[0], 12U * sizeof(double));
  }

  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  Omega = 6.2831853071795862 * input->Omega;
  OmegaDot = 6.2831853071795862 * input->OmegaDot;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input->sweepAngle * 3.1415926535897931 / 180.0,
                  input->coneAngle * 3.1415926535897931 / 180.0,
                  (input->rollAngle + 0.5 * (input->sectionProps.twist[0] +
    input->sectionProps.twist[1])) * 3.1415926535897931 / 180.0, lambda);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    K15_idx_0 = lambda[i + 24];
    H34_idx_3 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    Oel[i] = H34_idx_3 + K15_idx_0 * Omega;
    a_temp[i] = (input->CN2H[i] * 0.0 + input->CN2H[i + 3] * 0.0) + input->
      CN2H[i + 6] * 9.81;
    ODotel[i] = H34_idx_3 + K15_idx_0 * OmegaDot;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  H34_idx_0 = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input->xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &H34_idx_3);
    integrationFactor = H34_idx_3 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input->sectionProps.rhoA[0] + N_data[1] *
      input->sectionProps.rhoA[1];
    ycm = N_data[0] * input->sectionProps.ycm[0] + N_data[1] *
      input->sectionProps.ycm[1];
    zcm = N_data[0] * input->sectionProps.zcm[0] + N_data[1] *
      input->sectionProps.zcm[1];
    S14_idx_0 = N_data[0] * input->x.data[0] + N_data[1] * input->x.data[1];
    S12_idx_0 = N_data[0] * input->y.data[0] + N_data[1] * input->y.data[1];
    H34_idx_3 = N_data[0] * input->z.data[0] + N_data[1] * input->z.data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      K15_idx_0 = lambda[i + 12];
      S23_idx_0 = lambda[i + 24];
      posLocal[i] = (lambda[i] * S14_idx_0 + K15_idx_0 * S12_idx_0) + S23_idx_0 *
        H34_idx_3;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + K15_idx_0 * disMomentgp
                           [1]) + S23_idx_0 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // calculate static aerodynamic load
    // distributed/body force load calculations
    K15_idx_0 = (O2 * O2 + O3 * O3) * posLocal[0];
    H12_idx_0 = O2dot * posLocal[2];
    S24_idx_0 = O3dot * posLocal[1];
    S25_idx_0 = H34_idx_0 * posLocal[2];
    S26_idx_0 = O1 * O2 * posLocal[1];
    S13_idx_0 = rhoA * ((((K15_idx_0 - S26_idx_0) - S25_idx_0) + S24_idx_0) -
                        H12_idx_0) - disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    H34_idx_3 = O1dot * posLocal[2];
    S23_idx_0 = O3dot * posLocal[0];
    S35_idx_0 = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 *
      O3) - posLocal[0] * O1 * O2) + H34_idx_3) - S23_idx_0) - disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    S14_idx_0 = O2dot * posLocal[0];
    S12_idx_0 = O1dot * posLocal[1];
    S34_idx_0 = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - H34_idx_0 *
      posLocal[0]) - O2 * O3 * posLocal[1]) + S14_idx_0) - S12_idx_0) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    S14_idx_0 = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) -
      posLocal[1] * (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                          (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                         (S14_idx_0 - S12_idx_0)) - zcm * (H34_idx_3 - S23_idx_0))
      - disMomentgp[0];

    // This function is a general routine to calculate an element vector
    H34_idx_3 = rhoA * zcm * ((((K15_idx_0 - posLocal[1] * O1 * O2) - posLocal[2]
      * O1 * O3) - H12_idx_0) + S24_idx_0) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    S23_idx_0 = rhoA * ycm * ((((S25_idx_0 + S26_idx_0) - K15_idx_0) - S24_idx_0)
      + H12_idx_0) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    F1_data_idx_0 += S13_idx_0 * N_data[0] * integrationFactor;
    F2_data_idx_0 += S35_idx_0 * N_data[0] * integrationFactor;
    F3_data_idx_0 += S34_idx_0 * N_data[0] * integrationFactor;
    F4_data_idx_0 += S14_idx_0 * N_data[0] * integrationFactor;
    F5_data_idx_0 += H34_idx_3 * N_data[0] * integrationFactor;
    F6_data_idx_0 += S23_idx_0 * N_data[0] * integrationFactor;
    F1_data_idx_1 += S13_idx_0 * N_data[1] * integrationFactor;
    F2_data_idx_1 += S35_idx_0 * N_data[1] * integrationFactor;
    F3_data_idx_1 += S34_idx_0 * N_data[1] * integrationFactor;
    F4_data_idx_1 += S14_idx_0 * N_data[1] * integrationFactor;
    F5_data_idx_1 += H34_idx_3 * N_data[1] * integrationFactor;
    F6_data_idx_1 += S23_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  H34_idx_3 = Oel[0] * Oel[1];
  c_tmp = Oel[0] * Oel[0];
  b_c_tmp = c_tmp + Oel[2] * Oel[2];
  c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  S12_idx_0 = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  S23_idx_0 = elStorage->S23[0] * Oel[1] * Oel[2];
  S25_idx_0 = elStorage->S25[0] * H34_idx_3;
  S26_idx_0 = elStorage->S26[0] * H34_idx_3;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  O1 = elStorage->S36[0] * Oel[0] * Oel[2];
  S14_idx_0 = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel
    [0] * Oel[1];
  S24_idx_0 = elStorage->S24_1[0] * b_c_tmp + elStorage->S24_2[0] * Oel[1] *
    Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * c_tmp + elStorage->S34_2[0] * Oel[1] * Oel[2];
  O2 = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel[0] *
    Oel[1];
  O3 = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0] * Oel[0] *
    Oel[2];
  K15_idx_0 = elStorage->C34[0];
  C34_idx_0 = K15_idx_0 * Oel[0];
  H12_idx_0 = 0.5 * elStorage->C12[0] * ODotel[2];
  zcm = 0.5 * elStorage->C13[0] * ODotel[1];
  rhoA = 0.5 * elStorage->C23[0] * ODotel[0];
  O1dot = 0.5 * elStorage->C24[0] * ODotel[0];
  O2dot = 0.5 * elStorage->C25[0] * ODotel[2];
  O3dot = 0.5 * elStorage->C26[0] * ODotel[2];
  H34_idx_0 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_0 = 0.5 * elStorage->C35[0] * ODotel[1];
  H36_idx_0 = 0.5 * elStorage->C36[0] * ODotel[1];
  H14_idx_0 = 0.5 * (elStorage->C14_1[0] * ODotel[1] + elStorage->C14_2[0] *
                     ODotel[2]);
  H45_idx_0 = 0.5 * (elStorage->C45_1[0] * ODotel[2] + elStorage->C45_2[0] *
                     ODotel[1]);
  H46_idx_0 = 0.5 * (elStorage->C46_1[0] * ODotel[1] + elStorage->C46_2[0] *
                     ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * H34_idx_3;
  S26_idx_1 = elStorage->S26[1] * H34_idx_3;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * b_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * c_tmp + elStorage->S34_2[1] * Oel[1] * Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[1];
  C34_idx_1 = K15_idx_0 * Oel[0];
  H12_idx_1 = 0.5 * elStorage->C12[1] * ODotel[2];
  H13_idx_1 = 0.5 * elStorage->C13[1] * ODotel[1];
  H23_idx_1 = 0.5 * elStorage->C23[1] * ODotel[0];
  H24_idx_1 = 0.5 * elStorage->C24[1] * ODotel[0];
  H25_idx_1 = 0.5 * elStorage->C25[1] * ODotel[2];
  H26_idx_1 = 0.5 * elStorage->C26[1] * ODotel[2];
  H34_idx_1 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_1 = 0.5 * elStorage->C35[1] * ODotel[1];
  H36_idx_1 = 0.5 * elStorage->C36[1] * ODotel[1];
  H14_idx_1 = 0.5 * (elStorage->C14_1[1] * ODotel[1] + elStorage->C14_2[1] *
                     ODotel[2]);
  H45_idx_1 = 0.5 * (elStorage->C45_1[1] * ODotel[2] + elStorage->C45_2[1] *
                     ODotel[1]);
  H46_idx_1 = 0.5 * (elStorage->C46_1[1] * ODotel[1] + elStorage->C46_2[1] *
                     ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * H34_idx_3;
  S26_idx_2 = elStorage->S26[2] * H34_idx_3;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * b_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * c_tmp + elStorage->S34_2[2] * Oel[1] * Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[2];
  C34_idx_2 = K15_idx_0 * Oel[0];
  H12_idx_2 = 0.5 * elStorage->C12[2] * ODotel[2];
  H13_idx_2 = 0.5 * elStorage->C13[2] * ODotel[1];
  H23_idx_2 = 0.5 * elStorage->C23[2] * ODotel[0];
  H24_idx_2 = 0.5 * elStorage->C24[2] * ODotel[0];
  H25_idx_2 = 0.5 * elStorage->C25[2] * ODotel[2];
  H26_idx_2 = 0.5 * elStorage->C26[2] * ODotel[2];
  H34_idx_2 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_2 = 0.5 * elStorage->C35[2] * ODotel[1];
  H36_idx_2 = 0.5 * elStorage->C36[2] * ODotel[1];
  H14_idx_2 = 0.5 * (elStorage->C14_1[2] * ODotel[1] + elStorage->C14_2[2] *
                     ODotel[2]);
  H45_idx_2 = 0.5 * (elStorage->C45_1[2] * ODotel[2] + elStorage->C45_2[2] *
                     ODotel[1]);
  H46_idx_2 = 0.5 * (elStorage->C46_1[2] * ODotel[1] + elStorage->C46_2[2] *
                     ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * H34_idx_3;
  S26_idx_3 = elStorage->S26[3] * H34_idx_3;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * b_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * c_tmp + elStorage->S34_2[3] * Oel[1] * Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[3];
  C34_idx_3 = K15_idx_0 * Oel[0];
  H12_idx_3 = 0.5 * elStorage->C12[3] * ODotel[2];
  H13_idx_3 = 0.5 * elStorage->C13[3] * ODotel[1];
  H23_idx_3 = 0.5 * elStorage->C23[3] * ODotel[0];
  H24_idx_3 = 0.5 * elStorage->C24[3] * ODotel[0];
  H25_idx_3 = 0.5 * elStorage->C25[3] * ODotel[2];
  H26_idx_3 = 0.5 * elStorage->C26[3] * ODotel[2];
  H34_idx_3 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_3 = 0.5 * elStorage->C35[3] * ODotel[1];
  H36_idx_3 = 0.5 * elStorage->C36[3] * ODotel[1];
  H14_idx_3 = 0.5 * (elStorage->C14_1[3] * ODotel[1] + elStorage->C14_2[3] *
                     ODotel[2]);
  H45_idx_3 = 0.5 * (elStorage->C45_1[3] * ODotel[2] + elStorage->C45_2[3] *
                     ODotel[1]);
  H46_idx_3 = 0.5 * (elStorage->C46_1[3] * ODotel[1] + elStorage->C46_2[3] *
                     ODotel[2]);

  // compile stiffness matrix without rotational effects
  ktemp2[0] = elStorage->K11[0];
  ktemp2[24] = elStorage->K12[0];
  ktemp2[48] = elStorage->K13[0];
  ktemp2[72] = elStorage->K14[0];
  ktemp2[96] = elStorage->K15[0];
  ktemp2[120] = elStorage->K16[0];
  ktemp2[2] = elStorage->K12[0];
  ktemp2[26] = elStorage->K22[0];
  ktemp2[50] = elStorage->K23[0];
  ktemp2[74] = elStorage->K24[0];
  ktemp2[98] = elStorage->K25[0];
  ktemp2[122] = elStorage->K26[0];
  ktemp2[4] = elStorage->K13[0];
  ktemp2[28] = elStorage->K23[0];
  ktemp2[52] = elStorage->K33[0];
  ktemp2[76] = elStorage->K34[0];
  ktemp2[100] = elStorage->K35[0];
  ktemp2[124] = elStorage->K36[0];
  ktemp2[6] = elStorage->K13[0];
  ktemp2[30] = elStorage->K24[0];
  ktemp2[54] = elStorage->K34[0];
  ktemp2[78] = elStorage->K44[0];
  ktemp2[102] = elStorage->K45[0];
  ktemp2[126] = elStorage->K46[0];
  ktemp2[8] = elStorage->K15[0];
  ktemp2[32] = elStorage->K25[0];
  ktemp2[56] = elStorage->K35[0];
  ktemp2[80] = elStorage->K45[0];
  ktemp2[104] = elStorage->K55[0];
  ktemp2[128] = elStorage->K56[0];
  ktemp2[10] = elStorage->K16[0];
  ktemp2[34] = elStorage->K26[0];
  ktemp2[58] = elStorage->K36[0];
  ktemp2[82] = elStorage->K46[0];
  ktemp2[106] = elStorage->K56[0];
  ktemp2[130] = elStorage->K66[0];
  ktemp2[1] = elStorage->K11[1];
  ktemp2[25] = elStorage->K12[1];
  ktemp2[49] = elStorage->K13[1];
  ktemp2[73] = elStorage->K14[1];
  ktemp2[97] = elStorage->K15[1];
  ktemp2[121] = elStorage->K16[1];
  ktemp2[3] = elStorage->K12[2];
  ktemp2[27] = elStorage->K22[1];
  ktemp2[51] = elStorage->K23[1];
  ktemp2[75] = elStorage->K24[1];
  ktemp2[99] = elStorage->K25[1];
  ktemp2[123] = elStorage->K26[1];
  ktemp2[5] = elStorage->K13[2];
  ktemp2[29] = elStorage->K23[2];
  ktemp2[53] = elStorage->K33[1];
  ktemp2[77] = elStorage->K34[1];
  ktemp2[101] = elStorage->K35[1];
  ktemp2[125] = elStorage->K36[1];
  ktemp2[7] = elStorage->K13[2];
  ktemp2[31] = elStorage->K24[2];
  ktemp2[55] = elStorage->K34[2];
  ktemp2[79] = elStorage->K44[1];
  ktemp2[103] = elStorage->K45[1];
  ktemp2[127] = elStorage->K46[1];
  ktemp2[9] = elStorage->K15[2];
  ktemp2[33] = elStorage->K25[2];
  ktemp2[57] = elStorage->K35[2];
  ktemp2[81] = elStorage->K45[2];
  ktemp2[105] = elStorage->K55[1];
  ktemp2[129] = elStorage->K56[1];
  ktemp2[11] = elStorage->K16[2];
  ktemp2[35] = elStorage->K26[2];
  ktemp2[59] = elStorage->K36[2];
  ktemp2[83] = elStorage->K46[2];
  ktemp2[107] = elStorage->K56[2];
  ktemp2[131] = elStorage->K66[1];
  ktemp2[12] = elStorage->K11[2];
  ktemp2[36] = elStorage->K12[2];
  ktemp2[60] = elStorage->K13[2];
  ktemp2[84] = elStorage->K14[2];
  ktemp2[108] = elStorage->K15[2];
  ktemp2[132] = elStorage->K16[2];
  ktemp2[14] = elStorage->K12[1];
  ktemp2[38] = elStorage->K22[2];
  ktemp2[62] = elStorage->K23[2];
  ktemp2[86] = elStorage->K24[2];
  ktemp2[110] = elStorage->K25[2];
  ktemp2[134] = elStorage->K26[2];
  ktemp2[16] = elStorage->K13[1];
  ktemp2[40] = elStorage->K23[1];
  ktemp2[64] = elStorage->K33[2];
  ktemp2[88] = elStorage->K34[2];
  ktemp2[112] = elStorage->K35[2];
  ktemp2[136] = elStorage->K36[2];
  ktemp2[18] = elStorage->K13[1];
  ktemp2[42] = elStorage->K24[1];
  ktemp2[66] = elStorage->K34[1];
  ktemp2[90] = elStorage->K44[2];
  ktemp2[114] = elStorage->K45[2];
  ktemp2[138] = elStorage->K46[2];
  ktemp2[20] = elStorage->K15[1];
  ktemp2[44] = elStorage->K25[1];
  ktemp2[68] = elStorage->K35[1];
  ktemp2[92] = elStorage->K45[1];
  ktemp2[116] = elStorage->K55[2];
  ktemp2[140] = elStorage->K56[2];
  ktemp2[22] = elStorage->K16[1];
  ktemp2[46] = elStorage->K26[1];
  ktemp2[70] = elStorage->K36[1];
  ktemp2[94] = elStorage->K46[1];
  ktemp2[118] = elStorage->K56[1];
  ktemp2[142] = elStorage->K66[2];
  ktemp2[13] = elStorage->K11[3];
  ktemp2[37] = elStorage->K12[3];
  ktemp2[61] = elStorage->K13[3];
  ktemp2[85] = elStorage->K14[3];
  ktemp2[109] = elStorage->K15[3];
  ktemp2[133] = elStorage->K16[3];
  ktemp2[15] = elStorage->K12[3];
  ktemp2[39] = elStorage->K22[3];
  ktemp2[63] = elStorage->K23[3];
  ktemp2[87] = elStorage->K24[3];
  ktemp2[111] = elStorage->K25[3];
  ktemp2[135] = elStorage->K26[3];
  ktemp2[17] = elStorage->K13[3];
  ktemp2[41] = elStorage->K23[3];
  ktemp2[65] = elStorage->K33[3];
  ktemp2[89] = elStorage->K34[3];
  ktemp2[113] = elStorage->K35[3];
  ktemp2[137] = elStorage->K36[3];
  ktemp2[19] = elStorage->K13[3];
  ktemp2[43] = elStorage->K24[3];
  ktemp2[67] = elStorage->K34[3];
  ktemp2[91] = elStorage->K44[3];
  ktemp2[115] = elStorage->K45[3];
  ktemp2[139] = elStorage->K46[3];
  ktemp2[21] = elStorage->K15[3];
  ktemp2[45] = elStorage->K25[3];
  ktemp2[69] = elStorage->K35[3];
  ktemp2[93] = elStorage->K45[3];
  ktemp2[117] = elStorage->K55[3];
  ktemp2[141] = elStorage->K56[3];
  ktemp2[23] = elStorage->K16[3];
  ktemp2[47] = elStorage->K26[3];
  ktemp2[71] = elStorage->K36[3];
  ktemp2[95] = elStorage->K46[3];
  ktemp2[119] = elStorage->K56[3];
  ktemp2[143] = elStorage->K66[3];
  mapMatrixNonSym(ktemp2, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  c_c_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * c_c_tmp;
  ycm = elStorage->K16[0] + elStorage->S16[0] * c_c_tmp;
  integrationFactor = elStorage->K56[0] + elStorage->S56[0] * c_c_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * c_c_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * c_c_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * c_c_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * c_c_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * c_c_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * c_c_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * c_c_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * c_c_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * c_c_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  ktemp2[0] = elStorage->K11[0] + elStorage->S11[0] * c_c_tmp;
  ktemp2[24] = (elStorage->K12[0] + S12_idx_0) + H12_idx_0;
  ktemp2[48] = (elStorage->K13[0] + S13_idx_0) + zcm;
  ktemp2_tmp = elStorage->K14[0] + S14_idx_0;
  ktemp2[72] = ktemp2_tmp + H14_idx_0;
  ktemp2[96] = K15_idx_0;
  ktemp2[120] = ycm;
  ktemp2[2] = (elStorage->K12[0] + S12_idx_0) - H12_idx_0;
  ktemp2[26] = elStorage->K22[0] + elStorage->S22[0] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[0] + S23_idx_0;
  ktemp2[50] = b_ktemp2_tmp + rhoA;
  c_ktemp2_tmp = elStorage->K24[0] + S24_idx_0;
  ktemp2[74] = c_ktemp2_tmp + O1dot;
  d_ktemp2_tmp = elStorage->K25[0] + S25_idx_0;
  ktemp2[98] = d_ktemp2_tmp + O2dot;
  e_ktemp2_tmp = elStorage->K26[0] + S26_idx_0;
  ktemp2[122] = e_ktemp2_tmp + O3dot;
  ktemp2[4] = (elStorage->K13[0] + S13_idx_0) - zcm;
  ktemp2[28] = b_ktemp2_tmp - rhoA;
  ktemp2[52] = elStorage->K33[0] + elStorage->S33[0] * c_tmp;
  b_ktemp2_tmp = elStorage->K34[0] + S34_idx_0;
  ktemp2[76] = b_ktemp2_tmp + H34_idx_0;
  f_ktemp2_tmp = elStorage->K35[0] + S35_idx_0;
  ktemp2[100] = f_ktemp2_tmp + H35_idx_0;
  g_ktemp2_tmp = elStorage->K36[0] + O1;
  ktemp2[124] = g_ktemp2_tmp + H36_idx_0;
  ktemp2[6] = ktemp2_tmp - H14_idx_0;
  ktemp2[30] = c_ktemp2_tmp - O1dot;
  ktemp2[54] = b_ktemp2_tmp - H34_idx_0;
  ktemp2[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * b_c_tmp +
    elStorage->S44_2[0] * c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[0] + O2;
  ktemp2[102] = ktemp2_tmp + H45_idx_0;
  b_ktemp2_tmp = elStorage->K46[0] + O3;
  ktemp2[126] = b_ktemp2_tmp + H46_idx_0;
  ktemp2[8] = K15_idx_0;
  ktemp2[32] = d_ktemp2_tmp - O2dot;
  ktemp2[56] = f_ktemp2_tmp - H35_idx_0;
  ktemp2[80] = ktemp2_tmp - H45_idx_0;
  ktemp2[104] = elStorage->K55[0] + elStorage->S55[0] * c_c_tmp;
  ktemp2[128] = integrationFactor;
  ktemp2[10] = ycm;
  ktemp2[34] = e_ktemp2_tmp - O3dot;
  ktemp2[58] = g_ktemp2_tmp - H36_idx_0;
  ktemp2[82] = b_ktemp2_tmp - H46_idx_0;
  ktemp2[106] = integrationFactor;
  ktemp2[130] = elStorage->K66[0] + elStorage->S66[0] * c_c_tmp;
  ktemp2[1] = elStorage->K11[1] + elStorage->S11[1] * c_c_tmp;
  ktemp2[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  ktemp2[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  ktemp2_tmp = elStorage->K14[1] + S14_idx_1;
  ktemp2[73] = ktemp2_tmp + H14_idx_1;
  ktemp2[97] = K15_idx_1;
  ktemp2[121] = K16_idx_1;
  ktemp2[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  ktemp2[27] = elStorage->K22[1] + elStorage->S22[1] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[1] + S23_idx_1;
  ktemp2[51] = b_ktemp2_tmp + H23_idx_1;
  c_ktemp2_tmp = elStorage->K24[1] + S24_idx_1;
  ktemp2[75] = c_ktemp2_tmp + H24_idx_1;
  d_ktemp2_tmp = elStorage->K25[1] + S25_idx_1;
  ktemp2[99] = d_ktemp2_tmp + H25_idx_1;
  e_ktemp2_tmp = elStorage->K26[1] + S26_idx_1;
  ktemp2[123] = e_ktemp2_tmp + H26_idx_1;
  ktemp2[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  f_ktemp2_tmp = elStorage->K23[2] + S23_idx_2;
  ktemp2[29] = f_ktemp2_tmp - H23_idx_2;
  ktemp2[53] = elStorage->K33[1] + elStorage->S33[1] * c_tmp;
  g_ktemp2_tmp = elStorage->K34[1] + S34_idx_1;
  ktemp2[77] = g_ktemp2_tmp + H34_idx_1;
  O1 = elStorage->K35[1] + S35_idx_1;
  ktemp2[101] = O1 + H35_idx_1;
  integrationFactor = elStorage->K36[1] + S36_idx_1;
  ktemp2[125] = integrationFactor + H36_idx_1;
  ycm = elStorage->K14[2] + S14_idx_2;
  ktemp2[7] = ycm - H14_idx_2;
  rhoA = elStorage->K24[2] + S24_idx_2;
  ktemp2[31] = rhoA - H24_idx_2;
  zcm = elStorage->K34[2] + S34_idx_2;
  ktemp2[55] = zcm - H34_idx_2;
  ktemp2[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * b_c_tmp +
    elStorage->S44_2[1] * c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  S34_idx_0 = elStorage->K45[1] + S45_idx_1;
  ktemp2[103] = S34_idx_0 + H45_idx_1;
  S35_idx_0 = elStorage->K46[1] + S46_idx_1;
  ktemp2[127] = S35_idx_0 + H46_idx_1;
  ktemp2[9] = K15_idx_2;
  S13_idx_0 = elStorage->K25[2] + S25_idx_2;
  ktemp2[33] = S13_idx_0 - H25_idx_2;
  S26_idx_0 = elStorage->K35[2] + S35_idx_2;
  ktemp2[57] = S26_idx_0 - H35_idx_2;
  S25_idx_0 = elStorage->K45[2] + S45_idx_2;
  ktemp2[81] = S25_idx_0 - H45_idx_2;
  ktemp2[105] = elStorage->K55[1] + elStorage->S55[1] * c_c_tmp;
  ktemp2[129] = K56_idx_1;
  ktemp2[11] = K16_idx_2;
  S24_idx_0 = elStorage->K26[2] + S26_idx_2;
  ktemp2[35] = S24_idx_0 - H26_idx_2;
  H12_idx_0 = elStorage->K36[2] + S36_idx_2;
  ktemp2[59] = H12_idx_0 - H36_idx_2;
  K15_idx_0 = elStorage->K46[2] + S46_idx_2;
  ktemp2[83] = K15_idx_0 - H46_idx_2;
  ktemp2[107] = K56_idx_2;
  ktemp2[131] = elStorage->K66[1] + elStorage->S66[1] * c_c_tmp;
  ktemp2[12] = elStorage->K11[2] + elStorage->S11[2] * c_c_tmp;
  ktemp2[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  ktemp2[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  ktemp2[84] = ycm + H14_idx_2;
  ktemp2[108] = K15_idx_2;
  ktemp2[132] = K16_idx_2;
  ktemp2[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  ktemp2[38] = elStorage->K22[2] + elStorage->S22[2] * b_c_tmp;
  ktemp2[62] = f_ktemp2_tmp + H23_idx_2;
  ktemp2[86] = rhoA + H24_idx_2;
  ktemp2[110] = S13_idx_0 + H25_idx_2;
  ktemp2[134] = S24_idx_0 + H26_idx_2;
  ktemp2[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  ktemp2[40] = b_ktemp2_tmp - H23_idx_1;
  ktemp2[64] = elStorage->K33[2] + elStorage->S33[2] * c_tmp;
  ktemp2[88] = zcm + H34_idx_2;
  ktemp2[112] = S26_idx_0 + H35_idx_2;
  ktemp2[136] = H12_idx_0 + H36_idx_2;
  ktemp2[18] = ktemp2_tmp - H14_idx_1;
  ktemp2[42] = c_ktemp2_tmp - H24_idx_1;
  ktemp2[66] = g_ktemp2_tmp - H34_idx_1;
  ktemp2[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * b_c_tmp +
    elStorage->S44_2[2] * c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  ktemp2[114] = S25_idx_0 + H45_idx_2;
  ktemp2[138] = K15_idx_0 + H46_idx_2;
  ktemp2[20] = K15_idx_1;
  ktemp2[44] = d_ktemp2_tmp - H25_idx_1;
  ktemp2[68] = O1 - H35_idx_1;
  ktemp2[92] = S34_idx_0 - H45_idx_1;
  ktemp2[116] = elStorage->K55[2] + elStorage->S55[2] * c_c_tmp;
  ktemp2[140] = K56_idx_2;
  ktemp2[22] = K16_idx_1;
  ktemp2[46] = e_ktemp2_tmp - H26_idx_1;
  ktemp2[70] = integrationFactor - H36_idx_1;
  ktemp2[94] = S35_idx_0 - H46_idx_1;
  ktemp2[118] = K56_idx_1;
  ktemp2[142] = elStorage->K66[2] + elStorage->S66[2] * c_c_tmp;
  ktemp2[13] = elStorage->K11[3] + elStorage->S11[3] * c_c_tmp;
  ktemp2[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  ktemp2[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  ktemp2_tmp = elStorage->K14[3] + S14_idx_3;
  ktemp2[85] = ktemp2_tmp + H14_idx_3;
  ktemp2[109] = K15_idx_3;
  ktemp2[133] = K16_idx_3;
  ktemp2[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  ktemp2[39] = elStorage->K22[3] + elStorage->S22[3] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[3] + S23_idx_3;
  ktemp2[63] = b_ktemp2_tmp + H23_idx_3;
  c_ktemp2_tmp = elStorage->K24[3] + S24_idx_3;
  ktemp2[87] = c_ktemp2_tmp + H24_idx_3;
  d_ktemp2_tmp = elStorage->K25[3] + S25_idx_3;
  ktemp2[111] = d_ktemp2_tmp + H25_idx_3;
  e_ktemp2_tmp = elStorage->K26[3] + S26_idx_3;
  ktemp2[135] = e_ktemp2_tmp + H26_idx_3;
  ktemp2[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  ktemp2[41] = b_ktemp2_tmp - H23_idx_3;
  ktemp2[65] = elStorage->K33[3] + elStorage->S33[3] * c_tmp;
  b_ktemp2_tmp = elStorage->K34[3] + S34_idx_3;
  ktemp2[89] = b_ktemp2_tmp + H34_idx_3;
  f_ktemp2_tmp = elStorage->K35[3] + S35_idx_3;
  ktemp2[113] = f_ktemp2_tmp + H35_idx_3;
  g_ktemp2_tmp = elStorage->K36[3] + S36_idx_3;
  ktemp2[137] = g_ktemp2_tmp + H36_idx_3;
  ktemp2[19] = ktemp2_tmp - H14_idx_3;
  ktemp2[43] = c_ktemp2_tmp - H24_idx_3;
  ktemp2[67] = b_ktemp2_tmp - H34_idx_3;
  ktemp2[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * b_c_tmp +
    elStorage->S44_2[3] * c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[3] + S45_idx_3;
  ktemp2[115] = ktemp2_tmp + H45_idx_3;
  b_ktemp2_tmp = elStorage->K46[3] + S46_idx_3;
  ktemp2[139] = b_ktemp2_tmp + H46_idx_3;
  ktemp2[21] = K15_idx_3;
  ktemp2[45] = d_ktemp2_tmp - H25_idx_3;
  ktemp2[69] = f_ktemp2_tmp - H35_idx_3;
  ktemp2[93] = ktemp2_tmp - H45_idx_3;
  ktemp2[117] = elStorage->K55[3] + elStorage->S55[3] * c_c_tmp;
  ktemp2[141] = K56_idx_3;
  ktemp2[23] = K16_idx_3;
  ktemp2[47] = e_ktemp2_tmp - H26_idx_3;
  ktemp2[71] = g_ktemp2_tmp - H36_idx_3;
  ktemp2[95] = b_ktemp2_tmp - H46_idx_3;
  ktemp2[119] = K56_idx_3;
  ktemp2[143] = elStorage->K66[3] + elStorage->S66[3] * c_c_tmp;
  mapMatrixNonSym(ktemp2, Khate);

  //  Declare type
  // compile Coriolis/damping matrix
  //  ktemp = [zm,C12,C13,C14,zm,zm;
  //      -C12',zm,C23,C24,C25,C26;
  //      -C13',-C23',C33,C34,C35,C36;
  //      -C14',-C24',C43,C44,C45,C46;
  //      zm,-C25',-C35',-C45',zm,zm;
  //      zm,-C26',-C36',-C46',zm,zm];
  std::memset(&ktemp2[0], 0, 144U * sizeof(double));

  //  Row 1
  ktemp2_tmp = elStorage->C12[0] * Oel[2];
  ktemp2[2] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C12[2] * Oel[2];
  ktemp2[3] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C13[0] * Oel[1];
  ktemp2[4] = c_ktemp2_tmp;
  d_ktemp2_tmp = elStorage->C13[2] * Oel[1];
  ktemp2[5] = d_ktemp2_tmp;
  e_ktemp2_tmp = elStorage->C14_1[0] * Oel[1] + elStorage->C14_2[0] * Oel[2];
  ktemp2[6] = e_ktemp2_tmp;
  f_ktemp2_tmp = elStorage->C14_1[2] * Oel[1] + elStorage->C14_2[2] * Oel[2];
  ktemp2[7] = f_ktemp2_tmp;

  //  Row 2
  g_ktemp2_tmp = elStorage->C12[1] * Oel[2];
  ktemp2[14] = g_ktemp2_tmp;
  O1 = elStorage->C12[3] * Oel[2];
  ktemp2[15] = O1;
  integrationFactor = elStorage->C13[1] * Oel[1];
  ktemp2[16] = integrationFactor;
  ycm = elStorage->C13[3] * Oel[1];
  ktemp2[17] = ycm;
  rhoA = elStorage->C14_1[1] * Oel[1] + elStorage->C14_2[1] * Oel[2];
  ktemp2[18] = rhoA;
  zcm = elStorage->C14_1[3] * Oel[1] + elStorage->C14_2[3] * Oel[2];
  ktemp2[19] = zcm;

  //  Row 3
  ktemp2[24] = -ktemp2_tmp;
  ktemp2[25] = -g_ktemp2_tmp;
  ktemp2_tmp = elStorage->C23[0] * Oel[0];
  ktemp2[28] = ktemp2_tmp;
  g_ktemp2_tmp = elStorage->C23[2] * Oel[0];
  ktemp2[29] = g_ktemp2_tmp;
  S34_idx_0 = elStorage->C24[0] * Oel[0];
  ktemp2[30] = S34_idx_0;
  S35_idx_0 = elStorage->C24[2] * Oel[0];
  ktemp2[31] = S35_idx_0;
  S13_idx_0 = elStorage->C25[0] * Oel[2];
  ktemp2[32] = S13_idx_0;
  S26_idx_0 = elStorage->C25[2] * Oel[2];
  ktemp2[33] = S26_idx_0;
  S25_idx_0 = elStorage->C26[0] * Oel[2];
  ktemp2[34] = S25_idx_0;
  S24_idx_0 = elStorage->C26[2] * Oel[2];
  ktemp2[35] = S24_idx_0;

  //  Row 4
  ktemp2[36] = -b_ktemp2_tmp;
  ktemp2[37] = -O1;
  b_ktemp2_tmp = elStorage->C23[1] * Oel[0];
  ktemp2[40] = b_ktemp2_tmp;
  O1 = elStorage->C23[3] * Oel[0];
  ktemp2[41] = O1;
  H12_idx_0 = elStorage->C24[1] * Oel[0];
  ktemp2[42] = H12_idx_0;
  K15_idx_0 = elStorage->C24[3] * Oel[0];
  ktemp2[43] = K15_idx_0;
  H34_idx_3 = elStorage->C25[1] * Oel[2];
  ktemp2[44] = H34_idx_3;
  S14_idx_0 = elStorage->C25[3] * Oel[2];
  ktemp2[45] = S14_idx_0;
  S12_idx_0 = elStorage->C26[1] * Oel[2];
  ktemp2[46] = S12_idx_0;
  S23_idx_0 = elStorage->C26[3] * Oel[2];
  ktemp2[47] = S23_idx_0;

  //  Row 5
  ktemp2[48] = -c_ktemp2_tmp;
  ktemp2[49] = -integrationFactor;
  ktemp2[50] = -ktemp2_tmp;
  ktemp2[51] = -b_ktemp2_tmp;
  ktemp2[52] = 0.0;
  ktemp2[53] = 0.0;
  ktemp2[54] = C34_idx_0;
  ktemp2[55] = C34_idx_2;
  ktemp2_tmp = elStorage->C35[0] * Oel[1];
  ktemp2[56] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C35[2] * Oel[1];
  ktemp2[57] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C36[0] * Oel[1];
  ktemp2[58] = c_ktemp2_tmp;
  integrationFactor = elStorage->C36[2] * Oel[1];
  ktemp2[59] = integrationFactor;

  //  Row 6
  ktemp2[60] = -d_ktemp2_tmp;
  ktemp2[61] = -ycm;
  ktemp2[62] = -g_ktemp2_tmp;
  ktemp2[63] = -O1;
  ktemp2[64] = 0.0;
  ktemp2[65] = 0.0;
  ktemp2[66] = C34_idx_1;
  ktemp2[67] = C34_idx_3;
  d_ktemp2_tmp = elStorage->C35[1] * Oel[1];
  ktemp2[68] = d_ktemp2_tmp;
  g_ktemp2_tmp = elStorage->C35[3] * Oel[1];
  ktemp2[69] = g_ktemp2_tmp;
  O1 = elStorage->C36[1] * Oel[1];
  ktemp2[70] = O1;
  ycm = elStorage->C36[3] * Oel[1];
  ktemp2[71] = ycm;

  //  Row 7
  ktemp2[72] = -e_ktemp2_tmp;
  ktemp2[73] = -rhoA;
  ktemp2[74] = -S34_idx_0;
  ktemp2[75] = -H12_idx_0;
  ktemp2[76] = -C34_idx_0;
  ktemp2[77] = -C34_idx_1;
  ktemp2[78] = 0.0;
  ktemp2[79] = 0.0;
  e_ktemp2_tmp = elStorage->C45_1[0] * Oel[2] + elStorage->C45_2[0] * Oel[1];
  ktemp2[80] = e_ktemp2_tmp;
  rhoA = elStorage->C45_1[2] * Oel[2] + elStorage->C45_2[2] * Oel[1];
  ktemp2[81] = rhoA;
  S34_idx_0 = elStorage->C46_1[0] * Oel[1] + elStorage->C46_2[0] * Oel[2];
  ktemp2[82] = S34_idx_0;
  H12_idx_0 = elStorage->C46_1[2] * Oel[1] + elStorage->C46_2[2] * Oel[2];
  ktemp2[83] = H12_idx_0;

  //  Row 8
  ktemp2[84] = -f_ktemp2_tmp;
  ktemp2[85] = -zcm;
  ktemp2[86] = -S35_idx_0;
  ktemp2[87] = -K15_idx_0;
  ktemp2[88] = -C34_idx_2;
  ktemp2[89] = -C34_idx_3;
  ktemp2[90] = 0.0;
  ktemp2[91] = 0.0;
  f_ktemp2_tmp = elStorage->C45_1[1] * Oel[2] + elStorage->C45_2[1] * Oel[1];
  ktemp2[92] = f_ktemp2_tmp;
  zcm = elStorage->C45_1[3] * Oel[2] + elStorage->C45_2[3] * Oel[1];
  ktemp2[93] = zcm;
  S35_idx_0 = elStorage->C46_1[1] * Oel[1] + elStorage->C46_2[1] * Oel[2];
  ktemp2[94] = S35_idx_0;
  K15_idx_0 = elStorage->C46_1[3] * Oel[1] + elStorage->C46_2[3] * Oel[2];
  ktemp2[95] = K15_idx_0;

  //  Row 9
  ktemp2[98] = -S13_idx_0;
  ktemp2[99] = -H34_idx_3;
  ktemp2[100] = -ktemp2_tmp;
  ktemp2[101] = -d_ktemp2_tmp;
  ktemp2[102] = -e_ktemp2_tmp;
  ktemp2[103] = -f_ktemp2_tmp;

  //  Row 10
  ktemp2[110] = -S26_idx_0;
  ktemp2[111] = -S14_idx_0;
  ktemp2[112] = -b_ktemp2_tmp;
  ktemp2[113] = -g_ktemp2_tmp;
  ktemp2[114] = -rhoA;
  ktemp2[115] = -zcm;

  //  Row 11
  ktemp2[122] = -S25_idx_0;
  ktemp2[123] = -S12_idx_0;
  ktemp2[124] = -c_ktemp2_tmp;
  ktemp2[125] = -O1;
  ktemp2[126] = -S34_idx_0;
  ktemp2[127] = -S35_idx_0;

  //  Row 12
  ktemp2[134] = -S24_idx_0;
  ktemp2[135] = -S23_idx_0;
  ktemp2[136] = -integrationFactor;
  ktemp2[137] = -ycm;
  ktemp2[138] = -H12_idx_0;
  ktemp2[139] = -K15_idx_0;
  mapMatrixNonSym(ktemp2, Ce);

  // compile mass matrix
  ktemp2[0] = elStorage->M11[0];
  ktemp2[24] = 0.0;
  ktemp2[48] = 0.0;
  ktemp2[72] = 0.0;
  ktemp2[96] = elStorage->M15[0];
  ktemp2[120] = elStorage->M16[0];
  ktemp2[2] = 0.0;
  ktemp2[26] = elStorage->M22[0];
  ktemp2[50] = 0.0;
  ktemp2[74] = elStorage->M24[0];
  ktemp2[98] = 0.0;
  ktemp2[122] = 0.0;
  ktemp2[4] = 0.0;
  ktemp2[28] = 0.0;
  ktemp2[52] = elStorage->M33[0];
  ktemp2[76] = elStorage->M34[0];
  ktemp2[100] = 0.0;
  ktemp2[124] = 0.0;
  ktemp2[6] = 0.0;
  ktemp2[30] = elStorage->M24[0];
  ktemp2[54] = elStorage->M34[0];
  ktemp2[78] = elStorage->M44[0];
  ktemp2[102] = 0.0;
  ktemp2[126] = 0.0;
  ktemp2[8] = elStorage->M15[0];
  ktemp2[32] = 0.0;
  ktemp2[56] = 0.0;
  ktemp2[80] = 0.0;
  ktemp2[104] = elStorage->M55[0];
  ktemp2[128] = elStorage->M56[0];
  ktemp2[10] = elStorage->M16[0];
  ktemp2[34] = 0.0;
  ktemp2[58] = 0.0;
  ktemp2[82] = 0.0;
  ktemp2[106] = elStorage->M56[0];
  ktemp2[130] = elStorage->M66[0];
  ktemp2[1] = elStorage->M11[1];
  ktemp2[25] = 0.0;
  ktemp2[49] = 0.0;
  ktemp2[73] = 0.0;
  ktemp2[97] = elStorage->M15[1];
  ktemp2[121] = elStorage->M16[1];
  ktemp2[3] = 0.0;
  ktemp2[27] = elStorage->M22[1];
  ktemp2[51] = 0.0;
  ktemp2[75] = elStorage->M24[1];
  ktemp2[99] = 0.0;
  ktemp2[123] = 0.0;
  ktemp2[5] = 0.0;
  ktemp2[29] = 0.0;
  ktemp2[53] = elStorage->M33[1];
  ktemp2[77] = elStorage->M34[1];
  ktemp2[101] = 0.0;
  ktemp2[125] = 0.0;
  ktemp2[7] = 0.0;
  ktemp2[31] = elStorage->M24[2];
  ktemp2[55] = elStorage->M34[2];
  ktemp2[79] = elStorage->M44[1];
  ktemp2[103] = 0.0;
  ktemp2[127] = 0.0;
  ktemp2[9] = elStorage->M15[2];
  ktemp2[33] = 0.0;
  ktemp2[57] = 0.0;
  ktemp2[81] = 0.0;
  ktemp2[105] = elStorage->M55[1];
  ktemp2[129] = elStorage->M56[1];
  ktemp2[11] = elStorage->M16[2];
  ktemp2[35] = 0.0;
  ktemp2[59] = 0.0;
  ktemp2[83] = 0.0;
  ktemp2[107] = elStorage->M56[2];
  ktemp2[131] = elStorage->M66[1];
  ktemp2[12] = elStorage->M11[2];
  ktemp2[36] = 0.0;
  ktemp2[60] = 0.0;
  ktemp2[84] = 0.0;
  ktemp2[108] = elStorage->M15[2];
  ktemp2[132] = elStorage->M16[2];
  ktemp2[14] = 0.0;
  ktemp2[38] = elStorage->M22[2];
  ktemp2[62] = 0.0;
  ktemp2[86] = elStorage->M24[2];
  ktemp2[110] = 0.0;
  ktemp2[134] = 0.0;
  ktemp2[16] = 0.0;
  ktemp2[40] = 0.0;
  ktemp2[64] = elStorage->M33[2];
  ktemp2[88] = elStorage->M34[2];
  ktemp2[112] = 0.0;
  ktemp2[136] = 0.0;
  ktemp2[18] = 0.0;
  ktemp2[42] = elStorage->M24[1];
  ktemp2[66] = elStorage->M34[1];
  ktemp2[90] = elStorage->M44[2];
  ktemp2[114] = 0.0;
  ktemp2[138] = 0.0;
  ktemp2[20] = elStorage->M15[1];
  ktemp2[44] = 0.0;
  ktemp2[68] = 0.0;
  ktemp2[92] = 0.0;
  ktemp2[116] = elStorage->M55[2];
  ktemp2[140] = elStorage->M56[2];
  ktemp2[22] = elStorage->M16[1];
  ktemp2[46] = 0.0;
  ktemp2[70] = 0.0;
  ktemp2[94] = 0.0;
  ktemp2[118] = elStorage->M56[1];
  ktemp2[142] = elStorage->M66[2];
  ktemp2[13] = elStorage->M11[3];
  ktemp2[37] = 0.0;
  ktemp2[61] = 0.0;
  ktemp2[85] = 0.0;
  ktemp2[109] = elStorage->M15[3];
  ktemp2[133] = elStorage->M16[3];
  ktemp2[15] = 0.0;
  ktemp2[39] = elStorage->M22[3];
  ktemp2[63] = 0.0;
  ktemp2[87] = elStorage->M24[3];
  ktemp2[111] = 0.0;
  ktemp2[135] = 0.0;
  ktemp2[17] = 0.0;
  ktemp2[41] = 0.0;
  ktemp2[65] = elStorage->M33[3];
  ktemp2[89] = elStorage->M34[3];
  ktemp2[113] = 0.0;
  ktemp2[137] = 0.0;
  ktemp2[19] = 0.0;
  ktemp2[43] = elStorage->M24[3];
  ktemp2[67] = elStorage->M34[3];
  ktemp2[91] = elStorage->M44[3];
  ktemp2[115] = 0.0;
  ktemp2[139] = 0.0;
  ktemp2[21] = elStorage->M15[3];
  ktemp2[45] = 0.0;
  ktemp2[69] = 0.0;
  ktemp2[93] = 0.0;
  ktemp2[117] = elStorage->M55[3];
  ktemp2[141] = elStorage->M56[3];
  ktemp2[23] = elStorage->M16[3];
  ktemp2[47] = 0.0;
  ktemp2[71] = 0.0;
  ktemp2[95] = 0.0;
  ktemp2[119] = elStorage->M56[3];
  ktemp2[143] = elStorage->M66[3];
  mapMatrixNonSym(ktemp2, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input->RayleighAlpha * Kenr[i] + input->RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      ktemp2[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(ktemp2, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Khate,
                  lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Khate);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fhate[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fhate[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fhate, Fe);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input->concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input->concMass[0];
    Me[13] += input->concMass[0];
    Me[26] += input->concMass[0];
    Me[39] += input->concMass[1];
    Me[52] += input->concMass[2];
    Me[65] += input->concMass[3];
    Me[78] += input->concMass[4];
    Me[91] += input->concMass[4];
    Me[104] += input->concMass[4];
    Me[117] += input->concMass[5];
    Me[130] += input->concMass[6];
    Me[143] += input->concMass[7];

    // modify Ce for concentrated mass
    H34_idx_3 = 2.0 * input->concMass[0] * Omega;
    Ce[12] -= H34_idx_3;
    Ce[1] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[0] * 0.0;
    Ce[24] += H34_idx_3;
    Ce[2] -= H34_idx_3;
    Ce[25] -= H34_idx_3;
    Ce[14] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[4] * Omega;
    Ce[90] -= H34_idx_3;
    Ce[79] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[4] * 0.0;
    Ce[102] += H34_idx_3;
    Ce[80] -= H34_idx_3;
    Ce[103] -= H34_idx_3;
    Ce[92] += H34_idx_3;

    // modify Ke for concentrated mass
    H34_idx_3 = Omega * Omega;
    S23_idx_0 = input->concMass[0] * H34_idx_3;
    Khate[0] -= S23_idx_0;
    S14_idx_0 = input->concMass[0] * 0.0 * 0.0;
    S12_idx_0 = input->concMass[0] * OmegaDot;
    Khate[12] = (Khate[12] + S14_idx_0) - S12_idx_0;
    Khate[1] = (Khate[1] + S14_idx_0) + S12_idx_0;
    S14_idx_0 = input->concMass[0] * 0.0 * Omega;
    Khate[24] = (Khate[24] + S14_idx_0) + input->concMass[0] * 0.0;
    Khate[2] = (Khate[2] + S14_idx_0) - input->concMass[0] * 0.0;
    Khate[25] = (Khate[25] + S14_idx_0) - input->concMass[0] * 0.0;
    Khate[14] = (Khate[14] + S14_idx_0) + input->concMass[0] * 0.0;
    Khate[13] -= S23_idx_0;
    Khate[26] -= input->concMass[0] * 0.0;
    S23_idx_0 = input->concMass[4] * H34_idx_3;
    Khate[78] -= S23_idx_0;
    S14_idx_0 = input->concMass[4] * 0.0 * 0.0;
    S12_idx_0 = input->concMass[4] * OmegaDot;
    Khate[90] = (Khate[90] + S14_idx_0) - S12_idx_0;
    Khate[79] = (Khate[79] + S14_idx_0) + S12_idx_0;
    S14_idx_0 = input->concMass[4] * 0.0 * Omega;
    Khate[102] = (Khate[102] + S14_idx_0) + input->concMass[4] * 0.0;
    Khate[80] = (Khate[80] + S14_idx_0) - input->concMass[4] * 0.0;
    Khate[103] = (Khate[103] + S14_idx_0) - input->concMass[4] * 0.0;
    Khate[92] = (Khate[92] + S14_idx_0) + input->concMass[4] * 0.0;
    Khate[91] -= S23_idx_0;
    Khate[104] -= input->concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    H34_idx_3 = Omega * Omega;
    S23_idx_0 = 0.0 * Omega * input->z.data[0];
    Fe[0] = ((Fe[0] + input->concMass[0] * ((input->x.data[0] * H34_idx_3 - 0.0 *
                input->y.data[0]) - S23_idx_0)) + input->concMass[0] *
             (input->y.data[0] * OmegaDot - input->z.data[0] * 0.0)) -
      input->concMass[0] * a_temp[0];
    Fe[1] = ((Fe[1] + input->concMass[0] * ((input->y.data[0] * H34_idx_3 -
                S23_idx_0) - 0.0 * input->x.data[0])) + input->concMass[0] *
             (input->z.data[0] * 0.0 - input->x.data[0] * OmegaDot)) -
      input->concMass[0] * a_temp[1];
    Fe[2] = ((Fe[2] + input->concMass[0] * ((input->z.data[0] * 0.0 - Omega *
                0.0 * input->x.data[0]) - Omega * 0.0 * input->y.data[0])) +
             input->concMass[0] * (input->x.data[0] * 0.0 - input->y.data[0] *
              0.0)) - input->concMass[0] * a_temp[2];
    S23_idx_0 = 0.0 * Omega * input->z.data[1];
    Fe[6] = ((Fe[6] + input->concMass[4] * ((input->x.data[1] * H34_idx_3 - 0.0 *
                input->y.data[1]) - S23_idx_0)) + input->concMass[4] *
             (input->y.data[1] * OmegaDot - input->z.data[1] * 0.0)) -
      input->concMass[4] * a_temp[0];
    Fe[7] = ((Fe[7] + input->concMass[4] * ((input->y.data[1] * H34_idx_3 -
                S23_idx_0) - 0.0 * input->x.data[1])) + input->concMass[4] *
             (input->z.data[1] * 0.0 - input->x.data[1] * OmegaDot)) -
      input->concMass[4] * a_temp[1];
    Fe[8] = ((Fe[8] + input->concMass[4] * ((input->z.data[1] * 0.0 - Omega *
                0.0 * input->x.data[1]) - Omega * 0.0 * input->y.data[1])) +
             input->concMass[4] * (input->x.data[1] * 0.0 - input->y.data[1] *
              0.0)) - input->concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  std::memset(&Fhate[0], 0, 12U * sizeof(double));
  std::memset(&FhatLessConc[0], 0, 12U * sizeof(double));
  if (g_strcmp(input->analysisType.data, input->analysisType.size)) {
    // calculate effective stiffness matrix and force vector for Dean integrator 
    for (i = 0; i < 12; i++) {
      FhatLessConc[i] = dispm1_data[i];
      b_Khate[i] = 2.0 * input->disp.data[i] - dispm1_data[i];
      Fhate[i] = -0.001 * dispm1_data[i] - 0.001 * input->disp.data[i];
    }

    for (i = 0; i < 12; i++) {
      K15_idx_0 = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        K15_idx_0 += Me[i + 12 * b_i] * b_Khate[b_i];
      }

      Ftemp_data[i] = Fe[i] * 2000.0 + K15_idx_0;
    }

    for (i = 0; i < 12; i++) {
      K15_idx_0 = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        K15_idx_0 += Khate[i + 12 * b_i] * Fhate[b_i];
      }

      b_Khate[i] = K15_idx_0;
    }

    for (i = 0; i < 12; i++) {
      K15_idx_0 = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        K15_idx_0 += Ce[i + 12 * b_i] * (1.0E+6 * FhatLessConc[b_i]);
      }

      Fhate[i] = (Ftemp_data[i] + b_Khate[i]) + K15_idx_0;
    }

    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));

    // ........................................................
    // ..........................................................
    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] * 0.001 + 1.0E+6 * Ce[i]) + Me[i];
    }

    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  if (h_strcmp(input->analysisType.data, input->analysisType.size)) {
    // calculate effective stiffness matrix and load vector for Newmark-Beta integrator 
    //      a1 = timeInt.a1;
    //      a2 = timeInt.a2;
    if (e_strcmp(input->iterationType)) {
      // considerations if newton raphson iteration is used
      if (0 <= dispddot_size_idx_1 - 1) {
        std::memcpy(&b_data[0], &dispddot_data[0], dispddot_size_idx_1 * sizeof
                    (double));
      }

      if (dispddot_size_idx_1 == 1) {
        for (i = 0; i < 12; i++) {
          K15_idx_0 = 0.0;
          for (b_i = 0; b_i < 12; b_i++) {
            K15_idx_0 += Me[i + 12 * b_i] * b_data[b_i];
          }

          Fhate[i] = K15_idx_0;
        }
      } else {
        mtimes(Me, b_data, Fhate);
      }

      if (0 <= dispdot_size_idx_1 - 1) {
        std::memcpy(&b_data[0], &dispdot_data[0], dispdot_size_idx_1 * sizeof
                    (double));
      }

      if (dispdot_size_idx_1 == 1) {
        for (i = 0; i < 12; i++) {
          K15_idx_0 = 0.0;
          for (b_i = 0; b_i < 12; b_i++) {
            K15_idx_0 += Ce[i + 12 * b_i] * b_data[b_i];
          }

          b_Khate[i] = K15_idx_0;
        }
      } else {
        mtimes(Ce, b_data, b_Khate);
      }

      std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
      mtimes(Khate, b_data, FhatLessConc);
      for (i = 0; i < 12; i++) {
        Fhate[i] = ((Fe[i] - Fhate[i]) - b_Khate[i]) - FhatLessConc[i];
      }
    } else {
      if (f_strcmp(input->iterationType)) {
        // considerations if direct iteration is used or linear analysis
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, b_Khate);
        for (i = 0; i < 12; i++) {
          Fhate[i] = (Fe[i] + Fhate[i]) + b_Khate[i];
        }
      }
    }

    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] + 1.0E+6 * Me[i]) + 1000.0 * Ce[i];
    }

    // ........................................................
    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));
    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  if (i_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&FhatLessConc[0], &Fe[0], 12U * sizeof(double));
  }

  if ((i_strcmp(input->analysisType.data, input->analysisType.size) || k_strcmp
       (input->analysisType.data, input->analysisType.size)) && e_strcmp
      (input->iterationType)) {
    // considerations for newton-raphson iteration
    std::memcpy(&b_data[0], &disp_iter_data[0], 12U * sizeof(double));
    mtimes(Khate, b_data, b_Khate);
    for (i = 0; i < 12; i++) {
      Fe[i] -= b_Khate[i];
    }
  }

  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Khate[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &Fe[0], 12U * sizeof(double));
  output->Me.size[0] = 1;
  output->Me.size[1] = 1;
  output->Me.data[0] = 0.0;
  output->Ce.size[0] = 1;
  output->Ce.size[1] = 1;
  output->Ce.data[0] = 0.0;
  if (i_strcmp(input->analysisType.data, input->analysisType.size) || j_strcmp
      (input->analysisType.data, input->analysisType.size)) {
    output->Me.size[0] = 12;
    output->Me.size[1] = 12;
    output->Ce.size[0] = 12;
    output->Ce.size[1] = 12;
    std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
    std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));
  }

  if (g_strcmp(input->analysisType.data, input->analysisType.size) || h_strcmp
      (input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&output->FhatLessConc[0], &FhatLessConc[0], 12U * sizeof(double));
  }

  // ------------------------------------------
}

//
// calculateTimoshenkoElementNL performs nonlinear element calculations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [output] = calculateTimoshenkoElementNL(input,elStorage)
//
//    This function performs nonlinear element calculations.
//
//       input:
//       input      = object containing element input
//       elStorage  = obect containing precalculated element data
//
//       output:
//       output     = object containing element data
// Arguments    : const double input_xloc[2]
//                const double input_sectionProps_twist[2]
//                const double input_sectionProps_rhoA[2]
//                const double input_sectionProps_EA[2]
//                const double input_sectionProps_zcm[2]
//                const double input_sectionProps_ycm[2]
//                double input_sweepAngle
//                double input_coneAngle
//                double input_rollAngle
//                const double input_concMass[8]
//                const double input_disp_data[]
//                const double input_x_data[]
//                const double input_y_data[]
//                const double input_z_data[]
//                const double input_CN2H[9]
//                double input_RayleighAlpha
//                double input_RayleighBeta
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void c_calculateTimoshenkoElementNL(const double input_xloc[2], const double
  input_sectionProps_twist[2], const double input_sectionProps_rhoA[2], const
  double input_sectionProps_EA[2], const double input_sectionProps_zcm[2], const
  double input_sectionProps_ycm[2], double input_sweepAngle, double
  input_coneAngle, double input_rollAngle, const double input_concMass[8], const
  double input_disp_data[], const double input_x_data[], const double
  input_y_data[], const double input_z_data[], const double input_CN2H[9],
  double input_RayleighAlpha, double input_RayleighBeta, const f_struct_T
  *elStorage, o_struct_T *output)
{
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double SS22_data_idx_0;
  double SS33_data_idx_0;
  double SS22_data_idx_1;
  double SS33_data_idx_1;
  double SS22_data_idx_2;
  double SS33_data_idx_2;
  double SS22_data_idx_3;
  double SS33_data_idx_3;
  double lambda[144];
  double b_data[12];
  double dispLocal[12];
  int i;
  double O1;
  double Oel[3];
  double K16_idx_0;
  double O2;
  double K15_idx_0;
  double O3;
  double a_temp[3];
  double O1dot;
  double ODotel[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double H14_idx_0;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double integrationFactor;
  double c_tmp;
  double b_c_tmp;
  double S12_idx_0;
  double S13_idx_0;
  double rhoA;
  double S23_idx_0;
  double S25_idx_0;
  double S14_idx_0;
  double S26_idx_0;
  double Faxial;
  double S35_idx_0;
  double S24_idx_0;
  double S34_idx_0;
  double C34_idx_0;
  double H12_idx_0;
  double H13_idx_0;
  double zcm;
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double ycm;
  double S25_idx_1;
  double S26_idx_1;
  double S35_idx_1;
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double disMomentgp[3];
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C34_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double posLocal[3];
  double H23_idx_1;
  double disLoadgpLocal[3];
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C34_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C34_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H34_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double H46_idx_3;
  double ktemp2[144];
  double Kenr[144];
  double c_c_tmp;
  double K15_idx_1;
  double K16_idx_1;
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double ktemp2_tmp;
  double b_ktemp2_tmp;
  double c_ktemp2_tmp;
  double d_ktemp2_tmp;
  double e_ktemp2_tmp;
  double Ke[144];
  double Ce[144];
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fel_data[12];
  int tmp_size[1];
  boolean_T concMassFlag;

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  // declare type
  // declare type
  // declare type
  // options for Dean integrator
  // --------------------------------------------
  // setting for modal analysis flag
  // setting for initial reduced order model calculations
  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  SS22_data_idx_0 = 0.0;
  SS33_data_idx_0 = 0.0;
  SS22_data_idx_1 = 0.0;
  SS33_data_idx_1 = 0.0;
  SS22_data_idx_2 = 0.0;
  SS33_data_idx_2 = 0.0;
  SS22_data_idx_3 = 0.0;
  SS33_data_idx_3 = 0.0;

  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input_sweepAngle * 3.1415926535897931 / 180.0, input_coneAngle
                  * 3.1415926535897931 / 180.0, (input_rollAngle + 0.5 *
    (input_sectionProps_twist[0] + input_sectionProps_twist[1])) *
                  3.1415926535897931 / 180.0, lambda);
  std::memcpy(&b_data[0], &input_disp_data[0], 12U * sizeof(double));
  mtimes(lambda, b_data, dispLocal);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    K16_idx_0 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    K15_idx_0 = lambda[i + 24];
    a_temp[i] = (input_CN2H[i] * 0.0 + input_CN2H[i + 3] * 0.0) + input_CN2H[i +
      6] * 9.81;
    ODotel[i] = K16_idx_0 + K15_idx_0 * 0.0;
    Oel[i] = K16_idx_0 + K15_idx_0 * 3.2898681336964524;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  H14_idx_0 = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input_xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &K15_idx_0);
    integrationFactor = K15_idx_0 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input_sectionProps_rhoA[0] + N_data[1] *
      input_sectionProps_rhoA[1];
    K15_idx_0 = p_N_x_data[0] * dispLocal[1] + p_N_x_data[1] * dispLocal[7];
    S14_idx_0 = p_N_x_data[0] * dispLocal[2] + p_N_x_data[1] * dispLocal[8];
    Faxial = (N_data[0] * input_sectionProps_EA[0] + N_data[1] *
              input_sectionProps_EA[1]) * (((p_N_x_data[0] * dispLocal[0] +
      p_N_x_data[1] * dispLocal[6]) + 0.5 * (K15_idx_0 * K15_idx_0)) + 0.5 *
      (S14_idx_0 * S14_idx_0));

    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    ycm = N_data[0] * input_sectionProps_ycm[0] + N_data[1] *
      input_sectionProps_ycm[1];
    zcm = N_data[0] * input_sectionProps_zcm[0] + N_data[1] *
      input_sectionProps_zcm[1];
    S14_idx_0 = N_data[0] * input_x_data[0] + N_data[1] * input_x_data[1];
    S12_idx_0 = N_data[0] * input_y_data[0] + N_data[1] * input_y_data[1];
    S23_idx_0 = N_data[0] * input_z_data[0] + N_data[1] * input_z_data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      K16_idx_0 = lambda[i + 12];
      K15_idx_0 = lambda[i + 24];
      posLocal[i] = (lambda[i] * S14_idx_0 + K16_idx_0 * S12_idx_0) + K15_idx_0 *
        S23_idx_0;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + K16_idx_0 * disMomentgp
                           [1]) + K15_idx_0 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // stress-stiffening/pre-stress calculations
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // calculate static aerodynamic load
    // distributed/body force load calculations
    H12_idx_0 = (O2 * O2 + O3 * O3) * posLocal[0];
    S24_idx_0 = O2dot * posLocal[2];
    S25_idx_0 = O3dot * posLocal[1];
    S26_idx_0 = H14_idx_0 * posLocal[2];
    S13_idx_0 = O1 * O2 * posLocal[1];
    S35_idx_0 = rhoA * ((((H12_idx_0 - S13_idx_0) - S26_idx_0) + S25_idx_0) -
                        S24_idx_0) - disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = O1dot * posLocal[2];
    S14_idx_0 = O3dot * posLocal[0];
    S34_idx_0 = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 *
      O3) - posLocal[0] * O1 * O2) + K15_idx_0) - S14_idx_0) - disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    S12_idx_0 = O2dot * posLocal[0];
    S23_idx_0 = O1dot * posLocal[1];
    H13_idx_0 = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - H14_idx_0 *
      posLocal[0]) - O2 * O3 * posLocal[1]) + S12_idx_0) - S23_idx_0) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    K16_idx_0 = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) -
      posLocal[1] * (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                          (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                         (S12_idx_0 - S23_idx_0)) - zcm * (K15_idx_0 - S14_idx_0))
      - disMomentgp[0];

    // This function is a general routine to calculate an element vector
    S23_idx_0 = rhoA * zcm * ((((H12_idx_0 - posLocal[1] * O1 * O2) - posLocal[2]
      * O1 * O3) - S24_idx_0) + S25_idx_0) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = rhoA * ycm * ((((S26_idx_0 + S13_idx_0) - H12_idx_0) - S25_idx_0)
      + S24_idx_0) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    S14_idx_0 = Faxial * p_N_x_data[0];
    S12_idx_0 = S14_idx_0 * p_N_x_data[0] * integrationFactor;
    SS22_data_idx_0 += S12_idx_0;
    SS33_data_idx_0 += S12_idx_0;
    S12_idx_0 = S14_idx_0 * p_N_x_data[1] * integrationFactor;
    SS22_data_idx_2 += S12_idx_0;
    SS33_data_idx_2 += S12_idx_0;
    F1_data_idx_0 += S35_idx_0 * N_data[0] * integrationFactor;
    F2_data_idx_0 += S34_idx_0 * N_data[0] * integrationFactor;
    F3_data_idx_0 += H13_idx_0 * N_data[0] * integrationFactor;
    F4_data_idx_0 += K16_idx_0 * N_data[0] * integrationFactor;
    F5_data_idx_0 += S23_idx_0 * N_data[0] * integrationFactor;
    F6_data_idx_0 += K15_idx_0 * N_data[0] * integrationFactor;
    S14_idx_0 = Faxial * p_N_x_data[1];
    S12_idx_0 = S14_idx_0 * p_N_x_data[0] * integrationFactor;
    SS22_data_idx_1 += S12_idx_0;
    SS33_data_idx_1 += S12_idx_0;
    S12_idx_0 = S14_idx_0 * p_N_x_data[1] * integrationFactor;
    SS22_data_idx_3 += S12_idx_0;
    SS33_data_idx_3 += S12_idx_0;
    F1_data_idx_1 += S35_idx_0 * N_data[1] * integrationFactor;
    F2_data_idx_1 += S34_idx_0 * N_data[1] * integrationFactor;
    F3_data_idx_1 += H13_idx_0 * N_data[1] * integrationFactor;
    F4_data_idx_1 += K16_idx_0 * N_data[1] * integrationFactor;
    F5_data_idx_1 += S23_idx_0 * N_data[1] * integrationFactor;
    F6_data_idx_1 += K15_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  K15_idx_0 = Oel[0] * Oel[1];
  c_tmp = Oel[0] * Oel[0];
  b_c_tmp = c_tmp + Oel[2] * Oel[2];
  c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  S12_idx_0 = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  S23_idx_0 = elStorage->S23[0] * Oel[1] * Oel[2];
  S25_idx_0 = elStorage->S25[0] * K15_idx_0;
  S26_idx_0 = elStorage->S26[0] * K15_idx_0;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  rhoA = elStorage->S36[0] * Oel[0] * Oel[2];
  S14_idx_0 = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel
    [0] * Oel[1];
  S24_idx_0 = elStorage->S24_1[0] * b_c_tmp + elStorage->S24_2[0] * Oel[1] *
    Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * c_tmp + elStorage->S34_2[0] * Oel[1] * Oel[2];
  Faxial = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel[0] *
    Oel[1];
  integrationFactor = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0]
    * Oel[0] * Oel[2];
  K16_idx_0 = elStorage->C34[0];
  C34_idx_0 = K16_idx_0 * Oel[0];
  H12_idx_0 = 0.5 * elStorage->C12[0] * ODotel[2];
  H13_idx_0 = 0.5 * elStorage->C13[0] * ODotel[1];
  zcm = 0.5 * elStorage->C23[0] * ODotel[0];
  O1 = 0.5 * elStorage->C24[0] * ODotel[0];
  O2 = 0.5 * elStorage->C25[0] * ODotel[2];
  O3 = 0.5 * elStorage->C26[0] * ODotel[2];
  O1dot = 0.5 * K16_idx_0 * ODotel[0];
  O2dot = 0.5 * elStorage->C35[0] * ODotel[1];
  O3dot = 0.5 * elStorage->C36[0] * ODotel[1];
  H14_idx_0 = 0.5 * (elStorage->C14_1[0] * ODotel[1] + elStorage->C14_2[0] *
                     ODotel[2]);
  H45_idx_0 = 0.5 * (elStorage->C45_1[0] * ODotel[2] + elStorage->C45_2[0] *
                     ODotel[1]);
  H46_idx_0 = 0.5 * (elStorage->C46_1[0] * ODotel[1] + elStorage->C46_2[0] *
                     ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * K15_idx_0;
  S26_idx_1 = elStorage->S26[1] * K15_idx_0;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * b_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * c_tmp + elStorage->S34_2[1] * Oel[1] * Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  K16_idx_0 = elStorage->C34[1];
  C34_idx_1 = K16_idx_0 * Oel[0];
  H12_idx_1 = 0.5 * elStorage->C12[1] * ODotel[2];
  H13_idx_1 = 0.5 * elStorage->C13[1] * ODotel[1];
  H23_idx_1 = 0.5 * elStorage->C23[1] * ODotel[0];
  H24_idx_1 = 0.5 * elStorage->C24[1] * ODotel[0];
  H25_idx_1 = 0.5 * elStorage->C25[1] * ODotel[2];
  H26_idx_1 = 0.5 * elStorage->C26[1] * ODotel[2];
  H34_idx_1 = 0.5 * K16_idx_0 * ODotel[0];
  H35_idx_1 = 0.5 * elStorage->C35[1] * ODotel[1];
  H36_idx_1 = 0.5 * elStorage->C36[1] * ODotel[1];
  H14_idx_1 = 0.5 * (elStorage->C14_1[1] * ODotel[1] + elStorage->C14_2[1] *
                     ODotel[2]);
  H45_idx_1 = 0.5 * (elStorage->C45_1[1] * ODotel[2] + elStorage->C45_2[1] *
                     ODotel[1]);
  H46_idx_1 = 0.5 * (elStorage->C46_1[1] * ODotel[1] + elStorage->C46_2[1] *
                     ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * K15_idx_0;
  S26_idx_2 = elStorage->S26[2] * K15_idx_0;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * b_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * c_tmp + elStorage->S34_2[2] * Oel[1] * Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  K16_idx_0 = elStorage->C34[2];
  C34_idx_2 = K16_idx_0 * Oel[0];
  H12_idx_2 = 0.5 * elStorage->C12[2] * ODotel[2];
  H13_idx_2 = 0.5 * elStorage->C13[2] * ODotel[1];
  H23_idx_2 = 0.5 * elStorage->C23[2] * ODotel[0];
  H24_idx_2 = 0.5 * elStorage->C24[2] * ODotel[0];
  H25_idx_2 = 0.5 * elStorage->C25[2] * ODotel[2];
  H26_idx_2 = 0.5 * elStorage->C26[2] * ODotel[2];
  H34_idx_2 = 0.5 * K16_idx_0 * ODotel[0];
  H35_idx_2 = 0.5 * elStorage->C35[2] * ODotel[1];
  H36_idx_2 = 0.5 * elStorage->C36[2] * ODotel[1];
  H14_idx_2 = 0.5 * (elStorage->C14_1[2] * ODotel[1] + elStorage->C14_2[2] *
                     ODotel[2]);
  H45_idx_2 = 0.5 * (elStorage->C45_1[2] * ODotel[2] + elStorage->C45_2[2] *
                     ODotel[1]);
  H46_idx_2 = 0.5 * (elStorage->C46_1[2] * ODotel[1] + elStorage->C46_2[2] *
                     ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * K15_idx_0;
  S26_idx_3 = elStorage->S26[3] * K15_idx_0;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * b_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * c_tmp + elStorage->S34_2[3] * Oel[1] * Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  K16_idx_0 = elStorage->C34[3];
  C34_idx_3 = K16_idx_0 * Oel[0];
  H12_idx_3 = 0.5 * elStorage->C12[3] * ODotel[2];
  H13_idx_3 = 0.5 * elStorage->C13[3] * ODotel[1];
  H23_idx_3 = 0.5 * elStorage->C23[3] * ODotel[0];
  H24_idx_3 = 0.5 * elStorage->C24[3] * ODotel[0];
  H25_idx_3 = 0.5 * elStorage->C25[3] * ODotel[2];
  H26_idx_3 = 0.5 * elStorage->C26[3] * ODotel[2];
  H34_idx_3 = 0.5 * K16_idx_0 * ODotel[0];
  H35_idx_3 = 0.5 * elStorage->C35[3] * ODotel[1];
  H36_idx_3 = 0.5 * elStorage->C36[3] * ODotel[1];
  H14_idx_3 = 0.5 * (elStorage->C14_1[3] * ODotel[1] + elStorage->C14_2[3] *
                     ODotel[2]);
  H45_idx_3 = 0.5 * (elStorage->C45_1[3] * ODotel[2] + elStorage->C45_2[3] *
                     ODotel[1]);
  H46_idx_3 = 0.5 * (elStorage->C46_1[3] * ODotel[1] + elStorage->C46_2[3] *
                     ODotel[2]);

  // compile stiffness matrix without rotational effects
  ktemp2[0] = elStorage->K11[0];
  ktemp2[24] = elStorage->K12[0];
  ktemp2[48] = elStorage->K13[0];
  ktemp2[72] = elStorage->K14[0];
  ktemp2[96] = elStorage->K15[0];
  ktemp2[120] = elStorage->K16[0];
  ktemp2[2] = elStorage->K12[0];
  ktemp2[26] = elStorage->K22[0];
  ktemp2[50] = elStorage->K23[0];
  ktemp2[74] = elStorage->K24[0];
  ktemp2[98] = elStorage->K25[0];
  ktemp2[122] = elStorage->K26[0];
  ktemp2[4] = elStorage->K13[0];
  ktemp2[28] = elStorage->K23[0];
  ktemp2[52] = elStorage->K33[0];
  ktemp2[76] = elStorage->K34[0];
  ktemp2[100] = elStorage->K35[0];
  ktemp2[124] = elStorage->K36[0];
  ktemp2[6] = elStorage->K13[0];
  ktemp2[30] = elStorage->K24[0];
  ktemp2[54] = elStorage->K34[0];
  ktemp2[78] = elStorage->K44[0];
  ktemp2[102] = elStorage->K45[0];
  ktemp2[126] = elStorage->K46[0];
  ktemp2[8] = elStorage->K15[0];
  ktemp2[32] = elStorage->K25[0];
  ktemp2[56] = elStorage->K35[0];
  ktemp2[80] = elStorage->K45[0];
  ktemp2[104] = elStorage->K55[0];
  ktemp2[128] = elStorage->K56[0];
  ktemp2[10] = elStorage->K16[0];
  ktemp2[34] = elStorage->K26[0];
  ktemp2[58] = elStorage->K36[0];
  ktemp2[82] = elStorage->K46[0];
  ktemp2[106] = elStorage->K56[0];
  ktemp2[130] = elStorage->K66[0];
  ktemp2[1] = elStorage->K11[1];
  ktemp2[25] = elStorage->K12[1];
  ktemp2[49] = elStorage->K13[1];
  ktemp2[73] = elStorage->K14[1];
  ktemp2[97] = elStorage->K15[1];
  ktemp2[121] = elStorage->K16[1];
  ktemp2[3] = elStorage->K12[2];
  ktemp2[27] = elStorage->K22[1];
  ktemp2[51] = elStorage->K23[1];
  ktemp2[75] = elStorage->K24[1];
  ktemp2[99] = elStorage->K25[1];
  ktemp2[123] = elStorage->K26[1];
  ktemp2[5] = elStorage->K13[2];
  ktemp2[29] = elStorage->K23[2];
  ktemp2[53] = elStorage->K33[1];
  ktemp2[77] = elStorage->K34[1];
  ktemp2[101] = elStorage->K35[1];
  ktemp2[125] = elStorage->K36[1];
  ktemp2[7] = elStorage->K13[2];
  ktemp2[31] = elStorage->K24[2];
  ktemp2[55] = elStorage->K34[2];
  ktemp2[79] = elStorage->K44[1];
  ktemp2[103] = elStorage->K45[1];
  ktemp2[127] = elStorage->K46[1];
  ktemp2[9] = elStorage->K15[2];
  ktemp2[33] = elStorage->K25[2];
  ktemp2[57] = elStorage->K35[2];
  ktemp2[81] = elStorage->K45[2];
  ktemp2[105] = elStorage->K55[1];
  ktemp2[129] = elStorage->K56[1];
  ktemp2[11] = elStorage->K16[2];
  ktemp2[35] = elStorage->K26[2];
  ktemp2[59] = elStorage->K36[2];
  ktemp2[83] = elStorage->K46[2];
  ktemp2[107] = elStorage->K56[2];
  ktemp2[131] = elStorage->K66[1];
  ktemp2[12] = elStorage->K11[2];
  ktemp2[36] = elStorage->K12[2];
  ktemp2[60] = elStorage->K13[2];
  ktemp2[84] = elStorage->K14[2];
  ktemp2[108] = elStorage->K15[2];
  ktemp2[132] = elStorage->K16[2];
  ktemp2[14] = elStorage->K12[1];
  ktemp2[38] = elStorage->K22[2];
  ktemp2[62] = elStorage->K23[2];
  ktemp2[86] = elStorage->K24[2];
  ktemp2[110] = elStorage->K25[2];
  ktemp2[134] = elStorage->K26[2];
  ktemp2[16] = elStorage->K13[1];
  ktemp2[40] = elStorage->K23[1];
  ktemp2[64] = elStorage->K33[2];
  ktemp2[88] = elStorage->K34[2];
  ktemp2[112] = elStorage->K35[2];
  ktemp2[136] = elStorage->K36[2];
  ktemp2[18] = elStorage->K13[1];
  ktemp2[42] = elStorage->K24[1];
  ktemp2[66] = elStorage->K34[1];
  ktemp2[90] = elStorage->K44[2];
  ktemp2[114] = elStorage->K45[2];
  ktemp2[138] = elStorage->K46[2];
  ktemp2[20] = elStorage->K15[1];
  ktemp2[44] = elStorage->K25[1];
  ktemp2[68] = elStorage->K35[1];
  ktemp2[92] = elStorage->K45[1];
  ktemp2[116] = elStorage->K55[2];
  ktemp2[140] = elStorage->K56[2];
  ktemp2[22] = elStorage->K16[1];
  ktemp2[46] = elStorage->K26[1];
  ktemp2[70] = elStorage->K36[1];
  ktemp2[94] = elStorage->K46[1];
  ktemp2[118] = elStorage->K56[1];
  ktemp2[142] = elStorage->K66[2];
  ktemp2[13] = elStorage->K11[3];
  ktemp2[37] = elStorage->K12[3];
  ktemp2[61] = elStorage->K13[3];
  ktemp2[85] = elStorage->K14[3];
  ktemp2[109] = elStorage->K15[3];
  ktemp2[133] = elStorage->K16[3];
  ktemp2[15] = elStorage->K12[3];
  ktemp2[39] = elStorage->K22[3];
  ktemp2[63] = elStorage->K23[3];
  ktemp2[87] = elStorage->K24[3];
  ktemp2[111] = elStorage->K25[3];
  ktemp2[135] = elStorage->K26[3];
  ktemp2[17] = elStorage->K13[3];
  ktemp2[41] = elStorage->K23[3];
  ktemp2[65] = elStorage->K33[3];
  ktemp2[89] = elStorage->K34[3];
  ktemp2[113] = elStorage->K35[3];
  ktemp2[137] = elStorage->K36[3];
  ktemp2[19] = elStorage->K13[3];
  ktemp2[43] = elStorage->K24[3];
  ktemp2[67] = elStorage->K34[3];
  ktemp2[91] = elStorage->K44[3];
  ktemp2[115] = elStorage->K45[3];
  ktemp2[139] = elStorage->K46[3];
  ktemp2[21] = elStorage->K15[3];
  ktemp2[45] = elStorage->K25[3];
  ktemp2[69] = elStorage->K35[3];
  ktemp2[93] = elStorage->K45[3];
  ktemp2[117] = elStorage->K55[3];
  ktemp2[141] = elStorage->K56[3];
  ktemp2[23] = elStorage->K16[3];
  ktemp2[47] = elStorage->K26[3];
  ktemp2[71] = elStorage->K36[3];
  ktemp2[95] = elStorage->K46[3];
  ktemp2[119] = elStorage->K56[3];
  ktemp2[143] = elStorage->K66[3];
  mapMatrixNonSym(ktemp2, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  c_c_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * c_c_tmp;
  K16_idx_0 = elStorage->K16[0] + elStorage->S16[0] * c_c_tmp;
  ycm = elStorage->K56[0] + elStorage->S56[0] * c_c_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * c_c_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * c_c_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * c_c_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * c_c_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * c_c_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * c_c_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * c_c_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * c_c_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * c_c_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  ktemp2[0] = elStorage->K11[0] + elStorage->S11[0] * c_c_tmp;
  ktemp2[24] = (elStorage->K12[0] + S12_idx_0) + H12_idx_0;
  ktemp2[48] = (elStorage->K13[0] + S13_idx_0) + H13_idx_0;
  ktemp2_tmp = elStorage->K14[0] + S14_idx_0;
  ktemp2[72] = ktemp2_tmp + H14_idx_0;
  ktemp2[96] = K15_idx_0;
  ktemp2[120] = K16_idx_0;
  ktemp2[2] = (elStorage->K12[0] + S12_idx_0) - H12_idx_0;
  ktemp2[26] = (elStorage->K22[0] + elStorage->S22[0] * b_c_tmp) +
    SS22_data_idx_0;
  b_ktemp2_tmp = elStorage->K23[0] + S23_idx_0;
  ktemp2[50] = b_ktemp2_tmp + zcm;
  c_ktemp2_tmp = elStorage->K24[0] + S24_idx_0;
  ktemp2[74] = c_ktemp2_tmp + O1;
  d_ktemp2_tmp = elStorage->K25[0] + S25_idx_0;
  ktemp2[98] = d_ktemp2_tmp + O2;
  e_ktemp2_tmp = elStorage->K26[0] + S26_idx_0;
  ktemp2[122] = e_ktemp2_tmp + O3;
  ktemp2[4] = (elStorage->K13[0] + S13_idx_0) - H13_idx_0;
  ktemp2[28] = b_ktemp2_tmp - zcm;
  ktemp2[52] = (elStorage->K33[0] + elStorage->S33[0] * c_tmp) + SS33_data_idx_0;
  b_ktemp2_tmp = elStorage->K34[0] + S34_idx_0;
  ktemp2[76] = b_ktemp2_tmp + O1dot;
  SS33_data_idx_0 = elStorage->K35[0] + S35_idx_0;
  ktemp2[100] = SS33_data_idx_0 + O2dot;
  SS22_data_idx_0 = elStorage->K36[0] + rhoA;
  ktemp2[124] = SS22_data_idx_0 + O3dot;
  ktemp2[6] = ktemp2_tmp - H14_idx_0;
  ktemp2[30] = c_ktemp2_tmp - O1;
  ktemp2[54] = b_ktemp2_tmp - O1dot;
  ktemp2[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * b_c_tmp +
    elStorage->S44_2[0] * c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[0] + Faxial;
  ktemp2[102] = ktemp2_tmp + H45_idx_0;
  b_ktemp2_tmp = elStorage->K46[0] + integrationFactor;
  ktemp2[126] = b_ktemp2_tmp + H46_idx_0;
  ktemp2[8] = K15_idx_0;
  ktemp2[32] = d_ktemp2_tmp - O2;
  ktemp2[56] = SS33_data_idx_0 - O2dot;
  ktemp2[80] = ktemp2_tmp - H45_idx_0;
  ktemp2[104] = elStorage->K55[0] + elStorage->S55[0] * c_c_tmp;
  ktemp2[128] = ycm;
  ktemp2[10] = K16_idx_0;
  ktemp2[34] = e_ktemp2_tmp - O3;
  ktemp2[58] = SS22_data_idx_0 - O3dot;
  ktemp2[82] = b_ktemp2_tmp - H46_idx_0;
  ktemp2[106] = ycm;
  ktemp2[130] = elStorage->K66[0] + elStorage->S66[0] * c_c_tmp;
  ktemp2[1] = elStorage->K11[1] + elStorage->S11[1] * c_c_tmp;
  ktemp2[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  ktemp2[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  ktemp2_tmp = elStorage->K14[1] + S14_idx_1;
  ktemp2[73] = ktemp2_tmp + H14_idx_1;
  ktemp2[97] = K15_idx_1;
  ktemp2[121] = K16_idx_1;
  ktemp2[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  ktemp2[27] = (elStorage->K22[1] + elStorage->S22[1] * b_c_tmp) +
    SS22_data_idx_1;
  b_ktemp2_tmp = elStorage->K23[1] + S23_idx_1;
  ktemp2[51] = b_ktemp2_tmp + H23_idx_1;
  c_ktemp2_tmp = elStorage->K24[1] + S24_idx_1;
  ktemp2[75] = c_ktemp2_tmp + H24_idx_1;
  d_ktemp2_tmp = elStorage->K25[1] + S25_idx_1;
  ktemp2[99] = d_ktemp2_tmp + H25_idx_1;
  e_ktemp2_tmp = elStorage->K26[1] + S26_idx_1;
  ktemp2[123] = e_ktemp2_tmp + H26_idx_1;
  ktemp2[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  SS33_data_idx_0 = elStorage->K23[2] + S23_idx_2;
  ktemp2[29] = SS33_data_idx_0 - H23_idx_2;
  ktemp2[53] = (elStorage->K33[1] + elStorage->S33[1] * c_tmp) + SS33_data_idx_1;
  SS22_data_idx_0 = elStorage->K34[1] + S34_idx_1;
  ktemp2[77] = SS22_data_idx_0 + H34_idx_1;
  Faxial = elStorage->K35[1] + S35_idx_1;
  ktemp2[101] = Faxial + H35_idx_1;
  rhoA = elStorage->K36[1] + S36_idx_1;
  ktemp2[125] = rhoA + H36_idx_1;
  ycm = elStorage->K14[2] + S14_idx_2;
  ktemp2[7] = ycm - H14_idx_2;
  zcm = elStorage->K24[2] + S24_idx_2;
  ktemp2[31] = zcm - H24_idx_2;
  H13_idx_0 = elStorage->K34[2] + S34_idx_2;
  ktemp2[55] = H13_idx_0 - H34_idx_2;
  ktemp2[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * b_c_tmp +
    elStorage->S44_2[1] * c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  S34_idx_0 = elStorage->K45[1] + S45_idx_1;
  ktemp2[103] = S34_idx_0 + H45_idx_1;
  S35_idx_0 = elStorage->K46[1] + S46_idx_1;
  ktemp2[127] = S35_idx_0 + H46_idx_1;
  ktemp2[9] = K15_idx_2;
  S13_idx_0 = elStorage->K25[2] + S25_idx_2;
  ktemp2[33] = S13_idx_0 - H25_idx_2;
  S26_idx_0 = elStorage->K35[2] + S35_idx_2;
  ktemp2[57] = S26_idx_0 - H35_idx_2;
  S25_idx_0 = elStorage->K45[2] + S45_idx_2;
  ktemp2[81] = S25_idx_0 - H45_idx_2;
  ktemp2[105] = elStorage->K55[1] + elStorage->S55[1] * c_c_tmp;
  ktemp2[129] = K56_idx_1;
  ktemp2[11] = K16_idx_2;
  S24_idx_0 = elStorage->K26[2] + S26_idx_2;
  ktemp2[35] = S24_idx_0 - H26_idx_2;
  H12_idx_0 = elStorage->K36[2] + S36_idx_2;
  ktemp2[59] = H12_idx_0 - H36_idx_2;
  K16_idx_0 = elStorage->K46[2] + S46_idx_2;
  ktemp2[83] = K16_idx_0 - H46_idx_2;
  ktemp2[107] = K56_idx_2;
  ktemp2[131] = elStorage->K66[1] + elStorage->S66[1] * c_c_tmp;
  ktemp2[12] = elStorage->K11[2] + elStorage->S11[2] * c_c_tmp;
  ktemp2[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  ktemp2[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  ktemp2[84] = ycm + H14_idx_2;
  ktemp2[108] = K15_idx_2;
  ktemp2[132] = K16_idx_2;
  ktemp2[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  ktemp2[38] = (elStorage->K22[2] + elStorage->S22[2] * b_c_tmp) +
    SS22_data_idx_2;
  ktemp2[62] = SS33_data_idx_0 + H23_idx_2;
  ktemp2[86] = zcm + H24_idx_2;
  ktemp2[110] = S13_idx_0 + H25_idx_2;
  ktemp2[134] = S24_idx_0 + H26_idx_2;
  ktemp2[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  ktemp2[40] = b_ktemp2_tmp - H23_idx_1;
  ktemp2[64] = (elStorage->K33[2] + elStorage->S33[2] * c_tmp) + SS33_data_idx_2;
  ktemp2[88] = H13_idx_0 + H34_idx_2;
  ktemp2[112] = S26_idx_0 + H35_idx_2;
  ktemp2[136] = H12_idx_0 + H36_idx_2;
  ktemp2[18] = ktemp2_tmp - H14_idx_1;
  ktemp2[42] = c_ktemp2_tmp - H24_idx_1;
  ktemp2[66] = SS22_data_idx_0 - H34_idx_1;
  ktemp2[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * b_c_tmp +
    elStorage->S44_2[2] * c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  ktemp2[114] = S25_idx_0 + H45_idx_2;
  ktemp2[138] = K16_idx_0 + H46_idx_2;
  ktemp2[20] = K15_idx_1;
  ktemp2[44] = d_ktemp2_tmp - H25_idx_1;
  ktemp2[68] = Faxial - H35_idx_1;
  ktemp2[92] = S34_idx_0 - H45_idx_1;
  ktemp2[116] = elStorage->K55[2] + elStorage->S55[2] * c_c_tmp;
  ktemp2[140] = K56_idx_2;
  ktemp2[22] = K16_idx_1;
  ktemp2[46] = e_ktemp2_tmp - H26_idx_1;
  ktemp2[70] = rhoA - H36_idx_1;
  ktemp2[94] = S35_idx_0 - H46_idx_1;
  ktemp2[118] = K56_idx_1;
  ktemp2[142] = elStorage->K66[2] + elStorage->S66[2] * c_c_tmp;
  ktemp2[13] = elStorage->K11[3] + elStorage->S11[3] * c_c_tmp;
  ktemp2[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  ktemp2[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  ktemp2_tmp = elStorage->K14[3] + S14_idx_3;
  ktemp2[85] = ktemp2_tmp + H14_idx_3;
  ktemp2[109] = K15_idx_3;
  ktemp2[133] = K16_idx_3;
  ktemp2[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  ktemp2[39] = (elStorage->K22[3] + elStorage->S22[3] * b_c_tmp) +
    SS22_data_idx_3;
  b_ktemp2_tmp = elStorage->K23[3] + S23_idx_3;
  ktemp2[63] = b_ktemp2_tmp + H23_idx_3;
  c_ktemp2_tmp = elStorage->K24[3] + S24_idx_3;
  ktemp2[87] = c_ktemp2_tmp + H24_idx_3;
  d_ktemp2_tmp = elStorage->K25[3] + S25_idx_3;
  ktemp2[111] = d_ktemp2_tmp + H25_idx_3;
  e_ktemp2_tmp = elStorage->K26[3] + S26_idx_3;
  ktemp2[135] = e_ktemp2_tmp + H26_idx_3;
  ktemp2[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  ktemp2[41] = b_ktemp2_tmp - H23_idx_3;
  ktemp2[65] = (elStorage->K33[3] + elStorage->S33[3] * c_tmp) + SS33_data_idx_3;
  b_ktemp2_tmp = elStorage->K34[3] + S34_idx_3;
  ktemp2[89] = b_ktemp2_tmp + H34_idx_3;
  SS33_data_idx_0 = elStorage->K35[3] + S35_idx_3;
  ktemp2[113] = SS33_data_idx_0 + H35_idx_3;
  SS22_data_idx_0 = elStorage->K36[3] + S36_idx_3;
  ktemp2[137] = SS22_data_idx_0 + H36_idx_3;
  ktemp2[19] = ktemp2_tmp - H14_idx_3;
  ktemp2[43] = c_ktemp2_tmp - H24_idx_3;
  ktemp2[67] = b_ktemp2_tmp - H34_idx_3;
  ktemp2[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * b_c_tmp +
    elStorage->S44_2[3] * c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[3] + S45_idx_3;
  ktemp2[115] = ktemp2_tmp + H45_idx_3;
  b_ktemp2_tmp = elStorage->K46[3] + S46_idx_3;
  ktemp2[139] = b_ktemp2_tmp + H46_idx_3;
  ktemp2[21] = K15_idx_3;
  ktemp2[45] = d_ktemp2_tmp - H25_idx_3;
  ktemp2[69] = SS33_data_idx_0 - H35_idx_3;
  ktemp2[93] = ktemp2_tmp - H45_idx_3;
  ktemp2[117] = elStorage->K55[3] + elStorage->S55[3] * c_c_tmp;
  ktemp2[141] = K56_idx_3;
  ktemp2[23] = K16_idx_3;
  ktemp2[47] = e_ktemp2_tmp - H26_idx_3;
  ktemp2[71] = SS22_data_idx_0 - H36_idx_3;
  ktemp2[95] = b_ktemp2_tmp - H46_idx_3;
  ktemp2[119] = K56_idx_3;
  ktemp2[143] = elStorage->K66[3] + elStorage->S66[3] * c_c_tmp;
  mapMatrixNonSym(ktemp2, Ke);

  //  Declare type
  // compile Coriolis/damping matrix
  //  ktemp = [zm,C12,C13,C14,zm,zm;
  //      -C12',zm,C23,C24,C25,C26;
  //      -C13',-C23',C33,C34,C35,C36;
  //      -C14',-C24',C43,C44,C45,C46;
  //      zm,-C25',-C35',-C45',zm,zm;
  //      zm,-C26',-C36',-C46',zm,zm];
  std::memset(&ktemp2[0], 0, 144U * sizeof(double));

  //  Row 1
  ktemp2_tmp = elStorage->C12[0] * Oel[2];
  ktemp2[2] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C12[2] * Oel[2];
  ktemp2[3] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C13[0] * Oel[1];
  ktemp2[4] = c_ktemp2_tmp;
  d_ktemp2_tmp = elStorage->C13[2] * Oel[1];
  ktemp2[5] = d_ktemp2_tmp;
  e_ktemp2_tmp = elStorage->C14_1[0] * Oel[1] + elStorage->C14_2[0] * Oel[2];
  ktemp2[6] = e_ktemp2_tmp;
  SS33_data_idx_0 = elStorage->C14_1[2] * Oel[1] + elStorage->C14_2[2] * Oel[2];
  ktemp2[7] = SS33_data_idx_0;

  //  Row 2
  SS22_data_idx_0 = elStorage->C12[1] * Oel[2];
  ktemp2[14] = SS22_data_idx_0;
  Faxial = elStorage->C12[3] * Oel[2];
  ktemp2[15] = Faxial;
  rhoA = elStorage->C13[1] * Oel[1];
  ktemp2[16] = rhoA;
  ycm = elStorage->C13[3] * Oel[1];
  ktemp2[17] = ycm;
  zcm = elStorage->C14_1[1] * Oel[1] + elStorage->C14_2[1] * Oel[2];
  ktemp2[18] = zcm;
  H13_idx_0 = elStorage->C14_1[3] * Oel[1] + elStorage->C14_2[3] * Oel[2];
  ktemp2[19] = H13_idx_0;

  //  Row 3
  ktemp2[24] = -ktemp2_tmp;
  ktemp2[25] = -SS22_data_idx_0;
  ktemp2_tmp = elStorage->C23[0] * Oel[0];
  ktemp2[28] = ktemp2_tmp;
  SS22_data_idx_0 = elStorage->C23[2] * Oel[0];
  ktemp2[29] = SS22_data_idx_0;
  S34_idx_0 = elStorage->C24[0] * Oel[0];
  ktemp2[30] = S34_idx_0;
  S35_idx_0 = elStorage->C24[2] * Oel[0];
  ktemp2[31] = S35_idx_0;
  S13_idx_0 = elStorage->C25[0] * Oel[2];
  ktemp2[32] = S13_idx_0;
  S26_idx_0 = elStorage->C25[2] * Oel[2];
  ktemp2[33] = S26_idx_0;
  S25_idx_0 = elStorage->C26[0] * Oel[2];
  ktemp2[34] = S25_idx_0;
  S24_idx_0 = elStorage->C26[2] * Oel[2];
  ktemp2[35] = S24_idx_0;

  //  Row 4
  ktemp2[36] = -b_ktemp2_tmp;
  ktemp2[37] = -Faxial;
  b_ktemp2_tmp = elStorage->C23[1] * Oel[0];
  ktemp2[40] = b_ktemp2_tmp;
  Faxial = elStorage->C23[3] * Oel[0];
  ktemp2[41] = Faxial;
  H12_idx_0 = elStorage->C24[1] * Oel[0];
  ktemp2[42] = H12_idx_0;
  K16_idx_0 = elStorage->C24[3] * Oel[0];
  ktemp2[43] = K16_idx_0;
  K15_idx_0 = elStorage->C25[1] * Oel[2];
  ktemp2[44] = K15_idx_0;
  S14_idx_0 = elStorage->C25[3] * Oel[2];
  ktemp2[45] = S14_idx_0;
  S12_idx_0 = elStorage->C26[1] * Oel[2];
  ktemp2[46] = S12_idx_0;
  S23_idx_0 = elStorage->C26[3] * Oel[2];
  ktemp2[47] = S23_idx_0;

  //  Row 5
  ktemp2[48] = -c_ktemp2_tmp;
  ktemp2[49] = -rhoA;
  ktemp2[50] = -ktemp2_tmp;
  ktemp2[51] = -b_ktemp2_tmp;
  ktemp2[52] = 0.0;
  ktemp2[53] = 0.0;
  ktemp2[54] = C34_idx_0;
  ktemp2[55] = C34_idx_2;
  ktemp2_tmp = elStorage->C35[0] * Oel[1];
  ktemp2[56] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C35[2] * Oel[1];
  ktemp2[57] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C36[0] * Oel[1];
  ktemp2[58] = c_ktemp2_tmp;
  rhoA = elStorage->C36[2] * Oel[1];
  ktemp2[59] = rhoA;

  //  Row 6
  ktemp2[60] = -d_ktemp2_tmp;
  ktemp2[61] = -ycm;
  ktemp2[62] = -SS22_data_idx_0;
  ktemp2[63] = -Faxial;
  ktemp2[64] = 0.0;
  ktemp2[65] = 0.0;
  ktemp2[66] = C34_idx_1;
  ktemp2[67] = C34_idx_3;
  d_ktemp2_tmp = elStorage->C35[1] * Oel[1];
  ktemp2[68] = d_ktemp2_tmp;
  SS22_data_idx_0 = elStorage->C35[3] * Oel[1];
  ktemp2[69] = SS22_data_idx_0;
  Faxial = elStorage->C36[1] * Oel[1];
  ktemp2[70] = Faxial;
  ycm = elStorage->C36[3] * Oel[1];
  ktemp2[71] = ycm;

  //  Row 7
  ktemp2[72] = -e_ktemp2_tmp;
  ktemp2[73] = -zcm;
  ktemp2[74] = -S34_idx_0;
  ktemp2[75] = -H12_idx_0;
  ktemp2[76] = -C34_idx_0;
  ktemp2[77] = -C34_idx_1;
  ktemp2[78] = 0.0;
  ktemp2[79] = 0.0;
  e_ktemp2_tmp = elStorage->C45_1[0] * Oel[2] + elStorage->C45_2[0] * Oel[1];
  ktemp2[80] = e_ktemp2_tmp;
  zcm = elStorage->C45_1[2] * Oel[2] + elStorage->C45_2[2] * Oel[1];
  ktemp2[81] = zcm;
  S34_idx_0 = elStorage->C46_1[0] * Oel[1] + elStorage->C46_2[0] * Oel[2];
  ktemp2[82] = S34_idx_0;
  H12_idx_0 = elStorage->C46_1[2] * Oel[1] + elStorage->C46_2[2] * Oel[2];
  ktemp2[83] = H12_idx_0;

  //  Row 8
  ktemp2[84] = -SS33_data_idx_0;
  ktemp2[85] = -H13_idx_0;
  ktemp2[86] = -S35_idx_0;
  ktemp2[87] = -K16_idx_0;
  ktemp2[88] = -C34_idx_2;
  ktemp2[89] = -C34_idx_3;
  ktemp2[90] = 0.0;
  ktemp2[91] = 0.0;
  SS33_data_idx_0 = elStorage->C45_1[1] * Oel[2] + elStorage->C45_2[1] * Oel[1];
  ktemp2[92] = SS33_data_idx_0;
  H13_idx_0 = elStorage->C45_1[3] * Oel[2] + elStorage->C45_2[3] * Oel[1];
  ktemp2[93] = H13_idx_0;
  S35_idx_0 = elStorage->C46_1[1] * Oel[1] + elStorage->C46_2[1] * Oel[2];
  ktemp2[94] = S35_idx_0;
  K16_idx_0 = elStorage->C46_1[3] * Oel[1] + elStorage->C46_2[3] * Oel[2];
  ktemp2[95] = K16_idx_0;

  //  Row 9
  ktemp2[98] = -S13_idx_0;
  ktemp2[99] = -K15_idx_0;
  ktemp2[100] = -ktemp2_tmp;
  ktemp2[101] = -d_ktemp2_tmp;
  ktemp2[102] = -e_ktemp2_tmp;
  ktemp2[103] = -SS33_data_idx_0;

  //  Row 10
  ktemp2[110] = -S26_idx_0;
  ktemp2[111] = -S14_idx_0;
  ktemp2[112] = -b_ktemp2_tmp;
  ktemp2[113] = -SS22_data_idx_0;
  ktemp2[114] = -zcm;
  ktemp2[115] = -H13_idx_0;

  //  Row 11
  ktemp2[122] = -S25_idx_0;
  ktemp2[123] = -S12_idx_0;
  ktemp2[124] = -c_ktemp2_tmp;
  ktemp2[125] = -Faxial;
  ktemp2[126] = -S34_idx_0;
  ktemp2[127] = -S35_idx_0;

  //  Row 12
  ktemp2[134] = -S24_idx_0;
  ktemp2[135] = -S23_idx_0;
  ktemp2[136] = -rhoA;
  ktemp2[137] = -ycm;
  ktemp2[138] = -H12_idx_0;
  ktemp2[139] = -K16_idx_0;
  mapMatrixNonSym(ktemp2, Ce);

  // compile mass matrix
  ktemp2[0] = elStorage->M11[0];
  ktemp2[24] = 0.0;
  ktemp2[48] = 0.0;
  ktemp2[72] = 0.0;
  ktemp2[96] = elStorage->M15[0];
  ktemp2[120] = elStorage->M16[0];
  ktemp2[2] = 0.0;
  ktemp2[26] = elStorage->M22[0];
  ktemp2[50] = 0.0;
  ktemp2[74] = elStorage->M24[0];
  ktemp2[98] = 0.0;
  ktemp2[122] = 0.0;
  ktemp2[4] = 0.0;
  ktemp2[28] = 0.0;
  ktemp2[52] = elStorage->M33[0];
  ktemp2[76] = elStorage->M34[0];
  ktemp2[100] = 0.0;
  ktemp2[124] = 0.0;
  ktemp2[6] = 0.0;
  ktemp2[30] = elStorage->M24[0];
  ktemp2[54] = elStorage->M34[0];
  ktemp2[78] = elStorage->M44[0];
  ktemp2[102] = 0.0;
  ktemp2[126] = 0.0;
  ktemp2[8] = elStorage->M15[0];
  ktemp2[32] = 0.0;
  ktemp2[56] = 0.0;
  ktemp2[80] = 0.0;
  ktemp2[104] = elStorage->M55[0];
  ktemp2[128] = elStorage->M56[0];
  ktemp2[10] = elStorage->M16[0];
  ktemp2[34] = 0.0;
  ktemp2[58] = 0.0;
  ktemp2[82] = 0.0;
  ktemp2[106] = elStorage->M56[0];
  ktemp2[130] = elStorage->M66[0];
  ktemp2[1] = elStorage->M11[1];
  ktemp2[25] = 0.0;
  ktemp2[49] = 0.0;
  ktemp2[73] = 0.0;
  ktemp2[97] = elStorage->M15[1];
  ktemp2[121] = elStorage->M16[1];
  ktemp2[3] = 0.0;
  ktemp2[27] = elStorage->M22[1];
  ktemp2[51] = 0.0;
  ktemp2[75] = elStorage->M24[1];
  ktemp2[99] = 0.0;
  ktemp2[123] = 0.0;
  ktemp2[5] = 0.0;
  ktemp2[29] = 0.0;
  ktemp2[53] = elStorage->M33[1];
  ktemp2[77] = elStorage->M34[1];
  ktemp2[101] = 0.0;
  ktemp2[125] = 0.0;
  ktemp2[7] = 0.0;
  ktemp2[31] = elStorage->M24[2];
  ktemp2[55] = elStorage->M34[2];
  ktemp2[79] = elStorage->M44[1];
  ktemp2[103] = 0.0;
  ktemp2[127] = 0.0;
  ktemp2[9] = elStorage->M15[2];
  ktemp2[33] = 0.0;
  ktemp2[57] = 0.0;
  ktemp2[81] = 0.0;
  ktemp2[105] = elStorage->M55[1];
  ktemp2[129] = elStorage->M56[1];
  ktemp2[11] = elStorage->M16[2];
  ktemp2[35] = 0.0;
  ktemp2[59] = 0.0;
  ktemp2[83] = 0.0;
  ktemp2[107] = elStorage->M56[2];
  ktemp2[131] = elStorage->M66[1];
  ktemp2[12] = elStorage->M11[2];
  ktemp2[36] = 0.0;
  ktemp2[60] = 0.0;
  ktemp2[84] = 0.0;
  ktemp2[108] = elStorage->M15[2];
  ktemp2[132] = elStorage->M16[2];
  ktemp2[14] = 0.0;
  ktemp2[38] = elStorage->M22[2];
  ktemp2[62] = 0.0;
  ktemp2[86] = elStorage->M24[2];
  ktemp2[110] = 0.0;
  ktemp2[134] = 0.0;
  ktemp2[16] = 0.0;
  ktemp2[40] = 0.0;
  ktemp2[64] = elStorage->M33[2];
  ktemp2[88] = elStorage->M34[2];
  ktemp2[112] = 0.0;
  ktemp2[136] = 0.0;
  ktemp2[18] = 0.0;
  ktemp2[42] = elStorage->M24[1];
  ktemp2[66] = elStorage->M34[1];
  ktemp2[90] = elStorage->M44[2];
  ktemp2[114] = 0.0;
  ktemp2[138] = 0.0;
  ktemp2[20] = elStorage->M15[1];
  ktemp2[44] = 0.0;
  ktemp2[68] = 0.0;
  ktemp2[92] = 0.0;
  ktemp2[116] = elStorage->M55[2];
  ktemp2[140] = elStorage->M56[2];
  ktemp2[22] = elStorage->M16[1];
  ktemp2[46] = 0.0;
  ktemp2[70] = 0.0;
  ktemp2[94] = 0.0;
  ktemp2[118] = elStorage->M56[1];
  ktemp2[142] = elStorage->M66[2];
  ktemp2[13] = elStorage->M11[3];
  ktemp2[37] = 0.0;
  ktemp2[61] = 0.0;
  ktemp2[85] = 0.0;
  ktemp2[109] = elStorage->M15[3];
  ktemp2[133] = elStorage->M16[3];
  ktemp2[15] = 0.0;
  ktemp2[39] = elStorage->M22[3];
  ktemp2[63] = 0.0;
  ktemp2[87] = elStorage->M24[3];
  ktemp2[111] = 0.0;
  ktemp2[135] = 0.0;
  ktemp2[17] = 0.0;
  ktemp2[41] = 0.0;
  ktemp2[65] = elStorage->M33[3];
  ktemp2[89] = elStorage->M34[3];
  ktemp2[113] = 0.0;
  ktemp2[137] = 0.0;
  ktemp2[19] = 0.0;
  ktemp2[43] = elStorage->M24[3];
  ktemp2[67] = elStorage->M34[3];
  ktemp2[91] = elStorage->M44[3];
  ktemp2[115] = 0.0;
  ktemp2[139] = 0.0;
  ktemp2[21] = elStorage->M15[3];
  ktemp2[45] = 0.0;
  ktemp2[69] = 0.0;
  ktemp2[93] = 0.0;
  ktemp2[117] = elStorage->M55[3];
  ktemp2[141] = elStorage->M56[3];
  ktemp2[23] = elStorage->M16[3];
  ktemp2[47] = 0.0;
  ktemp2[71] = 0.0;
  ktemp2[95] = 0.0;
  ktemp2[119] = elStorage->M56[3];
  ktemp2[143] = elStorage->M66[3];
  mapMatrixNonSym(ktemp2, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input_RayleighAlpha * Kenr[i] + input_RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      ktemp2[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(ktemp2, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ke, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Ke);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fel_data[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fel_data[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fel_data,
                  dispLocal);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input_concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input_concMass[0];
    Me[13] += input_concMass[0];
    Me[26] += input_concMass[0];
    Me[39] += input_concMass[1];
    Me[52] += input_concMass[2];
    Me[65] += input_concMass[3];
    Me[78] += input_concMass[4];
    Me[91] += input_concMass[4];
    Me[104] += input_concMass[4];
    Me[117] += input_concMass[5];
    Me[130] += input_concMass[6];
    Me[143] += input_concMass[7];

    // modify Ce for concentrated mass
    K15_idx_0 = 2.0 * input_concMass[0] * 3.2898681336964524;
    Ce[12] -= K15_idx_0;
    Ce[1] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[0] * 0.0;
    Ce[24] += K15_idx_0;
    Ce[2] -= K15_idx_0;
    Ce[25] -= K15_idx_0;
    Ce[14] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[4] * 3.2898681336964524;
    Ce[90] -= K15_idx_0;
    Ce[79] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[4] * 0.0;
    Ce[102] += K15_idx_0;
    Ce[80] -= K15_idx_0;
    Ce[103] -= K15_idx_0;
    Ce[92] += K15_idx_0;

    // modify Ke for concentrated mass
    Ke[0] -= input_concMass[0] * 10.823232337111378;
    K15_idx_0 = input_concMass[0] * 0.0 * 0.0;
    Ke[12] = (Ke[12] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[1] = (Ke[1] + K15_idx_0) + input_concMass[0] * 0.0;
    K15_idx_0 = input_concMass[0] * 0.0 * 3.2898681336964524;
    Ke[24] = (Ke[24] + K15_idx_0) + input_concMass[0] * 0.0;
    Ke[2] = (Ke[2] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[25] = (Ke[25] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[14] = (Ke[14] + K15_idx_0) + input_concMass[0] * 0.0;
    Ke[13] -= input_concMass[0] * 10.823232337111378;
    Ke[26] -= input_concMass[0] * 0.0;
    Ke[78] -= input_concMass[4] * 10.823232337111378;
    K15_idx_0 = input_concMass[4] * 0.0 * 0.0;
    Ke[90] = (Ke[90] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[79] = (Ke[79] + K15_idx_0) + input_concMass[4] * 0.0;
    K15_idx_0 = input_concMass[4] * 0.0 * 3.2898681336964524;
    Ke[102] = (Ke[102] + K15_idx_0) + input_concMass[4] * 0.0;
    Ke[80] = (Ke[80] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[103] = (Ke[103] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[92] = (Ke[92] + K15_idx_0) + input_concMass[4] * 0.0;
    Ke[91] -= input_concMass[4] * 10.823232337111378;
    Ke[104] -= input_concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    dispLocal[0] = ((dispLocal[0] + input_concMass[0] * ((input_x_data[0] *
      10.823232337111378 - 0.0 * input_y_data[0]) - 0.0 * input_z_data[0])) +
                    input_concMass[0] * (input_y_data[0] * 0.0 - input_z_data[0]
      * 0.0)) - input_concMass[0] * a_temp[0];
    K15_idx_0 = input_z_data[0] * 0.0 - input_x_data[0] * 0.0;
    dispLocal[1] = ((dispLocal[1] + input_concMass[0] * ((input_y_data[0] *
      10.823232337111378 - 0.0 * input_z_data[0]) - 0.0 * input_x_data[0])) +
                    input_concMass[0] * K15_idx_0) - input_concMass[0] * a_temp
      [1];
    dispLocal[2] = ((dispLocal[2] + input_concMass[0] * (K15_idx_0 - 0.0 *
      input_y_data[0])) + input_concMass[0] * (input_x_data[0] * 0.0 -
      input_y_data[0] * 0.0)) - input_concMass[0] * a_temp[2];
    dispLocal[6] = ((dispLocal[6] + input_concMass[4] * ((input_x_data[1] *
      10.823232337111378 - 0.0 * input_y_data[1]) - 0.0 * input_z_data[1])) +
                    input_concMass[4] * (input_y_data[1] * 0.0 - input_z_data[1]
      * 0.0)) - input_concMass[4] * a_temp[0];
    K15_idx_0 = input_z_data[1] * 0.0 - input_x_data[1] * 0.0;
    dispLocal[7] = ((dispLocal[7] + input_concMass[4] * ((input_y_data[1] *
      10.823232337111378 - 0.0 * input_z_data[1]) - 0.0 * input_x_data[1])) +
                    input_concMass[4] * K15_idx_0) - input_concMass[4] * a_temp
      [1];
    dispLocal[8] = ((dispLocal[8] + input_concMass[4] * (K15_idx_0 - 0.0 *
      input_y_data[1])) + input_concMass[4] * (input_x_data[1] * 0.0 -
      input_y_data[1] * 0.0)) - input_concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Ke[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &dispLocal[0], 12U * sizeof(double));
  output->Me.size[0] = 12;
  output->Me.size[1] = 12;
  output->Ce.size[0] = 12;
  output->Ce.size[1] = 12;
  std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
  std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));

  // ------------------------------------------
}

//
// calculateTimoshenkoElementNL performs nonlinear element calculations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [output] = calculateTimoshenkoElementNL(input,elStorage)
//
//    This function performs nonlinear element calculations.
//
//       input:
//       input      = object containing element input
//       elStorage  = obect containing precalculated element data
//
//       output:
//       output     = object containing element data
// Arguments    : const n_struct_T *input
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void calculateTimoshenkoElementNL(const n_struct_T *input, const f_struct_T
  *elStorage, o_struct_T *output)
{
  int dispdot_size_idx_1;
  double dispdot_data[12];
  int dispddot_size_idx_1;
  double dispddot_data[12];
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double Omega;
  double OmegaDot;
  double lambda[144];
  int i;
  double O1;
  double Oel[3];
  double K15_idx_0;
  double O2;
  double H34_idx_3;
  double O3;
  double O1dot;
  double ODotel[3];
  double a_temp[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double H34_idx_0;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double integrationFactor;
  double c_tmp;
  double b_c_tmp;
  double S12_idx_0;
  double S13_idx_0;
  double S23_idx_0;
  double S25_idx_0;
  double S26_idx_0;
  double S35_idx_0;
  double S14_idx_0;
  double S24_idx_0;
  double S34_idx_0;
  double C34_idx_0;
  double H12_idx_0;
  double zcm;
  double rhoA;
  double ycm;
  double H35_idx_0;
  double H36_idx_0;
  double H14_idx_0;
  double disMomentgp[3];
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double posLocal[3];
  double S35_idx_1;
  double disLoadgpLocal[3];
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C34_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double H23_idx_1;
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C34_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C34_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double H46_idx_3;
  double ktemp2[144];
  double Kenr[144];
  double c_c_tmp;
  double K15_idx_1;
  double K16_idx_1;
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double ktemp2_tmp;
  double b_ktemp2_tmp;
  double c_ktemp2_tmp;
  double d_ktemp2_tmp;
  double e_ktemp2_tmp;
  double f_ktemp2_tmp;
  double g_ktemp2_tmp;
  double Khate[144];
  double Ce[144];
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fhate[12];
  double Fe[12];
  int tmp_size[1];
  boolean_T concMassFlag;
  double FhatLessConc[12];
  double b_data[12];

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  // declare type
  dispdot_size_idx_1 = 1;
  dispdot_data[0] = 0.0;

  // declare type
  dispddot_size_idx_1 = 1;
  dispddot_data[0] = 0.0;

  // declare type
  // options for Dean integrator
  if (c_strcmp(input->analysisType)) {
    // options for newmark beta integrator
    dispdot_size_idx_1 = 12;
    dispddot_size_idx_1 = 12;
    std::memcpy(&dispdot_data[0], &input->dispdot.data[0], 12U * sizeof(double));
    std::memcpy(&dispddot_data[0], &input->dispddot.data[0], 12U * sizeof(double));
  }

  // --------------------------------------------
  // setting for modal analysis flag
  // setting for initial reduced order model calculations
  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  Omega = 6.2831853071795862 * input->Omega;
  OmegaDot = 6.2831853071795862 * input->OmegaDot;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input->sweepAngle * 3.1415926535897931 / 180.0,
                  input->coneAngle * 3.1415926535897931 / 180.0,
                  (input->rollAngle + 0.5 * (input->sectionProps.twist[0] +
    input->sectionProps.twist[1])) * 3.1415926535897931 / 180.0, lambda);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    K15_idx_0 = lambda[i + 24];
    H34_idx_3 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    Oel[i] = H34_idx_3 + K15_idx_0 * Omega;
    a_temp[i] = (input->CN2H[i] * 0.0 + input->CN2H[i + 3] * 0.0) + input->
      CN2H[i + 6] * 9.81;
    ODotel[i] = H34_idx_3 + K15_idx_0 * OmegaDot;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  H34_idx_0 = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input->xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &H34_idx_3);
    integrationFactor = H34_idx_3 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input->sectionProps.rhoA[0] + N_data[1] *
      input->sectionProps.rhoA[1];
    ycm = N_data[0] * input->sectionProps.ycm[0] + N_data[1] *
      input->sectionProps.ycm[1];
    zcm = N_data[0] * input->sectionProps.zcm[0] + N_data[1] *
      input->sectionProps.zcm[1];
    S14_idx_0 = N_data[0] * input->x.data[0] + N_data[1] * input->x.data[1];
    S12_idx_0 = N_data[0] * input->y.data[0] + N_data[1] * input->y.data[1];
    H34_idx_3 = N_data[0] * input->z.data[0] + N_data[1] * input->z.data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      K15_idx_0 = lambda[i + 12];
      S23_idx_0 = lambda[i + 24];
      posLocal[i] = (lambda[i] * S14_idx_0 + K15_idx_0 * S12_idx_0) + S23_idx_0 *
        H34_idx_3;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + K15_idx_0 * disMomentgp
                           [1]) + S23_idx_0 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // calculate static aerodynamic load
    // distributed/body force load calculations
    K15_idx_0 = (O2 * O2 + O3 * O3) * posLocal[0];
    H12_idx_0 = O2dot * posLocal[2];
    S24_idx_0 = O3dot * posLocal[1];
    S25_idx_0 = H34_idx_0 * posLocal[2];
    S26_idx_0 = O1 * O2 * posLocal[1];
    S13_idx_0 = rhoA * ((((K15_idx_0 - S26_idx_0) - S25_idx_0) + S24_idx_0) -
                        H12_idx_0) - disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    H34_idx_3 = O1dot * posLocal[2];
    S23_idx_0 = O3dot * posLocal[0];
    S35_idx_0 = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 *
      O3) - posLocal[0] * O1 * O2) + H34_idx_3) - S23_idx_0) - disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    S14_idx_0 = O2dot * posLocal[0];
    S12_idx_0 = O1dot * posLocal[1];
    S34_idx_0 = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - H34_idx_0 *
      posLocal[0]) - O2 * O3 * posLocal[1]) + S14_idx_0) - S12_idx_0) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    S14_idx_0 = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) -
      posLocal[1] * (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                          (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                         (S14_idx_0 - S12_idx_0)) - zcm * (H34_idx_3 - S23_idx_0))
      - disMomentgp[0];

    // This function is a general routine to calculate an element vector
    H34_idx_3 = rhoA * zcm * ((((K15_idx_0 - posLocal[1] * O1 * O2) - posLocal[2]
      * O1 * O3) - H12_idx_0) + S24_idx_0) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    S23_idx_0 = rhoA * ycm * ((((S25_idx_0 + S26_idx_0) - K15_idx_0) - S24_idx_0)
      + H12_idx_0) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    F1_data_idx_0 += S13_idx_0 * N_data[0] * integrationFactor;
    F2_data_idx_0 += S35_idx_0 * N_data[0] * integrationFactor;
    F3_data_idx_0 += S34_idx_0 * N_data[0] * integrationFactor;
    F4_data_idx_0 += S14_idx_0 * N_data[0] * integrationFactor;
    F5_data_idx_0 += H34_idx_3 * N_data[0] * integrationFactor;
    F6_data_idx_0 += S23_idx_0 * N_data[0] * integrationFactor;
    F1_data_idx_1 += S13_idx_0 * N_data[1] * integrationFactor;
    F2_data_idx_1 += S35_idx_0 * N_data[1] * integrationFactor;
    F3_data_idx_1 += S34_idx_0 * N_data[1] * integrationFactor;
    F4_data_idx_1 += S14_idx_0 * N_data[1] * integrationFactor;
    F5_data_idx_1 += H34_idx_3 * N_data[1] * integrationFactor;
    F6_data_idx_1 += S23_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  H34_idx_3 = Oel[0] * Oel[1];
  c_tmp = Oel[0] * Oel[0];
  b_c_tmp = c_tmp + Oel[2] * Oel[2];
  c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  S12_idx_0 = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  S23_idx_0 = elStorage->S23[0] * Oel[1] * Oel[2];
  S25_idx_0 = elStorage->S25[0] * H34_idx_3;
  S26_idx_0 = elStorage->S26[0] * H34_idx_3;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  O1 = elStorage->S36[0] * Oel[0] * Oel[2];
  S14_idx_0 = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel
    [0] * Oel[1];
  S24_idx_0 = elStorage->S24_1[0] * b_c_tmp + elStorage->S24_2[0] * Oel[1] *
    Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * c_tmp + elStorage->S34_2[0] * Oel[1] * Oel[2];
  O2 = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel[0] *
    Oel[1];
  O3 = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0] * Oel[0] *
    Oel[2];
  K15_idx_0 = elStorage->C34[0];
  C34_idx_0 = K15_idx_0 * Oel[0];
  H12_idx_0 = 0.5 * elStorage->C12[0] * ODotel[2];
  zcm = 0.5 * elStorage->C13[0] * ODotel[1];
  rhoA = 0.5 * elStorage->C23[0] * ODotel[0];
  O1dot = 0.5 * elStorage->C24[0] * ODotel[0];
  O2dot = 0.5 * elStorage->C25[0] * ODotel[2];
  O3dot = 0.5 * elStorage->C26[0] * ODotel[2];
  H34_idx_0 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_0 = 0.5 * elStorage->C35[0] * ODotel[1];
  H36_idx_0 = 0.5 * elStorage->C36[0] * ODotel[1];
  H14_idx_0 = 0.5 * (elStorage->C14_1[0] * ODotel[1] + elStorage->C14_2[0] *
                     ODotel[2]);
  H45_idx_0 = 0.5 * (elStorage->C45_1[0] * ODotel[2] + elStorage->C45_2[0] *
                     ODotel[1]);
  H46_idx_0 = 0.5 * (elStorage->C46_1[0] * ODotel[1] + elStorage->C46_2[0] *
                     ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * H34_idx_3;
  S26_idx_1 = elStorage->S26[1] * H34_idx_3;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * b_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * c_tmp + elStorage->S34_2[1] * Oel[1] * Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[1];
  C34_idx_1 = K15_idx_0 * Oel[0];
  H12_idx_1 = 0.5 * elStorage->C12[1] * ODotel[2];
  H13_idx_1 = 0.5 * elStorage->C13[1] * ODotel[1];
  H23_idx_1 = 0.5 * elStorage->C23[1] * ODotel[0];
  H24_idx_1 = 0.5 * elStorage->C24[1] * ODotel[0];
  H25_idx_1 = 0.5 * elStorage->C25[1] * ODotel[2];
  H26_idx_1 = 0.5 * elStorage->C26[1] * ODotel[2];
  H34_idx_1 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_1 = 0.5 * elStorage->C35[1] * ODotel[1];
  H36_idx_1 = 0.5 * elStorage->C36[1] * ODotel[1];
  H14_idx_1 = 0.5 * (elStorage->C14_1[1] * ODotel[1] + elStorage->C14_2[1] *
                     ODotel[2]);
  H45_idx_1 = 0.5 * (elStorage->C45_1[1] * ODotel[2] + elStorage->C45_2[1] *
                     ODotel[1]);
  H46_idx_1 = 0.5 * (elStorage->C46_1[1] * ODotel[1] + elStorage->C46_2[1] *
                     ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * H34_idx_3;
  S26_idx_2 = elStorage->S26[2] * H34_idx_3;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * b_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * c_tmp + elStorage->S34_2[2] * Oel[1] * Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[2];
  C34_idx_2 = K15_idx_0 * Oel[0];
  H12_idx_2 = 0.5 * elStorage->C12[2] * ODotel[2];
  H13_idx_2 = 0.5 * elStorage->C13[2] * ODotel[1];
  H23_idx_2 = 0.5 * elStorage->C23[2] * ODotel[0];
  H24_idx_2 = 0.5 * elStorage->C24[2] * ODotel[0];
  H25_idx_2 = 0.5 * elStorage->C25[2] * ODotel[2];
  H26_idx_2 = 0.5 * elStorage->C26[2] * ODotel[2];
  H34_idx_2 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_2 = 0.5 * elStorage->C35[2] * ODotel[1];
  H36_idx_2 = 0.5 * elStorage->C36[2] * ODotel[1];
  H14_idx_2 = 0.5 * (elStorage->C14_1[2] * ODotel[1] + elStorage->C14_2[2] *
                     ODotel[2]);
  H45_idx_2 = 0.5 * (elStorage->C45_1[2] * ODotel[2] + elStorage->C45_2[2] *
                     ODotel[1]);
  H46_idx_2 = 0.5 * (elStorage->C46_1[2] * ODotel[1] + elStorage->C46_2[2] *
                     ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * H34_idx_3;
  S26_idx_3 = elStorage->S26[3] * H34_idx_3;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * b_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * c_tmp + elStorage->S34_2[3] * Oel[1] * Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  K15_idx_0 = elStorage->C34[3];
  C34_idx_3 = K15_idx_0 * Oel[0];
  H12_idx_3 = 0.5 * elStorage->C12[3] * ODotel[2];
  H13_idx_3 = 0.5 * elStorage->C13[3] * ODotel[1];
  H23_idx_3 = 0.5 * elStorage->C23[3] * ODotel[0];
  H24_idx_3 = 0.5 * elStorage->C24[3] * ODotel[0];
  H25_idx_3 = 0.5 * elStorage->C25[3] * ODotel[2];
  H26_idx_3 = 0.5 * elStorage->C26[3] * ODotel[2];
  H34_idx_3 = 0.5 * K15_idx_0 * ODotel[0];
  H35_idx_3 = 0.5 * elStorage->C35[3] * ODotel[1];
  H36_idx_3 = 0.5 * elStorage->C36[3] * ODotel[1];
  H14_idx_3 = 0.5 * (elStorage->C14_1[3] * ODotel[1] + elStorage->C14_2[3] *
                     ODotel[2]);
  H45_idx_3 = 0.5 * (elStorage->C45_1[3] * ODotel[2] + elStorage->C45_2[3] *
                     ODotel[1]);
  H46_idx_3 = 0.5 * (elStorage->C46_1[3] * ODotel[1] + elStorage->C46_2[3] *
                     ODotel[2]);

  // compile stiffness matrix without rotational effects
  ktemp2[0] = elStorage->K11[0];
  ktemp2[24] = elStorage->K12[0];
  ktemp2[48] = elStorage->K13[0];
  ktemp2[72] = elStorage->K14[0];
  ktemp2[96] = elStorage->K15[0];
  ktemp2[120] = elStorage->K16[0];
  ktemp2[2] = elStorage->K12[0];
  ktemp2[26] = elStorage->K22[0];
  ktemp2[50] = elStorage->K23[0];
  ktemp2[74] = elStorage->K24[0];
  ktemp2[98] = elStorage->K25[0];
  ktemp2[122] = elStorage->K26[0];
  ktemp2[4] = elStorage->K13[0];
  ktemp2[28] = elStorage->K23[0];
  ktemp2[52] = elStorage->K33[0];
  ktemp2[76] = elStorage->K34[0];
  ktemp2[100] = elStorage->K35[0];
  ktemp2[124] = elStorage->K36[0];
  ktemp2[6] = elStorage->K13[0];
  ktemp2[30] = elStorage->K24[0];
  ktemp2[54] = elStorage->K34[0];
  ktemp2[78] = elStorage->K44[0];
  ktemp2[102] = elStorage->K45[0];
  ktemp2[126] = elStorage->K46[0];
  ktemp2[8] = elStorage->K15[0];
  ktemp2[32] = elStorage->K25[0];
  ktemp2[56] = elStorage->K35[0];
  ktemp2[80] = elStorage->K45[0];
  ktemp2[104] = elStorage->K55[0];
  ktemp2[128] = elStorage->K56[0];
  ktemp2[10] = elStorage->K16[0];
  ktemp2[34] = elStorage->K26[0];
  ktemp2[58] = elStorage->K36[0];
  ktemp2[82] = elStorage->K46[0];
  ktemp2[106] = elStorage->K56[0];
  ktemp2[130] = elStorage->K66[0];
  ktemp2[1] = elStorage->K11[1];
  ktemp2[25] = elStorage->K12[1];
  ktemp2[49] = elStorage->K13[1];
  ktemp2[73] = elStorage->K14[1];
  ktemp2[97] = elStorage->K15[1];
  ktemp2[121] = elStorage->K16[1];
  ktemp2[3] = elStorage->K12[2];
  ktemp2[27] = elStorage->K22[1];
  ktemp2[51] = elStorage->K23[1];
  ktemp2[75] = elStorage->K24[1];
  ktemp2[99] = elStorage->K25[1];
  ktemp2[123] = elStorage->K26[1];
  ktemp2[5] = elStorage->K13[2];
  ktemp2[29] = elStorage->K23[2];
  ktemp2[53] = elStorage->K33[1];
  ktemp2[77] = elStorage->K34[1];
  ktemp2[101] = elStorage->K35[1];
  ktemp2[125] = elStorage->K36[1];
  ktemp2[7] = elStorage->K13[2];
  ktemp2[31] = elStorage->K24[2];
  ktemp2[55] = elStorage->K34[2];
  ktemp2[79] = elStorage->K44[1];
  ktemp2[103] = elStorage->K45[1];
  ktemp2[127] = elStorage->K46[1];
  ktemp2[9] = elStorage->K15[2];
  ktemp2[33] = elStorage->K25[2];
  ktemp2[57] = elStorage->K35[2];
  ktemp2[81] = elStorage->K45[2];
  ktemp2[105] = elStorage->K55[1];
  ktemp2[129] = elStorage->K56[1];
  ktemp2[11] = elStorage->K16[2];
  ktemp2[35] = elStorage->K26[2];
  ktemp2[59] = elStorage->K36[2];
  ktemp2[83] = elStorage->K46[2];
  ktemp2[107] = elStorage->K56[2];
  ktemp2[131] = elStorage->K66[1];
  ktemp2[12] = elStorage->K11[2];
  ktemp2[36] = elStorage->K12[2];
  ktemp2[60] = elStorage->K13[2];
  ktemp2[84] = elStorage->K14[2];
  ktemp2[108] = elStorage->K15[2];
  ktemp2[132] = elStorage->K16[2];
  ktemp2[14] = elStorage->K12[1];
  ktemp2[38] = elStorage->K22[2];
  ktemp2[62] = elStorage->K23[2];
  ktemp2[86] = elStorage->K24[2];
  ktemp2[110] = elStorage->K25[2];
  ktemp2[134] = elStorage->K26[2];
  ktemp2[16] = elStorage->K13[1];
  ktemp2[40] = elStorage->K23[1];
  ktemp2[64] = elStorage->K33[2];
  ktemp2[88] = elStorage->K34[2];
  ktemp2[112] = elStorage->K35[2];
  ktemp2[136] = elStorage->K36[2];
  ktemp2[18] = elStorage->K13[1];
  ktemp2[42] = elStorage->K24[1];
  ktemp2[66] = elStorage->K34[1];
  ktemp2[90] = elStorage->K44[2];
  ktemp2[114] = elStorage->K45[2];
  ktemp2[138] = elStorage->K46[2];
  ktemp2[20] = elStorage->K15[1];
  ktemp2[44] = elStorage->K25[1];
  ktemp2[68] = elStorage->K35[1];
  ktemp2[92] = elStorage->K45[1];
  ktemp2[116] = elStorage->K55[2];
  ktemp2[140] = elStorage->K56[2];
  ktemp2[22] = elStorage->K16[1];
  ktemp2[46] = elStorage->K26[1];
  ktemp2[70] = elStorage->K36[1];
  ktemp2[94] = elStorage->K46[1];
  ktemp2[118] = elStorage->K56[1];
  ktemp2[142] = elStorage->K66[2];
  ktemp2[13] = elStorage->K11[3];
  ktemp2[37] = elStorage->K12[3];
  ktemp2[61] = elStorage->K13[3];
  ktemp2[85] = elStorage->K14[3];
  ktemp2[109] = elStorage->K15[3];
  ktemp2[133] = elStorage->K16[3];
  ktemp2[15] = elStorage->K12[3];
  ktemp2[39] = elStorage->K22[3];
  ktemp2[63] = elStorage->K23[3];
  ktemp2[87] = elStorage->K24[3];
  ktemp2[111] = elStorage->K25[3];
  ktemp2[135] = elStorage->K26[3];
  ktemp2[17] = elStorage->K13[3];
  ktemp2[41] = elStorage->K23[3];
  ktemp2[65] = elStorage->K33[3];
  ktemp2[89] = elStorage->K34[3];
  ktemp2[113] = elStorage->K35[3];
  ktemp2[137] = elStorage->K36[3];
  ktemp2[19] = elStorage->K13[3];
  ktemp2[43] = elStorage->K24[3];
  ktemp2[67] = elStorage->K34[3];
  ktemp2[91] = elStorage->K44[3];
  ktemp2[115] = elStorage->K45[3];
  ktemp2[139] = elStorage->K46[3];
  ktemp2[21] = elStorage->K15[3];
  ktemp2[45] = elStorage->K25[3];
  ktemp2[69] = elStorage->K35[3];
  ktemp2[93] = elStorage->K45[3];
  ktemp2[117] = elStorage->K55[3];
  ktemp2[141] = elStorage->K56[3];
  ktemp2[23] = elStorage->K16[3];
  ktemp2[47] = elStorage->K26[3];
  ktemp2[71] = elStorage->K36[3];
  ktemp2[95] = elStorage->K46[3];
  ktemp2[119] = elStorage->K56[3];
  ktemp2[143] = elStorage->K66[3];
  mapMatrixNonSym(ktemp2, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  c_c_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * c_c_tmp;
  ycm = elStorage->K16[0] + elStorage->S16[0] * c_c_tmp;
  integrationFactor = elStorage->K56[0] + elStorage->S56[0] * c_c_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * c_c_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * c_c_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * c_c_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * c_c_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * c_c_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * c_c_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * c_c_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * c_c_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * c_c_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  ktemp2[0] = elStorage->K11[0] + elStorage->S11[0] * c_c_tmp;
  ktemp2[24] = (elStorage->K12[0] + S12_idx_0) + H12_idx_0;
  ktemp2[48] = (elStorage->K13[0] + S13_idx_0) + zcm;
  ktemp2_tmp = elStorage->K14[0] + S14_idx_0;
  ktemp2[72] = ktemp2_tmp + H14_idx_0;
  ktemp2[96] = K15_idx_0;
  ktemp2[120] = ycm;
  ktemp2[2] = (elStorage->K12[0] + S12_idx_0) - H12_idx_0;
  ktemp2[26] = elStorage->K22[0] + elStorage->S22[0] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[0] + S23_idx_0;
  ktemp2[50] = b_ktemp2_tmp + rhoA;
  c_ktemp2_tmp = elStorage->K24[0] + S24_idx_0;
  ktemp2[74] = c_ktemp2_tmp + O1dot;
  d_ktemp2_tmp = elStorage->K25[0] + S25_idx_0;
  ktemp2[98] = d_ktemp2_tmp + O2dot;
  e_ktemp2_tmp = elStorage->K26[0] + S26_idx_0;
  ktemp2[122] = e_ktemp2_tmp + O3dot;
  ktemp2[4] = (elStorage->K13[0] + S13_idx_0) - zcm;
  ktemp2[28] = b_ktemp2_tmp - rhoA;
  ktemp2[52] = elStorage->K33[0] + elStorage->S33[0] * c_tmp;
  b_ktemp2_tmp = elStorage->K34[0] + S34_idx_0;
  ktemp2[76] = b_ktemp2_tmp + H34_idx_0;
  f_ktemp2_tmp = elStorage->K35[0] + S35_idx_0;
  ktemp2[100] = f_ktemp2_tmp + H35_idx_0;
  g_ktemp2_tmp = elStorage->K36[0] + O1;
  ktemp2[124] = g_ktemp2_tmp + H36_idx_0;
  ktemp2[6] = ktemp2_tmp - H14_idx_0;
  ktemp2[30] = c_ktemp2_tmp - O1dot;
  ktemp2[54] = b_ktemp2_tmp - H34_idx_0;
  ktemp2[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * b_c_tmp +
    elStorage->S44_2[0] * c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[0] + O2;
  ktemp2[102] = ktemp2_tmp + H45_idx_0;
  b_ktemp2_tmp = elStorage->K46[0] + O3;
  ktemp2[126] = b_ktemp2_tmp + H46_idx_0;
  ktemp2[8] = K15_idx_0;
  ktemp2[32] = d_ktemp2_tmp - O2dot;
  ktemp2[56] = f_ktemp2_tmp - H35_idx_0;
  ktemp2[80] = ktemp2_tmp - H45_idx_0;
  ktemp2[104] = elStorage->K55[0] + elStorage->S55[0] * c_c_tmp;
  ktemp2[128] = integrationFactor;
  ktemp2[10] = ycm;
  ktemp2[34] = e_ktemp2_tmp - O3dot;
  ktemp2[58] = g_ktemp2_tmp - H36_idx_0;
  ktemp2[82] = b_ktemp2_tmp - H46_idx_0;
  ktemp2[106] = integrationFactor;
  ktemp2[130] = elStorage->K66[0] + elStorage->S66[0] * c_c_tmp;
  ktemp2[1] = elStorage->K11[1] + elStorage->S11[1] * c_c_tmp;
  ktemp2[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  ktemp2[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  ktemp2_tmp = elStorage->K14[1] + S14_idx_1;
  ktemp2[73] = ktemp2_tmp + H14_idx_1;
  ktemp2[97] = K15_idx_1;
  ktemp2[121] = K16_idx_1;
  ktemp2[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  ktemp2[27] = elStorage->K22[1] + elStorage->S22[1] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[1] + S23_idx_1;
  ktemp2[51] = b_ktemp2_tmp + H23_idx_1;
  c_ktemp2_tmp = elStorage->K24[1] + S24_idx_1;
  ktemp2[75] = c_ktemp2_tmp + H24_idx_1;
  d_ktemp2_tmp = elStorage->K25[1] + S25_idx_1;
  ktemp2[99] = d_ktemp2_tmp + H25_idx_1;
  e_ktemp2_tmp = elStorage->K26[1] + S26_idx_1;
  ktemp2[123] = e_ktemp2_tmp + H26_idx_1;
  ktemp2[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  f_ktemp2_tmp = elStorage->K23[2] + S23_idx_2;
  ktemp2[29] = f_ktemp2_tmp - H23_idx_2;
  ktemp2[53] = elStorage->K33[1] + elStorage->S33[1] * c_tmp;
  g_ktemp2_tmp = elStorage->K34[1] + S34_idx_1;
  ktemp2[77] = g_ktemp2_tmp + H34_idx_1;
  O1 = elStorage->K35[1] + S35_idx_1;
  ktemp2[101] = O1 + H35_idx_1;
  integrationFactor = elStorage->K36[1] + S36_idx_1;
  ktemp2[125] = integrationFactor + H36_idx_1;
  ycm = elStorage->K14[2] + S14_idx_2;
  ktemp2[7] = ycm - H14_idx_2;
  rhoA = elStorage->K24[2] + S24_idx_2;
  ktemp2[31] = rhoA - H24_idx_2;
  zcm = elStorage->K34[2] + S34_idx_2;
  ktemp2[55] = zcm - H34_idx_2;
  ktemp2[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * b_c_tmp +
    elStorage->S44_2[1] * c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  S34_idx_0 = elStorage->K45[1] + S45_idx_1;
  ktemp2[103] = S34_idx_0 + H45_idx_1;
  S35_idx_0 = elStorage->K46[1] + S46_idx_1;
  ktemp2[127] = S35_idx_0 + H46_idx_1;
  ktemp2[9] = K15_idx_2;
  S13_idx_0 = elStorage->K25[2] + S25_idx_2;
  ktemp2[33] = S13_idx_0 - H25_idx_2;
  S26_idx_0 = elStorage->K35[2] + S35_idx_2;
  ktemp2[57] = S26_idx_0 - H35_idx_2;
  S25_idx_0 = elStorage->K45[2] + S45_idx_2;
  ktemp2[81] = S25_idx_0 - H45_idx_2;
  ktemp2[105] = elStorage->K55[1] + elStorage->S55[1] * c_c_tmp;
  ktemp2[129] = K56_idx_1;
  ktemp2[11] = K16_idx_2;
  S24_idx_0 = elStorage->K26[2] + S26_idx_2;
  ktemp2[35] = S24_idx_0 - H26_idx_2;
  H12_idx_0 = elStorage->K36[2] + S36_idx_2;
  ktemp2[59] = H12_idx_0 - H36_idx_2;
  K15_idx_0 = elStorage->K46[2] + S46_idx_2;
  ktemp2[83] = K15_idx_0 - H46_idx_2;
  ktemp2[107] = K56_idx_2;
  ktemp2[131] = elStorage->K66[1] + elStorage->S66[1] * c_c_tmp;
  ktemp2[12] = elStorage->K11[2] + elStorage->S11[2] * c_c_tmp;
  ktemp2[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  ktemp2[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  ktemp2[84] = ycm + H14_idx_2;
  ktemp2[108] = K15_idx_2;
  ktemp2[132] = K16_idx_2;
  ktemp2[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  ktemp2[38] = elStorage->K22[2] + elStorage->S22[2] * b_c_tmp;
  ktemp2[62] = f_ktemp2_tmp + H23_idx_2;
  ktemp2[86] = rhoA + H24_idx_2;
  ktemp2[110] = S13_idx_0 + H25_idx_2;
  ktemp2[134] = S24_idx_0 + H26_idx_2;
  ktemp2[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  ktemp2[40] = b_ktemp2_tmp - H23_idx_1;
  ktemp2[64] = elStorage->K33[2] + elStorage->S33[2] * c_tmp;
  ktemp2[88] = zcm + H34_idx_2;
  ktemp2[112] = S26_idx_0 + H35_idx_2;
  ktemp2[136] = H12_idx_0 + H36_idx_2;
  ktemp2[18] = ktemp2_tmp - H14_idx_1;
  ktemp2[42] = c_ktemp2_tmp - H24_idx_1;
  ktemp2[66] = g_ktemp2_tmp - H34_idx_1;
  ktemp2[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * b_c_tmp +
    elStorage->S44_2[2] * c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  ktemp2[114] = S25_idx_0 + H45_idx_2;
  ktemp2[138] = K15_idx_0 + H46_idx_2;
  ktemp2[20] = K15_idx_1;
  ktemp2[44] = d_ktemp2_tmp - H25_idx_1;
  ktemp2[68] = O1 - H35_idx_1;
  ktemp2[92] = S34_idx_0 - H45_idx_1;
  ktemp2[116] = elStorage->K55[2] + elStorage->S55[2] * c_c_tmp;
  ktemp2[140] = K56_idx_2;
  ktemp2[22] = K16_idx_1;
  ktemp2[46] = e_ktemp2_tmp - H26_idx_1;
  ktemp2[70] = integrationFactor - H36_idx_1;
  ktemp2[94] = S35_idx_0 - H46_idx_1;
  ktemp2[118] = K56_idx_1;
  ktemp2[142] = elStorage->K66[2] + elStorage->S66[2] * c_c_tmp;
  ktemp2[13] = elStorage->K11[3] + elStorage->S11[3] * c_c_tmp;
  ktemp2[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  ktemp2[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  ktemp2_tmp = elStorage->K14[3] + S14_idx_3;
  ktemp2[85] = ktemp2_tmp + H14_idx_3;
  ktemp2[109] = K15_idx_3;
  ktemp2[133] = K16_idx_3;
  ktemp2[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  ktemp2[39] = elStorage->K22[3] + elStorage->S22[3] * b_c_tmp;
  b_ktemp2_tmp = elStorage->K23[3] + S23_idx_3;
  ktemp2[63] = b_ktemp2_tmp + H23_idx_3;
  c_ktemp2_tmp = elStorage->K24[3] + S24_idx_3;
  ktemp2[87] = c_ktemp2_tmp + H24_idx_3;
  d_ktemp2_tmp = elStorage->K25[3] + S25_idx_3;
  ktemp2[111] = d_ktemp2_tmp + H25_idx_3;
  e_ktemp2_tmp = elStorage->K26[3] + S26_idx_3;
  ktemp2[135] = e_ktemp2_tmp + H26_idx_3;
  ktemp2[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  ktemp2[41] = b_ktemp2_tmp - H23_idx_3;
  ktemp2[65] = elStorage->K33[3] + elStorage->S33[3] * c_tmp;
  b_ktemp2_tmp = elStorage->K34[3] + S34_idx_3;
  ktemp2[89] = b_ktemp2_tmp + H34_idx_3;
  f_ktemp2_tmp = elStorage->K35[3] + S35_idx_3;
  ktemp2[113] = f_ktemp2_tmp + H35_idx_3;
  g_ktemp2_tmp = elStorage->K36[3] + S36_idx_3;
  ktemp2[137] = g_ktemp2_tmp + H36_idx_3;
  ktemp2[19] = ktemp2_tmp - H14_idx_3;
  ktemp2[43] = c_ktemp2_tmp - H24_idx_3;
  ktemp2[67] = b_ktemp2_tmp - H34_idx_3;
  ktemp2[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * b_c_tmp +
    elStorage->S44_2[3] * c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  ktemp2_tmp = elStorage->K45[3] + S45_idx_3;
  ktemp2[115] = ktemp2_tmp + H45_idx_3;
  b_ktemp2_tmp = elStorage->K46[3] + S46_idx_3;
  ktemp2[139] = b_ktemp2_tmp + H46_idx_3;
  ktemp2[21] = K15_idx_3;
  ktemp2[45] = d_ktemp2_tmp - H25_idx_3;
  ktemp2[69] = f_ktemp2_tmp - H35_idx_3;
  ktemp2[93] = ktemp2_tmp - H45_idx_3;
  ktemp2[117] = elStorage->K55[3] + elStorage->S55[3] * c_c_tmp;
  ktemp2[141] = K56_idx_3;
  ktemp2[23] = K16_idx_3;
  ktemp2[47] = e_ktemp2_tmp - H26_idx_3;
  ktemp2[71] = g_ktemp2_tmp - H36_idx_3;
  ktemp2[95] = b_ktemp2_tmp - H46_idx_3;
  ktemp2[119] = K56_idx_3;
  ktemp2[143] = elStorage->K66[3] + elStorage->S66[3] * c_c_tmp;
  mapMatrixNonSym(ktemp2, Khate);

  //  Declare type
  // compile Coriolis/damping matrix
  //  ktemp = [zm,C12,C13,C14,zm,zm;
  //      -C12',zm,C23,C24,C25,C26;
  //      -C13',-C23',C33,C34,C35,C36;
  //      -C14',-C24',C43,C44,C45,C46;
  //      zm,-C25',-C35',-C45',zm,zm;
  //      zm,-C26',-C36',-C46',zm,zm];
  std::memset(&ktemp2[0], 0, 144U * sizeof(double));

  //  Row 1
  ktemp2_tmp = elStorage->C12[0] * Oel[2];
  ktemp2[2] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C12[2] * Oel[2];
  ktemp2[3] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C13[0] * Oel[1];
  ktemp2[4] = c_ktemp2_tmp;
  d_ktemp2_tmp = elStorage->C13[2] * Oel[1];
  ktemp2[5] = d_ktemp2_tmp;
  e_ktemp2_tmp = elStorage->C14_1[0] * Oel[1] + elStorage->C14_2[0] * Oel[2];
  ktemp2[6] = e_ktemp2_tmp;
  f_ktemp2_tmp = elStorage->C14_1[2] * Oel[1] + elStorage->C14_2[2] * Oel[2];
  ktemp2[7] = f_ktemp2_tmp;

  //  Row 2
  g_ktemp2_tmp = elStorage->C12[1] * Oel[2];
  ktemp2[14] = g_ktemp2_tmp;
  O1 = elStorage->C12[3] * Oel[2];
  ktemp2[15] = O1;
  integrationFactor = elStorage->C13[1] * Oel[1];
  ktemp2[16] = integrationFactor;
  ycm = elStorage->C13[3] * Oel[1];
  ktemp2[17] = ycm;
  rhoA = elStorage->C14_1[1] * Oel[1] + elStorage->C14_2[1] * Oel[2];
  ktemp2[18] = rhoA;
  zcm = elStorage->C14_1[3] * Oel[1] + elStorage->C14_2[3] * Oel[2];
  ktemp2[19] = zcm;

  //  Row 3
  ktemp2[24] = -ktemp2_tmp;
  ktemp2[25] = -g_ktemp2_tmp;
  ktemp2_tmp = elStorage->C23[0] * Oel[0];
  ktemp2[28] = ktemp2_tmp;
  g_ktemp2_tmp = elStorage->C23[2] * Oel[0];
  ktemp2[29] = g_ktemp2_tmp;
  S34_idx_0 = elStorage->C24[0] * Oel[0];
  ktemp2[30] = S34_idx_0;
  S35_idx_0 = elStorage->C24[2] * Oel[0];
  ktemp2[31] = S35_idx_0;
  S13_idx_0 = elStorage->C25[0] * Oel[2];
  ktemp2[32] = S13_idx_0;
  S26_idx_0 = elStorage->C25[2] * Oel[2];
  ktemp2[33] = S26_idx_0;
  S25_idx_0 = elStorage->C26[0] * Oel[2];
  ktemp2[34] = S25_idx_0;
  S24_idx_0 = elStorage->C26[2] * Oel[2];
  ktemp2[35] = S24_idx_0;

  //  Row 4
  ktemp2[36] = -b_ktemp2_tmp;
  ktemp2[37] = -O1;
  b_ktemp2_tmp = elStorage->C23[1] * Oel[0];
  ktemp2[40] = b_ktemp2_tmp;
  O1 = elStorage->C23[3] * Oel[0];
  ktemp2[41] = O1;
  H12_idx_0 = elStorage->C24[1] * Oel[0];
  ktemp2[42] = H12_idx_0;
  K15_idx_0 = elStorage->C24[3] * Oel[0];
  ktemp2[43] = K15_idx_0;
  H34_idx_3 = elStorage->C25[1] * Oel[2];
  ktemp2[44] = H34_idx_3;
  S14_idx_0 = elStorage->C25[3] * Oel[2];
  ktemp2[45] = S14_idx_0;
  S12_idx_0 = elStorage->C26[1] * Oel[2];
  ktemp2[46] = S12_idx_0;
  S23_idx_0 = elStorage->C26[3] * Oel[2];
  ktemp2[47] = S23_idx_0;

  //  Row 5
  ktemp2[48] = -c_ktemp2_tmp;
  ktemp2[49] = -integrationFactor;
  ktemp2[50] = -ktemp2_tmp;
  ktemp2[51] = -b_ktemp2_tmp;
  ktemp2[52] = 0.0;
  ktemp2[53] = 0.0;
  ktemp2[54] = C34_idx_0;
  ktemp2[55] = C34_idx_2;
  ktemp2_tmp = elStorage->C35[0] * Oel[1];
  ktemp2[56] = ktemp2_tmp;
  b_ktemp2_tmp = elStorage->C35[2] * Oel[1];
  ktemp2[57] = b_ktemp2_tmp;
  c_ktemp2_tmp = elStorage->C36[0] * Oel[1];
  ktemp2[58] = c_ktemp2_tmp;
  integrationFactor = elStorage->C36[2] * Oel[1];
  ktemp2[59] = integrationFactor;

  //  Row 6
  ktemp2[60] = -d_ktemp2_tmp;
  ktemp2[61] = -ycm;
  ktemp2[62] = -g_ktemp2_tmp;
  ktemp2[63] = -O1;
  ktemp2[64] = 0.0;
  ktemp2[65] = 0.0;
  ktemp2[66] = C34_idx_1;
  ktemp2[67] = C34_idx_3;
  d_ktemp2_tmp = elStorage->C35[1] * Oel[1];
  ktemp2[68] = d_ktemp2_tmp;
  g_ktemp2_tmp = elStorage->C35[3] * Oel[1];
  ktemp2[69] = g_ktemp2_tmp;
  O1 = elStorage->C36[1] * Oel[1];
  ktemp2[70] = O1;
  ycm = elStorage->C36[3] * Oel[1];
  ktemp2[71] = ycm;

  //  Row 7
  ktemp2[72] = -e_ktemp2_tmp;
  ktemp2[73] = -rhoA;
  ktemp2[74] = -S34_idx_0;
  ktemp2[75] = -H12_idx_0;
  ktemp2[76] = -C34_idx_0;
  ktemp2[77] = -C34_idx_1;
  ktemp2[78] = 0.0;
  ktemp2[79] = 0.0;
  e_ktemp2_tmp = elStorage->C45_1[0] * Oel[2] + elStorage->C45_2[0] * Oel[1];
  ktemp2[80] = e_ktemp2_tmp;
  rhoA = elStorage->C45_1[2] * Oel[2] + elStorage->C45_2[2] * Oel[1];
  ktemp2[81] = rhoA;
  S34_idx_0 = elStorage->C46_1[0] * Oel[1] + elStorage->C46_2[0] * Oel[2];
  ktemp2[82] = S34_idx_0;
  H12_idx_0 = elStorage->C46_1[2] * Oel[1] + elStorage->C46_2[2] * Oel[2];
  ktemp2[83] = H12_idx_0;

  //  Row 8
  ktemp2[84] = -f_ktemp2_tmp;
  ktemp2[85] = -zcm;
  ktemp2[86] = -S35_idx_0;
  ktemp2[87] = -K15_idx_0;
  ktemp2[88] = -C34_idx_2;
  ktemp2[89] = -C34_idx_3;
  ktemp2[90] = 0.0;
  ktemp2[91] = 0.0;
  f_ktemp2_tmp = elStorage->C45_1[1] * Oel[2] + elStorage->C45_2[1] * Oel[1];
  ktemp2[92] = f_ktemp2_tmp;
  zcm = elStorage->C45_1[3] * Oel[2] + elStorage->C45_2[3] * Oel[1];
  ktemp2[93] = zcm;
  S35_idx_0 = elStorage->C46_1[1] * Oel[1] + elStorage->C46_2[1] * Oel[2];
  ktemp2[94] = S35_idx_0;
  K15_idx_0 = elStorage->C46_1[3] * Oel[1] + elStorage->C46_2[3] * Oel[2];
  ktemp2[95] = K15_idx_0;

  //  Row 9
  ktemp2[98] = -S13_idx_0;
  ktemp2[99] = -H34_idx_3;
  ktemp2[100] = -ktemp2_tmp;
  ktemp2[101] = -d_ktemp2_tmp;
  ktemp2[102] = -e_ktemp2_tmp;
  ktemp2[103] = -f_ktemp2_tmp;

  //  Row 10
  ktemp2[110] = -S26_idx_0;
  ktemp2[111] = -S14_idx_0;
  ktemp2[112] = -b_ktemp2_tmp;
  ktemp2[113] = -g_ktemp2_tmp;
  ktemp2[114] = -rhoA;
  ktemp2[115] = -zcm;

  //  Row 11
  ktemp2[122] = -S25_idx_0;
  ktemp2[123] = -S12_idx_0;
  ktemp2[124] = -c_ktemp2_tmp;
  ktemp2[125] = -O1;
  ktemp2[126] = -S34_idx_0;
  ktemp2[127] = -S35_idx_0;

  //  Row 12
  ktemp2[134] = -S24_idx_0;
  ktemp2[135] = -S23_idx_0;
  ktemp2[136] = -integrationFactor;
  ktemp2[137] = -ycm;
  ktemp2[138] = -H12_idx_0;
  ktemp2[139] = -K15_idx_0;
  mapMatrixNonSym(ktemp2, Ce);

  // compile mass matrix
  ktemp2[0] = elStorage->M11[0];
  ktemp2[24] = 0.0;
  ktemp2[48] = 0.0;
  ktemp2[72] = 0.0;
  ktemp2[96] = elStorage->M15[0];
  ktemp2[120] = elStorage->M16[0];
  ktemp2[2] = 0.0;
  ktemp2[26] = elStorage->M22[0];
  ktemp2[50] = 0.0;
  ktemp2[74] = elStorage->M24[0];
  ktemp2[98] = 0.0;
  ktemp2[122] = 0.0;
  ktemp2[4] = 0.0;
  ktemp2[28] = 0.0;
  ktemp2[52] = elStorage->M33[0];
  ktemp2[76] = elStorage->M34[0];
  ktemp2[100] = 0.0;
  ktemp2[124] = 0.0;
  ktemp2[6] = 0.0;
  ktemp2[30] = elStorage->M24[0];
  ktemp2[54] = elStorage->M34[0];
  ktemp2[78] = elStorage->M44[0];
  ktemp2[102] = 0.0;
  ktemp2[126] = 0.0;
  ktemp2[8] = elStorage->M15[0];
  ktemp2[32] = 0.0;
  ktemp2[56] = 0.0;
  ktemp2[80] = 0.0;
  ktemp2[104] = elStorage->M55[0];
  ktemp2[128] = elStorage->M56[0];
  ktemp2[10] = elStorage->M16[0];
  ktemp2[34] = 0.0;
  ktemp2[58] = 0.0;
  ktemp2[82] = 0.0;
  ktemp2[106] = elStorage->M56[0];
  ktemp2[130] = elStorage->M66[0];
  ktemp2[1] = elStorage->M11[1];
  ktemp2[25] = 0.0;
  ktemp2[49] = 0.0;
  ktemp2[73] = 0.0;
  ktemp2[97] = elStorage->M15[1];
  ktemp2[121] = elStorage->M16[1];
  ktemp2[3] = 0.0;
  ktemp2[27] = elStorage->M22[1];
  ktemp2[51] = 0.0;
  ktemp2[75] = elStorage->M24[1];
  ktemp2[99] = 0.0;
  ktemp2[123] = 0.0;
  ktemp2[5] = 0.0;
  ktemp2[29] = 0.0;
  ktemp2[53] = elStorage->M33[1];
  ktemp2[77] = elStorage->M34[1];
  ktemp2[101] = 0.0;
  ktemp2[125] = 0.0;
  ktemp2[7] = 0.0;
  ktemp2[31] = elStorage->M24[2];
  ktemp2[55] = elStorage->M34[2];
  ktemp2[79] = elStorage->M44[1];
  ktemp2[103] = 0.0;
  ktemp2[127] = 0.0;
  ktemp2[9] = elStorage->M15[2];
  ktemp2[33] = 0.0;
  ktemp2[57] = 0.0;
  ktemp2[81] = 0.0;
  ktemp2[105] = elStorage->M55[1];
  ktemp2[129] = elStorage->M56[1];
  ktemp2[11] = elStorage->M16[2];
  ktemp2[35] = 0.0;
  ktemp2[59] = 0.0;
  ktemp2[83] = 0.0;
  ktemp2[107] = elStorage->M56[2];
  ktemp2[131] = elStorage->M66[1];
  ktemp2[12] = elStorage->M11[2];
  ktemp2[36] = 0.0;
  ktemp2[60] = 0.0;
  ktemp2[84] = 0.0;
  ktemp2[108] = elStorage->M15[2];
  ktemp2[132] = elStorage->M16[2];
  ktemp2[14] = 0.0;
  ktemp2[38] = elStorage->M22[2];
  ktemp2[62] = 0.0;
  ktemp2[86] = elStorage->M24[2];
  ktemp2[110] = 0.0;
  ktemp2[134] = 0.0;
  ktemp2[16] = 0.0;
  ktemp2[40] = 0.0;
  ktemp2[64] = elStorage->M33[2];
  ktemp2[88] = elStorage->M34[2];
  ktemp2[112] = 0.0;
  ktemp2[136] = 0.0;
  ktemp2[18] = 0.0;
  ktemp2[42] = elStorage->M24[1];
  ktemp2[66] = elStorage->M34[1];
  ktemp2[90] = elStorage->M44[2];
  ktemp2[114] = 0.0;
  ktemp2[138] = 0.0;
  ktemp2[20] = elStorage->M15[1];
  ktemp2[44] = 0.0;
  ktemp2[68] = 0.0;
  ktemp2[92] = 0.0;
  ktemp2[116] = elStorage->M55[2];
  ktemp2[140] = elStorage->M56[2];
  ktemp2[22] = elStorage->M16[1];
  ktemp2[46] = 0.0;
  ktemp2[70] = 0.0;
  ktemp2[94] = 0.0;
  ktemp2[118] = elStorage->M56[1];
  ktemp2[142] = elStorage->M66[2];
  ktemp2[13] = elStorage->M11[3];
  ktemp2[37] = 0.0;
  ktemp2[61] = 0.0;
  ktemp2[85] = 0.0;
  ktemp2[109] = elStorage->M15[3];
  ktemp2[133] = elStorage->M16[3];
  ktemp2[15] = 0.0;
  ktemp2[39] = elStorage->M22[3];
  ktemp2[63] = 0.0;
  ktemp2[87] = elStorage->M24[3];
  ktemp2[111] = 0.0;
  ktemp2[135] = 0.0;
  ktemp2[17] = 0.0;
  ktemp2[41] = 0.0;
  ktemp2[65] = elStorage->M33[3];
  ktemp2[89] = elStorage->M34[3];
  ktemp2[113] = 0.0;
  ktemp2[137] = 0.0;
  ktemp2[19] = 0.0;
  ktemp2[43] = elStorage->M24[3];
  ktemp2[67] = elStorage->M34[3];
  ktemp2[91] = elStorage->M44[3];
  ktemp2[115] = 0.0;
  ktemp2[139] = 0.0;
  ktemp2[21] = elStorage->M15[3];
  ktemp2[45] = 0.0;
  ktemp2[69] = 0.0;
  ktemp2[93] = 0.0;
  ktemp2[117] = elStorage->M55[3];
  ktemp2[141] = elStorage->M56[3];
  ktemp2[23] = elStorage->M16[3];
  ktemp2[47] = 0.0;
  ktemp2[71] = 0.0;
  ktemp2[95] = 0.0;
  ktemp2[119] = elStorage->M56[3];
  ktemp2[143] = elStorage->M66[3];
  mapMatrixNonSym(ktemp2, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input->RayleighAlpha * Kenr[i] + input->RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      ktemp2[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(ktemp2, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce, lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Khate,
                  lambda);
  c_sparse_mtimes(lambda, lambda_d, lambda_colidx, lambda_rowidx, Khate);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fhate[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fhate[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fhate, Fe);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input->concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input->concMass[0];
    Me[13] += input->concMass[0];
    Me[26] += input->concMass[0];
    Me[39] += input->concMass[1];
    Me[52] += input->concMass[2];
    Me[65] += input->concMass[3];
    Me[78] += input->concMass[4];
    Me[91] += input->concMass[4];
    Me[104] += input->concMass[4];
    Me[117] += input->concMass[5];
    Me[130] += input->concMass[6];
    Me[143] += input->concMass[7];

    // modify Ce for concentrated mass
    H34_idx_3 = 2.0 * input->concMass[0] * Omega;
    Ce[12] -= H34_idx_3;
    Ce[1] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[0] * 0.0;
    Ce[24] += H34_idx_3;
    Ce[2] -= H34_idx_3;
    Ce[25] -= H34_idx_3;
    Ce[14] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[4] * Omega;
    Ce[90] -= H34_idx_3;
    Ce[79] += H34_idx_3;
    H34_idx_3 = 2.0 * input->concMass[4] * 0.0;
    Ce[102] += H34_idx_3;
    Ce[80] -= H34_idx_3;
    Ce[103] -= H34_idx_3;
    Ce[92] += H34_idx_3;

    // modify Ke for concentrated mass
    H34_idx_3 = Omega * Omega;
    S23_idx_0 = input->concMass[0] * H34_idx_3;
    Khate[0] -= S23_idx_0;
    S14_idx_0 = input->concMass[0] * 0.0 * 0.0;
    S12_idx_0 = input->concMass[0] * OmegaDot;
    Khate[12] = (Khate[12] + S14_idx_0) - S12_idx_0;
    Khate[1] = (Khate[1] + S14_idx_0) + S12_idx_0;
    S14_idx_0 = input->concMass[0] * 0.0 * Omega;
    Khate[24] = (Khate[24] + S14_idx_0) + input->concMass[0] * 0.0;
    Khate[2] = (Khate[2] + S14_idx_0) - input->concMass[0] * 0.0;
    Khate[25] = (Khate[25] + S14_idx_0) - input->concMass[0] * 0.0;
    Khate[14] = (Khate[14] + S14_idx_0) + input->concMass[0] * 0.0;
    Khate[13] -= S23_idx_0;
    Khate[26] -= input->concMass[0] * 0.0;
    S23_idx_0 = input->concMass[4] * H34_idx_3;
    Khate[78] -= S23_idx_0;
    S14_idx_0 = input->concMass[4] * 0.0 * 0.0;
    S12_idx_0 = input->concMass[4] * OmegaDot;
    Khate[90] = (Khate[90] + S14_idx_0) - S12_idx_0;
    Khate[79] = (Khate[79] + S14_idx_0) + S12_idx_0;
    S14_idx_0 = input->concMass[4] * 0.0 * Omega;
    Khate[102] = (Khate[102] + S14_idx_0) + input->concMass[4] * 0.0;
    Khate[80] = (Khate[80] + S14_idx_0) - input->concMass[4] * 0.0;
    Khate[103] = (Khate[103] + S14_idx_0) - input->concMass[4] * 0.0;
    Khate[92] = (Khate[92] + S14_idx_0) + input->concMass[4] * 0.0;
    Khate[91] -= S23_idx_0;
    Khate[104] -= input->concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    H34_idx_3 = Omega * Omega;
    S23_idx_0 = 0.0 * Omega * input->z.data[0];
    Fe[0] = ((Fe[0] + input->concMass[0] * ((input->x.data[0] * H34_idx_3 - 0.0 *
                input->y.data[0]) - S23_idx_0)) + input->concMass[0] *
             (input->y.data[0] * OmegaDot - input->z.data[0] * 0.0)) -
      input->concMass[0] * a_temp[0];
    Fe[1] = ((Fe[1] + input->concMass[0] * ((input->y.data[0] * H34_idx_3 -
                S23_idx_0) - 0.0 * input->x.data[0])) + input->concMass[0] *
             (input->z.data[0] * 0.0 - input->x.data[0] * OmegaDot)) -
      input->concMass[0] * a_temp[1];
    Fe[2] = ((Fe[2] + input->concMass[0] * ((input->z.data[0] * 0.0 - Omega *
                0.0 * input->x.data[0]) - Omega * 0.0 * input->y.data[0])) +
             input->concMass[0] * (input->x.data[0] * 0.0 - input->y.data[0] *
              0.0)) - input->concMass[0] * a_temp[2];
    S23_idx_0 = 0.0 * Omega * input->z.data[1];
    Fe[6] = ((Fe[6] + input->concMass[4] * ((input->x.data[1] * H34_idx_3 - 0.0 *
                input->y.data[1]) - S23_idx_0)) + input->concMass[4] *
             (input->y.data[1] * OmegaDot - input->z.data[1] * 0.0)) -
      input->concMass[4] * a_temp[0];
    Fe[7] = ((Fe[7] + input->concMass[4] * ((input->y.data[1] * H34_idx_3 -
                S23_idx_0) - 0.0 * input->x.data[1])) + input->concMass[4] *
             (input->z.data[1] * 0.0 - input->x.data[1] * OmegaDot)) -
      input->concMass[4] * a_temp[1];
    Fe[8] = ((Fe[8] + input->concMass[4] * ((input->z.data[1] * 0.0 - Omega *
                0.0 * input->x.data[1]) - Omega * 0.0 * input->y.data[1])) +
             input->concMass[4] * (input->x.data[1] * 0.0 - input->y.data[1] *
              0.0)) - input->concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  std::memset(&Fhate[0], 0, 12U * sizeof(double));
  std::memset(&FhatLessConc[0], 0, 12U * sizeof(double));
  if (c_strcmp(input->analysisType)) {
    // calculate effective stiffness matrix and load vector for Newmark-Beta integrator 
    //      a1 = timeInt.a1;
    //      a2 = timeInt.a2;
    if (e_strcmp(input->iterationType)) {
      // considerations if newton raphson iteration is used
      if (input->firstIteration) {
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, FhatLessConc);
        std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
        mtimes(Khate, b_data, Ftemp_data);
        for (i = 0; i < 12; i++) {
          Fhate[i] = ((Fe[i] + Fhate[i]) + FhatLessConc[i]) - Ftemp_data[i];
        }
      } else {
        if (0 <= dispddot_size_idx_1 - 1) {
          std::memcpy(&b_data[0], &dispddot_data[0], dispddot_size_idx_1 *
                      sizeof(double));
        }

        if (dispddot_size_idx_1 == 1) {
          for (i = 0; i < 12; i++) {
            K15_idx_0 = 0.0;
            for (b_i = 0; b_i < 12; b_i++) {
              K15_idx_0 += Me[i + 12 * b_i] * b_data[b_i];
            }

            Fhate[i] = K15_idx_0;
          }
        } else {
          mtimes(Me, b_data, Fhate);
        }

        if (0 <= dispdot_size_idx_1 - 1) {
          std::memcpy(&b_data[0], &dispdot_data[0], dispdot_size_idx_1 * sizeof
                      (double));
        }

        if (dispdot_size_idx_1 == 1) {
          for (i = 0; i < 12; i++) {
            K15_idx_0 = 0.0;
            for (b_i = 0; b_i < 12; b_i++) {
              K15_idx_0 += Ce[i + 12 * b_i] * b_data[b_i];
            }

            FhatLessConc[i] = K15_idx_0;
          }
        } else {
          mtimes(Ce, b_data, FhatLessConc);
        }

        std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
        mtimes(Khate, b_data, Ftemp_data);
        for (i = 0; i < 12; i++) {
          Fhate[i] = ((Fe[i] - Fhate[i]) - FhatLessConc[i]) - Ftemp_data[i];
        }
      }
    } else {
      if (f_strcmp(input->iterationType)) {
        // considerations if direct iteration is used or linear analysis
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, FhatLessConc);
        for (i = 0; i < 12; i++) {
          Fhate[i] = (Fe[i] + Fhate[i]) + FhatLessConc[i];
        }
      }
    }

    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] + 1.0E+6 * Me[i]) + 1000.0 * Ce[i];
    }

    // ........................................................
    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));
    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Khate[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &Fe[0], 12U * sizeof(double));
  output->Me.size[0] = 1;
  output->Me.size[1] = 1;
  output->Me.data[0] = 0.0;
  output->Ce.size[0] = 1;
  output->Ce.size[1] = 1;
  output->Ce.data[0] = 0.0;
  if (d_strcmp(input->analysisType)) {
    output->Me.size[0] = 12;
    output->Me.size[1] = 12;
    output->Ce.size[0] = 12;
    output->Ce.size[1] = 12;
    std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
    std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));
  }

  if (c_strcmp(input->analysisType)) {
    std::memcpy(&output->FhatLessConc[0], &FhatLessConc[0], 12U * sizeof(double));
  }

  // ------------------------------------------
}

//
// File trailer for calculateTimoshenkoElementNL.cpp
//
// [EOF]
//
