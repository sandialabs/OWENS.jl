//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateLoadVecFromDistForce.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "calculateLoadVecFromDistForce.h"
#include "calculateLambda.h"
#include "calculateShapeFunctions.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Definitions

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
//                double input_sweepAngle
//                double input_coneAngle
//                double input_rollAngle
//                const double input_extDistF2Node[2]
//                const double input_extDistF3Node[2]
//                const double input_extDistF4Node[2]
//                double output_Fe[12]
// Return Type  : void
//
void calculateLoadVecFromDistForce(const double input_xloc[2], const double
  input_sectionProps_twist[2], double input_sweepAngle, double input_coneAngle,
  double input_rollAngle, const double input_extDistF2Node[2], const double
  input_extDistF3Node[2], const double input_extDistF4Node[2], double output_Fe
  [12])
{
  double F1_idx_0;
  double F3_idx_0;
  double F2_idx_0;
  double F4_idx_0;
  double F5_idx_0;
  double F6_idx_0;
  double F1_idx_1;
  double F3_idx_1;
  double F2_idx_1;
  double F4_idx_1;
  double F5_idx_1;
  double F6_idx_1;
  int i;
  double N_data[2];
  int N_size[2];
  double unusedU0_data[2];
  int unusedU0_size[1];
  double Jac;
  double integrationFactor;
  double Ftemp[12];
  double bc;
  double extDistF2;
  double extDistF3;
  double b[12];
  double extDistF4;
  double F1_idx_1_tmp;
  double b_dv[144];
  emxArray_real_T *a_d;
  emxArray_int32_T *a_colidx;
  int b_i;
  emxArray_int32_T *a_rowidx;
  double b_dv1[144];
  double Fe[12];
  int acol;
  int apend;
  int nap;
  int ap;
  int Fe_tmp;

  // -------- assign input block ----------------
  // --------------------------------------------
  // number of gauss points for full integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_idx_0 = 0.0;
  F3_idx_0 = 0.0;
  F2_idx_0 = 0.0;
  F4_idx_0 = 0.0;
  F5_idx_0 = 0.0;
  F6_idx_0 = 0.0;
  F1_idx_1 = 0.0;
  F3_idx_1 = 0.0;
  F2_idx_1 = 0.0;
  F4_idx_1 = 0.0;
  F5_idx_1 = 0.0;
  F6_idx_1 = 0.0;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  // Integration loop
  for (i = 0; i < 4; i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[i], input_xloc, N_data, N_size, unusedU0_data,
      unusedU0_size, &Jac);
    integrationFactor = Jac * dv1[i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // distributed/body force load calculations
    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    bc = 0.0 * N_data[0] * integrationFactor;
    F1_idx_0 += bc;
    extDistF2 = N_data[0] * input_extDistF2Node[0] + N_data[1] *
      input_extDistF2Node[1];
    extDistF3 = N_data[0] * input_extDistF3Node[0] + N_data[1] *
      input_extDistF3Node[1];
    extDistF4 = N_data[0] * input_extDistF4Node[0] + N_data[1] *
      input_extDistF4Node[1];
    F1_idx_1_tmp = 0.0 * N_data[1] * integrationFactor;
    F1_idx_1 += F1_idx_1_tmp;

    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element vector
    // Element calculation functions---------------------------------
    F2_idx_0 += extDistF2 * N_data[0] * integrationFactor;
    F3_idx_0 += extDistF3 * N_data[0] * integrationFactor;
    F4_idx_0 += extDistF4 * N_data[0] * integrationFactor;
    F5_idx_0 += bc;
    F6_idx_0 += bc;
    F2_idx_1 += extDistF2 * N_data[1] * integrationFactor;
    F3_idx_1 += extDistF3 * N_data[1] * integrationFactor;
    F4_idx_1 += extDistF4 * N_data[1] * integrationFactor;
    F5_idx_1 += F1_idx_1_tmp;
    F6_idx_1 += F1_idx_1_tmp;
  }

  // END OF INTEGRATION LOOP
  // ---------------------------------------------
  // compile element force vector
  Ftemp[0] = F1_idx_0;
  Ftemp[2] = F2_idx_0;
  Ftemp[4] = F3_idx_0;
  Ftemp[6] = F4_idx_0;
  Ftemp[8] = F5_idx_0;
  Ftemp[10] = F6_idx_0;
  Ftemp[1] = F1_idx_1;
  Ftemp[3] = F2_idx_1;
  Ftemp[5] = F3_idx_1;
  Ftemp[7] = F4_idx_1;
  Ftemp[9] = F5_idx_1;
  Ftemp[11] = F6_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  std::memset(&b[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  //  %------------------------------------------------------------------------- 
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  calculateLambda(input_sweepAngle * 3.1415926535897931 / 180.0, input_coneAngle
                  * 3.1415926535897931 / 180.0, (input_rollAngle + 0.5 *
    (input_sectionProps_twist[0] + input_sectionProps_twist[1])) *
                  3.1415926535897931 / 180.0, b_dv);
  for (i = 0; i < 12; i++) {
    b[iv1[i] - 1] = Ftemp[i];
    for (b_i = 0; b_i < 12; b_i++) {
      b_dv1[b_i + 12 * i] = b_dv[i + 12 * b_i];
    }

    Fe[i] = 0.0;
  }

  emxInit_real_T(&a_d, 1);
  emxInit_int32_T(&a_colidx, 1);
  emxInit_int32_T(&a_rowidx, 1);
  sparse(b_dv1, a_d, a_colidx, a_rowidx);
  if (a_colidx->data[a_colidx->size[0] - 1] - 1 != 0) {
    for (acol = 0; acol < 12; acol++) {
      bc = b[acol];
      b_i = a_colidx->data[acol];
      apend = a_colidx->data[acol + 1];
      nap = apend - a_colidx->data[acol];
      if (nap >= 4) {
        apend = (apend - nap) + ((nap / 4) << 2);
        for (ap = b_i; ap <= apend - 1; ap += 4) {
          Fe_tmp = a_rowidx->data[ap - 1] - 1;
          Fe[Fe_tmp] += a_d->data[ap - 1] * bc;
          Fe[a_rowidx->data[ap] - 1] += a_d->data[ap] * bc;
          Fe_tmp = a_rowidx->data[ap + 1] - 1;
          Fe[Fe_tmp] += a_d->data[ap + 1] * bc;
          Fe_tmp = a_rowidx->data[ap + 2] - 1;
          Fe[Fe_tmp] += a_d->data[ap + 2] * bc;
        }

        nap = a_colidx->data[acol + 1] - 1;
        for (ap = apend; ap <= nap; ap++) {
          Fe_tmp = a_rowidx->data[ap - 1] - 1;
          Fe[Fe_tmp] += a_d->data[ap - 1] * bc;
        }
      } else {
        apend--;
        for (ap = b_i; ap <= apend; ap++) {
          Fe_tmp = a_rowidx->data[ap - 1] - 1;
          Fe[Fe_tmp] += a_d->data[ap - 1] * bc;
        }
      }
    }
  }

  emxFree_int32_T(&a_rowidx);
  emxFree_int32_T(&a_colidx);
  emxFree_real_T(&a_d);

  // ----- assign output block ----------------
  std::memcpy(&output_Fe[0], &Fe[0], 12U * sizeof(double));

  // ------------------------------------------
}

//
// File trailer for calculateLoadVecFromDistForce.cpp
//
// [EOF]
//
