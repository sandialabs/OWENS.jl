//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateStrainForElements.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "calculateStrainForElements.h"
#include "calculateLambda.h"
#include "calculateShapeFunctions.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// calculate strains
// Arguments    : double numEl
//                const emxArray_real_T *conn
//                const emxArray_struct_T *el_props
//                const emxArray_real_T *el_elLen
//                const emxArray_real_T *el_psi
//                const emxArray_real_T *el_theta
//                const emxArray_real_T *el_roll
//                const emxArray_real_T *displ
//                b_emxArray_struct_T *elStrain
// Return Type  : void
//
void calculateStrainForElements(double numEl, const emxArray_real_T *conn, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const emxArray_real_T *displ, b_emxArray_struct_T
  *elStrain)
{
  int i;
  int loop_ub;
  int b_i;
  int b_index;
  double input_xloc[2];
  double input_disp_data[12];
  int j;
  int k;
  double a[144];
  double b_data[12];
  double dispLocal[12];
  int aoffset;
  double theta_yNode_idx_0;
  double theta_yNode_idx_1;
  double theta_zNode_idx_0;
  double theta_zNode_idx_1;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double unusedU1;
  double theta_x_prime;
  i = elStrain->size[0] * elStrain->size[1];
  elStrain->size[0] = 1;
  loop_ub = static_cast<int>(numEl);
  elStrain->size[1] = loop_ub;
  emxEnsureCapacity_struct_T1(elStrain, i);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    elStrain->data[b_i] = r1;

    // Calculate Ke and Fe for element i
    b_index = 0;
    input_xloc[0] = 0.0;
    input_xloc[1] = el_elLen->data[b_i];
    std::memset(&input_disp_data[0], 0, 12U * sizeof(double));
    for (j = 0; j < 2; j++) {
      // define element coordinates and displacements associated with element
      for (k = 0; k < 6; k++) {
        input_disp_data[b_index] = displ->data[static_cast<int>(((conn->data[b_i
          + conn->size[0] * j] - 1.0) * 6.0 + (static_cast<double>(k) + 1.0))) -
          1];
        b_index++;
      }
    }

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
    // -------- assign input block ----------------
    // --------------------------------------------
    // number of gauss points for full integration
    // calculate quad points
    // Initialize element sub matrices and sub vectors
    // Sort displacement vector
    // Written for 2 node element with 6 dof per node
    calculateLambda(el_psi->data[b_i] * 3.1415926535897931 / 180.0,
                    el_theta->data[b_i] * 3.1415926535897931 / 180.0,
                    (el_roll->data[b_i] + 0.5 * (el_props->data[b_i].twist[0] +
      el_props->data[b_i].twist[1])) * 3.1415926535897931 / 180.0, a);
    for (i = 0; i < 12; i++) {
      b_data[i] = input_disp_data[i];
      dispLocal[i] = 0.0;
    }

    for (k = 0; k < 12; k++) {
      aoffset = k * 12;
      for (i = 0; i < 12; i++) {
        dispLocal[i] += b_data[k] * a[aoffset + i];
      }
    }

    // '
    theta_yNode_idx_0 = dispLocal[4];
    theta_yNode_idx_1 = dispLocal[10];
    theta_zNode_idx_0 = dispLocal[5];
    theta_zNode_idx_1 = dispLocal[11];

    // Integration loop
    // END OF INTEGRATION LOOP
    for (i = 0; i < 4; i++) {
      // Calculate shape functions at quad point i
      calculateShapeFunctions(dv[i], input_xloc, N_data, N_size, p_N_x_data,
        p_N_x_size, &unusedU1);

      // N1 = N;
      // N2 = N;
      // N3 = N;
      // N4 = N;
      // calculate displacement derivatives at quad point i
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      // This function interpolates a value using distinct values at valNode
      // and the corresponding shape function N.
      theta_x_prime = p_N_x_data[0] * dispLocal[3] + p_N_x_data[1] * dispLocal[9];
      elStrain->data[b_i].eps_xx_0[i] = p_N_x_data[0] * dispLocal[0] +
        p_N_x_data[1] * dispLocal[6];
      elStrain->data[b_i].eps_xx_z[i] = p_N_x_data[0] * theta_yNode_idx_0 +
        p_N_x_data[1] * theta_yNode_idx_1;
      elStrain->data[b_i].eps_xx_y[i] = -(p_N_x_data[0] * theta_zNode_idx_0 +
        p_N_x_data[1] * theta_zNode_idx_1);
      elStrain->data[b_i].gam_xz_0[i] = (N_data[0] * theta_yNode_idx_0 + N_data
        [1] * theta_yNode_idx_1) + (p_N_x_data[0] * dispLocal[2] + p_N_x_data[1]
        * dispLocal[8]);
      elStrain->data[b_i].gam_xz_y[i] = theta_x_prime;
      elStrain->data[b_i].gam_xy_0[i] = -(N_data[0] * theta_zNode_idx_0 +
        N_data[1] * theta_zNode_idx_1) + (p_N_x_data[0] * dispLocal[1] +
        p_N_x_data[1] * dispLocal[7]);
      elStrain->data[b_i].gam_xy_z[i] = -theta_x_prime;
    }
  }
}

//
// File trailer for calculateStrainForElements.cpp
//
// [EOF]
//
