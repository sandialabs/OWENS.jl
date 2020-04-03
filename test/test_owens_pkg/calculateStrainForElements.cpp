//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateStrainForElements.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 03-Apr-2020 15:56:19
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
//                c_emxArray_struct_T *elStrain
// Return Type  : void
//
void calculateStrainForElements(double numEl, const emxArray_real_T *conn, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const emxArray_real_T *displ, c_emxArray_struct_T
  *elStrain)
{
  int b_index;
  int loop_ub;
  int i;
  static const f_struct_T b_r = { { 0.0, 0.0, 0.0, 0.0 },// eps_xx_0
    { 0.0, 0.0, 0.0, 0.0 },            // eps_xx_z
    { 0.0, 0.0, 0.0, 0.0 },            // eps_xx_y
    { 0.0, 0.0, 0.0, 0.0 },            // gam_xz_0
    { 0.0, 0.0, 0.0, 0.0 },            // gam_xz_y
    { 0.0, 0.0, 0.0, 0.0 },            // gam_xy_0
    { 0.0, 0.0, 0.0, 0.0 }             // gam_xy_z
  };

  double elInput_xloc[2];
  double elInput_disp_data[12];
  int aoffset;
  int k;
  double lambda[144];
  double b_data[12];
  double dispLocal[12];
  double theta_yNode_idx_0;
  double theta_yNode_idx_1;
  double theta_zNode_idx_0;
  double theta_zNode_idx_1;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double theta_x_prime;
  b_index = elStrain->size[0] * elStrain->size[1];
  elStrain->size[0] = 1;
  loop_ub = static_cast<int>(numEl);
  elStrain->size[1] = loop_ub;
  emxEnsureCapacity_struct_T5(elStrain, b_index);
  for (i = 0; i < loop_ub; i++) {
    elStrain->data[i] = b_r;

    // Calculate Ke and Fe for element i
    b_index = 0;
    elInput_xloc[0] = 0.0;
    elInput_xloc[1] = el_elLen->data[i];
    std::memset(&elInput_disp_data[0], 0, 12U * sizeof(double));
    for (aoffset = 0; aoffset < 2; aoffset++) {
      // define element coordinates and displacements associated with element
      for (k = 0; k < 6; k++) {
        elInput_disp_data[b_index] = displ->data[static_cast<int>(((conn->data[i
          + conn->size[0] * aoffset] - 1.0) * 6.0 + (static_cast<double>(k) +
          1.0))) - 1];
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
    calculateLambda(el_psi->data[i] * 3.1415926535897931 / 180.0, el_theta->
                    data[i] * 3.1415926535897931 / 180.0, (el_roll->data[i] +
      0.5 * (el_props->data[i].twist[0] + el_props->data[i].twist[1])) *
                    3.1415926535897931 / 180.0, lambda);
    for (b_index = 0; b_index < 12; b_index++) {
      b_data[b_index] = elInput_disp_data[b_index];
      dispLocal[b_index] = 0.0;
    }

    for (k = 0; k < 12; k++) {
      aoffset = k * 12;
      for (b_index = 0; b_index < 12; b_index++) {
        dispLocal[b_index] += b_data[k] * lambda[aoffset + b_index];
      }
    }

    // '
    theta_yNode_idx_0 = dispLocal[4];
    theta_yNode_idx_1 = dispLocal[10];
    theta_zNode_idx_0 = dispLocal[5];
    theta_zNode_idx_1 = dispLocal[11];

    // Integration loop
    // END OF INTEGRATION LOOP
    for (b_index = 0; b_index < 4; b_index++) {
      // Calculate shape functions at quad point i
      calculateShapeFunctions(dv[b_index], elInput_xloc, N_data, N_size,
        p_N_x_data, p_N_x_size, &theta_x_prime);

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
      elStrain->data[i].eps_xx_0[b_index] = p_N_x_data[0] * dispLocal[0] +
        p_N_x_data[1] * dispLocal[6];
      elStrain->data[i].eps_xx_z[b_index] = p_N_x_data[0] * theta_yNode_idx_0 +
        p_N_x_data[1] * theta_yNode_idx_1;
      elStrain->data[i].eps_xx_y[b_index] = -(p_N_x_data[0] * theta_zNode_idx_0
        + p_N_x_data[1] * theta_zNode_idx_1);
      elStrain->data[i].gam_xz_0[b_index] = (N_data[0] * theta_yNode_idx_0 +
        N_data[1] * theta_yNode_idx_1) + (p_N_x_data[0] * dispLocal[2] +
        p_N_x_data[1] * dispLocal[8]);
      elStrain->data[i].gam_xz_y[b_index] = theta_x_prime;
      elStrain->data[i].gam_xy_0[b_index] = -(N_data[0] * theta_zNode_idx_0 +
        N_data[1] * theta_zNode_idx_1) + (p_N_x_data[0] * dispLocal[1] +
        p_N_x_data[1] * dispLocal[7]);
      elStrain->data[i].gam_xy_z[b_index] = -theta_x_prime;
    }
  }
}

//
// File trailer for calculateStrainForElements.cpp
//
// [EOF]
//
