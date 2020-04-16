//
// File: readElementData.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "readElementData.h"
#include "fileManager.h"
#include "getSplitLine.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <stdio.h>
#include <string.h>

// Variable Definitions
static const e_struct_T r = { { 0.0, 0.0 },// ac
  { 0.0, 0.0 },                        // twist
  { 0.0, 0.0 },                        // rhoA
  { 0.0, 0.0 },                        // EIyy
  { 0.0, 0.0 },                        // EIzz
  { 0.0, 0.0 },                        // GJ
  { 0.0, 0.0 },                        // EA
  { 0.0, 0.0 },                        // rhoIyy
  { 0.0, 0.0 },                        // rhoIzz
  { 0.0, 0.0 },                        // rhoJ
  { 0.0, 0.0 },                        // zcm
  { 0.0, 0.0 },                        // ycm
  { 0.0, 0.0 },                        // a
  { 0.0, 0.0 },                        // EIyz
  { 0.0, 0.0 },                        // alpha1
  { 0.0, 0.0 },                        // alpha2
  { 0.0, 0.0 },                        // alpha3
  { 0.0, 0.0 },                        // alpha4
  { 0.0, 0.0 },                        // alpha5
  { 0.0, 0.0 },                        // alpha6
  { 0.0, 0.0 },                        // rhoIyz
  { 0.0, 0.0 },                        // b
  { 0.0, 0.0 },                        // a0
  { 0.0, 0.0 }                         // aeroCenterOffset
};

// Function Definitions

//
// readElementData  reads element data
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [el] = readElementData(numElements,elfile,ortfile,bldfile
//
//    This function reads element data and stores data in the element data
//    object.
//
//       input:
//       numElements   = number of elements in structural mesh
//       elfile        = element data filename
//       ortfile       = element orientation filename
//       bldfile       = blade data filename
// Arguments    : double numElements
//                const emxArray_char_T *elfile
//                const emxArray_char_T *ortfile
//                const double bladeData_struct_nodeNum_data[]
//                const int bladeData_struct_nodeNum_size[1]
//                const double c_bladeData_struct_elementNum_d[]
//                const int c_bladeData_struct_elementNum_s[1]
//                const double bladeData_struct_remaining_data[]
//                const int bladeData_struct_remaining_size[2]
//                emxArray_struct_T *el_props
//                emxArray_real_T *el_elLen
//                emxArray_real_T *el_psi
//                emxArray_real_T *el_theta
//                emxArray_real_T *el_roll
//                emxArray_boolean_T *el_rotationalEffects
// Return Type  : void
//
void b_readElementData(double numElements, const emxArray_char_T *elfile, const
  emxArray_char_T *ortfile, const double bladeData_struct_nodeNum_data[], const
  int bladeData_struct_nodeNum_size[1], const double
  c_bladeData_struct_elementNum_d[], const int c_bladeData_struct_elementNum_s[1],
  const double bladeData_struct_remaining_data[], const int
  bladeData_struct_remaining_size[2], emxArray_struct_T *el_props,
  emxArray_real_T *el_elLen, emxArray_real_T *el_psi, emxArray_real_T *el_theta,
  emxArray_real_T *el_roll, emxArray_boolean_T *el_rotationalEffects)
{
  signed char fileid;
  int i;
  int loop_ub;
  int n;
  int idx;
  double ex;
  int k;
  boolean_T exitg1;
  double d;
  double y[2];
  boolean_T b_y;
  emxArray_real_T *temp;

  //       output:
  //       el            = element data object
  fileid = c_cfopen(elfile, "rb");

  // open element data file
  i = el_props->size[0] * el_props->size[1];
  el_props->size[0] = 1;
  loop_ub = static_cast<int>(numElements);
  el_props->size[1] = loop_ub;
  emxEnsureCapacity_struct_T(el_props, i);
  for (i = 0; i < loop_ub; i++) {
    el_props->data[i] = r;
  }

  for (n = 0; n < loop_ub; n++) {
    getSplitLine(static_cast<double>(fileid), el_elLen);

    // read element data
    getSplitLine(static_cast<double>(fileid), el_psi);

    // structural properties
    el_props->data[n].ac[0] = -(el_elLen->data[1] - 0.5);
    el_props->data[n].ac[1] = -(el_psi->data[1] - 0.5);
    el_props->data[n].twist[0] = el_elLen->data[2];
    el_props->data[n].twist[1] = el_psi->data[2];
    el_props->data[n].rhoA[0] = el_elLen->data[3];
    el_props->data[n].rhoA[1] = el_psi->data[3];
    el_props->data[n].EIyy[0] = el_elLen->data[4];
    el_props->data[n].EIyy[1] = el_psi->data[4];
    el_props->data[n].EIzz[0] = el_elLen->data[5];
    el_props->data[n].EIzz[1] = el_psi->data[5];
    y[0] = std::abs(el_props->data[n].EIyy[0] - el_props->data[n].EIzz[0]);
    y[1] = std::abs(el_props->data[n].EIyy[1] - el_props->data[n].EIzz[1]);
    b_y = true;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 2)) {
      if (!(y[k] < 0.001)) {
        b_y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (b_y) {
      el_props->data[n].EIzz[0] *= 1.0001;
      el_props->data[n].EIzz[1] *= 1.0001;
    }

    el_props->data[n].GJ[0] = el_elLen->data[6];
    el_props->data[n].GJ[1] = el_psi->data[6];
    el_props->data[n].EA[0] = el_elLen->data[7];
    el_props->data[n].EA[1] = el_psi->data[7];
    el_props->data[n].alpha1[0] = el_elLen->data[8];
    el_props->data[n].alpha1[1] = el_psi->data[8];
    el_props->data[n].rhoIyy[0] = el_elLen->data[9];
    el_props->data[n].rhoIyy[1] = el_psi->data[9];
    el_props->data[n].rhoIzz[0] = el_elLen->data[10];
    el_props->data[n].rhoIzz[1] = el_psi->data[10];
    el_props->data[n].rhoJ[0] = el_elLen->data[9] + el_elLen->data[10];
    el_props->data[n].rhoJ[1] = el_psi->data[9] + el_psi->data[10];
    el_props->data[n].zcm[0] = el_elLen->data[13];
    el_props->data[n].zcm[1] = el_psi->data[13];
    el_props->data[n].ycm[0] = el_elLen->data[14];
    el_props->data[n].ycm[1] = el_psi->data[14];
    el_props->data[n].a[0] = el_elLen->data[16];
    el_props->data[n].a[1] = el_psi->data[16];

    // coupling factors
    el_props->data[n].EIyz[0] = 0.0;
    el_props->data[n].alpha1[0] = 0.0;
    el_props->data[n].alpha2[0] = 0.0;
    el_props->data[n].alpha3[0] = 0.0;
    el_props->data[n].alpha4[0] = 0.0;
    el_props->data[n].alpha5[0] = 0.0;
    el_props->data[n].alpha6[0] = 0.0;
    el_props->data[n].rhoIyz[0] = 0.0;
    el_props->data[n].b[0] = 0.0;
    el_props->data[n].a0[0] = 6.2831853071795862;
    el_props->data[n].EIyz[1] = 0.0;
    el_props->data[n].alpha1[1] = 0.0;
    el_props->data[n].alpha2[1] = 0.0;
    el_props->data[n].alpha3[1] = 0.0;
    el_props->data[n].alpha4[1] = 0.0;
    el_props->data[n].alpha5[1] = 0.0;
    el_props->data[n].alpha6[1] = 0.0;
    el_props->data[n].rhoIyz[1] = 0.0;
    el_props->data[n].b[1] = 0.0;
    el_props->data[n].a0[1] = 6.2831853071795862;
  }

  // node number associated with blade section
  // element number associated with blade sectino
  // blade data
  n = bladeData_struct_nodeNum_size[0];
  if (bladeData_struct_nodeNum_size[0] <= 2) {
    if (bladeData_struct_nodeNum_size[0] == 1) {
      ex = bladeData_struct_nodeNum_data[0];
    } else if ((bladeData_struct_nodeNum_data[0] <
                bladeData_struct_nodeNum_data[1]) || (rtIsNaN
                (bladeData_struct_nodeNum_data[0]) && (!rtIsNaN
                 (bladeData_struct_nodeNum_data[1])))) {
      ex = bladeData_struct_nodeNum_data[1];
    } else {
      ex = bladeData_struct_nodeNum_data[0];
    }
  } else {
    if (!rtIsNaN(bladeData_struct_nodeNum_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= bladeData_struct_nodeNum_size[0])) {
        if (!rtIsNaN(bladeData_struct_nodeNum_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = bladeData_struct_nodeNum_data[0];
    } else {
      ex = bladeData_struct_nodeNum_data[idx - 1];
      i = idx + 1;
      for (k = i; k <= n; k++) {
        d = bladeData_struct_nodeNum_data[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }

  idx = static_cast<int>(ex);
  i = el_elLen->size[0];
  el_elLen->size[0] = idx;
  emxEnsureCapacity_real_T(el_elLen, i);
  for (i = 0; i < idx; i++) {
    el_elLen->data[i] = 0.0;
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (n = 0; n < i; n++) {
    el_elLen->data[static_cast<int>(bladeData_struct_nodeNum_data[n]) - 1] =
      bladeData_struct_remaining_data[n + bladeData_struct_remaining_size[0] * 9];

    // store chord of blade sections
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (n = 0; n < i; n++) {
    if (c_bladeData_struct_elementNum_d[n] != -1.0) {
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .b[0] = 0.5 * el_elLen->data[static_cast<int>
        (bladeData_struct_nodeNum_data[n]) - 1];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .b[1] = 0.5 * el_elLen->data[static_cast<int>
        (bladeData_struct_nodeNum_data[n + 1]) - 1];

      // element semi chord
      idx = n + bladeData_struct_remaining_size[0] * 11;
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        a0[0] = bladeData_struct_remaining_data[idx];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        a0[1] = bladeData_struct_remaining_data[idx + 1];

      // element lift curve slope (needed for flutter analysis)
      // convert "a" to semichord fraction aft of halfchord
      idx = static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1;

      // convert "ac" to semichord fraction foreward of halfchord
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .a[0] = ((el_props->data[idx].a[0] + 0.25 * (2.0 * el_props->data[idx]
        .b[0])) - el_props->data[idx].b[0]) / el_props->data[idx].b[0];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .a[1] = ((el_props->data[idx].a[1] + 0.25 * (2.0 * el_props->data[idx]
        .b[1])) - el_props->data[idx].b[1]) / el_props->data[idx].b[1];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        ac[0] = el_props->data[idx].ac[0] * 2.0;
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        ac[1] = el_props->data[idx].ac[1] * 2.0;

      // physical aero center offset from elastic axis
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        aeroCenterOffset[0] = el_props->data[idx].ac[0] * el_props->data[idx].b
        [0] - el_props->data[idx].a[0];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        aeroCenterOffset[1] = el_props->data[idx].ac[1] * el_props->data[idx].b
        [1] - el_props->data[idx].a[1];
    }
  }

  printf("%s\n", "EIyz, rhoIyz deactivated");
  fflush(stdout);
  cfclose(static_cast<double>(fileid));

  // close element file
  // read element orientation data
  fileid = c_cfopen(ortfile, "rb");
  i = el_elLen->size[0];
  el_elLen->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_elLen, i);
  i = el_psi->size[0];
  el_psi->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_psi, i);
  i = el_theta->size[0];
  el_theta->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_theta, i);
  i = el_roll->size[0];
  el_roll->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_roll, i);
  emxInit_real_T(&temp, 1);
  for (n = 0; n < loop_ub; n++) {
    getSplitLine(static_cast<double>(fileid), temp);
    el_elLen->data[n] = temp->data[4];
    el_psi->data[n] = temp->data[1];
    el_theta->data[n] = temp->data[2];
    el_roll->data[n] = temp->data[3];
  }

  emxFree_real_T(&temp);

  // store data in element object
  i = el_rotationalEffects->size[0] * el_rotationalEffects->size[1];
  el_rotationalEffects->size[0] = 1;
  el_rotationalEffects->size[1] = loop_ub;
  emxEnsureCapacity_boolean_T(el_rotationalEffects, i);
  for (i = 0; i < loop_ub; i++) {
    el_rotationalEffects->data[i] = true;
  }

  cfclose(static_cast<double>(fileid));

  // close ort file
}

//
// readElementData  reads element data
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [el] = readElementData(numElements,elfile,ortfile,bldfile
//
//    This function reads element data and stores data in the element data
//    object.
//
//       input:
//       numElements   = number of elements in structural mesh
//       elfile        = element data filename
//       ortfile       = element orientation filename
//       bldfile       = blade data filename
// Arguments    : double numElements
//                const emxArray_char_T *elfile
//                const emxArray_char_T *ortfile
//                const double bladeData_struct_nodeNum_data[]
//                const int bladeData_struct_nodeNum_size[1]
//                const double c_bladeData_struct_elementNum_d[]
//                const int c_bladeData_struct_elementNum_s[1]
//                const double bladeData_struct_remaining_data[]
//                const int bladeData_struct_remaining_size[2]
//                emxArray_struct_T *el_props
//                emxArray_real_T *el_elLen
//                emxArray_real_T *el_psi
//                emxArray_real_T *el_theta
//                emxArray_real_T *el_roll
//                emxArray_boolean_T *el_rotationalEffects
// Return Type  : void
//
void readElementData(double numElements, const emxArray_char_T *elfile, const
                     emxArray_char_T *ortfile, const double
                     bladeData_struct_nodeNum_data[], const int
                     bladeData_struct_nodeNum_size[1], const double
                     c_bladeData_struct_elementNum_d[], const int
                     c_bladeData_struct_elementNum_s[1], const double
                     bladeData_struct_remaining_data[], const int
                     bladeData_struct_remaining_size[2], emxArray_struct_T
                     *el_props, emxArray_real_T *el_elLen, emxArray_real_T
                     *el_psi, emxArray_real_T *el_theta, emxArray_real_T
                     *el_roll, emxArray_boolean_T *el_rotationalEffects)
{
  signed char fileid;
  int i;
  int loop_ub;
  int n;
  int idx;
  double ex;
  int k;
  boolean_T exitg1;
  double d;
  double y[2];
  boolean_T b_y;
  emxArray_real_T *temp;

  //       output:
  //       el            = element data object
  fileid = b_cfopen(elfile, "rb");

  // open element data file
  i = el_props->size[0] * el_props->size[1];
  el_props->size[0] = 1;
  loop_ub = static_cast<int>(numElements);
  el_props->size[1] = loop_ub;
  emxEnsureCapacity_struct_T(el_props, i);
  for (i = 0; i < loop_ub; i++) {
    el_props->data[i] = r;
  }

  for (n = 0; n < loop_ub; n++) {
    getSplitLine(static_cast<double>(fileid), el_elLen);

    // read element data
    getSplitLine(static_cast<double>(fileid), el_psi);

    // structural properties
    el_props->data[n].ac[0] = -(el_elLen->data[1] - 0.5);
    el_props->data[n].ac[1] = -(el_psi->data[1] - 0.5);
    el_props->data[n].twist[0] = el_elLen->data[2];
    el_props->data[n].twist[1] = el_psi->data[2];
    el_props->data[n].rhoA[0] = el_elLen->data[3];
    el_props->data[n].rhoA[1] = el_psi->data[3];
    el_props->data[n].EIyy[0] = el_elLen->data[4];
    el_props->data[n].EIyy[1] = el_psi->data[4];
    el_props->data[n].EIzz[0] = el_elLen->data[5];
    el_props->data[n].EIzz[1] = el_psi->data[5];
    y[0] = std::abs(el_props->data[n].EIyy[0] - el_props->data[n].EIzz[0]);
    y[1] = std::abs(el_props->data[n].EIyy[1] - el_props->data[n].EIzz[1]);
    b_y = true;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 2)) {
      if (!(y[k] < 0.001)) {
        b_y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (b_y) {
      el_props->data[n].EIzz[0] *= 1.0001;
      el_props->data[n].EIzz[1] *= 1.0001;
    }

    el_props->data[n].GJ[0] = el_elLen->data[6];
    el_props->data[n].GJ[1] = el_psi->data[6];
    el_props->data[n].EA[0] = el_elLen->data[7];
    el_props->data[n].EA[1] = el_psi->data[7];
    el_props->data[n].alpha1[0] = el_elLen->data[8];
    el_props->data[n].alpha1[1] = el_psi->data[8];
    el_props->data[n].rhoIyy[0] = el_elLen->data[9];
    el_props->data[n].rhoIyy[1] = el_psi->data[9];
    el_props->data[n].rhoIzz[0] = el_elLen->data[10];
    el_props->data[n].rhoIzz[1] = el_psi->data[10];
    el_props->data[n].rhoJ[0] = el_elLen->data[9] + el_elLen->data[10];
    el_props->data[n].rhoJ[1] = el_psi->data[9] + el_psi->data[10];
    el_props->data[n].zcm[0] = el_elLen->data[13];
    el_props->data[n].zcm[1] = el_psi->data[13];
    el_props->data[n].ycm[0] = el_elLen->data[14];
    el_props->data[n].ycm[1] = el_psi->data[14];
    el_props->data[n].a[0] = el_elLen->data[16];
    el_props->data[n].a[1] = el_psi->data[16];

    // coupling factors
    el_props->data[n].EIyz[0] = 0.0;
    el_props->data[n].alpha1[0] = 0.0;
    el_props->data[n].alpha2[0] = 0.0;
    el_props->data[n].alpha3[0] = 0.0;
    el_props->data[n].alpha4[0] = 0.0;
    el_props->data[n].alpha5[0] = 0.0;
    el_props->data[n].alpha6[0] = 0.0;
    el_props->data[n].rhoIyz[0] = 0.0;
    el_props->data[n].b[0] = 0.0;
    el_props->data[n].a0[0] = 6.2831853071795862;
    el_props->data[n].EIyz[1] = 0.0;
    el_props->data[n].alpha1[1] = 0.0;
    el_props->data[n].alpha2[1] = 0.0;
    el_props->data[n].alpha3[1] = 0.0;
    el_props->data[n].alpha4[1] = 0.0;
    el_props->data[n].alpha5[1] = 0.0;
    el_props->data[n].alpha6[1] = 0.0;
    el_props->data[n].rhoIyz[1] = 0.0;
    el_props->data[n].b[1] = 0.0;
    el_props->data[n].a0[1] = 6.2831853071795862;
  }

  // node number associated with blade section
  // element number associated with blade sectino
  // blade data
  n = bladeData_struct_nodeNum_size[0];
  if (bladeData_struct_nodeNum_size[0] <= 2) {
    if (bladeData_struct_nodeNum_size[0] == 1) {
      ex = bladeData_struct_nodeNum_data[0];
    } else if ((bladeData_struct_nodeNum_data[0] <
                bladeData_struct_nodeNum_data[1]) || (rtIsNaN
                (bladeData_struct_nodeNum_data[0]) && (!rtIsNaN
                 (bladeData_struct_nodeNum_data[1])))) {
      ex = bladeData_struct_nodeNum_data[1];
    } else {
      ex = bladeData_struct_nodeNum_data[0];
    }
  } else {
    if (!rtIsNaN(bladeData_struct_nodeNum_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= bladeData_struct_nodeNum_size[0])) {
        if (!rtIsNaN(bladeData_struct_nodeNum_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = bladeData_struct_nodeNum_data[0];
    } else {
      ex = bladeData_struct_nodeNum_data[idx - 1];
      i = idx + 1;
      for (k = i; k <= n; k++) {
        d = bladeData_struct_nodeNum_data[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }

  idx = static_cast<int>(ex);
  i = el_elLen->size[0];
  el_elLen->size[0] = idx;
  emxEnsureCapacity_real_T(el_elLen, i);
  for (i = 0; i < idx; i++) {
    el_elLen->data[i] = 0.0;
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (n = 0; n < i; n++) {
    el_elLen->data[static_cast<int>(bladeData_struct_nodeNum_data[n]) - 1] =
      bladeData_struct_remaining_data[n + bladeData_struct_remaining_size[0] * 9];

    // store chord of blade sections
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (n = 0; n < i; n++) {
    if (c_bladeData_struct_elementNum_d[n] != -1.0) {
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .b[0] = 0.5 * el_elLen->data[static_cast<int>
        (bladeData_struct_nodeNum_data[n]) - 1];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .b[1] = 0.5 * el_elLen->data[static_cast<int>
        (bladeData_struct_nodeNum_data[n + 1]) - 1];

      // element semi chord
      idx = n + bladeData_struct_remaining_size[0] * 11;
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        a0[0] = bladeData_struct_remaining_data[idx];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        a0[1] = bladeData_struct_remaining_data[idx + 1];

      // element lift curve slope (needed for flutter analysis)
      // convert "a" to semichord fraction aft of halfchord
      idx = static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1;

      // convert "ac" to semichord fraction foreward of halfchord
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .a[0] = ((el_props->data[idx].a[0] + 0.25 * (2.0 * el_props->data[idx]
        .b[0])) - el_props->data[idx].b[0]) / el_props->data[idx].b[0];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1]
        .a[1] = ((el_props->data[idx].a[1] + 0.25 * (2.0 * el_props->data[idx]
        .b[1])) - el_props->data[idx].b[1]) / el_props->data[idx].b[1];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        ac[0] = el_props->data[idx].ac[0] * 2.0;
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        ac[1] = el_props->data[idx].ac[1] * 2.0;

      // physical aero center offset from elastic axis
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        aeroCenterOffset[0] = el_props->data[idx].ac[0] * el_props->data[idx].b
        [0] - el_props->data[idx].a[0];
      el_props->data[static_cast<int>(c_bladeData_struct_elementNum_d[n]) - 1].
        aeroCenterOffset[1] = el_props->data[idx].ac[1] * el_props->data[idx].b
        [1] - el_props->data[idx].a[1];
    }
  }

  printf("%s\n", "EIyz, rhoIyz deactivated");
  fflush(stdout);
  cfclose(static_cast<double>(fileid));

  // close element file
  // read element orientation data
  fileid = b_cfopen(ortfile, "rb");
  i = el_elLen->size[0];
  el_elLen->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_elLen, i);
  i = el_psi->size[0];
  el_psi->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_psi, i);
  i = el_theta->size[0];
  el_theta->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_theta, i);
  i = el_roll->size[0];
  el_roll->size[0] = loop_ub;
  emxEnsureCapacity_real_T(el_roll, i);
  emxInit_real_T(&temp, 1);
  for (n = 0; n < loop_ub; n++) {
    getSplitLine(static_cast<double>(fileid), temp);
    el_elLen->data[n] = temp->data[4];
    el_psi->data[n] = temp->data[1];
    el_theta->data[n] = temp->data[2];
    el_roll->data[n] = temp->data[3];
  }

  emxFree_real_T(&temp);

  // store data in element object
  i = el_rotationalEffects->size[0] * el_rotationalEffects->size[1];
  el_rotationalEffects->size[0] = 1;
  el_rotationalEffects->size[1] = loop_ub;
  emxEnsureCapacity_boolean_T(el_rotationalEffects, i);
  for (i = 0; i < loop_ub; i++) {
    el_rotationalEffects->data[i] = true;
  }

  cfclose(static_cast<double>(fileid));

  // close ort file
}

//
// File trailer for readElementData.cpp
//
// [EOF]
//
