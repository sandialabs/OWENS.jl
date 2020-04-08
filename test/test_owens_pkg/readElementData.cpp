//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readElementData.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
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
static const d_struct_T r = { { 0.0, 0.0 },// ac
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
  emxArray_struct_T *sectionPropsArray;
  signed char fileid;
  int i;
  int loop_ub;
  emxArray_real_T *data1;
  emxArray_real_T *data2;
  int b_i;
  int n;
  int idx;
  double ex;
  int k;
  boolean_T exitg1;
  emxArray_real_T *chord;
  double d;
  double y[2];
  boolean_T b_y;
  emxArray_real_T *elLen;
  emxArray_real_T *psi;
  emxArray_real_T *theta;
  emxArray_real_T *roll;
  emxInit_struct_T(&sectionPropsArray, 2);

  //       output:
  //       el            = element data object
  fileid = c_cfopen(elfile, "rb");

  // open element data file
  i = sectionPropsArray->size[0] * sectionPropsArray->size[1];
  sectionPropsArray->size[0] = 1;
  loop_ub = static_cast<int>(numElements);
  sectionPropsArray->size[1] = loop_ub;
  emxEnsureCapacity_struct_T(sectionPropsArray, i);
  for (i = 0; i < loop_ub; i++) {
    sectionPropsArray->data[i] = r;
  }

  emxInit_real_T(&data1, 1);
  emxInit_real_T(&data2, 1);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    getSplitLine(static_cast<double>(fileid), data1);

    // read element data
    getSplitLine(static_cast<double>(fileid), data2);

    // structural properties
    sectionPropsArray->data[b_i].ac[0] = -(data1->data[1] - 0.5);
    sectionPropsArray->data[b_i].ac[1] = -(data2->data[1] - 0.5);
    sectionPropsArray->data[b_i].twist[0] = data1->data[2];
    sectionPropsArray->data[b_i].twist[1] = data2->data[2];
    sectionPropsArray->data[b_i].rhoA[0] = data1->data[3];
    sectionPropsArray->data[b_i].rhoA[1] = data2->data[3];
    sectionPropsArray->data[b_i].EIyy[0] = data1->data[4];
    sectionPropsArray->data[b_i].EIyy[1] = data2->data[4];
    sectionPropsArray->data[b_i].EIzz[0] = data1->data[5];
    sectionPropsArray->data[b_i].EIzz[1] = data2->data[5];
    y[0] = std::abs(sectionPropsArray->data[b_i].EIyy[0] -
                    sectionPropsArray->data[b_i].EIzz[0]);
    y[1] = std::abs(sectionPropsArray->data[b_i].EIyy[1] -
                    sectionPropsArray->data[b_i].EIzz[1]);
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
      sectionPropsArray->data[b_i].EIzz[0] *= 1.0001;
      sectionPropsArray->data[b_i].EIzz[1] *= 1.0001;
    }

    sectionPropsArray->data[b_i].GJ[0] = data1->data[6];
    sectionPropsArray->data[b_i].GJ[1] = data2->data[6];
    sectionPropsArray->data[b_i].EA[0] = data1->data[7];
    sectionPropsArray->data[b_i].EA[1] = data2->data[7];
    sectionPropsArray->data[b_i].alpha1[0] = data1->data[8];
    sectionPropsArray->data[b_i].alpha1[1] = data2->data[8];
    sectionPropsArray->data[b_i].rhoIyy[0] = data1->data[9];
    sectionPropsArray->data[b_i].rhoIyy[1] = data2->data[9];
    sectionPropsArray->data[b_i].rhoIzz[0] = data1->data[10];
    sectionPropsArray->data[b_i].rhoIzz[1] = data2->data[10];
    sectionPropsArray->data[b_i].rhoJ[0] = data1->data[9] + data1->data[10];
    sectionPropsArray->data[b_i].rhoJ[1] = data2->data[9] + data2->data[10];
    sectionPropsArray->data[b_i].zcm[0] = data1->data[13];
    sectionPropsArray->data[b_i].zcm[1] = data2->data[13];
    sectionPropsArray->data[b_i].ycm[0] = data1->data[14];
    sectionPropsArray->data[b_i].ycm[1] = data2->data[14];
    sectionPropsArray->data[b_i].a[0] = data1->data[16];
    sectionPropsArray->data[b_i].a[1] = data2->data[16];

    // coupling factors
    sectionPropsArray->data[b_i].EIyz[0] = 0.0;
    sectionPropsArray->data[b_i].alpha1[0] = 0.0;
    sectionPropsArray->data[b_i].alpha2[0] = 0.0;
    sectionPropsArray->data[b_i].alpha3[0] = 0.0;
    sectionPropsArray->data[b_i].alpha4[0] = 0.0;
    sectionPropsArray->data[b_i].alpha5[0] = 0.0;
    sectionPropsArray->data[b_i].alpha6[0] = 0.0;
    sectionPropsArray->data[b_i].rhoIyz[0] = 0.0;
    sectionPropsArray->data[b_i].b[0] = 0.0;
    sectionPropsArray->data[b_i].a0[0] = 6.2831853071795862;
    sectionPropsArray->data[b_i].EIyz[1] = 0.0;
    sectionPropsArray->data[b_i].alpha1[1] = 0.0;
    sectionPropsArray->data[b_i].alpha2[1] = 0.0;
    sectionPropsArray->data[b_i].alpha3[1] = 0.0;
    sectionPropsArray->data[b_i].alpha4[1] = 0.0;
    sectionPropsArray->data[b_i].alpha5[1] = 0.0;
    sectionPropsArray->data[b_i].alpha6[1] = 0.0;
    sectionPropsArray->data[b_i].rhoIyz[1] = 0.0;
    sectionPropsArray->data[b_i].b[1] = 0.0;
    sectionPropsArray->data[b_i].a0[1] = 6.2831853071795862;
  }

  emxFree_real_T(&data2);

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

  emxInit_real_T(&chord, 1);
  idx = static_cast<int>(ex);
  i = chord->size[0];
  chord->size[0] = idx;
  emxEnsureCapacity_real_T(chord, i);
  for (i = 0; i < idx; i++) {
    chord->data[i] = 0.0;
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (b_i = 0; b_i < i; b_i++) {
    chord->data[static_cast<int>(bladeData_struct_nodeNum_data[b_i]) - 1] =
      bladeData_struct_remaining_data[b_i + bladeData_struct_remaining_size[0] *
      9];

    // store chord of blade sections
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (b_i = 0; b_i < i; b_i++) {
    if (c_bladeData_struct_elementNum_d[b_i] != -1.0) {
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].b[0] = 0.5 * chord->data[
        static_cast<int>(bladeData_struct_nodeNum_data[b_i]) - 1];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].b[1] = 0.5 * chord->data[
        static_cast<int>(bladeData_struct_nodeNum_data[b_i + 1]) - 1];

      // element semi chord
      idx = b_i + bladeData_struct_remaining_size[0] * 11;
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a0[0] =
        bladeData_struct_remaining_data[idx];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a0[1] =
        bladeData_struct_remaining_data[idx + 1];

      // element lift curve slope (needed for flutter analysis)
      // convert "a" to semichord fraction aft of halfchord
      idx = static_cast<int>(c_bladeData_struct_elementNum_d[b_i]) - 1;

      // convert "ac" to semichord fraction foreward of halfchord
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a[0] =
        ((sectionPropsArray->data[idx].a[0] + 0.25 * (2.0 *
           sectionPropsArray->data[idx].b[0])) - sectionPropsArray->data[idx].b
         [0]) / sectionPropsArray->data[idx].b[0];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a[1] =
        ((sectionPropsArray->data[idx].a[1] + 0.25 * (2.0 *
           sectionPropsArray->data[idx].b[1])) - sectionPropsArray->data[idx].b
         [1]) / sectionPropsArray->data[idx].b[1];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].ac[0] =
        sectionPropsArray->data[idx].ac[0] * 2.0;
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].ac[1] =
        sectionPropsArray->data[idx].ac[1] * 2.0;

      // physical aero center offset from elastic axis
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].aeroCenterOffset[0] =
        sectionPropsArray->data[idx].ac[0] * sectionPropsArray->data[idx].b[0] -
        sectionPropsArray->data[idx].a[0];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].aeroCenterOffset[1] =
        sectionPropsArray->data[idx].ac[1] * sectionPropsArray->data[idx].b[1] -
        sectionPropsArray->data[idx].a[1];
    }
  }

  emxFree_real_T(&chord);
  emxInit_real_T(&elLen, 1);
  emxInit_real_T(&psi, 1);
  emxInit_real_T(&theta, 1);
  emxInit_real_T(&roll, 1);
  printf("%s\n", "EIyz, rhoIyz deactivated");
  fflush(stdout);
  cfclose(static_cast<double>(fileid));

  // close element file
  // read element orientation data
  fileid = c_cfopen(ortfile, "rb");
  i = elLen->size[0];
  elLen->size[0] = loop_ub;
  emxEnsureCapacity_real_T(elLen, i);
  i = psi->size[0];
  psi->size[0] = loop_ub;
  emxEnsureCapacity_real_T(psi, i);
  i = theta->size[0];
  theta->size[0] = loop_ub;
  emxEnsureCapacity_real_T(theta, i);
  i = roll->size[0];
  roll->size[0] = loop_ub;
  emxEnsureCapacity_real_T(roll, i);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    getSplitLine(static_cast<double>(fileid), data1);
    elLen->data[b_i] = data1->data[4];
    psi->data[b_i] = data1->data[1];
    theta->data[b_i] = data1->data[2];
    roll->data[b_i] = data1->data[3];
  }

  emxFree_real_T(&data1);

  // store data in element object
  i = el_props->size[0] * el_props->size[1];
  el_props->size[0] = 1;
  el_props->size[1] = sectionPropsArray->size[1];
  emxEnsureCapacity_struct_T(el_props, i);
  idx = sectionPropsArray->size[0] * sectionPropsArray->size[1];
  for (i = 0; i < idx; i++) {
    el_props->data[i] = sectionPropsArray->data[i];
  }

  emxFree_struct_T(&sectionPropsArray);
  i = el_elLen->size[0];
  el_elLen->size[0] = elLen->size[0];
  emxEnsureCapacity_real_T(el_elLen, i);
  idx = elLen->size[0];
  for (i = 0; i < idx; i++) {
    el_elLen->data[i] = elLen->data[i];
  }

  emxFree_real_T(&elLen);
  i = el_psi->size[0];
  el_psi->size[0] = psi->size[0];
  emxEnsureCapacity_real_T(el_psi, i);
  idx = psi->size[0];
  for (i = 0; i < idx; i++) {
    el_psi->data[i] = psi->data[i];
  }

  emxFree_real_T(&psi);
  i = el_theta->size[0];
  el_theta->size[0] = theta->size[0];
  emxEnsureCapacity_real_T(el_theta, i);
  idx = theta->size[0];
  for (i = 0; i < idx; i++) {
    el_theta->data[i] = theta->data[i];
  }

  emxFree_real_T(&theta);
  i = el_roll->size[0];
  el_roll->size[0] = roll->size[0];
  emxEnsureCapacity_real_T(el_roll, i);
  idx = roll->size[0];
  for (i = 0; i < idx; i++) {
    el_roll->data[i] = roll->data[i];
  }

  emxFree_real_T(&roll);
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
  emxArray_struct_T *sectionPropsArray;
  signed char fileid;
  int i;
  int loop_ub;
  emxArray_real_T *data1;
  emxArray_real_T *data2;
  int b_i;
  int n;
  int idx;
  double ex;
  int k;
  boolean_T exitg1;
  emxArray_real_T *chord;
  double d;
  double y[2];
  boolean_T b_y;
  emxArray_real_T *elLen;
  emxArray_real_T *psi;
  emxArray_real_T *theta;
  emxArray_real_T *roll;
  emxInit_struct_T(&sectionPropsArray, 2);

  //       output:
  //       el            = element data object
  fileid = b_cfopen(elfile, "rb");

  // open element data file
  i = sectionPropsArray->size[0] * sectionPropsArray->size[1];
  sectionPropsArray->size[0] = 1;
  loop_ub = static_cast<int>(numElements);
  sectionPropsArray->size[1] = loop_ub;
  emxEnsureCapacity_struct_T(sectionPropsArray, i);
  for (i = 0; i < loop_ub; i++) {
    sectionPropsArray->data[i] = r;
  }

  emxInit_real_T(&data1, 1);
  emxInit_real_T(&data2, 1);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    getSplitLine(static_cast<double>(fileid), data1);

    // read element data
    getSplitLine(static_cast<double>(fileid), data2);

    // structural properties
    sectionPropsArray->data[b_i].ac[0] = -(data1->data[1] - 0.5);
    sectionPropsArray->data[b_i].ac[1] = -(data2->data[1] - 0.5);
    sectionPropsArray->data[b_i].twist[0] = data1->data[2];
    sectionPropsArray->data[b_i].twist[1] = data2->data[2];
    sectionPropsArray->data[b_i].rhoA[0] = data1->data[3];
    sectionPropsArray->data[b_i].rhoA[1] = data2->data[3];
    sectionPropsArray->data[b_i].EIyy[0] = data1->data[4];
    sectionPropsArray->data[b_i].EIyy[1] = data2->data[4];
    sectionPropsArray->data[b_i].EIzz[0] = data1->data[5];
    sectionPropsArray->data[b_i].EIzz[1] = data2->data[5];
    y[0] = std::abs(sectionPropsArray->data[b_i].EIyy[0] -
                    sectionPropsArray->data[b_i].EIzz[0]);
    y[1] = std::abs(sectionPropsArray->data[b_i].EIyy[1] -
                    sectionPropsArray->data[b_i].EIzz[1]);
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
      sectionPropsArray->data[b_i].EIzz[0] *= 1.0001;
      sectionPropsArray->data[b_i].EIzz[1] *= 1.0001;
    }

    sectionPropsArray->data[b_i].GJ[0] = data1->data[6];
    sectionPropsArray->data[b_i].GJ[1] = data2->data[6];
    sectionPropsArray->data[b_i].EA[0] = data1->data[7];
    sectionPropsArray->data[b_i].EA[1] = data2->data[7];
    sectionPropsArray->data[b_i].alpha1[0] = data1->data[8];
    sectionPropsArray->data[b_i].alpha1[1] = data2->data[8];
    sectionPropsArray->data[b_i].rhoIyy[0] = data1->data[9];
    sectionPropsArray->data[b_i].rhoIyy[1] = data2->data[9];
    sectionPropsArray->data[b_i].rhoIzz[0] = data1->data[10];
    sectionPropsArray->data[b_i].rhoIzz[1] = data2->data[10];
    sectionPropsArray->data[b_i].rhoJ[0] = data1->data[9] + data1->data[10];
    sectionPropsArray->data[b_i].rhoJ[1] = data2->data[9] + data2->data[10];
    sectionPropsArray->data[b_i].zcm[0] = data1->data[13];
    sectionPropsArray->data[b_i].zcm[1] = data2->data[13];
    sectionPropsArray->data[b_i].ycm[0] = data1->data[14];
    sectionPropsArray->data[b_i].ycm[1] = data2->data[14];
    sectionPropsArray->data[b_i].a[0] = data1->data[16];
    sectionPropsArray->data[b_i].a[1] = data2->data[16];

    // coupling factors
    sectionPropsArray->data[b_i].EIyz[0] = 0.0;
    sectionPropsArray->data[b_i].alpha1[0] = 0.0;
    sectionPropsArray->data[b_i].alpha2[0] = 0.0;
    sectionPropsArray->data[b_i].alpha3[0] = 0.0;
    sectionPropsArray->data[b_i].alpha4[0] = 0.0;
    sectionPropsArray->data[b_i].alpha5[0] = 0.0;
    sectionPropsArray->data[b_i].alpha6[0] = 0.0;
    sectionPropsArray->data[b_i].rhoIyz[0] = 0.0;
    sectionPropsArray->data[b_i].b[0] = 0.0;
    sectionPropsArray->data[b_i].a0[0] = 6.2831853071795862;
    sectionPropsArray->data[b_i].EIyz[1] = 0.0;
    sectionPropsArray->data[b_i].alpha1[1] = 0.0;
    sectionPropsArray->data[b_i].alpha2[1] = 0.0;
    sectionPropsArray->data[b_i].alpha3[1] = 0.0;
    sectionPropsArray->data[b_i].alpha4[1] = 0.0;
    sectionPropsArray->data[b_i].alpha5[1] = 0.0;
    sectionPropsArray->data[b_i].alpha6[1] = 0.0;
    sectionPropsArray->data[b_i].rhoIyz[1] = 0.0;
    sectionPropsArray->data[b_i].b[1] = 0.0;
    sectionPropsArray->data[b_i].a0[1] = 6.2831853071795862;
  }

  emxFree_real_T(&data2);

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

  emxInit_real_T(&chord, 1);
  idx = static_cast<int>(ex);
  i = chord->size[0];
  chord->size[0] = idx;
  emxEnsureCapacity_real_T(chord, i);
  for (i = 0; i < idx; i++) {
    chord->data[i] = 0.0;
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (b_i = 0; b_i < i; b_i++) {
    chord->data[static_cast<int>(bladeData_struct_nodeNum_data[b_i]) - 1] =
      bladeData_struct_remaining_data[b_i + bladeData_struct_remaining_size[0] *
      9];

    // store chord of blade sections
  }

  i = c_bladeData_struct_elementNum_s[0];
  for (b_i = 0; b_i < i; b_i++) {
    if (c_bladeData_struct_elementNum_d[b_i] != -1.0) {
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].b[0] = 0.5 * chord->data[
        static_cast<int>(bladeData_struct_nodeNum_data[b_i]) - 1];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].b[1] = 0.5 * chord->data[
        static_cast<int>(bladeData_struct_nodeNum_data[b_i + 1]) - 1];

      // element semi chord
      idx = b_i + bladeData_struct_remaining_size[0] * 11;
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a0[0] =
        bladeData_struct_remaining_data[idx];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a0[1] =
        bladeData_struct_remaining_data[idx + 1];

      // element lift curve slope (needed for flutter analysis)
      // convert "a" to semichord fraction aft of halfchord
      idx = static_cast<int>(c_bladeData_struct_elementNum_d[b_i]) - 1;

      // convert "ac" to semichord fraction foreward of halfchord
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a[0] =
        ((sectionPropsArray->data[idx].a[0] + 0.25 * (2.0 *
           sectionPropsArray->data[idx].b[0])) - sectionPropsArray->data[idx].b
         [0]) / sectionPropsArray->data[idx].b[0];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].a[1] =
        ((sectionPropsArray->data[idx].a[1] + 0.25 * (2.0 *
           sectionPropsArray->data[idx].b[1])) - sectionPropsArray->data[idx].b
         [1]) / sectionPropsArray->data[idx].b[1];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].ac[0] =
        sectionPropsArray->data[idx].ac[0] * 2.0;
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].ac[1] =
        sectionPropsArray->data[idx].ac[1] * 2.0;

      // physical aero center offset from elastic axis
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].aeroCenterOffset[0] =
        sectionPropsArray->data[idx].ac[0] * sectionPropsArray->data[idx].b[0] -
        sectionPropsArray->data[idx].a[0];
      sectionPropsArray->data[static_cast<int>
        (c_bladeData_struct_elementNum_d[b_i]) - 1].aeroCenterOffset[1] =
        sectionPropsArray->data[idx].ac[1] * sectionPropsArray->data[idx].b[1] -
        sectionPropsArray->data[idx].a[1];
    }
  }

  emxFree_real_T(&chord);
  emxInit_real_T(&elLen, 1);
  emxInit_real_T(&psi, 1);
  emxInit_real_T(&theta, 1);
  emxInit_real_T(&roll, 1);
  printf("%s\n", "EIyz, rhoIyz deactivated");
  fflush(stdout);
  cfclose(static_cast<double>(fileid));

  // close element file
  // read element orientation data
  fileid = b_cfopen(ortfile, "rb");
  i = elLen->size[0];
  elLen->size[0] = loop_ub;
  emxEnsureCapacity_real_T(elLen, i);
  i = psi->size[0];
  psi->size[0] = loop_ub;
  emxEnsureCapacity_real_T(psi, i);
  i = theta->size[0];
  theta->size[0] = loop_ub;
  emxEnsureCapacity_real_T(theta, i);
  i = roll->size[0];
  roll->size[0] = loop_ub;
  emxEnsureCapacity_real_T(roll, i);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    getSplitLine(static_cast<double>(fileid), data1);
    elLen->data[b_i] = data1->data[4];
    psi->data[b_i] = data1->data[1];
    theta->data[b_i] = data1->data[2];
    roll->data[b_i] = data1->data[3];
  }

  emxFree_real_T(&data1);

  // store data in element object
  i = el_props->size[0] * el_props->size[1];
  el_props->size[0] = 1;
  el_props->size[1] = sectionPropsArray->size[1];
  emxEnsureCapacity_struct_T(el_props, i);
  idx = sectionPropsArray->size[0] * sectionPropsArray->size[1];
  for (i = 0; i < idx; i++) {
    el_props->data[i] = sectionPropsArray->data[i];
  }

  emxFree_struct_T(&sectionPropsArray);
  i = el_elLen->size[0];
  el_elLen->size[0] = elLen->size[0];
  emxEnsureCapacity_real_T(el_elLen, i);
  idx = elLen->size[0];
  for (i = 0; i < idx; i++) {
    el_elLen->data[i] = elLen->data[i];
  }

  emxFree_real_T(&elLen);
  i = el_psi->size[0];
  el_psi->size[0] = psi->size[0];
  emxEnsureCapacity_real_T(el_psi, i);
  idx = psi->size[0];
  for (i = 0; i < idx; i++) {
    el_psi->data[i] = psi->data[i];
  }

  emxFree_real_T(&psi);
  i = el_theta->size[0];
  el_theta->size[0] = theta->size[0];
  emxEnsureCapacity_real_T(el_theta, i);
  idx = theta->size[0];
  for (i = 0; i < idx; i++) {
    el_theta->data[i] = theta->data[i];
  }

  emxFree_real_T(&theta);
  i = el_roll->size[0];
  el_roll->size[0] = roll->size[0];
  emxEnsureCapacity_real_T(el_roll, i);
  idx = roll->size[0];
  for (i = 0; i < idx; i++) {
    el_roll->data[i] = roll->data[i];
  }

  emxFree_real_T(&roll);
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
