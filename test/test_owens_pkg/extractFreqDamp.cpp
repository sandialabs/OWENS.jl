//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: extractFreqDamp.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "extractFreqDamp.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Declarations
static void c_constructReducedDispVecFromEi(const emxArray_creal_T *vec1, const
  emxArray_real_T *reducedDOFList, double BC_numpBC, const emxArray_real_T
  *BC_pBC, emxArray_creal_T *vec1Red);

// Function Definitions

//
// This function takes the original mode shape and modifies it to
// account for boundary conditions
// Arguments    : const emxArray_creal_T *vec1
//                const emxArray_real_T *reducedDOFList
//                double BC_numpBC
//                const emxArray_real_T *BC_pBC
//                emxArray_creal_T *vec1Red
// Return Type  : void
//
static void c_constructReducedDispVecFromEi(const emxArray_creal_T *vec1, const
  emxArray_real_T *reducedDOFList, double BC_numpBC, const emxArray_real_T
  *BC_pBC, emxArray_creal_T *vec1Red)
{
  emxArray_real_T *bcdoflist;
  int i;
  int k;
  int b_i;
  double b_index;
  double a;
  boolean_T tf;
  boolean_T exitg1;
  double b;
  double absx;
  int exponent;
  emxInit_real_T(&bcdoflist, 2);
  i = bcdoflist->size[0] * bcdoflist->size[1];
  bcdoflist->size[0] = 1;
  k = static_cast<int>(BC_numpBC);
  bcdoflist->size[1] = k;
  emxEnsureCapacity_real_T(bcdoflist, i);

  // form pBC DOF list
  for (b_i = 0; b_i < k; b_i++) {
    bcdoflist->data[b_i] = 0.0;
    bcdoflist->data[b_i] = (BC_pBC->data[b_i] - 1.0) * 6.0 + BC_pBC->data[b_i +
      BC_pBC->size[0]];
  }

  b_index = 1.0;
  i = vec1Red->size[0];
  vec1Red->size[0] = reducedDOFList->size[1];
  emxEnsureCapacity_creal_T(vec1Red, i);
  k = reducedDOFList->size[1];
  for (i = 0; i < k; i++) {
    vec1Red->data[i].re = 0.0;
    vec1Red->data[i].im = 0.0;
  }

  i = reducedDOFList->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    a = reducedDOFList->data[b_i];
    tf = false;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= bcdoflist->size[1] - 1)) {
      b = bcdoflist->data[k];
      absx = std::abs(b / 2.0);
      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &exponent);
          absx = std::ldexp(1.0, exponent - 53);
        }
      } else {
        absx = rtNaN;
      }

      if ((std::abs(bcdoflist->data[k] - a) < absx) || (rtIsInf(a) && rtIsInf(b)
           && ((a > 0.0) == (bcdoflist->data[k] > 0.0)))) {
        tf = true;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (!tf) {
      vec1Red->data[b_i] = vec1->data[static_cast<int>((static_cast<double>
        (vec1->size[0]) / 2.0 + b_index)) - 1];
      b_index++;
    } else {
      vec1Red->data[b_i].re = 0.0;
      vec1Red->data[b_i].im = 0.0;
    }
  }

  emxFree_real_T(&bcdoflist);
}

//
// extractFreqDamp   extract frequency, damping, mode shapes from eigsolution
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freq,damp,phase1,phase2,sortedModes] = extractFreqDamp(val,vec,...
//                                          numDOFPerNode,jointTransform,...
//                                          reducedDOFList,BC,analysisType)
//
//    This function calculates the eigenvalues and vectors of a structural
//    dynamic system.
//
//    input:
//    val            = eigenvalue
//    vec            = eigenvector
//    numDOFPerNode  = number of degrees of freedom per node
//    jointTransform = joint transformation matrix from reduced to full DOF
//                     list
//    reducedDOFList = listing of reduced DOFs
//    BC             = boundary condition object containing boundary
//                     condition info
//    analysisType   = analysis type
//
//    output:
//    freq        = modal frequency
//    damp        = modal damping
//    phase1      = in phase mode shape (real part of mode shape)
//    phase2      = out of phase mode shape (imaginary part of mode shape)
//    sortedModes = total, complex mode shape
// Arguments    : const creal_T val
//                const emxArray_creal_T *vec
//                const emxArray_real_T *jointTransform
//                const emxArray_real_T *reducedDOFList
//                double BC_numpBC
//                const emxArray_real_T *BC_pBC
//                double *freq
//                double *damp
//                emxArray_real_T *phase1
//                emxArray_real_T *phase2
//                emxArray_creal_T *sortedModes
// Return Type  : void
//
void extractFreqDamp(const creal_T val, const emxArray_creal_T *vec, const
                     emxArray_real_T *jointTransform, const emxArray_real_T
                     *reducedDOFList, double BC_numpBC, const emxArray_real_T
                     *BC_pBC, double *freq, double *damp, emxArray_real_T
                     *phase1, emxArray_real_T *phase2, emxArray_creal_T
                     *sortedModes)
{
  double freq_tmp;
  emxArray_creal_T *b;
  emxArray_creal_T *b_jointTransform;
  int i;
  int nx;
  emxArray_creal_T *dispOrig;
  emxArray_int8_T *sortedModes0;
  int k;
  int i1;
  int b_i;
  int j;
  emxArray_real_T *varargin_1;
  double maxval[6];
  boolean_T exitg1;
  double max1;
  double u0;
  double max2;
  double maxOverall;
  boolean_T guard1 = false;
  freq_tmp = std::abs(val.im);
  *freq = freq_tmp / 6.2831853071795862;

  // damped frequency for state space system
  *damp = -val.re / freq_tmp;

  // modal damping
  if (freq_tmp < 0.0001) {
    // if imaginary part of eigenvalue is numeric zero treat as spring-mass system 
    *freq = std::sqrt(std::abs(val.re)) / 6.2831853071795862;
    *damp = 0.0;
  }

  emxInit_creal_T(&b, 1);
  emxInit_creal_T(&b_jointTransform, 2);

  // for all but automated flutter analysis
  //           [len,numModeShapes] = size(vec);
  c_constructReducedDispVecFromEi(vec, reducedDOFList, BC_numpBC, BC_pBC, b);

  // construct mode shape vector with boundary conditinos
  i = b_jointTransform->size[0] * b_jointTransform->size[1];
  b_jointTransform->size[0] = jointTransform->size[0];
  b_jointTransform->size[1] = jointTransform->size[1];
  emxEnsureCapacity_creal_T(b_jointTransform, i);
  nx = jointTransform->size[0] * jointTransform->size[1];
  for (i = 0; i < nx; i++) {
    b_jointTransform->data[i].re = jointTransform->data[i];
    b_jointTransform->data[i].im = 0.0;
  }

  emxInit_creal_T(&dispOrig, 1);
  i = dispOrig->size[0];
  dispOrig->size[0] = b_jointTransform->size[0];
  emxEnsureCapacity_creal_T(dispOrig, i);
  nx = b_jointTransform->size[0];
  for (i = 0; i < nx; i++) {
    dispOrig->data[i].re = 0.0;
    dispOrig->data[i].im = 0.0;
    k = b_jointTransform->size[1];
    for (i1 = 0; i1 < k; i1++) {
      dispOrig->data[i].re += b_jointTransform->data[i + b_jointTransform->size
        [0] * i1].re * b->data[i1].re - b_jointTransform->data[i +
        b_jointTransform->size[0] * i1].im * b->data[i1].im;
      dispOrig->data[i].im += b_jointTransform->data[i + b_jointTransform->size
        [0] * i1].re * b->data[i1].im + b_jointTransform->data[i +
        b_jointTransform->size[0] * i1].im * b->data[i1].re;
    }
  }

  emxFree_creal_T(&b_jointTransform);
  emxFree_creal_T(&b);
  emxInit_int8_T(&sortedModes0, 2);

  // transform from reduced DOF list to full DOF list
  i = sortedModes0->size[0] * sortedModes0->size[1];
  sortedModes0->size[0] = static_cast<int>((static_cast<double>(dispOrig->size[0])
    / 6.0));
  sortedModes0->size[1] = 6;
  emxEnsureCapacity_int8_T(sortedModes0, i);
  nx = static_cast<int>((static_cast<double>(dispOrig->size[0]) / 6.0)) * 6;
  for (i = 0; i < nx; i++) {
    sortedModes0->data[i] = 0;
  }

  i = sortedModes->size[0] * sortedModes->size[1];
  sortedModes->size[0] = sortedModes0->size[0];
  sortedModes->size[1] = 6;
  emxEnsureCapacity_creal_T(sortedModes, i);
  nx = sortedModes0->size[0] * sortedModes0->size[1];
  for (i = 0; i < nx; i++) {
    sortedModes->data[i].re = sortedModes0->data[i];
    sortedModes->data[i].im = 0.0;
  }

  emxFree_int8_T(&sortedModes0);
  i = static_cast<int>((static_cast<double>(dispOrig->size[0]) / 6.0));
  for (b_i = 0; b_i < i; b_i++) {
    // construct matrix of nodal DOF values from full DOF eigenvector
    for (j = 0; j < 6; j++) {
      sortedModes->data[b_i + sortedModes->size[0] * j] = dispOrig->data[
        static_cast<int>((static_cast<unsigned int>((b_i * 6)) + j))];
    }
  }

  emxFree_creal_T(&dispOrig);
  i = phase1->size[0] * phase1->size[1];
  phase1->size[0] = sortedModes->size[0];
  phase1->size[1] = 6;
  emxEnsureCapacity_real_T(phase1, i);
  nx = sortedModes->size[0] * sortedModes->size[1];
  for (i = 0; i < nx; i++) {
    phase1->data[i] = sortedModes->data[i].re;
  }

  // phase 1 is real part of modeshape (0 deg in phase)
  i = phase2->size[0] * phase2->size[1];
  phase2->size[0] = sortedModes->size[0];
  phase2->size[1] = 6;
  emxEnsureCapacity_real_T(phase2, i);
  nx = sortedModes->size[0] * sortedModes->size[1];
  for (i = 0; i < nx; i++) {
    phase2->data[i] = sortedModes->data[i].im;
  }

  emxInit_real_T(&varargin_1, 2);

  // phase 2 is imag part of modeshape (90 deg out of phase)
  nx = phase1->size[0] * 6;
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = phase1->size[0];
  varargin_1->size[1] = 6;
  emxEnsureCapacity_real_T(varargin_1, i);
  for (k = 0; k < nx; k++) {
    varargin_1->data[k] = std::abs(phase1->data[k]);
  }

  nx = varargin_1->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = varargin_1->data[varargin_1->size[0] * j];
    for (b_i = 2; b_i <= nx; b_i++) {
      freq_tmp = varargin_1->data[(b_i + varargin_1->size[0] * j) - 1];
      if ((!rtIsNaN(freq_tmp)) && (rtIsNaN(maxval[j]) || (maxval[j] < freq_tmp)))
      {
        maxval[j] = freq_tmp;
      }
    }
  }

  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k <= 6)) {
      if (!rtIsNaN(maxval[k - 1])) {
        nx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (nx == 0) {
    max1 = maxval[0];
  } else {
    freq_tmp = maxval[nx - 1];
    i = nx + 1;
    for (k = i; k < 7; k++) {
      u0 = maxval[k - 1];
      if (freq_tmp < u0) {
        freq_tmp = u0;
      }
    }

    max1 = freq_tmp;
  }

  // find maximum values for modeshape normalization
  nx = phase2->size[0] * 6;
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = phase2->size[0];
  varargin_1->size[1] = 6;
  emxEnsureCapacity_real_T(varargin_1, i);
  for (k = 0; k < nx; k++) {
    varargin_1->data[k] = std::abs(phase2->data[k]);
  }

  nx = varargin_1->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = varargin_1->data[varargin_1->size[0] * j];
    for (b_i = 2; b_i <= nx; b_i++) {
      freq_tmp = varargin_1->data[(b_i + varargin_1->size[0] * j) - 1];
      if ((!rtIsNaN(freq_tmp)) && (rtIsNaN(maxval[j]) || (maxval[j] < freq_tmp)))
      {
        maxval[j] = freq_tmp;
      }
    }
  }

  emxFree_real_T(&varargin_1);
  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k <= 6)) {
      if (!rtIsNaN(maxval[k - 1])) {
        nx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (nx == 0) {
    max2 = maxval[0];
  } else {
    freq_tmp = maxval[nx - 1];
    i = nx + 1;
    for (k = i; k < 7; k++) {
      u0 = maxval[k - 1];
      if (freq_tmp < u0) {
        freq_tmp = u0;
      }
    }

    max2 = freq_tmp;
  }

  if ((max1 > max2) || rtIsNaN(max2)) {
    maxOverall = max1;
  } else {
    maxOverall = max2;
  }

  if (maxOverall == 0.0) {
    maxOverall = 1.0;
  }

  nx = phase1->size[0] * phase1->size[1];
  i = phase1->size[0] * phase1->size[1];
  phase1->size[1] = 6;
  emxEnsureCapacity_real_T(phase1, i);
  for (i = 0; i < nx; i++) {
    phase1->data[i] /= maxOverall;
  }

  // normalize modeshapes
  nx = phase2->size[0] * phase2->size[1];
  i = phase2->size[0] * phase2->size[1];
  phase2->size[1] = 6;
  emxEnsureCapacity_real_T(phase2, i);
  for (i = 0; i < nx; i++) {
    phase2->data[i] /= maxOverall;
  }

  nx = phase1->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = phase1->data[phase1->size[0] * j];
    for (b_i = 2; b_i <= nx; b_i++) {
      freq_tmp = phase1->data[(b_i + phase1->size[0] * j) - 1];
      if ((!rtIsNaN(freq_tmp)) && (rtIsNaN(maxval[j]) || (maxval[j] > freq_tmp)))
      {
        maxval[j] = freq_tmp;
      }
    }
  }

  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k <= 6)) {
      if (!rtIsNaN(maxval[k - 1])) {
        nx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (nx == 0) {
    freq_tmp = maxval[0];
  } else {
    freq_tmp = maxval[nx - 1];
    i = nx + 1;
    for (k = i; k < 7; k++) {
      u0 = maxval[k - 1];
      if (freq_tmp > u0) {
        freq_tmp = u0;
      }
    }
  }

  guard1 = false;
  if (std::abs(freq_tmp + 1.0) < 0.0001) {
    guard1 = true;
  } else {
    nx = phase2->size[0];
    for (j = 0; j < 6; j++) {
      maxval[j] = phase2->data[phase2->size[0] * j];
      for (b_i = 2; b_i <= nx; b_i++) {
        freq_tmp = phase2->data[(b_i + phase2->size[0] * j) - 1];
        if ((!rtIsNaN(freq_tmp)) && (rtIsNaN(maxval[j]) || (maxval[j] > freq_tmp)))
        {
          maxval[j] = freq_tmp;
        }
      }
    }

    if (!rtIsNaN(maxval[0])) {
      nx = 1;
    } else {
      nx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= 6)) {
        if (!rtIsNaN(maxval[k - 1])) {
          nx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (nx == 0) {
      freq_tmp = maxval[0];
    } else {
      freq_tmp = maxval[nx - 1];
      i = nx + 1;
      for (k = i; k < 7; k++) {
        u0 = maxval[k - 1];
        if (freq_tmp > u0) {
          freq_tmp = u0;
        }
      }
    }

    if (std::abs(freq_tmp + 1.0) < 0.0001) {
      guard1 = true;
    }
  }

  if (guard1) {
    nx = phase1->size[0] * phase1->size[1];
    i = phase1->size[0] * phase1->size[1];
    phase1->size[1] = 6;
    emxEnsureCapacity_real_T(phase1, i);
    for (i = 0; i < nx; i++) {
      phase1->data[i] = -phase1->data[i];
    }

    nx = phase2->size[0] * phase2->size[1];
    i = phase2->size[0] * phase2->size[1];
    phase2->size[1] = 6;
    emxEnsureCapacity_real_T(phase2, i);
    for (i = 0; i < nx; i++) {
      phase2->data[i] = -phase2->data[i];
    }
  }
}

//
// File trailer for extractFreqDamp.cpp
//
// [EOF]
//
