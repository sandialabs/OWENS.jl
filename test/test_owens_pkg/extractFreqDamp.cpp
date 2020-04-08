//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: extractFreqDamp.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
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
  int loop_ub;
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
  loop_ub = static_cast<int>(BC_numpBC);
  bcdoflist->size[1] = loop_ub;
  emxEnsureCapacity_real_T(bcdoflist, i);

  // form pBC DOF list
  for (b_i = 0; b_i < loop_ub; b_i++) {
    bcdoflist->data[b_i] = 0.0;
    bcdoflist->data[b_i] = (BC_pBC->data[b_i] - 1.0) * 6.0 + BC_pBC->data[b_i +
      BC_pBC->size[0]];
  }

  b_index = 1.0;
  i = vec1Red->size[0];
  vec1Red->size[0] = reducedDOFList->size[1];
  emxEnsureCapacity_creal_T(vec1Red, i);
  loop_ub = reducedDOFList->size[1];
  for (i = 0; i < loop_ub; i++) {
    vec1Red->data[i].re = 0.0;
    vec1Red->data[i].im = 0.0;
  }

  i = reducedDOFList->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    a = reducedDOFList->data[b_i];
    tf = false;
    loop_ub = 0;
    exitg1 = false;
    while ((!exitg1) && (loop_ub <= bcdoflist->size[1] - 1)) {
      b = bcdoflist->data[loop_ub];
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

      if ((std::abs(bcdoflist->data[loop_ub] - a) < absx) || (rtIsInf(a) &&
           rtIsInf(b) && ((a > 0.0) == (bcdoflist->data[loop_ub] > 0.0)))) {
        tf = true;
        exitg1 = true;
      } else {
        loop_ub++;
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
  double max2;
  emxArray_creal_T *dispReduc;
  emxArray_creal_T *b_jointTransform;
  int i;
  int nx;
  emxArray_creal_T *dispOrig;
  emxArray_real_T *sortedModes0;
  int loop_ub;
  int j;
  double maxval[6];
  boolean_T exitg1;
  double max1;
  double d;
  boolean_T guard1 = false;
  max2 = std::abs(val.im);
  *freq = max2 / 6.2831853071795862;

  // damped frequency for state space system
  *damp = -val.re / max2;

  // modal damping
  if (max2 < 0.0001) {
    // if imaginary part of eigenvalue is numeric zero treat as spring-mass system 
    *freq = std::sqrt(std::abs(val.re)) / 6.2831853071795862;
    *damp = 0.0;
  }

  emxInit_creal_T(&dispReduc, 1);
  emxInit_creal_T(&b_jointTransform, 2);

  // for all but automated flutter analysis
  //           [len,numModeShapes] = size(vec);
  c_constructReducedDispVecFromEi(vec, reducedDOFList, BC_numpBC, BC_pBC,
    dispReduc);

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
    loop_ub = b_jointTransform->size[1];
    for (j = 0; j < loop_ub; j++) {
      dispOrig->data[i].re += b_jointTransform->data[i + b_jointTransform->size
        [0] * j].re * dispReduc->data[j].re - b_jointTransform->data[i +
        b_jointTransform->size[0] * j].im * dispReduc->data[j].im;
      dispOrig->data[i].im += b_jointTransform->data[i + b_jointTransform->size
        [0] * j].re * dispReduc->data[j].im + b_jointTransform->data[i +
        b_jointTransform->size[0] * j].im * dispReduc->data[j].re;
    }
  }

  emxFree_creal_T(&b_jointTransform);
  emxFree_creal_T(&dispReduc);
  emxInit_real_T(&sortedModes0, 2);

  // transform from reduced DOF list to full DOF list
  i = sortedModes0->size[0] * sortedModes0->size[1];
  sortedModes0->size[0] = static_cast<int>((static_cast<double>(dispOrig->size[0])
    / 6.0));
  sortedModes0->size[1] = 6;
  emxEnsureCapacity_real_T(sortedModes0, i);
  nx = static_cast<int>((static_cast<double>(dispOrig->size[0]) / 6.0)) * 6;
  for (i = 0; i < nx; i++) {
    sortedModes0->data[i] = 0.0;
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

  i = static_cast<int>((static_cast<double>(dispOrig->size[0]) / 6.0));
  for (loop_ub = 0; loop_ub < i; loop_ub++) {
    // construct matrix of nodal DOF values from full DOF eigenvector
    for (j = 0; j < 6; j++) {
      sortedModes->data[loop_ub + sortedModes->size[0] * j] = dispOrig->data[
        static_cast<int>((static_cast<unsigned int>((loop_ub * 6)) + j))];
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

  // phase 2 is imag part of modeshape (90 deg out of phase)
  nx = phase1->size[0] * 6;
  i = sortedModes0->size[0] * sortedModes0->size[1];
  sortedModes0->size[0] = phase1->size[0];
  sortedModes0->size[1] = 6;
  emxEnsureCapacity_real_T(sortedModes0, i);
  for (loop_ub = 0; loop_ub < nx; loop_ub++) {
    sortedModes0->data[loop_ub] = std::abs(phase1->data[loop_ub]);
  }

  nx = sortedModes0->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = sortedModes0->data[sortedModes0->size[0] * j];
    for (loop_ub = 2; loop_ub <= nx; loop_ub++) {
      max2 = sortedModes0->data[(loop_ub + sortedModes0->size[0] * j) - 1];
      if ((!rtIsNaN(max2)) && (rtIsNaN(maxval[j]) || (maxval[j] < max2))) {
        maxval[j] = max2;
      }
    }
  }

  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    loop_ub = 2;
    exitg1 = false;
    while ((!exitg1) && (loop_ub <= 6)) {
      if (!rtIsNaN(maxval[loop_ub - 1])) {
        nx = loop_ub;
        exitg1 = true;
      } else {
        loop_ub++;
      }
    }
  }

  if (nx == 0) {
    max1 = maxval[0];
  } else {
    max1 = maxval[nx - 1];
    i = nx + 1;
    for (loop_ub = i; loop_ub < 7; loop_ub++) {
      d = maxval[loop_ub - 1];
      if (max1 < d) {
        max1 = d;
      }
    }
  }

  // find maximum values for modeshape normalization
  nx = phase2->size[0] * 6;
  i = sortedModes0->size[0] * sortedModes0->size[1];
  sortedModes0->size[0] = phase2->size[0];
  sortedModes0->size[1] = 6;
  emxEnsureCapacity_real_T(sortedModes0, i);
  for (loop_ub = 0; loop_ub < nx; loop_ub++) {
    sortedModes0->data[loop_ub] = std::abs(phase2->data[loop_ub]);
  }

  nx = sortedModes0->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = sortedModes0->data[sortedModes0->size[0] * j];
    for (loop_ub = 2; loop_ub <= nx; loop_ub++) {
      max2 = sortedModes0->data[(loop_ub + sortedModes0->size[0] * j) - 1];
      if ((!rtIsNaN(max2)) && (rtIsNaN(maxval[j]) || (maxval[j] < max2))) {
        maxval[j] = max2;
      }
    }
  }

  emxFree_real_T(&sortedModes0);
  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    loop_ub = 2;
    exitg1 = false;
    while ((!exitg1) && (loop_ub <= 6)) {
      if (!rtIsNaN(maxval[loop_ub - 1])) {
        nx = loop_ub;
        exitg1 = true;
      } else {
        loop_ub++;
      }
    }
  }

  if (nx == 0) {
    max2 = maxval[0];
  } else {
    max2 = maxval[nx - 1];
    i = nx + 1;
    for (loop_ub = i; loop_ub < 7; loop_ub++) {
      d = maxval[loop_ub - 1];
      if (max2 < d) {
        max2 = d;
      }
    }
  }

  if ((max1 > max2) || rtIsNaN(max2)) {
    max2 = max1;
  }

  if (max2 == 0.0) {
    max2 = 1.0;
  }

  nx = phase1->size[0] * phase1->size[1];
  i = phase1->size[0] * phase1->size[1];
  phase1->size[1] = 6;
  emxEnsureCapacity_real_T(phase1, i);
  for (i = 0; i < nx; i++) {
    phase1->data[i] /= max2;
  }

  // normalize modeshapes
  nx = phase2->size[0] * phase2->size[1];
  i = phase2->size[0] * phase2->size[1];
  phase2->size[1] = 6;
  emxEnsureCapacity_real_T(phase2, i);
  for (i = 0; i < nx; i++) {
    phase2->data[i] /= max2;
  }

  nx = phase1->size[0];
  for (j = 0; j < 6; j++) {
    maxval[j] = phase1->data[phase1->size[0] * j];
    for (loop_ub = 2; loop_ub <= nx; loop_ub++) {
      max2 = phase1->data[(loop_ub + phase1->size[0] * j) - 1];
      if ((!rtIsNaN(max2)) && (rtIsNaN(maxval[j]) || (maxval[j] > max2))) {
        maxval[j] = max2;
      }
    }
  }

  if (!rtIsNaN(maxval[0])) {
    nx = 1;
  } else {
    nx = 0;
    loop_ub = 2;
    exitg1 = false;
    while ((!exitg1) && (loop_ub <= 6)) {
      if (!rtIsNaN(maxval[loop_ub - 1])) {
        nx = loop_ub;
        exitg1 = true;
      } else {
        loop_ub++;
      }
    }
  }

  if (nx == 0) {
    max2 = maxval[0];
  } else {
    max2 = maxval[nx - 1];
    i = nx + 1;
    for (loop_ub = i; loop_ub < 7; loop_ub++) {
      d = maxval[loop_ub - 1];
      if (max2 > d) {
        max2 = d;
      }
    }
  }

  guard1 = false;
  if (std::abs(max2 + 1.0) < 0.0001) {
    guard1 = true;
  } else {
    nx = phase2->size[0];
    for (j = 0; j < 6; j++) {
      maxval[j] = phase2->data[phase2->size[0] * j];
      for (loop_ub = 2; loop_ub <= nx; loop_ub++) {
        max2 = phase2->data[(loop_ub + phase2->size[0] * j) - 1];
        if ((!rtIsNaN(max2)) && (rtIsNaN(maxval[j]) || (maxval[j] > max2))) {
          maxval[j] = max2;
        }
      }
    }

    if (!rtIsNaN(maxval[0])) {
      nx = 1;
    } else {
      nx = 0;
      loop_ub = 2;
      exitg1 = false;
      while ((!exitg1) && (loop_ub <= 6)) {
        if (!rtIsNaN(maxval[loop_ub - 1])) {
          nx = loop_ub;
          exitg1 = true;
        } else {
          loop_ub++;
        }
      }
    }

    if (nx == 0) {
      max2 = maxval[0];
    } else {
      max2 = maxval[nx - 1];
      i = nx + 1;
      for (loop_ub = i; loop_ub < 7; loop_ub++) {
        d = maxval[loop_ub - 1];
        if (max2 > d) {
          max2 = d;
        }
      }
    }

    if (std::abs(max2 + 1.0) < 0.0001) {
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
