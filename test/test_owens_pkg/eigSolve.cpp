//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: eigSolve.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "eigSolve.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "schur.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "test_owens_rtwutil.h"
#include "xdlanv2.h"
#include "xzggev.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// eigSolve   Calculates eigenvalues and vectors of structural dynamics rep.
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [eigVec,eigVal = eigSolve(M,C,K,numModes,flag)
//
//    This function calculates the eigenvalues and vectors of a structural
//    dynamic system.
//
//    input:
//    M         = system mass matrix
//    C         = system damping matrix
//    K         = system stiffness matrix
//    numModes  = number of lower system modes to extract
//
//    output:
//    eigVal    = diagonal matrix of eigenvalues
//    eigVec    = matrix of eigenvectors as columns
//    flag      = directs type of eigensolve
//                 1 = all eigenvalues extracted
//                 2 = subset of eigenvalues extracted
// Arguments    : const emxArray_real_T *M
//                const emxArray_real_T *C
//                const emxArray_real_T *K
//                emxArray_creal_T *eigVec
//                emxArray_creal_T *eigVal
// Return Type  : void
//
void eigSolve(const emxArray_real_T *M, const emxArray_real_T *C, const
              emxArray_real_T *K, emxArray_creal_T *eigVec, emxArray_creal_T
              *eigVal)
{
  int lastcol;
  int input_sizes_idx_1;
  emxArray_int8_T *result;
  emxArray_real_T *varargin_2;
  int i;
  int input_sizes_idx_0;
  int sizes_idx_1;
  boolean_T empty_non_axis_sizes;
  int n;
  emxArray_real_T *b_result;
  emxArray_real_T *varargin_1;
  emxArray_real_T *A;
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  boolean_T exitg2;
  emxArray_creal_T *At;
  int exitg1;
  emxArray_creal_T *alpha1;
  emxArray_creal_T *beta1;
  double colnorm;
  double scale;
  double absxk;
  double t;

  //      if(meirovitchFlag)
  //          disp('Solving with Meirovitch approach ...');
  //          [eigVec,eigVal] = eigSolveMeirovitch(M,C,K,flag);
  //      else
  //          disp('Solving with standard state space approach ...');
  if ((M->size[0] == 0) || (M->size[1] == 0)) {
    input_sizes_idx_1 = 0;
  } else {
    lastcol = M->size[0];
    input_sizes_idx_1 = M->size[1];
    if (lastcol > input_sizes_idx_1) {
      input_sizes_idx_1 = lastcol;
    }
  }

  emxInit_int8_T(&result, 2);
  emxInit_real_T(&varargin_2, 2);

  // constructs state space form (with mass matrix inverted)
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = input_sizes_idx_1;
  varargin_2->size[1] = input_sizes_idx_1;
  emxEnsureCapacity_real_T(varargin_2, i);
  input_sizes_idx_0 = input_sizes_idx_1 * input_sizes_idx_1;
  for (i = 0; i < input_sizes_idx_0; i++) {
    varargin_2->data[i] = 0.0;
  }

  if (input_sizes_idx_1 > 0) {
    for (sizes_idx_1 = 0; sizes_idx_1 < input_sizes_idx_1; sizes_idx_1++) {
      varargin_2->data[sizes_idx_1 + varargin_2->size[0] * sizes_idx_1] = 1.0;
    }
  }

  if (input_sizes_idx_1 != 0) {
    input_sizes_idx_0 = input_sizes_idx_1;
  } else if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    input_sizes_idx_0 = varargin_2->size[0];
  } else {
    input_sizes_idx_0 = 0;
    if (varargin_2->size[0] > 0) {
      input_sizes_idx_0 = varargin_2->size[0];
    }
  }

  empty_non_axis_sizes = (input_sizes_idx_0 == 0);
  if ((!empty_non_axis_sizes) && (input_sizes_idx_1 == 0)) {
    input_sizes_idx_1 = 0;
  }

  if (empty_non_axis_sizes || ((varargin_2->size[0] != 0) && (varargin_2->size[1]
        != 0))) {
    sizes_idx_1 = varargin_2->size[1];
  } else {
    sizes_idx_1 = 0;
  }

  i = result->size[0] * result->size[1];
  result->size[0] = input_sizes_idx_0;
  result->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_int8_T(result, i);
  for (i = 0; i < input_sizes_idx_1; i++) {
    for (n = 0; n < input_sizes_idx_0; n++) {
      result->data[n + result->size[0] * i] = 0;
    }
  }

  for (i = 0; i < sizes_idx_1; i++) {
    for (n = 0; n < input_sizes_idx_0; n++) {
      result->data[n + result->size[0] * (i + input_sizes_idx_1)] = static_cast<
        signed char>(varargin_2->data[n + input_sizes_idx_0 * i]);
    }
  }

  emxInit_real_T(&b_result, 2);
  i = b_result->size[0] * b_result->size[1];
  b_result->size[0] = M->size[0];
  b_result->size[1] = M->size[1];
  emxEnsureCapacity_real_T(b_result, i);
  input_sizes_idx_0 = M->size[0] * M->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_result->data[i] = -M->data[i];
  }

  emxInit_real_T(&varargin_1, 2);
  mldivide(b_result, K, varargin_1);
  i = b_result->size[0] * b_result->size[1];
  b_result->size[0] = M->size[0];
  b_result->size[1] = M->size[1];
  emxEnsureCapacity_real_T(b_result, i);
  input_sizes_idx_0 = M->size[0] * M->size[1];
  for (i = 0; i < input_sizes_idx_0; i++) {
    b_result->data[i] = -M->data[i];
  }

  mldivide(b_result, C, varargin_2);
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    lastcol = varargin_1->size[0];
  } else if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    lastcol = varargin_2->size[0];
  } else {
    lastcol = varargin_1->size[0];
    if (lastcol <= 0) {
      lastcol = 0;
    }

    if (varargin_2->size[0] > lastcol) {
      lastcol = varargin_2->size[0];
    }
  }

  empty_non_axis_sizes = (lastcol == 0);
  if (empty_non_axis_sizes || ((varargin_1->size[0] != 0) && (varargin_1->size[1]
        != 0))) {
    input_sizes_idx_1 = varargin_1->size[1];
  } else {
    input_sizes_idx_1 = 0;
  }

  if (empty_non_axis_sizes || ((varargin_2->size[0] != 0) && (varargin_2->size[1]
        != 0))) {
    sizes_idx_1 = varargin_2->size[1];
  } else {
    sizes_idx_1 = 0;
  }

  i = b_result->size[0] * b_result->size[1];
  b_result->size[0] = lastcol;
  b_result->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_real_T(b_result, i);
  for (i = 0; i < input_sizes_idx_1; i++) {
    for (n = 0; n < lastcol; n++) {
      b_result->data[n + b_result->size[0] * i] = varargin_1->data[n + lastcol *
        i];
    }
  }

  for (i = 0; i < sizes_idx_1; i++) {
    for (n = 0; n < lastcol; n++) {
      b_result->data[n + b_result->size[0] * (i + input_sizes_idx_1)] =
        varargin_2->data[n + lastcol * i];
    }
  }

  if ((result->size[0] != 0) && (result->size[1] != 0)) {
    lastcol = result->size[1];
  } else if ((b_result->size[0] != 0) && (b_result->size[1] != 0)) {
    lastcol = b_result->size[1];
  } else {
    lastcol = result->size[1];
    if (lastcol <= 0) {
      lastcol = 0;
    }

    if (b_result->size[1] > lastcol) {
      lastcol = b_result->size[1];
    }
  }

  empty_non_axis_sizes = (lastcol == 0);
  if (empty_non_axis_sizes || ((result->size[0] != 0) && (result->size[1] != 0)))
  {
    input_sizes_idx_0 = result->size[0];
  } else {
    input_sizes_idx_0 = 0;
  }

  if (empty_non_axis_sizes || ((b_result->size[0] != 0) && (b_result->size[1] !=
        0))) {
    sizes_idx_1 = b_result->size[0];
  } else {
    sizes_idx_1 = 0;
  }

  emxInit_real_T(&A, 2);
  i = A->size[0] * A->size[1];
  A->size[0] = input_sizes_idx_0 + sizes_idx_1;
  A->size[1] = lastcol;
  emxEnsureCapacity_real_T(A, i);
  for (i = 0; i < lastcol; i++) {
    for (n = 0; n < input_sizes_idx_0; n++) {
      A->data[n + A->size[0] * i] = result->data[n + input_sizes_idx_0 * i];
    }
  }

  emxFree_int8_T(&result);
  for (i = 0; i < lastcol; i++) {
    for (n = 0; n < sizes_idx_1; n++) {
      A->data[(n + input_sizes_idx_0) + A->size[0] * i] = b_result->data[n +
        sizes_idx_1 * i];
    }
  }

  emxFree_real_T(&b_result);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    i = eigVec->size[0] * eigVec->size[1];
    eigVec->size[0] = A->size[0];
    eigVec->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(eigVec, i);
    input_sizes_idx_0 = A->size[0] * A->size[1];
    for (i = 0; i < input_sizes_idx_0; i++) {
      eigVec->data[i].re = A->data[i];
      eigVec->data[i].im = 0.0;
    }

    i = eigVal->size[0] * eigVal->size[1];
    eigVal->size[0] = A->size[0];
    eigVal->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(eigVal, i);
    input_sizes_idx_0 = A->size[0] * A->size[1];
    for (i = 0; i < input_sizes_idx_0; i++) {
      eigVal->data[i].re = A->data[i];
      eigVal->data[i].im = 0.0;
    }
  } else {
    input_sizes_idx_0 = A->size[0] * A->size[1];
    empty_non_axis_sizes = true;
    for (sizes_idx_1 = 0; sizes_idx_1 < input_sizes_idx_0; sizes_idx_1++) {
      if ((!empty_non_axis_sizes) || (rtIsInf(A->data[sizes_idx_1]) || rtIsNaN
           (A->data[sizes_idx_1]))) {
        empty_non_axis_sizes = false;
      }
    }

    if (!empty_non_axis_sizes) {
      if ((A->size[0] == 1) && (A->size[1] == 1)) {
        i = eigVec->size[0] * eigVec->size[1];
        eigVec->size[0] = 1;
        eigVec->size[1] = 1;
        emxEnsureCapacity_creal_T(eigVec, i);
        for (i = 0; i < 1; i++) {
          eigVec->data[0].re = rtNaN;
          eigVec->data[0].im = 0.0;
        }

        i = eigVal->size[0] * eigVal->size[1];
        eigVal->size[0] = A->size[0];
        eigVal->size[1] = A->size[1];
        emxEnsureCapacity_creal_T(eigVal, i);
        input_sizes_idx_0 = A->size[0] * A->size[1];
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVal->data[i].re = rtNaN;
          eigVal->data[i].im = 0.0;
        }
      } else {
        unnamed_idx_0 = static_cast<unsigned int>(A->size[0]);
        unnamed_idx_1 = static_cast<unsigned int>(A->size[1]);
        i = eigVec->size[0] * eigVec->size[1];
        eigVec->size[0] = static_cast<int>(unnamed_idx_0);
        eigVec->size[1] = static_cast<int>(unnamed_idx_1);
        emxEnsureCapacity_creal_T(eigVec, i);
        input_sizes_idx_0 = static_cast<int>(unnamed_idx_0) * static_cast<int>
          (unnamed_idx_1);
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVec->data[i].re = rtNaN;
          eigVec->data[i].im = 0.0;
        }

        unnamed_idx_0 = static_cast<unsigned int>(A->size[0]);
        unnamed_idx_1 = static_cast<unsigned int>(A->size[1]);
        i = eigVal->size[0] * eigVal->size[1];
        eigVal->size[0] = static_cast<int>(unnamed_idx_0);
        eigVal->size[1] = static_cast<int>(unnamed_idx_1);
        emxEnsureCapacity_creal_T(eigVal, i);
        input_sizes_idx_0 = static_cast<int>(unnamed_idx_0) * static_cast<int>
          (unnamed_idx_1);
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVal->data[i].re = 0.0;
          eigVal->data[i].im = 0.0;
        }

        i = static_cast<int>(unnamed_idx_0);
        for (sizes_idx_1 = 0; sizes_idx_1 < i; sizes_idx_1++) {
          eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re = rtNaN;
          eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im = 0.0;
        }
      }
    } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
      i = eigVec->size[0] * eigVec->size[1];
      eigVec->size[0] = 1;
      eigVec->size[1] = 1;
      emxEnsureCapacity_creal_T(eigVec, i);
      eigVec->data[0].re = 1.0;
      eigVec->data[0].im = 0.0;
      i = eigVal->size[0] * eigVal->size[1];
      eigVal->size[0] = A->size[0];
      eigVal->size[1] = A->size[1];
      emxEnsureCapacity_creal_T(eigVal, i);
      input_sizes_idx_0 = A->size[0] * A->size[1];
      for (i = 0; i < input_sizes_idx_0; i++) {
        eigVal->data[i].re = A->data[i];
        eigVal->data[i].im = 0.0;
      }
    } else {
      empty_non_axis_sizes = (A->size[0] == A->size[1]);
      if (empty_non_axis_sizes) {
        sizes_idx_1 = 0;
        exitg2 = false;
        while ((!exitg2) && (sizes_idx_1 <= A->size[1] - 1)) {
          input_sizes_idx_0 = 0;
          do {
            exitg1 = 0;
            if (input_sizes_idx_0 <= sizes_idx_1) {
              if (!(A->data[input_sizes_idx_0 + A->size[0] * sizes_idx_1] ==
                    A->data[sizes_idx_1 + A->size[0] * input_sizes_idx_0])) {
                empty_non_axis_sizes = false;
                exitg1 = 1;
              } else {
                input_sizes_idx_0++;
              }
            } else {
              sizes_idx_1++;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      }

      if (empty_non_axis_sizes) {
        schur(A, varargin_2, varargin_1);
        i = eigVec->size[0] * eigVec->size[1];
        eigVec->size[0] = varargin_2->size[0];
        eigVec->size[1] = varargin_2->size[1];
        emxEnsureCapacity_creal_T(eigVec, i);
        input_sizes_idx_0 = varargin_2->size[0] * varargin_2->size[1];
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVec->data[i].re = varargin_2->data[i];
          eigVec->data[i].im = 0.0;
        }

        n = varargin_1->size[0];
        for (sizes_idx_1 = 2; sizes_idx_1 <= n; sizes_idx_1++) {
          varargin_1->data[(sizes_idx_1 + varargin_1->size[0] * (sizes_idx_1 - 2))
            - 1] = 0.0;
          for (input_sizes_idx_0 = 0; input_sizes_idx_0 <= sizes_idx_1 - 2;
               input_sizes_idx_0++) {
            varargin_1->data[input_sizes_idx_0 + varargin_1->size[0] *
              (sizes_idx_1 - 1)] = 0.0;
          }
        }

        i = eigVal->size[0] * eigVal->size[1];
        eigVal->size[0] = varargin_1->size[0];
        eigVal->size[1] = varargin_1->size[1];
        emxEnsureCapacity_creal_T(eigVal, i);
        input_sizes_idx_0 = varargin_1->size[0] * varargin_1->size[1];
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVal->data[i].re = varargin_1->data[i];
          eigVal->data[i].im = 0.0;
        }
      } else {
        emxInit_creal_T(&At, 2);
        i = At->size[0] * At->size[1];
        At->size[0] = A->size[0];
        At->size[1] = A->size[1];
        emxEnsureCapacity_creal_T(At, i);
        input_sizes_idx_0 = A->size[0] * A->size[1];
        for (i = 0; i < input_sizes_idx_0; i++) {
          At->data[i].re = A->data[i];
          At->data[i].im = 0.0;
        }

        emxInit_creal_T(&alpha1, 1);
        emxInit_creal_T(&beta1, 1);
        xzggev(At, &input_sizes_idx_0, alpha1, beta1, eigVec);
        n = A->size[0];
        lastcol = (A->size[0] - 1) * A->size[0] + 1;
        emxFree_creal_T(&At);
        for (input_sizes_idx_1 = 1; n < 0 ? input_sizes_idx_1 >= lastcol :
             input_sizes_idx_1 <= lastcol; input_sizes_idx_1 += n) {
          colnorm = 0.0;
          if (n == 1) {
            colnorm = rt_hypotd_snf(eigVec->data[input_sizes_idx_1 - 1].re,
              eigVec->data[input_sizes_idx_1 - 1].im);
          } else {
            scale = 3.3121686421112381E-170;
            input_sizes_idx_0 = (input_sizes_idx_1 + n) - 1;
            for (sizes_idx_1 = input_sizes_idx_1; sizes_idx_1 <=
                 input_sizes_idx_0; sizes_idx_1++) {
              absxk = std::abs(eigVec->data[sizes_idx_1 - 1].re);
              if (absxk > scale) {
                t = scale / absxk;
                colnorm = colnorm * t * t + 1.0;
                scale = absxk;
              } else {
                t = absxk / scale;
                colnorm += t * t;
              }

              absxk = std::abs(eigVec->data[sizes_idx_1 - 1].im);
              if (absxk > scale) {
                t = scale / absxk;
                colnorm = colnorm * t * t + 1.0;
                scale = absxk;
              } else {
                t = absxk / scale;
                colnorm += t * t;
              }
            }

            colnorm = scale * std::sqrt(colnorm);
          }

          i = (input_sizes_idx_1 + n) - 1;
          for (sizes_idx_1 = input_sizes_idx_1; sizes_idx_1 <= i; sizes_idx_1++)
          {
            absxk = eigVec->data[sizes_idx_1 - 1].re;
            scale = eigVec->data[sizes_idx_1 - 1].im;
            if (scale == 0.0) {
              absxk /= colnorm;
              scale = 0.0;
            } else if (absxk == 0.0) {
              absxk = 0.0;
              scale /= colnorm;
            } else {
              absxk /= colnorm;
              scale /= colnorm;
            }

            eigVec->data[sizes_idx_1 - 1].re = absxk;
            eigVec->data[sizes_idx_1 - 1].im = scale;
          }
        }

        i = eigVal->size[0] * eigVal->size[1];
        eigVal->size[0] = alpha1->size[0];
        eigVal->size[1] = alpha1->size[0];
        emxEnsureCapacity_creal_T(eigVal, i);
        input_sizes_idx_0 = alpha1->size[0] * alpha1->size[0];
        for (i = 0; i < input_sizes_idx_0; i++) {
          eigVal->data[i].re = 0.0;
          eigVal->data[i].im = 0.0;
        }

        i = alpha1->size[0];
        for (sizes_idx_1 = 0; sizes_idx_1 < i; sizes_idx_1++) {
          if (beta1->data[sizes_idx_1].im == 0.0) {
            if (alpha1->data[sizes_idx_1].im == 0.0) {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                alpha1->data[sizes_idx_1].re / beta1->data[sizes_idx_1].re;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im = 0.0;
            } else if (alpha1->data[sizes_idx_1].re == 0.0) {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re = 0.0;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                alpha1->data[sizes_idx_1].im / beta1->data[sizes_idx_1].re;
            } else {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                alpha1->data[sizes_idx_1].re / beta1->data[sizes_idx_1].re;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                alpha1->data[sizes_idx_1].im / beta1->data[sizes_idx_1].re;
            }
          } else if (beta1->data[sizes_idx_1].re == 0.0) {
            if (alpha1->data[sizes_idx_1].re == 0.0) {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                alpha1->data[sizes_idx_1].im / beta1->data[sizes_idx_1].im;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im = 0.0;
            } else if (alpha1->data[sizes_idx_1].im == 0.0) {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re = 0.0;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                -(alpha1->data[sizes_idx_1].re / beta1->data[sizes_idx_1].im);
            } else {
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                alpha1->data[sizes_idx_1].im / beta1->data[sizes_idx_1].im;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                -(alpha1->data[sizes_idx_1].re / beta1->data[sizes_idx_1].im);
            }
          } else {
            t = std::abs(beta1->data[sizes_idx_1].re);
            scale = std::abs(beta1->data[sizes_idx_1].im);
            if (t > scale) {
              scale = beta1->data[sizes_idx_1].im / beta1->data[sizes_idx_1].re;
              absxk = beta1->data[sizes_idx_1].re + scale * beta1->
                data[sizes_idx_1].im;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                (alpha1->data[sizes_idx_1].re + scale * alpha1->data[sizes_idx_1]
                 .im) / absxk;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                (alpha1->data[sizes_idx_1].im - scale * alpha1->data[sizes_idx_1]
                 .re) / absxk;
            } else if (scale == t) {
              if (beta1->data[sizes_idx_1].re > 0.0) {
                scale = 0.5;
              } else {
                scale = -0.5;
              }

              if (beta1->data[sizes_idx_1].im > 0.0) {
                absxk = 0.5;
              } else {
                absxk = -0.5;
              }

              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                (alpha1->data[sizes_idx_1].re * scale + alpha1->data[sizes_idx_1]
                 .im * absxk) / t;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                (alpha1->data[sizes_idx_1].im * scale - alpha1->data[sizes_idx_1]
                 .re * absxk) / t;
            } else {
              scale = beta1->data[sizes_idx_1].re / beta1->data[sizes_idx_1].im;
              absxk = beta1->data[sizes_idx_1].im + scale * beta1->
                data[sizes_idx_1].re;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].re =
                (scale * alpha1->data[sizes_idx_1].re + alpha1->data[sizes_idx_1]
                 .im) / absxk;
              eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1].im =
                (scale * alpha1->data[sizes_idx_1].im - alpha1->data[sizes_idx_1]
                 .re) / absxk;
            }
          }
        }

        emxFree_creal_T(&beta1);
        emxFree_creal_T(&alpha1);
      }
    }
  }

  emxFree_real_T(&A);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&varargin_2);

  // full eigenvalue solve
  //      end
  //  function [eigVec,eigVal] = eigSolveMeirovitch(M,C,K,flag)
  //      	[len,dum] = size(M);
  //
  //      %Force symmetry and skewsymmetry
  //      K=0.5*(K+K');
  //      M=0.5*(M+M');
  //      for i=1:len
  //         for j=1:i
  //             if(i~=j)
  //             C(i,j) = -C(j,i);
  //             end
  //         end
  //      end
  //
  //      Mstar = [M, zeros(len);
  //               zeros(len), K];
  //
  //      Gstar = [C, K;-K,zeros(len)];
  //
  //      Rstar = chol(Mstar);
  //      invRstar = inv(Rstar);
  //
  //      Kstar = Gstar'*(inv(Mstar))*Gstar;
  //      Kstar = 0.5*(Kstar+Kstar');
  //
  //  % 	sysMat = [zeros(len), eye(len);
  //  % 			  -(M^-1)*K, -(M^-1)*C];
  //
  //      sysMat = invRstar'*Kstar*invRstar;
  //      %force symmetry of system matrix
  //      sysMat = 0.5*(sysMat'+sysMat);
  //
  //      if(flag==1)
  //          [eigVec,eigVal] = eig(sysMat);
  //      end
  //      if(flag==2)
  //           sysMat=sparse(sysMat);
  //          [eigVec,eigVal] = eigs(sysMat,20,'SM');
  //      end
  //
  //
  //      %sort and reconstruct typical vectors from meirovitch solve
  //      [len,dum]=size(eigVal);
  //      for i=1:len
  //         valtemp(i)=eigVal(i,i);
  //      end
  //
  //      [valtemp,map,posIndex] = bubbleSort(valtemp);
  //
  //      len = 20;
  //      for i=1:len
  //          vecNew(:,i)=eigVec(:,map(i));
  //      end
  //
  //
  //      for i=1:len
  //  %         eigVec(:,i) = inv(Rstar)*eigVec(:,i);
  //         vecNew(:,i) = invRstar*vecNew(:,i);
  //      end
  //
  //      eigVec=vecNew;
  //
  //      for i=1:len
  //      eigValNew(i,i)=valtemp(i);
  //      end
  //      eigVal = eigValNew;
  //
  //      for i=1:len/2
  //          index1=(i-1)*2+1;
  //          index2=(i-1)*2+2;
  //
  //          vec1 = eigVec(:,index1);
  //          vec2 = eigVec(:,index2);
  //
  //          temp1 = vec1 + 1i.*vec2;
  //          temp2 = vec1 - 1i.*vec2;
  //
  //          eigVec(:,index1) = temp1;
  //          eigVec(:,index2) = temp2;
  //      end
  //  end
}

//
// File trailer for eigSolve.cpp
//
// [EOF]
//
