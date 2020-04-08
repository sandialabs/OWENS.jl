//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: eigSolve.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "eigSolve.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "schur.h"
#include "sort.h"
#include "sortLE.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "test_owens_rtwutil.h"
#include "xzggev.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// ,numModes,flag)
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
  int n;
  int input_sizes_idx_1;
  emxArray_int8_T *result;
  emxArray_real_T *varargin_2;
  int i;
  int sizes_idx_1;
  int input_sizes_idx_0;
  boolean_T empty_non_axis_sizes;
  int coltop;
  emxArray_real_T *b_result;
  emxArray_real_T *T;
  emxArray_real_T *A;
  emxArray_creal_T *eigVec0;
  emxArray_creal_T *sorted_eigVal;
  emxArray_creal_T *alpha1;
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  boolean_T exitg2;
  emxArray_creal_T *At;
  int exitg1;
  emxArray_creal_T *beta1;
  emxArray_real_T *d;
  double colnorm;
  double scale;
  emxArray_int32_T *iidx;
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
    n = M->size[0];
    input_sizes_idx_1 = M->size[1];
    if (n > input_sizes_idx_1) {
      input_sizes_idx_1 = n;
    }
  }

  emxInit_int8_T(&result, 2);
  emxInit_real_T(&varargin_2, 2);

  //  eyelen = eye(len);
  //  zeroslen = zeros(len);
  //  if(flag == 1)
  // constructs state space form (with mass matrix inverted)
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = input_sizes_idx_1;
  varargin_2->size[1] = input_sizes_idx_1;
  emxEnsureCapacity_real_T(varargin_2, i);
  sizes_idx_1 = input_sizes_idx_1 * input_sizes_idx_1;
  for (i = 0; i < sizes_idx_1; i++) {
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
    for (coltop = 0; coltop < input_sizes_idx_0; coltop++) {
      result->data[coltop + result->size[0] * i] = 0;
    }
  }

  for (i = 0; i < sizes_idx_1; i++) {
    for (coltop = 0; coltop < input_sizes_idx_0; coltop++) {
      result->data[coltop + result->size[0] * (i + input_sizes_idx_1)] =
        static_cast<signed char>(varargin_2->data[coltop + input_sizes_idx_0 * i]);
    }
  }

  emxInit_real_T(&b_result, 2);
  i = b_result->size[0] * b_result->size[1];
  b_result->size[0] = M->size[0];
  b_result->size[1] = M->size[1];
  emxEnsureCapacity_real_T(b_result, i);
  sizes_idx_1 = M->size[0] * M->size[1];
  for (i = 0; i < sizes_idx_1; i++) {
    b_result->data[i] = -M->data[i];
  }

  emxInit_real_T(&T, 2);
  mldivide(b_result, K, T);
  i = b_result->size[0] * b_result->size[1];
  b_result->size[0] = M->size[0];
  b_result->size[1] = M->size[1];
  emxEnsureCapacity_real_T(b_result, i);
  sizes_idx_1 = M->size[0] * M->size[1];
  for (i = 0; i < sizes_idx_1; i++) {
    b_result->data[i] = -M->data[i];
  }

  mldivide(b_result, C, varargin_2);
  if ((T->size[0] != 0) && (T->size[1] != 0)) {
    n = T->size[0];
  } else if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    n = varargin_2->size[0];
  } else {
    n = T->size[0];
    if (n <= 0) {
      n = 0;
    }

    if (varargin_2->size[0] > n) {
      n = varargin_2->size[0];
    }
  }

  empty_non_axis_sizes = (n == 0);
  if (empty_non_axis_sizes || ((T->size[0] != 0) && (T->size[1] != 0))) {
    input_sizes_idx_1 = T->size[1];
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
  b_result->size[0] = n;
  b_result->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_real_T(b_result, i);
  for (i = 0; i < input_sizes_idx_1; i++) {
    for (coltop = 0; coltop < n; coltop++) {
      b_result->data[coltop + b_result->size[0] * i] = T->data[coltop + n * i];
    }
  }

  for (i = 0; i < sizes_idx_1; i++) {
    for (coltop = 0; coltop < n; coltop++) {
      b_result->data[coltop + b_result->size[0] * (i + input_sizes_idx_1)] =
        varargin_2->data[coltop + n * i];
    }
  }

  if ((result->size[0] != 0) && (result->size[1] != 0)) {
    n = result->size[1];
  } else if ((b_result->size[0] != 0) && (b_result->size[1] != 0)) {
    n = b_result->size[1];
  } else {
    n = result->size[1];
    if (n <= 0) {
      n = 0;
    }

    if (b_result->size[1] > n) {
      n = b_result->size[1];
    }
  }

  empty_non_axis_sizes = (n == 0);
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

  //  end
  //  if(flag == 2)
  //  sysMat2 = [zeroslen, eyelen;      %construct state space form and lets eigs invert mass matrix 
  //      -K, -C];
  //  sysMat2 = sparse(sysMat2);
  //
  //  MMat2 = [eyelen,zeroslen;zeroslen,M];
  //  MMat2 = sparse(MMat2);
  //  end
  //  if(flag==1)
  i = A->size[0] * A->size[1];
  A->size[0] = input_sizes_idx_0 + sizes_idx_1;
  A->size[1] = n;
  emxEnsureCapacity_real_T(A, i);
  for (i = 0; i < n; i++) {
    for (coltop = 0; coltop < input_sizes_idx_0; coltop++) {
      A->data[coltop + A->size[0] * i] = result->data[coltop + input_sizes_idx_0
        * i];
    }
  }

  emxFree_int8_T(&result);
  for (i = 0; i < n; i++) {
    for (coltop = 0; coltop < sizes_idx_1; coltop++) {
      A->data[(coltop + input_sizes_idx_0) + A->size[0] * i] = b_result->
        data[coltop + sizes_idx_1 * i];
    }
  }

  emxFree_real_T(&b_result);
  emxInit_creal_T(&eigVec0, 2);
  emxInit_creal_T(&sorted_eigVal, 1);
  emxInit_creal_T(&alpha1, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    i = eigVec0->size[0] * eigVec0->size[1];
    eigVec0->size[0] = A->size[0];
    eigVec0->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(eigVec0, i);
    sizes_idx_1 = A->size[0] * A->size[1];
    for (i = 0; i < sizes_idx_1; i++) {
      eigVec0->data[i].re = A->data[i];
      eigVec0->data[i].im = 0.0;
    }

    i = sorted_eigVal->size[0];
    sorted_eigVal->size[0] = A->size[0];
    emxEnsureCapacity_creal_T(sorted_eigVal, i);
    sizes_idx_1 = A->size[0];
    for (i = 0; i < sizes_idx_1; i++) {
      sorted_eigVal->data[i].re = 0.0;
      sorted_eigVal->data[i].im = 0.0;
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
        i = eigVec0->size[0] * eigVec0->size[1];
        eigVec0->size[0] = 1;
        eigVec0->size[1] = 1;
        emxEnsureCapacity_creal_T(eigVec0, i);
        for (i = 0; i < 1; i++) {
          eigVec0->data[0].re = rtNaN;
          eigVec0->data[0].im = 0.0;
        }

        i = sorted_eigVal->size[0];
        sorted_eigVal->size[0] = 1;
        emxEnsureCapacity_creal_T(sorted_eigVal, i);
        for (i = 0; i < 1; i++) {
          sorted_eigVal->data[0].re = rtNaN;
          sorted_eigVal->data[0].im = 0.0;
        }
      } else {
        unnamed_idx_0 = static_cast<unsigned int>(A->size[0]);
        unnamed_idx_1 = static_cast<unsigned int>(A->size[1]);
        i = eigVec0->size[0] * eigVec0->size[1];
        eigVec0->size[0] = static_cast<int>(unnamed_idx_0);
        eigVec0->size[1] = static_cast<int>(unnamed_idx_1);
        emxEnsureCapacity_creal_T(eigVec0, i);
        sizes_idx_1 = static_cast<int>(unnamed_idx_0) * static_cast<int>
          (unnamed_idx_1);
        for (i = 0; i < sizes_idx_1; i++) {
          eigVec0->data[i].re = rtNaN;
          eigVec0->data[i].im = 0.0;
        }

        i = sorted_eigVal->size[0];
        sorted_eigVal->size[0] = A->size[0];
        emxEnsureCapacity_creal_T(sorted_eigVal, i);
        sizes_idx_1 = A->size[0];
        for (i = 0; i < sizes_idx_1; i++) {
          sorted_eigVal->data[i].re = rtNaN;
          sorted_eigVal->data[i].im = 0.0;
        }
      }
    } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
      i = eigVec0->size[0] * eigVec0->size[1];
      eigVec0->size[0] = 1;
      eigVec0->size[1] = 1;
      emxEnsureCapacity_creal_T(eigVec0, i);
      eigVec0->data[0].re = 1.0;
      eigVec0->data[0].im = 0.0;
      i = sorted_eigVal->size[0];
      sorted_eigVal->size[0] = 1;
      emxEnsureCapacity_creal_T(sorted_eigVal, i);
      sorted_eigVal->data[0].re = A->data[0];
      sorted_eigVal->data[0].im = 0.0;
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
        schur(A, varargin_2, T);
        i = eigVec0->size[0] * eigVec0->size[1];
        eigVec0->size[0] = varargin_2->size[0];
        eigVec0->size[1] = varargin_2->size[1];
        emxEnsureCapacity_creal_T(eigVec0, i);
        sizes_idx_1 = varargin_2->size[0] * varargin_2->size[1];
        for (i = 0; i < sizes_idx_1; i++) {
          eigVec0->data[i].re = varargin_2->data[i];
          eigVec0->data[i].im = 0.0;
        }

        emxInit_real_T(&d, 1);
        n = T->size[0];
        i = d->size[0];
        d->size[0] = T->size[0];
        emxEnsureCapacity_real_T(d, i);
        for (sizes_idx_1 = 0; sizes_idx_1 < n; sizes_idx_1++) {
          d->data[sizes_idx_1] = T->data[sizes_idx_1 + T->size[0] * sizes_idx_1];
        }

        i = sorted_eigVal->size[0];
        sorted_eigVal->size[0] = d->size[0];
        emxEnsureCapacity_creal_T(sorted_eigVal, i);
        sizes_idx_1 = d->size[0];
        for (i = 0; i < sizes_idx_1; i++) {
          sorted_eigVal->data[i].re = d->data[i];
          sorted_eigVal->data[i].im = 0.0;
        }

        emxFree_real_T(&d);
      } else {
        emxInit_creal_T(&At, 2);
        i = At->size[0] * At->size[1];
        At->size[0] = A->size[0];
        At->size[1] = A->size[1];
        emxEnsureCapacity_creal_T(At, i);
        sizes_idx_1 = A->size[0] * A->size[1];
        for (i = 0; i < sizes_idx_1; i++) {
          At->data[i].re = A->data[i];
          At->data[i].im = 0.0;
        }

        emxInit_creal_T(&beta1, 1);
        xzggev(At, &input_sizes_idx_0, alpha1, beta1, eigVec0);
        n = A->size[0];
        input_sizes_idx_1 = (A->size[0] - 1) * A->size[0] + 1;
        emxFree_creal_T(&At);
        for (coltop = 1; n < 0 ? coltop >= input_sizes_idx_1 : coltop <=
             input_sizes_idx_1; coltop += n) {
          colnorm = 0.0;
          if (n == 1) {
            colnorm = rt_hypotd_snf(eigVec0->data[coltop - 1].re, eigVec0->
              data[coltop - 1].im);
          } else {
            scale = 3.3121686421112381E-170;
            input_sizes_idx_0 = (coltop + n) - 1;
            for (sizes_idx_1 = coltop; sizes_idx_1 <= input_sizes_idx_0;
                 sizes_idx_1++) {
              absxk = std::abs(eigVec0->data[sizes_idx_1 - 1].re);
              if (absxk > scale) {
                t = scale / absxk;
                colnorm = colnorm * t * t + 1.0;
                scale = absxk;
              } else {
                t = absxk / scale;
                colnorm += t * t;
              }

              absxk = std::abs(eigVec0->data[sizes_idx_1 - 1].im);
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

          i = (coltop + n) - 1;
          for (sizes_idx_1 = coltop; sizes_idx_1 <= i; sizes_idx_1++) {
            absxk = eigVec0->data[sizes_idx_1 - 1].re;
            scale = eigVec0->data[sizes_idx_1 - 1].im;
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

            eigVec0->data[sizes_idx_1 - 1].re = absxk;
            eigVec0->data[sizes_idx_1 - 1].im = scale;
          }
        }

        i = sorted_eigVal->size[0];
        sorted_eigVal->size[0] = alpha1->size[0];
        emxEnsureCapacity_creal_T(sorted_eigVal, i);
        sizes_idx_1 = alpha1->size[0];
        for (i = 0; i < sizes_idx_1; i++) {
          if (beta1->data[i].im == 0.0) {
            if (alpha1->data[i].im == 0.0) {
              sorted_eigVal->data[i].re = alpha1->data[i].re / beta1->data[i].re;
              sorted_eigVal->data[i].im = 0.0;
            } else if (alpha1->data[i].re == 0.0) {
              sorted_eigVal->data[i].re = 0.0;
              sorted_eigVal->data[i].im = alpha1->data[i].im / beta1->data[i].re;
            } else {
              sorted_eigVal->data[i].re = alpha1->data[i].re / beta1->data[i].re;
              sorted_eigVal->data[i].im = alpha1->data[i].im / beta1->data[i].re;
            }
          } else if (beta1->data[i].re == 0.0) {
            if (alpha1->data[i].re == 0.0) {
              sorted_eigVal->data[i].re = alpha1->data[i].im / beta1->data[i].im;
              sorted_eigVal->data[i].im = 0.0;
            } else if (alpha1->data[i].im == 0.0) {
              sorted_eigVal->data[i].re = 0.0;
              sorted_eigVal->data[i].im = -(alpha1->data[i].re / beta1->data[i].
                im);
            } else {
              sorted_eigVal->data[i].re = alpha1->data[i].im / beta1->data[i].im;
              sorted_eigVal->data[i].im = -(alpha1->data[i].re / beta1->data[i].
                im);
            }
          } else {
            t = std::abs(beta1->data[i].re);
            scale = std::abs(beta1->data[i].im);
            if (t > scale) {
              scale = beta1->data[i].im / beta1->data[i].re;
              absxk = beta1->data[i].re + scale * beta1->data[i].im;
              sorted_eigVal->data[i].re = (alpha1->data[i].re + scale *
                alpha1->data[i].im) / absxk;
              sorted_eigVal->data[i].im = (alpha1->data[i].im - scale *
                alpha1->data[i].re) / absxk;
            } else if (scale == t) {
              if (beta1->data[i].re > 0.0) {
                scale = 0.5;
              } else {
                scale = -0.5;
              }

              if (beta1->data[i].im > 0.0) {
                absxk = 0.5;
              } else {
                absxk = -0.5;
              }

              sorted_eigVal->data[i].re = (alpha1->data[i].re * scale +
                alpha1->data[i].im * absxk) / t;
              sorted_eigVal->data[i].im = (alpha1->data[i].im * scale -
                alpha1->data[i].re * absxk) / t;
            } else {
              scale = beta1->data[i].re / beta1->data[i].im;
              absxk = beta1->data[i].im + scale * beta1->data[i].re;
              sorted_eigVal->data[i].re = (scale * alpha1->data[i].re +
                alpha1->data[i].im) / absxk;
              sorted_eigVal->data[i].im = (scale * alpha1->data[i].im -
                alpha1->data[i].re) / absxk;
            }
          }
        }

        emxFree_creal_T(&beta1);
      }
    }
  }

  emxFree_real_T(&T);
  emxFree_real_T(&A);
  emxFree_real_T(&varargin_2);

  // full eigenvalue solve
  //  test0 = sysMat*eigVec0 - eigVec0*diag(eigVal0);
  //  output = real(test0(1,1));
  i = alpha1->size[0];
  alpha1->size[0] = sorted_eigVal->size[0];
  emxEnsureCapacity_creal_T(alpha1, i);
  sizes_idx_1 = sorted_eigVal->size[0];
  for (i = 0; i < sizes_idx_1; i++) {
    alpha1->data[i] = sorted_eigVal->data[i];
  }

  emxFree_creal_T(&sorted_eigVal);
  emxInit_int32_T(&iidx, 1);
  b_sort(alpha1, iidx);
  sizes_idx_1 = eigVec0->size[0];
  i = eigVec->size[0] * eigVec->size[1];
  eigVec->size[0] = eigVec0->size[0];
  eigVec->size[1] = iidx->size[0];
  emxEnsureCapacity_creal_T(eigVec, i);
  input_sizes_idx_0 = iidx->size[0];
  for (i = 0; i < input_sizes_idx_0; i++) {
    for (coltop = 0; coltop < sizes_idx_1; coltop++) {
      eigVec->data[coltop + eigVec->size[0] * i] = eigVec0->data[coltop +
        eigVec0->size[0] * (iidx->data[i] - 1)];
    }
  }

  emxFree_int32_T(&iidx);
  emxFree_creal_T(&eigVec0);
  input_sizes_idx_0 = alpha1->size[0];
  i = eigVal->size[0] * eigVal->size[1];
  eigVal->size[0] = alpha1->size[0];
  eigVal->size[1] = alpha1->size[0];
  emxEnsureCapacity_creal_T(eigVal, i);
  sizes_idx_1 = alpha1->size[0] * alpha1->size[0];
  for (i = 0; i < sizes_idx_1; i++) {
    eigVal->data[i].re = 0.0;
    eigVal->data[i].im = 0.0;
  }

  for (sizes_idx_1 = 0; sizes_idx_1 < input_sizes_idx_0; sizes_idx_1++) {
    eigVal->data[sizes_idx_1 + eigVal->size[0] * sizes_idx_1] = alpha1->
      data[sizes_idx_1];
  }

  emxFree_creal_T(&alpha1);

  //  return to diag
  //  test0 = sysMat*eigVec - eigVec*eigVal;
  //  output = max(max(real(test0)));
  //  fprintf('%e\n',output)
  //  end
  //  if(flag==2)
  //      error('Cannot compute subset of modes for eiganvalue solve when deployed, eigs not supported for compilition! Run with full set of modes.') 
  //  [eigVec,eigVal] = eigs(sysMat2,MMat2,828,'SM');  %subest of modes for eigenvalue solve 
  //  test = sysMat2*eigVec - MMat2*eigVec*eigVal;
  //  output = max(max(real(test)));
  //  fprintf('%e\n',output)
  //  end
  //  if(flag==3)
  //      sysMat=inv(M)*K;                      %eigenvalue solve on spring mass system only 
  //      [eigVec,eigVal] = eig(sysMat);
  //      [eigVec] = sortEigOutput(diag(eigVal),eigVec,numModes); %eigenvalues/vectors sorted in ascending frequency before returning 
  //  end
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
