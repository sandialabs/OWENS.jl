//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mtimes1.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// Arguments    : const double a[144]
//                const emxArray_real_T *b_d
//                const emxArray_int32_T *b_colidx
//                const emxArray_int32_T *b_rowidx
//                double c[144]
// Return Type  : void
//
void b_sparse_mtimes(const double a[144], const emxArray_real_T *b_d, const
                     emxArray_int32_T *b_colidx, const emxArray_int32_T
                     *b_rowidx, double c[144])
{
  int ccol;
  int coff;
  int bpend;
  int i;
  int bp;
  int aoff;
  double bd;
  int cidx;
  int aidx;
  std::memset(&c[0], 0, 144U * sizeof(double));
  if (b_colidx->data[b_colidx->size[0] - 1] - 1 != 0) {
    for (ccol = 0; ccol < 12; ccol++) {
      coff = ccol * 12 + 1;
      bpend = b_colidx->data[ccol + 1] - 1;
      i = b_colidx->data[ccol];
      for (bp = i; bp <= bpend; bp++) {
        aoff = (b_rowidx->data[bp - 1] - 1) * 12;
        bd = b_d->data[bp - 1];
        c[coff - 1] += a[aoff] * bd;
        c[coff] += a[aoff + 1] * bd;
        c[coff + 1] += a[aoff + 2] * bd;
        c[coff + 2] += a[aoff + 3] * bd;
        cidx = coff + 4;
        aidx = aoff + 4;
        c[cidx - 1] += a[aoff + 4] * bd;
        c[cidx] += a[aidx + 1] * bd;
        c[cidx + 1] += a[aidx + 2] * bd;
        c[cidx + 2] += a[aidx + 3] * bd;
        cidx = coff + 8;
        aidx = aoff + 8;
        c[cidx - 1] += a[aoff + 8] * bd;
        c[cidx] += a[aidx + 1] * bd;
        c[cidx + 1] += a[aidx + 2] * bd;
        c[cidx + 2] += a[aidx + 3] * bd;
      }
    }
  }
}

//
// Arguments    : const emxArray_real_T *a_d
//                const emxArray_int32_T *a_colidx
//                const emxArray_int32_T *a_rowidx
//                const double b_data[]
//                double c[12]
// Return Type  : void
//
void c_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
                     *a_colidx, const emxArray_int32_T *a_rowidx, const double
                     b_data[], double c[12])
{
  int acol;
  double bc;
  int i;
  int apend;
  int nap;
  int ap;
  int c_tmp;
  std::memset(&c[0], 0, 12U * sizeof(double));
  if (a_colidx->data[a_colidx->size[0] - 1] - 1 != 0) {
    for (acol = 0; acol < 12; acol++) {
      bc = b_data[acol];
      i = a_colidx->data[acol];
      apend = a_colidx->data[acol + 1];
      nap = apend - a_colidx->data[acol];
      if (nap >= 4) {
        apend = (apend - nap) + ((nap / 4) << 2);
        for (ap = i; ap <= apend - 1; ap += 4) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c[c_tmp] += a_d->data[ap - 1] * bc;
          c[a_rowidx->data[ap] - 1] += a_d->data[ap] * bc;
          c_tmp = a_rowidx->data[ap + 1] - 1;
          c[c_tmp] += a_d->data[ap + 1] * bc;
          c_tmp = a_rowidx->data[ap + 2] - 1;
          c[c_tmp] += a_d->data[ap + 2] * bc;
        }

        nap = a_colidx->data[acol + 1] - 1;
        for (ap = apend; ap <= nap; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c[c_tmp] += a_d->data[ap - 1] * bc;
        }
      } else {
        apend--;
        for (ap = i; ap <= apend; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c[c_tmp] += a_d->data[ap - 1] * bc;
        }
      }
    }
  }
}

//
// Arguments    : const emxArray_real_T *a_d
//                const emxArray_int32_T *a_colidx
//                const emxArray_int32_T *a_rowidx
//                const double b[144]
//                double c[144]
// Return Type  : void
//
void sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, const double b[144], double
                   c[144])
{
  int j;
  int coff;
  int acol;
  double bc;
  int i;
  int apend;
  int nap;
  int ap;
  int c_tmp;
  std::memset(&c[0], 0, 144U * sizeof(double));
  if (a_colidx->data[a_colidx->size[0] - 1] - 1 != 0) {
    for (j = 0; j < 12; j++) {
      coff = j * 12 - 1;
      for (acol = 0; acol < 12; acol++) {
        bc = b[acol + 12 * j];
        i = a_colidx->data[acol];
        apend = a_colidx->data[acol + 1];
        nap = apend - a_colidx->data[acol];
        if (nap >= 4) {
          apend = (apend - nap) + ((nap / 4) << 2);
          for (ap = i; ap <= apend - 1; ap += 4) {
            c_tmp = a_rowidx->data[ap - 1] + coff;
            c[c_tmp] += a_d->data[ap - 1] * bc;
            c_tmp = a_rowidx->data[ap] + coff;
            c[c_tmp] += a_d->data[ap] * bc;
            c_tmp = a_rowidx->data[ap + 1] + coff;
            c[c_tmp] += a_d->data[ap + 1] * bc;
            c_tmp = a_rowidx->data[ap + 2] + coff;
            c[c_tmp] += a_d->data[ap + 2] * bc;
          }

          nap = a_colidx->data[acol + 1] - 1;
          for (ap = apend; ap <= nap; ap++) {
            c_tmp = (a_rowidx->data[ap - 1] + 12 * j) - 1;
            c[c_tmp] += a_d->data[ap - 1] * bc;
          }
        } else {
          apend--;
          for (ap = i; ap <= apend; ap++) {
            c_tmp = a_rowidx->data[ap - 1] + coff;
            c[c_tmp] += a_d->data[ap - 1] * bc;
          }
        }
      }
    }
  }
}

//
// File trailer for mtimes1.cpp
//
// [EOF]
//
