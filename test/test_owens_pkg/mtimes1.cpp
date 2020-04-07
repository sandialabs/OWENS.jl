//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mtimes1.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "mtimes1.h"
#include "introsort.h"
#include "rt_nonfinite.h"
#include "sparse1.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_real_T *a_d
//                const emxArray_int32_T *a_colidx
//                const emxArray_int32_T *a_rowidx
//                const double b[144]
//                double c[144]
// Return Type  : void
//
void b_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
                     *a_colidx, const emxArray_int32_T *a_rowidx, const double
                     b[144], double c[144])
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
// Arguments    : const double a[144]
//                const emxArray_real_T *b_d
//                const emxArray_int32_T *b_colidx
//                const emxArray_int32_T *b_rowidx
//                double c[144]
// Return Type  : void
//
void c_sparse_mtimes(const double a[144], const emxArray_real_T *b_d, const
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
void d_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
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
//                int a_m
//                const emxArray_real_T *b_d
//                const emxArray_int32_T *b_colidx
//                const emxArray_int32_T *b_rowidx
//                int b_n
//                emxArray_real_T *c_d
//                emxArray_int32_T *c_colidx
//                emxArray_int32_T *c_rowidx
//                int *c_m
//                int *c_n
// Return Type  : void
//
void e_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
                     *a_colidx, const emxArray_int32_T *a_rowidx, int a_m, const
                     emxArray_real_T *b_d, const emxArray_int32_T *b_colidx,
                     const emxArray_int32_T *b_rowidx, int b_n, emxArray_real_T *
                     c_d, emxArray_int32_T *c_colidx, emxArray_int32_T *c_rowidx,
                     int *c_m, int *c_n)
{
  emxArray_int32_T *ccolidx;
  int i;
  int blen;
  emxArray_int32_T *flag;
  int cnnz;
  int j;
  int exitg1;
  int bcidx;
  int cstart;
  int cmax;
  int aend;
  emxArray_real_T *wd;
  int i1;
  int pb;
  boolean_T needSort;
  int pcstart;
  double bd;
  emxInit_int32_T(&ccolidx, 1);
  i = ccolidx->size[0];
  ccolidx->size[0] = b_colidx->size[0];
  emxEnsureCapacity_int32_T(ccolidx, i);
  blen = b_colidx->size[0];
  for (i = 0; i < blen; i++) {
    ccolidx->data[i] = 0;
  }

  emxInit_int32_T(&flag, 1);
  i = flag->size[0];
  flag->size[0] = a_m;
  emxEnsureCapacity_int32_T(flag, i);
  for (i = 0; i < a_m; i++) {
    flag->data[i] = 0;
  }

  cnnz = 0;
  j = 0;
  do {
    exitg1 = 0;
    if (j <= b_n - 1) {
      bcidx = b_colidx->data[j];
      cstart = cnnz;
      cmax = cnnz + a_m;
      ccolidx->data[j] = cnnz + 1;
      while ((bcidx < b_colidx->data[j + 1]) && (cnnz <= cmax)) {
        blen = b_rowidx->data[bcidx - 1];
        aend = a_colidx->data[blen] - 1;
        i = a_colidx->data[blen - 1];
        for (blen = i; blen <= aend; blen++) {
          i1 = a_rowidx->data[blen - 1] - 1;
          if (flag->data[i1] != j + 1) {
            flag->data[i1] = j + 1;
            cnnz++;
          }
        }

        bcidx++;
      }

      if (cnnz < cstart) {
        exitg1 = 1;
      } else {
        j++;
      }
    } else {
      ccolidx->data[b_n] = cnnz + 1;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  sparse_sparse(a_m, b_n, cnnz, c_d, c_colidx, c_rowidx, c_m, c_n);
  i = c_colidx->size[0];
  c_colidx->size[0] = ccolidx->size[0];
  emxEnsureCapacity_int32_T(c_colidx, i);
  blen = ccolidx->size[0];
  for (i = 0; i < blen; i++) {
    c_colidx->data[i] = ccolidx->data[i];
  }

  emxInit_real_T(&wd, 1);
  i = wd->size[0];
  wd->size[0] = a_m;
  emxEnsureCapacity_real_T(wd, i);
  i = flag->size[0];
  flag->size[0] = a_m;
  emxEnsureCapacity_int32_T(flag, i);
  for (i = 0; i < a_m; i++) {
    flag->data[i] = 0;
  }

  pb = 0;
  cnnz = -1;
  for (j = 0; j < b_n; j++) {
    aend = j + 1;
    needSort = false;
    pcstart = cnnz + 2;
    blen = (b_colidx->data[aend] - pb) - 1;
    if (blen != 0) {
      if (blen == 1) {
        cstart = a_colidx->data[b_rowidx->data[pb]] - 1;
        i = a_colidx->data[b_rowidx->data[pb] - 1];
        for (cmax = i; cmax <= cstart; cmax++) {
          cnnz++;
          blen = a_rowidx->data[cmax - 1];
          c_rowidx->data[cnnz] = blen;
          wd->data[blen - 1] = a_d->data[cmax - 1] * b_d->data[pb];
        }

        pb++;
      } else {
        cstart = a_colidx->data[b_rowidx->data[pb]] - 1;
        i = a_colidx->data[b_rowidx->data[pb] - 1];
        for (cmax = i; cmax <= cstart; cmax++) {
          cnnz++;
          blen = a_rowidx->data[cmax - 1];
          bcidx = blen - 1;
          flag->data[bcidx] = cnnz + 1;
          c_rowidx->data[cnnz] = blen;
          wd->data[bcidx] = a_d->data[cmax - 1] * b_d->data[pb];
        }

        for (pb++; pb + 1 < b_colidx->data[aend]; pb++) {
          bd = b_d->data[pb];
          cstart = a_colidx->data[b_rowidx->data[pb]] - 1;
          i = a_colidx->data[b_rowidx->data[pb] - 1];
          for (cmax = i; cmax <= cstart; cmax++) {
            i1 = a_rowidx->data[cmax - 1];
            blen = i1 - 1;
            if (flag->data[blen] < pcstart) {
              cnnz++;
              flag->data[blen] = cnnz + 1;
              c_rowidx->data[cnnz] = i1;
              wd->data[blen] = a_d->data[cmax - 1] * bd;
              needSort = true;
            } else {
              wd->data[blen] += a_d->data[cmax - 1] * bd;
            }
          }
        }
      }
    }

    cmax = ccolidx->data[aend] - 1;
    pcstart = ccolidx->data[aend - 1];
    if (needSort) {
      introsort(c_rowidx, ccolidx->data[aend - 1], ccolidx->data[aend] - 1);
    }

    for (cstart = pcstart; cstart <= cmax; cstart++) {
      c_d->data[cstart - 1] = wd->data[c_rowidx->data[cstart - 1] - 1];
    }
  }

  emxFree_int32_T(&flag);
  emxFree_real_T(&wd);
  blen = 1;
  i = ccolidx->size[0];
  emxFree_int32_T(&ccolidx);
  for (cstart = 0; cstart <= i - 2; cstart++) {
    cmax = c_colidx->data[cstart];
    c_colidx->data[cstart] = blen;
    while (cmax < c_colidx->data[cstart + 1]) {
      bcidx = c_rowidx->data[cmax - 1];
      bd = c_d->data[cmax - 1];
      cmax++;
      if (bd != 0.0) {
        c_d->data[blen - 1] = bd;
        c_rowidx->data[blen - 1] = bcidx;
        blen++;
      }
    }
  }

  c_colidx->data[c_colidx->size[0] - 1] = blen;
}

//
// Arguments    : const emxArray_real_T *a_d
//                const emxArray_int32_T *a_colidx
//                const emxArray_int32_T *a_rowidx
//                const emxArray_real_T *b_d
//                const emxArray_int32_T *b_colidx
//                const emxArray_int32_T *b_rowidx
//                emxArray_real_T *c_d
//                emxArray_int32_T *c_colidx
//                emxArray_int32_T *c_rowidx
// Return Type  : void
//
void sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d,
                   const emxArray_int32_T *b_colidx, const emxArray_int32_T
                   *b_rowidx, emxArray_real_T *c_d, emxArray_int32_T *c_colidx,
                   emxArray_int32_T *c_rowidx)
{
  int i;
  int b_i;
  int cnnz;
  int flag[12];
  int j;
  int exitg1;
  int bcidx;
  int cstart;
  int cmax;
  int paend;
  int i1;
  int pb;
  boolean_T needSort;
  int pa;
  double bd;
  double wd[12];
  i = c_colidx->size[0];
  c_colidx->size[0] = b_colidx->size[0];
  emxEnsureCapacity_int32_T(c_colidx, i);
  b_i = b_colidx->size[0];
  for (i = 0; i < b_i; i++) {
    c_colidx->data[i] = 0;
  }

  for (b_i = 0; b_i < 12; b_i++) {
    flag[b_i] = 0;
  }

  cnnz = 0;
  j = 0;
  do {
    exitg1 = 0;
    if (j < 12) {
      bcidx = b_colidx->data[j];
      cstart = cnnz;
      cmax = cnnz + 12;
      c_colidx->data[j] = cnnz + 1;
      while ((bcidx < b_colidx->data[j + 1]) && (cnnz <= cmax)) {
        b_i = b_rowidx->data[bcidx - 1];
        paend = a_colidx->data[b_i] - 1;
        i = a_colidx->data[b_i - 1];
        for (b_i = i; b_i <= paend; b_i++) {
          i1 = a_rowidx->data[b_i - 1] - 1;
          if (flag[i1] != j + 1) {
            flag[i1] = j + 1;
            cnnz++;
          }
        }

        bcidx++;
      }

      if (cnnz < cstart) {
        exitg1 = 1;
      } else {
        j++;
      }
    } else {
      c_colidx->data[12] = cnnz + 1;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (cnnz >= 1) {
    b_i = cnnz;
  } else {
    b_i = 1;
  }

  i = c_d->size[0];
  c_d->size[0] = b_i;
  emxEnsureCapacity_real_T(c_d, i);
  for (i = 0; i < b_i; i++) {
    c_d->data[i] = 0.0;
  }

  i = c_rowidx->size[0];
  c_rowidx->size[0] = b_i;
  emxEnsureCapacity_int32_T(c_rowidx, i);
  for (i = 0; i < b_i; i++) {
    c_rowidx->data[i] = 0;
  }

  for (b_i = 0; b_i < 12; b_i++) {
    flag[b_i] = 0;
  }

  pb = 0;
  cnnz = -1;
  for (j = 0; j < 12; j++) {
    needSort = false;
    cstart = cnnz + 2;
    cmax = b_colidx->data[j + 1];
    b_i = (cmax - pb) - 1;
    if (b_i != 0) {
      if (b_i == 1) {
        paend = a_colidx->data[b_rowidx->data[pb]] - 1;
        i = a_colidx->data[b_rowidx->data[pb] - 1];
        for (pa = i; pa <= paend; pa++) {
          cnnz++;
          b_i = a_rowidx->data[pa - 1];
          c_rowidx->data[cnnz] = b_i;
          wd[b_i - 1] = a_d->data[pa - 1] * b_d->data[pb];
        }

        pb++;
      } else {
        paend = a_colidx->data[b_rowidx->data[pb]] - 1;
        i = a_colidx->data[b_rowidx->data[pb] - 1];
        for (pa = i; pa <= paend; pa++) {
          cnnz++;
          b_i = a_rowidx->data[pa - 1];
          bcidx = b_i - 1;
          flag[bcidx] = cnnz + 1;
          c_rowidx->data[cnnz] = b_i;
          wd[bcidx] = a_d->data[pa - 1] * b_d->data[pb];
        }

        for (pb++; pb + 1 < cmax; pb++) {
          bd = b_d->data[pb];
          paend = a_colidx->data[b_rowidx->data[pb]] - 1;
          i = a_colidx->data[b_rowidx->data[pb] - 1];
          for (pa = i; pa <= paend; pa++) {
            i1 = a_rowidx->data[pa - 1];
            b_i = i1 - 1;
            if (flag[b_i] < cstart) {
              cnnz++;
              flag[b_i] = cnnz + 1;
              c_rowidx->data[cnnz] = i1;
              wd[b_i] = a_d->data[pa - 1] * bd;
              needSort = true;
            } else {
              wd[b_i] += a_d->data[pa - 1] * bd;
            }
          }
        }
      }
    }

    b_i = c_colidx->data[j + 1] - 1;
    i = c_colidx->data[j];
    if (needSort) {
      introsort(c_rowidx, c_colidx->data[j], c_colidx->data[j + 1] - 1);
    }

    for (cmax = i; cmax <= b_i; cmax++) {
      c_d->data[cmax - 1] = wd[c_rowidx->data[cmax - 1] - 1];
    }
  }

  b_i = 1;
  i = c_colidx->size[0];
  for (cmax = 0; cmax <= i - 2; cmax++) {
    bcidx = c_colidx->data[cmax];
    c_colidx->data[cmax] = b_i;
    while (bcidx < c_colidx->data[cmax + 1]) {
      cstart = c_rowidx->data[bcidx - 1];
      bd = c_d->data[bcidx - 1];
      bcidx++;
      if (bd != 0.0) {
        c_d->data[b_i - 1] = bd;
        c_rowidx->data[b_i - 1] = cstart;
        b_i++;
      }
    }
  }

  c_colidx->data[c_colidx->size[0] - 1] = b_i;
}

//
// File trailer for mtimes1.cpp
//
// [EOF]
//
