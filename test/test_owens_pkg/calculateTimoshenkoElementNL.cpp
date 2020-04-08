//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateTimoshenkoElementNL.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//

// Include Files
#include "calculateTimoshenkoElementNL.h"
#include "calculateLambda.h"
#include "calculateShapeFunctions.h"
#include "fillIn.h"
#include "find.h"
#include "introsort.h"
#include "mtimes.h"
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "sparse1.h"
#include "strcmp.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Declarations
static void b_mapMatrixNonSym(const double Ktemp_data[], double Kel[144]);
static void mapMatrixNonSym(const double Ktemp[144], double Kel[144]);

// Function Definitions

//
// ----- function to form total stifness matrix and transform to desired
//  DOF mapping
// Arguments    : const double Ktemp_data[]
//                double Kel[144]
// Return Type  : void
//
static void b_mapMatrixNonSym(const double Ktemp_data[], double Kel[144])
{
  emxArray_real_T *y_d;
  emxArray_int32_T *y_colidx;
  emxArray_int8_T *y_rowidx;
  int numalloc;
  int cstart;
  int i;
  int bcidx;
  emxArray_real_T *a_d;
  emxArray_int32_T *a_colidx;
  emxArray_int32_T *a_rowidx;
  emxArray_int32_T *ccolidx;
  int cnnz;
  int flag[12];
  int j;
  int exitg1;
  int cmax;
  coder_internal_sparse c;
  int pcstart;
  int i1;
  int pb;
  boolean_T needSort;
  coder_internal_sparse expl_temp;
  int blen_tmp;
  double wd[12];
  double bd;
  emxInit_real_T(&y_d, 1);
  emxInit_int32_T(&y_colidx, 1);
  emxInit_int8_T(&y_rowidx, 1);

  // map to FEA numbering
  numalloc = 0;
  for (cstart = 0; cstart < 144; cstart++) {
    if (Ktemp_data[cstart] != 0.0) {
      numalloc++;
    }
  }

  if (numalloc < 1) {
    numalloc = 1;
  }

  i = y_d->size[0];
  y_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(y_d, i);
  for (i = 0; i < numalloc; i++) {
    y_d->data[i] = 0.0;
  }

  i = y_colidx->size[0];
  y_colidx->size[0] = 13;
  emxEnsureCapacity_int32_T(y_colidx, i);
  for (i = 0; i < 13; i++) {
    y_colidx->data[i] = 0;
  }

  y_colidx->data[0] = 1;
  i = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int8_T(y_rowidx, i);
  for (i = 0; i < numalloc; i++) {
    y_rowidx->data[i] = 0;
  }

  y_rowidx->data[0] = 1;
  numalloc = 0;
  for (bcidx = 0; bcidx < 12; bcidx++) {
    for (cstart = 0; cstart < 12; cstart++) {
      i = cstart + 12 * bcidx;
      if (Ktemp_data[i] != 0.0) {
        y_rowidx->data[numalloc] = static_cast<signed char>((cstart + 1));
        y_d->data[numalloc] = Ktemp_data[i];
        numalloc++;
      }
    }

    y_colidx->data[bcidx + 1] = numalloc + 1;
  }

  emxInit_real_T(&a_d, 1);
  emxInit_int32_T(&a_colidx, 1);
  emxInit_int32_T(&a_rowidx, 1);
  emxInit_int32_T(&ccolidx, 1);
  b_sparse(a_d, a_colidx, a_rowidx);
  i = ccolidx->size[0];
  ccolidx->size[0] = 13;
  emxEnsureCapacity_int32_T(ccolidx, i);
  for (i = 0; i < 13; i++) {
    ccolidx->data[i] = 0;
  }

  for (numalloc = 0; numalloc < 12; numalloc++) {
    flag[numalloc] = 0;
  }

  cnnz = 0;
  j = 0;
  do {
    exitg1 = 0;
    if (j < 12) {
      bcidx = y_colidx->data[j];
      cstart = cnnz;
      cmax = cnnz + 12;
      ccolidx->data[j] = cnnz + 1;
      while ((bcidx < y_colidx->data[j + 1]) && (cnnz <= cmax)) {
        numalloc = y_rowidx->data[bcidx - 1];
        pcstart = a_colidx->data[numalloc] - 1;
        i = a_colidx->data[numalloc - 1];
        for (numalloc = i; numalloc <= pcstart; numalloc++) {
          i1 = a_rowidx->data[numalloc - 1] - 1;
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
      ccolidx->data[12] = cnnz + 1;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  c_emxInitStruct_coder_internal_(&c);
  sparse_spallocLike(cnnz, c.d, c.colidx, c.rowidx);
  i = c.colidx->size[0];
  c.colidx->size[0] = ccolidx->size[0];
  emxEnsureCapacity_int32_T(c.colidx, i);
  numalloc = ccolidx->size[0];
  for (i = 0; i < numalloc; i++) {
    c.colidx->data[i] = ccolidx->data[i];
  }

  for (numalloc = 0; numalloc < 12; numalloc++) {
    flag[numalloc] = 0;
  }

  pb = 0;
  cnnz = -1;
  for (j = 0; j < 12; j++) {
    needSort = false;
    pcstart = cnnz + 2;
    blen_tmp = y_colidx->data[j + 1];
    numalloc = (blen_tmp - pb) - 1;
    if (numalloc != 0) {
      if (numalloc == 1) {
        cstart = a_colidx->data[y_rowidx->data[pb]] - 1;
        i = a_colidx->data[y_rowidx->data[pb] - 1];
        for (cmax = i; cmax <= cstart; cmax++) {
          cnnz++;
          i1 = a_rowidx->data[cmax - 1];
          c.rowidx->data[cnnz] = i1;
          wd[i1 - 1] = a_d->data[cmax - 1] * y_d->data[pb];
        }

        pb++;
      } else {
        cstart = a_colidx->data[y_rowidx->data[pb]] - 1;
        i = a_colidx->data[y_rowidx->data[pb] - 1];
        for (cmax = i; cmax <= cstart; cmax++) {
          cnnz++;
          numalloc = a_rowidx->data[cmax - 1];
          bcidx = numalloc - 1;
          flag[bcidx] = cnnz + 1;
          c.rowidx->data[cnnz] = numalloc;
          wd[bcidx] = a_d->data[cmax - 1] * y_d->data[pb];
        }

        for (pb++; pb + 1 < blen_tmp; pb++) {
          bd = y_d->data[pb];
          cstart = a_colidx->data[y_rowidx->data[pb]] - 1;
          i = a_colidx->data[y_rowidx->data[pb] - 1];
          for (cmax = i; cmax <= cstart; cmax++) {
            i1 = a_rowidx->data[cmax - 1];
            numalloc = i1 - 1;
            if (flag[numalloc] < pcstart) {
              cnnz++;
              flag[numalloc] = cnnz + 1;
              c.rowidx->data[cnnz] = i1;
              wd[numalloc] = a_d->data[cmax - 1] * bd;
              needSort = true;
            } else {
              wd[numalloc] += a_d->data[cmax - 1] * bd;
            }
          }
        }
      }
    }

    numalloc = ccolidx->data[j + 1] - 1;
    i = ccolidx->data[j];
    if (needSort) {
      introsort(c.rowidx, ccolidx->data[j], ccolidx->data[j + 1] - 1);
    }

    for (cstart = i; cstart <= numalloc; cstart++) {
      c.d->data[cstart - 1] = wd[c.rowidx->data[cstart - 1] - 1];
    }
  }

  emxFree_int32_T(&ccolidx);
  emxFree_int8_T(&y_rowidx);
  c_emxInitStruct_coder_internal_(&expl_temp);
  sparse_fillIn(&c);
  c_sparse(a_d, a_colidx, a_rowidx);
  sparse_mtimes(c.d, c.colidx, c.rowidx, a_d, a_colidx, a_rowidx, &expl_temp);
  i = y_d->size[0];
  y_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(y_d, i);
  numalloc = expl_temp.d->size[0];
  emxFree_int32_T(&a_rowidx);
  emxFree_real_T(&a_d);
  c_emxFreeStruct_coder_internal_(&c);
  for (i = 0; i < numalloc; i++) {
    y_d->data[i] = expl_temp.d->data[i];
  }

  i = y_colidx->size[0];
  y_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(y_colidx, i);
  numalloc = expl_temp.colidx->size[0];
  for (i = 0; i < numalloc; i++) {
    y_colidx->data[i] = expl_temp.colidx->data[i];
  }

  i = a_colidx->size[0];
  a_colidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(a_colidx, i);
  numalloc = expl_temp.rowidx->size[0];
  for (i = 0; i < numalloc; i++) {
    a_colidx->data[i] = expl_temp.rowidx->data[i];
  }

  c_emxFreeStruct_coder_internal_(&expl_temp);
  std::memset(&Kel[0], 0, 144U * sizeof(double));
  for (numalloc = 0; numalloc < 12; numalloc++) {
    bcidx = y_colidx->data[numalloc + 1] - 1;
    i = y_colidx->data[numalloc];
    for (cstart = i; cstart <= bcidx; cstart++) {
      Kel[(a_colidx->data[cstart - 1] + 12 * numalloc) - 1] = y_d->data[cstart -
        1];
    }
  }

  emxFree_int32_T(&a_colidx);
  emxFree_int32_T(&y_colidx);
  emxFree_real_T(&y_d);

  // declare map
  //  map = [1, 7, 2, 8, 3, 9,...
  //        4, 10, 5, 11, 6, 12];
  //
  //  %map to FEA numbering
  //  for i=1:a
  //      I=map(i);
  //      for j=1:a
  //          J=map(j);
  //          Kel(I,J) = Ktemp(i,j);
  //      end
  //  end
}

//
// ----- function to form total stifness matrix and transform to desired
//  DOF mapping
// Arguments    : const double Ktemp[144]
//                double Kel[144]
// Return Type  : void
//
static void mapMatrixNonSym(const double Ktemp[144], double Kel[144])
{
  emxArray_real_T *this_d;
  emxArray_int32_T *this_colidx;
  emxArray_int32_T *this_rowidx;
  emxArray_real_T *t1_d;
  emxArray_int32_T *t1_colidx;
  emxArray_int32_T *t1_rowidx;
  coder_internal_sparse expl_temp;
  int i;
  int loop_ub;
  int cend;
  int idx;
  emxInit_real_T(&this_d, 1);
  emxInit_int32_T(&this_colidx, 1);
  emxInit_int32_T(&this_rowidx, 1);
  emxInit_real_T(&t1_d, 1);
  emxInit_int32_T(&t1_colidx, 1);
  emxInit_int32_T(&t1_rowidx, 1);
  c_emxInitStruct_coder_internal_(&expl_temp);

  // map to FEA numbering
  b_sparse(t1_d, t1_colidx, t1_rowidx);
  sparse(Ktemp, this_d, this_colidx, this_rowidx);
  sparse_mtimes(t1_d, t1_colidx, t1_rowidx, this_d, this_colidx, this_rowidx,
                &expl_temp);
  i = this_d->size[0];
  this_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(this_d, i);
  loop_ub = expl_temp.d->size[0];
  for (i = 0; i < loop_ub; i++) {
    this_d->data[i] = expl_temp.d->data[i];
  }

  i = this_colidx->size[0];
  this_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(this_colidx, i);
  loop_ub = expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    this_colidx->data[i] = expl_temp.colidx->data[i];
  }

  i = this_rowidx->size[0];
  this_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(this_rowidx, i);
  loop_ub = expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    this_rowidx->data[i] = expl_temp.rowidx->data[i];
  }

  c_sparse(t1_d, t1_colidx, t1_rowidx);
  sparse_mtimes(this_d, this_colidx, this_rowidx, t1_d, t1_colidx, t1_rowidx,
                &expl_temp);
  i = this_d->size[0];
  this_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(this_d, i);
  loop_ub = expl_temp.d->size[0];
  emxFree_int32_T(&t1_rowidx);
  emxFree_int32_T(&t1_colidx);
  emxFree_real_T(&t1_d);
  for (i = 0; i < loop_ub; i++) {
    this_d->data[i] = expl_temp.d->data[i];
  }

  i = this_colidx->size[0];
  this_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(this_colidx, i);
  loop_ub = expl_temp.colidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    this_colidx->data[i] = expl_temp.colidx->data[i];
  }

  i = this_rowidx->size[0];
  this_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(this_rowidx, i);
  loop_ub = expl_temp.rowidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    this_rowidx->data[i] = expl_temp.rowidx->data[i];
  }

  c_emxFreeStruct_coder_internal_(&expl_temp);
  std::memset(&Kel[0], 0, 144U * sizeof(double));
  for (loop_ub = 0; loop_ub < 12; loop_ub++) {
    cend = this_colidx->data[loop_ub + 1] - 1;
    i = this_colidx->data[loop_ub];
    for (idx = i; idx <= cend; idx++) {
      Kel[(this_rowidx->data[idx - 1] + 12 * loop_ub) - 1] = this_d->data[idx -
        1];
    }
  }

  emxFree_int32_T(&this_rowidx);
  emxFree_int32_T(&this_colidx);
  emxFree_real_T(&this_d);

  // declare map
  //  map = [1, 7, 2, 8, 3, 9,...
  //        4, 10, 5, 11, 6, 12];
  //
  //  %map to FEA numbering
  //  for i=1:a
  //      I=map(i);
  //      for j=1:a
  //          J=map(j);
  //          Kel(I,J) = Ktemp(i,j);
  //      end
  //  end
}

//
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
// Arguments    : const p_struct_T *input
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void b_calculateTimoshenkoElementNL(const p_struct_T *input, const f_struct_T
  *elStorage, o_struct_T *output)
{
  double disp_iter_data[12];
  double dispm1_data[12];
  int dispdot_size_idx_1;
  double dispdot_data[12];
  int dispddot_size_idx_1;
  double dispddot_data[12];
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double Omega;
  double OmegaDot;
  double lambda[144];
  int i;
  double O1;
  double Oel[3];
  double c_tmp;
  double O2;
  double H46_idx_3;
  double O3;
  double O1dot;
  double ODotel[3];
  double a_temp[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double S26_idx_0;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double integrationFactor;
  double b_c_tmp;
  double c_c_tmp;
  double S13_idx_0;
  double S35_idx_0;
  double S36_idx_0;
  double ycm;
  double S34_idx_0;
  double S45_idx_0;
  double S46_idx_0;
  double C12_idx_0;
  double K15_idx_1;
  double C13_idx_0;
  double K15_idx_0;
  double rhoA;
  double C23_idx_0;
  double K16_idx_0;
  double zcm;
  double C24_idx_0;
  double K56_idx_0;
  double C25_idx_0;
  double f_tmp;
  double C26_idx_0;
  double disMomentgp[3];
  double b_f_tmp;
  double C34_idx_0;
  double c_f_tmp;
  double C35_idx_0;
  double d_f_tmp;
  double C36_idx_0;
  double e_f_tmp;
  double posLocal[3];
  double f;
  double disLoadgpLocal[3];
  double C14_idx_0;
  double b_f;
  double c_f;
  double C45_idx_0;
  double C46_idx_0;
  double H13_idx_0;
  double H23_idx_0;
  double H24_idx_0;
  double H25_idx_0;
  double H26_idx_0;
  double H34_idx_0;
  double H35_idx_0;
  double H36_idx_0;
  double H14_idx_0;
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double S35_idx_1;
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C12_idx_1;
  double C13_idx_1;
  double C23_idx_1;
  double C24_idx_1;
  double C25_idx_1;
  double C26_idx_1;
  double C34_idx_1;
  double C35_idx_1;
  double C36_idx_1;
  double C14_idx_1;
  double C45_idx_1;
  double C46_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double H23_idx_1;
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C12_idx_2;
  double C13_idx_2;
  double C23_idx_2;
  double C24_idx_2;
  double C25_idx_2;
  double C26_idx_2;
  double C34_idx_2;
  double C35_idx_2;
  double C36_idx_2;
  double C14_idx_2;
  double C45_idx_2;
  double C46_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C12_idx_3;
  double C13_idx_3;
  double C23_idx_3;
  double C24_idx_3;
  double C25_idx_3;
  double C26_idx_3;
  double C34_idx_3;
  double C35_idx_3;
  double C36_idx_3;
  double C14_idx_3;
  double C45_idx_3;
  double C46_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H34_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double b_elStorage[144];
  double Kenr[144];
  double K16_idx_1;
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double Khate[144];
  double reshapes_f3_data[24];
  double reshapes_f4_data[24];
  double reshapes_f1[24];
  double reshapes_f2[24];
  double reshapes_f5[24];
  double reshapes_f6[24];
  double Ce[144];
  int reshapes_f1_data_tmp;
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fhate[12];
  double Fe[12];
  int tmp_size[1];
  boolean_T concMassFlag;
  double FhatLessConc[12];
  double b_Khate[12];
  double b_data[12];

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  std::memcpy(&disp_iter_data[0], &input->displ_iter.data[0], 12U * sizeof
              (double));
  dispm1_data[0] = 0.0;

  // declare type
  dispdot_size_idx_1 = 1;
  dispdot_data[0] = 0.0;

  // declare type
  dispddot_size_idx_1 = 1;
  dispddot_data[0] = 0.0;

  // declare type
  // options for Dean integrator
  if (g_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&dispm1_data[0], &input->dispm1.data[0], 12U * sizeof(double));
  } else {
    if (h_strcmp(input->analysisType.data, input->analysisType.size)) {
      // options for newmark beta integrator
      dispdot_size_idx_1 = 12;
      dispddot_size_idx_1 = 12;
      std::memcpy(&dispdot_data[0], &input->dispdot.data[0], 12U * sizeof(double));
      std::memcpy(&dispddot_data[0], &input->dispddot.data[0], 12U * sizeof
                  (double));
    }
  }

  // --------------------------------------------
  // setting for modal analysis flag
  if (i_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&disp_iter_data[0], &input->disp.data[0], 12U * sizeof(double));
  }

  // setting for initial reduced order model calculations
  if (j_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&disp_iter_data[0], &input->disp.data[0], 12U * sizeof(double));
  }

  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  Omega = 6.2831853071795862 * input->Omega;
  OmegaDot = 6.2831853071795862 * input->OmegaDot;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input->sweepAngle * 3.1415926535897931 / 180.0,
                  input->coneAngle * 3.1415926535897931 / 180.0,
                  (input->rollAngle + 0.5 * (input->sectionProps.twist[0] +
    input->sectionProps.twist[1])) * 3.1415926535897931 / 180.0, lambda);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    c_tmp = lambda[i + 24];
    H46_idx_3 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    Oel[i] = H46_idx_3 + c_tmp * Omega;
    a_temp[i] = (input->CN2H[i] * 0.0 + input->CN2H[i + 3] * 0.0) + input->
      CN2H[i + 6] * 9.81;
    ODotel[i] = H46_idx_3 + c_tmp * OmegaDot;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  S26_idx_0 = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input->xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &H46_idx_3);
    integrationFactor = H46_idx_3 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input->sectionProps.rhoA[0] + N_data[1] *
      input->sectionProps.rhoA[1];
    ycm = N_data[0] * input->sectionProps.ycm[0] + N_data[1] *
      input->sectionProps.ycm[1];
    zcm = N_data[0] * input->sectionProps.zcm[0] + N_data[1] *
      input->sectionProps.zcm[1];
    K15_idx_0 = N_data[0] * input->x.data[0] + N_data[1] * input->x.data[1];
    K16_idx_0 = N_data[0] * input->y.data[0] + N_data[1] * input->y.data[1];
    H46_idx_3 = N_data[0] * input->z.data[0] + N_data[1] * input->z.data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      c_tmp = lambda[i + 12];
      K15_idx_1 = lambda[i + 24];
      posLocal[i] = (lambda[i] * K15_idx_0 + c_tmp * K16_idx_0) + K15_idx_1 *
        H46_idx_3;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + c_tmp * disMomentgp[1])
        + K15_idx_1 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // calculate static aerodynamic load
    // distributed/body force load calculations
    f_tmp = (O2 * O2 + O3 * O3) * posLocal[0];
    b_f_tmp = O2dot * posLocal[2];
    c_f_tmp = O3dot * posLocal[1];
    d_f_tmp = S26_idx_0 * posLocal[2];
    e_f_tmp = O1 * O2 * posLocal[1];
    f = rhoA * ((((f_tmp - e_f_tmp) - d_f_tmp) + c_f_tmp) - b_f_tmp) -
      disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    H46_idx_3 = O1dot * posLocal[2];
    K56_idx_0 = O3dot * posLocal[0];
    b_f = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 * O3)
                    - posLocal[0] * O1 * O2) + H46_idx_3) - K56_idx_0) -
      disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = O2dot * posLocal[0];
    K16_idx_0 = O1dot * posLocal[1];
    c_f = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - S26_idx_0 * posLocal[0])
                    - O2 * O3 * posLocal[1]) + K15_idx_0) - K16_idx_0) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) -
      posLocal[1] * (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                          (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                         (K15_idx_0 - K16_idx_0)) - zcm * (H46_idx_3 - K56_idx_0))
      - disMomentgp[0];

    // This function is a general routine to calculate an element vector
    H46_idx_3 = rhoA * zcm * ((((f_tmp - posLocal[1] * O1 * O2) - posLocal[2] *
      O1 * O3) - b_f_tmp) + c_f_tmp) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    K56_idx_0 = rhoA * ycm * ((((d_f_tmp + e_f_tmp) - f_tmp) - c_f_tmp) +
      b_f_tmp) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    F1_data_idx_0 += f * N_data[0] * integrationFactor;
    F2_data_idx_0 += b_f * N_data[0] * integrationFactor;
    F3_data_idx_0 += c_f * N_data[0] * integrationFactor;
    F4_data_idx_0 += K15_idx_0 * N_data[0] * integrationFactor;
    F5_data_idx_0 += H46_idx_3 * N_data[0] * integrationFactor;
    F6_data_idx_0 += K56_idx_0 * N_data[0] * integrationFactor;
    F1_data_idx_1 += f * N_data[1] * integrationFactor;
    F2_data_idx_1 += b_f * N_data[1] * integrationFactor;
    F3_data_idx_1 += c_f * N_data[1] * integrationFactor;
    F4_data_idx_1 += K15_idx_0 * N_data[1] * integrationFactor;
    F5_data_idx_1 += H46_idx_3 * N_data[1] * integrationFactor;
    F6_data_idx_1 += K56_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  H46_idx_3 = Oel[0] * Oel[1];
  b_c_tmp = Oel[0] * Oel[0];
  c_c_tmp = b_c_tmp + Oel[2] * Oel[2];
  b_c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  integrationFactor = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  O1 = elStorage->S23[0] * Oel[1] * Oel[2];
  O2dot = elStorage->S25[0] * H46_idx_3;
  S26_idx_0 = elStorage->S26[0] * H46_idx_3;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  S36_idx_0 = elStorage->S36[0] * Oel[0] * Oel[2];
  ycm = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel[0] *
    Oel[1];
  O3 = elStorage->S24_1[0] * c_c_tmp + elStorage->S24_2[0] * Oel[1] * Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * b_c_tmp + elStorage->S34_2[0] * Oel[1] *
    Oel[2];
  S45_idx_0 = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel
    [0] * Oel[1];
  S46_idx_0 = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[0];
  C12_idx_0 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[0];
  C13_idx_0 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[0];
  C23_idx_0 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[0];
  C24_idx_0 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[0];
  C25_idx_0 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[0];
  C26_idx_0 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[0];
  C34_idx_0 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[0];
  C35_idx_0 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[0];
  C36_idx_0 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[0];
  f = elStorage->C14_2[0];
  C14_idx_0 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[0];
  c_f = elStorage->C45_2[0];
  C45_idx_0 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[0];
  rhoA = elStorage->C46_2[0];
  C46_idx_0 = zcm * Oel[1] + rhoA * Oel[2];
  O2 = 0.5 * c_tmp * ODotel[2];
  H13_idx_0 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_0 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_0 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_0 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_0 = 0.5 * f_tmp * ODotel[2];
  H34_idx_0 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_0 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_0 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_0 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_0 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_0 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * H46_idx_3;
  S26_idx_1 = elStorage->S26[1] * H46_idx_3;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * c_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * b_c_tmp + elStorage->S34_2[1] * Oel[1] *
    Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[1];
  C12_idx_1 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[1];
  C13_idx_1 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[1];
  C23_idx_1 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[1];
  C24_idx_1 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[1];
  C25_idx_1 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[1];
  C26_idx_1 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[1];
  C34_idx_1 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[1];
  C35_idx_1 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[1];
  C36_idx_1 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[1];
  f = elStorage->C14_2[1];
  C14_idx_1 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[1];
  c_f = elStorage->C45_2[1];
  C45_idx_1 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[1];
  rhoA = elStorage->C46_2[1];
  C46_idx_1 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_1 = 0.5 * c_tmp * ODotel[2];
  H13_idx_1 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_1 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_1 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_1 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_1 = 0.5 * f_tmp * ODotel[2];
  H34_idx_1 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_1 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_1 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_1 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_1 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_1 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * H46_idx_3;
  S26_idx_2 = elStorage->S26[2] * H46_idx_3;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * c_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * b_c_tmp + elStorage->S34_2[2] * Oel[1] *
    Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[2];
  C12_idx_2 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[2];
  C13_idx_2 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[2];
  C23_idx_2 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[2];
  C24_idx_2 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[2];
  C25_idx_2 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[2];
  C26_idx_2 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[2];
  C34_idx_2 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[2];
  C35_idx_2 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[2];
  C36_idx_2 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[2];
  f = elStorage->C14_2[2];
  C14_idx_2 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[2];
  c_f = elStorage->C45_2[2];
  C45_idx_2 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[2];
  rhoA = elStorage->C46_2[2];
  C46_idx_2 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_2 = 0.5 * c_tmp * ODotel[2];
  H13_idx_2 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_2 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_2 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_2 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_2 = 0.5 * f_tmp * ODotel[2];
  H34_idx_2 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_2 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_2 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_2 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_2 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_2 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * H46_idx_3;
  S26_idx_3 = elStorage->S26[3] * H46_idx_3;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * c_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * b_c_tmp + elStorage->S34_2[3] * Oel[1] *
    Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[3];
  C12_idx_3 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[3];
  C13_idx_3 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[3];
  C23_idx_3 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[3];
  C24_idx_3 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[3];
  C25_idx_3 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[3];
  C26_idx_3 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[3];
  C34_idx_3 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[3];
  C35_idx_3 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[3];
  C36_idx_3 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[3];
  f = elStorage->C14_2[3];
  C14_idx_3 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[3];
  c_f = elStorage->C45_2[3];
  C45_idx_3 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[3];
  rhoA = elStorage->C46_2[3];
  C46_idx_3 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_3 = 0.5 * c_tmp * ODotel[2];
  H13_idx_3 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_3 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_3 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_3 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_3 = 0.5 * f_tmp * ODotel[2];
  H34_idx_3 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_3 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_3 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_3 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_3 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_3 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);

  // compile stiffness matrix without rotational effects
  b_elStorage[0] = elStorage->K11[0];
  b_elStorage[24] = elStorage->K12[0];
  b_elStorage[48] = elStorage->K13[0];
  b_elStorage[72] = elStorage->K14[0];
  b_elStorage[96] = elStorage->K15[0];
  b_elStorage[120] = elStorage->K16[0];
  b_elStorage[2] = elStorage->K12[0];
  b_elStorage[26] = elStorage->K22[0];
  b_elStorage[50] = elStorage->K23[0];
  b_elStorage[74] = elStorage->K24[0];
  b_elStorage[98] = elStorage->K25[0];
  b_elStorage[122] = elStorage->K26[0];
  b_elStorage[4] = elStorage->K13[0];
  b_elStorage[28] = elStorage->K23[0];
  b_elStorage[52] = elStorage->K33[0];
  b_elStorage[76] = elStorage->K34[0];
  b_elStorage[100] = elStorage->K35[0];
  b_elStorage[124] = elStorage->K36[0];
  b_elStorage[6] = elStorage->K13[0];
  b_elStorage[30] = elStorage->K24[0];
  b_elStorage[54] = elStorage->K34[0];
  b_elStorage[78] = elStorage->K44[0];
  b_elStorage[102] = elStorage->K45[0];
  b_elStorage[126] = elStorage->K46[0];
  b_elStorage[8] = elStorage->K15[0];
  b_elStorage[32] = elStorage->K25[0];
  b_elStorage[56] = elStorage->K35[0];
  b_elStorage[80] = elStorage->K45[0];
  b_elStorage[104] = elStorage->K55[0];
  b_elStorage[128] = elStorage->K56[0];
  b_elStorage[10] = elStorage->K16[0];
  b_elStorage[34] = elStorage->K26[0];
  b_elStorage[58] = elStorage->K36[0];
  b_elStorage[82] = elStorage->K46[0];
  b_elStorage[106] = elStorage->K56[0];
  b_elStorage[130] = elStorage->K66[0];
  b_elStorage[1] = elStorage->K11[1];
  b_elStorage[25] = elStorage->K12[1];
  b_elStorage[49] = elStorage->K13[1];
  b_elStorage[73] = elStorage->K14[1];
  b_elStorage[97] = elStorage->K15[1];
  b_elStorage[121] = elStorage->K16[1];
  b_elStorage[3] = elStorage->K12[2];
  b_elStorage[27] = elStorage->K22[1];
  b_elStorage[51] = elStorage->K23[1];
  b_elStorage[75] = elStorage->K24[1];
  b_elStorage[99] = elStorage->K25[1];
  b_elStorage[123] = elStorage->K26[1];
  b_elStorage[5] = elStorage->K13[2];
  b_elStorage[29] = elStorage->K23[2];
  b_elStorage[53] = elStorage->K33[1];
  b_elStorage[77] = elStorage->K34[1];
  b_elStorage[101] = elStorage->K35[1];
  b_elStorage[125] = elStorage->K36[1];
  b_elStorage[7] = elStorage->K13[2];
  b_elStorage[31] = elStorage->K24[2];
  b_elStorage[55] = elStorage->K34[2];
  b_elStorage[79] = elStorage->K44[1];
  b_elStorage[103] = elStorage->K45[1];
  b_elStorage[127] = elStorage->K46[1];
  b_elStorage[9] = elStorage->K15[2];
  b_elStorage[33] = elStorage->K25[2];
  b_elStorage[57] = elStorage->K35[2];
  b_elStorage[81] = elStorage->K45[2];
  b_elStorage[105] = elStorage->K55[1];
  b_elStorage[129] = elStorage->K56[1];
  b_elStorage[11] = elStorage->K16[2];
  b_elStorage[35] = elStorage->K26[2];
  b_elStorage[59] = elStorage->K36[2];
  b_elStorage[83] = elStorage->K46[2];
  b_elStorage[107] = elStorage->K56[2];
  b_elStorage[131] = elStorage->K66[1];
  b_elStorage[12] = elStorage->K11[2];
  b_elStorage[36] = elStorage->K12[2];
  b_elStorage[60] = elStorage->K13[2];
  b_elStorage[84] = elStorage->K14[2];
  b_elStorage[108] = elStorage->K15[2];
  b_elStorage[132] = elStorage->K16[2];
  b_elStorage[14] = elStorage->K12[1];
  b_elStorage[38] = elStorage->K22[2];
  b_elStorage[62] = elStorage->K23[2];
  b_elStorage[86] = elStorage->K24[2];
  b_elStorage[110] = elStorage->K25[2];
  b_elStorage[134] = elStorage->K26[2];
  b_elStorage[16] = elStorage->K13[1];
  b_elStorage[40] = elStorage->K23[1];
  b_elStorage[64] = elStorage->K33[2];
  b_elStorage[88] = elStorage->K34[2];
  b_elStorage[112] = elStorage->K35[2];
  b_elStorage[136] = elStorage->K36[2];
  b_elStorage[18] = elStorage->K13[1];
  b_elStorage[42] = elStorage->K24[1];
  b_elStorage[66] = elStorage->K34[1];
  b_elStorage[90] = elStorage->K44[2];
  b_elStorage[114] = elStorage->K45[2];
  b_elStorage[138] = elStorage->K46[2];
  b_elStorage[20] = elStorage->K15[1];
  b_elStorage[44] = elStorage->K25[1];
  b_elStorage[68] = elStorage->K35[1];
  b_elStorage[92] = elStorage->K45[1];
  b_elStorage[116] = elStorage->K55[2];
  b_elStorage[140] = elStorage->K56[2];
  b_elStorage[22] = elStorage->K16[1];
  b_elStorage[46] = elStorage->K26[1];
  b_elStorage[70] = elStorage->K36[1];
  b_elStorage[94] = elStorage->K46[1];
  b_elStorage[118] = elStorage->K56[1];
  b_elStorage[142] = elStorage->K66[2];
  b_elStorage[13] = elStorage->K11[3];
  b_elStorage[37] = elStorage->K12[3];
  b_elStorage[61] = elStorage->K13[3];
  b_elStorage[85] = elStorage->K14[3];
  b_elStorage[109] = elStorage->K15[3];
  b_elStorage[133] = elStorage->K16[3];
  b_elStorage[15] = elStorage->K12[3];
  b_elStorage[39] = elStorage->K22[3];
  b_elStorage[63] = elStorage->K23[3];
  b_elStorage[87] = elStorage->K24[3];
  b_elStorage[111] = elStorage->K25[3];
  b_elStorage[135] = elStorage->K26[3];
  b_elStorage[17] = elStorage->K13[3];
  b_elStorage[41] = elStorage->K23[3];
  b_elStorage[65] = elStorage->K33[3];
  b_elStorage[89] = elStorage->K34[3];
  b_elStorage[113] = elStorage->K35[3];
  b_elStorage[137] = elStorage->K36[3];
  b_elStorage[19] = elStorage->K13[3];
  b_elStorage[43] = elStorage->K24[3];
  b_elStorage[67] = elStorage->K34[3];
  b_elStorage[91] = elStorage->K44[3];
  b_elStorage[115] = elStorage->K45[3];
  b_elStorage[139] = elStorage->K46[3];
  b_elStorage[21] = elStorage->K15[3];
  b_elStorage[45] = elStorage->K25[3];
  b_elStorage[69] = elStorage->K35[3];
  b_elStorage[93] = elStorage->K45[3];
  b_elStorage[117] = elStorage->K55[3];
  b_elStorage[141] = elStorage->K56[3];
  b_elStorage[23] = elStorage->K16[3];
  b_elStorage[47] = elStorage->K26[3];
  b_elStorage[71] = elStorage->K36[3];
  b_elStorage[95] = elStorage->K46[3];
  b_elStorage[119] = elStorage->K56[3];
  b_elStorage[143] = elStorage->K66[3];
  mapMatrixNonSym(b_elStorage, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  c_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * c_tmp;
  K16_idx_0 = elStorage->K16[0] + elStorage->S16[0] * c_tmp;
  K56_idx_0 = elStorage->K56[0] + elStorage->S56[0] * c_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * c_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * c_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * c_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * c_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * c_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * c_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * c_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * c_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * c_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  b_elStorage[0] = elStorage->K11[0] + elStorage->S11[0] * c_tmp;
  b_elStorage[24] = (elStorage->K12[0] + integrationFactor) + O2;
  b_elStorage[48] = (elStorage->K13[0] + S13_idx_0) + H13_idx_0;
  O3dot = elStorage->K14[0] + ycm;
  b_elStorage[72] = O3dot + H14_idx_0;
  b_elStorage[96] = K15_idx_0;
  b_elStorage[120] = K16_idx_0;
  b_elStorage[2] = (elStorage->K12[0] + integrationFactor) - O2;
  b_elStorage[26] = elStorage->K22[0] + elStorage->S22[0] * c_c_tmp;
  O1dot = elStorage->K23[0] + O1;
  b_elStorage[50] = O1dot + H23_idx_0;
  O3 += elStorage->K24[0];
  b_elStorage[74] = O3 + H24_idx_0;
  O2 = elStorage->K25[0] + O2dot;
  b_elStorage[98] = O2 + H25_idx_0;
  O1 = elStorage->K26[0] + S26_idx_0;
  b_elStorage[122] = O1 + H26_idx_0;
  b_elStorage[4] = (elStorage->K13[0] + S13_idx_0) - H13_idx_0;
  b_elStorage[28] = O1dot - H23_idx_0;
  b_elStorage[52] = elStorage->K33[0] + elStorage->S33[0] * b_c_tmp;
  O1dot = elStorage->K34[0] + S34_idx_0;
  b_elStorage[76] = O1dot + H34_idx_0;
  integrationFactor = elStorage->K35[0] + S35_idx_0;
  b_elStorage[100] = integrationFactor + H35_idx_0;
  ycm = elStorage->K36[0] + S36_idx_0;
  b_elStorage[124] = ycm + H36_idx_0;
  b_elStorage[6] = O3dot - H14_idx_0;
  b_elStorage[30] = O3 - H24_idx_0;
  b_elStorage[54] = O1dot - H34_idx_0;
  b_elStorage[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * c_c_tmp +
    elStorage->S44_2[0] * b_c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  O3dot = elStorage->K45[0] + S45_idx_0;
  b_elStorage[102] = O3dot + H45_idx_0;
  O1dot = elStorage->K46[0] + S46_idx_0;
  b_elStorage[126] = O1dot + H46_idx_0;
  b_elStorage[8] = K15_idx_0;
  b_elStorage[32] = O2 - H25_idx_0;
  b_elStorage[56] = integrationFactor - H35_idx_0;
  b_elStorage[80] = O3dot - H45_idx_0;
  b_elStorage[104] = elStorage->K55[0] + elStorage->S55[0] * c_tmp;
  b_elStorage[128] = K56_idx_0;
  b_elStorage[10] = K16_idx_0;
  b_elStorage[34] = O1 - H26_idx_0;
  b_elStorage[58] = ycm - H36_idx_0;
  b_elStorage[82] = O1dot - H46_idx_0;
  b_elStorage[106] = K56_idx_0;
  b_elStorage[130] = elStorage->K66[0] + elStorage->S66[0] * c_tmp;
  b_elStorage[1] = elStorage->K11[1] + elStorage->S11[1] * c_tmp;
  b_elStorage[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  b_elStorage[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  O3dot = elStorage->K14[1] + S14_idx_1;
  b_elStorage[73] = O3dot + H14_idx_1;
  b_elStorage[97] = K15_idx_1;
  b_elStorage[121] = K16_idx_1;
  b_elStorage[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  b_elStorage[27] = elStorage->K22[1] + elStorage->S22[1] * c_c_tmp;
  O1dot = elStorage->K23[1] + S23_idx_1;
  b_elStorage[51] = O1dot + H23_idx_1;
  O3 = elStorage->K24[1] + S24_idx_1;
  b_elStorage[75] = O3 + H24_idx_1;
  O2 = elStorage->K25[1] + S25_idx_1;
  b_elStorage[99] = O2 + H25_idx_1;
  O1 = elStorage->K26[1] + S26_idx_1;
  b_elStorage[123] = O1 + H26_idx_1;
  b_elStorage[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  integrationFactor = elStorage->K23[2] + S23_idx_2;
  b_elStorage[29] = integrationFactor - H23_idx_2;
  b_elStorage[53] = elStorage->K33[1] + elStorage->S33[1] * b_c_tmp;
  ycm = elStorage->K34[1] + S34_idx_1;
  b_elStorage[77] = ycm + H34_idx_1;
  rhoA = elStorage->K35[1] + S35_idx_1;
  b_elStorage[101] = rhoA + H35_idx_1;
  zcm = elStorage->K36[1] + S36_idx_1;
  b_elStorage[125] = zcm + H36_idx_1;
  c_f = elStorage->K14[2] + S14_idx_2;
  b_elStorage[7] = c_f - H14_idx_2;
  b_f = elStorage->K24[2] + S24_idx_2;
  b_elStorage[31] = b_f - H24_idx_2;
  f = elStorage->K34[2] + S34_idx_2;
  b_elStorage[55] = f - H34_idx_2;
  b_elStorage[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * c_c_tmp +
    elStorage->S44_2[1] * b_c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  e_f_tmp = elStorage->K45[1] + S45_idx_1;
  b_elStorage[103] = e_f_tmp + H45_idx_1;
  d_f_tmp = elStorage->K46[1] + S46_idx_1;
  b_elStorage[127] = d_f_tmp + H46_idx_1;
  b_elStorage[9] = K15_idx_2;
  c_f_tmp = elStorage->K25[2] + S25_idx_2;
  b_elStorage[33] = c_f_tmp - H25_idx_2;
  b_f_tmp = elStorage->K35[2] + S35_idx_2;
  b_elStorage[57] = b_f_tmp - H35_idx_2;
  f_tmp = elStorage->K45[2] + S45_idx_2;
  b_elStorage[81] = f_tmp - H45_idx_2;
  b_elStorage[105] = elStorage->K55[1] + elStorage->S55[1] * c_tmp;
  b_elStorage[129] = K56_idx_1;
  b_elStorage[11] = K16_idx_2;
  K56_idx_0 = elStorage->K26[2] + S26_idx_2;
  b_elStorage[35] = K56_idx_0 - H26_idx_2;
  K16_idx_0 = elStorage->K36[2] + S36_idx_2;
  b_elStorage[59] = K16_idx_0 - H36_idx_2;
  K15_idx_0 = elStorage->K46[2] + S46_idx_2;
  b_elStorage[83] = K15_idx_0 - H46_idx_2;
  b_elStorage[107] = K56_idx_2;
  b_elStorage[131] = elStorage->K66[1] + elStorage->S66[1] * c_tmp;
  b_elStorage[12] = elStorage->K11[2] + elStorage->S11[2] * c_tmp;
  b_elStorage[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  b_elStorage[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  b_elStorage[84] = c_f + H14_idx_2;
  b_elStorage[108] = K15_idx_2;
  b_elStorage[132] = K16_idx_2;
  b_elStorage[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  b_elStorage[38] = elStorage->K22[2] + elStorage->S22[2] * c_c_tmp;
  b_elStorage[62] = integrationFactor + H23_idx_2;
  b_elStorage[86] = b_f + H24_idx_2;
  b_elStorage[110] = c_f_tmp + H25_idx_2;
  b_elStorage[134] = K56_idx_0 + H26_idx_2;
  b_elStorage[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  b_elStorage[40] = O1dot - H23_idx_1;
  b_elStorage[64] = elStorage->K33[2] + elStorage->S33[2] * b_c_tmp;
  b_elStorage[88] = f + H34_idx_2;
  b_elStorage[112] = b_f_tmp + H35_idx_2;
  b_elStorage[136] = K16_idx_0 + H36_idx_2;
  b_elStorage[18] = O3dot - H14_idx_1;
  b_elStorage[42] = O3 - H24_idx_1;
  b_elStorage[66] = ycm - H34_idx_1;
  b_elStorage[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * c_c_tmp +
    elStorage->S44_2[2] * b_c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  b_elStorage[114] = f_tmp + H45_idx_2;
  b_elStorage[138] = K15_idx_0 + H46_idx_2;
  b_elStorage[20] = K15_idx_1;
  b_elStorage[44] = O2 - H25_idx_1;
  b_elStorage[68] = rhoA - H35_idx_1;
  b_elStorage[92] = e_f_tmp - H45_idx_1;
  b_elStorage[116] = elStorage->K55[2] + elStorage->S55[2] * c_tmp;
  b_elStorage[140] = K56_idx_2;
  b_elStorage[22] = K16_idx_1;
  b_elStorage[46] = O1 - H26_idx_1;
  b_elStorage[70] = zcm - H36_idx_1;
  b_elStorage[94] = d_f_tmp - H46_idx_1;
  b_elStorage[118] = K56_idx_1;
  b_elStorage[142] = elStorage->K66[2] + elStorage->S66[2] * c_tmp;
  b_elStorage[13] = elStorage->K11[3] + elStorage->S11[3] * c_tmp;
  b_elStorage[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  b_elStorage[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  O3dot = elStorage->K14[3] + S14_idx_3;
  b_elStorage[85] = O3dot + H14_idx_3;
  b_elStorage[109] = K15_idx_3;
  b_elStorage[133] = K16_idx_3;
  b_elStorage[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  b_elStorage[39] = elStorage->K22[3] + elStorage->S22[3] * c_c_tmp;
  O1dot = elStorage->K23[3] + S23_idx_3;
  b_elStorage[63] = O1dot + H23_idx_3;
  O3 = elStorage->K24[3] + S24_idx_3;
  b_elStorage[87] = O3 + H24_idx_3;
  O2 = elStorage->K25[3] + S25_idx_3;
  b_elStorage[111] = O2 + H25_idx_3;
  O1 = elStorage->K26[3] + S26_idx_3;
  b_elStorage[135] = O1 + H26_idx_3;
  b_elStorage[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  b_elStorage[41] = O1dot - H23_idx_3;
  b_elStorage[65] = elStorage->K33[3] + elStorage->S33[3] * b_c_tmp;
  O1dot = elStorage->K34[3] + S34_idx_3;
  b_elStorage[89] = O1dot + H34_idx_3;
  integrationFactor = elStorage->K35[3] + S35_idx_3;
  b_elStorage[113] = integrationFactor + H35_idx_3;
  ycm = elStorage->K36[3] + S36_idx_3;
  b_elStorage[137] = ycm + H36_idx_3;
  b_elStorage[19] = O3dot - H14_idx_3;
  b_elStorage[43] = O3 - H24_idx_3;
  b_elStorage[67] = O1dot - H34_idx_3;
  b_elStorage[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * c_c_tmp +
    elStorage->S44_2[3] * b_c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  O3dot = elStorage->K45[3] + S45_idx_3;
  b_elStorage[115] = O3dot + H45_idx_3;
  O1dot = elStorage->K46[3] + S46_idx_3;
  b_elStorage[139] = O1dot + H46_idx_3;
  b_elStorage[21] = K15_idx_3;
  b_elStorage[45] = O2 - H25_idx_3;
  b_elStorage[69] = integrationFactor - H35_idx_3;
  b_elStorage[93] = O3dot - H45_idx_3;
  b_elStorage[117] = elStorage->K55[3] + elStorage->S55[3] * c_tmp;
  b_elStorage[141] = K56_idx_3;
  b_elStorage[23] = K16_idx_3;
  b_elStorage[47] = O1 - H26_idx_3;
  b_elStorage[71] = ycm - H36_idx_3;
  b_elStorage[95] = O1dot - H46_idx_3;
  b_elStorage[119] = K56_idx_3;
  b_elStorage[143] = elStorage->K66[3] + elStorage->S66[3] * c_tmp;
  mapMatrixNonSym(b_elStorage, Khate);

  //  Declare type
  // compile Coriolis/damping matrix
  reshapes_f3_data[0] = -C13_idx_0;
  reshapes_f3_data[4] = -C23_idx_0;
  reshapes_f3_data[8] = 0.0;
  reshapes_f3_data[12] = C34_idx_0;
  reshapes_f3_data[16] = C35_idx_0;
  reshapes_f3_data[20] = C36_idx_0;
  reshapes_f4_data[0] = -C14_idx_0;
  reshapes_f4_data[4] = -C24_idx_0;
  reshapes_f4_data[8] = -C34_idx_0;
  reshapes_f4_data[12] = 0.0;
  reshapes_f4_data[16] = C45_idx_0;
  reshapes_f4_data[20] = C46_idx_0;
  reshapes_f1[0] = 0.0;
  reshapes_f1[4] = C12_idx_0;
  reshapes_f1[8] = C13_idx_0;
  reshapes_f1[12] = C14_idx_0;
  reshapes_f1[16] = 0.0;
  reshapes_f1[20] = 0.0;
  reshapes_f2[0] = -C12_idx_0;
  reshapes_f2[4] = 0.0;
  reshapes_f2[8] = C23_idx_0;
  reshapes_f2[12] = C24_idx_0;
  reshapes_f2[16] = C25_idx_0;
  reshapes_f2[20] = C26_idx_0;
  reshapes_f5[0] = 0.0;
  reshapes_f5[4] = -C25_idx_0;
  reshapes_f5[8] = -C35_idx_0;
  reshapes_f5[12] = -C45_idx_0;
  reshapes_f5[16] = 0.0;
  reshapes_f5[20] = 0.0;
  reshapes_f6[0] = 0.0;
  reshapes_f6[4] = -C26_idx_0;
  reshapes_f6[8] = -C36_idx_0;
  reshapes_f6[12] = -C46_idx_0;
  reshapes_f6[16] = 0.0;
  reshapes_f6[20] = 0.0;
  reshapes_f3_data[1] = -C13_idx_2;
  reshapes_f3_data[5] = -C23_idx_2;
  reshapes_f3_data[9] = 0.0;
  reshapes_f3_data[13] = C34_idx_1;
  reshapes_f3_data[17] = C35_idx_1;
  reshapes_f3_data[21] = C36_idx_1;
  reshapes_f4_data[1] = -C14_idx_2;
  reshapes_f4_data[5] = -C24_idx_2;
  reshapes_f4_data[9] = -C34_idx_2;
  reshapes_f4_data[13] = 0.0;
  reshapes_f4_data[17] = C45_idx_1;
  reshapes_f4_data[21] = C46_idx_1;
  reshapes_f1[1] = 0.0;
  reshapes_f1[5] = C12_idx_1;
  reshapes_f1[9] = C13_idx_1;
  reshapes_f1[13] = C14_idx_1;
  reshapes_f1[17] = 0.0;
  reshapes_f1[21] = 0.0;
  reshapes_f2[1] = -C12_idx_2;
  reshapes_f2[5] = 0.0;
  reshapes_f2[9] = C23_idx_1;
  reshapes_f2[13] = C24_idx_1;
  reshapes_f2[17] = C25_idx_1;
  reshapes_f2[21] = C26_idx_1;
  reshapes_f5[1] = 0.0;
  reshapes_f5[5] = -C25_idx_2;
  reshapes_f5[9] = -C35_idx_2;
  reshapes_f5[13] = -C45_idx_2;
  reshapes_f5[17] = 0.0;
  reshapes_f5[21] = 0.0;
  reshapes_f6[1] = 0.0;
  reshapes_f6[5] = -C26_idx_2;
  reshapes_f6[9] = -C36_idx_2;
  reshapes_f6[13] = -C46_idx_2;
  reshapes_f6[17] = 0.0;
  reshapes_f6[21] = 0.0;
  reshapes_f3_data[2] = -C13_idx_1;
  reshapes_f3_data[6] = -C23_idx_1;
  reshapes_f3_data[10] = 0.0;
  reshapes_f3_data[14] = C34_idx_2;
  reshapes_f3_data[18] = C35_idx_2;
  reshapes_f3_data[22] = C36_idx_2;
  reshapes_f4_data[2] = -C14_idx_1;
  reshapes_f4_data[6] = -C24_idx_1;
  reshapes_f4_data[10] = -C34_idx_1;
  reshapes_f4_data[14] = 0.0;
  reshapes_f4_data[18] = C45_idx_2;
  reshapes_f4_data[22] = C46_idx_2;
  reshapes_f1[2] = 0.0;
  reshapes_f1[6] = C12_idx_2;
  reshapes_f1[10] = C13_idx_2;
  reshapes_f1[14] = C14_idx_2;
  reshapes_f1[18] = 0.0;
  reshapes_f1[22] = 0.0;
  reshapes_f2[2] = -C12_idx_1;
  reshapes_f2[6] = 0.0;
  reshapes_f2[10] = C23_idx_2;
  reshapes_f2[14] = C24_idx_2;
  reshapes_f2[18] = C25_idx_2;
  reshapes_f2[22] = C26_idx_2;
  reshapes_f5[2] = 0.0;
  reshapes_f5[6] = -C25_idx_1;
  reshapes_f5[10] = -C35_idx_1;
  reshapes_f5[14] = -C45_idx_1;
  reshapes_f5[18] = 0.0;
  reshapes_f5[22] = 0.0;
  reshapes_f6[2] = 0.0;
  reshapes_f6[6] = -C26_idx_1;
  reshapes_f6[10] = -C36_idx_1;
  reshapes_f6[14] = -C46_idx_1;
  reshapes_f6[18] = 0.0;
  reshapes_f6[22] = 0.0;
  reshapes_f3_data[3] = -C13_idx_3;
  reshapes_f3_data[7] = -C23_idx_3;
  reshapes_f3_data[11] = 0.0;
  reshapes_f3_data[15] = C34_idx_3;
  reshapes_f3_data[19] = C35_idx_3;
  reshapes_f3_data[23] = C36_idx_3;
  reshapes_f4_data[3] = -C14_idx_3;
  reshapes_f4_data[7] = -C24_idx_3;
  reshapes_f4_data[11] = -C34_idx_3;
  reshapes_f4_data[15] = 0.0;
  reshapes_f4_data[19] = C45_idx_3;
  reshapes_f4_data[23] = C46_idx_3;
  reshapes_f1[3] = 0.0;
  reshapes_f1[7] = C12_idx_3;
  reshapes_f1[11] = C13_idx_3;
  reshapes_f1[15] = C14_idx_3;
  reshapes_f1[19] = 0.0;
  reshapes_f1[23] = 0.0;
  reshapes_f2[3] = -C12_idx_3;
  reshapes_f2[7] = 0.0;
  reshapes_f2[11] = C23_idx_3;
  reshapes_f2[15] = C24_idx_3;
  reshapes_f2[19] = C25_idx_3;
  reshapes_f2[23] = C26_idx_3;
  reshapes_f5[3] = 0.0;
  reshapes_f5[7] = -C25_idx_3;
  reshapes_f5[11] = -C35_idx_3;
  reshapes_f5[15] = -C45_idx_3;
  reshapes_f5[19] = 0.0;
  reshapes_f5[23] = 0.0;
  reshapes_f6[3] = 0.0;
  reshapes_f6[7] = -C26_idx_3;
  reshapes_f6[11] = -C36_idx_3;
  reshapes_f6[15] = -C46_idx_3;
  reshapes_f6[19] = 0.0;
  reshapes_f6[23] = 0.0;
  for (i = 0; i < 12; i++) {
    b_i = i << 1;
    b_elStorage[12 * i] = reshapes_f1[b_i];
    b_elStorage[12 * i + 2] = reshapes_f2[b_i];
    b_elStorage[12 * i + 4] = reshapes_f3_data[2 * i];
    b_elStorage[12 * i + 6] = reshapes_f4_data[2 * i];
    b_elStorage[12 * i + 8] = reshapes_f5[b_i];
    b_elStorage[12 * i + 10] = reshapes_f6[b_i];
    b_i++;
    b_elStorage[12 * i + 1] = reshapes_f1[b_i];
    b_elStorage[12 * i + 3] = reshapes_f2[b_i];
    reshapes_f1_data_tmp = 2 * i + 1;
    b_elStorage[12 * i + 5] = reshapes_f3_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 7] = reshapes_f4_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 9] = reshapes_f5[b_i];
    b_elStorage[12 * i + 11] = reshapes_f6[b_i];
  }

  b_mapMatrixNonSym(b_elStorage, Ce);

  // compile mass matrix
  b_elStorage[0] = elStorage->M11[0];
  b_elStorage[24] = 0.0;
  b_elStorage[48] = 0.0;
  b_elStorage[72] = 0.0;
  b_elStorage[96] = elStorage->M15[0];
  b_elStorage[120] = elStorage->M16[0];
  b_elStorage[2] = 0.0;
  b_elStorage[26] = elStorage->M22[0];
  b_elStorage[50] = 0.0;
  b_elStorage[74] = elStorage->M24[0];
  b_elStorage[98] = 0.0;
  b_elStorage[122] = 0.0;
  b_elStorage[4] = 0.0;
  b_elStorage[28] = 0.0;
  b_elStorage[52] = elStorage->M33[0];
  b_elStorage[76] = elStorage->M34[0];
  b_elStorage[100] = 0.0;
  b_elStorage[124] = 0.0;
  b_elStorage[6] = 0.0;
  b_elStorage[30] = elStorage->M24[0];
  b_elStorage[54] = elStorage->M34[0];
  b_elStorage[78] = elStorage->M44[0];
  b_elStorage[102] = 0.0;
  b_elStorage[126] = 0.0;
  b_elStorage[8] = elStorage->M15[0];
  b_elStorage[32] = 0.0;
  b_elStorage[56] = 0.0;
  b_elStorage[80] = 0.0;
  b_elStorage[104] = elStorage->M55[0];
  b_elStorage[128] = elStorage->M56[0];
  b_elStorage[10] = elStorage->M16[0];
  b_elStorage[34] = 0.0;
  b_elStorage[58] = 0.0;
  b_elStorage[82] = 0.0;
  b_elStorage[106] = elStorage->M56[0];
  b_elStorage[130] = elStorage->M66[0];
  b_elStorage[1] = elStorage->M11[1];
  b_elStorage[25] = 0.0;
  b_elStorage[49] = 0.0;
  b_elStorage[73] = 0.0;
  b_elStorage[97] = elStorage->M15[1];
  b_elStorage[121] = elStorage->M16[1];
  b_elStorage[3] = 0.0;
  b_elStorage[27] = elStorage->M22[1];
  b_elStorage[51] = 0.0;
  b_elStorage[75] = elStorage->M24[1];
  b_elStorage[99] = 0.0;
  b_elStorage[123] = 0.0;
  b_elStorage[5] = 0.0;
  b_elStorage[29] = 0.0;
  b_elStorage[53] = elStorage->M33[1];
  b_elStorage[77] = elStorage->M34[1];
  b_elStorage[101] = 0.0;
  b_elStorage[125] = 0.0;
  b_elStorage[7] = 0.0;
  b_elStorage[31] = elStorage->M24[2];
  b_elStorage[55] = elStorage->M34[2];
  b_elStorage[79] = elStorage->M44[1];
  b_elStorage[103] = 0.0;
  b_elStorage[127] = 0.0;
  b_elStorage[9] = elStorage->M15[2];
  b_elStorage[33] = 0.0;
  b_elStorage[57] = 0.0;
  b_elStorage[81] = 0.0;
  b_elStorage[105] = elStorage->M55[1];
  b_elStorage[129] = elStorage->M56[1];
  b_elStorage[11] = elStorage->M16[2];
  b_elStorage[35] = 0.0;
  b_elStorage[59] = 0.0;
  b_elStorage[83] = 0.0;
  b_elStorage[107] = elStorage->M56[2];
  b_elStorage[131] = elStorage->M66[1];
  b_elStorage[12] = elStorage->M11[2];
  b_elStorage[36] = 0.0;
  b_elStorage[60] = 0.0;
  b_elStorage[84] = 0.0;
  b_elStorage[108] = elStorage->M15[2];
  b_elStorage[132] = elStorage->M16[2];
  b_elStorage[14] = 0.0;
  b_elStorage[38] = elStorage->M22[2];
  b_elStorage[62] = 0.0;
  b_elStorage[86] = elStorage->M24[2];
  b_elStorage[110] = 0.0;
  b_elStorage[134] = 0.0;
  b_elStorage[16] = 0.0;
  b_elStorage[40] = 0.0;
  b_elStorage[64] = elStorage->M33[2];
  b_elStorage[88] = elStorage->M34[2];
  b_elStorage[112] = 0.0;
  b_elStorage[136] = 0.0;
  b_elStorage[18] = 0.0;
  b_elStorage[42] = elStorage->M24[1];
  b_elStorage[66] = elStorage->M34[1];
  b_elStorage[90] = elStorage->M44[2];
  b_elStorage[114] = 0.0;
  b_elStorage[138] = 0.0;
  b_elStorage[20] = elStorage->M15[1];
  b_elStorage[44] = 0.0;
  b_elStorage[68] = 0.0;
  b_elStorage[92] = 0.0;
  b_elStorage[116] = elStorage->M55[2];
  b_elStorage[140] = elStorage->M56[2];
  b_elStorage[22] = elStorage->M16[1];
  b_elStorage[46] = 0.0;
  b_elStorage[70] = 0.0;
  b_elStorage[94] = 0.0;
  b_elStorage[118] = elStorage->M56[1];
  b_elStorage[142] = elStorage->M66[2];
  b_elStorage[13] = elStorage->M11[3];
  b_elStorage[37] = 0.0;
  b_elStorage[61] = 0.0;
  b_elStorage[85] = 0.0;
  b_elStorage[109] = elStorage->M15[3];
  b_elStorage[133] = elStorage->M16[3];
  b_elStorage[15] = 0.0;
  b_elStorage[39] = elStorage->M22[3];
  b_elStorage[63] = 0.0;
  b_elStorage[87] = elStorage->M24[3];
  b_elStorage[111] = 0.0;
  b_elStorage[135] = 0.0;
  b_elStorage[17] = 0.0;
  b_elStorage[41] = 0.0;
  b_elStorage[65] = elStorage->M33[3];
  b_elStorage[89] = elStorage->M34[3];
  b_elStorage[113] = 0.0;
  b_elStorage[137] = 0.0;
  b_elStorage[19] = 0.0;
  b_elStorage[43] = elStorage->M24[3];
  b_elStorage[67] = elStorage->M34[3];
  b_elStorage[91] = elStorage->M44[3];
  b_elStorage[115] = 0.0;
  b_elStorage[139] = 0.0;
  b_elStorage[21] = elStorage->M15[3];
  b_elStorage[45] = 0.0;
  b_elStorage[69] = 0.0;
  b_elStorage[93] = 0.0;
  b_elStorage[117] = elStorage->M55[3];
  b_elStorage[141] = elStorage->M56[3];
  b_elStorage[23] = elStorage->M16[3];
  b_elStorage[47] = 0.0;
  b_elStorage[71] = 0.0;
  b_elStorage[95] = 0.0;
  b_elStorage[119] = elStorage->M56[3];
  b_elStorage[143] = elStorage->M66[3];
  mapMatrixNonSym(b_elStorage, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input->RayleighAlpha * Kenr[i] + input->RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      b_elStorage[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(b_elStorage, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Khate,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Khate);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fhate[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fhate[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fhate, Fe);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input->concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input->concMass[0];
    Me[13] += input->concMass[0];
    Me[26] += input->concMass[0];
    Me[39] += input->concMass[1];
    Me[52] += input->concMass[2];
    Me[65] += input->concMass[3];
    Me[78] += input->concMass[4];
    Me[91] += input->concMass[4];
    Me[104] += input->concMass[4];
    Me[117] += input->concMass[5];
    Me[130] += input->concMass[6];
    Me[143] += input->concMass[7];

    // modify Ce for concentrated mass
    H46_idx_3 = 2.0 * input->concMass[0] * Omega;
    Ce[12] -= H46_idx_3;
    Ce[1] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[0] * 0.0;
    Ce[24] += H46_idx_3;
    Ce[2] -= H46_idx_3;
    Ce[25] -= H46_idx_3;
    Ce[14] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[4] * Omega;
    Ce[90] -= H46_idx_3;
    Ce[79] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[4] * 0.0;
    Ce[102] += H46_idx_3;
    Ce[80] -= H46_idx_3;
    Ce[103] -= H46_idx_3;
    Ce[92] += H46_idx_3;

    // modify Ke for concentrated mass
    H46_idx_3 = Omega * Omega;
    K56_idx_0 = input->concMass[0] * H46_idx_3;
    Khate[0] -= K56_idx_0;
    K15_idx_0 = input->concMass[0] * 0.0 * 0.0;
    K16_idx_0 = input->concMass[0] * OmegaDot;
    Khate[12] = (Khate[12] + K15_idx_0) - K16_idx_0;
    Khate[1] = (Khate[1] + K15_idx_0) + K16_idx_0;
    K15_idx_0 = input->concMass[0] * 0.0 * Omega;
    Khate[24] = (Khate[24] + K15_idx_0) + input->concMass[0] * 0.0;
    Khate[2] = (Khate[2] + K15_idx_0) - input->concMass[0] * 0.0;
    Khate[25] = (Khate[25] + K15_idx_0) - input->concMass[0] * 0.0;
    Khate[14] = (Khate[14] + K15_idx_0) + input->concMass[0] * 0.0;
    Khate[13] -= K56_idx_0;
    Khate[26] -= input->concMass[0] * 0.0;
    K56_idx_0 = input->concMass[4] * H46_idx_3;
    Khate[78] -= K56_idx_0;
    K15_idx_0 = input->concMass[4] * 0.0 * 0.0;
    K16_idx_0 = input->concMass[4] * OmegaDot;
    Khate[90] = (Khate[90] + K15_idx_0) - K16_idx_0;
    Khate[79] = (Khate[79] + K15_idx_0) + K16_idx_0;
    K15_idx_0 = input->concMass[4] * 0.0 * Omega;
    Khate[102] = (Khate[102] + K15_idx_0) + input->concMass[4] * 0.0;
    Khate[80] = (Khate[80] + K15_idx_0) - input->concMass[4] * 0.0;
    Khate[103] = (Khate[103] + K15_idx_0) - input->concMass[4] * 0.0;
    Khate[92] = (Khate[92] + K15_idx_0) + input->concMass[4] * 0.0;
    Khate[91] -= K56_idx_0;
    Khate[104] -= input->concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    H46_idx_3 = Omega * Omega;
    K56_idx_0 = 0.0 * Omega * input->z.data[0];
    Fe[0] = ((Fe[0] + input->concMass[0] * ((input->x.data[0] * H46_idx_3 - 0.0 *
                input->y.data[0]) - K56_idx_0)) + input->concMass[0] *
             (input->y.data[0] * OmegaDot - input->z.data[0] * 0.0)) -
      input->concMass[0] * a_temp[0];
    Fe[1] = ((Fe[1] + input->concMass[0] * ((input->y.data[0] * H46_idx_3 -
                K56_idx_0) - 0.0 * input->x.data[0])) + input->concMass[0] *
             (input->z.data[0] * 0.0 - input->x.data[0] * OmegaDot)) -
      input->concMass[0] * a_temp[1];
    Fe[2] = ((Fe[2] + input->concMass[0] * ((input->z.data[0] * 0.0 - Omega *
                0.0 * input->x.data[0]) - Omega * 0.0 * input->y.data[0])) +
             input->concMass[0] * (input->x.data[0] * 0.0 - input->y.data[0] *
              0.0)) - input->concMass[0] * a_temp[2];
    K56_idx_0 = 0.0 * Omega * input->z.data[1];
    Fe[6] = ((Fe[6] + input->concMass[4] * ((input->x.data[1] * H46_idx_3 - 0.0 *
                input->y.data[1]) - K56_idx_0)) + input->concMass[4] *
             (input->y.data[1] * OmegaDot - input->z.data[1] * 0.0)) -
      input->concMass[4] * a_temp[0];
    Fe[7] = ((Fe[7] + input->concMass[4] * ((input->y.data[1] * H46_idx_3 -
                K56_idx_0) - 0.0 * input->x.data[1])) + input->concMass[4] *
             (input->z.data[1] * 0.0 - input->x.data[1] * OmegaDot)) -
      input->concMass[4] * a_temp[1];
    Fe[8] = ((Fe[8] + input->concMass[4] * ((input->z.data[1] * 0.0 - Omega *
                0.0 * input->x.data[1]) - Omega * 0.0 * input->y.data[1])) +
             input->concMass[4] * (input->x.data[1] * 0.0 - input->y.data[1] *
              0.0)) - input->concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  std::memset(&Fhate[0], 0, 12U * sizeof(double));
  std::memset(&FhatLessConc[0], 0, 12U * sizeof(double));
  if (g_strcmp(input->analysisType.data, input->analysisType.size)) {
    // calculate effective stiffness matrix and force vector for Dean integrator 
    for (i = 0; i < 12; i++) {
      FhatLessConc[i] = dispm1_data[i];
      b_Khate[i] = 2.0 * input->disp.data[i] - dispm1_data[i];
      Fhate[i] = -0.001 * dispm1_data[i] - 0.001 * input->disp.data[i];
    }

    for (i = 0; i < 12; i++) {
      c_tmp = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        c_tmp += Me[i + 12 * b_i] * b_Khate[b_i];
      }

      Ftemp_data[i] = Fe[i] * 2000.0 + c_tmp;
    }

    for (i = 0; i < 12; i++) {
      c_tmp = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        c_tmp += Khate[i + 12 * b_i] * Fhate[b_i];
      }

      b_Khate[i] = c_tmp;
    }

    for (i = 0; i < 12; i++) {
      c_tmp = 0.0;
      for (b_i = 0; b_i < 12; b_i++) {
        c_tmp += Ce[i + 12 * b_i] * (1.0E+6 * FhatLessConc[b_i]);
      }

      Fhate[i] = (Ftemp_data[i] + b_Khate[i]) + c_tmp;
    }

    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));

    // ........................................................
    // ..........................................................
    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] * 0.001 + 1.0E+6 * Ce[i]) + Me[i];
    }

    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  if (h_strcmp(input->analysisType.data, input->analysisType.size)) {
    // calculate effective stiffness matrix and load vector for Newmark-Beta integrator 
    //      a1 = timeInt.a1;
    //      a2 = timeInt.a2;
    if (e_strcmp(input->iterationType)) {
      // considerations if newton raphson iteration is used
      if (0 <= dispddot_size_idx_1 - 1) {
        std::memcpy(&b_data[0], &dispddot_data[0], dispddot_size_idx_1 * sizeof
                    (double));
      }

      if (dispddot_size_idx_1 == 1) {
        for (i = 0; i < 12; i++) {
          c_tmp = 0.0;
          for (b_i = 0; b_i < 12; b_i++) {
            c_tmp += Me[i + 12 * b_i] * b_data[b_i];
          }

          Fhate[i] = c_tmp;
        }
      } else {
        mtimes(Me, b_data, Fhate);
      }

      if (0 <= dispdot_size_idx_1 - 1) {
        std::memcpy(&b_data[0], &dispdot_data[0], dispdot_size_idx_1 * sizeof
                    (double));
      }

      if (dispdot_size_idx_1 == 1) {
        for (i = 0; i < 12; i++) {
          c_tmp = 0.0;
          for (b_i = 0; b_i < 12; b_i++) {
            c_tmp += Ce[i + 12 * b_i] * b_data[b_i];
          }

          b_Khate[i] = c_tmp;
        }
      } else {
        mtimes(Ce, b_data, b_Khate);
      }

      std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
      mtimes(Khate, b_data, FhatLessConc);
      for (i = 0; i < 12; i++) {
        Fhate[i] = ((Fe[i] - Fhate[i]) - b_Khate[i]) - FhatLessConc[i];
      }
    } else {
      if (f_strcmp(input->iterationType)) {
        // considerations if direct iteration is used or linear analysis
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, b_Khate);
        for (i = 0; i < 12; i++) {
          Fhate[i] = (Fe[i] + Fhate[i]) + b_Khate[i];
        }
      }
    }

    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] + 1.0E+6 * Me[i]) + 1000.0 * Ce[i];
    }

    // ........................................................
    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));
    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  if (i_strcmp(input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&FhatLessConc[0], &Fe[0], 12U * sizeof(double));
  }

  if ((i_strcmp(input->analysisType.data, input->analysisType.size) || k_strcmp
       (input->analysisType.data, input->analysisType.size)) && e_strcmp
      (input->iterationType)) {
    // considerations for newton-raphson iteration
    std::memcpy(&b_data[0], &disp_iter_data[0], 12U * sizeof(double));
    mtimes(Khate, b_data, b_Khate);
    for (i = 0; i < 12; i++) {
      Fe[i] -= b_Khate[i];
    }
  }

  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Khate[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &Fe[0], 12U * sizeof(double));
  output->Me.size[0] = 1;
  output->Me.size[1] = 1;
  output->Me.data[0] = 0.0;
  output->Ce.size[0] = 1;
  output->Ce.size[1] = 1;
  output->Ce.data[0] = 0.0;
  if (i_strcmp(input->analysisType.data, input->analysisType.size) || j_strcmp
      (input->analysisType.data, input->analysisType.size)) {
    output->Me.size[0] = 12;
    output->Me.size[1] = 12;
    output->Ce.size[0] = 12;
    output->Ce.size[1] = 12;
    std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
    std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));
  }

  if (g_strcmp(input->analysisType.data, input->analysisType.size) || h_strcmp
      (input->analysisType.data, input->analysisType.size)) {
    std::memcpy(&output->FhatLessConc[0], &FhatLessConc[0], 12U * sizeof(double));
  }

  // ------------------------------------------
}

//
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
// Arguments    : const double input_xloc[2]
//                const double input_sectionProps_twist[2]
//                const double input_sectionProps_rhoA[2]
//                const double input_sectionProps_EA[2]
//                const double input_sectionProps_zcm[2]
//                const double input_sectionProps_ycm[2]
//                double input_sweepAngle
//                double input_coneAngle
//                double input_rollAngle
//                const double input_concMass[8]
//                const double input_disp_data[]
//                const double input_x_data[]
//                const double input_y_data[]
//                const double input_z_data[]
//                const double input_CN2H[9]
//                double input_RayleighAlpha
//                double input_RayleighBeta
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void c_calculateTimoshenkoElementNL(const double input_xloc[2], const double
  input_sectionProps_twist[2], const double input_sectionProps_rhoA[2], const
  double input_sectionProps_EA[2], const double input_sectionProps_zcm[2], const
  double input_sectionProps_ycm[2], double input_sweepAngle, double
  input_coneAngle, double input_rollAngle, const double input_concMass[8], const
  double input_disp_data[], const double input_x_data[], const double
  input_y_data[], const double input_z_data[], const double input_CN2H[9],
  double input_RayleighAlpha, double input_RayleighBeta, const f_struct_T
  *elStorage, o_struct_T *output)
{
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double SS22_data_idx_0;
  double SS33_data_idx_0;
  double SS22_data_idx_1;
  double SS33_data_idx_1;
  double SS22_data_idx_2;
  double SS33_data_idx_2;
  double SS22_data_idx_3;
  double SS33_data_idx_3;
  double lambda[144];
  double b_data[12];
  double dispLocal[12];
  int i;
  double O1;
  double Oel[3];
  double K15_idx_1;
  double O2;
  double K16_idx_1;
  double O3;
  double a_temp[3];
  double O1dot;
  double ODotel[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double f_tmp;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double K15_idx_0;
  double integrationFactor;
  double c_tmp;
  double b_c_tmp;
  double rhoA;
  double S13_idx_0;
  double Faxial;
  double K16_idx_0;
  double S35_idx_0;
  double S36_idx_0;
  double ycm;
  double S34_idx_0;
  double S45_idx_0;
  double S46_idx_0;
  double C12_idx_0;
  double C13_idx_0;
  double C23_idx_0;
  double K56_idx_0;
  double C24_idx_0;
  double b_f_tmp;
  double C25_idx_0;
  double f;
  double C26_idx_0;
  double c_f_tmp;
  double C34_idx_0;
  double d_f_tmp;
  double C35_idx_0;
  double e_f_tmp;
  double C36_idx_0;
  double zcm;
  double f_f_tmp;
  double g_f_tmp;
  double C14_idx_0;
  double b_f;
  double c_f;
  double disMomentgp[3];
  double C45_idx_0;
  double d_f;
  double C46_idx_0;
  double H13_idx_0;
  double H23_idx_0;
  double posLocal[3];
  double H24_idx_0;
  double disLoadgpLocal[3];
  double H25_idx_0;
  double H26_idx_0;
  double H34_idx_0;
  double H35_idx_0;
  double H36_idx_0;
  double H14_idx_0;
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double S35_idx_1;
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C12_idx_1;
  double C13_idx_1;
  double C23_idx_1;
  double C24_idx_1;
  double C25_idx_1;
  double C26_idx_1;
  double C34_idx_1;
  double C35_idx_1;
  double C36_idx_1;
  double C14_idx_1;
  double C45_idx_1;
  double C46_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double H23_idx_1;
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C12_idx_2;
  double C13_idx_2;
  double C23_idx_2;
  double C24_idx_2;
  double C25_idx_2;
  double C26_idx_2;
  double C34_idx_2;
  double C35_idx_2;
  double C36_idx_2;
  double C14_idx_2;
  double C45_idx_2;
  double C46_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C12_idx_3;
  double C13_idx_3;
  double C23_idx_3;
  double C24_idx_3;
  double C25_idx_3;
  double C26_idx_3;
  double C34_idx_3;
  double C35_idx_3;
  double C36_idx_3;
  double C14_idx_3;
  double C45_idx_3;
  double C46_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H34_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double b_elStorage[144];
  double Kenr[144];
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double Ke[144];
  double reshapes_f3_data[24];
  double reshapes_f4_data[24];
  double reshapes_f1[24];
  double reshapes_f2[24];
  double reshapes_f5[24];
  double reshapes_f6[24];
  double Ce[144];
  int reshapes_f1_data_tmp;
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fel_data[12];
  int tmp_size[1];
  boolean_T concMassFlag;

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  // declare type
  // declare type
  // declare type
  // options for Dean integrator
  // --------------------------------------------
  // setting for modal analysis flag
  // setting for initial reduced order model calculations
  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  SS22_data_idx_0 = 0.0;
  SS33_data_idx_0 = 0.0;
  SS22_data_idx_1 = 0.0;
  SS33_data_idx_1 = 0.0;
  SS22_data_idx_2 = 0.0;
  SS33_data_idx_2 = 0.0;
  SS22_data_idx_3 = 0.0;
  SS33_data_idx_3 = 0.0;

  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input_sweepAngle * 3.1415926535897931 / 180.0, input_coneAngle
                  * 3.1415926535897931 / 180.0, (input_rollAngle + 0.5 *
    (input_sectionProps_twist[0] + input_sectionProps_twist[1])) *
                  3.1415926535897931 / 180.0, lambda);
  std::memcpy(&b_data[0], &input_disp_data[0], 12U * sizeof(double));
  mtimes(lambda, b_data, dispLocal);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    K15_idx_1 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    K16_idx_1 = lambda[i + 24];
    a_temp[i] = (input_CN2H[i] * 0.0 + input_CN2H[i + 3] * 0.0) + input_CN2H[i +
      6] * 9.81;
    ODotel[i] = K15_idx_1 + K16_idx_1 * 0.0;
    Oel[i] = K15_idx_1 + K16_idx_1 * 3.2898681336964524;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  f_tmp = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input_xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &K15_idx_0);
    integrationFactor = K15_idx_0 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input_sectionProps_rhoA[0] + N_data[1] *
      input_sectionProps_rhoA[1];
    K15_idx_0 = p_N_x_data[0] * dispLocal[1] + p_N_x_data[1] * dispLocal[7];
    K16_idx_0 = p_N_x_data[0] * dispLocal[2] + p_N_x_data[1] * dispLocal[8];
    Faxial = (N_data[0] * input_sectionProps_EA[0] + N_data[1] *
              input_sectionProps_EA[1]) * (((p_N_x_data[0] * dispLocal[0] +
      p_N_x_data[1] * dispLocal[6]) + 0.5 * (K15_idx_0 * K15_idx_0)) + 0.5 *
      (K16_idx_0 * K16_idx_0));

    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    ycm = N_data[0] * input_sectionProps_ycm[0] + N_data[1] *
      input_sectionProps_ycm[1];
    zcm = N_data[0] * input_sectionProps_zcm[0] + N_data[1] *
      input_sectionProps_zcm[1];
    K15_idx_0 = N_data[0] * input_x_data[0] + N_data[1] * input_x_data[1];
    K16_idx_0 = N_data[0] * input_y_data[0] + N_data[1] * input_y_data[1];
    K56_idx_0 = N_data[0] * input_z_data[0] + N_data[1] * input_z_data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      K15_idx_1 = lambda[i + 12];
      K16_idx_1 = lambda[i + 24];
      posLocal[i] = (lambda[i] * K15_idx_0 + K15_idx_1 * K16_idx_0) + K16_idx_1 *
        K56_idx_0;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + K15_idx_1 * disMomentgp
                           [1]) + K16_idx_1 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // stress-stiffening/pre-stress calculations
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // calculate static aerodynamic load
    // distributed/body force load calculations
    c_f_tmp = (O2 * O2 + O3 * O3) * posLocal[0];
    d_f_tmp = O2dot * posLocal[2];
    e_f_tmp = O3dot * posLocal[1];
    f_f_tmp = f_tmp * posLocal[2];
    g_f_tmp = O1 * O2 * posLocal[1];
    b_f = rhoA * ((((c_f_tmp - g_f_tmp) - f_f_tmp) + e_f_tmp) - d_f_tmp) -
      disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = O1dot * posLocal[2];
    K16_idx_0 = O3dot * posLocal[0];
    c_f = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 * O3)
                    - posLocal[0] * O1 * O2) + K15_idx_0) - K16_idx_0) -
      disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    K56_idx_0 = O2dot * posLocal[0];
    b_f_tmp = O1dot * posLocal[1];
    d_f = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - f_tmp * posLocal[0]) -
                    O2 * O3 * posLocal[1]) + K56_idx_0) - b_f_tmp) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    f = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) - posLocal[1] *
                   (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                  (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                 (K56_idx_0 - b_f_tmp)) - zcm * (K15_idx_0 - K16_idx_0)) -
      disMomentgp[0];

    // This function is a general routine to calculate an element vector
    b_f_tmp = rhoA * zcm * ((((c_f_tmp - posLocal[1] * O1 * O2) - posLocal[2] *
      O1 * O3) - d_f_tmp) + e_f_tmp) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = rhoA * ycm * ((((f_f_tmp + g_f_tmp) - c_f_tmp) - e_f_tmp) +
      d_f_tmp) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    K16_idx_0 = Faxial * p_N_x_data[0];
    K56_idx_0 = K16_idx_0 * p_N_x_data[0] * integrationFactor;
    SS22_data_idx_0 += K56_idx_0;
    SS33_data_idx_0 += K56_idx_0;
    K56_idx_0 = K16_idx_0 * p_N_x_data[1] * integrationFactor;
    SS22_data_idx_2 += K56_idx_0;
    SS33_data_idx_2 += K56_idx_0;
    F1_data_idx_0 += b_f * N_data[0] * integrationFactor;
    F2_data_idx_0 += c_f * N_data[0] * integrationFactor;
    F3_data_idx_0 += d_f * N_data[0] * integrationFactor;
    F4_data_idx_0 += f * N_data[0] * integrationFactor;
    F5_data_idx_0 += b_f_tmp * N_data[0] * integrationFactor;
    F6_data_idx_0 += K15_idx_0 * N_data[0] * integrationFactor;
    K16_idx_0 = Faxial * p_N_x_data[1];
    K56_idx_0 = K16_idx_0 * p_N_x_data[0] * integrationFactor;
    SS22_data_idx_1 += K56_idx_0;
    SS33_data_idx_1 += K56_idx_0;
    K56_idx_0 = K16_idx_0 * p_N_x_data[1] * integrationFactor;
    SS22_data_idx_3 += K56_idx_0;
    SS33_data_idx_3 += K56_idx_0;
    F1_data_idx_1 += b_f * N_data[1] * integrationFactor;
    F2_data_idx_1 += c_f * N_data[1] * integrationFactor;
    F3_data_idx_1 += d_f * N_data[1] * integrationFactor;
    F4_data_idx_1 += f * N_data[1] * integrationFactor;
    F5_data_idx_1 += b_f_tmp * N_data[1] * integrationFactor;
    F6_data_idx_1 += K15_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  K15_idx_0 = Oel[0] * Oel[1];
  c_tmp = Oel[0] * Oel[0];
  b_c_tmp = c_tmp + Oel[2] * Oel[2];
  c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  rhoA = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  Faxial = elStorage->S23[0] * Oel[1] * Oel[2];
  O3 = elStorage->S25[0] * K15_idx_0;
  O2dot = elStorage->S26[0] * K15_idx_0;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  S36_idx_0 = elStorage->S36[0] * Oel[0] * Oel[2];
  ycm = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel[0] *
    Oel[1];
  O2 = elStorage->S24_1[0] * b_c_tmp + elStorage->S24_2[0] * Oel[1] * Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * c_tmp + elStorage->S34_2[0] * Oel[1] * Oel[2];
  S45_idx_0 = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel
    [0] * Oel[1];
  S46_idx_0 = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0] * Oel
    [0] * Oel[2];
  K15_idx_1 = elStorage->C12[0];
  C12_idx_0 = K15_idx_1 * Oel[2];
  K16_idx_1 = elStorage->C13[0];
  C13_idx_0 = K16_idx_1 * Oel[1];
  K16_idx_0 = elStorage->C23[0];
  C23_idx_0 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C24[0];
  C24_idx_0 = K56_idx_0 * Oel[0];
  b_f_tmp = elStorage->C25[0];
  C25_idx_0 = b_f_tmp * Oel[2];
  f = elStorage->C26[0];
  C26_idx_0 = f * Oel[2];
  c_f_tmp = elStorage->C34[0];
  C34_idx_0 = c_f_tmp * Oel[0];
  d_f_tmp = elStorage->C35[0];
  C35_idx_0 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C36[0];
  C36_idx_0 = e_f_tmp * Oel[1];
  f_f_tmp = elStorage->C14_1[0];
  g_f_tmp = elStorage->C14_2[0];
  C14_idx_0 = f_f_tmp * Oel[1] + g_f_tmp * Oel[2];
  b_f = elStorage->C45_1[0];
  c_f = elStorage->C45_2[0];
  C45_idx_0 = b_f * Oel[2] + c_f * Oel[1];
  d_f = elStorage->C46_1[0];
  zcm = elStorage->C46_2[0];
  C46_idx_0 = d_f * Oel[1] + zcm * Oel[2];
  integrationFactor = 0.5 * K15_idx_1 * ODotel[2];
  H13_idx_0 = 0.5 * K16_idx_1 * ODotel[1];
  H23_idx_0 = 0.5 * K16_idx_0 * ODotel[0];
  H24_idx_0 = 0.5 * K56_idx_0 * ODotel[0];
  H25_idx_0 = 0.5 * b_f_tmp * ODotel[2];
  H26_idx_0 = 0.5 * f * ODotel[2];
  H34_idx_0 = 0.5 * c_f_tmp * ODotel[0];
  H35_idx_0 = 0.5 * d_f_tmp * ODotel[1];
  H36_idx_0 = 0.5 * e_f_tmp * ODotel[1];
  H14_idx_0 = 0.5 * (f_f_tmp * ODotel[1] + g_f_tmp * ODotel[2]);
  H45_idx_0 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_0 = 0.5 * (d_f * ODotel[1] + zcm * ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * K15_idx_0;
  S26_idx_1 = elStorage->S26[1] * K15_idx_0;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * b_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * c_tmp + elStorage->S34_2[1] * Oel[1] * Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  K15_idx_1 = elStorage->C12[1];
  C12_idx_1 = K15_idx_1 * Oel[2];
  K16_idx_1 = elStorage->C13[1];
  C13_idx_1 = K16_idx_1 * Oel[1];
  K16_idx_0 = elStorage->C23[1];
  C23_idx_1 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C24[1];
  C24_idx_1 = K56_idx_0 * Oel[0];
  b_f_tmp = elStorage->C25[1];
  C25_idx_1 = b_f_tmp * Oel[2];
  f = elStorage->C26[1];
  C26_idx_1 = f * Oel[2];
  c_f_tmp = elStorage->C34[1];
  C34_idx_1 = c_f_tmp * Oel[0];
  d_f_tmp = elStorage->C35[1];
  C35_idx_1 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C36[1];
  C36_idx_1 = e_f_tmp * Oel[1];
  f_f_tmp = elStorage->C14_1[1];
  g_f_tmp = elStorage->C14_2[1];
  C14_idx_1 = f_f_tmp * Oel[1] + g_f_tmp * Oel[2];
  b_f = elStorage->C45_1[1];
  c_f = elStorage->C45_2[1];
  C45_idx_1 = b_f * Oel[2] + c_f * Oel[1];
  d_f = elStorage->C46_1[1];
  zcm = elStorage->C46_2[1];
  C46_idx_1 = d_f * Oel[1] + zcm * Oel[2];
  H12_idx_1 = 0.5 * K15_idx_1 * ODotel[2];
  H13_idx_1 = 0.5 * K16_idx_1 * ODotel[1];
  H23_idx_1 = 0.5 * K16_idx_0 * ODotel[0];
  H24_idx_1 = 0.5 * K56_idx_0 * ODotel[0];
  H25_idx_1 = 0.5 * b_f_tmp * ODotel[2];
  H26_idx_1 = 0.5 * f * ODotel[2];
  H34_idx_1 = 0.5 * c_f_tmp * ODotel[0];
  H35_idx_1 = 0.5 * d_f_tmp * ODotel[1];
  H36_idx_1 = 0.5 * e_f_tmp * ODotel[1];
  H14_idx_1 = 0.5 * (f_f_tmp * ODotel[1] + g_f_tmp * ODotel[2]);
  H45_idx_1 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_1 = 0.5 * (d_f * ODotel[1] + zcm * ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * K15_idx_0;
  S26_idx_2 = elStorage->S26[2] * K15_idx_0;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * b_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * c_tmp + elStorage->S34_2[2] * Oel[1] * Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  K15_idx_1 = elStorage->C12[2];
  C12_idx_2 = K15_idx_1 * Oel[2];
  K16_idx_1 = elStorage->C13[2];
  C13_idx_2 = K16_idx_1 * Oel[1];
  K16_idx_0 = elStorage->C23[2];
  C23_idx_2 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C24[2];
  C24_idx_2 = K56_idx_0 * Oel[0];
  b_f_tmp = elStorage->C25[2];
  C25_idx_2 = b_f_tmp * Oel[2];
  f = elStorage->C26[2];
  C26_idx_2 = f * Oel[2];
  c_f_tmp = elStorage->C34[2];
  C34_idx_2 = c_f_tmp * Oel[0];
  d_f_tmp = elStorage->C35[2];
  C35_idx_2 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C36[2];
  C36_idx_2 = e_f_tmp * Oel[1];
  f_f_tmp = elStorage->C14_1[2];
  g_f_tmp = elStorage->C14_2[2];
  C14_idx_2 = f_f_tmp * Oel[1] + g_f_tmp * Oel[2];
  b_f = elStorage->C45_1[2];
  c_f = elStorage->C45_2[2];
  C45_idx_2 = b_f * Oel[2] + c_f * Oel[1];
  d_f = elStorage->C46_1[2];
  zcm = elStorage->C46_2[2];
  C46_idx_2 = d_f * Oel[1] + zcm * Oel[2];
  H12_idx_2 = 0.5 * K15_idx_1 * ODotel[2];
  H13_idx_2 = 0.5 * K16_idx_1 * ODotel[1];
  H23_idx_2 = 0.5 * K16_idx_0 * ODotel[0];
  H24_idx_2 = 0.5 * K56_idx_0 * ODotel[0];
  H25_idx_2 = 0.5 * b_f_tmp * ODotel[2];
  H26_idx_2 = 0.5 * f * ODotel[2];
  H34_idx_2 = 0.5 * c_f_tmp * ODotel[0];
  H35_idx_2 = 0.5 * d_f_tmp * ODotel[1];
  H36_idx_2 = 0.5 * e_f_tmp * ODotel[1];
  H14_idx_2 = 0.5 * (f_f_tmp * ODotel[1] + g_f_tmp * ODotel[2]);
  H45_idx_2 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_2 = 0.5 * (d_f * ODotel[1] + zcm * ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * K15_idx_0;
  S26_idx_3 = elStorage->S26[3] * K15_idx_0;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * b_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * c_tmp + elStorage->S34_2[3] * Oel[1] * Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  K15_idx_1 = elStorage->C12[3];
  C12_idx_3 = K15_idx_1 * Oel[2];
  K16_idx_1 = elStorage->C13[3];
  C13_idx_3 = K16_idx_1 * Oel[1];
  K16_idx_0 = elStorage->C23[3];
  C23_idx_3 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C24[3];
  C24_idx_3 = K56_idx_0 * Oel[0];
  b_f_tmp = elStorage->C25[3];
  C25_idx_3 = b_f_tmp * Oel[2];
  f = elStorage->C26[3];
  C26_idx_3 = f * Oel[2];
  c_f_tmp = elStorage->C34[3];
  C34_idx_3 = c_f_tmp * Oel[0];
  d_f_tmp = elStorage->C35[3];
  C35_idx_3 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C36[3];
  C36_idx_3 = e_f_tmp * Oel[1];
  f_f_tmp = elStorage->C14_1[3];
  g_f_tmp = elStorage->C14_2[3];
  C14_idx_3 = f_f_tmp * Oel[1] + g_f_tmp * Oel[2];
  b_f = elStorage->C45_1[3];
  c_f = elStorage->C45_2[3];
  C45_idx_3 = b_f * Oel[2] + c_f * Oel[1];
  d_f = elStorage->C46_1[3];
  zcm = elStorage->C46_2[3];
  C46_idx_3 = d_f * Oel[1] + zcm * Oel[2];
  H12_idx_3 = 0.5 * K15_idx_1 * ODotel[2];
  H13_idx_3 = 0.5 * K16_idx_1 * ODotel[1];
  H23_idx_3 = 0.5 * K16_idx_0 * ODotel[0];
  H24_idx_3 = 0.5 * K56_idx_0 * ODotel[0];
  H25_idx_3 = 0.5 * b_f_tmp * ODotel[2];
  H26_idx_3 = 0.5 * f * ODotel[2];
  H34_idx_3 = 0.5 * c_f_tmp * ODotel[0];
  H35_idx_3 = 0.5 * d_f_tmp * ODotel[1];
  H36_idx_3 = 0.5 * e_f_tmp * ODotel[1];
  H14_idx_3 = 0.5 * (f_f_tmp * ODotel[1] + g_f_tmp * ODotel[2]);
  H45_idx_3 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  O3dot = 0.5 * (d_f * ODotel[1] + zcm * ODotel[2]);

  // compile stiffness matrix without rotational effects
  b_elStorage[0] = elStorage->K11[0];
  b_elStorage[24] = elStorage->K12[0];
  b_elStorage[48] = elStorage->K13[0];
  b_elStorage[72] = elStorage->K14[0];
  b_elStorage[96] = elStorage->K15[0];
  b_elStorage[120] = elStorage->K16[0];
  b_elStorage[2] = elStorage->K12[0];
  b_elStorage[26] = elStorage->K22[0];
  b_elStorage[50] = elStorage->K23[0];
  b_elStorage[74] = elStorage->K24[0];
  b_elStorage[98] = elStorage->K25[0];
  b_elStorage[122] = elStorage->K26[0];
  b_elStorage[4] = elStorage->K13[0];
  b_elStorage[28] = elStorage->K23[0];
  b_elStorage[52] = elStorage->K33[0];
  b_elStorage[76] = elStorage->K34[0];
  b_elStorage[100] = elStorage->K35[0];
  b_elStorage[124] = elStorage->K36[0];
  b_elStorage[6] = elStorage->K13[0];
  b_elStorage[30] = elStorage->K24[0];
  b_elStorage[54] = elStorage->K34[0];
  b_elStorage[78] = elStorage->K44[0];
  b_elStorage[102] = elStorage->K45[0];
  b_elStorage[126] = elStorage->K46[0];
  b_elStorage[8] = elStorage->K15[0];
  b_elStorage[32] = elStorage->K25[0];
  b_elStorage[56] = elStorage->K35[0];
  b_elStorage[80] = elStorage->K45[0];
  b_elStorage[104] = elStorage->K55[0];
  b_elStorage[128] = elStorage->K56[0];
  b_elStorage[10] = elStorage->K16[0];
  b_elStorage[34] = elStorage->K26[0];
  b_elStorage[58] = elStorage->K36[0];
  b_elStorage[82] = elStorage->K46[0];
  b_elStorage[106] = elStorage->K56[0];
  b_elStorage[130] = elStorage->K66[0];
  b_elStorage[1] = elStorage->K11[1];
  b_elStorage[25] = elStorage->K12[1];
  b_elStorage[49] = elStorage->K13[1];
  b_elStorage[73] = elStorage->K14[1];
  b_elStorage[97] = elStorage->K15[1];
  b_elStorage[121] = elStorage->K16[1];
  b_elStorage[3] = elStorage->K12[2];
  b_elStorage[27] = elStorage->K22[1];
  b_elStorage[51] = elStorage->K23[1];
  b_elStorage[75] = elStorage->K24[1];
  b_elStorage[99] = elStorage->K25[1];
  b_elStorage[123] = elStorage->K26[1];
  b_elStorage[5] = elStorage->K13[2];
  b_elStorage[29] = elStorage->K23[2];
  b_elStorage[53] = elStorage->K33[1];
  b_elStorage[77] = elStorage->K34[1];
  b_elStorage[101] = elStorage->K35[1];
  b_elStorage[125] = elStorage->K36[1];
  b_elStorage[7] = elStorage->K13[2];
  b_elStorage[31] = elStorage->K24[2];
  b_elStorage[55] = elStorage->K34[2];
  b_elStorage[79] = elStorage->K44[1];
  b_elStorage[103] = elStorage->K45[1];
  b_elStorage[127] = elStorage->K46[1];
  b_elStorage[9] = elStorage->K15[2];
  b_elStorage[33] = elStorage->K25[2];
  b_elStorage[57] = elStorage->K35[2];
  b_elStorage[81] = elStorage->K45[2];
  b_elStorage[105] = elStorage->K55[1];
  b_elStorage[129] = elStorage->K56[1];
  b_elStorage[11] = elStorage->K16[2];
  b_elStorage[35] = elStorage->K26[2];
  b_elStorage[59] = elStorage->K36[2];
  b_elStorage[83] = elStorage->K46[2];
  b_elStorage[107] = elStorage->K56[2];
  b_elStorage[131] = elStorage->K66[1];
  b_elStorage[12] = elStorage->K11[2];
  b_elStorage[36] = elStorage->K12[2];
  b_elStorage[60] = elStorage->K13[2];
  b_elStorage[84] = elStorage->K14[2];
  b_elStorage[108] = elStorage->K15[2];
  b_elStorage[132] = elStorage->K16[2];
  b_elStorage[14] = elStorage->K12[1];
  b_elStorage[38] = elStorage->K22[2];
  b_elStorage[62] = elStorage->K23[2];
  b_elStorage[86] = elStorage->K24[2];
  b_elStorage[110] = elStorage->K25[2];
  b_elStorage[134] = elStorage->K26[2];
  b_elStorage[16] = elStorage->K13[1];
  b_elStorage[40] = elStorage->K23[1];
  b_elStorage[64] = elStorage->K33[2];
  b_elStorage[88] = elStorage->K34[2];
  b_elStorage[112] = elStorage->K35[2];
  b_elStorage[136] = elStorage->K36[2];
  b_elStorage[18] = elStorage->K13[1];
  b_elStorage[42] = elStorage->K24[1];
  b_elStorage[66] = elStorage->K34[1];
  b_elStorage[90] = elStorage->K44[2];
  b_elStorage[114] = elStorage->K45[2];
  b_elStorage[138] = elStorage->K46[2];
  b_elStorage[20] = elStorage->K15[1];
  b_elStorage[44] = elStorage->K25[1];
  b_elStorage[68] = elStorage->K35[1];
  b_elStorage[92] = elStorage->K45[1];
  b_elStorage[116] = elStorage->K55[2];
  b_elStorage[140] = elStorage->K56[2];
  b_elStorage[22] = elStorage->K16[1];
  b_elStorage[46] = elStorage->K26[1];
  b_elStorage[70] = elStorage->K36[1];
  b_elStorage[94] = elStorage->K46[1];
  b_elStorage[118] = elStorage->K56[1];
  b_elStorage[142] = elStorage->K66[2];
  b_elStorage[13] = elStorage->K11[3];
  b_elStorage[37] = elStorage->K12[3];
  b_elStorage[61] = elStorage->K13[3];
  b_elStorage[85] = elStorage->K14[3];
  b_elStorage[109] = elStorage->K15[3];
  b_elStorage[133] = elStorage->K16[3];
  b_elStorage[15] = elStorage->K12[3];
  b_elStorage[39] = elStorage->K22[3];
  b_elStorage[63] = elStorage->K23[3];
  b_elStorage[87] = elStorage->K24[3];
  b_elStorage[111] = elStorage->K25[3];
  b_elStorage[135] = elStorage->K26[3];
  b_elStorage[17] = elStorage->K13[3];
  b_elStorage[41] = elStorage->K23[3];
  b_elStorage[65] = elStorage->K33[3];
  b_elStorage[89] = elStorage->K34[3];
  b_elStorage[113] = elStorage->K35[3];
  b_elStorage[137] = elStorage->K36[3];
  b_elStorage[19] = elStorage->K13[3];
  b_elStorage[43] = elStorage->K24[3];
  b_elStorage[67] = elStorage->K34[3];
  b_elStorage[91] = elStorage->K44[3];
  b_elStorage[115] = elStorage->K45[3];
  b_elStorage[139] = elStorage->K46[3];
  b_elStorage[21] = elStorage->K15[3];
  b_elStorage[45] = elStorage->K25[3];
  b_elStorage[69] = elStorage->K35[3];
  b_elStorage[93] = elStorage->K45[3];
  b_elStorage[117] = elStorage->K55[3];
  b_elStorage[141] = elStorage->K56[3];
  b_elStorage[23] = elStorage->K16[3];
  b_elStorage[47] = elStorage->K26[3];
  b_elStorage[71] = elStorage->K36[3];
  b_elStorage[95] = elStorage->K46[3];
  b_elStorage[119] = elStorage->K56[3];
  b_elStorage[143] = elStorage->K66[3];
  mapMatrixNonSym(b_elStorage, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  f_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * f_tmp;
  K16_idx_0 = elStorage->K16[0] + elStorage->S16[0] * f_tmp;
  K56_idx_0 = elStorage->K56[0] + elStorage->S56[0] * f_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * f_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * f_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * f_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * f_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * f_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * f_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * f_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * f_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * f_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  b_elStorage[0] = elStorage->K11[0] + elStorage->S11[0] * f_tmp;
  b_elStorage[24] = (elStorage->K12[0] + rhoA) + integrationFactor;
  b_elStorage[48] = (elStorage->K13[0] + S13_idx_0) + H13_idx_0;
  O1dot = elStorage->K14[0] + ycm;
  b_elStorage[72] = O1dot + H14_idx_0;
  b_elStorage[96] = K15_idx_0;
  b_elStorage[120] = K16_idx_0;
  b_elStorage[2] = (elStorage->K12[0] + rhoA) - integrationFactor;
  b_elStorage[26] = (elStorage->K22[0] + elStorage->S22[0] * b_c_tmp) +
    SS22_data_idx_0;
  O1 = elStorage->K23[0] + Faxial;
  b_elStorage[50] = O1 + H23_idx_0;
  integrationFactor = elStorage->K24[0] + O2;
  b_elStorage[74] = integrationFactor + H24_idx_0;
  Faxial = elStorage->K25[0] + O3;
  b_elStorage[98] = Faxial + H25_idx_0;
  rhoA = elStorage->K26[0] + O2dot;
  b_elStorage[122] = rhoA + H26_idx_0;
  b_elStorage[4] = (elStorage->K13[0] + S13_idx_0) - H13_idx_0;
  b_elStorage[28] = O1 - H23_idx_0;
  b_elStorage[52] = (elStorage->K33[0] + elStorage->S33[0] * c_tmp) +
    SS33_data_idx_0;
  O1 = elStorage->K34[0] + S34_idx_0;
  b_elStorage[76] = O1 + H34_idx_0;
  ycm = elStorage->K35[0] + S35_idx_0;
  b_elStorage[100] = ycm + H35_idx_0;
  zcm = elStorage->K36[0] + S36_idx_0;
  b_elStorage[124] = zcm + H36_idx_0;
  b_elStorage[6] = O1dot - H14_idx_0;
  b_elStorage[30] = integrationFactor - H24_idx_0;
  b_elStorage[54] = O1 - H34_idx_0;
  b_elStorage[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * b_c_tmp +
    elStorage->S44_2[0] * c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  O1dot = elStorage->K45[0] + S45_idx_0;
  b_elStorage[102] = O1dot + H45_idx_0;
  O1 = elStorage->K46[0] + S46_idx_0;
  b_elStorage[126] = O1 + H46_idx_0;
  b_elStorage[8] = K15_idx_0;
  b_elStorage[32] = Faxial - H25_idx_0;
  b_elStorage[56] = ycm - H35_idx_0;
  b_elStorage[80] = O1dot - H45_idx_0;
  b_elStorage[104] = elStorage->K55[0] + elStorage->S55[0] * f_tmp;
  b_elStorage[128] = K56_idx_0;
  b_elStorage[10] = K16_idx_0;
  b_elStorage[34] = rhoA - H26_idx_0;
  b_elStorage[58] = zcm - H36_idx_0;
  b_elStorage[82] = O1 - H46_idx_0;
  b_elStorage[106] = K56_idx_0;
  b_elStorage[130] = elStorage->K66[0] + elStorage->S66[0] * f_tmp;
  b_elStorage[1] = elStorage->K11[1] + elStorage->S11[1] * f_tmp;
  b_elStorage[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  b_elStorage[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  O1dot = elStorage->K14[1] + S14_idx_1;
  b_elStorage[73] = O1dot + H14_idx_1;
  b_elStorage[97] = K15_idx_1;
  b_elStorage[121] = K16_idx_1;
  b_elStorage[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  b_elStorage[27] = (elStorage->K22[1] + elStorage->S22[1] * b_c_tmp) +
    SS22_data_idx_1;
  O1 = elStorage->K23[1] + S23_idx_1;
  b_elStorage[51] = O1 + H23_idx_1;
  integrationFactor = elStorage->K24[1] + S24_idx_1;
  b_elStorage[75] = integrationFactor + H24_idx_1;
  Faxial = elStorage->K25[1] + S25_idx_1;
  b_elStorage[99] = Faxial + H25_idx_1;
  rhoA = elStorage->K26[1] + S26_idx_1;
  b_elStorage[123] = rhoA + H26_idx_1;
  b_elStorage[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  ycm = elStorage->K23[2] + S23_idx_2;
  b_elStorage[29] = ycm - H23_idx_2;
  b_elStorage[53] = (elStorage->K33[1] + elStorage->S33[1] * c_tmp) +
    SS33_data_idx_1;
  zcm = elStorage->K34[1] + S34_idx_1;
  b_elStorage[77] = zcm + H34_idx_1;
  d_f = elStorage->K35[1] + S35_idx_1;
  b_elStorage[101] = d_f + H35_idx_1;
  c_f = elStorage->K36[1] + S36_idx_1;
  b_elStorage[125] = c_f + H36_idx_1;
  b_f = elStorage->K14[2] + S14_idx_2;
  b_elStorage[7] = b_f - H14_idx_2;
  g_f_tmp = elStorage->K24[2] + S24_idx_2;
  b_elStorage[31] = g_f_tmp - H24_idx_2;
  f_f_tmp = elStorage->K34[2] + S34_idx_2;
  b_elStorage[55] = f_f_tmp - H34_idx_2;
  b_elStorage[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * b_c_tmp +
    elStorage->S44_2[1] * c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  e_f_tmp = elStorage->K45[1] + S45_idx_1;
  b_elStorage[103] = e_f_tmp + H45_idx_1;
  d_f_tmp = elStorage->K46[1] + S46_idx_1;
  b_elStorage[127] = d_f_tmp + H46_idx_1;
  b_elStorage[9] = K15_idx_2;
  c_f_tmp = elStorage->K25[2] + S25_idx_2;
  b_elStorage[33] = c_f_tmp - H25_idx_2;
  f = elStorage->K35[2] + S35_idx_2;
  b_elStorage[57] = f - H35_idx_2;
  b_f_tmp = elStorage->K45[2] + S45_idx_2;
  b_elStorage[81] = b_f_tmp - H45_idx_2;
  b_elStorage[105] = elStorage->K55[1] + elStorage->S55[1] * f_tmp;
  b_elStorage[129] = K56_idx_1;
  b_elStorage[11] = K16_idx_2;
  K56_idx_0 = elStorage->K26[2] + S26_idx_2;
  b_elStorage[35] = K56_idx_0 - H26_idx_2;
  K16_idx_0 = elStorage->K36[2] + S36_idx_2;
  b_elStorage[59] = K16_idx_0 - H36_idx_2;
  K15_idx_0 = elStorage->K46[2] + S46_idx_2;
  b_elStorage[83] = K15_idx_0 - H46_idx_2;
  b_elStorage[107] = K56_idx_2;
  b_elStorage[131] = elStorage->K66[1] + elStorage->S66[1] * f_tmp;
  b_elStorage[12] = elStorage->K11[2] + elStorage->S11[2] * f_tmp;
  b_elStorage[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  b_elStorage[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  b_elStorage[84] = b_f + H14_idx_2;
  b_elStorage[108] = K15_idx_2;
  b_elStorage[132] = K16_idx_2;
  b_elStorage[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  b_elStorage[38] = (elStorage->K22[2] + elStorage->S22[2] * b_c_tmp) +
    SS22_data_idx_2;
  b_elStorage[62] = ycm + H23_idx_2;
  b_elStorage[86] = g_f_tmp + H24_idx_2;
  b_elStorage[110] = c_f_tmp + H25_idx_2;
  b_elStorage[134] = K56_idx_0 + H26_idx_2;
  b_elStorage[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  b_elStorage[40] = O1 - H23_idx_1;
  b_elStorage[64] = (elStorage->K33[2] + elStorage->S33[2] * c_tmp) +
    SS33_data_idx_2;
  b_elStorage[88] = f_f_tmp + H34_idx_2;
  b_elStorage[112] = f + H35_idx_2;
  b_elStorage[136] = K16_idx_0 + H36_idx_2;
  b_elStorage[18] = O1dot - H14_idx_1;
  b_elStorage[42] = integrationFactor - H24_idx_1;
  b_elStorage[66] = zcm - H34_idx_1;
  b_elStorage[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * b_c_tmp +
    elStorage->S44_2[2] * c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  b_elStorage[114] = b_f_tmp + H45_idx_2;
  b_elStorage[138] = K15_idx_0 + H46_idx_2;
  b_elStorage[20] = K15_idx_1;
  b_elStorage[44] = Faxial - H25_idx_1;
  b_elStorage[68] = d_f - H35_idx_1;
  b_elStorage[92] = e_f_tmp - H45_idx_1;
  b_elStorage[116] = elStorage->K55[2] + elStorage->S55[2] * f_tmp;
  b_elStorage[140] = K56_idx_2;
  b_elStorage[22] = K16_idx_1;
  b_elStorage[46] = rhoA - H26_idx_1;
  b_elStorage[70] = c_f - H36_idx_1;
  b_elStorage[94] = d_f_tmp - H46_idx_1;
  b_elStorage[118] = K56_idx_1;
  b_elStorage[142] = elStorage->K66[2] + elStorage->S66[2] * f_tmp;
  b_elStorage[13] = elStorage->K11[3] + elStorage->S11[3] * f_tmp;
  b_elStorage[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  b_elStorage[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  O1dot = elStorage->K14[3] + S14_idx_3;
  b_elStorage[85] = O1dot + H14_idx_3;
  b_elStorage[109] = K15_idx_3;
  b_elStorage[133] = K16_idx_3;
  b_elStorage[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  b_elStorage[39] = (elStorage->K22[3] + elStorage->S22[3] * b_c_tmp) +
    SS22_data_idx_3;
  O1 = elStorage->K23[3] + S23_idx_3;
  b_elStorage[63] = O1 + H23_idx_3;
  integrationFactor = elStorage->K24[3] + S24_idx_3;
  b_elStorage[87] = integrationFactor + H24_idx_3;
  Faxial = elStorage->K25[3] + S25_idx_3;
  b_elStorage[111] = Faxial + H25_idx_3;
  rhoA = elStorage->K26[3] + S26_idx_3;
  b_elStorage[135] = rhoA + H26_idx_3;
  b_elStorage[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  b_elStorage[41] = O1 - H23_idx_3;
  b_elStorage[65] = (elStorage->K33[3] + elStorage->S33[3] * c_tmp) +
    SS33_data_idx_3;
  O1 = elStorage->K34[3] + S34_idx_3;
  b_elStorage[89] = O1 + H34_idx_3;
  ycm = elStorage->K35[3] + S35_idx_3;
  b_elStorage[113] = ycm + H35_idx_3;
  zcm = elStorage->K36[3] + S36_idx_3;
  b_elStorage[137] = zcm + H36_idx_3;
  b_elStorage[19] = O1dot - H14_idx_3;
  b_elStorage[43] = integrationFactor - H24_idx_3;
  b_elStorage[67] = O1 - H34_idx_3;
  b_elStorage[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * b_c_tmp +
    elStorage->S44_2[3] * c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  O1dot = elStorage->K45[3] + S45_idx_3;
  b_elStorage[115] = O1dot + H45_idx_3;
  O1 = elStorage->K46[3] + S46_idx_3;
  b_elStorage[139] = O1 + O3dot;
  b_elStorage[21] = K15_idx_3;
  b_elStorage[45] = Faxial - H25_idx_3;
  b_elStorage[69] = ycm - H35_idx_3;
  b_elStorage[93] = O1dot - H45_idx_3;
  b_elStorage[117] = elStorage->K55[3] + elStorage->S55[3] * f_tmp;
  b_elStorage[141] = K56_idx_3;
  b_elStorage[23] = K16_idx_3;
  b_elStorage[47] = rhoA - H26_idx_3;
  b_elStorage[71] = zcm - H36_idx_3;
  b_elStorage[95] = O1 - O3dot;
  b_elStorage[119] = K56_idx_3;
  b_elStorage[143] = elStorage->K66[3] + elStorage->S66[3] * f_tmp;
  mapMatrixNonSym(b_elStorage, Ke);

  //  Declare type
  // compile Coriolis/damping matrix
  reshapes_f3_data[0] = -C13_idx_0;
  reshapes_f3_data[4] = -C23_idx_0;
  reshapes_f3_data[8] = 0.0;
  reshapes_f3_data[12] = C34_idx_0;
  reshapes_f3_data[16] = C35_idx_0;
  reshapes_f3_data[20] = C36_idx_0;
  reshapes_f4_data[0] = -C14_idx_0;
  reshapes_f4_data[4] = -C24_idx_0;
  reshapes_f4_data[8] = -C34_idx_0;
  reshapes_f4_data[12] = 0.0;
  reshapes_f4_data[16] = C45_idx_0;
  reshapes_f4_data[20] = C46_idx_0;
  reshapes_f1[0] = 0.0;
  reshapes_f1[4] = C12_idx_0;
  reshapes_f1[8] = C13_idx_0;
  reshapes_f1[12] = C14_idx_0;
  reshapes_f1[16] = 0.0;
  reshapes_f1[20] = 0.0;
  reshapes_f2[0] = -C12_idx_0;
  reshapes_f2[4] = 0.0;
  reshapes_f2[8] = C23_idx_0;
  reshapes_f2[12] = C24_idx_0;
  reshapes_f2[16] = C25_idx_0;
  reshapes_f2[20] = C26_idx_0;
  reshapes_f5[0] = 0.0;
  reshapes_f5[4] = -C25_idx_0;
  reshapes_f5[8] = -C35_idx_0;
  reshapes_f5[12] = -C45_idx_0;
  reshapes_f5[16] = 0.0;
  reshapes_f5[20] = 0.0;
  reshapes_f6[0] = 0.0;
  reshapes_f6[4] = -C26_idx_0;
  reshapes_f6[8] = -C36_idx_0;
  reshapes_f6[12] = -C46_idx_0;
  reshapes_f6[16] = 0.0;
  reshapes_f6[20] = 0.0;
  reshapes_f3_data[1] = -C13_idx_2;
  reshapes_f3_data[5] = -C23_idx_2;
  reshapes_f3_data[9] = 0.0;
  reshapes_f3_data[13] = C34_idx_1;
  reshapes_f3_data[17] = C35_idx_1;
  reshapes_f3_data[21] = C36_idx_1;
  reshapes_f4_data[1] = -C14_idx_2;
  reshapes_f4_data[5] = -C24_idx_2;
  reshapes_f4_data[9] = -C34_idx_2;
  reshapes_f4_data[13] = 0.0;
  reshapes_f4_data[17] = C45_idx_1;
  reshapes_f4_data[21] = C46_idx_1;
  reshapes_f1[1] = 0.0;
  reshapes_f1[5] = C12_idx_1;
  reshapes_f1[9] = C13_idx_1;
  reshapes_f1[13] = C14_idx_1;
  reshapes_f1[17] = 0.0;
  reshapes_f1[21] = 0.0;
  reshapes_f2[1] = -C12_idx_2;
  reshapes_f2[5] = 0.0;
  reshapes_f2[9] = C23_idx_1;
  reshapes_f2[13] = C24_idx_1;
  reshapes_f2[17] = C25_idx_1;
  reshapes_f2[21] = C26_idx_1;
  reshapes_f5[1] = 0.0;
  reshapes_f5[5] = -C25_idx_2;
  reshapes_f5[9] = -C35_idx_2;
  reshapes_f5[13] = -C45_idx_2;
  reshapes_f5[17] = 0.0;
  reshapes_f5[21] = 0.0;
  reshapes_f6[1] = 0.0;
  reshapes_f6[5] = -C26_idx_2;
  reshapes_f6[9] = -C36_idx_2;
  reshapes_f6[13] = -C46_idx_2;
  reshapes_f6[17] = 0.0;
  reshapes_f6[21] = 0.0;
  reshapes_f3_data[2] = -C13_idx_1;
  reshapes_f3_data[6] = -C23_idx_1;
  reshapes_f3_data[10] = 0.0;
  reshapes_f3_data[14] = C34_idx_2;
  reshapes_f3_data[18] = C35_idx_2;
  reshapes_f3_data[22] = C36_idx_2;
  reshapes_f4_data[2] = -C14_idx_1;
  reshapes_f4_data[6] = -C24_idx_1;
  reshapes_f4_data[10] = -C34_idx_1;
  reshapes_f4_data[14] = 0.0;
  reshapes_f4_data[18] = C45_idx_2;
  reshapes_f4_data[22] = C46_idx_2;
  reshapes_f1[2] = 0.0;
  reshapes_f1[6] = C12_idx_2;
  reshapes_f1[10] = C13_idx_2;
  reshapes_f1[14] = C14_idx_2;
  reshapes_f1[18] = 0.0;
  reshapes_f1[22] = 0.0;
  reshapes_f2[2] = -C12_idx_1;
  reshapes_f2[6] = 0.0;
  reshapes_f2[10] = C23_idx_2;
  reshapes_f2[14] = C24_idx_2;
  reshapes_f2[18] = C25_idx_2;
  reshapes_f2[22] = C26_idx_2;
  reshapes_f5[2] = 0.0;
  reshapes_f5[6] = -C25_idx_1;
  reshapes_f5[10] = -C35_idx_1;
  reshapes_f5[14] = -C45_idx_1;
  reshapes_f5[18] = 0.0;
  reshapes_f5[22] = 0.0;
  reshapes_f6[2] = 0.0;
  reshapes_f6[6] = -C26_idx_1;
  reshapes_f6[10] = -C36_idx_1;
  reshapes_f6[14] = -C46_idx_1;
  reshapes_f6[18] = 0.0;
  reshapes_f6[22] = 0.0;
  reshapes_f3_data[3] = -C13_idx_3;
  reshapes_f3_data[7] = -C23_idx_3;
  reshapes_f3_data[11] = 0.0;
  reshapes_f3_data[15] = C34_idx_3;
  reshapes_f3_data[19] = C35_idx_3;
  reshapes_f3_data[23] = C36_idx_3;
  reshapes_f4_data[3] = -C14_idx_3;
  reshapes_f4_data[7] = -C24_idx_3;
  reshapes_f4_data[11] = -C34_idx_3;
  reshapes_f4_data[15] = 0.0;
  reshapes_f4_data[19] = C45_idx_3;
  reshapes_f4_data[23] = C46_idx_3;
  reshapes_f1[3] = 0.0;
  reshapes_f1[7] = C12_idx_3;
  reshapes_f1[11] = C13_idx_3;
  reshapes_f1[15] = C14_idx_3;
  reshapes_f1[19] = 0.0;
  reshapes_f1[23] = 0.0;
  reshapes_f2[3] = -C12_idx_3;
  reshapes_f2[7] = 0.0;
  reshapes_f2[11] = C23_idx_3;
  reshapes_f2[15] = C24_idx_3;
  reshapes_f2[19] = C25_idx_3;
  reshapes_f2[23] = C26_idx_3;
  reshapes_f5[3] = 0.0;
  reshapes_f5[7] = -C25_idx_3;
  reshapes_f5[11] = -C35_idx_3;
  reshapes_f5[15] = -C45_idx_3;
  reshapes_f5[19] = 0.0;
  reshapes_f5[23] = 0.0;
  reshapes_f6[3] = 0.0;
  reshapes_f6[7] = -C26_idx_3;
  reshapes_f6[11] = -C36_idx_3;
  reshapes_f6[15] = -C46_idx_3;
  reshapes_f6[19] = 0.0;
  reshapes_f6[23] = 0.0;
  for (i = 0; i < 12; i++) {
    b_i = i << 1;
    b_elStorage[12 * i] = reshapes_f1[b_i];
    b_elStorage[12 * i + 2] = reshapes_f2[b_i];
    b_elStorage[12 * i + 4] = reshapes_f3_data[2 * i];
    b_elStorage[12 * i + 6] = reshapes_f4_data[2 * i];
    b_elStorage[12 * i + 8] = reshapes_f5[b_i];
    b_elStorage[12 * i + 10] = reshapes_f6[b_i];
    b_i++;
    b_elStorage[12 * i + 1] = reshapes_f1[b_i];
    b_elStorage[12 * i + 3] = reshapes_f2[b_i];
    reshapes_f1_data_tmp = 2 * i + 1;
    b_elStorage[12 * i + 5] = reshapes_f3_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 7] = reshapes_f4_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 9] = reshapes_f5[b_i];
    b_elStorage[12 * i + 11] = reshapes_f6[b_i];
  }

  b_mapMatrixNonSym(b_elStorage, Ce);

  // compile mass matrix
  b_elStorage[0] = elStorage->M11[0];
  b_elStorage[24] = 0.0;
  b_elStorage[48] = 0.0;
  b_elStorage[72] = 0.0;
  b_elStorage[96] = elStorage->M15[0];
  b_elStorage[120] = elStorage->M16[0];
  b_elStorage[2] = 0.0;
  b_elStorage[26] = elStorage->M22[0];
  b_elStorage[50] = 0.0;
  b_elStorage[74] = elStorage->M24[0];
  b_elStorage[98] = 0.0;
  b_elStorage[122] = 0.0;
  b_elStorage[4] = 0.0;
  b_elStorage[28] = 0.0;
  b_elStorage[52] = elStorage->M33[0];
  b_elStorage[76] = elStorage->M34[0];
  b_elStorage[100] = 0.0;
  b_elStorage[124] = 0.0;
  b_elStorage[6] = 0.0;
  b_elStorage[30] = elStorage->M24[0];
  b_elStorage[54] = elStorage->M34[0];
  b_elStorage[78] = elStorage->M44[0];
  b_elStorage[102] = 0.0;
  b_elStorage[126] = 0.0;
  b_elStorage[8] = elStorage->M15[0];
  b_elStorage[32] = 0.0;
  b_elStorage[56] = 0.0;
  b_elStorage[80] = 0.0;
  b_elStorage[104] = elStorage->M55[0];
  b_elStorage[128] = elStorage->M56[0];
  b_elStorage[10] = elStorage->M16[0];
  b_elStorage[34] = 0.0;
  b_elStorage[58] = 0.0;
  b_elStorage[82] = 0.0;
  b_elStorage[106] = elStorage->M56[0];
  b_elStorage[130] = elStorage->M66[0];
  b_elStorage[1] = elStorage->M11[1];
  b_elStorage[25] = 0.0;
  b_elStorage[49] = 0.0;
  b_elStorage[73] = 0.0;
  b_elStorage[97] = elStorage->M15[1];
  b_elStorage[121] = elStorage->M16[1];
  b_elStorage[3] = 0.0;
  b_elStorage[27] = elStorage->M22[1];
  b_elStorage[51] = 0.0;
  b_elStorage[75] = elStorage->M24[1];
  b_elStorage[99] = 0.0;
  b_elStorage[123] = 0.0;
  b_elStorage[5] = 0.0;
  b_elStorage[29] = 0.0;
  b_elStorage[53] = elStorage->M33[1];
  b_elStorage[77] = elStorage->M34[1];
  b_elStorage[101] = 0.0;
  b_elStorage[125] = 0.0;
  b_elStorage[7] = 0.0;
  b_elStorage[31] = elStorage->M24[2];
  b_elStorage[55] = elStorage->M34[2];
  b_elStorage[79] = elStorage->M44[1];
  b_elStorage[103] = 0.0;
  b_elStorage[127] = 0.0;
  b_elStorage[9] = elStorage->M15[2];
  b_elStorage[33] = 0.0;
  b_elStorage[57] = 0.0;
  b_elStorage[81] = 0.0;
  b_elStorage[105] = elStorage->M55[1];
  b_elStorage[129] = elStorage->M56[1];
  b_elStorage[11] = elStorage->M16[2];
  b_elStorage[35] = 0.0;
  b_elStorage[59] = 0.0;
  b_elStorage[83] = 0.0;
  b_elStorage[107] = elStorage->M56[2];
  b_elStorage[131] = elStorage->M66[1];
  b_elStorage[12] = elStorage->M11[2];
  b_elStorage[36] = 0.0;
  b_elStorage[60] = 0.0;
  b_elStorage[84] = 0.0;
  b_elStorage[108] = elStorage->M15[2];
  b_elStorage[132] = elStorage->M16[2];
  b_elStorage[14] = 0.0;
  b_elStorage[38] = elStorage->M22[2];
  b_elStorage[62] = 0.0;
  b_elStorage[86] = elStorage->M24[2];
  b_elStorage[110] = 0.0;
  b_elStorage[134] = 0.0;
  b_elStorage[16] = 0.0;
  b_elStorage[40] = 0.0;
  b_elStorage[64] = elStorage->M33[2];
  b_elStorage[88] = elStorage->M34[2];
  b_elStorage[112] = 0.0;
  b_elStorage[136] = 0.0;
  b_elStorage[18] = 0.0;
  b_elStorage[42] = elStorage->M24[1];
  b_elStorage[66] = elStorage->M34[1];
  b_elStorage[90] = elStorage->M44[2];
  b_elStorage[114] = 0.0;
  b_elStorage[138] = 0.0;
  b_elStorage[20] = elStorage->M15[1];
  b_elStorage[44] = 0.0;
  b_elStorage[68] = 0.0;
  b_elStorage[92] = 0.0;
  b_elStorage[116] = elStorage->M55[2];
  b_elStorage[140] = elStorage->M56[2];
  b_elStorage[22] = elStorage->M16[1];
  b_elStorage[46] = 0.0;
  b_elStorage[70] = 0.0;
  b_elStorage[94] = 0.0;
  b_elStorage[118] = elStorage->M56[1];
  b_elStorage[142] = elStorage->M66[2];
  b_elStorage[13] = elStorage->M11[3];
  b_elStorage[37] = 0.0;
  b_elStorage[61] = 0.0;
  b_elStorage[85] = 0.0;
  b_elStorage[109] = elStorage->M15[3];
  b_elStorage[133] = elStorage->M16[3];
  b_elStorage[15] = 0.0;
  b_elStorage[39] = elStorage->M22[3];
  b_elStorage[63] = 0.0;
  b_elStorage[87] = elStorage->M24[3];
  b_elStorage[111] = 0.0;
  b_elStorage[135] = 0.0;
  b_elStorage[17] = 0.0;
  b_elStorage[41] = 0.0;
  b_elStorage[65] = elStorage->M33[3];
  b_elStorage[89] = elStorage->M34[3];
  b_elStorage[113] = 0.0;
  b_elStorage[137] = 0.0;
  b_elStorage[19] = 0.0;
  b_elStorage[43] = elStorage->M24[3];
  b_elStorage[67] = elStorage->M34[3];
  b_elStorage[91] = elStorage->M44[3];
  b_elStorage[115] = 0.0;
  b_elStorage[139] = 0.0;
  b_elStorage[21] = elStorage->M15[3];
  b_elStorage[45] = 0.0;
  b_elStorage[69] = 0.0;
  b_elStorage[93] = 0.0;
  b_elStorage[117] = elStorage->M55[3];
  b_elStorage[141] = elStorage->M56[3];
  b_elStorage[23] = elStorage->M16[3];
  b_elStorage[47] = 0.0;
  b_elStorage[71] = 0.0;
  b_elStorage[95] = 0.0;
  b_elStorage[119] = elStorage->M56[3];
  b_elStorage[143] = elStorage->M66[3];
  mapMatrixNonSym(b_elStorage, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input_RayleighAlpha * Kenr[i] + input_RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      b_elStorage[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(b_elStorage, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ke,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Ke);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fel_data[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fel_data[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fel_data,
                  dispLocal);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input_concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input_concMass[0];
    Me[13] += input_concMass[0];
    Me[26] += input_concMass[0];
    Me[39] += input_concMass[1];
    Me[52] += input_concMass[2];
    Me[65] += input_concMass[3];
    Me[78] += input_concMass[4];
    Me[91] += input_concMass[4];
    Me[104] += input_concMass[4];
    Me[117] += input_concMass[5];
    Me[130] += input_concMass[6];
    Me[143] += input_concMass[7];

    // modify Ce for concentrated mass
    K15_idx_0 = 2.0 * input_concMass[0] * 3.2898681336964524;
    Ce[12] -= K15_idx_0;
    Ce[1] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[0] * 0.0;
    Ce[24] += K15_idx_0;
    Ce[2] -= K15_idx_0;
    Ce[25] -= K15_idx_0;
    Ce[14] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[4] * 3.2898681336964524;
    Ce[90] -= K15_idx_0;
    Ce[79] += K15_idx_0;
    K15_idx_0 = 2.0 * input_concMass[4] * 0.0;
    Ce[102] += K15_idx_0;
    Ce[80] -= K15_idx_0;
    Ce[103] -= K15_idx_0;
    Ce[92] += K15_idx_0;

    // modify Ke for concentrated mass
    Ke[0] -= input_concMass[0] * 10.823232337111378;
    K15_idx_0 = input_concMass[0] * 0.0 * 0.0;
    Ke[12] = (Ke[12] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[1] = (Ke[1] + K15_idx_0) + input_concMass[0] * 0.0;
    K15_idx_0 = input_concMass[0] * 0.0 * 3.2898681336964524;
    Ke[24] = (Ke[24] + K15_idx_0) + input_concMass[0] * 0.0;
    Ke[2] = (Ke[2] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[25] = (Ke[25] + K15_idx_0) - input_concMass[0] * 0.0;
    Ke[14] = (Ke[14] + K15_idx_0) + input_concMass[0] * 0.0;
    Ke[13] -= input_concMass[0] * 10.823232337111378;
    Ke[26] -= input_concMass[0] * 0.0;
    Ke[78] -= input_concMass[4] * 10.823232337111378;
    K15_idx_0 = input_concMass[4] * 0.0 * 0.0;
    Ke[90] = (Ke[90] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[79] = (Ke[79] + K15_idx_0) + input_concMass[4] * 0.0;
    K15_idx_0 = input_concMass[4] * 0.0 * 3.2898681336964524;
    Ke[102] = (Ke[102] + K15_idx_0) + input_concMass[4] * 0.0;
    Ke[80] = (Ke[80] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[103] = (Ke[103] + K15_idx_0) - input_concMass[4] * 0.0;
    Ke[92] = (Ke[92] + K15_idx_0) + input_concMass[4] * 0.0;
    Ke[91] -= input_concMass[4] * 10.823232337111378;
    Ke[104] -= input_concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    dispLocal[0] = ((dispLocal[0] + input_concMass[0] * ((input_x_data[0] *
      10.823232337111378 - 0.0 * input_y_data[0]) - 0.0 * input_z_data[0])) +
                    input_concMass[0] * (input_y_data[0] * 0.0 - input_z_data[0]
      * 0.0)) - input_concMass[0] * a_temp[0];
    K15_idx_0 = input_z_data[0] * 0.0 - input_x_data[0] * 0.0;
    dispLocal[1] = ((dispLocal[1] + input_concMass[0] * ((input_y_data[0] *
      10.823232337111378 - 0.0 * input_z_data[0]) - 0.0 * input_x_data[0])) +
                    input_concMass[0] * K15_idx_0) - input_concMass[0] * a_temp
      [1];
    dispLocal[2] = ((dispLocal[2] + input_concMass[0] * (K15_idx_0 - 0.0 *
      input_y_data[0])) + input_concMass[0] * (input_x_data[0] * 0.0 -
      input_y_data[0] * 0.0)) - input_concMass[0] * a_temp[2];
    dispLocal[6] = ((dispLocal[6] + input_concMass[4] * ((input_x_data[1] *
      10.823232337111378 - 0.0 * input_y_data[1]) - 0.0 * input_z_data[1])) +
                    input_concMass[4] * (input_y_data[1] * 0.0 - input_z_data[1]
      * 0.0)) - input_concMass[4] * a_temp[0];
    K15_idx_0 = input_z_data[1] * 0.0 - input_x_data[1] * 0.0;
    dispLocal[7] = ((dispLocal[7] + input_concMass[4] * ((input_y_data[1] *
      10.823232337111378 - 0.0 * input_z_data[1]) - 0.0 * input_x_data[1])) +
                    input_concMass[4] * K15_idx_0) - input_concMass[4] * a_temp
      [1];
    dispLocal[8] = ((dispLocal[8] + input_concMass[4] * (K15_idx_0 - 0.0 *
      input_y_data[1])) + input_concMass[4] * (input_x_data[1] * 0.0 -
      input_y_data[1] * 0.0)) - input_concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Ke[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &dispLocal[0], 12U * sizeof(double));
  output->Me.size[0] = 12;
  output->Me.size[1] = 12;
  output->Ce.size[0] = 12;
  output->Ce.size[1] = 12;
  std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
  std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));

  // ------------------------------------------
}

//
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
// Arguments    : const n_struct_T *input
//                const f_struct_T *elStorage
//                o_struct_T *output
// Return Type  : void
//
void calculateTimoshenkoElementNL(const n_struct_T *input, const f_struct_T
  *elStorage, o_struct_T *output)
{
  int dispdot_size_idx_1;
  double dispdot_data[12];
  int dispddot_size_idx_1;
  double dispddot_data[12];
  double F1_data_idx_0;
  double F3_data_idx_0;
  double F2_data_idx_0;
  double F4_data_idx_0;
  double F5_data_idx_0;
  double F6_data_idx_0;
  double F1_data_idx_1;
  double F3_data_idx_1;
  double F2_data_idx_1;
  double F4_data_idx_1;
  double F5_data_idx_1;
  double F6_data_idx_1;
  double Omega;
  double OmegaDot;
  double lambda[144];
  int i;
  double O1;
  double Oel[3];
  double c_tmp;
  double O2;
  double H46_idx_3;
  double O3;
  double O1dot;
  double ODotel[3];
  double a_temp[3];
  double O2dot;
  double O3dot;
  double b_dv[9];
  double S26_idx_0;
  int b_i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double integrationFactor;
  double b_c_tmp;
  double c_c_tmp;
  double S13_idx_0;
  double S35_idx_0;
  double S36_idx_0;
  double ycm;
  double S34_idx_0;
  double S45_idx_0;
  double S46_idx_0;
  double C12_idx_0;
  double K15_idx_1;
  double C13_idx_0;
  double K15_idx_0;
  double rhoA;
  double C23_idx_0;
  double K16_idx_0;
  double zcm;
  double C24_idx_0;
  double K56_idx_0;
  double C25_idx_0;
  double f_tmp;
  double C26_idx_0;
  double disMomentgp[3];
  double b_f_tmp;
  double C34_idx_0;
  double c_f_tmp;
  double C35_idx_0;
  double d_f_tmp;
  double C36_idx_0;
  double e_f_tmp;
  double posLocal[3];
  double f;
  double disLoadgpLocal[3];
  double C14_idx_0;
  double b_f;
  double c_f;
  double C45_idx_0;
  double C46_idx_0;
  double H13_idx_0;
  double H23_idx_0;
  double H24_idx_0;
  double H25_idx_0;
  double H26_idx_0;
  double H34_idx_0;
  double H35_idx_0;
  double H36_idx_0;
  double H14_idx_0;
  double H45_idx_0;
  double H46_idx_0;
  double S12_idx_1;
  double S13_idx_1;
  double S23_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double S35_idx_1;
  double S36_idx_1;
  double S14_idx_1;
  double S24_idx_1;
  double S34_idx_1;
  double S45_idx_1;
  double S46_idx_1;
  double C12_idx_1;
  double C13_idx_1;
  double C23_idx_1;
  double C24_idx_1;
  double C25_idx_1;
  double C26_idx_1;
  double C34_idx_1;
  double C35_idx_1;
  double C36_idx_1;
  double C14_idx_1;
  double C45_idx_1;
  double C46_idx_1;
  double H12_idx_1;
  double H13_idx_1;
  double H23_idx_1;
  double H24_idx_1;
  double H25_idx_1;
  double H26_idx_1;
  double H34_idx_1;
  double H35_idx_1;
  double H36_idx_1;
  double H14_idx_1;
  double H45_idx_1;
  double H46_idx_1;
  double S12_idx_2;
  double S13_idx_2;
  double S23_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S14_idx_2;
  double S24_idx_2;
  double S34_idx_2;
  double S45_idx_2;
  double S46_idx_2;
  double C12_idx_2;
  double C13_idx_2;
  double C23_idx_2;
  double C24_idx_2;
  double C25_idx_2;
  double C26_idx_2;
  double C34_idx_2;
  double C35_idx_2;
  double C36_idx_2;
  double C14_idx_2;
  double C45_idx_2;
  double C46_idx_2;
  double H12_idx_2;
  double H13_idx_2;
  double H23_idx_2;
  double H24_idx_2;
  double H25_idx_2;
  double H26_idx_2;
  double H34_idx_2;
  double H35_idx_2;
  double H36_idx_2;
  double H14_idx_2;
  double H45_idx_2;
  double H46_idx_2;
  double S12_idx_3;
  double S13_idx_3;
  double S23_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S14_idx_3;
  double S24_idx_3;
  double S34_idx_3;
  double S45_idx_3;
  double S46_idx_3;
  double C12_idx_3;
  double C13_idx_3;
  double C23_idx_3;
  double C24_idx_3;
  double C25_idx_3;
  double C26_idx_3;
  double C34_idx_3;
  double C35_idx_3;
  double C36_idx_3;
  double C14_idx_3;
  double C45_idx_3;
  double C46_idx_3;
  double H12_idx_3;
  double H13_idx_3;
  double H23_idx_3;
  double H24_idx_3;
  double H25_idx_3;
  double H26_idx_3;
  double H34_idx_3;
  double H35_idx_3;
  double H36_idx_3;
  double H14_idx_3;
  double H45_idx_3;
  double b_elStorage[144];
  double Kenr[144];
  double K16_idx_1;
  double K56_idx_1;
  double K15_idx_2;
  double K16_idx_2;
  double K56_idx_2;
  double K15_idx_3;
  double K16_idx_3;
  double K56_idx_3;
  double Khate[144];
  double reshapes_f3_data[24];
  double reshapes_f4_data[24];
  double reshapes_f1[24];
  double reshapes_f2[24];
  double reshapes_f5[24];
  double reshapes_f6[24];
  double Ce[144];
  int reshapes_f1_data_tmp;
  double Me[144];
  emxArray_real_T *lambda_d;
  emxArray_int32_T *lambda_colidx;
  emxArray_int32_T *lambda_rowidx;
  emxArray_real_T *lambdaTran_d;
  emxArray_int32_T *lambdaTran_colidx;
  emxArray_int32_T *lambdaTran_rowidx;
  double Ftemp_data[12];
  double Fhate[12];
  double Fe[12];
  int tmp_size[1];
  boolean_T concMassFlag;
  double FhatLessConc[12];
  double b_data[12];

  // -------- assign input block ----------------
  //  modalFlag      = input.modalFlag;
  // initialize CN2H to identity for static or modal analysis
  // declare type
  dispdot_size_idx_1 = 1;
  dispdot_data[0] = 0.0;

  // declare type
  dispddot_size_idx_1 = 1;
  dispddot_data[0] = 0.0;

  // declare type
  // options for Dean integrator
  if (c_strcmp(input->analysisType)) {
    // options for newmark beta integrator
    dispdot_size_idx_1 = 12;
    dispddot_size_idx_1 = 12;
    std::memcpy(&dispdot_data[0], &input->dispdot.data[0], 12U * sizeof(double));
    std::memcpy(&dispddot_data[0], &input->dispddot.data[0], 12U * sizeof(double));
  }

  // --------------------------------------------
  // setting for modal analysis flag
  // setting for initial reduced order model calculations
  // settings if aeroelastic analysis is active
  // Not used, but must be declared
  // number of gauss points for full integration
  // number of gauss points for reduced integration
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  F1_data_idx_0 = 0.0;
  F3_data_idx_0 = 0.0;
  F2_data_idx_0 = 0.0;
  F4_data_idx_0 = 0.0;
  F5_data_idx_0 = 0.0;
  F6_data_idx_0 = 0.0;
  F1_data_idx_1 = 0.0;
  F3_data_idx_1 = 0.0;
  F2_data_idx_1 = 0.0;
  F4_data_idx_1 = 0.0;
  F5_data_idx_1 = 0.0;
  F6_data_idx_1 = 0.0;

  // initialize pre-stress (stress stiffening matrices)
  // initialize nonlinear element matrices, only used if (useDisp)
  // initialize aeroelastic matrices only used if aeroElasticOn, but must declare type 
  // Convert frequencies from Hz to radians
  Omega = 6.2831853071795862 * input->Omega;
  OmegaDot = 6.2831853071795862 * input->OmegaDot;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input->sweepAngle * 3.1415926535897931 / 180.0,
                  input->coneAngle * 3.1415926535897931 / 180.0,
                  (input->rollAngle + 0.5 * (input->sectionProps.twist[0] +
    input->sectionProps.twist[1])) * 3.1415926535897931 / 180.0, lambda);

  //      theta_xNode = [dispLocal(4)  dispLocal(10)];
  //      theta_yNode = [dispLocal(5)  dispLocal(11)];
  //      theta_zNode = [dispLocal(6)  dispLocal(12)];
  for (i = 0; i < 3; i++) {
    c_tmp = lambda[i + 24];
    H46_idx_3 = lambda[i] * 0.0 + lambda[i + 12] * 0.0;
    Oel[i] = H46_idx_3 + c_tmp * Omega;
    a_temp[i] = (input->CN2H[i] * 0.0 + input->CN2H[i + 3] * 0.0) + input->
      CN2H[i + 6] * 9.81;
    ODotel[i] = H46_idx_3 + c_tmp * OmegaDot;
  }

  O1 = Oel[0];
  O2 = Oel[1];
  O3 = Oel[2];
  O1dot = ODotel[0];
  O2dot = ODotel[1];
  O3dot = ODotel[2];

  // gravitational acceleration [m/s^2]
  // acceleration of body in hub frame (from platform rigid body motion)
  // accelerations in inertial frame
  // Integration loop
  b_dv[0] = 0.0;
  b_dv[4] = 0.0;
  b_dv[7] = -0.0;
  b_dv[5] = 0.0;
  b_dv[8] = 0.0;
  S26_idx_0 = O1 * O3;
  for (b_i = 0; b_i < 4; b_i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[b_i], input->xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &H46_idx_3);
    integrationFactor = H46_idx_3 * dv1[b_i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    //  Only used if (useDisp || preStress)
    // mass center offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Not used but must declare type
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoA = N_data[0] * input->sectionProps.rhoA[0] + N_data[1] *
      input->sectionProps.rhoA[1];
    ycm = N_data[0] * input->sectionProps.ycm[0] + N_data[1] *
      input->sectionProps.ycm[1];
    zcm = N_data[0] * input->sectionProps.zcm[0] + N_data[1] *
      input->sectionProps.zcm[1];
    K15_idx_0 = N_data[0] * input->x.data[0] + N_data[1] * input->x.data[1];
    K16_idx_0 = N_data[0] * input->y.data[0] + N_data[1] * input->y.data[1];
    H46_idx_3 = N_data[0] * input->z.data[0] + N_data[1] * input->z.data[1];

    // let these loads be defined in the inertial frame
    disMomentgp[0] = rhoA * a_temp[0];
    disMomentgp[1] = rhoA * a_temp[1];
    disMomentgp[2] = rhoA * a_temp[2];
    for (i = 0; i < 3; i++) {
      c_tmp = lambda[i + 12];
      K15_idx_1 = lambda[i + 24];
      posLocal[i] = (lambda[i] * K15_idx_0 + c_tmp * K16_idx_0) + K15_idx_1 *
        H46_idx_3;
      disLoadgpLocal[i] = (lambda[i] * disMomentgp[0] + c_tmp * disMomentgp[1])
        + K15_idx_1 * disMomentgp[2];
    }

    b_dv[3] = -zcm;
    b_dv[6] = ycm;
    b_dv[1] = zcm;
    b_dv[2] = -ycm;
    for (i = 0; i < 3; i++) {
      disMomentgp[i] = (b_dv[i] * disLoadgpLocal[0] + b_dv[i + 3] *
                        disLoadgpLocal[1]) + b_dv[i + 6] * disLoadgpLocal[2];
    }

    // calculate static aerodynamic load
    // distributed/body force load calculations
    f_tmp = (O2 * O2 + O3 * O3) * posLocal[0];
    b_f_tmp = O2dot * posLocal[2];
    c_f_tmp = O3dot * posLocal[1];
    d_f_tmp = S26_idx_0 * posLocal[2];
    e_f_tmp = O1 * O2 * posLocal[1];
    f = rhoA * ((((f_tmp - e_f_tmp) - d_f_tmp) + c_f_tmp) - b_f_tmp) -
      disLoadgpLocal[0];

    // This function is a general routine to calculate an element vector
    H46_idx_3 = O1dot * posLocal[2];
    K56_idx_0 = O3dot * posLocal[0];
    b_f = rhoA * (((((O1 * O1 + O3 * O3) * posLocal[1] - posLocal[2] * O2 * O3)
                    - posLocal[0] * O1 * O2) + H46_idx_3) - K56_idx_0) -
      disLoadgpLocal[1];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = O2dot * posLocal[0];
    K16_idx_0 = O1dot * posLocal[1];
    c_f = rhoA * (((((O1 * O1 + O2 * O2) * posLocal[2] - S26_idx_0 * posLocal[0])
                    - O2 * O3 * posLocal[1]) + K15_idx_0) - K16_idx_0) -
      disLoadgpLocal[2];

    // This function is a general routine to calculate an element vector
    K15_idx_0 = rhoA * ((((posLocal[0] * (O1 * O2 * zcm - ycm * O1 * O3) -
      posLocal[1] * (ycm * O2 * O3 + zcm * (O1 * O1 + O3 * O3))) + posLocal[2] *
                          (ycm * (O1 * O1 + O2 * O2) + zcm * O2 * O3)) + ycm *
                         (K15_idx_0 - K16_idx_0)) - zcm * (H46_idx_3 - K56_idx_0))
      - disMomentgp[0];

    // This function is a general routine to calculate an element vector
    H46_idx_3 = rhoA * zcm * ((((f_tmp - posLocal[1] * O1 * O2) - posLocal[2] *
      O1 * O3) - b_f_tmp) + c_f_tmp) - disMomentgp[1];

    // This function is a general routine to calculate an element vector
    K56_idx_0 = rhoA * ycm * ((((d_f_tmp + e_f_tmp) - f_tmp) - c_f_tmp) +
      b_f_tmp) - disMomentgp[2];

    // This function is a general routine to calculate an element vector
    F1_data_idx_0 += f * N_data[0] * integrationFactor;
    F2_data_idx_0 += b_f * N_data[0] * integrationFactor;
    F3_data_idx_0 += c_f * N_data[0] * integrationFactor;
    F4_data_idx_0 += K15_idx_0 * N_data[0] * integrationFactor;
    F5_data_idx_0 += H46_idx_3 * N_data[0] * integrationFactor;
    F6_data_idx_0 += K56_idx_0 * N_data[0] * integrationFactor;
    F1_data_idx_1 += f * N_data[1] * integrationFactor;
    F2_data_idx_1 += b_f * N_data[1] * integrationFactor;
    F3_data_idx_1 += c_f * N_data[1] * integrationFactor;
    F4_data_idx_1 += K15_idx_0 * N_data[1] * integrationFactor;
    F5_data_idx_1 += H46_idx_3 * N_data[1] * integrationFactor;
    F6_data_idx_1 += K56_idx_0 * N_data[1] * integrationFactor;
  }

  // END OF INTEGRATION LOOP
  // Integration loop
  // Calculate shape functions at quad point i
  // ..... interpolate for value at quad point .....
  // END OF REDUCED INTEGRATION LOOP
  // unpack stored element stiffness data
  //  Only used if (useDisp)
  // unpack stored element mass data
  // unpack and scale stored element spin softening data
  H46_idx_3 = Oel[0] * Oel[1];
  b_c_tmp = Oel[0] * Oel[0];
  c_c_tmp = b_c_tmp + Oel[2] * Oel[2];
  b_c_tmp += Oel[1] * Oel[1];

  // unpack and scale stored element Corilois data
  // unpack and scale stored element Circulatory data
  integrationFactor = elStorage->S12[0] * Oel[0] * Oel[1];
  S13_idx_0 = elStorage->S13[0] * Oel[0] * Oel[2];
  O1 = elStorage->S23[0] * Oel[1] * Oel[2];
  O2dot = elStorage->S25[0] * H46_idx_3;
  S26_idx_0 = elStorage->S26[0] * H46_idx_3;
  S35_idx_0 = elStorage->S35[0] * Oel[0] * Oel[2];
  S36_idx_0 = elStorage->S36[0] * Oel[0] * Oel[2];
  ycm = elStorage->S14_1[0] * Oel[0] * Oel[2] + elStorage->S14_2[0] * Oel[0] *
    Oel[1];
  O3 = elStorage->S24_1[0] * c_c_tmp + elStorage->S24_2[0] * Oel[1] * Oel[2];
  S34_idx_0 = elStorage->S34_1[0] * b_c_tmp + elStorage->S34_2[0] * Oel[1] *
    Oel[2];
  S45_idx_0 = elStorage->S45_1[0] * Oel[0] * Oel[2] + elStorage->S45_2[0] * Oel
    [0] * Oel[1];
  S46_idx_0 = elStorage->S46_1[0] * Oel[0] * Oel[1] + elStorage->S46_2[0] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[0];
  C12_idx_0 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[0];
  C13_idx_0 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[0];
  C23_idx_0 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[0];
  C24_idx_0 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[0];
  C25_idx_0 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[0];
  C26_idx_0 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[0];
  C34_idx_0 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[0];
  C35_idx_0 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[0];
  C36_idx_0 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[0];
  f = elStorage->C14_2[0];
  C14_idx_0 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[0];
  c_f = elStorage->C45_2[0];
  C45_idx_0 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[0];
  rhoA = elStorage->C46_2[0];
  C46_idx_0 = zcm * Oel[1] + rhoA * Oel[2];
  O2 = 0.5 * c_tmp * ODotel[2];
  H13_idx_0 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_0 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_0 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_0 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_0 = 0.5 * f_tmp * ODotel[2];
  H34_idx_0 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_0 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_0 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_0 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_0 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_0 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_1 = elStorage->S12[1] * Oel[0] * Oel[1];
  S13_idx_1 = elStorage->S13[1] * Oel[0] * Oel[2];
  S23_idx_1 = elStorage->S23[1] * Oel[1] * Oel[2];
  S25_idx_1 = elStorage->S25[1] * H46_idx_3;
  S26_idx_1 = elStorage->S26[1] * H46_idx_3;
  S35_idx_1 = elStorage->S35[1] * Oel[0] * Oel[2];
  S36_idx_1 = elStorage->S36[1] * Oel[0] * Oel[2];
  S14_idx_1 = elStorage->S14_1[1] * Oel[0] * Oel[2] + elStorage->S14_2[1] * Oel
    [0] * Oel[1];
  S24_idx_1 = elStorage->S24_1[1] * c_c_tmp + elStorage->S24_2[1] * Oel[1] *
    Oel[2];
  S34_idx_1 = elStorage->S34_1[1] * b_c_tmp + elStorage->S34_2[1] * Oel[1] *
    Oel[2];
  S45_idx_1 = elStorage->S45_1[1] * Oel[0] * Oel[2] + elStorage->S45_2[1] * Oel
    [0] * Oel[1];
  S46_idx_1 = elStorage->S46_1[1] * Oel[0] * Oel[1] + elStorage->S46_2[1] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[1];
  C12_idx_1 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[1];
  C13_idx_1 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[1];
  C23_idx_1 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[1];
  C24_idx_1 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[1];
  C25_idx_1 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[1];
  C26_idx_1 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[1];
  C34_idx_1 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[1];
  C35_idx_1 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[1];
  C36_idx_1 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[1];
  f = elStorage->C14_2[1];
  C14_idx_1 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[1];
  c_f = elStorage->C45_2[1];
  C45_idx_1 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[1];
  rhoA = elStorage->C46_2[1];
  C46_idx_1 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_1 = 0.5 * c_tmp * ODotel[2];
  H13_idx_1 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_1 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_1 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_1 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_1 = 0.5 * f_tmp * ODotel[2];
  H34_idx_1 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_1 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_1 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_1 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_1 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_1 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_2 = elStorage->S12[2] * Oel[0] * Oel[1];
  S13_idx_2 = elStorage->S13[2] * Oel[0] * Oel[2];
  S23_idx_2 = elStorage->S23[2] * Oel[1] * Oel[2];
  S25_idx_2 = elStorage->S25[2] * H46_idx_3;
  S26_idx_2 = elStorage->S26[2] * H46_idx_3;
  S35_idx_2 = elStorage->S35[2] * Oel[0] * Oel[2];
  S36_idx_2 = elStorage->S36[2] * Oel[0] * Oel[2];
  S14_idx_2 = elStorage->S14_1[2] * Oel[0] * Oel[2] + elStorage->S14_2[2] * Oel
    [0] * Oel[1];
  S24_idx_2 = elStorage->S24_1[2] * c_c_tmp + elStorage->S24_2[2] * Oel[1] *
    Oel[2];
  S34_idx_2 = elStorage->S34_1[2] * b_c_tmp + elStorage->S34_2[2] * Oel[1] *
    Oel[2];
  S45_idx_2 = elStorage->S45_1[2] * Oel[0] * Oel[2] + elStorage->S45_2[2] * Oel
    [0] * Oel[1];
  S46_idx_2 = elStorage->S46_1[2] * Oel[0] * Oel[1] + elStorage->S46_2[2] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[2];
  C12_idx_2 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[2];
  C13_idx_2 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[2];
  C23_idx_2 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[2];
  C24_idx_2 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[2];
  C25_idx_2 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[2];
  C26_idx_2 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[2];
  C34_idx_2 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[2];
  C35_idx_2 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[2];
  C36_idx_2 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[2];
  f = elStorage->C14_2[2];
  C14_idx_2 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[2];
  c_f = elStorage->C45_2[2];
  C45_idx_2 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[2];
  rhoA = elStorage->C46_2[2];
  C46_idx_2 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_2 = 0.5 * c_tmp * ODotel[2];
  H13_idx_2 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_2 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_2 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_2 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_2 = 0.5 * f_tmp * ODotel[2];
  H34_idx_2 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_2 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_2 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_2 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_2 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_2 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);
  S12_idx_3 = elStorage->S12[3] * Oel[0] * Oel[1];
  S13_idx_3 = elStorage->S13[3] * Oel[0] * Oel[2];
  S23_idx_3 = elStorage->S23[3] * Oel[1] * Oel[2];
  S25_idx_3 = elStorage->S25[3] * H46_idx_3;
  S26_idx_3 = elStorage->S26[3] * H46_idx_3;
  S35_idx_3 = elStorage->S35[3] * Oel[0] * Oel[2];
  S36_idx_3 = elStorage->S36[3] * Oel[0] * Oel[2];
  S14_idx_3 = elStorage->S14_1[3] * Oel[0] * Oel[2] + elStorage->S14_2[3] * Oel
    [0] * Oel[1];
  S24_idx_3 = elStorage->S24_1[3] * c_c_tmp + elStorage->S24_2[3] * Oel[1] *
    Oel[2];
  S34_idx_3 = elStorage->S34_1[3] * b_c_tmp + elStorage->S34_2[3] * Oel[1] *
    Oel[2];
  S45_idx_3 = elStorage->S45_1[3] * Oel[0] * Oel[2] + elStorage->S45_2[3] * Oel
    [0] * Oel[1];
  S46_idx_3 = elStorage->S46_1[3] * Oel[0] * Oel[1] + elStorage->S46_2[3] * Oel
    [0] * Oel[2];
  c_tmp = elStorage->C12[3];
  C12_idx_3 = c_tmp * Oel[2];
  K15_idx_1 = elStorage->C13[3];
  C13_idx_3 = K15_idx_1 * Oel[1];
  K15_idx_0 = elStorage->C23[3];
  C23_idx_3 = K15_idx_0 * Oel[0];
  K16_idx_0 = elStorage->C24[3];
  C24_idx_3 = K16_idx_0 * Oel[0];
  K56_idx_0 = elStorage->C25[3];
  C25_idx_3 = K56_idx_0 * Oel[2];
  f_tmp = elStorage->C26[3];
  C26_idx_3 = f_tmp * Oel[2];
  b_f_tmp = elStorage->C34[3];
  C34_idx_3 = b_f_tmp * Oel[0];
  c_f_tmp = elStorage->C35[3];
  C35_idx_3 = c_f_tmp * Oel[1];
  d_f_tmp = elStorage->C36[3];
  C36_idx_3 = d_f_tmp * Oel[1];
  e_f_tmp = elStorage->C14_1[3];
  f = elStorage->C14_2[3];
  C14_idx_3 = e_f_tmp * Oel[1] + f * Oel[2];
  b_f = elStorage->C45_1[3];
  c_f = elStorage->C45_2[3];
  C45_idx_3 = b_f * Oel[2] + c_f * Oel[1];
  zcm = elStorage->C46_1[3];
  rhoA = elStorage->C46_2[3];
  C46_idx_3 = zcm * Oel[1] + rhoA * Oel[2];
  H12_idx_3 = 0.5 * c_tmp * ODotel[2];
  H13_idx_3 = 0.5 * K15_idx_1 * ODotel[1];
  H23_idx_3 = 0.5 * K15_idx_0 * ODotel[0];
  H24_idx_3 = 0.5 * K16_idx_0 * ODotel[0];
  H25_idx_3 = 0.5 * K56_idx_0 * ODotel[2];
  H26_idx_3 = 0.5 * f_tmp * ODotel[2];
  H34_idx_3 = 0.5 * b_f_tmp * ODotel[0];
  H35_idx_3 = 0.5 * c_f_tmp * ODotel[1];
  H36_idx_3 = 0.5 * d_f_tmp * ODotel[1];
  H14_idx_3 = 0.5 * (e_f_tmp * ODotel[1] + f * ODotel[2]);
  H45_idx_3 = 0.5 * (b_f * ODotel[2] + c_f * ODotel[1]);
  H46_idx_3 = 0.5 * (zcm * ODotel[1] + rhoA * ODotel[2]);

  // compile stiffness matrix without rotational effects
  b_elStorage[0] = elStorage->K11[0];
  b_elStorage[24] = elStorage->K12[0];
  b_elStorage[48] = elStorage->K13[0];
  b_elStorage[72] = elStorage->K14[0];
  b_elStorage[96] = elStorage->K15[0];
  b_elStorage[120] = elStorage->K16[0];
  b_elStorage[2] = elStorage->K12[0];
  b_elStorage[26] = elStorage->K22[0];
  b_elStorage[50] = elStorage->K23[0];
  b_elStorage[74] = elStorage->K24[0];
  b_elStorage[98] = elStorage->K25[0];
  b_elStorage[122] = elStorage->K26[0];
  b_elStorage[4] = elStorage->K13[0];
  b_elStorage[28] = elStorage->K23[0];
  b_elStorage[52] = elStorage->K33[0];
  b_elStorage[76] = elStorage->K34[0];
  b_elStorage[100] = elStorage->K35[0];
  b_elStorage[124] = elStorage->K36[0];
  b_elStorage[6] = elStorage->K13[0];
  b_elStorage[30] = elStorage->K24[0];
  b_elStorage[54] = elStorage->K34[0];
  b_elStorage[78] = elStorage->K44[0];
  b_elStorage[102] = elStorage->K45[0];
  b_elStorage[126] = elStorage->K46[0];
  b_elStorage[8] = elStorage->K15[0];
  b_elStorage[32] = elStorage->K25[0];
  b_elStorage[56] = elStorage->K35[0];
  b_elStorage[80] = elStorage->K45[0];
  b_elStorage[104] = elStorage->K55[0];
  b_elStorage[128] = elStorage->K56[0];
  b_elStorage[10] = elStorage->K16[0];
  b_elStorage[34] = elStorage->K26[0];
  b_elStorage[58] = elStorage->K36[0];
  b_elStorage[82] = elStorage->K46[0];
  b_elStorage[106] = elStorage->K56[0];
  b_elStorage[130] = elStorage->K66[0];
  b_elStorage[1] = elStorage->K11[1];
  b_elStorage[25] = elStorage->K12[1];
  b_elStorage[49] = elStorage->K13[1];
  b_elStorage[73] = elStorage->K14[1];
  b_elStorage[97] = elStorage->K15[1];
  b_elStorage[121] = elStorage->K16[1];
  b_elStorage[3] = elStorage->K12[2];
  b_elStorage[27] = elStorage->K22[1];
  b_elStorage[51] = elStorage->K23[1];
  b_elStorage[75] = elStorage->K24[1];
  b_elStorage[99] = elStorage->K25[1];
  b_elStorage[123] = elStorage->K26[1];
  b_elStorage[5] = elStorage->K13[2];
  b_elStorage[29] = elStorage->K23[2];
  b_elStorage[53] = elStorage->K33[1];
  b_elStorage[77] = elStorage->K34[1];
  b_elStorage[101] = elStorage->K35[1];
  b_elStorage[125] = elStorage->K36[1];
  b_elStorage[7] = elStorage->K13[2];
  b_elStorage[31] = elStorage->K24[2];
  b_elStorage[55] = elStorage->K34[2];
  b_elStorage[79] = elStorage->K44[1];
  b_elStorage[103] = elStorage->K45[1];
  b_elStorage[127] = elStorage->K46[1];
  b_elStorage[9] = elStorage->K15[2];
  b_elStorage[33] = elStorage->K25[2];
  b_elStorage[57] = elStorage->K35[2];
  b_elStorage[81] = elStorage->K45[2];
  b_elStorage[105] = elStorage->K55[1];
  b_elStorage[129] = elStorage->K56[1];
  b_elStorage[11] = elStorage->K16[2];
  b_elStorage[35] = elStorage->K26[2];
  b_elStorage[59] = elStorage->K36[2];
  b_elStorage[83] = elStorage->K46[2];
  b_elStorage[107] = elStorage->K56[2];
  b_elStorage[131] = elStorage->K66[1];
  b_elStorage[12] = elStorage->K11[2];
  b_elStorage[36] = elStorage->K12[2];
  b_elStorage[60] = elStorage->K13[2];
  b_elStorage[84] = elStorage->K14[2];
  b_elStorage[108] = elStorage->K15[2];
  b_elStorage[132] = elStorage->K16[2];
  b_elStorage[14] = elStorage->K12[1];
  b_elStorage[38] = elStorage->K22[2];
  b_elStorage[62] = elStorage->K23[2];
  b_elStorage[86] = elStorage->K24[2];
  b_elStorage[110] = elStorage->K25[2];
  b_elStorage[134] = elStorage->K26[2];
  b_elStorage[16] = elStorage->K13[1];
  b_elStorage[40] = elStorage->K23[1];
  b_elStorage[64] = elStorage->K33[2];
  b_elStorage[88] = elStorage->K34[2];
  b_elStorage[112] = elStorage->K35[2];
  b_elStorage[136] = elStorage->K36[2];
  b_elStorage[18] = elStorage->K13[1];
  b_elStorage[42] = elStorage->K24[1];
  b_elStorage[66] = elStorage->K34[1];
  b_elStorage[90] = elStorage->K44[2];
  b_elStorage[114] = elStorage->K45[2];
  b_elStorage[138] = elStorage->K46[2];
  b_elStorage[20] = elStorage->K15[1];
  b_elStorage[44] = elStorage->K25[1];
  b_elStorage[68] = elStorage->K35[1];
  b_elStorage[92] = elStorage->K45[1];
  b_elStorage[116] = elStorage->K55[2];
  b_elStorage[140] = elStorage->K56[2];
  b_elStorage[22] = elStorage->K16[1];
  b_elStorage[46] = elStorage->K26[1];
  b_elStorage[70] = elStorage->K36[1];
  b_elStorage[94] = elStorage->K46[1];
  b_elStorage[118] = elStorage->K56[1];
  b_elStorage[142] = elStorage->K66[2];
  b_elStorage[13] = elStorage->K11[3];
  b_elStorage[37] = elStorage->K12[3];
  b_elStorage[61] = elStorage->K13[3];
  b_elStorage[85] = elStorage->K14[3];
  b_elStorage[109] = elStorage->K15[3];
  b_elStorage[133] = elStorage->K16[3];
  b_elStorage[15] = elStorage->K12[3];
  b_elStorage[39] = elStorage->K22[3];
  b_elStorage[63] = elStorage->K23[3];
  b_elStorage[87] = elStorage->K24[3];
  b_elStorage[111] = elStorage->K25[3];
  b_elStorage[135] = elStorage->K26[3];
  b_elStorage[17] = elStorage->K13[3];
  b_elStorage[41] = elStorage->K23[3];
  b_elStorage[65] = elStorage->K33[3];
  b_elStorage[89] = elStorage->K34[3];
  b_elStorage[113] = elStorage->K35[3];
  b_elStorage[137] = elStorage->K36[3];
  b_elStorage[19] = elStorage->K13[3];
  b_elStorage[43] = elStorage->K24[3];
  b_elStorage[67] = elStorage->K34[3];
  b_elStorage[91] = elStorage->K44[3];
  b_elStorage[115] = elStorage->K45[3];
  b_elStorage[139] = elStorage->K46[3];
  b_elStorage[21] = elStorage->K15[3];
  b_elStorage[45] = elStorage->K25[3];
  b_elStorage[69] = elStorage->K35[3];
  b_elStorage[93] = elStorage->K45[3];
  b_elStorage[117] = elStorage->K55[3];
  b_elStorage[141] = elStorage->K56[3];
  b_elStorage[23] = elStorage->K16[3];
  b_elStorage[47] = elStorage->K26[3];
  b_elStorage[71] = elStorage->K36[3];
  b_elStorage[95] = elStorage->K46[3];
  b_elStorage[119] = elStorage->K56[3];
  b_elStorage[143] = elStorage->K66[3];
  mapMatrixNonSym(b_elStorage, Kenr);

  // add spin softening and circulatory effects to stiffness marix
  c_tmp = Oel[1] * Oel[1] + Oel[2] * Oel[2];
  K15_idx_0 = elStorage->K15[0] + elStorage->S15[0] * c_tmp;
  K16_idx_0 = elStorage->K16[0] + elStorage->S16[0] * c_tmp;
  K56_idx_0 = elStorage->K56[0] + elStorage->S56[0] * c_tmp;
  K15_idx_1 = elStorage->K15[1] + elStorage->S15[1] * c_tmp;
  K16_idx_1 = elStorage->K16[1] + elStorage->S16[1] * c_tmp;
  K56_idx_1 = elStorage->K56[1] + elStorage->S56[1] * c_tmp;
  K15_idx_2 = elStorage->K15[2] + elStorage->S15[2] * c_tmp;
  K16_idx_2 = elStorage->K16[2] + elStorage->S16[2] * c_tmp;
  K56_idx_2 = elStorage->K56[2] + elStorage->S56[2] * c_tmp;
  K15_idx_3 = elStorage->K15[3] + elStorage->S15[3] * c_tmp;
  K16_idx_3 = elStorage->K16[3] + elStorage->S16[3] * c_tmp;
  K56_idx_3 = elStorage->K56[3] + elStorage->S56[3] * c_tmp;

  // ---------------------------------------------
  // compile stiffness matrix with rotational effects
  b_elStorage[0] = elStorage->K11[0] + elStorage->S11[0] * c_tmp;
  b_elStorage[24] = (elStorage->K12[0] + integrationFactor) + O2;
  b_elStorage[48] = (elStorage->K13[0] + S13_idx_0) + H13_idx_0;
  O3dot = elStorage->K14[0] + ycm;
  b_elStorage[72] = O3dot + H14_idx_0;
  b_elStorage[96] = K15_idx_0;
  b_elStorage[120] = K16_idx_0;
  b_elStorage[2] = (elStorage->K12[0] + integrationFactor) - O2;
  b_elStorage[26] = elStorage->K22[0] + elStorage->S22[0] * c_c_tmp;
  O1dot = elStorage->K23[0] + O1;
  b_elStorage[50] = O1dot + H23_idx_0;
  O3 += elStorage->K24[0];
  b_elStorage[74] = O3 + H24_idx_0;
  O2 = elStorage->K25[0] + O2dot;
  b_elStorage[98] = O2 + H25_idx_0;
  O1 = elStorage->K26[0] + S26_idx_0;
  b_elStorage[122] = O1 + H26_idx_0;
  b_elStorage[4] = (elStorage->K13[0] + S13_idx_0) - H13_idx_0;
  b_elStorage[28] = O1dot - H23_idx_0;
  b_elStorage[52] = elStorage->K33[0] + elStorage->S33[0] * b_c_tmp;
  O1dot = elStorage->K34[0] + S34_idx_0;
  b_elStorage[76] = O1dot + H34_idx_0;
  integrationFactor = elStorage->K35[0] + S35_idx_0;
  b_elStorage[100] = integrationFactor + H35_idx_0;
  ycm = elStorage->K36[0] + S36_idx_0;
  b_elStorage[124] = ycm + H36_idx_0;
  b_elStorage[6] = O3dot - H14_idx_0;
  b_elStorage[30] = O3 - H24_idx_0;
  b_elStorage[54] = O1dot - H34_idx_0;
  b_elStorage[78] = elStorage->K44[0] + ((elStorage->S44_1[0] * c_c_tmp +
    elStorage->S44_2[0] * b_c_tmp) + elStorage->S44_3[0] * Oel[1] * Oel[2]);
  O3dot = elStorage->K45[0] + S45_idx_0;
  b_elStorage[102] = O3dot + H45_idx_0;
  O1dot = elStorage->K46[0] + S46_idx_0;
  b_elStorage[126] = O1dot + H46_idx_0;
  b_elStorage[8] = K15_idx_0;
  b_elStorage[32] = O2 - H25_idx_0;
  b_elStorage[56] = integrationFactor - H35_idx_0;
  b_elStorage[80] = O3dot - H45_idx_0;
  b_elStorage[104] = elStorage->K55[0] + elStorage->S55[0] * c_tmp;
  b_elStorage[128] = K56_idx_0;
  b_elStorage[10] = K16_idx_0;
  b_elStorage[34] = O1 - H26_idx_0;
  b_elStorage[58] = ycm - H36_idx_0;
  b_elStorage[82] = O1dot - H46_idx_0;
  b_elStorage[106] = K56_idx_0;
  b_elStorage[130] = elStorage->K66[0] + elStorage->S66[0] * c_tmp;
  b_elStorage[1] = elStorage->K11[1] + elStorage->S11[1] * c_tmp;
  b_elStorage[25] = (elStorage->K12[1] + S12_idx_1) + H12_idx_1;
  b_elStorage[49] = (elStorage->K13[1] + S13_idx_1) + H13_idx_1;
  O3dot = elStorage->K14[1] + S14_idx_1;
  b_elStorage[73] = O3dot + H14_idx_1;
  b_elStorage[97] = K15_idx_1;
  b_elStorage[121] = K16_idx_1;
  b_elStorage[3] = (elStorage->K12[2] + S12_idx_2) - H12_idx_2;
  b_elStorage[27] = elStorage->K22[1] + elStorage->S22[1] * c_c_tmp;
  O1dot = elStorage->K23[1] + S23_idx_1;
  b_elStorage[51] = O1dot + H23_idx_1;
  O3 = elStorage->K24[1] + S24_idx_1;
  b_elStorage[75] = O3 + H24_idx_1;
  O2 = elStorage->K25[1] + S25_idx_1;
  b_elStorage[99] = O2 + H25_idx_1;
  O1 = elStorage->K26[1] + S26_idx_1;
  b_elStorage[123] = O1 + H26_idx_1;
  b_elStorage[5] = (elStorage->K13[2] + S13_idx_2) - H13_idx_2;
  integrationFactor = elStorage->K23[2] + S23_idx_2;
  b_elStorage[29] = integrationFactor - H23_idx_2;
  b_elStorage[53] = elStorage->K33[1] + elStorage->S33[1] * b_c_tmp;
  ycm = elStorage->K34[1] + S34_idx_1;
  b_elStorage[77] = ycm + H34_idx_1;
  rhoA = elStorage->K35[1] + S35_idx_1;
  b_elStorage[101] = rhoA + H35_idx_1;
  zcm = elStorage->K36[1] + S36_idx_1;
  b_elStorage[125] = zcm + H36_idx_1;
  c_f = elStorage->K14[2] + S14_idx_2;
  b_elStorage[7] = c_f - H14_idx_2;
  b_f = elStorage->K24[2] + S24_idx_2;
  b_elStorage[31] = b_f - H24_idx_2;
  f = elStorage->K34[2] + S34_idx_2;
  b_elStorage[55] = f - H34_idx_2;
  b_elStorage[79] = elStorage->K44[1] + ((elStorage->S44_1[1] * c_c_tmp +
    elStorage->S44_2[1] * b_c_tmp) + elStorage->S44_3[1] * Oel[1] * Oel[2]);
  e_f_tmp = elStorage->K45[1] + S45_idx_1;
  b_elStorage[103] = e_f_tmp + H45_idx_1;
  d_f_tmp = elStorage->K46[1] + S46_idx_1;
  b_elStorage[127] = d_f_tmp + H46_idx_1;
  b_elStorage[9] = K15_idx_2;
  c_f_tmp = elStorage->K25[2] + S25_idx_2;
  b_elStorage[33] = c_f_tmp - H25_idx_2;
  b_f_tmp = elStorage->K35[2] + S35_idx_2;
  b_elStorage[57] = b_f_tmp - H35_idx_2;
  f_tmp = elStorage->K45[2] + S45_idx_2;
  b_elStorage[81] = f_tmp - H45_idx_2;
  b_elStorage[105] = elStorage->K55[1] + elStorage->S55[1] * c_tmp;
  b_elStorage[129] = K56_idx_1;
  b_elStorage[11] = K16_idx_2;
  K56_idx_0 = elStorage->K26[2] + S26_idx_2;
  b_elStorage[35] = K56_idx_0 - H26_idx_2;
  K16_idx_0 = elStorage->K36[2] + S36_idx_2;
  b_elStorage[59] = K16_idx_0 - H36_idx_2;
  K15_idx_0 = elStorage->K46[2] + S46_idx_2;
  b_elStorage[83] = K15_idx_0 - H46_idx_2;
  b_elStorage[107] = K56_idx_2;
  b_elStorage[131] = elStorage->K66[1] + elStorage->S66[1] * c_tmp;
  b_elStorage[12] = elStorage->K11[2] + elStorage->S11[2] * c_tmp;
  b_elStorage[36] = (elStorage->K12[2] + S12_idx_2) + H12_idx_2;
  b_elStorage[60] = (elStorage->K13[2] + S13_idx_2) + H13_idx_2;
  b_elStorage[84] = c_f + H14_idx_2;
  b_elStorage[108] = K15_idx_2;
  b_elStorage[132] = K16_idx_2;
  b_elStorage[14] = (elStorage->K12[1] + S12_idx_1) - H12_idx_1;
  b_elStorage[38] = elStorage->K22[2] + elStorage->S22[2] * c_c_tmp;
  b_elStorage[62] = integrationFactor + H23_idx_2;
  b_elStorage[86] = b_f + H24_idx_2;
  b_elStorage[110] = c_f_tmp + H25_idx_2;
  b_elStorage[134] = K56_idx_0 + H26_idx_2;
  b_elStorage[16] = (elStorage->K13[1] + S13_idx_1) - H13_idx_1;
  b_elStorage[40] = O1dot - H23_idx_1;
  b_elStorage[64] = elStorage->K33[2] + elStorage->S33[2] * b_c_tmp;
  b_elStorage[88] = f + H34_idx_2;
  b_elStorage[112] = b_f_tmp + H35_idx_2;
  b_elStorage[136] = K16_idx_0 + H36_idx_2;
  b_elStorage[18] = O3dot - H14_idx_1;
  b_elStorage[42] = O3 - H24_idx_1;
  b_elStorage[66] = ycm - H34_idx_1;
  b_elStorage[90] = elStorage->K44[2] + ((elStorage->S44_1[2] * c_c_tmp +
    elStorage->S44_2[2] * b_c_tmp) + elStorage->S44_3[2] * Oel[1] * Oel[2]);
  b_elStorage[114] = f_tmp + H45_idx_2;
  b_elStorage[138] = K15_idx_0 + H46_idx_2;
  b_elStorage[20] = K15_idx_1;
  b_elStorage[44] = O2 - H25_idx_1;
  b_elStorage[68] = rhoA - H35_idx_1;
  b_elStorage[92] = e_f_tmp - H45_idx_1;
  b_elStorage[116] = elStorage->K55[2] + elStorage->S55[2] * c_tmp;
  b_elStorage[140] = K56_idx_2;
  b_elStorage[22] = K16_idx_1;
  b_elStorage[46] = O1 - H26_idx_1;
  b_elStorage[70] = zcm - H36_idx_1;
  b_elStorage[94] = d_f_tmp - H46_idx_1;
  b_elStorage[118] = K56_idx_1;
  b_elStorage[142] = elStorage->K66[2] + elStorage->S66[2] * c_tmp;
  b_elStorage[13] = elStorage->K11[3] + elStorage->S11[3] * c_tmp;
  b_elStorage[37] = (elStorage->K12[3] + S12_idx_3) + H12_idx_3;
  b_elStorage[61] = (elStorage->K13[3] + S13_idx_3) + H13_idx_3;
  O3dot = elStorage->K14[3] + S14_idx_3;
  b_elStorage[85] = O3dot + H14_idx_3;
  b_elStorage[109] = K15_idx_3;
  b_elStorage[133] = K16_idx_3;
  b_elStorage[15] = (elStorage->K12[3] + S12_idx_3) - H12_idx_3;
  b_elStorage[39] = elStorage->K22[3] + elStorage->S22[3] * c_c_tmp;
  O1dot = elStorage->K23[3] + S23_idx_3;
  b_elStorage[63] = O1dot + H23_idx_3;
  O3 = elStorage->K24[3] + S24_idx_3;
  b_elStorage[87] = O3 + H24_idx_3;
  O2 = elStorage->K25[3] + S25_idx_3;
  b_elStorage[111] = O2 + H25_idx_3;
  O1 = elStorage->K26[3] + S26_idx_3;
  b_elStorage[135] = O1 + H26_idx_3;
  b_elStorage[17] = (elStorage->K13[3] + S13_idx_3) - H13_idx_3;
  b_elStorage[41] = O1dot - H23_idx_3;
  b_elStorage[65] = elStorage->K33[3] + elStorage->S33[3] * b_c_tmp;
  O1dot = elStorage->K34[3] + S34_idx_3;
  b_elStorage[89] = O1dot + H34_idx_3;
  integrationFactor = elStorage->K35[3] + S35_idx_3;
  b_elStorage[113] = integrationFactor + H35_idx_3;
  ycm = elStorage->K36[3] + S36_idx_3;
  b_elStorage[137] = ycm + H36_idx_3;
  b_elStorage[19] = O3dot - H14_idx_3;
  b_elStorage[43] = O3 - H24_idx_3;
  b_elStorage[67] = O1dot - H34_idx_3;
  b_elStorage[91] = elStorage->K44[3] + ((elStorage->S44_1[3] * c_c_tmp +
    elStorage->S44_2[3] * b_c_tmp) + elStorage->S44_3[3] * Oel[1] * Oel[2]);
  O3dot = elStorage->K45[3] + S45_idx_3;
  b_elStorage[115] = O3dot + H45_idx_3;
  O1dot = elStorage->K46[3] + S46_idx_3;
  b_elStorage[139] = O1dot + H46_idx_3;
  b_elStorage[21] = K15_idx_3;
  b_elStorage[45] = O2 - H25_idx_3;
  b_elStorage[69] = integrationFactor - H35_idx_3;
  b_elStorage[93] = O3dot - H45_idx_3;
  b_elStorage[117] = elStorage->K55[3] + elStorage->S55[3] * c_tmp;
  b_elStorage[141] = K56_idx_3;
  b_elStorage[23] = K16_idx_3;
  b_elStorage[47] = O1 - H26_idx_3;
  b_elStorage[71] = ycm - H36_idx_3;
  b_elStorage[95] = O1dot - H46_idx_3;
  b_elStorage[119] = K56_idx_3;
  b_elStorage[143] = elStorage->K66[3] + elStorage->S66[3] * c_tmp;
  mapMatrixNonSym(b_elStorage, Khate);

  //  Declare type
  // compile Coriolis/damping matrix
  reshapes_f3_data[0] = -C13_idx_0;
  reshapes_f3_data[4] = -C23_idx_0;
  reshapes_f3_data[8] = 0.0;
  reshapes_f3_data[12] = C34_idx_0;
  reshapes_f3_data[16] = C35_idx_0;
  reshapes_f3_data[20] = C36_idx_0;
  reshapes_f4_data[0] = -C14_idx_0;
  reshapes_f4_data[4] = -C24_idx_0;
  reshapes_f4_data[8] = -C34_idx_0;
  reshapes_f4_data[12] = 0.0;
  reshapes_f4_data[16] = C45_idx_0;
  reshapes_f4_data[20] = C46_idx_0;
  reshapes_f1[0] = 0.0;
  reshapes_f1[4] = C12_idx_0;
  reshapes_f1[8] = C13_idx_0;
  reshapes_f1[12] = C14_idx_0;
  reshapes_f1[16] = 0.0;
  reshapes_f1[20] = 0.0;
  reshapes_f2[0] = -C12_idx_0;
  reshapes_f2[4] = 0.0;
  reshapes_f2[8] = C23_idx_0;
  reshapes_f2[12] = C24_idx_0;
  reshapes_f2[16] = C25_idx_0;
  reshapes_f2[20] = C26_idx_0;
  reshapes_f5[0] = 0.0;
  reshapes_f5[4] = -C25_idx_0;
  reshapes_f5[8] = -C35_idx_0;
  reshapes_f5[12] = -C45_idx_0;
  reshapes_f5[16] = 0.0;
  reshapes_f5[20] = 0.0;
  reshapes_f6[0] = 0.0;
  reshapes_f6[4] = -C26_idx_0;
  reshapes_f6[8] = -C36_idx_0;
  reshapes_f6[12] = -C46_idx_0;
  reshapes_f6[16] = 0.0;
  reshapes_f6[20] = 0.0;
  reshapes_f3_data[1] = -C13_idx_2;
  reshapes_f3_data[5] = -C23_idx_2;
  reshapes_f3_data[9] = 0.0;
  reshapes_f3_data[13] = C34_idx_1;
  reshapes_f3_data[17] = C35_idx_1;
  reshapes_f3_data[21] = C36_idx_1;
  reshapes_f4_data[1] = -C14_idx_2;
  reshapes_f4_data[5] = -C24_idx_2;
  reshapes_f4_data[9] = -C34_idx_2;
  reshapes_f4_data[13] = 0.0;
  reshapes_f4_data[17] = C45_idx_1;
  reshapes_f4_data[21] = C46_idx_1;
  reshapes_f1[1] = 0.0;
  reshapes_f1[5] = C12_idx_1;
  reshapes_f1[9] = C13_idx_1;
  reshapes_f1[13] = C14_idx_1;
  reshapes_f1[17] = 0.0;
  reshapes_f1[21] = 0.0;
  reshapes_f2[1] = -C12_idx_2;
  reshapes_f2[5] = 0.0;
  reshapes_f2[9] = C23_idx_1;
  reshapes_f2[13] = C24_idx_1;
  reshapes_f2[17] = C25_idx_1;
  reshapes_f2[21] = C26_idx_1;
  reshapes_f5[1] = 0.0;
  reshapes_f5[5] = -C25_idx_2;
  reshapes_f5[9] = -C35_idx_2;
  reshapes_f5[13] = -C45_idx_2;
  reshapes_f5[17] = 0.0;
  reshapes_f5[21] = 0.0;
  reshapes_f6[1] = 0.0;
  reshapes_f6[5] = -C26_idx_2;
  reshapes_f6[9] = -C36_idx_2;
  reshapes_f6[13] = -C46_idx_2;
  reshapes_f6[17] = 0.0;
  reshapes_f6[21] = 0.0;
  reshapes_f3_data[2] = -C13_idx_1;
  reshapes_f3_data[6] = -C23_idx_1;
  reshapes_f3_data[10] = 0.0;
  reshapes_f3_data[14] = C34_idx_2;
  reshapes_f3_data[18] = C35_idx_2;
  reshapes_f3_data[22] = C36_idx_2;
  reshapes_f4_data[2] = -C14_idx_1;
  reshapes_f4_data[6] = -C24_idx_1;
  reshapes_f4_data[10] = -C34_idx_1;
  reshapes_f4_data[14] = 0.0;
  reshapes_f4_data[18] = C45_idx_2;
  reshapes_f4_data[22] = C46_idx_2;
  reshapes_f1[2] = 0.0;
  reshapes_f1[6] = C12_idx_2;
  reshapes_f1[10] = C13_idx_2;
  reshapes_f1[14] = C14_idx_2;
  reshapes_f1[18] = 0.0;
  reshapes_f1[22] = 0.0;
  reshapes_f2[2] = -C12_idx_1;
  reshapes_f2[6] = 0.0;
  reshapes_f2[10] = C23_idx_2;
  reshapes_f2[14] = C24_idx_2;
  reshapes_f2[18] = C25_idx_2;
  reshapes_f2[22] = C26_idx_2;
  reshapes_f5[2] = 0.0;
  reshapes_f5[6] = -C25_idx_1;
  reshapes_f5[10] = -C35_idx_1;
  reshapes_f5[14] = -C45_idx_1;
  reshapes_f5[18] = 0.0;
  reshapes_f5[22] = 0.0;
  reshapes_f6[2] = 0.0;
  reshapes_f6[6] = -C26_idx_1;
  reshapes_f6[10] = -C36_idx_1;
  reshapes_f6[14] = -C46_idx_1;
  reshapes_f6[18] = 0.0;
  reshapes_f6[22] = 0.0;
  reshapes_f3_data[3] = -C13_idx_3;
  reshapes_f3_data[7] = -C23_idx_3;
  reshapes_f3_data[11] = 0.0;
  reshapes_f3_data[15] = C34_idx_3;
  reshapes_f3_data[19] = C35_idx_3;
  reshapes_f3_data[23] = C36_idx_3;
  reshapes_f4_data[3] = -C14_idx_3;
  reshapes_f4_data[7] = -C24_idx_3;
  reshapes_f4_data[11] = -C34_idx_3;
  reshapes_f4_data[15] = 0.0;
  reshapes_f4_data[19] = C45_idx_3;
  reshapes_f4_data[23] = C46_idx_3;
  reshapes_f1[3] = 0.0;
  reshapes_f1[7] = C12_idx_3;
  reshapes_f1[11] = C13_idx_3;
  reshapes_f1[15] = C14_idx_3;
  reshapes_f1[19] = 0.0;
  reshapes_f1[23] = 0.0;
  reshapes_f2[3] = -C12_idx_3;
  reshapes_f2[7] = 0.0;
  reshapes_f2[11] = C23_idx_3;
  reshapes_f2[15] = C24_idx_3;
  reshapes_f2[19] = C25_idx_3;
  reshapes_f2[23] = C26_idx_3;
  reshapes_f5[3] = 0.0;
  reshapes_f5[7] = -C25_idx_3;
  reshapes_f5[11] = -C35_idx_3;
  reshapes_f5[15] = -C45_idx_3;
  reshapes_f5[19] = 0.0;
  reshapes_f5[23] = 0.0;
  reshapes_f6[3] = 0.0;
  reshapes_f6[7] = -C26_idx_3;
  reshapes_f6[11] = -C36_idx_3;
  reshapes_f6[15] = -C46_idx_3;
  reshapes_f6[19] = 0.0;
  reshapes_f6[23] = 0.0;
  for (i = 0; i < 12; i++) {
    b_i = i << 1;
    b_elStorage[12 * i] = reshapes_f1[b_i];
    b_elStorage[12 * i + 2] = reshapes_f2[b_i];
    b_elStorage[12 * i + 4] = reshapes_f3_data[2 * i];
    b_elStorage[12 * i + 6] = reshapes_f4_data[2 * i];
    b_elStorage[12 * i + 8] = reshapes_f5[b_i];
    b_elStorage[12 * i + 10] = reshapes_f6[b_i];
    b_i++;
    b_elStorage[12 * i + 1] = reshapes_f1[b_i];
    b_elStorage[12 * i + 3] = reshapes_f2[b_i];
    reshapes_f1_data_tmp = 2 * i + 1;
    b_elStorage[12 * i + 5] = reshapes_f3_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 7] = reshapes_f4_data[reshapes_f1_data_tmp];
    b_elStorage[12 * i + 9] = reshapes_f5[b_i];
    b_elStorage[12 * i + 11] = reshapes_f6[b_i];
  }

  b_mapMatrixNonSym(b_elStorage, Ce);

  // compile mass matrix
  b_elStorage[0] = elStorage->M11[0];
  b_elStorage[24] = 0.0;
  b_elStorage[48] = 0.0;
  b_elStorage[72] = 0.0;
  b_elStorage[96] = elStorage->M15[0];
  b_elStorage[120] = elStorage->M16[0];
  b_elStorage[2] = 0.0;
  b_elStorage[26] = elStorage->M22[0];
  b_elStorage[50] = 0.0;
  b_elStorage[74] = elStorage->M24[0];
  b_elStorage[98] = 0.0;
  b_elStorage[122] = 0.0;
  b_elStorage[4] = 0.0;
  b_elStorage[28] = 0.0;
  b_elStorage[52] = elStorage->M33[0];
  b_elStorage[76] = elStorage->M34[0];
  b_elStorage[100] = 0.0;
  b_elStorage[124] = 0.0;
  b_elStorage[6] = 0.0;
  b_elStorage[30] = elStorage->M24[0];
  b_elStorage[54] = elStorage->M34[0];
  b_elStorage[78] = elStorage->M44[0];
  b_elStorage[102] = 0.0;
  b_elStorage[126] = 0.0;
  b_elStorage[8] = elStorage->M15[0];
  b_elStorage[32] = 0.0;
  b_elStorage[56] = 0.0;
  b_elStorage[80] = 0.0;
  b_elStorage[104] = elStorage->M55[0];
  b_elStorage[128] = elStorage->M56[0];
  b_elStorage[10] = elStorage->M16[0];
  b_elStorage[34] = 0.0;
  b_elStorage[58] = 0.0;
  b_elStorage[82] = 0.0;
  b_elStorage[106] = elStorage->M56[0];
  b_elStorage[130] = elStorage->M66[0];
  b_elStorage[1] = elStorage->M11[1];
  b_elStorage[25] = 0.0;
  b_elStorage[49] = 0.0;
  b_elStorage[73] = 0.0;
  b_elStorage[97] = elStorage->M15[1];
  b_elStorage[121] = elStorage->M16[1];
  b_elStorage[3] = 0.0;
  b_elStorage[27] = elStorage->M22[1];
  b_elStorage[51] = 0.0;
  b_elStorage[75] = elStorage->M24[1];
  b_elStorage[99] = 0.0;
  b_elStorage[123] = 0.0;
  b_elStorage[5] = 0.0;
  b_elStorage[29] = 0.0;
  b_elStorage[53] = elStorage->M33[1];
  b_elStorage[77] = elStorage->M34[1];
  b_elStorage[101] = 0.0;
  b_elStorage[125] = 0.0;
  b_elStorage[7] = 0.0;
  b_elStorage[31] = elStorage->M24[2];
  b_elStorage[55] = elStorage->M34[2];
  b_elStorage[79] = elStorage->M44[1];
  b_elStorage[103] = 0.0;
  b_elStorage[127] = 0.0;
  b_elStorage[9] = elStorage->M15[2];
  b_elStorage[33] = 0.0;
  b_elStorage[57] = 0.0;
  b_elStorage[81] = 0.0;
  b_elStorage[105] = elStorage->M55[1];
  b_elStorage[129] = elStorage->M56[1];
  b_elStorage[11] = elStorage->M16[2];
  b_elStorage[35] = 0.0;
  b_elStorage[59] = 0.0;
  b_elStorage[83] = 0.0;
  b_elStorage[107] = elStorage->M56[2];
  b_elStorage[131] = elStorage->M66[1];
  b_elStorage[12] = elStorage->M11[2];
  b_elStorage[36] = 0.0;
  b_elStorage[60] = 0.0;
  b_elStorage[84] = 0.0;
  b_elStorage[108] = elStorage->M15[2];
  b_elStorage[132] = elStorage->M16[2];
  b_elStorage[14] = 0.0;
  b_elStorage[38] = elStorage->M22[2];
  b_elStorage[62] = 0.0;
  b_elStorage[86] = elStorage->M24[2];
  b_elStorage[110] = 0.0;
  b_elStorage[134] = 0.0;
  b_elStorage[16] = 0.0;
  b_elStorage[40] = 0.0;
  b_elStorage[64] = elStorage->M33[2];
  b_elStorage[88] = elStorage->M34[2];
  b_elStorage[112] = 0.0;
  b_elStorage[136] = 0.0;
  b_elStorage[18] = 0.0;
  b_elStorage[42] = elStorage->M24[1];
  b_elStorage[66] = elStorage->M34[1];
  b_elStorage[90] = elStorage->M44[2];
  b_elStorage[114] = 0.0;
  b_elStorage[138] = 0.0;
  b_elStorage[20] = elStorage->M15[1];
  b_elStorage[44] = 0.0;
  b_elStorage[68] = 0.0;
  b_elStorage[92] = 0.0;
  b_elStorage[116] = elStorage->M55[2];
  b_elStorage[140] = elStorage->M56[2];
  b_elStorage[22] = elStorage->M16[1];
  b_elStorage[46] = 0.0;
  b_elStorage[70] = 0.0;
  b_elStorage[94] = 0.0;
  b_elStorage[118] = elStorage->M56[1];
  b_elStorage[142] = elStorage->M66[2];
  b_elStorage[13] = elStorage->M11[3];
  b_elStorage[37] = 0.0;
  b_elStorage[61] = 0.0;
  b_elStorage[85] = 0.0;
  b_elStorage[109] = elStorage->M15[3];
  b_elStorage[133] = elStorage->M16[3];
  b_elStorage[15] = 0.0;
  b_elStorage[39] = elStorage->M22[3];
  b_elStorage[63] = 0.0;
  b_elStorage[87] = elStorage->M24[3];
  b_elStorage[111] = 0.0;
  b_elStorage[135] = 0.0;
  b_elStorage[17] = 0.0;
  b_elStorage[41] = 0.0;
  b_elStorage[65] = elStorage->M33[3];
  b_elStorage[89] = elStorage->M34[3];
  b_elStorage[113] = 0.0;
  b_elStorage[137] = 0.0;
  b_elStorage[19] = 0.0;
  b_elStorage[43] = elStorage->M24[3];
  b_elStorage[67] = elStorage->M34[3];
  b_elStorage[91] = elStorage->M44[3];
  b_elStorage[115] = 0.0;
  b_elStorage[139] = 0.0;
  b_elStorage[21] = elStorage->M15[3];
  b_elStorage[45] = 0.0;
  b_elStorage[69] = 0.0;
  b_elStorage[93] = 0.0;
  b_elStorage[117] = elStorage->M55[3];
  b_elStorage[141] = elStorage->M56[3];
  b_elStorage[23] = elStorage->M16[3];
  b_elStorage[47] = 0.0;
  b_elStorage[71] = 0.0;
  b_elStorage[95] = 0.0;
  b_elStorage[119] = elStorage->M56[3];
  b_elStorage[143] = elStorage->M66[3];
  mapMatrixNonSym(b_elStorage, Me);

  // account for rayleigh damping
  for (i = 0; i < 144; i++) {
    Ce[i] += input->RayleighAlpha * Kenr[i] + input->RayleighBeta * Me[i];
  }

  emxInit_real_T(&lambda_d, 1);
  emxInit_int32_T(&lambda_colidx, 1);
  emxInit_int32_T(&lambda_rowidx, 1);

  // compile element force vector
  //  transform matrices for sweep
  //  Note,a negative sweep angle, will sweep away from the direction of
  //  positive rotation
  sparse(lambda, lambda_d, lambda_colidx, lambda_rowidx);
  for (i = 0; i < 12; i++) {
    for (b_i = 0; b_i < 12; b_i++) {
      b_elStorage[b_i + 12 * i] = lambda[i + 12 * b_i];
    }
  }

  emxInit_real_T(&lambdaTran_d, 1);
  emxInit_int32_T(&lambdaTran_colidx, 1);
  emxInit_int32_T(&lambdaTran_rowidx, 1);
  sparse(b_elStorage, lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Me,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Me);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Ce,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Ce);
  b_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Khate,
                  b_elStorage);
  c_sparse_mtimes(b_elStorage, lambda_d, lambda_colidx, lambda_rowidx, Khate);
  Ftemp_data[0] = F1_data_idx_0;
  Ftemp_data[2] = F2_data_idx_0;
  Ftemp_data[4] = F3_data_idx_0;
  Ftemp_data[6] = F4_data_idx_0;
  Ftemp_data[8] = F5_data_idx_0;
  Ftemp_data[10] = F6_data_idx_0;
  Ftemp_data[1] = F1_data_idx_1;
  Ftemp_data[3] = F2_data_idx_1;
  Ftemp_data[5] = F3_data_idx_1;
  Ftemp_data[7] = F4_data_idx_1;
  Ftemp_data[9] = F5_data_idx_1;
  Ftemp_data[11] = F6_data_idx_1;

  // ----- function to form total force vector and transform to desired
  //  DOF mapping
  emxFree_int32_T(&lambda_rowidx);
  emxFree_int32_T(&lambda_colidx);
  emxFree_real_T(&lambda_d);
  std::memset(&Fhate[0], 0, 12U * sizeof(double));

  //
  //  %declare map
  for (b_i = 0; b_i < 12; b_i++) {
    Fhate[iv1[b_i] - 1] = Ftemp_data[b_i];
  }

  //  %------------------------------------------------------------------------- 
  d_sparse_mtimes(lambdaTran_d, lambdaTran_colidx, lambdaTran_rowidx, Fhate, Fe);

  //
  // concentrated mass
  // NOTE: Concentrated mass terms would modify 4,5,6 and 10,11,12 entries
  //  if some ycm or zcm offset from the node was accounted for in concentrated mass terms 
  d_eml_find(input->concMass, p_N_x_size, tmp_size);
  concMassFlag = (tmp_size[0] != 0);
  emxFree_int32_T(&lambdaTran_rowidx);
  emxFree_int32_T(&lambdaTran_colidx);
  emxFree_real_T(&lambdaTran_d);
  if (concMassFlag) {
    // modify Me for concentrated mass
    Me[0] += input->concMass[0];
    Me[13] += input->concMass[0];
    Me[26] += input->concMass[0];
    Me[39] += input->concMass[1];
    Me[52] += input->concMass[2];
    Me[65] += input->concMass[3];
    Me[78] += input->concMass[4];
    Me[91] += input->concMass[4];
    Me[104] += input->concMass[4];
    Me[117] += input->concMass[5];
    Me[130] += input->concMass[6];
    Me[143] += input->concMass[7];

    // modify Ce for concentrated mass
    H46_idx_3 = 2.0 * input->concMass[0] * Omega;
    Ce[12] -= H46_idx_3;
    Ce[1] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[0] * 0.0;
    Ce[24] += H46_idx_3;
    Ce[2] -= H46_idx_3;
    Ce[25] -= H46_idx_3;
    Ce[14] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[4] * Omega;
    Ce[90] -= H46_idx_3;
    Ce[79] += H46_idx_3;
    H46_idx_3 = 2.0 * input->concMass[4] * 0.0;
    Ce[102] += H46_idx_3;
    Ce[80] -= H46_idx_3;
    Ce[103] -= H46_idx_3;
    Ce[92] += H46_idx_3;

    // modify Ke for concentrated mass
    H46_idx_3 = Omega * Omega;
    K56_idx_0 = input->concMass[0] * H46_idx_3;
    Khate[0] -= K56_idx_0;
    K15_idx_0 = input->concMass[0] * 0.0 * 0.0;
    K16_idx_0 = input->concMass[0] * OmegaDot;
    Khate[12] = (Khate[12] + K15_idx_0) - K16_idx_0;
    Khate[1] = (Khate[1] + K15_idx_0) + K16_idx_0;
    K15_idx_0 = input->concMass[0] * 0.0 * Omega;
    Khate[24] = (Khate[24] + K15_idx_0) + input->concMass[0] * 0.0;
    Khate[2] = (Khate[2] + K15_idx_0) - input->concMass[0] * 0.0;
    Khate[25] = (Khate[25] + K15_idx_0) - input->concMass[0] * 0.0;
    Khate[14] = (Khate[14] + K15_idx_0) + input->concMass[0] * 0.0;
    Khate[13] -= K56_idx_0;
    Khate[26] -= input->concMass[0] * 0.0;
    K56_idx_0 = input->concMass[4] * H46_idx_3;
    Khate[78] -= K56_idx_0;
    K15_idx_0 = input->concMass[4] * 0.0 * 0.0;
    K16_idx_0 = input->concMass[4] * OmegaDot;
    Khate[90] = (Khate[90] + K15_idx_0) - K16_idx_0;
    Khate[79] = (Khate[79] + K15_idx_0) + K16_idx_0;
    K15_idx_0 = input->concMass[4] * 0.0 * Omega;
    Khate[102] = (Khate[102] + K15_idx_0) + input->concMass[4] * 0.0;
    Khate[80] = (Khate[80] + K15_idx_0) - input->concMass[4] * 0.0;
    Khate[103] = (Khate[103] + K15_idx_0) - input->concMass[4] * 0.0;
    Khate[92] = (Khate[92] + K15_idx_0) + input->concMass[4] * 0.0;
    Khate[91] -= K56_idx_0;
    Khate[104] -= input->concMass[4] * 0.0;
  }

  // modify Fe for  concentrated load
  if (concMassFlag) {
    H46_idx_3 = Omega * Omega;
    K56_idx_0 = 0.0 * Omega * input->z.data[0];
    Fe[0] = ((Fe[0] + input->concMass[0] * ((input->x.data[0] * H46_idx_3 - 0.0 *
                input->y.data[0]) - K56_idx_0)) + input->concMass[0] *
             (input->y.data[0] * OmegaDot - input->z.data[0] * 0.0)) -
      input->concMass[0] * a_temp[0];
    Fe[1] = ((Fe[1] + input->concMass[0] * ((input->y.data[0] * H46_idx_3 -
                K56_idx_0) - 0.0 * input->x.data[0])) + input->concMass[0] *
             (input->z.data[0] * 0.0 - input->x.data[0] * OmegaDot)) -
      input->concMass[0] * a_temp[1];
    Fe[2] = ((Fe[2] + input->concMass[0] * ((input->z.data[0] * 0.0 - Omega *
                0.0 * input->x.data[0]) - Omega * 0.0 * input->y.data[0])) +
             input->concMass[0] * (input->x.data[0] * 0.0 - input->y.data[0] *
              0.0)) - input->concMass[0] * a_temp[2];
    K56_idx_0 = 0.0 * Omega * input->z.data[1];
    Fe[6] = ((Fe[6] + input->concMass[4] * ((input->x.data[1] * H46_idx_3 - 0.0 *
                input->y.data[1]) - K56_idx_0)) + input->concMass[4] *
             (input->y.data[1] * OmegaDot - input->z.data[1] * 0.0)) -
      input->concMass[4] * a_temp[0];
    Fe[7] = ((Fe[7] + input->concMass[4] * ((input->y.data[1] * H46_idx_3 -
                K56_idx_0) - 0.0 * input->x.data[1])) + input->concMass[4] *
             (input->z.data[1] * 0.0 - input->x.data[1] * OmegaDot)) -
      input->concMass[4] * a_temp[1];
    Fe[8] = ((Fe[8] + input->concMass[4] * ((input->z.data[1] * 0.0 - Omega *
                0.0 * input->x.data[1]) - Omega * 0.0 * input->y.data[1])) +
             input->concMass[4] * (input->x.data[1] * 0.0 - input->y.data[1] *
              0.0)) - input->concMass[4] * a_temp[2];
  }

  //
  //  Declare Types
  std::memset(&Fhate[0], 0, 12U * sizeof(double));
  std::memset(&FhatLessConc[0], 0, 12U * sizeof(double));
  if (c_strcmp(input->analysisType)) {
    // calculate effective stiffness matrix and load vector for Newmark-Beta integrator 
    //      a1 = timeInt.a1;
    //      a2 = timeInt.a2;
    if (e_strcmp(input->iterationType)) {
      // considerations if newton raphson iteration is used
      if (input->firstIteration) {
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, FhatLessConc);
        std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
        mtimes(Khate, b_data, Ftemp_data);
        for (i = 0; i < 12; i++) {
          Fhate[i] = ((Fe[i] + Fhate[i]) + FhatLessConc[i]) - Ftemp_data[i];
        }
      } else {
        if (0 <= dispddot_size_idx_1 - 1) {
          std::memcpy(&b_data[0], &dispddot_data[0], dispddot_size_idx_1 *
                      sizeof(double));
        }

        if (dispddot_size_idx_1 == 1) {
          for (i = 0; i < 12; i++) {
            c_tmp = 0.0;
            for (b_i = 0; b_i < 12; b_i++) {
              c_tmp += Me[i + 12 * b_i] * b_data[b_i];
            }

            Fhate[i] = c_tmp;
          }
        } else {
          mtimes(Me, b_data, Fhate);
        }

        if (0 <= dispdot_size_idx_1 - 1) {
          std::memcpy(&b_data[0], &dispdot_data[0], dispdot_size_idx_1 * sizeof
                      (double));
        }

        if (dispdot_size_idx_1 == 1) {
          for (i = 0; i < 12; i++) {
            c_tmp = 0.0;
            for (b_i = 0; b_i < 12; b_i++) {
              c_tmp += Ce[i + 12 * b_i] * b_data[b_i];
            }

            FhatLessConc[i] = c_tmp;
          }
        } else {
          mtimes(Ce, b_data, FhatLessConc);
        }

        std::memcpy(&b_data[0], &input->disp.data[0], 12U * sizeof(double));
        mtimes(Khate, b_data, Ftemp_data);
        for (i = 0; i < 12; i++) {
          Fhate[i] = ((Fe[i] - Fhate[i]) - FhatLessConc[i]) - Ftemp_data[i];
        }
      }
    } else {
      if (f_strcmp(input->iterationType)) {
        // considerations if direct iteration is used or linear analysis
        for (i = 0; i < 12; i++) {
          b_data[i] = (1.0E+6 * input->disp.data[i] + 2000.0 * dispdot_data[i])
            + dispddot_data[i];
        }

        mtimes(Me, b_data, Fhate);
        for (i = 0; i < 12; i++) {
          b_data[i] = (1000.0 * input->disp.data[i] + dispdot_data[i]) + 0.0 *
            dispddot_data[i];
        }

        mtimes(Ce, b_data, FhatLessConc);
        for (i = 0; i < 12; i++) {
          Fhate[i] = (Fe[i] + Fhate[i]) + FhatLessConc[i];
        }
      }
    }

    for (i = 0; i < 144; i++) {
      Khate[i] = (Khate[i] + 1.0E+6 * Me[i]) + 1000.0 * Ce[i];
    }

    // ........................................................
    std::memcpy(&FhatLessConc[0], &Fhate[0], 12U * sizeof(double));
    std::memcpy(&Fe[0], &Fhate[0], 12U * sizeof(double));
  }

  // ----- assign output block ----------------
  std::memset(&output->FhatLessConc[0], 0, 12U * sizeof(double));
  std::memcpy(&output->Ke[0], &Khate[0], 144U * sizeof(double));
  std::memcpy(&output->Fe[0], &Fe[0], 12U * sizeof(double));
  output->Me.size[0] = 1;
  output->Me.size[1] = 1;
  output->Me.data[0] = 0.0;
  output->Ce.size[0] = 1;
  output->Ce.size[1] = 1;
  output->Ce.data[0] = 0.0;
  if (d_strcmp(input->analysisType)) {
    output->Me.size[0] = 12;
    output->Me.size[1] = 12;
    output->Ce.size[0] = 12;
    output->Ce.size[1] = 12;
    std::memcpy(&output->Me.data[0], &Me[0], 144U * sizeof(double));
    std::memcpy(&output->Ce.data[0], &Ce[0], 144U * sizeof(double));
  }

  if (c_strcmp(input->analysisType)) {
    std::memcpy(&output->FhatLessConc[0], &FhatLessConc[0], 12U * sizeof(double));
  }

  // ------------------------------------------
}

//
// File trailer for calculateTimoshenkoElementNL.cpp
//
// [EOF]
//
