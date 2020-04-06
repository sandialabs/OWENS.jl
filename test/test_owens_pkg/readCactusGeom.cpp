//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readCactusGeom.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "readCactusGeom.h"
#include "fileManager.h"
#include "myfgetl.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "test_owens_rtwutil.h"
#include <string.h>

// Function Declarations
static void processLine(double fid, emxArray_real_T *data);

// Function Definitions

//
// Arguments    : double fid
//                emxArray_real_T *data
// Return Type  : void
//
static void processLine(double fid, emxArray_real_T *data)
{
  emxArray_char_T *line;
  emxArray_boolean_T *delimiter_log_idx;
  int i;
  int loop_ub;
  double delimiter_idx_size;
  int n;
  emxArray_uint32_T *delimiter_idx;
  unsigned int j;
  int i1;
  emxArray_char_T *b_line;
  unsigned int u;
  int i2;
  int i3;
  creal_T dc;
  emxInit_char_T(&line, 2);
  emxInit_boolean_T(&delimiter_log_idx, 2);
  myfgetl(fid, line);

  //  Find where all of the delimiters are
  i = delimiter_log_idx->size[0] * delimiter_log_idx->size[1];
  delimiter_log_idx->size[0] = line->size[0];
  delimiter_log_idx->size[1] = line->size[1];
  emxEnsureCapacity_boolean_T(delimiter_log_idx, i);
  loop_ub = line->size[0] * line->size[1];
  for (i = 0; i < loop_ub; i++) {
    delimiter_log_idx->data[i] = (line->data[i] == ' ');
  }

  //  Reduce delimiter index to just two consecutive spaces
  delimiter_idx_size = 0.0;

  // Figure out the delimiter array size
  if ((delimiter_log_idx->size[0] == 0) || (delimiter_log_idx->size[1] == 0)) {
    n = -3;
  } else {
    n = delimiter_log_idx->size[1] - 3;
  }

  for (loop_ub = 0; loop_ub <= n; loop_ub++) {
    if (delimiter_log_idx->data[loop_ub] && delimiter_log_idx->data[loop_ub + 1]
        && (!delimiter_log_idx->data[loop_ub + 2])) {
      delimiter_idx_size++;
    }
  }

  emxInit_uint32_T(&delimiter_idx, 2);
  i = delimiter_idx->size[0] * delimiter_idx->size[1];
  delimiter_idx->size[0] = 1;
  loop_ub = static_cast<int>(delimiter_idx_size);
  delimiter_idx->size[1] = loop_ub;
  emxEnsureCapacity_uint32_T(delimiter_idx, i);
  for (i = 0; i < loop_ub; i++) {
    delimiter_idx->data[i] = 0U;
  }

  // Initialize variable scope
  j = 1U;
  if ((delimiter_log_idx->size[0] == 0) || (delimiter_log_idx->size[1] == 0)) {
    n = -3;
  } else {
    n = delimiter_log_idx->size[1] - 3;
  }

  for (loop_ub = 0; loop_ub <= n; loop_ub++) {
    if (delimiter_log_idx->data[loop_ub] && delimiter_log_idx->data[loop_ub + 1]
        && (!delimiter_log_idx->data[loop_ub + 2])) {
      delimiter_idx->data[static_cast<int>(j) - 1] = static_cast<unsigned int>
        ((loop_ub + 2));
      j++;
    }
  }

  emxFree_boolean_T(&delimiter_log_idx);
  if ((line->size[0] == 0) || (line->size[1] == 0)) {
    n = 0;
  } else {
    n = line->size[1];
  }

  i = delimiter_idx->size[1];
  i1 = delimiter_idx->size[0] * delimiter_idx->size[1];
  delimiter_idx->size[1]++;
  emxEnsureCapacity_uint32_T(delimiter_idx, i1);
  delimiter_idx->data[i] = n + 1U;

  // They all have a prefix in the .geom file, so don't include a 0 index and thus don't grab the prefix.  The logic above also skips the type and blade lines, so no need to parse for those either. 
  i = data->size[0];
  data->size[0] = delimiter_idx->size[1] - 1;
  emxEnsureCapacity_real_T(data, i);
  loop_ub = delimiter_idx->size[1];
  for (i = 0; i <= loop_ub - 2; i++) {
    data->data[i] = 0.0;
  }

  //  Extract the data from the beginning to the last delimiter
  i = delimiter_idx->size[1];
  emxInit_char_T(&b_line, 2);
  for (n = 0; n <= i - 2; n++) {
    j = delimiter_idx->data[n];
    u = delimiter_idx->data[n + 1];
    if (static_cast<double>(j) + 1.0 > static_cast<double>(u) - 1.0) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(j);
      i2 = static_cast<int>((static_cast<double>(u) - 1.0));
    }

    i3 = b_line->size[0] * b_line->size[1];
    b_line->size[0] = 1;
    loop_ub = i2 - i1;
    b_line->size[1] = loop_ub;
    emxEnsureCapacity_char_T(b_line, i3);
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_line->data[i2] = line->data[i1 + i2];
    }

    dc = c_str2double(b_line);
    data->data[n] = dc.re;
  }

  emxFree_char_T(&b_line);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);
}

//
// Arguments    : const emxArray_char_T *geomfn
//                int *cactusGeom_NBlade
//                int *cactusGeom_NStrut
//                emxArray_real_T *cactusGeom_RotN
//                emxArray_real_T *cactusGeom_RotP
//                emxArray_real_T *cactusGeom_RefAR
//                emxArray_real_T *cactusGeom_RefR
//                d_emxArray_struct_T *cactusGeom_blade
//                f_emxArray_struct_T *cactusGeom_strut
// Return Type  : void
//
void readCactusGeom(const emxArray_char_T *geomfn, int *cactusGeom_NBlade, int
                    *cactusGeom_NStrut, emxArray_real_T *cactusGeom_RotN,
                    emxArray_real_T *cactusGeom_RotP, emxArray_real_T
                    *cactusGeom_RefAR, emxArray_real_T *cactusGeom_RefR,
                    d_emxArray_struct_T *cactusGeom_blade, f_emxArray_struct_T
                    *cactusGeom_strut)
{
  emxArray_int32_T *NBlade;
  emxArray_real_T *FlipN;
  signed char fileid;
  int i;
  int loop_ub;
  emxArray_int32_T *NStrut;
  double d;
  int i1;
  emxArray_real_T *NElem;
  emxArray_int8_T *blade_QCx;
  emxArray_int8_T *blade_QCy;
  emxArray_int8_T *blade_QCz;
  emxArray_int8_T *blade_tx;
  emxArray_int8_T *blade_ty;
  emxArray_int8_T *blade_tz;
  emxArray_int8_T *blade_CtoR;
  int b_loop_ub;
  emxArray_int8_T *blade_PEx;
  emxArray_int8_T *blade_PEy;
  emxArray_int8_T *blade_PEz;
  emxArray_int8_T *blade_tEx;
  emxArray_int8_T *blade_tEy;
  emxArray_int8_T *blade_tEz;
  emxArray_int8_T *blade_nEx;
  emxArray_int8_T *blade_nEy;
  emxArray_int8_T *blade_nEz;
  emxArray_int8_T *blade_sEx;
  emxArray_int8_T *blade_sEy;
  emxArray_int8_T *blade_sEz;
  emxArray_int8_T *blade_ECtoR;
  emxArray_int8_T *blade_EAreaR;
  emxArray_int8_T *blade_iSect;
  k_struct_T expl_temp;
  int c_loop_ub;
  int d_loop_ub;
  int e_loop_ub;
  int f_loop_ub;
  int g_loop_ub;
  int h_loop_ub;
  int i_loop_ub;
  int j_loop_ub;
  int k_loop_ub;
  int l_loop_ub;
  int m_loop_ub;
  int n_loop_ub;
  int o_loop_ub;
  int p_loop_ub;
  int q_loop_ub;
  int r_loop_ub;
  int s_loop_ub;
  int t_loop_ub;
  int u_loop_ub;
  int v_loop_ub;
  int w_loop_ub;
  int x_loop_ub;
  m_struct_T b_expl_temp;
  int y_loop_ub;
  int ab_loop_ub;
  int bb_loop_ub;
  int cb_loop_ub;
  int db_loop_ub;
  int eb_loop_ub;
  int fb_loop_ub;
  int gb_loop_ub;
  int hb_loop_ub;
  int ib_loop_ub;
  int jb_loop_ub;
  int kb_loop_ub;
  emxInit_int32_T(&NBlade, 1);
  emxInit_real_T(&FlipN, 1);
  fileid = c_cfopen(geomfn, "rb");
  processLine(static_cast<double>(fileid), FlipN);
  i = NBlade->size[0];
  NBlade->size[0] = FlipN->size[0];
  emxEnsureCapacity_int32_T(NBlade, i);
  loop_ub = FlipN->size[0];
  for (i = 0; i < loop_ub; i++) {
    d = rt_roundd_snf(FlipN->data[i]);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i1 = static_cast<int>(d);
      } else {
        i1 = MIN_int32_T;
      }
    } else if (d >= 2.147483648E+9) {
      i1 = MAX_int32_T;
    } else {
      i1 = 0;
    }

    NBlade->data[i] = i1;
  }

  emxInit_int32_T(&NStrut, 1);
  *cactusGeom_NBlade = NBlade->data[0];
  processLine(static_cast<double>(fileid), FlipN);
  i = NStrut->size[0];
  NStrut->size[0] = FlipN->size[0];
  emxEnsureCapacity_int32_T(NStrut, i);
  loop_ub = FlipN->size[0];
  for (i = 0; i < loop_ub; i++) {
    d = rt_roundd_snf(FlipN->data[i]);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i1 = static_cast<int>(d);
      } else {
        i1 = MIN_int32_T;
      }
    } else if (d >= 2.147483648E+9) {
      i1 = MAX_int32_T;
    } else {
      i1 = 0;
    }

    NStrut->data[i] = i1;
  }

  emxInit_real_T(&NElem, 1);
  emxInit_int8_T(&blade_QCx, 1);
  *cactusGeom_NStrut = NStrut->data[0];
  processLine(static_cast<double>(fileid), cactusGeom_RotN);
  processLine(static_cast<double>(fileid), cactusGeom_RotP);
  processLine(static_cast<double>(fileid), cactusGeom_RefAR);
  processLine(static_cast<double>(fileid), cactusGeom_RefR);
  b_myfgetl(static_cast<double>(fileid));

  //  skip a line
  b_myfgetl(static_cast<double>(fileid));

  //  skip a line
  processLine(static_cast<double>(fileid), NElem);
  loop_ub = static_cast<int>((NElem->data[0] + 1.0));
  i = blade_QCx->size[0];
  blade_QCx->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCx, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCx->data[i] = 0;
  }

  emxInit_int8_T(&blade_QCy, 1);
  i = blade_QCy->size[0];
  blade_QCy->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCy, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCy->data[i] = 0;
  }

  emxInit_int8_T(&blade_QCz, 1);
  i = blade_QCz->size[0];
  blade_QCz->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCz, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCz->data[i] = 0;
  }

  emxInit_int8_T(&blade_tx, 1);
  i = blade_tx->size[0];
  blade_tx->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_tx, i);
  for (i = 0; i < loop_ub; i++) {
    blade_tx->data[i] = 0;
  }

  emxInit_int8_T(&blade_ty, 1);
  i = blade_ty->size[0];
  blade_ty->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_ty, i);
  for (i = 0; i < loop_ub; i++) {
    blade_ty->data[i] = 0;
  }

  emxInit_int8_T(&blade_tz, 1);
  i = blade_tz->size[0];
  blade_tz->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_tz, i);
  for (i = 0; i < loop_ub; i++) {
    blade_tz->data[i] = 0;
  }

  emxInit_int8_T(&blade_CtoR, 1);
  b_loop_ub = static_cast<int>(NElem->data[0]);
  i = blade_CtoR->size[0];
  blade_CtoR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_CtoR, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_CtoR->data[i] = 0;
  }

  emxInit_int8_T(&blade_PEx, 1);
  i = blade_PEx->size[0];
  blade_PEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEx->data[i] = 0;
  }

  emxInit_int8_T(&blade_PEy, 1);
  i = blade_PEy->size[0];
  blade_PEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEy->data[i] = 0;
  }

  emxInit_int8_T(&blade_PEz, 1);
  i = blade_PEz->size[0];
  blade_PEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEz->data[i] = 0;
  }

  emxInit_int8_T(&blade_tEx, 1);
  i = blade_tEx->size[0];
  blade_tEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tEx->data[i] = 0;
  }

  emxInit_int8_T(&blade_tEy, 1);
  i = blade_tEy->size[0];
  blade_tEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tEy->data[i] = 0;
  }

  emxInit_int8_T(&blade_tEz, 1);
  i = blade_tEz->size[0];
  blade_tEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tEz->data[i] = 0;
  }

  emxInit_int8_T(&blade_nEx, 1);
  i = blade_nEx->size[0];
  blade_nEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_nEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_nEx->data[i] = 0;
  }

  emxInit_int8_T(&blade_nEy, 1);
  i = blade_nEy->size[0];
  blade_nEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_nEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_nEy->data[i] = 0;
  }

  emxInit_int8_T(&blade_nEz, 1);
  i = blade_nEz->size[0];
  blade_nEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_nEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_nEz->data[i] = 0;
  }

  emxInit_int8_T(&blade_sEx, 1);
  i = blade_sEx->size[0];
  blade_sEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_sEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_sEx->data[i] = 0;
  }

  emxInit_int8_T(&blade_sEy, 1);
  i = blade_sEy->size[0];
  blade_sEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_sEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_sEy->data[i] = 0;
  }

  emxInit_int8_T(&blade_sEz, 1);
  i = blade_sEz->size[0];
  blade_sEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_sEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_sEz->data[i] = 0;
  }

  emxInit_int8_T(&blade_ECtoR, 1);
  i = blade_ECtoR->size[0];
  blade_ECtoR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_ECtoR, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_ECtoR->data[i] = 0;
  }

  emxInit_int8_T(&blade_EAreaR, 1);
  i = blade_EAreaR->size[0];
  blade_EAreaR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_EAreaR, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_EAreaR->data[i] = 0;
  }

  emxInit_int8_T(&blade_iSect, 1);
  i = blade_iSect->size[0];
  blade_iSect->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_iSect, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_iSect->data[i] = 0;
  }

  emxInitStruct_struct_T5(&expl_temp);
  i = expl_temp.iSect->size[0];
  expl_temp.iSect->size[0] = blade_iSect->size[0];
  emxEnsureCapacity_real_T(expl_temp.iSect, i);
  b_loop_ub = blade_iSect->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.iSect->data[i] = blade_iSect->data[i];
  }

  emxFree_int8_T(&blade_iSect);
  i = expl_temp.EAreaR->size[0];
  expl_temp.EAreaR->size[0] = blade_EAreaR->size[0];
  emxEnsureCapacity_real_T(expl_temp.EAreaR, i);
  b_loop_ub = blade_EAreaR->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.EAreaR->data[i] = blade_EAreaR->data[i];
  }

  emxFree_int8_T(&blade_EAreaR);
  i = expl_temp.ECtoR->size[0];
  expl_temp.ECtoR->size[0] = blade_ECtoR->size[0];
  emxEnsureCapacity_real_T(expl_temp.ECtoR, i);
  b_loop_ub = blade_ECtoR->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.ECtoR->data[i] = blade_ECtoR->data[i];
  }

  emxFree_int8_T(&blade_ECtoR);
  i = expl_temp.sEz->size[0];
  expl_temp.sEz->size[0] = blade_sEz->size[0];
  emxEnsureCapacity_real_T(expl_temp.sEz, i);
  b_loop_ub = blade_sEz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.sEz->data[i] = blade_sEz->data[i];
  }

  emxFree_int8_T(&blade_sEz);
  i = expl_temp.sEy->size[0];
  expl_temp.sEy->size[0] = blade_sEy->size[0];
  emxEnsureCapacity_real_T(expl_temp.sEy, i);
  b_loop_ub = blade_sEy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.sEy->data[i] = blade_sEy->data[i];
  }

  emxFree_int8_T(&blade_sEy);
  i = expl_temp.sEx->size[0];
  expl_temp.sEx->size[0] = blade_sEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.sEx, i);
  b_loop_ub = blade_sEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.sEx->data[i] = blade_sEx->data[i];
  }

  emxFree_int8_T(&blade_sEx);
  i = expl_temp.nEz->size[0];
  expl_temp.nEz->size[0] = blade_nEz->size[0];
  emxEnsureCapacity_real_T(expl_temp.nEz, i);
  b_loop_ub = blade_nEz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.nEz->data[i] = blade_nEz->data[i];
  }

  emxFree_int8_T(&blade_nEz);
  i = expl_temp.nEy->size[0];
  expl_temp.nEy->size[0] = blade_nEy->size[0];
  emxEnsureCapacity_real_T(expl_temp.nEy, i);
  b_loop_ub = blade_nEy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.nEy->data[i] = blade_nEy->data[i];
  }

  emxFree_int8_T(&blade_nEy);
  i = expl_temp.nEx->size[0];
  expl_temp.nEx->size[0] = blade_nEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.nEx, i);
  b_loop_ub = blade_nEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.nEx->data[i] = blade_nEx->data[i];
  }

  emxFree_int8_T(&blade_nEx);
  i = expl_temp.tEz->size[0];
  expl_temp.tEz->size[0] = blade_tEz->size[0];
  emxEnsureCapacity_real_T(expl_temp.tEz, i);
  b_loop_ub = blade_tEz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tEz->data[i] = blade_tEz->data[i];
  }

  emxFree_int8_T(&blade_tEz);
  i = expl_temp.tEy->size[0];
  expl_temp.tEy->size[0] = blade_tEy->size[0];
  emxEnsureCapacity_real_T(expl_temp.tEy, i);
  b_loop_ub = blade_tEy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tEy->data[i] = blade_tEy->data[i];
  }

  i = expl_temp.tEx->size[0];
  expl_temp.tEx->size[0] = blade_tEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.tEx, i);
  b_loop_ub = blade_tEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tEx->data[i] = blade_tEx->data[i];
  }

  i = expl_temp.PEz->size[0];
  expl_temp.PEz->size[0] = blade_PEz->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEz, i);
  b_loop_ub = blade_PEz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEz->data[i] = blade_PEz->data[i];
  }

  i = expl_temp.PEy->size[0];
  expl_temp.PEy->size[0] = blade_PEy->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEy, i);
  b_loop_ub = blade_PEy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEy->data[i] = blade_PEy->data[i];
  }

  i = expl_temp.PEx->size[0];
  expl_temp.PEx->size[0] = blade_PEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEx, i);
  b_loop_ub = blade_PEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEx->data[i] = blade_PEx->data[i];
  }

  i = expl_temp.CtoR->size[0];
  expl_temp.CtoR->size[0] = blade_CtoR->size[0];
  emxEnsureCapacity_real_T(expl_temp.CtoR, i);
  b_loop_ub = blade_CtoR->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.CtoR->data[i] = blade_CtoR->data[i];
  }

  i = expl_temp.tz->size[0];
  expl_temp.tz->size[0] = blade_tz->size[0];
  emxEnsureCapacity_real_T(expl_temp.tz, i);
  b_loop_ub = blade_tz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tz->data[i] = blade_tz->data[i];
  }

  i = expl_temp.ty->size[0];
  expl_temp.ty->size[0] = blade_ty->size[0];
  emxEnsureCapacity_real_T(expl_temp.ty, i);
  b_loop_ub = blade_ty->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.ty->data[i] = blade_ty->data[i];
  }

  i = expl_temp.tx->size[0];
  expl_temp.tx->size[0] = blade_tx->size[0];
  emxEnsureCapacity_real_T(expl_temp.tx, i);
  b_loop_ub = blade_tx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tx->data[i] = blade_tx->data[i];
  }

  i = expl_temp.QCz->size[0];
  expl_temp.QCz->size[0] = blade_QCz->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCz, i);
  b_loop_ub = blade_QCz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCz->data[i] = blade_QCz->data[i];
  }

  i = expl_temp.QCy->size[0];
  expl_temp.QCy->size[0] = blade_QCy->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCy, i);
  b_loop_ub = blade_QCy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCy->data[i] = blade_QCy->data[i];
  }

  i = expl_temp.QCx->size[0];
  expl_temp.QCx->size[0] = blade_QCx->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCx, i);
  b_loop_ub = blade_QCx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCx->data[i] = blade_QCx->data[i];
  }

  expl_temp.FlipN = 0.0;
  expl_temp.NElem = 0.0;
  repmat(&expl_temp, NBlade->data[0], cactusGeom_blade);
  i = NBlade->data[0];
  emxFreeStruct_struct_T4(&expl_temp);
  if (0 <= NBlade->data[0] - 1) {
    if (1.0 > NElem->data[0] + 1.0) {
      loop_ub = 0;
    }

    c_loop_ub = loop_ub;
    if (1.0 > NElem->data[0] + 1.0) {
      d_loop_ub = 0;
    } else {
      d_loop_ub = static_cast<int>((NElem->data[0] + 1.0));
    }

    if (1.0 > NElem->data[0] + 1.0) {
      e_loop_ub = 0;
    } else {
      e_loop_ub = static_cast<int>((NElem->data[0] + 1.0));
    }

    if (1.0 > NElem->data[0] + 1.0) {
      f_loop_ub = 0;
    } else {
      f_loop_ub = static_cast<int>((NElem->data[0] + 1.0));
    }

    if (1.0 > NElem->data[0] + 1.0) {
      g_loop_ub = 0;
    } else {
      g_loop_ub = static_cast<int>((NElem->data[0] + 1.0));
    }

    if (1.0 > NElem->data[0] + 1.0) {
      h_loop_ub = 0;
    } else {
      h_loop_ub = static_cast<int>((NElem->data[0] + 1.0));
    }

    if (1.0 > NElem->data[0]) {
      i_loop_ub = 0;
    } else {
      i_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      j_loop_ub = 0;
    } else {
      j_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      k_loop_ub = 0;
    } else {
      k_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      l_loop_ub = 0;
    } else {
      l_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      m_loop_ub = 0;
    } else {
      m_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      n_loop_ub = 0;
    } else {
      n_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      o_loop_ub = 0;
    } else {
      o_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      p_loop_ub = 0;
    } else {
      p_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      q_loop_ub = 0;
    } else {
      q_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      r_loop_ub = 0;
    } else {
      r_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      s_loop_ub = 0;
    } else {
      s_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      t_loop_ub = 0;
    } else {
      t_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      u_loop_ub = 0;
    } else {
      u_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      v_loop_ub = 0;
    } else {
      v_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      w_loop_ub = 0;
    } else {
      w_loop_ub = static_cast<int>(NElem->data[0]);
    }

    if (1.0 > NElem->data[0]) {
      x_loop_ub = 0;
    } else {
      x_loop_ub = static_cast<int>(NElem->data[0]);
    }
  }

  emxFree_int32_T(&NBlade);
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    if (b_loop_ub + 1 != 1) {
      b_myfgetl(static_cast<double>(fileid));

      //  skip a line
      b_myfgetl(static_cast<double>(fileid));

      //  skip a line
    }

    cactusGeom_blade->data[b_loop_ub].NElem = NElem->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_blade->data[b_loop_ub].FlipN = FlipN->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].QCx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < d_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].QCy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < e_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].QCz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < f_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].tx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < g_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].ty->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < h_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].tz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < i_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].CtoR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < j_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].PEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < k_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].PEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < l_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].PEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < m_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].tEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < n_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].tEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < o_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].tEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < p_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].nEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < q_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].nEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < r_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].nEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < s_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].sEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < t_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].sEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < u_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].sEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < v_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].ECtoR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < w_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].EAreaR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < x_loop_ub; i1++) {
      cactusGeom_blade->data[b_loop_ub].iSect->data[i1] = FlipN->data[i1];
    }
  }

  b_myfgetl(static_cast<double>(fileid));

  //  skip line
  processLine(static_cast<double>(fileid), NElem);
  loop_ub = static_cast<int>((NElem->data[0] + 1.0));
  i = blade_QCx->size[0];
  blade_QCx->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCx, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCx->data[i] = 0;
  }

  i = blade_QCy->size[0];
  blade_QCy->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCy, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCy->data[i] = 0;
  }

  i = blade_QCz->size[0];
  blade_QCz->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_QCz, i);
  for (i = 0; i < loop_ub; i++) {
    blade_QCz->data[i] = 0;
  }

  i = blade_tx->size[0];
  blade_tx->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(blade_tx, i);
  for (i = 0; i < loop_ub; i++) {
    blade_tx->data[i] = 0;
  }

  b_loop_ub = static_cast<int>(NElem->data[0]);
  i = blade_ty->size[0];
  blade_ty->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_ty, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_ty->data[i] = 0;
  }

  i = blade_tz->size[0];
  blade_tz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tz->data[i] = 0;
  }

  i = blade_CtoR->size[0];
  blade_CtoR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_CtoR, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_CtoR->data[i] = 0;
  }

  i = blade_PEx->size[0];
  blade_PEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEx->data[i] = 0;
  }

  i = blade_PEy->size[0];
  blade_PEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEy->data[i] = 0;
  }

  i = blade_PEz->size[0];
  blade_PEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_PEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_PEz->data[i] = 0;
  }

  i = blade_tEx->size[0];
  blade_tEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tEx->data[i] = 0;
  }

  i = blade_tEy->size[0];
  blade_tEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(blade_tEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    blade_tEy->data[i] = 0;
  }

  emxInitStruct_struct_T6(&b_expl_temp);
  b_expl_temp.EIndE = 0.0;
  b_expl_temp.BIndE = 0.0;
  b_expl_temp.EIndS = 0.0;
  b_expl_temp.BIndS = 0.0;
  i = b_expl_temp.EAreaR->size[0];
  b_expl_temp.EAreaR->size[0] = blade_tEy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.EAreaR, i);
  c_loop_ub = blade_tEy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.EAreaR->data[i] = blade_tEy->data[i];
  }

  emxFree_int8_T(&blade_tEy);
  i = b_expl_temp.ECtoR->size[0];
  b_expl_temp.ECtoR->size[0] = blade_tEx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.ECtoR, i);
  c_loop_ub = blade_tEx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.ECtoR->data[i] = blade_tEx->data[i];
  }

  emxFree_int8_T(&blade_tEx);
  i = b_expl_temp.sEz->size[0];
  b_expl_temp.sEz->size[0] = blade_PEz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEz, i);
  c_loop_ub = blade_PEz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEz->data[i] = blade_PEz->data[i];
  }

  emxFree_int8_T(&blade_PEz);
  i = b_expl_temp.sEy->size[0];
  b_expl_temp.sEy->size[0] = blade_PEy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEy, i);
  c_loop_ub = blade_PEy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEy->data[i] = blade_PEy->data[i];
  }

  emxFree_int8_T(&blade_PEy);
  i = b_expl_temp.sEx->size[0];
  b_expl_temp.sEx->size[0] = blade_PEx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEx, i);
  c_loop_ub = blade_PEx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEx->data[i] = blade_PEx->data[i];
  }

  emxFree_int8_T(&blade_PEx);
  i = b_expl_temp.PEz->size[0];
  b_expl_temp.PEz->size[0] = blade_CtoR->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEz, i);
  c_loop_ub = blade_CtoR->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEz->data[i] = blade_CtoR->data[i];
  }

  emxFree_int8_T(&blade_CtoR);
  i = b_expl_temp.PEy->size[0];
  b_expl_temp.PEy->size[0] = blade_tz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEy, i);
  c_loop_ub = blade_tz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEy->data[i] = blade_tz->data[i];
  }

  emxFree_int8_T(&blade_tz);
  i = b_expl_temp.PEx->size[0];
  b_expl_temp.PEx->size[0] = blade_ty->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEx, i);
  c_loop_ub = blade_ty->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEx->data[i] = blade_ty->data[i];
  }

  emxFree_int8_T(&blade_ty);
  i = b_expl_temp.CtoR->size[0];
  b_expl_temp.CtoR->size[0] = blade_tx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.CtoR, i);
  c_loop_ub = blade_tx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.CtoR->data[i] = blade_tx->data[i];
  }

  emxFree_int8_T(&blade_tx);
  i = b_expl_temp.MCz->size[0];
  b_expl_temp.MCz->size[0] = blade_QCz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCz, i);
  c_loop_ub = blade_QCz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCz->data[i] = blade_QCz->data[i];
  }

  emxFree_int8_T(&blade_QCz);
  i = b_expl_temp.MCy->size[0];
  b_expl_temp.MCy->size[0] = blade_QCy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCy, i);
  c_loop_ub = blade_QCy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCy->data[i] = blade_QCy->data[i];
  }

  emxFree_int8_T(&blade_QCy);
  i = b_expl_temp.MCx->size[0];
  b_expl_temp.MCx->size[0] = blade_QCx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCx, i);
  c_loop_ub = blade_QCx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCx->data[i] = blade_QCx->data[i];
  }

  emxFree_int8_T(&blade_QCx);
  b_expl_temp.TtoC = 0.0;
  b_expl_temp.NElem = 0.0;
  b_repmat(&b_expl_temp, NStrut->data[0], cactusGeom_strut);
  i = NStrut->data[0];
  emxFreeStruct_struct_T6(&b_expl_temp);
  if (0 <= NStrut->data[0] - 1) {
    if (1.0 > NElem->data[0] + 1.0) {
      y_loop_ub = 0;
    } else {
      y_loop_ub = loop_ub;
    }

    if (1.0 > NElem->data[0] + 1.0) {
      ab_loop_ub = 0;
    } else {
      ab_loop_ub = loop_ub;
    }

    if (1.0 > NElem->data[0] + 1.0) {
      bb_loop_ub = 0;
    } else {
      bb_loop_ub = loop_ub;
    }

    if (1.0 > NElem->data[0] + 1.0) {
      loop_ub = 0;
    }

    cb_loop_ub = loop_ub;
    if (1.0 > NElem->data[0]) {
      db_loop_ub = 0;
    } else {
      db_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      eb_loop_ub = 0;
    } else {
      eb_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      fb_loop_ub = 0;
    } else {
      fb_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      gb_loop_ub = 0;
    } else {
      gb_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      hb_loop_ub = 0;
    } else {
      hb_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      ib_loop_ub = 0;
    } else {
      ib_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      jb_loop_ub = 0;
    } else {
      jb_loop_ub = b_loop_ub;
    }

    if (1.0 > NElem->data[0]) {
      b_loop_ub = 0;
    }

    kb_loop_ub = b_loop_ub;
  }

  emxFree_int32_T(&NStrut);
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    if (b_loop_ub + 1 != 1) {
      b_myfgetl(static_cast<double>(fileid));

      //  skip line
      b_myfgetl(static_cast<double>(fileid));

      //  skip line
    }

    cactusGeom_strut->data[b_loop_ub].NElem = NElem->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_strut->data[b_loop_ub].TtoC = FlipN->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < y_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].MCx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < ab_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].MCy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < bb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].MCz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < cb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].CtoR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < db_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].PEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < eb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].PEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < fb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].PEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < gb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].sEx->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < hb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].sEy->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < ib_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].sEz->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < jb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].ECtoR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    for (i1 = 0; i1 < kb_loop_ub; i1++) {
      cactusGeom_strut->data[b_loop_ub].EAreaR->data[i1] = FlipN->data[i1];
    }

    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_strut->data[b_loop_ub].BIndS = FlipN->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_strut->data[b_loop_ub].EIndS = FlipN->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_strut->data[b_loop_ub].BIndE = FlipN->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_strut->data[b_loop_ub].EIndE = FlipN->data[0];
  }

  emxFree_real_T(&FlipN);
  emxFree_real_T(&NElem);
  cfclose(static_cast<double>(fileid));
}

//
// File trailer for readCactusGeom.cpp
//
// [EOF]
//
