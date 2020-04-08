//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readCactusGeom.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
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
  int b_i;
  emxArray_uint32_T *delimiter_idx;
  unsigned int j;
  int i1;
  emxArray_char_T *b_line;
  int k;
  unsigned int u;
  unsigned int u1;
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
    loop_ub = -3;
  } else {
    loop_ub = delimiter_log_idx->size[1] - 3;
  }

  for (b_i = 0; b_i <= loop_ub; b_i++) {
    if (delimiter_log_idx->data[b_i] && delimiter_log_idx->data[b_i + 1] &&
        (!delimiter_log_idx->data[b_i + 2])) {
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
    loop_ub = -3;
  } else {
    loop_ub = delimiter_log_idx->size[1] - 3;
  }

  for (b_i = 0; b_i <= loop_ub; b_i++) {
    if (delimiter_log_idx->data[b_i] && delimiter_log_idx->data[b_i + 1] &&
        (!delimiter_log_idx->data[b_i + 2])) {
      delimiter_idx->data[static_cast<int>(j) - 1] = static_cast<unsigned int>
        ((b_i + 2));
      j++;
    }
  }

  emxFree_boolean_T(&delimiter_log_idx);
  if ((line->size[0] == 0) || (line->size[1] == 0)) {
    loop_ub = 0;
  } else {
    loop_ub = line->size[1];
  }

  i = delimiter_idx->size[1];
  i1 = delimiter_idx->size[0] * delimiter_idx->size[1];
  delimiter_idx->size[1]++;
  emxEnsureCapacity_uint32_T(delimiter_idx, i1);
  delimiter_idx->data[i] = loop_ub + 1U;

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
  for (k = 0; k <= i - 2; k++) {
    u = delimiter_idx->data[k];
    u1 = delimiter_idx->data[k + 1];
    if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(u);
      i2 = static_cast<int>((static_cast<double>(u1) - 1.0));
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
    data->data[k] = dc.re;
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
  emxArray_real_T *b_r;
  signed char fileid;
  int i;
  int loop_ub;
  emxArray_int32_T *NStrut;
  double d;
  int b_loop_ub;
  emxArray_real_T *NElem;
  emxArray_int8_T *blade_QCx;
  emxArray_int8_T *blade_QCy;
  emxArray_int8_T *blade_QCz;
  emxArray_int8_T *blade_tx;
  emxArray_int8_T *blade_ty;
  emxArray_int8_T *blade_tz;
  emxArray_int8_T *blade_CtoR;
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
  emxArray_real_T *FlipN;
  int c_loop_ub;
  emxArray_real_T *QCx;
  emxArray_real_T *QCy;
  int d_loop_ub;
  emxArray_real_T *QCz;
  emxArray_real_T *tx;
  emxArray_real_T *ty;
  int e_loop_ub;
  emxArray_real_T *tz;
  emxArray_real_T *CtoR;
  emxArray_real_T *PEx;
  int f_loop_ub;
  emxArray_real_T *PEy;
  emxArray_real_T *PEz;
  emxArray_real_T *tEx;
  int g_loop_ub;
  emxArray_real_T *tEy;
  emxArray_real_T *tEz;
  emxArray_real_T *nEx;
  int h_loop_ub;
  emxArray_real_T *nEy;
  emxArray_real_T *nEz;
  emxArray_real_T *sEx;
  int i_loop_ub;
  emxArray_real_T *sEy;
  emxArray_real_T *sEz;
  emxArray_real_T *ECtoR;
  int j_loop_ub;
  emxArray_real_T *EAreaR;
  emxArray_real_T *iSect;
  int b_i;
  int k_loop_ub;
  int l_loop_ub;
  int m_loop_ub;
  int n_loop_ub;
  int o_loop_ub;
  int p_loop_ub;
  emxArray_int8_T *strut_MCx;
  int q_loop_ub;
  int r_loop_ub;
  int s_loop_ub;
  emxArray_int8_T *strut_MCy;
  int t_loop_ub;
  int u_loop_ub;
  emxArray_int8_T *strut_MCz;
  int v_loop_ub;
  int w_loop_ub;
  emxArray_int8_T *strut_CtoR;
  int x_loop_ub;
  emxArray_int8_T *strut_PEx;
  emxArray_int8_T *strut_PEy;
  emxArray_int8_T *strut_PEz;
  emxArray_int8_T *strut_sEx;
  emxArray_int8_T *strut_sEy;
  emxArray_int8_T *strut_sEz;
  emxArray_int8_T *strut_ECtoR;
  emxArray_int8_T *strut_EAreaR;
  m_struct_T b_expl_temp;
  int y_loop_ub;
  emxArray_real_T *TtoC;
  emxArray_real_T *MCx;
  emxArray_real_T *MCy;
  int ab_loop_ub;
  emxArray_real_T *MCz;
  emxArray_real_T *BIndS;
  emxArray_real_T *EIndS;
  int bb_loop_ub;
  emxArray_real_T *BIndE;
  emxArray_real_T *EIndE;
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
  emxInit_real_T(&b_r, 1);
  fileid = c_cfopen(geomfn, "rb");
  processLine(static_cast<double>(fileid), b_r);
  i = NBlade->size[0];
  NBlade->size[0] = b_r->size[0];
  emxEnsureCapacity_int32_T(NBlade, i);
  loop_ub = b_r->size[0];
  for (i = 0; i < loop_ub; i++) {
    d = rt_roundd_snf(b_r->data[i]);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        b_loop_ub = static_cast<int>(d);
      } else {
        b_loop_ub = MIN_int32_T;
      }
    } else if (d >= 2.147483648E+9) {
      b_loop_ub = MAX_int32_T;
    } else {
      b_loop_ub = 0;
    }

    NBlade->data[i] = b_loop_ub;
  }

  emxInit_int32_T(&NStrut, 1);
  *cactusGeom_NBlade = NBlade->data[0];
  processLine(static_cast<double>(fileid), b_r);
  i = NStrut->size[0];
  NStrut->size[0] = b_r->size[0];
  emxEnsureCapacity_int32_T(NStrut, i);
  loop_ub = b_r->size[0];
  for (i = 0; i < loop_ub; i++) {
    d = rt_roundd_snf(b_r->data[i]);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        b_loop_ub = static_cast<int>(d);
      } else {
        b_loop_ub = MIN_int32_T;
      }
    } else if (d >= 2.147483648E+9) {
      b_loop_ub = MAX_int32_T;
    } else {
      b_loop_ub = 0;
    }

    NStrut->data[i] = b_loop_ub;
  }

  emxFree_real_T(&b_r);
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

  emxFree_int8_T(&blade_tEy);
  i = expl_temp.tEx->size[0];
  expl_temp.tEx->size[0] = blade_tEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.tEx, i);
  b_loop_ub = blade_tEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tEx->data[i] = blade_tEx->data[i];
  }

  emxFree_int8_T(&blade_tEx);
  i = expl_temp.PEz->size[0];
  expl_temp.PEz->size[0] = blade_PEz->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEz, i);
  b_loop_ub = blade_PEz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEz->data[i] = blade_PEz->data[i];
  }

  emxFree_int8_T(&blade_PEz);
  i = expl_temp.PEy->size[0];
  expl_temp.PEy->size[0] = blade_PEy->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEy, i);
  b_loop_ub = blade_PEy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEy->data[i] = blade_PEy->data[i];
  }

  emxFree_int8_T(&blade_PEy);
  i = expl_temp.PEx->size[0];
  expl_temp.PEx->size[0] = blade_PEx->size[0];
  emxEnsureCapacity_real_T(expl_temp.PEx, i);
  b_loop_ub = blade_PEx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.PEx->data[i] = blade_PEx->data[i];
  }

  emxFree_int8_T(&blade_PEx);
  i = expl_temp.CtoR->size[0];
  expl_temp.CtoR->size[0] = blade_CtoR->size[0];
  emxEnsureCapacity_real_T(expl_temp.CtoR, i);
  b_loop_ub = blade_CtoR->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.CtoR->data[i] = blade_CtoR->data[i];
  }

  emxFree_int8_T(&blade_CtoR);
  i = expl_temp.tz->size[0];
  expl_temp.tz->size[0] = blade_tz->size[0];
  emxEnsureCapacity_real_T(expl_temp.tz, i);
  b_loop_ub = blade_tz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tz->data[i] = blade_tz->data[i];
  }

  emxFree_int8_T(&blade_tz);
  i = expl_temp.ty->size[0];
  expl_temp.ty->size[0] = blade_ty->size[0];
  emxEnsureCapacity_real_T(expl_temp.ty, i);
  b_loop_ub = blade_ty->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.ty->data[i] = blade_ty->data[i];
  }

  emxFree_int8_T(&blade_ty);
  i = expl_temp.tx->size[0];
  expl_temp.tx->size[0] = blade_tx->size[0];
  emxEnsureCapacity_real_T(expl_temp.tx, i);
  b_loop_ub = blade_tx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.tx->data[i] = blade_tx->data[i];
  }

  emxFree_int8_T(&blade_tx);
  i = expl_temp.QCz->size[0];
  expl_temp.QCz->size[0] = blade_QCz->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCz, i);
  b_loop_ub = blade_QCz->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCz->data[i] = blade_QCz->data[i];
  }

  emxFree_int8_T(&blade_QCz);
  i = expl_temp.QCy->size[0];
  expl_temp.QCy->size[0] = blade_QCy->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCy, i);
  b_loop_ub = blade_QCy->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCy->data[i] = blade_QCy->data[i];
  }

  emxFree_int8_T(&blade_QCy);
  i = expl_temp.QCx->size[0];
  expl_temp.QCx->size[0] = blade_QCx->size[0];
  emxEnsureCapacity_real_T(expl_temp.QCx, i);
  b_loop_ub = blade_QCx->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    expl_temp.QCx->data[i] = blade_QCx->data[i];
  }

  emxFree_int8_T(&blade_QCx);
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
  emxInit_real_T(&FlipN, 1);
  emxInit_real_T(&QCx, 1);
  emxInit_real_T(&QCy, 1);
  emxInit_real_T(&QCz, 1);
  emxInit_real_T(&tx, 1);
  emxInit_real_T(&ty, 1);
  emxInit_real_T(&tz, 1);
  emxInit_real_T(&CtoR, 1);
  emxInit_real_T(&PEx, 1);
  emxInit_real_T(&PEy, 1);
  emxInit_real_T(&PEz, 1);
  emxInit_real_T(&tEx, 1);
  emxInit_real_T(&tEy, 1);
  emxInit_real_T(&tEz, 1);
  emxInit_real_T(&nEx, 1);
  emxInit_real_T(&nEy, 1);
  emxInit_real_T(&nEz, 1);
  emxInit_real_T(&sEx, 1);
  emxInit_real_T(&sEy, 1);
  emxInit_real_T(&sEz, 1);
  emxInit_real_T(&ECtoR, 1);
  emxInit_real_T(&EAreaR, 1);
  emxInit_real_T(&iSect, 1);
  for (b_i = 0; b_i < i; b_i++) {
    if (b_i + 1 != 1) {
      b_myfgetl(static_cast<double>(fileid));

      //  skip a line
      b_myfgetl(static_cast<double>(fileid));

      //  skip a line
    }

    cactusGeom_blade->data[b_i].NElem = NElem->data[0];
    processLine(static_cast<double>(fileid), FlipN);
    cactusGeom_blade->data[b_i].FlipN = FlipN->data[0];
    processLine(static_cast<double>(fileid), QCx);
    for (b_loop_ub = 0; b_loop_ub < c_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].QCx->data[b_loop_ub] = QCx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), QCy);
    for (b_loop_ub = 0; b_loop_ub < d_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].QCy->data[b_loop_ub] = QCy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), QCz);
    for (b_loop_ub = 0; b_loop_ub < e_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].QCz->data[b_loop_ub] = QCz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), tx);
    for (b_loop_ub = 0; b_loop_ub < f_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].tx->data[b_loop_ub] = tx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), ty);
    for (b_loop_ub = 0; b_loop_ub < g_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].ty->data[b_loop_ub] = ty->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), tz);
    for (b_loop_ub = 0; b_loop_ub < h_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].tz->data[b_loop_ub] = tz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), CtoR);
    for (b_loop_ub = 0; b_loop_ub < i_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].CtoR->data[b_loop_ub] = CtoR->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEx);
    for (b_loop_ub = 0; b_loop_ub < j_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].PEx->data[b_loop_ub] = PEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEy);
    for (b_loop_ub = 0; b_loop_ub < k_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].PEy->data[b_loop_ub] = PEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEz);
    for (b_loop_ub = 0; b_loop_ub < l_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].PEz->data[b_loop_ub] = PEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), tEx);
    for (b_loop_ub = 0; b_loop_ub < m_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].tEx->data[b_loop_ub] = tEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), tEy);
    for (b_loop_ub = 0; b_loop_ub < n_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].tEy->data[b_loop_ub] = tEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), tEz);
    for (b_loop_ub = 0; b_loop_ub < o_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].tEz->data[b_loop_ub] = tEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), nEx);
    for (b_loop_ub = 0; b_loop_ub < p_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].nEx->data[b_loop_ub] = nEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), nEy);
    for (b_loop_ub = 0; b_loop_ub < q_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].nEy->data[b_loop_ub] = nEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), nEz);
    for (b_loop_ub = 0; b_loop_ub < r_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].nEz->data[b_loop_ub] = nEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEx);
    for (b_loop_ub = 0; b_loop_ub < s_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].sEx->data[b_loop_ub] = sEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEy);
    for (b_loop_ub = 0; b_loop_ub < t_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].sEy->data[b_loop_ub] = sEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEz);
    for (b_loop_ub = 0; b_loop_ub < u_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].sEz->data[b_loop_ub] = sEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), ECtoR);
    for (b_loop_ub = 0; b_loop_ub < v_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].ECtoR->data[b_loop_ub] = ECtoR->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), EAreaR);
    for (b_loop_ub = 0; b_loop_ub < w_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].EAreaR->data[b_loop_ub] = EAreaR->
        data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), iSect);
    for (b_loop_ub = 0; b_loop_ub < x_loop_ub; b_loop_ub++) {
      cactusGeom_blade->data[b_i].iSect->data[b_loop_ub] = iSect->data[b_loop_ub];
    }
  }

  emxFree_real_T(&iSect);
  emxFree_real_T(&nEz);
  emxFree_real_T(&nEy);
  emxFree_real_T(&nEx);
  emxFree_real_T(&tEz);
  emxFree_real_T(&tEy);
  emxFree_real_T(&tEx);
  emxFree_real_T(&tz);
  emxFree_real_T(&ty);
  emxFree_real_T(&tx);
  emxFree_real_T(&QCz);
  emxFree_real_T(&QCy);
  emxFree_real_T(&QCx);
  emxFree_real_T(&FlipN);
  emxInit_int8_T(&strut_MCx, 1);
  b_myfgetl(static_cast<double>(fileid));

  //  skip line
  processLine(static_cast<double>(fileid), NElem);
  loop_ub = static_cast<int>((NElem->data[0] + 1.0));
  i = strut_MCx->size[0];
  strut_MCx->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(strut_MCx, i);
  for (i = 0; i < loop_ub; i++) {
    strut_MCx->data[i] = 0;
  }

  emxInit_int8_T(&strut_MCy, 1);
  i = strut_MCy->size[0];
  strut_MCy->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(strut_MCy, i);
  for (i = 0; i < loop_ub; i++) {
    strut_MCy->data[i] = 0;
  }

  emxInit_int8_T(&strut_MCz, 1);
  i = strut_MCz->size[0];
  strut_MCz->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(strut_MCz, i);
  for (i = 0; i < loop_ub; i++) {
    strut_MCz->data[i] = 0;
  }

  emxInit_int8_T(&strut_CtoR, 1);
  i = strut_CtoR->size[0];
  strut_CtoR->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(strut_CtoR, i);
  for (i = 0; i < loop_ub; i++) {
    strut_CtoR->data[i] = 0;
  }

  emxInit_int8_T(&strut_PEx, 1);
  b_loop_ub = static_cast<int>(NElem->data[0]);
  i = strut_PEx->size[0];
  strut_PEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_PEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_PEx->data[i] = 0;
  }

  emxInit_int8_T(&strut_PEy, 1);
  i = strut_PEy->size[0];
  strut_PEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_PEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_PEy->data[i] = 0;
  }

  emxInit_int8_T(&strut_PEz, 1);
  i = strut_PEz->size[0];
  strut_PEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_PEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_PEz->data[i] = 0;
  }

  emxInit_int8_T(&strut_sEx, 1);
  i = strut_sEx->size[0];
  strut_sEx->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_sEx, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_sEx->data[i] = 0;
  }

  emxInit_int8_T(&strut_sEy, 1);
  i = strut_sEy->size[0];
  strut_sEy->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_sEy, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_sEy->data[i] = 0;
  }

  emxInit_int8_T(&strut_sEz, 1);
  i = strut_sEz->size[0];
  strut_sEz->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_sEz, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_sEz->data[i] = 0;
  }

  emxInit_int8_T(&strut_ECtoR, 1);
  i = strut_ECtoR->size[0];
  strut_ECtoR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_ECtoR, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_ECtoR->data[i] = 0;
  }

  emxInit_int8_T(&strut_EAreaR, 1);
  i = strut_EAreaR->size[0];
  strut_EAreaR->size[0] = b_loop_ub;
  emxEnsureCapacity_int8_T(strut_EAreaR, i);
  for (i = 0; i < b_loop_ub; i++) {
    strut_EAreaR->data[i] = 0;
  }

  emxInitStruct_struct_T6(&b_expl_temp);
  b_expl_temp.EIndE = 0.0;
  b_expl_temp.BIndE = 0.0;
  b_expl_temp.EIndS = 0.0;
  b_expl_temp.BIndS = 0.0;
  i = b_expl_temp.EAreaR->size[0];
  b_expl_temp.EAreaR->size[0] = strut_EAreaR->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.EAreaR, i);
  c_loop_ub = strut_EAreaR->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.EAreaR->data[i] = strut_EAreaR->data[i];
  }

  emxFree_int8_T(&strut_EAreaR);
  i = b_expl_temp.ECtoR->size[0];
  b_expl_temp.ECtoR->size[0] = strut_ECtoR->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.ECtoR, i);
  c_loop_ub = strut_ECtoR->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.ECtoR->data[i] = strut_ECtoR->data[i];
  }

  emxFree_int8_T(&strut_ECtoR);
  i = b_expl_temp.sEz->size[0];
  b_expl_temp.sEz->size[0] = strut_sEz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEz, i);
  c_loop_ub = strut_sEz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEz->data[i] = strut_sEz->data[i];
  }

  emxFree_int8_T(&strut_sEz);
  i = b_expl_temp.sEy->size[0];
  b_expl_temp.sEy->size[0] = strut_sEy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEy, i);
  c_loop_ub = strut_sEy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEy->data[i] = strut_sEy->data[i];
  }

  emxFree_int8_T(&strut_sEy);
  i = b_expl_temp.sEx->size[0];
  b_expl_temp.sEx->size[0] = strut_sEx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.sEx, i);
  c_loop_ub = strut_sEx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.sEx->data[i] = strut_sEx->data[i];
  }

  emxFree_int8_T(&strut_sEx);
  i = b_expl_temp.PEz->size[0];
  b_expl_temp.PEz->size[0] = strut_PEz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEz, i);
  c_loop_ub = strut_PEz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEz->data[i] = strut_PEz->data[i];
  }

  emxFree_int8_T(&strut_PEz);
  i = b_expl_temp.PEy->size[0];
  b_expl_temp.PEy->size[0] = strut_PEy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEy, i);
  c_loop_ub = strut_PEy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEy->data[i] = strut_PEy->data[i];
  }

  emxFree_int8_T(&strut_PEy);
  i = b_expl_temp.PEx->size[0];
  b_expl_temp.PEx->size[0] = strut_PEx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.PEx, i);
  c_loop_ub = strut_PEx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.PEx->data[i] = strut_PEx->data[i];
  }

  emxFree_int8_T(&strut_PEx);
  i = b_expl_temp.CtoR->size[0];
  b_expl_temp.CtoR->size[0] = strut_CtoR->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.CtoR, i);
  c_loop_ub = strut_CtoR->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.CtoR->data[i] = strut_CtoR->data[i];
  }

  emxFree_int8_T(&strut_CtoR);
  i = b_expl_temp.MCz->size[0];
  b_expl_temp.MCz->size[0] = strut_MCz->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCz, i);
  c_loop_ub = strut_MCz->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCz->data[i] = strut_MCz->data[i];
  }

  emxFree_int8_T(&strut_MCz);
  i = b_expl_temp.MCy->size[0];
  b_expl_temp.MCy->size[0] = strut_MCy->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCy, i);
  c_loop_ub = strut_MCy->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCy->data[i] = strut_MCy->data[i];
  }

  emxFree_int8_T(&strut_MCy);
  i = b_expl_temp.MCx->size[0];
  b_expl_temp.MCx->size[0] = strut_MCx->size[0];
  emxEnsureCapacity_real_T(b_expl_temp.MCx, i);
  c_loop_ub = strut_MCx->size[0];
  for (i = 0; i < c_loop_ub; i++) {
    b_expl_temp.MCx->data[i] = strut_MCx->data[i];
  }

  emxFree_int8_T(&strut_MCx);
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
  emxInit_real_T(&TtoC, 1);
  emxInit_real_T(&MCx, 1);
  emxInit_real_T(&MCy, 1);
  emxInit_real_T(&MCz, 1);
  emxInit_real_T(&BIndS, 1);
  emxInit_real_T(&EIndS, 1);
  emxInit_real_T(&BIndE, 1);
  emxInit_real_T(&EIndE, 1);
  for (b_i = 0; b_i < i; b_i++) {
    if (b_i + 1 != 1) {
      b_myfgetl(static_cast<double>(fileid));

      //  skip line
      b_myfgetl(static_cast<double>(fileid));

      //  skip line
    }

    cactusGeom_strut->data[b_i].NElem = NElem->data[0];
    processLine(static_cast<double>(fileid), TtoC);
    cactusGeom_strut->data[b_i].TtoC = TtoC->data[0];
    processLine(static_cast<double>(fileid), MCx);
    for (b_loop_ub = 0; b_loop_ub < y_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].MCx->data[b_loop_ub] = MCx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), MCy);
    for (b_loop_ub = 0; b_loop_ub < ab_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].MCy->data[b_loop_ub] = MCy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), MCz);
    for (b_loop_ub = 0; b_loop_ub < bb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].MCz->data[b_loop_ub] = MCz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), CtoR);
    for (b_loop_ub = 0; b_loop_ub < cb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].CtoR->data[b_loop_ub] = CtoR->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEx);
    for (b_loop_ub = 0; b_loop_ub < db_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].PEx->data[b_loop_ub] = PEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEy);
    for (b_loop_ub = 0; b_loop_ub < eb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].PEy->data[b_loop_ub] = PEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), PEz);
    for (b_loop_ub = 0; b_loop_ub < fb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].PEz->data[b_loop_ub] = PEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEx);
    for (b_loop_ub = 0; b_loop_ub < gb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].sEx->data[b_loop_ub] = sEx->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEy);
    for (b_loop_ub = 0; b_loop_ub < hb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].sEy->data[b_loop_ub] = sEy->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), sEz);
    for (b_loop_ub = 0; b_loop_ub < ib_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].sEz->data[b_loop_ub] = sEz->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), ECtoR);
    for (b_loop_ub = 0; b_loop_ub < jb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].ECtoR->data[b_loop_ub] = ECtoR->data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), EAreaR);
    for (b_loop_ub = 0; b_loop_ub < kb_loop_ub; b_loop_ub++) {
      cactusGeom_strut->data[b_i].EAreaR->data[b_loop_ub] = EAreaR->
        data[b_loop_ub];
    }

    processLine(static_cast<double>(fileid), BIndS);
    cactusGeom_strut->data[b_i].BIndS = BIndS->data[0];
    processLine(static_cast<double>(fileid), EIndS);
    cactusGeom_strut->data[b_i].EIndS = EIndS->data[0];
    processLine(static_cast<double>(fileid), BIndE);
    cactusGeom_strut->data[b_i].BIndE = BIndE->data[0];
    processLine(static_cast<double>(fileid), EIndE);
    cactusGeom_strut->data[b_i].EIndE = EIndE->data[0];
  }

  emxFree_real_T(&EIndE);
  emxFree_real_T(&BIndE);
  emxFree_real_T(&EIndS);
  emxFree_real_T(&BIndS);
  emxFree_real_T(&MCz);
  emxFree_real_T(&MCy);
  emxFree_real_T(&MCx);
  emxFree_real_T(&TtoC);
  emxFree_real_T(&EAreaR);
  emxFree_real_T(&ECtoR);
  emxFree_real_T(&sEz);
  emxFree_real_T(&sEy);
  emxFree_real_T(&sEx);
  emxFree_real_T(&PEz);
  emxFree_real_T(&PEy);
  emxFree_real_T(&PEx);
  emxFree_real_T(&CtoR);
  emxFree_real_T(&NElem);
  cfclose(static_cast<double>(fileid));
}

//
// File trailer for readCactusGeom.cpp
//
// [EOF]
//
