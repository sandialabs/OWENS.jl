//
// File: importCactusFile.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "importCactusFile.h"
#include "feof.h"
#include "fileManager.h"
#include "myfgetl.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_char_T *filename
//                emxArray_real_T *data
// Return Type  : void
//
void b_importCactusFile(const emxArray_char_T *filename, emxArray_real_T *data)
{
  boolean_T b_bool;
  int kstr;
  signed char fileid;
  int fid;
  int exitg1;
  int i;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  static double b_data[44044];
  double j;
  emxArray_char_T *line;
  emxArray_uint32_T *delimiter_idx;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  emxArray_char_T *b_line;
  double b_x;
  int loop_ub;
  int nx;
  int idx;
  int i1;
  boolean_T exitg2;
  unsigned int u;
  unsigned int u1;
  creal_T dc;
  b_bool = false;
  if (filename->size[1] == 3) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (filename->data[kstr] != b_cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = c_cfopen(filename, "rb");
    fid = fileid;
  }

  for (i = 0; i < 44044; i++) {
    b_data[i] = rtNaN;
  }

  // TODO: dont make this hard coded
  //  skip header
  b_myfgetl(static_cast<double>(fid));
  j = 0.0;
  emxInit_char_T(&line, 2);
  emxInit_uint32_T(&delimiter_idx, 2);
  emxInit_boolean_T(&x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_char_T(&b_line, 2);
  do {
    exitg1 = 0;
    b_x = b_feof(static_cast<double>(fid));
    if (!(b_x != 0.0)) {
      j++;
      myfgetl(static_cast<double>(fid), line);

      //  Find where all of the delimiters are
      i = x->size[0] * x->size[1];
      x->size[0] = line->size[1];
      x->size[1] = line->size[0];
      emxEnsureCapacity_boolean_T(x, i);
      loop_ub = line->size[0];
      for (i = 0; i < loop_ub; i++) {
        kstr = line->size[1];
        for (i1 = 0; i1 < kstr; i1++) {
          x->data[i1 + x->size[0] * i] = (line->data[i + line->size[0] * i1] ==
            ',');
        }
      }

      nx = x->size[0] * x->size[1];
      idx = 0;
      i = ii->size[0];
      ii->size[0] = nx;
      emxEnsureCapacity_int32_T(ii, i);
      kstr = 0;
      exitg2 = false;
      while ((!exitg2) && (kstr <= nx - 1)) {
        if (x->data[kstr]) {
          idx++;
          ii->data[idx - 1] = kstr + 1;
          if (idx >= nx) {
            exitg2 = true;
          } else {
            kstr++;
          }
        } else {
          kstr++;
        }
      }

      if (nx == 1) {
        if (idx == 0) {
          ii->size[0] = 0;
        }
      } else {
        i = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity_int32_T(ii, i);
      }

      if ((line->size[0] == 0) || (line->size[1] == 0)) {
        kstr = 0;
      } else {
        kstr = line->size[1];
      }

      i = delimiter_idx->size[0] * delimiter_idx->size[1];
      delimiter_idx->size[0] = 1;
      delimiter_idx->size[1] = ii->size[0] + 2;
      emxEnsureCapacity_uint32_T(delimiter_idx, i);
      delimiter_idx->data[0] = 0U;
      loop_ub = ii->size[0];
      for (i = 0; i < loop_ub; i++) {
        delimiter_idx->data[i + 1] = static_cast<unsigned int>(ii->data[i]);
      }

      delimiter_idx->data[ii->size[0] + 1] = kstr + 1U;

      //  Extract the data from the beginning to the last delimiter
      i = delimiter_idx->size[1];
      for (kstr = 0; kstr <= i - 2; kstr++) {
        u = delimiter_idx->data[kstr];
        u1 = delimiter_idx->data[kstr + 1];
        if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
          i1 = 0;
          nx = 0;
        } else {
          i1 = static_cast<int>(u);
          nx = static_cast<int>((static_cast<double>(u1) - 1.0));
        }

        idx = b_line->size[0] * b_line->size[1];
        b_line->size[0] = 1;
        loop_ub = nx - i1;
        b_line->size[1] = loop_ub;
        emxEnsureCapacity_char_T(b_line, idx);
        for (nx = 0; nx < loop_ub; nx++) {
          b_line->data[nx] = line->data[i1 + nx];
        }

        dc = c_str2double(b_line);
        b_data[(static_cast<int>(j) + 2002 * kstr) - 1] = dc.re;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_line);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);
  cfclose(static_cast<double>(fid));
  if (1.0 > j - 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>((j - 1.0));
  }

  i = data->size[0] * data->size[1];
  data->size[0] = loop_ub;
  data->size[1] = 22;
  emxEnsureCapacity_real_T(data, i);
  for (i = 0; i < 22; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      data->data[i1 + data->size[0] * i] = b_data[i1 + 2002 * i];
    }
  }

  // trim the excess off
}

//
// Arguments    : const emxArray_char_T *filename
//                double data_data[]
//                int data_size[2]
// Return Type  : void
//
void importCactusFile(const emxArray_char_T *filename, double data_data[], int
                      data_size[2])
{
  boolean_T b_bool;
  int kstr;
  signed char fileid;
  int fid;
  int exitg1;
  int i;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  double data[960];
  double j;
  emxArray_char_T *line;
  emxArray_uint32_T *delimiter_idx;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  emxArray_char_T *b_line;
  double b_x;
  int loop_ub;
  int nx;
  int idx;
  int i1;
  boolean_T exitg2;
  unsigned int u;
  unsigned int u1;
  creal_T dc;
  b_bool = false;
  if (filename->size[1] == 3) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (filename->data[kstr] != b_cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = c_cfopen(filename, "rb");
    fid = fileid;
  }

  for (i = 0; i < 960; i++) {
    data[i] = rtNaN;
  }

  // TODO: dont make this hard coded
  //  skip header
  j = 0.0;
  emxInit_char_T(&line, 2);
  emxInit_uint32_T(&delimiter_idx, 2);
  emxInit_boolean_T(&x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_char_T(&b_line, 2);
  do {
    exitg1 = 0;
    b_x = b_feof(static_cast<double>(fid));
    if (!(b_x != 0.0)) {
      j++;
      myfgetl(static_cast<double>(fid), line);

      //  Find where all of the delimiters are
      i = x->size[0] * x->size[1];
      x->size[0] = line->size[1];
      x->size[1] = line->size[0];
      emxEnsureCapacity_boolean_T(x, i);
      loop_ub = line->size[0];
      for (i = 0; i < loop_ub; i++) {
        kstr = line->size[1];
        for (i1 = 0; i1 < kstr; i1++) {
          x->data[i1 + x->size[0] * i] = (line->data[i + line->size[0] * i1] ==
            '\x09');
        }
      }

      nx = x->size[0] * x->size[1];
      idx = 0;
      i = ii->size[0];
      ii->size[0] = nx;
      emxEnsureCapacity_int32_T(ii, i);
      kstr = 0;
      exitg2 = false;
      while ((!exitg2) && (kstr <= nx - 1)) {
        if (x->data[kstr]) {
          idx++;
          ii->data[idx - 1] = kstr + 1;
          if (idx >= nx) {
            exitg2 = true;
          } else {
            kstr++;
          }
        } else {
          kstr++;
        }
      }

      if (nx == 1) {
        if (idx == 0) {
          ii->size[0] = 0;
        }
      } else {
        i = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity_int32_T(ii, i);
      }

      if ((line->size[0] == 0) || (line->size[1] == 0)) {
        kstr = 0;
      } else {
        kstr = line->size[1];
      }

      i = delimiter_idx->size[0] * delimiter_idx->size[1];
      delimiter_idx->size[0] = 1;
      delimiter_idx->size[1] = ii->size[0] + 2;
      emxEnsureCapacity_uint32_T(delimiter_idx, i);
      delimiter_idx->data[0] = 0U;
      loop_ub = ii->size[0];
      for (i = 0; i < loop_ub; i++) {
        delimiter_idx->data[i + 1] = static_cast<unsigned int>(ii->data[i]);
      }

      delimiter_idx->data[ii->size[0] + 1] = kstr + 1U;

      //  Extract the data from the beginning to the last delimiter
      i = delimiter_idx->size[1];
      for (kstr = 0; kstr <= i - 2; kstr++) {
        u = delimiter_idx->data[kstr];
        u1 = delimiter_idx->data[kstr + 1];
        if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
          i1 = 0;
          nx = 0;
        } else {
          i1 = static_cast<int>(u);
          nx = static_cast<int>((static_cast<double>(u1) - 1.0));
        }

        idx = b_line->size[0] * b_line->size[1];
        b_line->size[0] = 1;
        loop_ub = nx - i1;
        b_line->size[1] = loop_ub;
        emxEnsureCapacity_char_T(b_line, idx);
        for (nx = 0; nx < loop_ub; nx++) {
          b_line->data[nx] = line->data[i1 + nx];
        }

        dc = c_str2double(b_line);
        data[(static_cast<int>(j) + 60 * kstr) - 1] = dc.re;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_line);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);
  cfclose(static_cast<double>(fid));
  if (1.0 > j - 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>((j - 1.0));
  }

  data_size[0] = loop_ub;
  data_size[1] = 16;
  for (i = 0; i < 16; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      data_data[i1 + loop_ub * i] = data[i1 + 60 * i];
    }
  }

  // trim the excess off
}

//
// File trailer for importCactusFile.cpp
//
// [EOF]
//
