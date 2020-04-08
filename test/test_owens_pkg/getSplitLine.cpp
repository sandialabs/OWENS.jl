//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: getSplitLine.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "getSplitLine.h"
#include "myfgetl.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : double fid
//                emxArray_real_T *data
// Return Type  : void
//
void getSplitLine(double fid, emxArray_real_T *data)
{
  emxArray_char_T *line;
  emxArray_boolean_T *x;
  int i;
  int loop_ub;
  emxArray_int32_T *ii;
  int b_ii;
  int nx;
  int idx;
  boolean_T exitg1;
  emxArray_uint32_T *delimiter_idx;
  emxArray_char_T *b_line;
  int k;
  unsigned int u;
  unsigned int u1;
  creal_T dc;
  emxInit_char_T(&line, 2);
  emxInit_boolean_T(&x, 2);
  myfgetl(fid, line);

  //  Find where all of the delimiters are
  i = x->size[0] * x->size[1];
  x->size[0] = line->size[1];
  x->size[1] = line->size[0];
  emxEnsureCapacity_boolean_T(x, i);
  loop_ub = line->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_ii = line->size[1];
    for (idx = 0; idx < b_ii; idx++) {
      x->data[idx + x->size[0] * i] = (line->data[i + line->size[0] * idx] ==
        '\x09');
    }
  }

  emxInit_int32_T(&ii, 1);
  nx = x->size[0] * x->size[1];
  idx = 0;
  i = ii->size[0];
  ii->size[0] = nx;
  emxEnsureCapacity_int32_T(ii, i);
  b_ii = 0;
  exitg1 = false;
  while ((!exitg1) && (b_ii <= nx - 1)) {
    if (x->data[b_ii]) {
      idx++;
      ii->data[idx - 1] = b_ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        b_ii++;
      }
    } else {
      b_ii++;
    }
  }

  emxFree_boolean_T(&x);
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
    b_ii = 0;
  } else {
    b_ii = line->size[1];
  }

  emxInit_uint32_T(&delimiter_idx, 2);
  i = delimiter_idx->size[0] * delimiter_idx->size[1];
  delimiter_idx->size[0] = 1;
  delimiter_idx->size[1] = ii->size[0] + 2;
  emxEnsureCapacity_uint32_T(delimiter_idx, i);
  delimiter_idx->data[0] = 0U;
  loop_ub = ii->size[0];
  for (i = 0; i < loop_ub; i++) {
    delimiter_idx->data[i + 1] = static_cast<unsigned int>(ii->data[i]);
  }

  delimiter_idx->data[ii->size[0] + 1] = b_ii + 1U;

  //  % Reduce index, getting rid of duplicate delimiters if spaces
  //  if delim == ' '
  //      use_idx = true(1,length(delimiter_idx));
  //      for i = 1:length(delimiter_idx)
  //          for j = 1:length(delimiter_idx)-i
  //              if delimiter_idx(i)+j == delimiter_idx(i+j)
  //                  use_idx(i) = false;
  //              else
  //                  break
  //              end
  //          end
  //      end
  //      delimiter_idx = delimiter_idx(use_idx);
  //  end
  //  Extract the data from the beginning to the last delimiter
  i = data->size[0];
  data->size[0] = delimiter_idx->size[1] - 1;
  emxEnsureCapacity_real_T(data, i);
  loop_ub = delimiter_idx->size[1];
  emxFree_int32_T(&ii);
  for (i = 0; i <= loop_ub - 2; i++) {
    data->data[i] = 0.0;
  }

  i = delimiter_idx->size[1];
  emxInit_char_T(&b_line, 2);
  for (k = 0; k <= i - 2; k++) {
    u = delimiter_idx->data[k];
    u1 = delimiter_idx->data[k + 1];
    if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
      idx = 0;
      b_ii = 0;
    } else {
      idx = static_cast<int>(u);
      b_ii = static_cast<int>((static_cast<double>(u1) - 1.0));
    }

    nx = b_line->size[0] * b_line->size[1];
    b_line->size[0] = 1;
    loop_ub = b_ii - idx;
    b_line->size[1] = loop_ub;
    emxEnsureCapacity_char_T(b_line, nx);
    for (b_ii = 0; b_ii < loop_ub; b_ii++) {
      b_line->data[b_ii] = line->data[idx + b_ii];
    }

    dc = c_str2double(b_line);
    data->data[k] = dc.re;
  }

  emxFree_char_T(&b_line);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);
}

//
// File trailer for getSplitLine.cpp
//
// [EOF]
//
