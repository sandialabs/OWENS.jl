//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: repmat.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "repmat.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const k_struct_T *a
//                int varargin_2
//                f_emxArray_struct_T *b
// Return Type  : void
//
void b_repmat(const k_struct_T *a, int varargin_2, f_emxArray_struct_T *b)
{
  int i;
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = varargin_2;
  emxEnsureCapacity_struct_T2(b, i);
  for (i = 0; i < varargin_2; i++) {
    emxCopyStruct_struct_T1(&b->data[i], a);
  }
}

//
// Arguments    : const j_struct_T a
//                int varargin_2
//                e_emxArray_struct_T *b
// Return Type  : void
//
void c_repmat(const j_struct_T a, int varargin_2, e_emxArray_struct_T *b)
{
  int i;
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = varargin_2;
  emxEnsureCapacity_struct_T3(b, i);
  for (i = 0; i < varargin_2; i++) {
    emxCopyStruct_struct_T2(&b->data[i], &a);
  }
}

//
// Arguments    : const i_struct_T *a
//                int varargin_2
//                d_emxArray_struct_T *b
// Return Type  : void
//
void repmat(const i_struct_T *a, int varargin_2, d_emxArray_struct_T *b)
{
  int i;
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = varargin_2;
  emxEnsureCapacity_struct_T1(b, i);
  for (i = 0; i < varargin_2; i++) {
    emxCopyStruct_struct_T(&b->data[i], a);
  }
}

//
// File trailer for repmat.cpp
//
// [EOF]
//
