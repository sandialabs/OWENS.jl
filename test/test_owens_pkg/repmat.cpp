//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: repmat.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 03-Apr-2020 15:56:19
//

// Include Files
#include "repmat.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const m_struct_T *a
//                int varargin_2
//                f_emxArray_struct_T *b
// Return Type  : void
//
void b_repmat(const m_struct_T *a, int varargin_2, f_emxArray_struct_T *b)
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
// Arguments    : const l_struct_T a
//                int varargin_2
//                e_emxArray_struct_T *b
// Return Type  : void
//
void c_repmat(const l_struct_T a, int varargin_2, e_emxArray_struct_T *b)
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
// Arguments    : const k_struct_T *a
//                int varargin_2
//                d_emxArray_struct_T *b
// Return Type  : void
//
void repmat(const k_struct_T *a, int varargin_2, d_emxArray_struct_T *b)
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
