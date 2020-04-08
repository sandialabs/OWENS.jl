//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: strcmp.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "strcmp.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const char a[3]
// Return Type  : boolean_T
//
boolean_T b_strcmp(const char a[3])
{
  int ret;
  static const char b[3] = { 'R', 'O', 'M' };

  ret = memcmp(&a[0], &b[0], 3);
  return ret == 0;
}

//
// Arguments    : const char a[3]
// Return Type  : boolean_T
//
boolean_T c_strcmp(const char a[3])
{
  int ret;
  static const char b[3] = { 'T', 'N', 'B' };

  ret = memcmp(&a[0], &b[0], 3);
  return ret == 0;
}

//
// Arguments    : const char a[3]
// Return Type  : boolean_T
//
boolean_T d_strcmp(const char a[3])
{
  int ret;
  static const char b[3] = { 'R', 'M', '0' };

  ret = memcmp(&a[0], &b[0], 3);
  return ret == 0;
}

//
// Arguments    : const char a[2]
// Return Type  : boolean_T
//
boolean_T e_strcmp(const char a[2])
{
  int ret;
  static const char b[2] = { 'N', 'R' };

  ret = memcmp(&a[0], &b[0], 2);
  return ret == 0;
}

//
// Arguments    : const char a[2]
// Return Type  : boolean_T
//
boolean_T f_strcmp(const char a[2])
{
  int ret;
  static const char b[2] = { 'D', 'I' };

  ret = memcmp(&a[0], &b[0], 2);
  return ret == 0;
}

//
// Arguments    : const char a_data[]
//                const int a_size[2]
// Return Type  : boolean_T
//
boolean_T g_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char b_cv[2] = { 'T', 'D' };

  b_bool = false;
  if (a_size[1] == 2) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 2) {
        if (a_data[kstr] != b_cv[kstr]) {
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

  return b_bool;
}

//
// Arguments    : const char a_data[]
//                const int a_size[2]
// Return Type  : boolean_T
//
boolean_T h_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char b_cv[3] = { 'T', 'N', 'B' };

  b_bool = false;
  if (a_size[1] == 3) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (a_data[kstr] != b_cv[kstr]) {
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

  return b_bool;
}

//
// Arguments    : const char a_data[]
//                const int a_size[2]
// Return Type  : boolean_T
//
boolean_T i_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  b_bool = false;
  if (a_size[1] == 1) {
    b_bool = ((!(a_data[0] != 'M')) || b_bool);
  }

  return b_bool;
}

//
// Arguments    : const char a_data[]
//                const int a_size[2]
// Return Type  : boolean_T
//
boolean_T j_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char b_cv[3] = { 'R', 'M', '0' };

  b_bool = false;
  if (a_size[1] == 3) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (a_data[kstr] != b_cv[kstr]) {
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

  return b_bool;
}

//
// Arguments    : const char a_data[]
//                const int a_size[2]
// Return Type  : boolean_T
//
boolean_T k_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  b_bool = false;
  if (a_size[1] == 1) {
    b_bool = ((!(a_data[0] != 'S')) || b_bool);
  }

  return b_bool;
}

//
// File trailer for strcmp.cpp
//
// [EOF]
//
