//
// File: readInitCond.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "readInitCond.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// readInitCond reads initial conditions
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [initCond] = readInitCond(filename)
//
//    This function reads initial conditions from file
//
//    input:
//    filename      = string containing file name for initial conditions file
//
//    output:
//    initCond      = array containing initial conditions
// Arguments    : const emxArray_char_T *filename
// Return Type  : void
//
void readInitCond(const emxArray_char_T *filename)
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  // initialize intial condition to null
  b_bool = false;
  if ((filename->size[0] == 1) && (filename->size[1] == 3)) {
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

  if (!b_bool) {
    b_cfopen(filename, "rb");
  }

  // open initial  conditions file
}

//
// File trailer for readInitCond.cpp
//
// [EOF]
//
