//
// File: mtimes.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// Arguments    : const double A[144]
//                const double B_data[]
//                double C[12]
// Return Type  : void
//
void mtimes(const double A[144], const double B_data[], double C[12])
{
  int k;
  int aoffset;
  int i;
  std::memset(&C[0], 0, 12U * sizeof(double));
  for (k = 0; k < 12; k++) {
    aoffset = k * 12;
    for (i = 0; i < 12; i++) {
      C[i] += B_data[k] * A[aoffset + i];
    }
  }
}

//
// File trailer for mtimes.cpp
//
// [EOF]
//
