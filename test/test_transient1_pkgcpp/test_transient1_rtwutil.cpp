//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: test_transient1_rtwutil.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "test_transient1_rtwutil.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : double u
// Return Type  : double
//
double rt_roundd_snf(double u)
{
  double y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

//
// File trailer for test_transient1_rtwutil.cpp
//
// [EOF]
//
