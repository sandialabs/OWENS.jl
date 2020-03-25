/*
 * File: test_transient1_rtwutil.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "test_transient1_rtwutil.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <math.h>
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : double u
 * Return Type  : double
 */
double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * File trailer for test_transient1_rtwutil.c
 *
 * [EOF]
 */
