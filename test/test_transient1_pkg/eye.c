/*
 * File: eye.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "eye.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : double b_I[9]
 * Return Type  : void
 */
void eye(double b_I[9])
{
  memset(&b_I[0], 0, 9U * sizeof(double));
  b_I[0] = 1.0;
  b_I[4] = 1.0;
  b_I[8] = 1.0;
}

/*
 * File trailer for eye.c
 *
 * [EOF]
 */
