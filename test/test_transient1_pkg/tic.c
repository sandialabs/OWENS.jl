/*
 * File: tic.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "tic.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "timeKeeper.h"
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void tic(void)
{
  struct timespec b_timespec;
  clock_gettime(CLOCK_MONOTONIC, &b_timespec);
  timeKeeper();
}

/*
 * File trailer for tic.c
 *
 * [EOF]
 */
