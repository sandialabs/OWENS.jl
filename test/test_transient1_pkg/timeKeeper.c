/*
 * File: timeKeeper.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "timeKeeper.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <string.h>

/* Variable Definitions */
static boolean_T savedTime_not_empty;

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void savedTime_not_empty_init(void)
{
  savedTime_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void timeKeeper(void)
{
  struct timespec b_timespec;
  if (!savedTime_not_empty) {
    clock_gettime(CLOCK_MONOTONIC, &b_timespec);
    savedTime_not_empty = true;
  }
}

/*
 * File trailer for timeKeeper.c
 *
 * [EOF]
 */
