/*
 * File: test_transient1_initialize.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "test_transient1_initialize.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_data.h"
#include "timeKeeper.h"
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void test_transient1_initialize(void)
{
  rt_InitInfAndNaN();
  savedTime_not_empty_init();
  filedata_init();
  isInitialized_test_transient1 = true;
}

/*
 * File trailer for test_transient1_initialize.c
 *
 * [EOF]
 */
