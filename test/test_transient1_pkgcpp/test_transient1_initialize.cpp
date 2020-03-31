//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: test_transient1_initialize.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "test_transient1_initialize.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_data.h"
#include "timeKeeper.h"
#include <string.h>

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void test_transient1_initialize()
{
  rt_InitInfAndNaN();
  savedTime_not_empty_init();
  filedata_init();
  isInitialized_test_transient1 = true;
}

//
// File trailer for test_transient1_initialize.cpp
//
// [EOF]
//
