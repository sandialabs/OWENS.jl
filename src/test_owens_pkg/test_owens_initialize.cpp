//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: test_owens_initialize.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "test_owens_initialize.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "timeKeeper.h"
#include <string.h>

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void test_owens_initialize()
{
  rt_InitInfAndNaN();
  savedTime_not_empty_init();
  filedata_init();
  isInitialized_test_owens = true;
}

//
// File trailer for test_owens_initialize.cpp
//
// [EOF]
//
