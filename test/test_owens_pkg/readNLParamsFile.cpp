//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readNLParamsFile.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "readNLParamsFile.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// readNLParamsFile   reads file for nonlinear iteration/load stepping
//  **********************************************************************
//  *                   Part of the SNL OWENS toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [nlParams] = readNLParamsFile(inputfile)
//
//    This function reads a nonlinear parameter for nonlinear iteration and
//    load stepping.
//
//       input:
//       inputfile    = string containing filename of main input file
//
//       output:
//       nlParams     = struct containing nonlinear parameter information.
// Arguments    : char nlParams_iterationType[2]
//                boolean_T *c_nlParams_adaptiveLoadStepping
//                double *nlParams_tolerance
//                double *nlParams_maxIterations
//                double *nlParams_maxNumLoadSteps
//                double *nlParams_minLoadStepDelta
//                double *nlParams_minLoadStep
//                double *nlParams_prescribedLoadStep
// Return Type  : void
//
void readNLParamsFile(char nlParams_iterationType[2], boolean_T
                      *c_nlParams_adaptiveLoadStepping, double
                      *nlParams_tolerance, double *nlParams_maxIterations,
                      double *nlParams_maxNumLoadSteps, double
                      *nlParams_minLoadStepDelta, double *nlParams_minLoadStep,
                      double *nlParams_prescribedLoadStep)
{
  // get main file prefix
  // create .nl file name string
  cfopen("./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.nl",
         "rb");

  // attempt to open file
  // assign parameters to nlParams file
  nlParams_iterationType[0] = 'N';
  nlParams_iterationType[1] = 'R';
  *c_nlParams_adaptiveLoadStepping = true;
  *nlParams_tolerance = 1.0E-6;
  *nlParams_maxIterations = 50.0;
  *nlParams_maxNumLoadSteps = 20.0;
  *nlParams_minLoadStepDelta = 0.05;
  *nlParams_minLoadStep = 0.05;
  *nlParams_prescribedLoadStep = 1.0;

  //  not used but must be declared
}

//
// File trailer for readNLParamsFile.cpp
//
// [EOF]
//
