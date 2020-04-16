//
// File: readNLParamsFile.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
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
  signed char fileid;

  // get main file prefix
  // create .nl file name string
  fileid = cfopen("./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.nl",
                  "rb");

  // attempt to open file
  if (fileid == -1) {
    // default parameters if file cannot be opened (doesnt exist)
    // uses adaptive load stepping by default
    nlParams_iterationType[0] = 'N';
    nlParams_iterationType[1] = 'R';
  } else {
    // if file can be opened read it
    //          iterationType = fscanf(fid,'%c',2); fgetl(fid); %read iteration type 'NR' = Newton Raphson, 'DI' = Direct Iteration 
    //          tolerance = fscanf(fid,'%f',1); fgetl(fid);
    //          maxIterations    = fscanf(fid,'%i',1); fgetl(fid); %read in maximum iterations allowed per time step 
    //          temp  = fscanf(fid,'%i',1);
    //          if(temp == 0) % if temp = 0 adaptive load stepping, read in params 
    //              fgetl(fid);
    //              adaptiveLoadSteppingFlag = true;
    //              maxNumLoadSteps     = fscanf(fid,'%i',1); fgetl(fid);
    //              minLoadStep      = fscanf(fid,'%f',1); fgetl(fid);
    //              minLoadStepDelta= fscanf(fid,'%f',1); fgetl(fid);
    //              fclose(fid); %close file
    //          elseif(temp>0) %if temp > 0, prescribed load stepping profile, read in params 
    //              adaptiveLoadSteppingFlag = false;
    //              prescribedLoadStep = zeros(temp,1);
    //              for i=1:temp
    //                  prescribedLoadStep(i) = fscanf(fid,'%f',1);
    //              end
    //              fgetl(fid);
    //              maxNumLoadSteps = 1e6;
    //              fclose(fid); %close file
    //          else
    //              fclose(fid);
    //              error('Load stepping parameter in .nl file not recognized. Exiting.'); 
    //          end
  }

  // assign parameters to nlParams file
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
