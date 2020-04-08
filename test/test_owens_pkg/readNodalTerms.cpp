//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readNodalTerms.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "readNodalTerms.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// readNodalTerms reads concentrated nodal terms file
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [nodalTerms] = readNodalTerms(filename)
//
//    This function reads the nodal terms file and stores data in the nodal
//    terms object.
//
//       input:
//       filename      = string containing nodal terms filename
//
//       output:
//       nodalTerms    = object containing concentrated nodal data
// Arguments    : const emxArray_char_T *filename
// Return Type  : void
//
void readNodalTerms(const emxArray_char_T *filename)
{
  // initialize stiffness, load, and mass arrays to null
  b_cfopen(filename, "rb");

  // open nodal terms file
  //      index = 1;
  // store concentrated nodal term data in nodalTerms object
}

//
// File trailer for readNodalTerms.cpp
//
// [EOF]
//
