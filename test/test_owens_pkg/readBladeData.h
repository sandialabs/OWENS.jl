//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readBladeData.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//
#ifndef READBLADEDATA_H
#define READBLADEDATA_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void readBladeData(const emxArray_char_T *filename, double
  *bladeData_numBlades, double bladeData_bladeNum_data[], int
  bladeData_bladeNum_size[1], double bladeData_h_data[], int bladeData_h_size[1],
  double bladeData_nodeNum_data[], int bladeData_nodeNum_size[1], double
  bladeData_elementNum_data[], int bladeData_elementNum_size[1], double
  bladeData_remaining_data[], int bladeData_remaining_size[2]);

#endif

//
// File trailer for readBladeData.h
//
// [EOF]
//
