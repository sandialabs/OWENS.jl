//
// File: owens.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//
#ifndef OWENS_H
#define OWENS_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_owens(emxArray_real_T *freq, emxArray_real_T *damp);
extern void c_owens(double freq_data[], int freq_size[2], double damp_data[],
                    int damp_size[2]);
extern void owens(const double varargin_7[2], double *freq, double *damp);

#endif

//
// File trailer for owens.h
//
// [EOF]
//
