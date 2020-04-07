//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: assemblyMatrixOnly.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//
#ifndef ASSEMBLYMATRIXONLY_H
#define ASSEMBLYMATRIXONLY_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void assemblyMatrixOnly(const double Ke[144], const double conn[2],
  emxArray_real_T *Kg);
extern void b_assemblyMatrixOnly(const double Ke_data[], const int Ke_size[2],
  const double conn[2], emxArray_real_T *Kg);

#endif

//
// File trailer for assemblyMatrixOnly.h
//
// [EOF]
//
