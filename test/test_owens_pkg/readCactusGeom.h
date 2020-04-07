//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readCactusGeom.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//
#ifndef READCACTUSGEOM_H
#define READCACTUSGEOM_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void readCactusGeom(const emxArray_char_T *geomfn, int *cactusGeom_NBlade,
  int *cactusGeom_NStrut, emxArray_real_T *cactusGeom_RotN, emxArray_real_T
  *cactusGeom_RotP, emxArray_real_T *cactusGeom_RefAR, emxArray_real_T
  *cactusGeom_RefR, d_emxArray_struct_T *cactusGeom_blade, f_emxArray_struct_T
  *cactusGeom_strut);

#endif

//
// File trailer for readCactusGeom.h
//
// [EOF]
//
