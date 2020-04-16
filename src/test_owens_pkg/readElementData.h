//
// File: readElementData.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//
#ifndef READELEMENTDATA_H
#define READELEMENTDATA_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_readElementData(double numElements, const emxArray_char_T *elfile,
  const emxArray_char_T *ortfile, const double bladeData_struct_nodeNum_data[],
  const int bladeData_struct_nodeNum_size[1], const double
  c_bladeData_struct_elementNum_d[], const int c_bladeData_struct_elementNum_s[1],
  const double bladeData_struct_remaining_data[], const int
  bladeData_struct_remaining_size[2], emxArray_struct_T *el_props,
  emxArray_real_T *el_elLen, emxArray_real_T *el_psi, emxArray_real_T *el_theta,
  emxArray_real_T *el_roll, emxArray_boolean_T *el_rotationalEffects);
extern void readElementData(double numElements, const emxArray_char_T *elfile,
  const emxArray_char_T *ortfile, const double bladeData_struct_nodeNum_data[],
  const int bladeData_struct_nodeNum_size[1], const double
  c_bladeData_struct_elementNum_d[], const int c_bladeData_struct_elementNum_s[1],
  const double bladeData_struct_remaining_data[], const int
  bladeData_struct_remaining_size[2], emxArray_struct_T *el_props,
  emxArray_real_T *el_elLen, emxArray_real_T *el_psi, emxArray_real_T *el_theta,
  emxArray_real_T *el_roll, emxArray_boolean_T *el_rotationalEffects);

#endif

//
// File trailer for readElementData.h
//
// [EOF]
//
