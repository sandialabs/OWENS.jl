//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readMesh.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//
#ifndef READMESH_H
#define READMESH_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_readMesh(const emxArray_char_T *filename, emxArray_real_T
  *mesh_nodeNum, double *mesh_numEl, double *mesh_numNodes, emxArray_real_T
  *mesh_x, emxArray_real_T *mesh_y, emxArray_real_T *mesh_z, emxArray_real_T
  *mesh_elNum, emxArray_real_T *mesh_conn);
extern void readMesh(const emxArray_char_T *filename, emxArray_real_T
                     *mesh_nodeNum, double *mesh_numEl, double *mesh_numNodes,
                     emxArray_real_T *mesh_x, emxArray_real_T *mesh_y,
                     emxArray_real_T *mesh_z, emxArray_real_T *mesh_elNum,
                     emxArray_real_T *mesh_conn);

#endif

//
// File trailer for readMesh.h
//
// [EOF]
//
