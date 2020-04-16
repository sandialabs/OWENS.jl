//
// File: linearAnalysisModal.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//
#ifndef LINEARANALYSISMODAL_H
#define LINEARANALYSISMODAL_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern void b_linearAnalysisModal(const char model_analysisType[2], boolean_T
  model_aeroElasticOn, double model_guessFreq, double model_RayleighAlpha,
  double model_RayleighBeta, double model_BC_numpBC, const emxArray_real_T
  *model_BC_pBC, const emxArray_real_T *model_BC_map, const emxArray_real_T
  *model_joint, const char model_outFilename_data[], const int
  model_outFilename_size[2], const emxArray_real_T *model_jointTransform, const
  emxArray_real_T *model_reducedDOFList, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const i_struct_T el, const
  emxArray_real_T *displ, const c_emxArray_struct_T *elStorage, emxArray_real_T *
  freq, emxArray_real_T *damp, emxArray_real_T *phase1, emxArray_real_T *phase2,
  emxArray_real_T *imagCompSign);
extern void linearAnalysisModal(double model_RayleighAlpha, double
  model_RayleighBeta, double model_BC_numpBC, const emxArray_real_T
  *model_BC_pBC, const emxArray_real_T *model_BC_map, const emxArray_real_T
  *model_joint, const char model_outFilename_data[], const int
  model_outFilename_size[2], const emxArray_real_T *model_jointTransform, const
  emxArray_real_T *model_reducedDOFList, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const i_struct_T el, const
  emxArray_real_T *displ, const c_emxArray_struct_T *elStorage, emxArray_real_T *
  freq, emxArray_real_T *damp, emxArray_real_T *phase1, emxArray_real_T *phase2,
  emxArray_real_T *imagCompSign);

#endif

//
// File trailer for linearAnalysisModal.h
//
// [EOF]
//
