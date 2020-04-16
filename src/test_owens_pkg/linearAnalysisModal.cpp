//
// File: linearAnalysisModal.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "linearAnalysisModal.h"
#include "ConcMassAssociatedWithElement.h"
#include "applyBCModal.h"
#include "assemblyMatrixOnly.h"
#include "calculateTimoshenkoElementNL.h"
#include "eigSolve.h"
#include "extractFreqDamp.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "staticAnalysis.h"
#include "strcmp.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "tic.h"
#include "writeOutput.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// linearAnalysisModal performs modal analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freq,damp,phase1,phase2,imagCompSign] = linearAnalysisModal(model,mesh,
//                                                  el,displ,Omega,elStorage)
//
//    This function performs a modal analysis of a structural dynamics
//    system and returns freq, damping, and mode shapes
//
//    input:
//    model          = object containing model information
//    mesh           = object containing mesh information
//    el             = object containing element information
//    displ          = displacement vector for use in pre-stressed analysis
//    Omega          = rotor speed (Hz)
//    elStorage      = previously calculated element system matrices
//
//
//    output:
//    freq         = modal frequencies
//    damp         = modal damping ratios
//    phase1       = in phase mode shapes (real part of mode shape)
//    phase2       = out of phase mode shapes (imaginary part of mode shape)
//    imagCompSign = sign of imaginary component of eigenvalues
// Arguments    : const char model_analysisType[2]
//                boolean_T model_aeroElasticOn
//                double model_guessFreq
//                double model_RayleighAlpha
//                double model_RayleighBeta
//                double model_BC_numpBC
//                const emxArray_real_T *model_BC_pBC
//                const emxArray_real_T *model_BC_map
//                const emxArray_real_T *model_joint
//                const char model_outFilename_data[]
//                const int model_outFilename_size[2]
//                const emxArray_real_T *model_jointTransform
//                const emxArray_real_T *model_reducedDOFList
//                double mesh_numEl
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                const i_struct_T el
//                const emxArray_real_T *displ
//                const c_emxArray_struct_T *elStorage
//                emxArray_real_T *freq
//                emxArray_real_T *damp
//                emxArray_real_T *phase1
//                emxArray_real_T *phase2
//                emxArray_real_T *imagCompSign
// Return Type  : void
//
void b_linearAnalysisModal(const char model_analysisType[2], boolean_T
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
  emxArray_real_T *imagCompSign)
{
  emxArray_real_T *Kg;
  double totalNumDOF;
  int i;
  int i1;
  int loop_ub_tmp;
  emxArray_real_T *Mg;
  emxArray_real_T *Cg;
  int b_i;
  double elInput_freq;
  emxArray_creal_T *eigVec;
  emxArray_creal_T *eigVal;
  double eldisp_data[12];
  emxArray_real_T *b_r;
  emxArray_real_T *b_r1;
  int b_index;
  double elInput_xloc[2];
  int j;
  double b_mesh_conn[2];
  double elInput_concMass[8];
  double elInput_sectionProps_ac[2];
  double elInput_sectionProps_twist[2];
  double elInput_sectionProps_rhoA[2];
  double elInput_sectionProps_EA[2];
  double elInput_sectionProps_zcm[2];
  double elInput_CN2H[9];
  double elInput_sectionProps_ycm[2];
  double elInput_sectionProps_a[2];
  double elInput_sectionProps_b[2];
  double elInput_sectionProps_a0[2];
  double elx_data[2];
  double elInput_y_data[2];
  double elInput_z_data[2];
  p_struct_T expl_temp;
  double elOutput_Ke[144];
  int elOutput_Me_size[2];
  double elOutput_Me_data[144];
  int elOutput_Ce_size[2];
  double elOutput_Ce_data[144];
  emxArray_real_T *b_r2;
  emxArray_real_T *r3;
  emxArray_creal_T *r4;
  emxArray_creal_T *b_eigVec;
  emxArray_real_T *b_freq;
  emxArray_char_T b_model_outFilename_data;
  signed char fileid;
  emxArray_real_T *b_damp;
  emxArray_real_T *b_imagCompSign;
  emxInit_real_T(&Kg, 2);
  tic();

  // extract mesh information from mesh object
  // extract element order from model
  // extract boundary  condition object from model
  // do initialization
  totalNumDOF = static_cast<double>(mesh_x->size[0]) * 6.0;
  i = Kg->size[0] * Kg->size[1];
  i1 = static_cast<int>(totalNumDOF);
  Kg->size[0] = i1;
  Kg->size[1] = i1;
  emxEnsureCapacity_real_T(Kg, i);
  loop_ub_tmp = i1 * i1;
  for (i = 0; i < loop_ub_tmp; i++) {
    Kg->data[i] = 0.0;
  }

  emxInit_real_T(&Mg, 2);
  i = Mg->size[0] * Mg->size[1];
  Mg->size[0] = i1;
  Mg->size[1] = i1;
  emxEnsureCapacity_real_T(Mg, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    Mg->data[i] = 0.0;
  }

  emxInit_real_T(&Cg, 2);
  i = Cg->size[0] * Cg->size[1];
  Cg->size[0] = i1;
  Cg->size[1] = i1;
  emxEnsureCapacity_real_T(Cg, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    Cg->data[i] = 0.0;
  }

  // extract concentrated nodal terms from model
  i = static_cast<int>(mesh_numEl);
  if (0 <= i - 1) {
    if (model_aeroElasticOn) {
      elInput_freq = model_guessFreq * 2.0 * 3.1415926535897931;

      // set guess frequency if aeroelastic analysis
    } else {
      elInput_freq = 0.0;

      // Declare variable on all execution paths
    }
  }

  for (b_i = 0; b_i < i; b_i++) {
    // element loop
    std::memset(&eldisp_data[0], 0, 12U * sizeof(double));

    // Calculate Ke and Fe for element i
    b_index = 0;

    // define element input object flags and element properties from el object
    elInput_xloc[0] = 0.0;
    elInput_xloc[1] = el.elLen->data[b_i];

    // initialize element coordinate list
    // retrieve concentrated nodal terms associated with element
    for (j = 0; j < 2; j++) {
      elInput_sectionProps_ac[j] = el.props->data[b_i].ac[j];
      elInput_sectionProps_twist[j] = el.props->data[b_i].twist[j];
      elInput_sectionProps_rhoA[j] = el.props->data[b_i].rhoA[j];
      elInput_sectionProps_EA[j] = el.props->data[b_i].EA[j];
      elInput_sectionProps_zcm[j] = el.props->data[b_i].zcm[j];
      elInput_sectionProps_ycm[j] = el.props->data[b_i].ycm[j];
      elInput_sectionProps_a[j] = el.props->data[b_i].a[j];
      elInput_sectionProps_b[j] = el.props->data[b_i].b[j];
      elInput_sectionProps_a0[j] = el.props->data[b_i].a0[j];

      // define element coordinates and displacements associated with element
      loop_ub_tmp = static_cast<int>(mesh_conn->data[b_i + mesh_conn->size[0] *
        j]) - 1;
      elx_data[j] = mesh_x->data[loop_ub_tmp];
      elInput_y_data[j] = mesh_y->data[loop_ub_tmp];
      elInput_z_data[j] = mesh_z->data[loop_ub_tmp];
      for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
        eldisp_data[b_index] = displ->data[static_cast<int>(((mesh_conn->
          data[b_i + mesh_conn->size[0] * j] - 1.0) * 6.0 + (static_cast<double>
          (loop_ub_tmp) + 1.0))) - 1];
        b_index++;
      }

      b_mesh_conn[j] = mesh_conn->data[b_i + mesh_conn->size[0] * j];
    }

    ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);

    // assign concentrated nodal terms and coordinates to element input
    // object
    // activate or deactivate rotational effects for element
    // set aeroelastic flag
    std::memset(&elInput_CN2H[0], 0, 9U * sizeof(double));
    elInput_CN2H[0] = 1.0;
    elInput_CN2H[4] = 1.0;
    elInput_CN2H[8] = 1.0;
    d_calculateTimoshenkoElementNL(elInput_xloc, elInput_sectionProps_ac,
      elInput_sectionProps_twist, elInput_sectionProps_rhoA,
      elInput_sectionProps_EA, elInput_sectionProps_zcm,
      elInput_sectionProps_ycm, elInput_sectionProps_a, elInput_sectionProps_b,
      elInput_sectionProps_a0, el.psi->data[b_i], el.theta->data[b_i],
      el.roll->data[b_i], elInput_concMass, eldisp_data, elx_data,
      elInput_y_data, elInput_z_data, 0.11833333333333333, model_aeroElasticOn,
      elInput_CN2H, model_RayleighAlpha, model_RayleighBeta, elInput_freq,
      &elStorage->data[b_i], &expl_temp);
    std::memcpy(&elOutput_Ke[0], &expl_temp.Ke[0], 144U * sizeof(double));
    elOutput_Me_size[0] = expl_temp.Me.size[0];
    elOutput_Me_size[1] = expl_temp.Me.size[1];
    loop_ub_tmp = expl_temp.Me.size[0] * expl_temp.Me.size[1];
    if (0 <= loop_ub_tmp - 1) {
      std::memcpy(&elOutput_Me_data[0], &expl_temp.Me.data[0], loop_ub_tmp *
                  sizeof(double));
    }

    elOutput_Ce_size[0] = expl_temp.Ce.size[0];
    elOutput_Ce_size[1] = expl_temp.Ce.size[1];
    loop_ub_tmp = expl_temp.Ce.size[0] * expl_temp.Ce.size[1];
    if (0 <= loop_ub_tmp - 1) {
      std::memcpy(&elOutput_Ce_data[0], &expl_temp.Ce.data[0], loop_ub_tmp *
                  sizeof(double));
    }

    // do element calculation
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    assemblyMatrixOnly(elOutput_Ke, b_mesh_conn, Kg);

    // assemble element into global stiffness matrix
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    b_assemblyMatrixOnly(elOutput_Me_data, elOutput_Me_size, b_mesh_conn, Mg);

    // assemble element into global mass matrix
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    b_assemblyMatrixOnly(elOutput_Ce_data, elOutput_Ce_size, b_mesh_conn, Cg);

    // assemble element into global damping matrix
  }

  emxInit_creal_T(&eigVec, 2);
  emxInit_creal_T(&eigVal, 2);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_r1, 2);

  // apply general 6x6  mass, damping, and stiffness matrices to nodes
  // ----------------------------------------------------------------------
  // APPLY CONSTRAINT
  b_applyConstraints(Kg, model_jointTransform);

  // modify global matrices for joint constraints using joint transform
  b_applyConstraints(Mg, model_jointTransform);
  b_applyConstraints(Cg, model_jointTransform);

  // APPLY BOUNDARY CONDITIONS
  // apply boundary conditions to global matrices
  applyBCModal(Mg, model_BC_numpBC, model_BC_map, b_r);
  applyBCModal(Cg, model_BC_numpBC, model_BC_map, b_r1);
  applyBCModal(Kg, model_BC_numpBC, model_BC_map, Mg);
  eigSolve(b_r, b_r1, Mg, eigVec, eigVal);

  // ,... %eigensolve of global system
  // model.numModesToExtract,solveFlag);
  // save eigVectors eigVec %save eigenvector for later use (if needed) TODO: Doesn't appear to be used 
  // extract frequency, damping, mode shapes from eigenvalues and vectors
  i = phase1->size[0] * phase1->size[1] * phase1->size[2];
  phase1->size[0] = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0));
  phase1->size[1] = 6;
  phase1->size[2] = eigVal->size[1];
  emxEnsureCapacity_real_T(phase1, i);
  loop_ub_tmp = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0)) *
    6 * eigVal->size[1];
  emxFree_real_T(&b_r1);
  emxFree_real_T(&b_r);
  emxFree_real_T(&Cg);
  emxFree_real_T(&Mg);
  emxFree_real_T(&Kg);
  for (i = 0; i < loop_ub_tmp; i++) {
    phase1->data[i] = 0.0;
  }

  i = phase2->size[0] * phase2->size[1] * phase2->size[2];
  phase2->size[0] = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0));
  phase2->size[1] = 6;
  phase2->size[2] = eigVal->size[1];
  emxEnsureCapacity_real_T(phase2, i);
  loop_ub_tmp = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0)) *
    6 * eigVal->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    phase2->data[i] = 0.0;
  }

  i = eigVal->size[1];
  i1 = freq->size[0] * freq->size[1];
  freq->size[0] = 1;
  freq->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(freq, i1);
  i1 = damp->size[0] * damp->size[1];
  damp->size[0] = 1;
  damp->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(damp, i1);
  i1 = imagCompSign->size[0] * imagCompSign->size[1];
  imagCompSign->size[0] = 1;
  imagCompSign->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(imagCompSign, i1);
  emxInit_real_T(&b_r2, 2);
  emxInit_real_T(&r3, 2);
  emxInit_creal_T(&r4, 2);
  emxInit_creal_T(&b_eigVec, 1);
  for (b_i = 0; b_i < i; b_i++) {
    loop_ub_tmp = eigVec->size[0];
    i1 = b_eigVec->size[0];
    b_eigVec->size[0] = eigVec->size[0];
    emxEnsureCapacity_creal_T(b_eigVec, i1);
    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      b_eigVec->data[i1] = eigVec->data[i1 + eigVec->size[0] * b_i];
    }

    extractFreqDamp(eigVal->data[b_i + eigVal->size[0] * b_i], b_eigVec,
                    model_jointTransform, model_reducedDOFList, model_BC_numpBC,
                    model_BC_pBC, &freq->data[b_i], &damp->data[b_i], b_r2, r3,
                    r4);
    loop_ub_tmp = b_r2->size[0];
    b_index = r3->size[0];
    for (i1 = 0; i1 < 6; i1++) {
      for (j = 0; j < loop_ub_tmp; j++) {
        phase1->data[(j + phase1->size[0] * i1) + phase1->size[0] * 6 * b_i] =
          b_r2->data[j + b_r2->size[0] * i1];
      }

      for (j = 0; j < b_index; j++) {
        phase2->data[(j + phase2->size[0] * i1) + phase2->size[0] * 6 * b_i] =
          r3->data[j + r3->size[0] * i1];
      }
    }

    totalNumDOF = eigVal->data[b_i + eigVal->size[0] * b_i].im;
    if (eigVal->data[b_i + eigVal->size[0] * b_i].im < 0.0) {
      totalNumDOF = -1.0;
    } else if (eigVal->data[b_i + eigVal->size[0] * b_i].im > 0.0) {
      totalNumDOF = 1.0;
    } else {
      if (eigVal->data[b_i + eigVal->size[0] * b_i].im == 0.0) {
        totalNumDOF = 0.0;
      }
    }

    imagCompSign->data[b_i] = totalNumDOF;
  }

  emxFree_creal_T(&b_eigVec);
  emxFree_creal_T(&r4);
  emxFree_real_T(&r3);
  emxFree_real_T(&b_r2);
  emxFree_creal_T(&eigVal);
  emxFree_creal_T(&eigVec);

  //  %write output
  //  t_modal = toc;
  //  disp('Elapsed time for modal analysis(s):');
  //  disp(t_modal);
  if (!n_strcmp(model_analysisType)) {
    emxInit_real_T(&b_freq, 2);
    b_model_outFilename_data.data = const_cast<char *>(&model_outFilename_data[0]);
    b_model_outFilename_data.size = const_cast<int *>(&model_outFilename_size[0]);
    b_model_outFilename_data.allocatedSize = -1;
    b_model_outFilename_data.numDimensions = 2;
    b_model_outFilename_data.canFreeData = false;
    fileid = c_cfopen(&b_model_outFilename_data, "wb");
    i = b_freq->size[0] * b_freq->size[1];
    b_freq->size[0] = 1;
    b_freq->size[1] = freq->size[1];
    emxEnsureCapacity_real_T(b_freq, i);
    loop_ub_tmp = freq->size[0] * freq->size[1] - 1;
    for (i = 0; i <= loop_ub_tmp; i++) {
      b_freq->data[i] = freq->data[i];
    }

    emxInit_real_T(&b_damp, 2);
    i = b_damp->size[0] * b_damp->size[1];
    b_damp->size[0] = 1;
    b_damp->size[1] = damp->size[1];
    emxEnsureCapacity_real_T(b_damp, i);
    loop_ub_tmp = damp->size[0] * damp->size[1] - 1;
    for (i = 0; i <= loop_ub_tmp; i++) {
      b_damp->data[i] = damp->data[i];
    }

    emxInit_real_T(&b_imagCompSign, 2);
    i = b_imagCompSign->size[0] * b_imagCompSign->size[1];
    b_imagCompSign->size[0] = 1;
    b_imagCompSign->size[1] = imagCompSign->size[1];
    emxEnsureCapacity_real_T(b_imagCompSign, i);
    loop_ub_tmp = imagCompSign->size[0] * imagCompSign->size[1] - 1;
    for (i = 0; i <= loop_ub_tmp; i++) {
      b_imagCompSign->data[i] = imagCompSign->data[i];
    }

    writeOutput(b_freq, b_damp, phase1, phase2, b_imagCompSign, static_cast<
                double>(fileid), freq, damp, imagCompSign);
    cfclose(static_cast<double>(fileid));
    emxFree_real_T(&b_imagCompSign);
    emxFree_real_T(&b_damp);
    emxFree_real_T(&b_freq);
  }
}

//
// linearAnalysisModal performs modal analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freq,damp,phase1,phase2,imagCompSign] = linearAnalysisModal(model,mesh,
//                                                  el,displ,Omega,elStorage)
//
//    This function performs a modal analysis of a structural dynamics
//    system and returns freq, damping, and mode shapes
//
//    input:
//    model          = object containing model information
//    mesh           = object containing mesh information
//    el             = object containing element information
//    displ          = displacement vector for use in pre-stressed analysis
//    Omega          = rotor speed (Hz)
//    elStorage      = previously calculated element system matrices
//
//
//    output:
//    freq         = modal frequencies
//    damp         = modal damping ratios
//    phase1       = in phase mode shapes (real part of mode shape)
//    phase2       = out of phase mode shapes (imaginary part of mode shape)
//    imagCompSign = sign of imaginary component of eigenvalues
// Arguments    : double model_RayleighAlpha
//                double model_RayleighBeta
//                double model_BC_numpBC
//                const emxArray_real_T *model_BC_pBC
//                const emxArray_real_T *model_BC_map
//                const emxArray_real_T *model_joint
//                const char model_outFilename_data[]
//                const int model_outFilename_size[2]
//                const emxArray_real_T *model_jointTransform
//                const emxArray_real_T *model_reducedDOFList
//                double mesh_numEl
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                const i_struct_T el
//                const emxArray_real_T *displ
//                const c_emxArray_struct_T *elStorage
//                emxArray_real_T *freq
//                emxArray_real_T *damp
//                emxArray_real_T *phase1
//                emxArray_real_T *phase2
//                emxArray_real_T *imagCompSign
// Return Type  : void
//
void linearAnalysisModal(double model_RayleighAlpha, double model_RayleighBeta,
  double model_BC_numpBC, const emxArray_real_T *model_BC_pBC, const
  emxArray_real_T *model_BC_map, const emxArray_real_T *model_joint, const char
  model_outFilename_data[], const int model_outFilename_size[2], const
  emxArray_real_T *model_jointTransform, const emxArray_real_T
  *model_reducedDOFList, double mesh_numEl, const emxArray_real_T *mesh_x, const
  emxArray_real_T *mesh_y, const emxArray_real_T *mesh_z, const emxArray_real_T *
  mesh_conn, const i_struct_T el, const emxArray_real_T *displ, const
  c_emxArray_struct_T *elStorage, emxArray_real_T *freq, emxArray_real_T *damp,
  emxArray_real_T *phase1, emxArray_real_T *phase2, emxArray_real_T
  *imagCompSign)
{
  emxArray_real_T *Kg;
  double totalNumDOF;
  int i;
  int i1;
  int loop_ub_tmp;
  emxArray_real_T *Mg;
  emxArray_real_T *Cg;
  int b_i;
  emxArray_creal_T *eigVec;
  emxArray_creal_T *eigVal;
  double eldisp_data[12];
  emxArray_real_T *b_r;
  emxArray_real_T *b_r1;
  int b_index;
  double elInput_xloc[2];
  int j;
  double b_mesh_conn[2];
  double elInput_concMass[8];
  double elInput_sectionProps_ac[2];
  double elInput_sectionProps_twist[2];
  double elInput_sectionProps_rhoA[2];
  double elInput_sectionProps_EA[2];
  double elInput_sectionProps_zcm[2];
  double elInput_CN2H[9];
  double elInput_sectionProps_ycm[2];
  double elInput_sectionProps_a[2];
  double elInput_sectionProps_b[2];
  double elInput_sectionProps_a0[2];
  double elx_data[2];
  double elInput_y_data[2];
  double elInput_z_data[2];
  p_struct_T expl_temp;
  double elOutput_Ke[144];
  int elOutput_Me_size[2];
  double elOutput_Me_data[144];
  int elOutput_Ce_size[2];
  double elOutput_Ce_data[144];
  emxArray_real_T *b_r2;
  emxArray_real_T *r3;
  emxArray_creal_T *r4;
  emxArray_creal_T *b_eigVec;
  emxArray_real_T *b_freq;
  emxArray_char_T b_model_outFilename_data;
  signed char fileid;
  emxArray_real_T *b_damp;
  emxArray_real_T *b_imagCompSign;
  emxInit_real_T(&Kg, 2);
  tic();

  // extract mesh information from mesh object
  // extract element order from model
  // extract boundary  condition object from model
  // do initialization
  totalNumDOF = static_cast<double>(mesh_x->size[0]) * 6.0;
  i = Kg->size[0] * Kg->size[1];
  i1 = static_cast<int>(totalNumDOF);
  Kg->size[0] = i1;
  Kg->size[1] = i1;
  emxEnsureCapacity_real_T(Kg, i);
  loop_ub_tmp = i1 * i1;
  for (i = 0; i < loop_ub_tmp; i++) {
    Kg->data[i] = 0.0;
  }

  emxInit_real_T(&Mg, 2);
  i = Mg->size[0] * Mg->size[1];
  Mg->size[0] = i1;
  Mg->size[1] = i1;
  emxEnsureCapacity_real_T(Mg, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    Mg->data[i] = 0.0;
  }

  emxInit_real_T(&Cg, 2);
  i = Cg->size[0] * Cg->size[1];
  Cg->size[0] = i1;
  Cg->size[1] = i1;
  emxEnsureCapacity_real_T(Cg, i);
  for (i = 0; i < loop_ub_tmp; i++) {
    Cg->data[i] = 0.0;
  }

  // extract concentrated nodal terms from model
  i = static_cast<int>(mesh_numEl);
  for (b_i = 0; b_i < i; b_i++) {
    // element loop
    std::memset(&eldisp_data[0], 0, 12U * sizeof(double));

    // Calculate Ke and Fe for element i
    b_index = 0;

    // define element input object flags and element properties from el object
    elInput_xloc[0] = 0.0;
    elInput_xloc[1] = el.elLen->data[b_i];

    // initialize element coordinate list
    // retrieve concentrated nodal terms associated with element
    for (j = 0; j < 2; j++) {
      elInput_sectionProps_ac[j] = el.props->data[b_i].ac[j];
      elInput_sectionProps_twist[j] = el.props->data[b_i].twist[j];
      elInput_sectionProps_rhoA[j] = el.props->data[b_i].rhoA[j];
      elInput_sectionProps_EA[j] = el.props->data[b_i].EA[j];
      elInput_sectionProps_zcm[j] = el.props->data[b_i].zcm[j];
      elInput_sectionProps_ycm[j] = el.props->data[b_i].ycm[j];
      elInput_sectionProps_a[j] = el.props->data[b_i].a[j];
      elInput_sectionProps_b[j] = el.props->data[b_i].b[j];
      elInput_sectionProps_a0[j] = el.props->data[b_i].a0[j];

      // define element coordinates and displacements associated with element
      totalNumDOF = mesh_conn->data[b_i + mesh_conn->size[0] * j];
      loop_ub_tmp = static_cast<int>(totalNumDOF) - 1;
      elx_data[j] = mesh_x->data[loop_ub_tmp];
      elInput_y_data[j] = mesh_y->data[loop_ub_tmp];
      elInput_z_data[j] = mesh_z->data[loop_ub_tmp];
      for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
        eldisp_data[b_index] = 0.0;
        b_index++;
      }

      b_mesh_conn[j] = totalNumDOF;
    }

    ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);

    // assign concentrated nodal terms and coordinates to element input
    // object
    // activate or deactivate rotational effects for element
    // set aeroelastic flag
    std::memset(&elInput_CN2H[0], 0, 9U * sizeof(double));
    elInput_CN2H[0] = 1.0;
    elInput_CN2H[4] = 1.0;
    elInput_CN2H[8] = 1.0;

    // Declare variable on all execution paths
    d_calculateTimoshenkoElementNL(elInput_xloc, elInput_sectionProps_ac,
      elInput_sectionProps_twist, elInput_sectionProps_rhoA,
      elInput_sectionProps_EA, elInput_sectionProps_zcm,
      elInput_sectionProps_ycm, elInput_sectionProps_a, elInput_sectionProps_b,
      elInput_sectionProps_a0, el.psi->data[b_i], el.theta->data[b_i],
      el.roll->data[b_i], elInput_concMass, eldisp_data, elx_data,
      elInput_y_data, elInput_z_data, 0.52359877559829882, false, elInput_CN2H,
      model_RayleighAlpha, model_RayleighBeta, 0.0, &elStorage->data[b_i],
      &expl_temp);
    std::memcpy(&elOutput_Ke[0], &expl_temp.Ke[0], 144U * sizeof(double));
    elOutput_Me_size[0] = expl_temp.Me.size[0];
    elOutput_Me_size[1] = expl_temp.Me.size[1];
    loop_ub_tmp = expl_temp.Me.size[0] * expl_temp.Me.size[1];
    if (0 <= loop_ub_tmp - 1) {
      std::memcpy(&elOutput_Me_data[0], &expl_temp.Me.data[0], loop_ub_tmp *
                  sizeof(double));
    }

    elOutput_Ce_size[0] = expl_temp.Ce.size[0];
    elOutput_Ce_size[1] = expl_temp.Ce.size[1];
    loop_ub_tmp = expl_temp.Ce.size[0] * expl_temp.Ce.size[1];
    if (0 <= loop_ub_tmp - 1) {
      std::memcpy(&elOutput_Ce_data[0], &expl_temp.Ce.data[0], loop_ub_tmp *
                  sizeof(double));
    }

    // do element calculation
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    assemblyMatrixOnly(elOutput_Ke, b_mesh_conn, Kg);

    // assemble element into global stiffness matrix
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    b_assemblyMatrixOnly(elOutput_Me_data, elOutput_Me_size, b_mesh_conn, Mg);

    // assemble element into global mass matrix
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    b_assemblyMatrixOnly(elOutput_Ce_data, elOutput_Ce_size, b_mesh_conn, Cg);

    // assemble element into global damping matrix
  }

  emxInit_creal_T(&eigVec, 2);
  emxInit_creal_T(&eigVal, 2);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_r1, 2);

  // apply general 6x6  mass, damping, and stiffness matrices to nodes
  // ----------------------------------------------------------------------
  // APPLY CONSTRAINT
  b_applyConstraints(Kg, model_jointTransform);

  // modify global matrices for joint constraints using joint transform
  b_applyConstraints(Mg, model_jointTransform);
  b_applyConstraints(Cg, model_jointTransform);

  // APPLY BOUNDARY CONDITIONS
  // apply boundary conditions to global matrices
  applyBCModal(Mg, model_BC_numpBC, model_BC_map, b_r);
  applyBCModal(Cg, model_BC_numpBC, model_BC_map, b_r1);
  applyBCModal(Kg, model_BC_numpBC, model_BC_map, Mg);
  eigSolve(b_r, b_r1, Mg, eigVec, eigVal);

  // ,... %eigensolve of global system
  // model.numModesToExtract,solveFlag);
  // save eigVectors eigVec %save eigenvector for later use (if needed) TODO: Doesn't appear to be used 
  // extract frequency, damping, mode shapes from eigenvalues and vectors
  i = phase1->size[0] * phase1->size[1] * phase1->size[2];
  phase1->size[0] = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0));
  phase1->size[1] = 6;
  phase1->size[2] = eigVal->size[1];
  emxEnsureCapacity_real_T(phase1, i);
  loop_ub_tmp = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0)) *
    6 * eigVal->size[1];
  emxFree_real_T(&b_r1);
  emxFree_real_T(&b_r);
  emxFree_real_T(&Cg);
  emxFree_real_T(&Mg);
  emxFree_real_T(&Kg);
  for (i = 0; i < loop_ub_tmp; i++) {
    phase1->data[i] = 0.0;
  }

  i = phase2->size[0] * phase2->size[1] * phase2->size[2];
  phase2->size[0] = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0));
  phase2->size[1] = 6;
  phase2->size[2] = eigVal->size[1];
  emxEnsureCapacity_real_T(phase2, i);
  loop_ub_tmp = static_cast<int>((static_cast<double>(displ->size[0]) / 6.0)) *
    6 * eigVal->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    phase2->data[i] = 0.0;
  }

  i = eigVal->size[1];
  i1 = freq->size[0] * freq->size[1];
  freq->size[0] = 1;
  freq->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(freq, i1);
  i1 = damp->size[0] * damp->size[1];
  damp->size[0] = 1;
  damp->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(damp, i1);
  i1 = imagCompSign->size[0] * imagCompSign->size[1];
  imagCompSign->size[0] = 1;
  imagCompSign->size[1] = eigVal->size[1];
  emxEnsureCapacity_real_T(imagCompSign, i1);
  emxInit_real_T(&b_r2, 2);
  emxInit_real_T(&r3, 2);
  emxInit_creal_T(&r4, 2);
  emxInit_creal_T(&b_eigVec, 1);
  for (b_i = 0; b_i < i; b_i++) {
    loop_ub_tmp = eigVec->size[0];
    i1 = b_eigVec->size[0];
    b_eigVec->size[0] = eigVec->size[0];
    emxEnsureCapacity_creal_T(b_eigVec, i1);
    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      b_eigVec->data[i1] = eigVec->data[i1 + eigVec->size[0] * b_i];
    }

    extractFreqDamp(eigVal->data[b_i + eigVal->size[0] * b_i], b_eigVec,
                    model_jointTransform, model_reducedDOFList, model_BC_numpBC,
                    model_BC_pBC, &freq->data[b_i], &damp->data[b_i], b_r2, r3,
                    r4);
    loop_ub_tmp = b_r2->size[0];
    b_index = r3->size[0];
    for (i1 = 0; i1 < 6; i1++) {
      for (j = 0; j < loop_ub_tmp; j++) {
        phase1->data[(j + phase1->size[0] * i1) + phase1->size[0] * 6 * b_i] =
          b_r2->data[j + b_r2->size[0] * i1];
      }

      for (j = 0; j < b_index; j++) {
        phase2->data[(j + phase2->size[0] * i1) + phase2->size[0] * 6 * b_i] =
          r3->data[j + r3->size[0] * i1];
      }
    }

    totalNumDOF = eigVal->data[b_i + eigVal->size[0] * b_i].im;
    if (eigVal->data[b_i + eigVal->size[0] * b_i].im < 0.0) {
      totalNumDOF = -1.0;
    } else if (eigVal->data[b_i + eigVal->size[0] * b_i].im > 0.0) {
      totalNumDOF = 1.0;
    } else {
      if (eigVal->data[b_i + eigVal->size[0] * b_i].im == 0.0) {
        totalNumDOF = 0.0;
      }
    }

    imagCompSign->data[b_i] = totalNumDOF;
  }

  emxFree_creal_T(&b_eigVec);
  emxFree_creal_T(&r4);
  emxFree_real_T(&r3);
  emxFree_real_T(&b_r2);
  emxFree_creal_T(&eigVal);
  emxFree_creal_T(&eigVec);
  emxInit_real_T(&b_freq, 2);

  //  %write output
  //  t_modal = toc;
  //  disp('Elapsed time for modal analysis(s):');
  //  disp(t_modal);
  b_model_outFilename_data.data = const_cast<char *>(&model_outFilename_data[0]);
  b_model_outFilename_data.size = const_cast<int *>(&model_outFilename_size[0]);
  b_model_outFilename_data.allocatedSize = -1;
  b_model_outFilename_data.numDimensions = 2;
  b_model_outFilename_data.canFreeData = false;
  fileid = c_cfopen(&b_model_outFilename_data, "wb");
  i = b_freq->size[0] * b_freq->size[1];
  b_freq->size[0] = 1;
  b_freq->size[1] = freq->size[1];
  emxEnsureCapacity_real_T(b_freq, i);
  loop_ub_tmp = freq->size[0] * freq->size[1] - 1;
  for (i = 0; i <= loop_ub_tmp; i++) {
    b_freq->data[i] = freq->data[i];
  }

  emxInit_real_T(&b_damp, 2);
  i = b_damp->size[0] * b_damp->size[1];
  b_damp->size[0] = 1;
  b_damp->size[1] = damp->size[1];
  emxEnsureCapacity_real_T(b_damp, i);
  loop_ub_tmp = damp->size[0] * damp->size[1] - 1;
  for (i = 0; i <= loop_ub_tmp; i++) {
    b_damp->data[i] = damp->data[i];
  }

  emxInit_real_T(&b_imagCompSign, 2);
  i = b_imagCompSign->size[0] * b_imagCompSign->size[1];
  b_imagCompSign->size[0] = 1;
  b_imagCompSign->size[1] = imagCompSign->size[1];
  emxEnsureCapacity_real_T(b_imagCompSign, i);
  loop_ub_tmp = imagCompSign->size[0] * imagCompSign->size[1] - 1;
  for (i = 0; i <= loop_ub_tmp; i++) {
    b_imagCompSign->data[i] = imagCompSign->data[i];
  }

  writeOutput(b_freq, b_damp, phase1, phase2, b_imagCompSign, static_cast<double>
              (fileid), freq, damp, imagCompSign);
  cfclose(static_cast<double>(fileid));
  emxFree_real_T(&b_imagCompSign);
  emxFree_real_T(&b_damp);
  emxFree_real_T(&b_freq);
}

//
// File trailer for linearAnalysisModal.cpp
//
// [EOF]
//
