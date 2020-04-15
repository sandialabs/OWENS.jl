//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: modalExecAuto.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "modalExecAuto.h"
#include "fminbnd.h"
#include "initialElementCalculations.h"
#include "linearAnalysisModal.h"
#include "rt_nonfinite.h"
#include "staticAnalysis.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <stdio.h>
#include <string.h>

// Function Definitions

//
// modalExec  Executive function for automated flutter/modal analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freq,damp]=modalExecAuto(model,mesh,el,displ,omegaArray,OmegaStart)
//
//    This function executes modal analysis.
//
//    input:
//    model          = object containing model information
//    mesh           = object containing mesh information
//    el             = object containing element information
//    displ          = displacement vector for use in pre-stressed analysis
//    Omega          = rotor speed (Hz)
//    OmegaStart     = rotor speed (Hz) from previous analysis if stepping
//                     through various rotor speeds, may be useful in load
//                     stepping
//
//    output:
//    freq           = modal frequencies of system
//    damp           = modal damping ratios of system
// Arguments    : const char model_analysisType[2]
//                const char model_nlParams_iterationType[2]
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
//                const h_struct_T mesh
//                const i_struct_T el
//                const emxArray_real_T *displ
//                double freq_data[]
//                int freq_size[2]
//                double damp_data[]
//                int damp_size[2]
// Return Type  : void
//
void modalExecAuto(const char model_analysisType[2], const char
                   model_nlParams_iterationType[2], double model_RayleighAlpha,
                   double model_RayleighBeta, double model_BC_numpBC, const
                   emxArray_real_T *model_BC_pBC, const emxArray_real_T
                   *model_BC_map, const emxArray_real_T *model_joint, const char
                   model_outFilename_data[], const int model_outFilename_size[2],
                   const emxArray_real_T *model_jointTransform, const
                   emxArray_real_T *model_reducedDOFList, const h_struct_T mesh,
                   const i_struct_T el, const emxArray_real_T *displ, double
                   freq_data[], int freq_size[2], double damp_data[], int
                   damp_size[2])
{
  c_emxArray_struct_T *elStorage;
  emxArray_real_T *b_displ;
  int i;
  int loop_ub;
  b_emxArray_struct_T *unusedU0;
  boolean_T staticAnalysisSuccessful;
  emxArray_real_T *freqOrig;
  emxArray_real_T *unusedU1;
  emxArray_real_T *unusedU2;
  emxArray_real_T *unusedU3;
  emxArray_real_T *unusedU4;
  emxArray_real_T *unusedU7;
  double model_guessFreq;
  double fval;
  double exitflag;
  double expl_temp;
  double b_expl_temp;
  emxInit_struct_T2(&elStorage, 2);
  emxInit_real_T(&b_displ, 1);
  initialElementCalculations(model_joint, el.props, el.elLen, el.psi, el.theta,
    el.roll, mesh.numEl, mesh.x, mesh.y, mesh.z, mesh.conn, elStorage);

  // performs initial element calculations
  //  [structureMass,structureMOI,structureMassCenter]=calculateStructureMassProps(elStorage); %calculate mass properties of structure 
  freq_size[0] = 1;
  freq_size[1] = 4;
  damp_size[0] = 1;
  damp_size[1] = 4;
  freq_data[0] = 0.0;
  damp_data[0] = 0.0;
  freq_data[1] = 0.0;
  damp_data[1] = 0.0;
  freq_data[2] = 0.0;
  damp_data[2] = 0.0;
  freq_data[3] = 0.0;
  damp_data[3] = 0.0;

  // loops over rotor speeds of interest
  // Do nonlinear iteration if needed
  i = b_displ->size[0];
  b_displ->size[0] = displ->size[0];
  emxEnsureCapacity_real_T(b_displ, i);
  loop_ub = displ->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_displ->data[i] = displ->data[i];
  }

  emxInit_struct_T1(&unusedU0, 2);
  staticAnalysis(model_nlParams_iterationType, model_RayleighAlpha,
                 model_RayleighBeta, model_BC_numpBC, model_BC_pBC, model_joint,
                 model_jointTransform, mesh.numEl, mesh.x, mesh.y, mesh.z,
                 mesh.conn, el, b_displ, elStorage, unusedU0,
                 &staticAnalysisSuccessful);

  // performs static analysis about specified operating condition
  //          save displnl displ
  emxFree_struct_T1(&unusedU0);
  if (staticAnalysisSuccessful) {
    emxInit_real_T(&freqOrig, 2);
    emxInit_real_T(&unusedU1, 2);
    emxInit_real_T(&unusedU2, 3);
    emxInit_real_T(&unusedU3, 3);
    emxInit_real_T(&unusedU4, 2);
    emxInit_real_T(&unusedU7, 2);

    // get nominal frequencies without aeroelastic effects
    // for the first rotor speed perform a modal analysis with initial guess frequency of zero to obtain frequency estimate 
    b_linearAnalysisModal(model_analysisType, false, 0.0, model_RayleighAlpha,
                          model_RayleighBeta, model_BC_numpBC, model_BC_pBC,
                          model_BC_map, model_joint, model_outFilename_data,
                          model_outFilename_size, model_jointTransform,
                          model_reducedDOFList, mesh.numEl, mesh.x, mesh.y,
                          mesh.z, mesh.conn, el, b_displ, elStorage, freqOrig,
                          unusedU1, unusedU2, unusedU3, unusedU4);

    // If zero, then the output freq is always the input freq and no need to iterate 
    //                  [finalfreq,fval,exitflag,output] = fzero(rootfun,freqOrig(i),options); 
    fminbnd(2.0, model_analysisType, true, model_RayleighAlpha,
            model_RayleighBeta, model_BC_numpBC, model_BC_pBC, model_BC_map,
            model_joint, model_outFilename_data, model_outFilename_size,
            model_jointTransform, model_reducedDOFList, mesh, el, b_displ,
            elStorage, freqOrig->data[1] * 2.0, &model_guessFreq, &fval,
            &exitflag, &expl_temp, &b_expl_temp);
    b_linearAnalysisModal(model_analysisType, true, model_guessFreq,
                          model_RayleighAlpha, model_RayleighBeta,
                          model_BC_numpBC, model_BC_pBC, model_BC_map,
                          model_joint, model_outFilename_data,
                          model_outFilename_size, model_jointTransform,
                          model_reducedDOFList, mesh.numEl, mesh.x, mesh.y,
                          mesh.z, mesh.conn, el, b_displ, elStorage, unusedU1,
                          unusedU4, unusedU2, unusedU3, unusedU7);

    // performs modal analysis during iteration
    freq_data[1] = unusedU1->data[1];

    // if converged, store and move on to next frequency
    damp_data[1] = unusedU4->data[1];

    //              model.guessFreq = freqOrig(i); %do p-k iteration for flutter analysis of modes of interest 
    //              converged = false;
    //              while(~converged) %TODO: this is terribly inefficient to calculate all of the modes, but only converge on one at a time, need to use some sort of ND root finder. 
    //                  [freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); %performs modal analysis during iteration 
    //
    printf("%f\n", 1.0);
    fflush(stdout);
    printf("%f\n", 2.0);
    fflush(stdout);

    //                  if(abs(freq(i) - model.guessFreq)<tol)  %check if modal frequency is converged 
    //                      converged = true;
    //                      %             i
    //
    //                      convergedFreq(j,i) = freq(i);    %if converged, store and move on to next frequency 
    //                      convergedDamp(j,i) = damp(i);
    //                  else
    //                      model.guessFreq = 0.5*(freq(i) + model.guessFreq); %if not converged select another guess frequency 
    //                  end
    //              end
    // If zero, then the output freq is always the input freq and no need to iterate 
    //                  [finalfreq,fval,exitflag,output] = fzero(rootfun,freqOrig(i),options); 
    fminbnd(4.0, model_analysisType, true, model_RayleighAlpha,
            model_RayleighBeta, model_BC_numpBC, model_BC_pBC, model_BC_map,
            model_joint, model_outFilename_data, model_outFilename_size,
            model_jointTransform, model_reducedDOFList, mesh, el, b_displ,
            elStorage, freqOrig->data[3] * 2.0, &model_guessFreq, &fval,
            &exitflag, &expl_temp, &b_expl_temp);
    b_linearAnalysisModal(model_analysisType, true, model_guessFreq,
                          model_RayleighAlpha, model_RayleighBeta,
                          model_BC_numpBC, model_BC_pBC, model_BC_map,
                          model_joint, model_outFilename_data,
                          model_outFilename_size, model_jointTransform,
                          model_reducedDOFList, mesh.numEl, mesh.x, mesh.y,
                          mesh.z, mesh.conn, el, b_displ, elStorage, unusedU1,
                          unusedU4, unusedU2, unusedU3, unusedU7);

    // performs modal analysis during iteration
    freq_data[3] = unusedU1->data[3];

    // if converged, store and move on to next frequency
    damp_data[3] = unusedU4->data[3];

    //              model.guessFreq = freqOrig(i); %do p-k iteration for flutter analysis of modes of interest 
    //              converged = false;
    //              while(~converged) %TODO: this is terribly inefficient to calculate all of the modes, but only converge on one at a time, need to use some sort of ND root finder. 
    //                  [freq,damp,~,~,~] = linearAnalysisModal(model,mesh,el,displ,Omega,elStorage); %performs modal analysis during iteration 
    //
    printf("%f\n", 1.0);
    fflush(stdout);
    printf("%f\n", 4.0);
    fflush(stdout);

    //                  if(abs(freq(i) - model.guessFreq)<tol)  %check if modal frequency is converged 
    //                      converged = true;
    //                      %             i
    //
    //                      convergedFreq(j,i) = freq(i);    %if converged, store and move on to next frequency 
    //                      convergedDamp(j,i) = damp(i);
    //                  else
    //                      model.guessFreq = 0.5*(freq(i) + model.guessFreq); %if not converged select another guess frequency 
    //                  end
    //              end
    emxFree_real_T(&unusedU7);
    emxFree_real_T(&unusedU4);
    emxFree_real_T(&unusedU3);
    emxFree_real_T(&unusedU2);
    emxFree_real_T(&unusedU1);
    emxFree_real_T(&freqOrig);
  }

  emxFree_real_T(&b_displ);
  emxFree_struct_T2(&elStorage);

  //  save flutterRun freq damp omegaArray %save frequency, damping, and rotor speed array to .mat file 
}

//
// File trailer for modalExecAuto.cpp
//
// [EOF]
//
