//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: staticAnalysis.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "staticAnalysis.h"
#include "ConcMassAssociatedWithElement.h"
#include "applyBC.h"
#include "assembly.h"
#include "calculateStrainForElements.h"
#include "calculateTimoshenkoElementNL.h"
#include "mldivide.h"
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "sparse1.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "tic.h"
#include "toc.h"
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <string.h>

// Function Declarations
static void applyConstraintsVec(emxArray_real_T *Fg, const emxArray_real_T
  *transMatrix);

// Function Definitions

//
// This function transforms a vector by the transformation matrix to
// enforce joint constraints
// Arguments    : emxArray_real_T *Fg
//                const emxArray_real_T *transMatrix
// Return Type  : void
//
static void applyConstraintsVec(emxArray_real_T *Fg, const emxArray_real_T
  *transMatrix)
{
  emxArray_real_T *transMatrix_d;
  emxArray_int32_T *transMatrix_colidx;
  emxArray_int32_T *transMatrix_rowidx;
  emxArray_real_T *a_d;
  emxArray_int32_T *a_colidx;
  emxArray_int32_T *a_rowidx;
  int apend;
  int nap;
  int apend1;
  int a_n;
  int i;
  int acol;
  double bc;
  int ap;
  emxInit_real_T(&transMatrix_d, 1);
  emxInit_int32_T(&transMatrix_colidx, 1);
  emxInit_int32_T(&transMatrix_rowidx, 1);
  emxInit_real_T(&a_d, 1);
  emxInit_int32_T(&a_colidx, 1);
  emxInit_int32_T(&a_rowidx, 1);
  d_sparse(transMatrix, transMatrix_d, transMatrix_colidx, transMatrix_rowidx,
           &apend, &nap);
  sparse_ctranspose(transMatrix_d, transMatrix_colidx, transMatrix_rowidx, apend,
                    nap, a_d, a_colidx, a_rowidx, &apend1, &a_n);
  i = transMatrix_d->size[0];
  transMatrix_d->size[0] = Fg->size[0];
  emxEnsureCapacity_real_T(transMatrix_d, i);
  apend = Fg->size[0];
  emxFree_int32_T(&transMatrix_rowidx);
  emxFree_int32_T(&transMatrix_colidx);
  for (i = 0; i < apend; i++) {
    transMatrix_d->data[i] = Fg->data[i];
  }

  i = Fg->size[0];
  Fg->size[0] = apend1;
  emxEnsureCapacity_real_T(Fg, i);
  for (i = 0; i < apend1; i++) {
    Fg->data[i] = 0.0;
  }

  if ((a_n != 0) && (apend1 != 0) && (a_colidx->data[a_colidx->size[0] - 1] - 1
       != 0)) {
    for (acol = 0; acol < a_n; acol++) {
      bc = transMatrix_d->data[acol];
      i = a_colidx->data[acol];
      apend = a_colidx->data[acol + 1];
      nap = apend - a_colidx->data[acol];
      if (nap >= 4) {
        apend1 = (apend - nap) + ((nap / 4) << 2);
        for (ap = i; ap <= apend1 - 1; ap += 4) {
          nap = a_rowidx->data[ap - 1] - 1;
          Fg->data[nap] += a_d->data[ap - 1] * bc;
          Fg->data[a_rowidx->data[ap] - 1] += a_d->data[ap] * bc;
          nap = a_rowidx->data[ap + 1] - 1;
          Fg->data[nap] += a_d->data[ap + 1] * bc;
          nap = a_rowidx->data[ap + 2] - 1;
          Fg->data[nap] += a_d->data[ap + 2] * bc;
        }

        nap = a_colidx->data[acol + 1] - 1;
        for (ap = apend1; ap <= nap; ap++) {
          i = a_rowidx->data[ap - 1] - 1;
          Fg->data[i] += a_d->data[ap - 1] * bc;
        }
      } else {
        apend--;
        for (ap = i; ap <= apend; ap++) {
          nap = a_rowidx->data[ap - 1] - 1;
          Fg->data[nap] += a_d->data[ap - 1] * bc;
        }
      }
    }
  }

  emxFree_int32_T(&a_rowidx);
  emxFree_int32_T(&a_colidx);
  emxFree_real_T(&a_d);
  emxFree_real_T(&transMatrix_d);
}

//
// This function transforms a matrix by the transformation matrix to
// enforce joint constraints
// Arguments    : emxArray_real_T *Kg
//                const emxArray_real_T *transMatrix
// Return Type  : void
//
void b_applyConstraints(emxArray_real_T *Kg, const emxArray_real_T *transMatrix)
{
  emxArray_real_T *transMatrix_d;
  emxArray_int32_T *transMatrix_colidx;
  emxArray_int32_T *transMatrix_rowidx;
  emxArray_real_T *Kg_d;
  emxArray_int32_T *Kg_colidx;
  emxArray_int32_T *Kg_rowidx;
  emxArray_real_T *t0_d;
  emxArray_int32_T *t0_colidx;
  emxArray_int32_T *t0_rowidx;
  emxArray_real_T *t1_d;
  emxArray_int32_T *t1_colidx;
  emxArray_int32_T *t1_rowidx;
  int transMatrix_m;
  int transMatrix_n;
  int cend;
  int Kg_n;
  int expl_temp;
  emxInit_real_T(&transMatrix_d, 1);
  emxInit_int32_T(&transMatrix_colidx, 1);
  emxInit_int32_T(&transMatrix_rowidx, 1);
  emxInit_real_T(&Kg_d, 1);
  emxInit_int32_T(&Kg_colidx, 1);
  emxInit_int32_T(&Kg_rowidx, 1);
  emxInit_real_T(&t0_d, 1);
  emxInit_int32_T(&t0_colidx, 1);
  emxInit_int32_T(&t0_rowidx, 1);
  emxInit_real_T(&t1_d, 1);
  emxInit_int32_T(&t1_colidx, 1);
  emxInit_int32_T(&t1_rowidx, 1);
  d_sparse(transMatrix, transMatrix_d, transMatrix_colidx, transMatrix_rowidx,
           &transMatrix_m, &transMatrix_n);
  d_sparse(Kg, Kg_d, Kg_colidx, Kg_rowidx, &cend, &Kg_n);
  sparse_ctranspose(transMatrix_d, transMatrix_colidx, transMatrix_rowidx,
                    transMatrix_m, transMatrix_n, t0_d, t0_colidx, t0_rowidx,
                    &cend, &expl_temp);
  g_sparse_mtimes(t0_d, t0_colidx, t0_rowidx, cend, Kg_d, Kg_colidx, Kg_rowidx,
                  Kg_n, t1_d, t1_colidx, t1_rowidx, &transMatrix_m, &expl_temp);
  g_sparse_mtimes(t1_d, t1_colidx, t1_rowidx, transMatrix_m, transMatrix_d,
                  transMatrix_colidx, transMatrix_rowidx, transMatrix_n, Kg_d,
                  Kg_colidx, Kg_rowidx, &cend, &Kg_n);
  transMatrix_n = Kg->size[0] * Kg->size[1];
  Kg->size[0] = cend;
  Kg->size[1] = Kg_n;
  emxEnsureCapacity_real_T(Kg, transMatrix_n);
  emxFree_int32_T(&t1_rowidx);
  emxFree_int32_T(&t1_colidx);
  emxFree_real_T(&t1_d);
  emxFree_int32_T(&t0_rowidx);
  emxFree_int32_T(&t0_colidx);
  emxFree_real_T(&t0_d);
  emxFree_int32_T(&transMatrix_rowidx);
  emxFree_int32_T(&transMatrix_colidx);
  emxFree_real_T(&transMatrix_d);
  for (transMatrix_n = 0; transMatrix_n < Kg_n; transMatrix_n++) {
    for (transMatrix_m = 0; transMatrix_m < cend; transMatrix_m++) {
      Kg->data[transMatrix_m + Kg->size[0] * transMatrix_n] = 0.0;
    }
  }

  for (transMatrix_m = 0; transMatrix_m < Kg_n; transMatrix_m++) {
    cend = Kg_colidx->data[transMatrix_m + 1] - 1;
    transMatrix_n = Kg_colidx->data[transMatrix_m];
    for (expl_temp = transMatrix_n; expl_temp <= cend; expl_temp++) {
      Kg->data[(Kg_rowidx->data[expl_temp - 1] + Kg->size[0] * transMatrix_m) -
        1] = Kg_d->data[expl_temp - 1];
    }
  }

  emxFree_int32_T(&Kg_rowidx);
  emxFree_int32_T(&Kg_colidx);
  emxFree_real_T(&Kg_d);
}

//
// staticAnalysis performs static analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [displ,staticAnalysisSuccessful]=staticAnalysis(model,mesh,el,displ,...
//                                     Omega,OmegaStart,elStorag
//
//    This function performs a static analysis and returns displacement
//    values and a flag denoting successful/unsuccessful analysis
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
//    elStorage      = previously calculated element system matrices
//
//
//    output:
//    displ                    = vector of displacemetns
//    staticAnalysisSuccessful = boolean flag denoting successful static
//                               analysis
// Arguments    : const char model_nlParams_iterationType[2]
//                double model_RayleighAlpha
//                double model_RayleighBeta
//                double model_BC_numpBC
//                const emxArray_real_T *model_BC_pBC
//                const emxArray_real_T *model_joint
//                const emxArray_real_T *model_jointTransform
//                double mesh_numEl
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                const i_struct_T el
//                emxArray_real_T *displ
//                const c_emxArray_struct_T *elStorage
//                b_emxArray_struct_T *elStrain
//                boolean_T *staticAnalysisSuccessful
// Return Type  : void
//
void staticAnalysis(const char model_nlParams_iterationType[2], double
                    model_RayleighAlpha, double model_RayleighBeta, double
                    model_BC_numpBC, const emxArray_real_T *model_BC_pBC, const
                    emxArray_real_T *model_joint, const emxArray_real_T
                    *model_jointTransform, double mesh_numEl, const
                    emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y,
                    const emxArray_real_T *mesh_z, const emxArray_real_T
                    *mesh_conn, const i_struct_T el, emxArray_real_T *displ,
                    const c_emxArray_struct_T *elStorage, b_emxArray_struct_T
                    *elStrain, boolean_T *staticAnalysisSuccessful)
{
  emxArray_real_T *BC_pBC;
  double BC_numpBC;
  int i;
  int msgId;
  emxArray_real_T *displPrev;
  double totalNumDOF;
  int loadStepCount;
  emxArray_real_T *dispOld;
  boolean_T staticAnalysisComplete;
  double loadStepPrev;
  double loadStep;
  char iterationType_data[2];
  emxArray_real_T *Kg;
  emxArray_real_T *Fg;
  emxArray_real_T *dispEval;
  double uNorm;
  int iterationCount;
  int b_index;
  int b_i;
  int inner;
  char elInput_iterationType_data[2];
  double eldisp_data[12];
  double elInput_xloc[2];
  double b_mesh_conn[2];
  double elInput_concMass[8];
  double elInput_sectionProps_ac[2];
  double elInput_sectionProps_twist[2];
  static const signed char b_iv[18] = { 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 0, 0,
    0, 0, 0, 0 };

  double elInput_sectionProps_rhoA[2];
  double elInput_sectionProps_EA[2];
  double elInput_sectionProps_zcm[2];
  boolean_T b_bool;
  double elInput_sectionProps_ycm[2];
  double elInput_sectionProps_a[2];
  double elInput_sectionProps_b[2];
  double elInput_sectionProps_a0[2];
  double elx_data[2];
  double elInput_y_data[2];
  double elInput_z_data[2];
  static const double elInput_CN2H[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0 };

  p_struct_T expl_temp;
  int exitg1;
  double elOutput_Ke[144];
  double elOutput_Fe[12];
  static const char b_cv[2] = { 'N', 'R' };

  int k;
  static const char cv1[2] = { 'D', 'I' };

  double denom;
  double a;
  emxInit_real_T(&BC_pBC, 2);
  tic();

  // extract mesh information from mesh object
  // extract element order from model
  BC_numpBC = model_BC_numpBC;
  i = BC_pBC->size[0] * BC_pBC->size[1];
  BC_pBC->size[0] = model_BC_pBC->size[0];
  BC_pBC->size[1] = 3;
  emxEnsureCapacity_real_T(BC_pBC, i);
  msgId = model_BC_pBC->size[0] * model_BC_pBC->size[1];
  for (i = 0; i < msgId; i++) {
    BC_pBC->data[i] = model_BC_pBC->data[i];
  }

  emxInit_real_T(&displPrev, 1);

  // extract boundary  condition object from model
  // do initialization
  totalNumDOF = static_cast<double>(mesh_x->size[0]) * 6.0;

  // extract concentrated nodal terms from model
  // extract extra copy of concentrated nodal terms from model
  // load stepping paramters
  loadStepCount = 1;
  i = displPrev->size[0];
  displPrev->size[0] = displ->size[0];
  emxEnsureCapacity_real_T(displPrev, i);
  msgId = displ->size[0];
  for (i = 0; i < msgId; i++) {
    displPrev->data[i] = displ->data[i];
  }

  emxInit_real_T(&dispOld, 1);

  // copy of initial displacement vector
  staticAnalysisComplete = false;

  // initialize to false
  loadStepPrev = 0.0;
  loadStep = 1.0;
  i = dispOld->size[0];
  dispOld->size[0] = displ->size[0];
  emxEnsureCapacity_real_T(dispOld, i);
  msgId = displ->size[0];
  for (i = 0; i < msgId; i++) {
    dispOld->data[i] = displ->data[i];
  }

  // initialize dispOld, first iteration logic below for accurate calculations
  // .........................................................................
  *staticAnalysisSuccessful = false;

  // initialize variable
  iterationType_data[0] = 'N';
  iterationType_data[1] = 'A';

  // Initialize iterationType so it is in scope.  Is set below.
  emxInit_real_T(&Kg, 2);
  emxInit_real_T(&Fg, 1);
  emxInit_real_T(&dispEval, 1);
  while ((!staticAnalysisComplete) && (loadStepCount < 20)) {
    //  staticAnalysisSuccessful = false; %initialize staticAnalysisSuccessful flag 
    uNorm = 1.0;

    // initialize norm for convergence check
    iterationCount = 0;

    // initialize iteration count
    while ((uNorm > 1.0E-6) && (iterationCount < 50)) {
      // iteration loop (convergence tolerance of 1.0e-6)
      msgId = static_cast<int>(totalNumDOF);
      i = Kg->size[0] * Kg->size[1];
      Kg->size[0] = msgId;
      Kg->size[1] = msgId;
      emxEnsureCapacity_real_T(Kg, i);
      b_index = msgId * msgId;
      for (i = 0; i < b_index; i++) {
        Kg->data[i] = 0.0;
      }

      // initialize global stiffness matrix
      i = Fg->size[0];
      Fg->size[0] = msgId;
      emxEnsureCapacity_real_T(Fg, i);
      for (i = 0; i < msgId; i++) {
        Fg->data[i] = 0.0;
      }

      // initialize global force vector
      i = static_cast<int>(mesh_numEl);
      if (0 <= i - 1) {
        if (iterationCount >= 1) {
          // option for acceleration of iterative procedure (gamma = 0 or gamma=0.5 are typical) 
          inner = dispEval->size[0];
          dispEval->size[0] = dispOld->size[0];
          emxEnsureCapacity_real_T(dispEval, inner);
          msgId = dispOld->size[0];
          for (inner = 0; inner < msgId; inner++) {
            dispEval->data[inner] = dispOld->data[inner] * 0.0 + displ->
              data[inner];
          }
        } else {
          inner = dispEval->size[0];
          dispEval->size[0] = displ->size[0];
          emxEnsureCapacity_real_T(dispEval, inner);
          msgId = displ->size[0];
          for (inner = 0; inner < msgId; inner++) {
            dispEval->data[inner] = displ->data[inner];
          }
        }

        iterationType_data[0] = model_nlParams_iterationType[0];
        iterationType_data[1] = model_nlParams_iterationType[1];
      }

      for (b_i = 0; b_i < i; b_i++) {
        // Calculate Ke and Fe for element i
        b_index = 0;

        // define element input object flags and element properties from el object 
        // define nonlinear iteration type
        elInput_iterationType_data[0] = model_nlParams_iterationType[0];
        elInput_iterationType_data[1] = model_nlParams_iterationType[1];
        std::memset(&eldisp_data[0], 0, 12U * sizeof(double));

        //  %Initialize variable, but isn't used
        elInput_xloc[0] = 0.0;
        elInput_xloc[1] = el.elLen->data[b_i];

        // initialize element coordinate list
        // retrieve concentrated nodal terms associated with element
        for (inner = 0; inner < 2; inner++) {
          elInput_sectionProps_ac[inner] = el.props->data[b_i].ac[inner];
          elInput_sectionProps_twist[inner] = el.props->data[b_i].twist[inner];
          elInput_sectionProps_rhoA[inner] = el.props->data[b_i].rhoA[inner];
          elInput_sectionProps_EA[inner] = el.props->data[b_i].EA[inner];
          elInput_sectionProps_zcm[inner] = el.props->data[b_i].zcm[inner];
          elInput_sectionProps_ycm[inner] = el.props->data[b_i].ycm[inner];
          elInput_sectionProps_a[inner] = el.props->data[b_i].a[inner];
          elInput_sectionProps_b[inner] = el.props->data[b_i].b[inner];
          elInput_sectionProps_a0[inner] = el.props->data[b_i].a0[inner];

          // define element coordinates and displacements associated with element 
          msgId = static_cast<int>(mesh_conn->data[b_i + mesh_conn->size[0] *
            inner]) - 1;
          elx_data[inner] = mesh_x->data[msgId];
          elInput_y_data[inner] = mesh_y->data[msgId];
          elInput_z_data[inner] = mesh_z->data[msgId];
          for (k = 0; k < 6; k++) {
            eldisp_data[b_index] = dispEval->data[static_cast<int>
              (((mesh_conn->data[b_i + mesh_conn->size[0] * inner] - 1.0) * 6.0
                + (static_cast<double>(k) + 1.0))) - 1];
            b_index++;
          }

          b_mesh_conn[inner] = mesh_conn->data[b_i + mesh_conn->size[0] * inner];
        }

        ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);

        // assign concentrated nodal terms and coordinates to element input
        // object
        // activate or deactivate rotational effects for element
        // deactivate flutter type input
        // turn off for static calculation
        // not used for static analysis but in general the calculate element function will look for this data in the input struct 
        c_calculateTimoshenkoElementNL(elInput_iterationType_data, elInput_xloc,
          elInput_sectionProps_ac, elInput_sectionProps_twist,
          elInput_sectionProps_rhoA, elInput_sectionProps_EA,
          elInput_sectionProps_zcm, elInput_sectionProps_ycm,
          elInput_sectionProps_a, elInput_sectionProps_b,
          elInput_sectionProps_a0, el.psi->data[b_i], el.theta->data[b_i],
          el.roll->data[b_i], elInput_concMass, eldisp_data, elx_data,
          elInput_y_data, elInput_z_data, elInput_CN2H, model_RayleighAlpha,
          model_RayleighBeta, loadStep, &elStorage->data[b_i], &expl_temp);
        std::memcpy(&elOutput_Ke[0], &expl_temp.Ke[0], 144U * sizeof(double));
        std::memcpy(&elOutput_Fe[0], &expl_temp.Fe[0], 12U * sizeof(double));

        // do element calculation
        b_mesh_conn[0] = mesh_conn->data[b_i];
        b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
        assembly(elOutput_Ke, elOutput_Fe, b_mesh_conn, Kg, Fg);

        // assemble element i into global system
      }

      // get arbitrary external loads from externalForcingStatic() function
      // ----------------------------------------------------------------------
      // apply general 6x6  mass, damping, and stiffness matrices to nodes
      // APPLY BOUNDARY CONDITIONS
      b_applyConstraints(Kg, model_jointTransform);

      // modify global stiffness matrix for joint constraints using joint transform 
      applyConstraintsVec(Fg, model_jointTransform);

      // modify global force vector for joint constraints using joint transform
      if (BC_numpBC == 0.0) {
        BC_numpBC = 6.0;
        i = BC_pBC->size[0] * BC_pBC->size[1];
        BC_pBC->size[0] = 6;
        BC_pBC->size[1] = 3;
        emxEnsureCapacity_real_T(BC_pBC, i);
        for (i = 0; i < 18; i++) {
          BC_pBC->data[i] = b_iv[i];
        }
      }

      applyBC(Kg, Fg, BC_numpBC, BC_pBC);

      // apply boundary conditions to global stiffness matrix and force vector
      i = dispOld->size[0];
      dispOld->size[0] = displ->size[0];
      emxEnsureCapacity_real_T(dispOld, i);
      msgId = displ->size[0];
      for (i = 0; i < msgId; i++) {
        dispOld->data[i] = displ->data[i];
      }

      // assign displacement vector from previous iteration
      b_bool = false;
      msgId = 0;
      do {
        exitg1 = 0;
        if (msgId < 2) {
          if (iterationType_data[msgId] != b_cv[msgId]) {
            exitg1 = 1;
          } else {
            msgId++;
          }
        } else {
          b_bool = true;
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (b_bool) {
        // system solve, norm calculation for newton-raphson iteration
        b_mldivide(Kg, Fg);
        if ((model_jointTransform->size[1] == 1) || (Fg->size[0] == 1)) {
          i = dispEval->size[0];
          dispEval->size[0] = model_jointTransform->size[0];
          emxEnsureCapacity_real_T(dispEval, i);
          msgId = model_jointTransform->size[0];
          for (i = 0; i < msgId; i++) {
            dispEval->data[i] = 0.0;
            b_index = model_jointTransform->size[1];
            for (inner = 0; inner < b_index; inner++) {
              dispEval->data[i] += model_jointTransform->data[i +
                model_jointTransform->size[0] * inner] * Fg->data[inner];
            }
          }
        } else {
          b_index = model_jointTransform->size[0];
          inner = model_jointTransform->size[1];
          i = dispEval->size[0];
          dispEval->size[0] = model_jointTransform->size[0];
          emxEnsureCapacity_real_T(dispEval, i);
          for (b_i = 0; b_i < b_index; b_i++) {
            dispEval->data[b_i] = 0.0;
          }

          for (k = 0; k < inner; k++) {
            msgId = k * b_index;
            for (b_i = 0; b_i < b_index; b_i++) {
              dispEval->data[b_i] += Fg->data[k] * model_jointTransform->
                data[msgId + b_i];
            }
          }
        }

        msgId = displ->size[0];
        for (i = 0; i < msgId; i++) {
          displ->data[i] += dispEval->data[i];
        }

        i = dispEval->size[0];
        dispEval->size[0] = displ->size[0];
        emxEnsureCapacity_real_T(dispEval, i);
        msgId = displ->size[0];
        for (i = 0; i < msgId; i++) {
          dispEval->data[i] = displ->data[i] - dispEval->data[i];
        }

        // Calculates a relative norm between two vectors: u and uPrev
        uNorm = 0.0;
        denom = 0.0;
        i = displ->size[0];
        for (b_i = 0; b_i < i; b_i++) {
          a = displ->data[b_i] - dispEval->data[b_i];
          uNorm += a * a;
          denom += displ->data[b_i] * displ->data[b_i];
        }

        uNorm = std::sqrt(uNorm / denom);
      } else {
        b_bool = false;
        msgId = 0;
        do {
          exitg1 = 0;
          if (msgId < 2) {
            if (iterationType_data[msgId] != cv1[msgId]) {
              exitg1 = 1;
            } else {
              msgId++;
            }
          } else {
            b_bool = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);

        if (b_bool) {
          // system solve, norm calculation for direct iteration
          i = dispEval->size[0];
          dispEval->size[0] = displ->size[0];
          emxEnsureCapacity_real_T(dispEval, i);
          msgId = displ->size[0];
          for (i = 0; i < msgId; i++) {
            dispEval->data[i] = displ->data[i];
          }

          b_mldivide(Kg, Fg);
          if ((model_jointTransform->size[1] == 1) || (Fg->size[0] == 1)) {
            i = displ->size[0];
            displ->size[0] = model_jointTransform->size[0];
            emxEnsureCapacity_real_T(displ, i);
            msgId = model_jointTransform->size[0];
            for (i = 0; i < msgId; i++) {
              displ->data[i] = 0.0;
              b_index = model_jointTransform->size[1];
              for (inner = 0; inner < b_index; inner++) {
                displ->data[i] += model_jointTransform->data[i +
                  model_jointTransform->size[0] * inner] * Fg->data[inner];
              }
            }
          } else {
            b_index = model_jointTransform->size[0];
            inner = model_jointTransform->size[1];
            i = displ->size[0];
            displ->size[0] = model_jointTransform->size[0];
            emxEnsureCapacity_real_T(displ, i);
            for (b_i = 0; b_i < b_index; b_i++) {
              displ->data[b_i] = 0.0;
            }

            for (k = 0; k < inner; k++) {
              msgId = k * b_index;
              for (b_i = 0; b_i < b_index; b_i++) {
                displ->data[b_i] += Fg->data[k] * model_jointTransform->
                  data[msgId + b_i];
              }
            }
          }

          // Calculates a relative norm between two vectors: u and uPrev
          uNorm = 0.0;
          denom = 0.0;
          i = displ->size[0];
          for (b_i = 0; b_i < i; b_i++) {
            a = displ->data[b_i] - dispEval->data[b_i];
            uNorm += a * a;
            denom += displ->data[b_i] * displ->data[b_i];
          }

          uNorm = std::sqrt(uNorm / denom);
          msgId = displ->size[0];
          for (i = 0; i < msgId; i++) {
            displ->data[i] = 0.5 * displ->data[i] + 0.5 * dispEval->data[i];
          }
        } else {
          // system solve for linear case
          b_mldivide(Kg, Fg);
          if ((model_jointTransform->size[1] == 1) || (Fg->size[0] == 1)) {
            i = displ->size[0];
            displ->size[0] = model_jointTransform->size[0];
            emxEnsureCapacity_real_T(displ, i);
            msgId = model_jointTransform->size[0];
            for (i = 0; i < msgId; i++) {
              displ->data[i] = 0.0;
              b_index = model_jointTransform->size[1];
              for (inner = 0; inner < b_index; inner++) {
                displ->data[i] += model_jointTransform->data[i +
                  model_jointTransform->size[0] * inner] * Fg->data[inner];
              }
            }
          } else {
            b_index = model_jointTransform->size[0];
            inner = model_jointTransform->size[1];
            i = displ->size[0];
            displ->size[0] = model_jointTransform->size[0];
            emxEnsureCapacity_real_T(displ, i);
            for (b_i = 0; b_i < b_index; b_i++) {
              displ->data[b_i] = 0.0;
            }

            for (k = 0; k < inner; k++) {
              msgId = k * b_index;
              for (b_i = 0; b_i < b_index; b_i++) {
                displ->data[b_i] += Fg->data[k] * model_jointTransform->
                  data[msgId + b_i];
              }
            }
          }

          uNorm = 0.0;
        }
      }

      iterationCount++;

      // increment iteration count
    }

    loadStepCount++;

    // increment load step count
    // update load step whether adaptive or prescribed
    // updateLoadStep    updates load step for static nonlinear analysis
    //  **********************************************************************
    //  *                   Part of the SNL OWENS Toolkit                    *
    //  * Developed by Sandia National Laboratories Wind Energy Technologies *
    //  *             See license.txt for disclaimer information             *
    //  **********************************************************************
    //    [loadStep,loadStepPrev,displ,displCopy,staticAnalysisSuccessful,...
    //     staticAnalysisComplete] = updateLoadStep(iterationCount,...
    //     loadStepParams,loadStep,loadStepPrev,loadStepCount,displCopy,displ)
    //
    //    This function updates the load stepping parameter whether through means 
    //    of adaptive loadstepping or a specified load step profile.
    //
    //       input:
    //       iterationCount      = number of iterations for current load step
    //       loadStepParams      = struct containing load step parameters
    //       loadStep            = load step value for current load step
    //       loadStepPrev        = load step value for previous load st ep
    //       loadStepCount       = number of load steps performed up to this
    //                             point
    //       displPrev           = converged displacement vector form previous
    //                             load step
    //       displ               = displacement vector at current load step
    //
    //       output:
    //       loadStep            = new load step value
    //       loadStepPrev        = load step value for previous load step
    //       displ               = most up to date displacement vector in load
    //                             stepping procedure
    //       displPrev           = displacement vector at end of previous load
    //                             step
    //       staticAnalysisSuccessful = boolean flag, true if load
    //                                  step was completed successfully
    //       staticAnalysisComplete   = boolean flag, true if analysis is complete 
    //
    //
    // initialize variable
    // for adaptive load stepping option
    // check if maximum number of load steps has been exceeded.
    // calculate new loadstep adaptively
    uNorm = loadStep;

    // This function performs updates a loadstep adaptively.
    if (iterationCount >= 50) {
      // check for exceeding max iterations
      msgId = 1;

      // corresponds to a message  saying load step was unsuccessful and is being reduced 
      *staticAnalysisSuccessful = false;
    } else {
      if (loadStep == 1.0) {
        msgId = 2;

        // corresponds to a message saying loadstep was successful and load stepping is finished 
      } else {
        msgId = 3;

        // corresponds to a message saying loads tep was successful and analysis is proceeding to next load step 
      }

      *staticAnalysisSuccessful = true;
      staticAnalysisComplete = (loadStep == 1.0);
      loadStepPrev = loadStep;

      // update previous load step and end loadstep
      uNorm = 1.0;
    }

    // make copy of load step for later checks
    uNorm = 0.5 * (uNorm + loadStepPrev);

    // update load step
    if (std::abs(uNorm - loadStepPrev) < 0.05) {
      // enforces delta load step is not below the minimum specified value
      if (uNorm < loadStepPrev) {
        uNorm -= 0.05;
      } else {
        uNorm += 0.05;
      }
    }

    if (uNorm < 0.05) {
      // check that load step is not below the minimum specified load step
      uNorm = 0.05;
    } else {
      if (uNorm > 1.0) {
        // if load step has extended beyond 1.0, set to 1.0
        uNorm = 1.0;
      }
    }

    // print loadstepping message to command line
    if (msgId == 1) {
      printf("Max iterations exceeded for loadstep = %f.\t\t Reducing load step to %f.\n",
             loadStep, uNorm);
      fflush(stdout);
    } else if (msgId == 2) {
      printf("Nonlinear iteration successful for loadstep = %f.\t Nonlinear static analysis complete.\n",
             uNorm);
      fflush(stdout);
    } else {
      printf("Nonlinear iteration successful for loadstep = %f.\t Increasing loadstep size to %f.\n",
             loadStep, uNorm);
      fflush(stdout);
    }

    loadStep = uNorm;
    if (*staticAnalysisSuccessful) {
      i = displPrev->size[0];
      displPrev->size[0] = displ->size[0];
      emxEnsureCapacity_real_T(displPrev, i);
      msgId = displ->size[0];
      for (i = 0; i < msgId; i++) {
        displPrev->data[i] = displ->data[i];
      }

      // update displacementPrev variable if previous load step was successful.
    } else {
      i = displ->size[0];
      displ->size[0] = displPrev->size[0];
      emxEnsureCapacity_real_T(displ, i);
      msgId = displPrev->size[0];
      for (i = 0; i < msgId; i++) {
        displ->data[i] = displPrev->data[i];
      }

      // reset to displ vector to that at start of load step if load step was unsuccessful 
    }
  }

  emxFree_real_T(&dispEval);
  emxFree_real_T(&Fg);
  emxFree_real_T(&Kg);
  emxFree_real_T(&dispOld);
  emxFree_real_T(&displPrev);
  emxFree_real_T(&BC_pBC);
  toc();

  //      reactionNodeNumber = 1; %place holder for nodal reaction force
  //      [FReaction] = calculateReactionForceAtNode(reactionNodeNumber,model,mesh,el,... 
  //          elStorage,[],[],displ,[],Omega,0,[]);
  calculateStrainForElements(mesh_numEl, mesh_conn, el.props, el.elLen, el.psi,
    el.theta, el.roll, displ, true, elStrain);
}

//
// File trailer for staticAnalysis.cpp
//
// [EOF]
//
