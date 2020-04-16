//
// File: structuralDynamicsTransient.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "structuralDynamicsTransient.h"
#include "ConcMassAssociatedWithElement.h"
#include "applyBC.h"
#include "assembly.h"
#include "calculateReactionForceAtNode.h"
#include "calculateStrainForElements.h"
#include "calculateTimoshenkoElementNL.h"
#include "mldivide.h"
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "strcmp.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Declarations
static void applyConstraints(emxArray_real_T *Kg, const emxArray_real_T
  *transMatrix);

// Function Definitions

//
// This function transforms a matrix by the transformation matrix to
// enforce joint constraints
//  Kg = transMatrix'*(Kg*transMatrix);
// Arguments    : emxArray_real_T *Kg
//                const emxArray_real_T *transMatrix
// Return Type  : void
//
static void applyConstraints(emxArray_real_T *Kg, const emxArray_real_T
  *transMatrix)
{
  emxArray_real_T *b_transMatrix;
  int i;
  int loop_ub;
  emxArray_real_T *this_d;
  int cend;
  emxArray_int32_T *this_colidx;
  int t9_m;
  emxArray_int32_T *this_rowidx;
  emxArray_real_T *t8_d;
  emxArray_int32_T *t8_colidx;
  emxArray_int32_T *t8_rowidx;
  emxArray_real_T *t9_d;
  emxArray_int32_T *t9_colidx;
  emxArray_int32_T *t9_rowidx;
  int this_m;
  int this_n;
  emxInit_real_T(&b_transMatrix, 2);
  i = b_transMatrix->size[0] * b_transMatrix->size[1];
  b_transMatrix->size[0] = transMatrix->size[1];
  b_transMatrix->size[1] = transMatrix->size[0];
  emxEnsureCapacity_real_T(b_transMatrix, i);
  loop_ub = transMatrix->size[0];
  for (i = 0; i < loop_ub; i++) {
    cend = transMatrix->size[1];
    for (t9_m = 0; t9_m < cend; t9_m++) {
      b_transMatrix->data[t9_m + b_transMatrix->size[0] * i] = transMatrix->
        data[i + transMatrix->size[0] * t9_m];
    }
  }

  emxInit_real_T(&this_d, 1);
  emxInit_int32_T(&this_colidx, 1);
  emxInit_int32_T(&this_rowidx, 1);
  emxInit_real_T(&t8_d, 1);
  emxInit_int32_T(&t8_colidx, 1);
  emxInit_int32_T(&t8_rowidx, 1);
  emxInit_real_T(&t9_d, 1);
  emxInit_int32_T(&t9_colidx, 1);
  emxInit_int32_T(&t9_rowidx, 1);
  d_sparse(b_transMatrix, t8_d, t8_colidx, t8_rowidx, &cend, &loop_ub);
  d_sparse(Kg, this_d, this_colidx, this_rowidx, &this_m, &this_n);
  g_sparse_mtimes(t8_d, t8_colidx, t8_rowidx, cend, this_d, this_colidx,
                  this_rowidx, this_n, t9_d, t9_colidx, t9_rowidx, &t9_m,
                  &loop_ub);
  d_sparse(transMatrix, t8_d, t8_colidx, t8_rowidx, &cend, &loop_ub);
  g_sparse_mtimes(t9_d, t9_colidx, t9_rowidx, t9_m, t8_d, t8_colidx, t8_rowidx,
                  loop_ub, this_d, this_colidx, this_rowidx, &this_m, &this_n);
  i = Kg->size[0] * Kg->size[1];
  Kg->size[0] = this_m;
  Kg->size[1] = this_n;
  emxEnsureCapacity_real_T(Kg, i);
  emxFree_real_T(&b_transMatrix);
  emxFree_int32_T(&t9_rowidx);
  emxFree_int32_T(&t9_colidx);
  emxFree_real_T(&t9_d);
  emxFree_int32_T(&t8_rowidx);
  emxFree_int32_T(&t8_colidx);
  emxFree_real_T(&t8_d);
  for (i = 0; i < this_n; i++) {
    for (t9_m = 0; t9_m < this_m; t9_m++) {
      Kg->data[t9_m + Kg->size[0] * i] = 0.0;
    }
  }

  for (loop_ub = 0; loop_ub < this_n; loop_ub++) {
    cend = this_colidx->data[loop_ub + 1] - 1;
    i = this_colidx->data[loop_ub];
    for (t9_m = i; t9_m <= cend; t9_m++) {
      Kg->data[(this_rowidx->data[t9_m - 1] + Kg->size[0] * loop_ub) - 1] =
        this_d->data[t9_m - 1];
    }
  }

  emxFree_int32_T(&this_rowidx);
  emxFree_int32_T(&this_colidx);
  emxFree_real_T(&this_d);
}

//
// structuralDynamicsTransient perform transient analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,...
//                              dispData,Omega,OmegaDot,time,delta_t,...
//                              elStorage,Fexternal,Fdof,CN2H,rbData)
//
//    This function performs transient structural dynamics analysis.
//
//    input:
//    model      = object containing model data
//    mesh       = object containing mesh data
//    el         = object containing element data
//    dispData   = object containing displacement data
//    Omega      = rotor speed (Hz)
//    OmegaDot   = rotor acceleratin (Hz)
//    time       = current simulation time
//    delta_t    = time step size
//    elStorage  = object containing stored element data
//    Fexternal  = vector containing external force values
//    Fdof       = vector containing global DOF numbering associated with
//                 external force values
//    CN2H       = transformation matrix from inertial frame to hub frame
//    rbData     = vector containing rigid body displacement, velocity, and
//                 acceleration
//
//    output:
//    dispOut       = object containing displacement data at end of time step
//    FReaction_sp1 = vector containing reaction force at turbine base at
//                    end of time step
// Arguments    : const char model_analysisType[3]
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
//                const j_struct_T dispData
//                double Omega
//                double OmegaDot
//                const c_emxArray_struct_T *elStorage
//                const emxArray_real_T *Fexternal
//                const emxArray_real_T *Fdof
//                const double CN2H[9]
//                k_struct_T *dispOut
//                double FReaction_sp1[6]
// Return Type  : void
//
void structuralDynamicsTransient(const char model_analysisType[3], double
  model_RayleighAlpha, double model_RayleighBeta, double model_BC_numpBC, const
  emxArray_real_T *model_BC_pBC, const emxArray_real_T *model_joint, const
  emxArray_real_T *model_jointTransform, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const i_struct_T el, const
  j_struct_T dispData, double Omega, double OmegaDot, const c_emxArray_struct_T *
  elStorage, const emxArray_real_T *Fexternal, const emxArray_real_T *Fdof,
  const double CN2H[9], k_struct_T *dispOut, double FReaction_sp1[6])
{
  double totalNumDOF;
  double elz_data[4];
  double eldisp_data[12];
  double eldispdot_data[12];
  double eldispddot_data[12];
  double eldispiter_data[12];
  emxArray_real_T *disp_s;
  emxArray_real_T *dispdot_s;
  emxArray_real_T *dispddot_s;
  emxArray_real_T *solution;
  int i;
  int aoffset;
  emxArray_real_T *Kg;
  emxArray_real_T *Fg;
  emxArray_real_T *a;
  int m;
  boolean_T elInput_firstIteration;
  int b_i;
  int j;
  double b_mesh_conn[2];
  double elInput_concMass[8];
  double elx_data[2];
  o_struct_T expl_temp;
  double ely_data[2];
  int k;
  int inner;
  b_struct_T b_expl_temp;
  p_struct_T c_expl_temp;
  double elOutput_Ke[144];
  double elOutput_Fe[12];

  // -------- get model information -----------
  totalNumDOF = static_cast<double>(mesh_x->size[0]) * 6.0;

  //  [~,numReducedDOF]=size(model.jointTransform);
  // -----------------------------------------
  // initialize displacements, tolerance, uNorm, iteration count for nonlinear
  // iteration
  elz_data[0] = 0.0;
  elz_data[1] = 0.0;
  std::memset(&eldisp_data[0], 0, 12U * sizeof(double));
  std::memset(&eldispdot_data[0], 0, 12U * sizeof(double));
  std::memset(&eldispddot_data[0], 0, 12U * sizeof(double));
  std::memset(&eldispiter_data[0], 0, 12U * sizeof(double));

  //  if(model.nlOn)
  //       iterationType = model.nlParams.iterationType;
  //  else
  //       iterationType = 'LINEAR';
  //  end
  emxInit_real_T(&disp_s, 1);
  emxInit_real_T(&dispdot_s, 1);
  emxInit_real_T(&dispddot_s, 1);
  emxInit_real_T(&solution, 1);
  if (c_strcmp(model_analysisType)) {
    // ------ newmark integration parameters ---------
    i = disp_s->size[0];
    disp_s->size[0] = dispData.displ_s->size[0];
    emxEnsureCapacity_real_T(disp_s, i);
    aoffset = dispData.displ_s->size[0];
    for (i = 0; i < aoffset; i++) {
      disp_s->data[i] = dispData.displ_s->data[i];
    }

    i = dispdot_s->size[0];
    dispdot_s->size[0] = dispData.displdot_s->size[0];
    emxEnsureCapacity_real_T(dispdot_s, i);
    aoffset = dispData.displdot_s->size[0];
    for (i = 0; i < aoffset; i++) {
      dispdot_s->data[i] = dispData.displdot_s->data[i];
    }

    i = dispddot_s->size[0];
    dispddot_s->size[0] = dispData.displddot_s->size[0];
    emxEnsureCapacity_real_T(dispddot_s, i);
    aoffset = dispData.displddot_s->size[0];
    for (i = 0; i < aoffset; i++) {
      dispddot_s->data[i] = dispData.displddot_s->data[i];
    }

    i = solution->size[0];
    solution->size[0] = dispData.displ_s->size[0];
    emxEnsureCapacity_real_T(solution, i);
    aoffset = dispData.displ_s->size[0];
    for (i = 0; i < aoffset; i++) {
      solution->data[i] = dispData.displ_s->data[i];
    }
  }

  // -----------------------------------------------
  //  Initialize elInput, and DO NOT redundantly re-assign the memory in the
  //  while and for loops below.
  //  %Is not used for this model type, but must be declared
  // Try to force matlab to write the C++ code so as to allocate the memory here and not inside the convergence loop. 
  emxInit_real_T(&Kg, 2);
  emxInit_real_T(&Fg, 1);
  emxInit_real_T(&a, 2);

  // iteration loop
  // ------- intitialization -----------------
  i = Kg->size[0] * Kg->size[1];
  aoffset = static_cast<int>(totalNumDOF);
  Kg->size[0] = aoffset;
  Kg->size[1] = aoffset;
  emxEnsureCapacity_real_T(Kg, i);
  m = aoffset * aoffset;
  for (i = 0; i < m; i++) {
    Kg->data[i] = 0.0;
  }

  // initialize global stiffness and force vector
  i = Fg->size[0];
  Fg->size[0] = aoffset;
  emxEnsureCapacity_real_T(Fg, i);
  for (i = 0; i < aoffset; i++) {
    Fg->data[i] = 0.0;
  }

  // -------------------------------------------
  // ---- element  calculation and assembly ----------------------------------
  i = static_cast<int>(mesh_numEl);
  if (0 <= i - 1) {
    elInput_firstIteration = true;
  }

  for (b_i = 0; b_i < i; b_i++) {
    // Calculate Ke and Fe for element i
    totalNumDOF = 1.0;

    // initialize element data
    // get concentrated terms associated with elemetn
    for (j = 0; j < 2; j++) {
      // get element cooridnates
      aoffset = static_cast<int>(mesh_conn->data[b_i + mesh_conn->size[0] * j])
        - 1;
      elx_data[j] = mesh_x->data[aoffset];
      ely_data[j] = mesh_y->data[aoffset];
      elz_data[j] = mesh_z->data[aoffset];

      // get element nodal displacements at s and s-1 time step
      for (k = 0; k < 6; k++) {
        //                  if(strcmp(analysisType,'TD'))
        //                      eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k); 
        //                      eldisp_sm1(index) = disp_sm1((conn(i,j)-1)*numDOFPerNode + k); 
        //                      eldispiter(index) = displ_iter((conn(i,j)-1)*numDOFPerNode + k); 
        //                  end
        if (c_strcmp(model_analysisType)) {
          m = static_cast<int>(((mesh_conn->data[b_i + mesh_conn->size[0] * j] -
            1.0) * 6.0 + (static_cast<double>(k) + 1.0))) - 1;
          inner = static_cast<int>(totalNumDOF) - 1;
          eldispiter_data[inner] = solution->data[m];
          eldisp_data[inner] = disp_s->data[m];
          eldispdot_data[inner] = dispdot_s->data[m];
          eldispddot_data[inner] = dispddot_s->data[m];
        }

        totalNumDOF++;
      }

      b_mesh_conn[j] = mesh_conn->data[b_i + mesh_conn->size[0] * j];
    }

    ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);

    //  specific to 'TD', but must be declared
    //  specific to 'TNB' , but must be declared
    std::memcpy(&expl_temp.CN2H[0], &CN2H[0], 9U * sizeof(double));
    expl_temp.RayleighBeta = model_RayleighBeta;
    expl_temp.RayleighAlpha = model_RayleighAlpha;
    expl_temp.gravityOn = true;
    expl_temp.airDensity = 0.0;
    expl_temp.aeroForceOn = false;
    expl_temp.aeroElasticOn = false;
    expl_temp.freq = 0.0;
    expl_temp.iterationType[0] = 'D';
    expl_temp.iterationType[1] = 'I';
    expl_temp.preStress = false;
    expl_temp.useDisp = false;
    expl_temp.displ_iter.size[0] = 1;
    expl_temp.displ_iter.size[1] = 12;
    expl_temp.OmegaDot = OmegaDot;
    expl_temp.Omega = Omega;
    expl_temp.omegaDotVec[0] = 0.0;
    expl_temp.omegaVec[0] = 0.0;
    expl_temp.accelVec[0] = 0.0;
    expl_temp.omegaDotVec[1] = 0.0;
    expl_temp.omegaVec[1] = 0.0;
    expl_temp.accelVec[1] = 0.0;
    expl_temp.omegaDotVec[2] = 0.0;
    expl_temp.omegaVec[2] = 0.0;
    expl_temp.accelVec[2] = 0.0;
    expl_temp.z.size[0] = 2;
    expl_temp.z.size[1] = 2;
    expl_temp.z.data[0] = elz_data[0];
    expl_temp.z.data[1] = elz_data[1];
    expl_temp.z.data[2] = 0.0;
    expl_temp.z.data[3] = 0.0;
    expl_temp.y.size[0] = 2;
    expl_temp.x.size[0] = 2;
    expl_temp.y.data[0] = ely_data[0];
    expl_temp.x.data[0] = elx_data[0];
    expl_temp.y.data[1] = ely_data[1];
    expl_temp.x.data[1] = elx_data[1];
    expl_temp.dispddot.size[0] = 1;
    expl_temp.dispddot.size[1] = 12;
    expl_temp.dispdot.size[0] = 1;
    expl_temp.dispdot.size[1] = 12;
    expl_temp.dispm1.size[0] = 1;
    expl_temp.dispm1.size[1] = 12;
    expl_temp.disp.size[0] = 1;
    expl_temp.disp.size[1] = 12;
    for (inner = 0; inner < 12; inner++) {
      expl_temp.displ_iter.data[inner] = eldispiter_data[inner];
      expl_temp.dispddot.data[inner] = eldispddot_data[inner];
      expl_temp.dispdot.data[inner] = eldispdot_data[inner];
      expl_temp.dispm1.data[inner] = 0.0;
      expl_temp.disp.data[inner] = eldisp_data[inner];
      expl_temp.concLoad[inner] = 0.0;
      expl_temp.concStiff[inner] = 0.0;
    }

    std::memcpy(&expl_temp.concMass[0], &elInput_concMass[0], 8U * sizeof(double));
    expl_temp.firstIteration = elInput_firstIteration;
    expl_temp.aeroSweepAngle = 0.0;
    expl_temp.rollAngle = el.roll->data[b_i];
    expl_temp.coneAngle = el.theta->data[b_i];
    expl_temp.sweepAngle = el.psi->data[b_i];
    expl_temp.sectionProps = el.props->data[b_i];
    expl_temp.xloc[0] = 0.0;
    expl_temp.xloc[1] = el.elLen->data[b_i];
    expl_temp.timeInt.delta_t = 0.002;
    expl_temp.timeInt.a1 = 0.001;
    expl_temp.timeInt.a2 = 0.001;
    expl_temp.timeInt.a3 = 1.0E+6;
    expl_temp.timeInt.a4 = 2000.0;
    expl_temp.timeInt.a5 = 1.0;
    expl_temp.timeInt.a6 = 1000.0;
    expl_temp.timeInt.a7 = 1.0;
    expl_temp.timeInt.a8 = 0.0;
    expl_temp.modalFlag = true;
    expl_temp.elementOrder = 1.0;
    expl_temp.analysisType[0] = model_analysisType[0];
    expl_temp.analysisType[1] = model_analysisType[1];
    expl_temp.analysisType[2] = model_analysisType[2];
    calculateTimoshenkoElementNL(&expl_temp, &elStorage->data[b_i], &c_expl_temp);
    std::memcpy(&elOutput_Ke[0], &c_expl_temp.Ke[0], 144U * sizeof(double));
    std::memcpy(&elOutput_Fe[0], &c_expl_temp.Fe[0], 12U * sizeof(double));

    // calculate timoshenko element
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    assembly(elOutput_Ke, elOutput_Fe, b_mesh_conn, Kg, Fg);

    // assemble element stiffness matrix and force vector
    //          Erestotal = Erestotal + elOutput.Eres;
    // ................................................
  }

  // ------- end element calculation and assembly ------------------
  //     %%
  // ----------------------------------------------------------------------
  //     %%
  // Apply external loads to structure
  if ((Fexternal->size[0] == 0) || (Fexternal->size[1] == 0)) {
    aoffset = 0;
  } else {
    aoffset = Fexternal->size[1];
  }

  for (b_i = 0; b_i < aoffset; b_i++) {
    if (c_strcmp(model_analysisType)) {
      Fg->data[static_cast<int>(Fdof->data[b_i]) - 1] += Fexternal->data[b_i];
    }
  }

  // ------ apply constraints on system -----------------------------------
  applyConstraints(Kg, model_jointTransform);

  // This function transforms a vector by the transformation matrix to
  // enforce joint constraints
  i = a->size[0] * a->size[1];
  a->size[0] = model_jointTransform->size[1];
  a->size[1] = model_jointTransform->size[0];
  emxEnsureCapacity_real_T(a, i);
  aoffset = model_jointTransform->size[0];
  for (i = 0; i < aoffset; i++) {
    m = model_jointTransform->size[1];
    for (inner = 0; inner < m; inner++) {
      a->data[inner + a->size[0] * i] = model_jointTransform->data[i +
        model_jointTransform->size[0] * inner];
    }
  }

  i = solution->size[0];
  solution->size[0] = Fg->size[0];
  emxEnsureCapacity_real_T(solution, i);
  aoffset = Fg->size[0];
  for (i = 0; i < aoffset; i++) {
    solution->data[i] = Fg->data[i];
  }

  if ((a->size[1] == 1) || (Fg->size[0] == 1)) {
    i = solution->size[0];
    solution->size[0] = a->size[0];
    emxEnsureCapacity_real_T(solution, i);
    aoffset = a->size[0];
    for (i = 0; i < aoffset; i++) {
      solution->data[i] = 0.0;
      m = a->size[1];
      for (inner = 0; inner < m; inner++) {
        solution->data[i] += a->data[i + a->size[0] * inner] * Fg->data[inner];
      }
    }

    i = Fg->size[0];
    Fg->size[0] = solution->size[0];
    emxEnsureCapacity_real_T(Fg, i);
    aoffset = solution->size[0];
    for (i = 0; i < aoffset; i++) {
      Fg->data[i] = solution->data[i];
    }
  } else {
    m = a->size[0];
    inner = a->size[1];
    i = Fg->size[0];
    Fg->size[0] = a->size[0];
    emxEnsureCapacity_real_T(Fg, i);
    for (b_i = 0; b_i < m; b_i++) {
      Fg->data[b_i] = 0.0;
    }

    for (k = 0; k < inner; k++) {
      aoffset = k * m;
      for (b_i = 0; b_i < m; b_i++) {
        Fg->data[b_i] += solution->data[k] * a->data[aoffset + b_i];
      }
    }
  }

  // ----------------------------------------------------------------------
  //     %%
  // Apply BCs to global system
  applyBC(Kg, Fg, model_BC_numpBC, model_BC_pBC);
  b_mldivide(Kg, Fg);

  // solve for displacements
  if ((model_jointTransform->size[1] == 1) || (Fg->size[0] == 1)) {
    i = solution->size[0];
    solution->size[0] = model_jointTransform->size[0];
    emxEnsureCapacity_real_T(solution, i);
    aoffset = model_jointTransform->size[0];
    for (i = 0; i < aoffset; i++) {
      solution->data[i] = 0.0;
      m = model_jointTransform->size[1];
      for (inner = 0; inner < m; inner++) {
        solution->data[i] += model_jointTransform->data[i +
          model_jointTransform->size[0] * inner] * Fg->data[inner];
      }
    }
  } else {
    m = model_jointTransform->size[0];
    inner = model_jointTransform->size[1];
    i = solution->size[0];
    solution->size[0] = model_jointTransform->size[0];
    emxEnsureCapacity_real_T(solution, i);
    for (b_i = 0; b_i < m; b_i++) {
      solution->data[b_i] = 0.0;
    }

    for (k = 0; k < inner; k++) {
      aoffset = k * m;
      for (b_i = 0; b_i < m; b_i++) {
        solution->data[b_i] += Fg->data[k] * model_jointTransform->data[aoffset
          + b_i];
      }
    }
  }

  // transform to full dof listing
  emxFree_real_T(&a);
  emxFree_real_T(&Fg);
  emxFree_real_T(&Kg);
  emxFree_real_T(&dispddot_s);
  emxFree_real_T(&dispdot_s);
  emxFree_real_T(&disp_s);

  // Calculate reaction at turbine base (hardwired to node number 1)
  b_expl_temp.a8 = 0.0;
  b_expl_temp.a7 = 1.0;
  b_expl_temp.a6 = 1000.0;
  b_expl_temp.a5 = 1.0;
  b_expl_temp.a4 = 2000.0;
  b_expl_temp.a3 = 1.0E+6;
  b_expl_temp.a2 = 0.001;
  b_expl_temp.a1 = 0.001;
  b_expl_temp.delta_t = 0.002;
  calculateReactionForceAtNode(model_analysisType, model_RayleighAlpha,
    model_RayleighBeta, model_joint, mesh_x, mesh_y, mesh_z, mesh_conn, el.props,
    el.elLen, el.psi, el.theta, el.roll, elStorage, &b_expl_temp, dispData,
    solution, Omega, OmegaDot, CN2H, FReaction_sp1);

  // Calculate strain
  calculateStrainForElements(mesh_numEl, mesh_conn, el.props, el.elLen, el.psi,
    el.theta, el.roll, solution, false, dispOut->elStrain);
  i = dispOut->displ_sp1->size[0];
  dispOut->displ_sp1->size[0] = solution->size[0];
  emxEnsureCapacity_real_T(dispOut->displ_sp1, i);
  aoffset = solution->size[0];
  for (i = 0; i < aoffset; i++) {
    dispOut->displ_sp1->data[i] = solution->data[i];
  }

  // store displacement vector in dispOut
  //  Specific to TNB, but must be declared
  aoffset = solution->size[0];
  for (i = 0; i < aoffset; i++) {
    solution->data[i] = (1.0E+6 * (solution->data[i] - dispData.displ_s->data[i])
                         - 2000.0 * dispData.displdot_s->data[i]) -
      dispData.displddot_s->data[i];
  }

  // store velocity vector in dispOut
  i = dispOut->displddot_sp1->size[0];
  dispOut->displddot_sp1->size[0] = solution->size[0];
  emxEnsureCapacity_real_T(dispOut->displddot_sp1, i);
  aoffset = solution->size[0];
  for (i = 0; i < aoffset; i++) {
    dispOut->displddot_sp1->data[i] = solution->data[i];
  }

  i = dispOut->displdot_sp1->size[0];
  dispOut->displdot_sp1->size[0] = dispData.displdot_s->size[0];
  emxEnsureCapacity_real_T(dispOut->displdot_sp1, i);
  aoffset = dispData.displdot_s->size[0];
  for (i = 0; i < aoffset; i++) {
    dispOut->displdot_sp1->data[i] = (dispData.displdot_s->data[i] + 0.001 *
      dispData.displddot_s->data[i]) + 0.001 * solution->data[i];
  }

  emxFree_real_T(&solution);

  // store acceleration vector in dispOut
}

//
// File trailer for structuralDynamicsTransient.cpp
//
// [EOF]
//
