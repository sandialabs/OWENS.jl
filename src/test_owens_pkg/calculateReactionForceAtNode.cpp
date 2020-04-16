//
// File: calculateReactionForceAtNode.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "calculateReactionForceAtNode.h"
#include "elementPostProcess.h"
#include "findElementsAssociatedWithNodeNumber.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// calculateReactionForceAtNode calculates reaction force at a node
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [cummulativeForce] = calculateReactionForceAtNode(nodeNum,model,mesh,...
//     el,elStorage,timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
//
//    This function calculates the reaction force at a node by post
//    processing all element associated with a node through connectivity or
//    joint constraints.
//
//    input:
//    nodeNum    = node number joint constraints are desired at
//    model      = object containing model data
//    mesh       = object containing mesh data
//    elStorage  = object containing stored element data
//    el         = object containing element data
//    timeInt    = object containing time integration parameters
//    dispData   = object containing displacement data
//    displ_iter = converged displacement solution
//    rbData     = vector containing rigid body displacement, velocity, and
//                 acceleration
//    Omega      = rotor speed (Hz)
//    OmegaDot   = rotor acceleratin (Hz)
//    CN2H       = transformation matrix from inertial frame to hub frame
//
//    output:
//    cummulativeForce  = vector containing reaction force at nodeNum
// Arguments    : const char model_analysisType[3]
//                double model_RayleighAlpha
//                double model_RayleighBeta
//                const emxArray_real_T *model_joint
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                const emxArray_struct_T *el_props
//                const emxArray_real_T *el_elLen
//                const emxArray_real_T *el_psi
//                const emxArray_real_T *el_theta
//                const emxArray_real_T *el_roll
//                const c_emxArray_struct_T *elStorage
//                const b_struct_T *timeInt
//                const j_struct_T dispData
//                const emxArray_real_T *displ_iter
//                double Omega
//                double OmegaDot
//                const double CN2H[9]
//                double cummulativeForce[6]
// Return Type  : void
//
void calculateReactionForceAtNode(const char model_analysisType[3], double
  model_RayleighAlpha, double model_RayleighBeta, const emxArray_real_T
  *model_joint, const emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y,
  const emxArray_real_T *mesh_z, const emxArray_real_T *mesh_conn, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const c_emxArray_struct_T *elStorage, const
  b_struct_T *timeInt, const j_struct_T dispData, const emxArray_real_T
  *displ_iter, double Omega, double OmegaDot, const double CN2H[9], double
  cummulativeForce[6])
{
  boolean_T b;
  int i;
  emxArray_real_T *elList;
  emxArray_real_T *elListLocalNodeNumbers;
  int b_i;
  double b_dv[9];
  double Fpp[12];
  double b_elListLocalNodeNumbers;
  int i1;
  b = false;

  // get connectivity list
  for (i = 0; i < 6; i++) {
    cummulativeForce[i] = 0.0;
  }

  emxInit_real_T(&elList, 2);
  emxInit_real_T(&elListLocalNodeNumbers, 2);

  // initialize force at node
  // find elements associated with nodeNum due to mesh connectivity or
  // constraints
  c_findElementsAssociatedWithNod(1.0, mesh_conn, model_joint, elList,
    elListLocalNodeNumbers);

  // process elements for nodal reaction forces and compile to find total
  // reaction force at specified node
  b_i = elList->size[1];
  for (i = 0; i < b_i; i++) {
    if (!b) {
      std::memset(&b_dv[0], 0, 9U * sizeof(double));
      b = true;
    }

    elementPostProcess(elList->data[i], model_analysisType, model_RayleighAlpha,
                       model_RayleighBeta, model_joint, mesh_x, mesh_y, mesh_z,
                       mesh_conn, el_props, el_elLen, el_psi, el_theta, el_roll,
                       elStorage, timeInt, dispData, displ_iter, b_dv, Omega,
                       OmegaDot, CN2H, Fpp);
    b_elListLocalNodeNumbers = (elListLocalNodeNumbers->data[i] - 1.0) * 6.0;
    for (i1 = 0; i1 < 6; i1++) {
      cummulativeForce[i1] += Fpp[static_cast<int>((b_elListLocalNodeNumbers + (
        static_cast<double>(i1) + 1.0))) - 1];
    }
  }

  emxFree_real_T(&elListLocalNodeNumbers);
  emxFree_real_T(&elList);
}

//
// File trailer for calculateReactionForceAtNode.cpp
//
// [EOF]
//
