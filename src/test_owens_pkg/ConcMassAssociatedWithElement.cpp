//
// File: ConcMassAssociatedWithElement.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//

// Include Files
#include "ConcMassAssociatedWithElement.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Definitions

//
// ConcMassAssociatedWithElement gets concentrated terms associated w/ el
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [mass,stiff,load,modJoint,modNodalMassTerms,...
//     modNodalStiffnessTerms,modNodalLoads] = ...
//      ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,...
//      nodalStiffnessTerms,nodalLoads)
//
//    This function compiles concentrated mass, stiffness, and load
//    associated with a node from both ndl and joint files. The mod*
//    variables are passed back with these terms removed to prevent
//    duplicate application of shared nodal terms between elements
//
//       input:
//       conn                = connectivity list for element
//       joint               = joint array for nodal terms
//       nodalMassTerms      = listing of concentrated nodal mass terms
//       nodalStiffnessTerms = listing of concentrated nodal stiffness terms
//       nodalLoads          = listing of concentrated nodal loads terms
//
//
//       output:
//       mass                = array of concentrated mass associated with element
//       stiff               = array of concentrated stiffness associated with
//                             element
//       load                = array of concentrated loads associated with element
//       modJoint            = modified joint object removing nodal terms that
//                             have/will be applied to the element calculations
//       modNodalMassTerms   = modified nodal mass object removing nodal terms that
//                             have/will be applied to the element calculations
//       modalStiffnessTerms = modified nodal stiffness object removing nodal terms that
//                             have/will be applied to the element calculations
//       modNodalLoads       = modified nodal loads object removing nodal terms that
//                             have/will be applied to the element calculations
// Arguments    : const double conn[2]
//                const emxArray_real_T *joint
//                double mass[8]
// Return Type  : void
//
void ConcMassAssociatedWithElement(const double conn[2], const emxArray_real_T
  *joint, double mass[8])
{
  double node1;
  double node2;
  double mass1;
  double mass2;
  emxArray_boolean_T *node1flag;
  emxArray_boolean_T *node2flag;
  int i;
  int loop_ub;
  int i1;
  double absx;
  int exponent;
  int b_exponent;
  node1 = conn[0];

  // define node #1 and node #2
  node2 = conn[1];
  mass1 = 0.0;

  // initialize concentrated mass amd moi for nodes
  mass2 = 0.0;

  // initialize concentrated stifness for nodes
  // initialize concentrated loads/moments
  // create copies of joint, and nodal mass, stiffness, loads arrays
  // get number of joints in model
  emxInit_boolean_T(&node1flag, 1);
  emxInit_boolean_T(&node2flag, 1);
  if (joint->size[0] > 0) {
    i = joint->size[0] - 1;
    loop_ub = joint->size[0];
    i1 = node1flag->size[0];
    node1flag->size[0] = loop_ub;
    emxEnsureCapacity_boolean_T(node1flag, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      node1flag->data[i1] = false;
    }

    for (loop_ub = 0; loop_ub <= i; loop_ub++) {
      absx = std::abs(node1 / 2.0);
      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &exponent);
          absx = std::ldexp(1.0, exponent - 53);
        }
      } else {
        absx = rtNaN;
      }

      if ((std::abs(node1 - joint->data[loop_ub + joint->size[0]]) < absx) ||
          (rtIsInf(joint->data[loop_ub + joint->size[0]]) && rtIsInf(node1) &&
           ((joint->data[loop_ub + joint->size[0]] > 0.0) == (node1 > 0.0)))) {
        node1flag->data[loop_ub] = true;
      }
    }

    // see if nodes are associated with a joint constraint as a master node
    i = joint->size[0] - 1;
    loop_ub = joint->size[0];
    i1 = node2flag->size[0];
    node2flag->size[0] = loop_ub;
    emxEnsureCapacity_boolean_T(node2flag, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      node2flag->data[i1] = false;
    }

    for (loop_ub = 0; loop_ub <= i; loop_ub++) {
      absx = std::abs(node2 / 2.0);
      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &b_exponent);
          absx = std::ldexp(1.0, b_exponent - 53);
        }
      } else {
        absx = rtNaN;
      }

      if ((std::abs(node2 - joint->data[loop_ub + joint->size[0]]) < absx) ||
          (rtIsInf(joint->data[loop_ub + joint->size[0]]) && rtIsInf(node2) &&
           ((joint->data[loop_ub + joint->size[0]] > 0.0) == (node2 > 0.0)))) {
        node2flag->data[loop_ub] = true;
      }
    }
  } else {
    i = node1flag->size[0];
    node1flag->size[0] = 1;
    emxEnsureCapacity_boolean_T(node1flag, i);
    node1flag->data[0] = false;
    i = node2flag->size[0];
    node2flag->size[0] = 1;
    emxEnsureCapacity_boolean_T(node2flag, i);
    node2flag->data[0] = false;
  }

  i = joint->size[0];
  for (loop_ub = 0; loop_ub < i; loop_ub++) {
    // if nodes are associated with joint constraint, use (if any) mass and stiffness specification from the joint file 
    if (node1flag->data[loop_ub]) {
      mass1 += joint->data[loop_ub + joint->size[0] * 4];

      //              stiff1x = stiff1x + joint(i,6);
      //              stiff1y = stiff1y + joint(i,6);
      //              stiff1z = stiff1z + joint(i,6);
      //              stiff1mx = stiff1mx + joint(i,6);
      //              stiff1my = stiff1my + joint(i,6);
      //              stiff1mz = stiff1mz + joint(i,6);
      //              modJoint(i,5)=0.0;
      //              modJoint(i,6)=0.0;
    }

    if (node2flag->data[loop_ub]) {
      mass2 += joint->data[loop_ub + joint->size[0] * 4];

      //              stiff2x = stiff2x + joint(i,6);
      //              stiff2y = stiff2y + joint(i,6);
      //              stiff2z = stiff2z + joint(i,6);
      //              stiff2mx = stiff2mx + joint(i,6);
      //              stiff2my = stiff2my + joint(i,6);
      //              stiff2mz = stiff2mz + joint(i,6);
      //              modJoint(i,5)=0.0;
      //              modJoint(i,6)=0.0;
    }
  }

  emxFree_boolean_T(&node2flag);
  emxFree_boolean_T(&node1flag);

  // apply concentrated mass/stiffness from NDL file
  // compile nodal concentrated terms into mass, stiffness, and load arrays
  mass[0] = mass1;
  mass[4] = mass2;
  mass[1] = 0.0;
  mass[2] = 0.0;
  mass[3] = 0.0;
  mass[5] = 0.0;
  mass[6] = 0.0;
  mass[7] = 0.0;
}

//
// File trailer for ConcMassAssociatedWithElement.cpp
//
// [EOF]
//
