//
// File: calculateBCMap.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "calculateBCMap.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Definitions

//
// calculateBCMap   calculates a boundary condition map
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [bcMap] = calculateBCMap(numpBC,pBC,numDofPerNode,reducedDofList)
//
//    This function creates a boundary condition map between full and reduced
//    dof listing as a result of constraints.
//
//       input:
//       numpBC            = number of boundary conditions
//       pBC               = array of boundary  condition data
//       numDofPerNode     = number of degrees of freedom per node
//       reducedDofList    = array of reduced DOF numbering
//
//       output:
//       elStorage         = map for boundary conditions between full and
//                           reduced dof list
// Arguments    : double numpBC
//                const emxArray_real_T *pBC
//                const emxArray_real_T *reducedDofList
//                emxArray_real_T *bcMap
// Return Type  : void
//
void calculateBCMap(double numpBC, const emxArray_real_T *pBC, const
                    emxArray_real_T *reducedDofList, emxArray_real_T *bcMap)
{
  emxArray_real_T *constrainedDof;
  int i;
  int b_i;
  double b_index;
  double a;
  boolean_T tf;
  int k;
  boolean_T exitg1;
  double absx;
  int exponent;
  emxInit_real_T(&constrainedDof, 1);
  i = static_cast<int>(numpBC);
  b_i = constrainedDof->size[0];
  constrainedDof->size[0] = i;
  emxEnsureCapacity_real_T(constrainedDof, b_i);
  for (b_i = 0; b_i < i; b_i++) {
    constrainedDof->data[b_i] = (pBC->data[b_i] - 1.0) * 6.0 + pBC->data[b_i +
      pBC->size[0]];

    // creates an array of constrained DOFs
  }

  sort(constrainedDof);
  i = bcMap->size[0];
  bcMap->size[0] = reducedDofList->size[1];
  emxEnsureCapacity_real_T(bcMap, i);
  b_index = 1.0;
  i = reducedDofList->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    a = reducedDofList->data[b_i];
    tf = false;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= constrainedDof->size[0] - 1)) {
      absx = std::abs(constrainedDof->data[k] / 2.0);
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

      if ((std::abs(constrainedDof->data[k] - a) < absx) || (rtIsInf(a) &&
           rtIsInf(constrainedDof->data[k]) && ((a > 0.0) ==
            (constrainedDof->data[k] > 0.0)))) {
        tf = true;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (tf) {
      // searches reduced DOF for constrained DOFs
      bcMap->data[b_i] = -1.0;
    } else {
      bcMap->data[b_i] = b_index;
      b_index++;
    }
  }

  emxFree_real_T(&constrainedDof);
}

//
// File trailer for calculateBCMap.cpp
//
// [EOF]
//
