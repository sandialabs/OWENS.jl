//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: applyBCModal.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//

// Include Files
#include "applyBCModal.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// applyBCModal Applies boundary conditions to system for modal analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
//
//    This function applies boundary conditions to a system matrix for modal
//    analysis
//
//       input:
//       K             = assembled global system matrix
//       BC            = struct of boundary condition information
//       numDofPerNode = number of degrees of freedom per node
// Arguments    : const emxArray_real_T *K
//                double numpBC
//                const emxArray_real_T *bcMap
//                emxArray_real_T *Knew
// Return Type  : void
//
void applyBCModal(const emxArray_real_T *K, double numpBC, const emxArray_real_T
                  *bcMap, emxArray_real_T *Knew)
{
  emxArray_boolean_T *x;
  unsigned int b_index;
  int i;
  int vlen;
  int nz;
  int k;
  emxArray_uint32_T *indVec;
  emxInit_boolean_T(&x, 1);

  //       output:
  //       K             = global system matrix with boundary conditions
  //       dofVector     = reduced DOF vector after imposing BCs
  //  indVec = zeros(numEq-numpBC,1); %can't pre-initialize...puts zeros in map
  //  that causes problems when creating Knew
  b_index = 1U;
  i = x->size[0];
  x->size[0] = bcMap->size[0];
  emxEnsureCapacity_boolean_T(x, i);
  vlen = bcMap->size[0];
  for (i = 0; i < vlen; i++) {
    x->data[i] = (bcMap->data[i] != -1.0);
  }

  vlen = x->size[0];
  if (x->size[0] == 0) {
    nz = 0;
  } else {
    nz = x->data[0];
    for (k = 2; k <= vlen; k++) {
      nz += x->data[k - 1];
    }
  }

  emxFree_boolean_T(&x);
  emxInit_uint32_T(&indVec, 2);
  i = indVec->size[0] * indVec->size[1];
  indVec->size[0] = 1;
  indVec->size[1] = nz;
  emxEnsureCapacity_uint32_T(indVec, i);
  for (i = 0; i < nz; i++) {
    indVec->data[i] = 0U;
  }

  i = bcMap->size[0];
  for (vlen = 0; vlen < i; vlen++) {
    if (bcMap->data[vlen] != -1.0) {
      indVec->data[static_cast<int>(b_index) - 1] = static_cast<unsigned int>
        ((vlen + 1));
      b_index++;
    }
  }

  if (numpBC > 0.0) {
    i = Knew->size[0] * Knew->size[1];
    Knew->size[0] = indVec->size[1];
    Knew->size[1] = indVec->size[1];
    emxEnsureCapacity_real_T(Knew, i);
    vlen = indVec->size[1];
    for (i = 0; i < vlen; i++) {
      k = indVec->size[1];
      for (nz = 0; nz < k; nz++) {
        Knew->data[nz + Knew->size[0] * i] = K->data[(static_cast<int>
          (indVec->data[nz]) + K->size[0] * (static_cast<int>(indVec->data[i]) -
          1)) - 1];
      }
    }
  } else {
    i = Knew->size[0] * Knew->size[1];
    Knew->size[0] = K->size[0];
    Knew->size[1] = K->size[1];
    emxEnsureCapacity_real_T(Knew, i);
    vlen = K->size[0] * K->size[1];
    for (i = 0; i < vlen; i++) {
      Knew->data[i] = K->data[i];
    }
  }

  emxFree_uint32_T(&indVec);
}

//
// File trailer for applyBCModal.cpp
//
// [EOF]
//
