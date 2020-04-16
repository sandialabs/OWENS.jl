//
// File: applyBC.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "applyBC.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// applyBC Applies boundary conditions to system for static analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Kg,Fg] = applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode)
//
//    This function applies boundary conditions to the stiffness matrix and
//    load vector for a static analysis.
//
//       input:
//       Kg            = assembled global stiffness matrix
//       Fg            = assembled global load vector
//       BC            = struct of boundary condition information
//       u             = global displacement vector
//       iterationType = for nonlinear analysis, not used in BLAST
//       numDofPerNode = number of degrees of freedom per node
// Arguments    : emxArray_real_T *Kg
//                emxArray_real_T *Fg
//                double BC_numpBC
//                const emxArray_real_T *BC_pBC
// Return Type  : void
//
void applyBC(emxArray_real_T *Kg, emxArray_real_T *Fg, double BC_numpBC, const
             emxArray_real_T *BC_pBC)
{
  int numEq;
  int i;
  int b_i;
  double specVal;
  double eqNumber;
  int j;
  int i1;

  //       output:
  //       Kg            = global stiffness matrix with boundary conditions
  //       Fg            = global load vector with boundary condition
  numEq = Kg->size[0] - 1;

  // APPLY BCs FOR PRIMARY VARIABLE
  if (BC_numpBC > 0.0) {
    i = BC_pBC->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      specVal = BC_pBC->data[b_i + BC_pBC->size[0] * 2];
      eqNumber = (BC_pBC->data[b_i] - 1.0) * 6.0 + BC_pBC->data[b_i +
        BC_pBC->size[0]];
      for (j = 0; j <= numEq; j++) {
        i1 = static_cast<int>(eqNumber) - 1;
        Kg->data[i1 + Kg->size[0] * j] = 0.0;
        Fg->data[j] -= Kg->data[j + Kg->size[0] * i1] * specVal;
        Kg->data[j + Kg->size[0] * i1] = 0.0;
      }

      i1 = static_cast<int>(eqNumber) - 1;
      Fg->data[i1] = specVal;
      Kg->data[i1 + Kg->size[0] * i1] = 1.0;
    }
  }

  // APPLY BCs FOR SECONDARY VARIABLE
  //  if(BC.numsBC > 0) % This does not appear to be used
  //      sBC = BC.sBC;
  //      [numsBC,~] = size(sBC);
  //
  //      for i=1:numsBC
  //          nodeNumber = sBC(i,1);
  //          dofNumber = sBC(i,2);
  //          specVal =  sBC(i,3);
  //
  //          eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
  //
  //          Fg(eqNumber) = Fg(eqNumber) + specVal;
  //
  //      end
  //  end
}

//
// File trailer for applyBC.cpp
//
// [EOF]
//
