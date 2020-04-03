//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: assembly.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 03-Apr-2020 15:56:19
//

// Include Files
#include "assembly.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// assembly Assembles element matrices into global system of equations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Kg,Fg] = assembly(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg)
//
//    This function assembles the element matrix and load vector into the
//    global system of equations
//
//       input:
//       Ke             = element matrix
//       Fe             = element vector
//       conn           = element connectivity
//       numNodesPerEl  = number of nodes per element
//       numDofPerNode  = number of degrees of freedom per node
//       Kg             = global system matrix
//       Fg             = global load vector
// Arguments    : const double Ke[144]
//                const double Fe[12]
//                const double conn[2]
//                emxArray_real_T *Kg
//                emxArray_real_T *Fg
// Return Type  : void
//
void assembly(const double Ke[144], const double Fe[12], const double conn[2],
              emxArray_real_T *Kg, emxArray_real_T *Fg)
{
  int count;
  double dofList_data[12];
  int i;
  int j;
  int b_i;

  //       output:
  //       Kg             = global system matrix with assembled element
  //       Fg             = global load vector with assembled element
  count = 0;
  std::memset(&dofList_data[0], 0, 12U * sizeof(double));
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 6; j++) {
      dofList_data[count] = (conn[i] - 1.0) * 6.0 + (static_cast<double>(j) +
        1.0);
      count++;
    }
  }

  // Assemble element i into global system
  for (j = 0; j < 12; j++) {
    count = static_cast<int>(dofList_data[j]) - 1;
    Fg->data[count] += Fe[j];
    for (i = 0; i < 12; i++) {
      b_i = static_cast<int>(dofList_data[i]) - 1;
      Kg->data[count + Kg->size[0] * b_i] += Ke[j + 12 * i];
    }
  }
}

//
// File trailer for assembly.cpp
//
// [EOF]
//
