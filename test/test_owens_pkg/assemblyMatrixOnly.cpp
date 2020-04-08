//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: assemblyMatrixOnly.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "assemblyMatrixOnly.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// assemblyMatrixOnly Assembles element matrices into global sys of equations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Kg] = assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)
//
//    This function assembles the element matrix into the
//    global system of equations
//
//       input:
//       Ke             = element matrix
//       conn           = element connectivity
//       numNodesPerEl  = number of nodes per element
//       numDofPerNode  = number of degrees of freedom per node
//       Kg             = global system matrix
// Arguments    : const double Ke[144]
//                const double conn[2]
//                emxArray_real_T *Kg
// Return Type  : void
//
void assemblyMatrixOnly(const double Ke[144], const double conn[2],
  emxArray_real_T *Kg)
{
  int count;
  double dofList_data[12];
  int i;
  int b_i;
  int j;
  int i1;
  double tmp_data[144];
  int b_tmp_data[12];
  int i2;
  int c_tmp_data[12];

  //       output:
  //       Kg             = global system matrix with assembled element
  count = 0;
  std::memset(&dofList_data[0], 0, 12U * sizeof(double));
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 6; j++) {
      dofList_data[count] = (conn[i] - 1.0) * 6.0 + (static_cast<double>(j) +
        1.0);
      count++;
    }
  }

  for (b_i = 0; b_i < 12; b_i++) {
    for (i1 = 0; i1 < 12; i1++) {
      tmp_data[i1 + 12 * b_i] = Kg->data[(static_cast<int>(dofList_data[i1]) +
        Kg->size[0] * (static_cast<int>(dofList_data[b_i]) - 1)) - 1];
    }

    i1 = static_cast<int>(dofList_data[b_i]) - 1;
    b_tmp_data[b_i] = i1;
    c_tmp_data[b_i] = i1;
  }

  for (b_i = 0; b_i < 12; b_i++) {
    for (i1 = 0; i1 < 12; i1++) {
      i2 = i1 + 12 * b_i;
      Kg->data[b_tmp_data[i1] + Kg->size[0] * c_tmp_data[b_i]] = tmp_data[i2] +
        Ke[i2];
    }
  }

  //  numDOFPerEl = length(dofList);
  //  %Assemble element i into global system
  //          for j=1:numDOFPerEl
  //              J = dofList(j);
  //              for m=1:numDOFPerEl
  //                  M = dofList(m);
  //                  Kg(J,M) = Kg(J,M) + Ke(j,m);
  //              end
  //          end
}

//
// assemblyMatrixOnly Assembles element matrices into global sys of equations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Kg] = assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg)
//
//    This function assembles the element matrix into the
//    global system of equations
//
//       input:
//       Ke             = element matrix
//       conn           = element connectivity
//       numNodesPerEl  = number of nodes per element
//       numDofPerNode  = number of degrees of freedom per node
//       Kg             = global system matrix
// Arguments    : const double Ke_data[]
//                const int Ke_size[2]
//                const double conn[2]
//                emxArray_real_T *Kg
// Return Type  : void
//
void b_assemblyMatrixOnly(const double Ke_data[], const int Ke_size[2], const
  double conn[2], emxArray_real_T *Kg)
{
  int count;
  double dofList_data[12];
  int i;
  int b_i;
  int j;
  int i1;
  double tmp_data[144];
  int b_tmp_data[12];
  int c_tmp_data[12];

  //       output:
  //       Kg             = global system matrix with assembled element
  count = 0;
  std::memset(&dofList_data[0], 0, 12U * sizeof(double));
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 6; j++) {
      dofList_data[count] = (conn[i] - 1.0) * 6.0 + (static_cast<double>(j) +
        1.0);
      count++;
    }
  }

  for (b_i = 0; b_i < 12; b_i++) {
    for (i1 = 0; i1 < 12; i1++) {
      tmp_data[i1 + 12 * b_i] = Kg->data[(static_cast<int>(dofList_data[i1]) +
        Kg->size[0] * (static_cast<int>(dofList_data[b_i]) - 1)) - 1];
    }

    i1 = static_cast<int>(dofList_data[b_i]) - 1;
    b_tmp_data[b_i] = i1;
    c_tmp_data[b_i] = i1;
  }

  for (b_i = 0; b_i < 12; b_i++) {
    for (i1 = 0; i1 < 12; i1++) {
      Kg->data[b_tmp_data[i1] + Kg->size[0] * c_tmp_data[b_i]] = tmp_data[i1 +
        12 * b_i] + Ke_data[i1 + Ke_size[0] * b_i];
    }
  }

  //  numDOFPerEl = length(dofList);
  //  %Assemble element i into global system
  //          for j=1:numDOFPerEl
  //              J = dofList(j);
  //              for m=1:numDOFPerEl
  //                  M = dofList(m);
  //                  Kg(J,M) = Kg(J,M) + Ke(j,m);
  //              end
  //          end
}

//
// File trailer for assemblyMatrixOnly.cpp
//
// [EOF]
//
