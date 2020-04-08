//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readBladeData.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "readBladeData.h"
#include "importCactusFile.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// readBladeDAta reads blade data
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [bladeData] = readBladeData(filename)
//
//    This function reads blade data from file
//
//    input:
//    filename      = string containing file name for blade data file
//
//    output:
//    bladeData     = object containing blade data
// Arguments    : const emxArray_char_T *filename
//                double *bladeData_numBlades
//                double bladeData_bladeNum_data[]
//                int bladeData_bladeNum_size[1]
//                double bladeData_h_data[]
//                int bladeData_h_size[1]
//                double bladeData_nodeNum_data[]
//                int bladeData_nodeNum_size[1]
//                double bladeData_elementNum_data[]
//                int bladeData_elementNum_size[1]
//                double bladeData_remaining_data[]
//                int bladeData_remaining_size[2]
// Return Type  : void
//
void readBladeData(const emxArray_char_T *filename, double *bladeData_numBlades,
                   double bladeData_bladeNum_data[], int
                   bladeData_bladeNum_size[1], double bladeData_h_data[], int
                   bladeData_h_size[1], double bladeData_nodeNum_data[], int
                   bladeData_nodeNum_size[1], double bladeData_elementNum_data[],
                   int bladeData_elementNum_size[1], double
                   bladeData_remaining_data[], int bladeData_remaining_size[2])
{
  double a_data[960];
  int a_size[2];
  int i;
  int idx;
  int k;
  boolean_T exitg1;
  double ex;
  int i1;
  double d;
  int strutStartIndex;
  int b_i;
  importCactusFile(filename, a_data, a_size);
  i = a_size[0];
  if (a_size[0] <= 2) {
    if (a_size[0] == 1) {
      *bladeData_numBlades = a_data[0];
    } else if ((a_data[0] < a_data[1]) || (rtIsNaN(a_data[0]) && (!rtIsNaN
                 (a_data[1])))) {
      *bladeData_numBlades = a_data[1];
    } else {
      *bladeData_numBlades = a_data[0];
    }
  } else {
    if (!rtIsNaN(a_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= a_size[0])) {
        if (!rtIsNaN(a_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      *bladeData_numBlades = a_data[0];
    } else {
      ex = a_data[idx - 1];
      i1 = idx + 1;
      for (k = i1; k <= i; k++) {
        d = a_data[k - 1];
        if (ex < d) {
          ex = d;
        }
      }

      *bladeData_numBlades = ex;
    }
  }

  //      numStruts = min(bladeNum);
  //      if(numStruts>0)
  //          numStruts = 0;
  //      else
  //          numStruts = abs(numStruts);
  //      end
  strutStartIndex = -1;
  b_i = 0;
  exitg1 = false;
  while ((!exitg1) && (b_i <= a_size[0] - 1)) {
    if (rtIsNaN(a_data[b_i + a_size[0] * 15])) {
      strutStartIndex = b_i;
      exitg1 = true;
    } else {
      b_i++;
    }
  }

  if (strutStartIndex + 1 == 0) {
    strutStartIndex = a_size[0];
  } else {
    //          strutDataBlock = a(strutStartIndex:end,:);
    //          [strutEntries, ~] = size(strutDataBlock);
    //          numNodesPerStrut = strutEntries/numStruts;
    //          numElPerStrut = numNodesPerStrut - 1;
  }

  if (1 > strutStartIndex) {
    idx = -1;
  } else {
    idx = strutStartIndex - 1;
  }

  // assign data to bladeData object %TODO: Should not be loading this file in multiple times 
  bladeData_bladeNum_size[0] = idx + 1;
  bladeData_h_size[0] = idx + 1;
  bladeData_nodeNum_size[0] = idx + 1;
  bladeData_elementNum_size[0] = idx + 1;
  if (0 <= idx) {
    std::memcpy(&bladeData_bladeNum_data[0], &a_data[0], (idx + 1) * sizeof
                (double));
  }

  for (i = 0; i <= idx; i++) {
    bladeData_h_data[i] = a_data[i + a_size[0]];
    bladeData_nodeNum_data[i] = a_data[i + a_size[0] * 2];
    bladeData_elementNum_data[i] = a_data[i + a_size[0] * 3];
  }

  bladeData_remaining_size[0] = idx + 1;
  bladeData_remaining_size[1] = 12;
  for (i = 0; i < 12; i++) {
    for (i1 = 0; i1 <= idx; i1++) {
      bladeData_remaining_data[i1 + bladeData_remaining_size[0] * i] = a_data[i1
        + a_size[0] * (i + 4)];
    }
  }
}

//
// File trailer for readBladeData.cpp
//
// [EOF]
//
