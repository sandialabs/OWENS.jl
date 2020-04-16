//
// File: readBladeData.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
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
  double d;
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
      *bladeData_numBlades = a_data[idx - 1];
      idx++;
      for (k = idx; k <= i; k++) {
        d = a_data[k - 1];
        if (*bladeData_numBlades < d) {
          *bladeData_numBlades = d;
        }
      }
    }
  }

  //      numStruts = min(bladeNum);
  //      if(numStruts>0)
  //          numStruts = 0;
  //      else
  //          numStruts = abs(numStruts);
  //      end
  idx = -1;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= a_size[0] - 1)) {
    if (rtIsNaN(a_data[k + a_size[0] * 15])) {
      idx = k;
      exitg1 = true;
    } else {
      k++;
    }
  }

  if (idx + 1 == 0) {
    idx = a_size[0];
  } else {
    //          strutDataBlock = a(strutStartIndex:end,:);
    //          [strutEntries, ~] = size(strutDataBlock);
    //          numNodesPerStrut = strutEntries/numStruts;
    //          numElPerStrut = numNodesPerStrut - 1;
  }

  if (1 > idx) {
    k = -1;
  } else {
    k = idx - 1;
  }

  // assign data to bladeData object %TODO: Should not be loading this file in multiple times 
  bladeData_bladeNum_size[0] = k + 1;
  bladeData_h_size[0] = k + 1;
  bladeData_nodeNum_size[0] = k + 1;
  bladeData_elementNum_size[0] = k + 1;
  if (0 <= k) {
    std::memcpy(&bladeData_bladeNum_data[0], &a_data[0], (k + 1) * sizeof(double));
  }

  for (i = 0; i <= k; i++) {
    bladeData_h_data[i] = a_data[i + a_size[0]];
    bladeData_nodeNum_data[i] = a_data[i + a_size[0] * 2];
    bladeData_elementNum_data[i] = a_data[i + a_size[0] * 3];
  }

  bladeData_remaining_size[0] = k + 1;
  bladeData_remaining_size[1] = 12;
  for (i = 0; i < 12; i++) {
    for (idx = 0; idx <= k; idx++) {
      bladeData_remaining_data[idx + bladeData_remaining_size[0] * i] =
        a_data[idx + a_size[0] * (i + 4)];
    }
  }
}

//
// File trailer for readBladeData.cpp
//
// [EOF]
//
