//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: writeOutput.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "writeOutput.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <stdio.h>
#include <string.h>

// Function Definitions

//
// writeOutput writes output from modal analysis
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freqSorted,dampSorted,imagCompSignSorted] = writeOutput(freq,damp,...
//                                                 phase1,phase2,...
//                                                 imagComponentSign,fid)
//
//    This function writes an output file for modal analysis.
//
//       input:
//       freq               = array of modal frequencies
//       damp               = array of modal damping ratios
//       phase1             = array of in phase mode shapes
//       phase2             = array of out of phase mode shapes
//       imagComponentSign  = array of sign of imaginary components
//       fid                = file identifier for output
//
//       output:
//       freqSorted         = array of sorted(by frequency) modal frequencies
//       dampSorted         = array of sorted(by frequency) modal damping ratios
//       imagCompSignSorted = array of sorted(by frequency) of imaginarycomponentSign array
// Arguments    : const emxArray_real_T *freq
//                const emxArray_real_T *damp
//                const emxArray_real_T *phase1
//                const emxArray_real_T *phase2
//                const emxArray_real_T *imagComponentSign
//                double fid
//                emxArray_real_T *freqSorted
//                emxArray_real_T *dampSorted
//                emxArray_real_T *imagCompSignSorted
// Return Type  : void
//
void writeOutput(const emxArray_real_T *freq, const emxArray_real_T *damp, const
                 emxArray_real_T *phase1, const emxArray_real_T *phase2, const
                 emxArray_real_T *imagComponentSign, double fid, emxArray_real_T
                 *freqSorted, emxArray_real_T *dampSorted, emxArray_real_T
                 *imagCompSignSorted)
{
  emxArray_real_T *b_freq;
  int u0;
  int loop_ub;
  int u1;
  boolean_T swapped;
  emxArray_real_T *map;
  double b_index;
  FILE * b_NULL;
  FILE * c_NULL;
  FILE * d_NULL;
  FILE * e_NULL;
  FILE * f_NULL;
  FILE * filestar;
  int i;
  FILE * g_NULL;
  FILE * h_NULL;
  FILE * i_NULL;
  int i1;
  int b_u1;
  FILE * j_NULL;
  emxInit_real_T(&b_freq, 2);

  // gets number of nodes for mode shape printing
  // bubbleSort Sorts the vector A in ascending order
  //  **********************************************************************
  //  *                   Part of the SNL OWENS Toolkit                    *
  //  * Developed by Sandia National Laboratories Wind Energy Technologies *
  //  *             See license.txt for disclaimer information             *
  //  **********************************************************************
  //    [A,origMap,posIndex] = bubbleSort(A)
  //
  //    This function accepts a vector A, sorts the vector in ascending order,
  //    outputting the sorted vector, a map to the original ordering, and the
  //    index at which a positive value first occurs.
  //
  //       input:
  //       A        = vector to be sorted
  //
  //       output:
  //       A        = sorted vector
  //       origMap  = map of sorted vector to original ordering
  //       posIndex = index at which positive value first occurs
  u0 = b_freq->size[0] * b_freq->size[1];
  b_freq->size[0] = freq->size[0];
  b_freq->size[1] = freq->size[1];
  emxEnsureCapacity_real_T(b_freq, u0);
  loop_ub = freq->size[0] * freq->size[1];
  for (u0 = 0; u0 < loop_ub; u0++) {
    b_freq->data[u0] = freq->data[u0];
  }

  if ((freq->size[0] == 0) || (freq->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = freq->size[0];
    u1 = freq->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  swapped = true;
  emxInit_real_T(&map, 2);
  if (u1 < 1) {
    map->size[0] = 1;
    map->size[1] = 0;
  } else {
    u0 = map->size[0] * map->size[1];
    map->size[0] = 1;
    loop_ub = u1 - 1;
    map->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(map, u0);
    for (u0 = 0; u0 <= loop_ub; u0++) {
      map->data[u0] = static_cast<double>(u0) + 1.0;
    }
  }

  while (swapped) {
    swapped = false;
    for (loop_ub = 0; loop_ub <= u1 - 2; loop_ub++) {
      if (b_freq->data[loop_ub + 1] < b_freq->data[loop_ub]) {
        b_index = b_freq->data[loop_ub];
        b_freq->data[loop_ub] = b_freq->data[loop_ub + 1];
        b_freq->data[loop_ub + 1] = b_index;
        swapped = true;
        b_index = map->data[loop_ub];
        map->data[loop_ub] = map->data[loop_ub + 1];
        map->data[loop_ub + 1] = b_index;
      }
    }
  }

  //      posIndex = length(A)/2+1;
  //      [posIndex] = findPositiveCrossOver(A);
  // sorts frequency
  if ((b_freq->size[0] == 0) || (b_freq->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = b_freq->size[0];
    u1 = b_freq->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  u0 = dampSorted->size[0] * dampSorted->size[1];
  dampSorted->size[0] = 1;
  dampSorted->size[1] = u1;
  emxEnsureCapacity_real_T(dampSorted, u0);
  for (u0 = 0; u0 < u1; u0++) {
    dampSorted->data[u0] = 0.0;
  }

  if ((b_freq->size[0] == 0) || (b_freq->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = b_freq->size[0];
    u1 = b_freq->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  u0 = freqSorted->size[0] * freqSorted->size[1];
  freqSorted->size[0] = 1;
  freqSorted->size[1] = u1;
  emxEnsureCapacity_real_T(freqSorted, u0);
  for (u0 = 0; u0 < u1; u0++) {
    freqSorted->data[u0] = 0.0;
  }

  if ((b_freq->size[0] == 0) || (b_freq->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = b_freq->size[0];
    u1 = b_freq->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  u0 = imagCompSignSorted->size[0] * imagCompSignSorted->size[1];
  imagCompSignSorted->size[0] = 1;
  imagCompSignSorted->size[1] = u1;
  emxEnsureCapacity_real_T(imagCompSignSorted, u0);
  for (u0 = 0; u0 < u1; u0++) {
    imagCompSignSorted->data[u0] = 0.0;
  }

  b_index = 1.0;
  if ((b_freq->size[0] == 0) || (b_freq->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = b_freq->size[0];
    u1 = b_freq->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  if (0 <= u1 - 1) {
    b_NULL = NULL;
    c_NULL = NULL;
    d_NULL = NULL;
    e_NULL = NULL;
    f_NULL = NULL;
    i = phase1->size[0] - 1;
    g_NULL = NULL;
    h_NULL = NULL;
    i_NULL = NULL;
    i1 = phase1->size[0] - 1;
    if ((b_freq->size[0] == 0) || (b_freq->size[1] == 0)) {
      b_u1 = 0;
    } else {
      u0 = b_freq->size[0];
      b_u1 = b_freq->size[1];
      if (u0 > b_u1) {
        b_u1 = u0;
      }
    }
  }

  for (loop_ub = 0; loop_ub < u1; loop_ub++) {
    // prints mode frequency, damping and in/out of phase mode shapes
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "MODE # %0.0f \n\n", b_index);
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == c_NULL)) {
      fprintf(filestar, "Frequency: %e: \n", b_freq->data[loop_ub]);
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == d_NULL)) {
      fprintf(filestar, "Damping %e: \n", damp->data[static_cast<int>(map->
               data[loop_ub]) - 1]);
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == e_NULL)) {
      fprintf(filestar, "0 deg Mode Shape:\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == f_NULL)) {
      fprintf(filestar,
              "U_x          U_y          U_z          theta_x     theta_y     theta_z \n");
      if (swapped) {
        fflush(filestar);
      }
    }

    u0 = static_cast<int>(map->data[loop_ub]) - 1;
    dampSorted->data[loop_ub] = damp->data[u0];
    freqSorted->data[loop_ub] = b_freq->data[loop_ub];
    imagCompSignSorted->data[loop_ub] = imagComponentSign->data[u0];
    for (u0 = 0; u0 <= i; u0++) {
      j_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == j_NULL)) {
        fprintf(filestar, "%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",
                phase1->data[u0 + phase1->size[0] * 6 * (static_cast<int>
                 (map->data[loop_ub]) - 1)], phase1->data[(u0 + phase1->size[0])
                + phase1->size[0] * 6 * (static_cast<int>(map->data[loop_ub]) -
                 1)], phase1->data[(u0 + phase1->size[0] * 2) + phase1->size[0] *
                6 * (static_cast<int>(map->data[loop_ub]) - 1)], phase1->data
                [(u0 + phase1->size[0] * 3) + phase1->size[0] * 6 * (
                 static_cast<int>(map->data[loop_ub]) - 1)], phase1->data[(u0 +
                 phase1->size[0] * 4) + phase1->size[0] * 6 * (static_cast<int>
                 (map->data[loop_ub]) - 1)], phase1->data[(u0 + phase1->size[0] *
                 5) + phase1->size[0] * 6 * (static_cast<int>(map->data[loop_ub])
                 - 1)]);
        if (swapped) {
          fflush(filestar);
        }
      }

      j_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == j_NULL)) {
        fprintf(filestar, "\n");
        if (swapped) {
          fflush(filestar);
        }
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == g_NULL)) {
      fprintf(filestar, "\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == h_NULL)) {
      fprintf(filestar, "90 deg Mode Shape:\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == i_NULL)) {
      fprintf(filestar,
              "U_x          U_y          U_z          theta_x     theta_y     theta_z \n");
      if (swapped) {
        fflush(filestar);
      }
    }

    for (u0 = 0; u0 <= i1; u0++) {
      j_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == j_NULL)) {
        fprintf(filestar, "%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",
                phase2->data[u0 + phase2->size[0] * 6 * (static_cast<int>
                 (map->data[loop_ub]) - 1)], phase2->data[(u0 + phase2->size[0])
                + phase2->size[0] * 6 * (static_cast<int>(map->data[loop_ub]) -
                 1)], phase2->data[(u0 + phase2->size[0] * 2) + phase2->size[0] *
                6 * (static_cast<int>(map->data[loop_ub]) - 1)], phase2->data
                [(u0 + phase2->size[0] * 3) + phase2->size[0] * 6 * (
                 static_cast<int>(map->data[loop_ub]) - 1)], phase2->data[(u0 +
                 phase2->size[0] * 4) + phase2->size[0] * 6 * (static_cast<int>
                 (map->data[loop_ub]) - 1)], phase2->data[(u0 + phase2->size[0] *
                 5) + phase2->size[0] * 6 * (static_cast<int>(map->data[loop_ub])
                 - 1)]);
        if (swapped) {
          fflush(filestar);
        }
      }

      j_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == j_NULL)) {
        fprintf(filestar, "\n");
        if (swapped) {
          fflush(filestar);
        }
      }
    }

    if (loop_ub + 1 < b_u1) {
      j_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == j_NULL)) {
        fprintf(filestar, "\n\n");
        if (swapped) {
          fflush(filestar);
        }
      }
    }

    b_index++;
  }

  emxFree_real_T(&map);
  emxFree_real_T(&b_freq);
}

//
// File trailer for writeOutput.cpp
//
// [EOF]
//
