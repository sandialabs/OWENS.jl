//
// File: writeOutput.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
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
  int l;
  int i;
  int loop_ub;
  boolean_T swapped;
  emxArray_real_T *map;
  int b_i;
  double b_index;
  double d;
  FILE * b_NULL;
  FILE * filestar;
  emxInit_real_T(&b_freq, 2);
  l = phase1->size[0] - 1;

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
  i = b_freq->size[0] * b_freq->size[1];
  b_freq->size[0] = 1;
  b_freq->size[1] = freq->size[1];
  emxEnsureCapacity_real_T(b_freq, i);
  loop_ub = freq->size[0] * freq->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_freq->data[i] = freq->data[i];
  }

  swapped = true;
  emxInit_real_T(&map, 2);
  if (freq->size[1] < 1) {
    map->size[0] = 1;
    map->size[1] = 0;
  } else {
    i = map->size[0] * map->size[1];
    map->size[0] = 1;
    map->size[1] = static_cast<int>((static_cast<double>(freq->size[1]) - 1.0))
      + 1;
    emxEnsureCapacity_real_T(map, i);
    loop_ub = static_cast<int>((static_cast<double>(freq->size[1]) - 1.0));
    for (i = 0; i <= loop_ub; i++) {
      map->data[i] = static_cast<double>(i) + 1.0;
    }
  }

  while (swapped) {
    swapped = false;
    i = freq->size[1];
    for (b_i = 0; b_i <= i - 2; b_i++) {
      b_index = b_freq->data[b_i];
      d = b_freq->data[b_i + 1];
      if (d < b_index) {
        b_freq->data[b_i] = d;
        b_freq->data[b_i + 1] = b_index;
        swapped = true;
        b_index = map->data[b_i];
        map->data[b_i] = map->data[b_i + 1];
        map->data[b_i + 1] = b_index;
      }
    }
  }

  //      posIndex = length(A)/2+1;
  //      [posIndex] = findPositiveCrossOver(A);
  // sorts frequency
  i = dampSorted->size[0] * dampSorted->size[1];
  dampSorted->size[0] = 1;
  dampSorted->size[1] = b_freq->size[1];
  emxEnsureCapacity_real_T(dampSorted, i);
  loop_ub = b_freq->size[1];
  for (i = 0; i < loop_ub; i++) {
    dampSorted->data[i] = 0.0;
  }

  i = freqSorted->size[0] * freqSorted->size[1];
  freqSorted->size[0] = 1;
  freqSorted->size[1] = b_freq->size[1];
  emxEnsureCapacity_real_T(freqSorted, i);
  loop_ub = b_freq->size[1];
  for (i = 0; i < loop_ub; i++) {
    freqSorted->data[i] = 0.0;
  }

  i = imagCompSignSorted->size[0] * imagCompSignSorted->size[1];
  imagCompSignSorted->size[0] = 1;
  imagCompSignSorted->size[1] = b_freq->size[1];
  emxEnsureCapacity_real_T(imagCompSignSorted, i);
  loop_ub = b_freq->size[1];
  for (i = 0; i < loop_ub; i++) {
    imagCompSignSorted->data[i] = 0.0;
  }

  b_index = 1.0;
  i = b_freq->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    // prints mode frequency, damping and in/out of phase mode shapes
    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "MODE # %0.0f \n\n", b_index);
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "Frequency: %e: \n", b_freq->data[b_i]);
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "Damping %e: \n", damp->data[static_cast<int>(map->
               data[b_i]) - 1]);
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "0 deg Mode Shape:\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar,
              "U_x          U_y          U_z          theta_x     theta_y     theta_z \n");
      if (swapped) {
        fflush(filestar);
      }
    }

    loop_ub = static_cast<int>(map->data[b_i]) - 1;
    dampSorted->data[b_i] = damp->data[loop_ub];
    freqSorted->data[b_i] = b_freq->data[b_i];
    imagCompSignSorted->data[b_i] = imagComponentSign->data[loop_ub];
    for (loop_ub = 0; loop_ub <= l; loop_ub++) {
      b_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",
                phase1->data[loop_ub + phase1->size[0] * 6 * (static_cast<int>
                 (map->data[b_i]) - 1)], phase1->data[(loop_ub + phase1->size[0])
                + phase1->size[0] * 6 * (static_cast<int>(map->data[b_i]) - 1)],
                phase1->data[(loop_ub + phase1->size[0] * 2) + phase1->size[0] *
                6 * (static_cast<int>(map->data[b_i]) - 1)], phase1->data
                [(loop_ub + phase1->size[0] * 3) + phase1->size[0] * 6 * (
                 static_cast<int>(map->data[b_i]) - 1)], phase1->data[(loop_ub +
                 phase1->size[0] * 4) + phase1->size[0] * 6 * (static_cast<int>
                 (map->data[b_i]) - 1)], phase1->data[(loop_ub + phase1->size[0]
                 * 5) + phase1->size[0] * 6 * (static_cast<int>(map->data[b_i])
                 - 1)]);
        if (swapped) {
          fflush(filestar);
        }
      }

      b_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "\n");
        if (swapped) {
          fflush(filestar);
        }
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "90 deg Mode Shape:\n");
      if (swapped) {
        fflush(filestar);
      }
    }

    b_NULL = NULL;
    getfilestar(fid, &filestar, &swapped);
    if (!(filestar == b_NULL)) {
      fprintf(filestar,
              "U_x          U_y          U_z          theta_x     theta_y     theta_z \n");
      if (swapped) {
        fflush(filestar);
      }
    }

    for (loop_ub = 0; loop_ub <= l; loop_ub++) {
      b_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t%8.6f \t",
                phase2->data[loop_ub + phase2->size[0] * 6 * (static_cast<int>
                 (map->data[b_i]) - 1)], phase2->data[(loop_ub + phase2->size[0])
                + phase2->size[0] * 6 * (static_cast<int>(map->data[b_i]) - 1)],
                phase2->data[(loop_ub + phase2->size[0] * 2) + phase2->size[0] *
                6 * (static_cast<int>(map->data[b_i]) - 1)], phase2->data
                [(loop_ub + phase2->size[0] * 3) + phase2->size[0] * 6 * (
                 static_cast<int>(map->data[b_i]) - 1)], phase2->data[(loop_ub +
                 phase2->size[0] * 4) + phase2->size[0] * 6 * (static_cast<int>
                 (map->data[b_i]) - 1)], phase2->data[(loop_ub + phase2->size[0]
                 * 5) + phase2->size[0] * 6 * (static_cast<int>(map->data[b_i])
                 - 1)]);
        if (swapped) {
          fflush(filestar);
        }
      }

      b_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "\n");
        if (swapped) {
          fflush(filestar);
        }
      }
    }

    if (b_i + 1 < b_freq->size[1]) {
      b_NULL = NULL;
      getfilestar(fid, &filestar, &swapped);
      if (!(filestar == b_NULL)) {
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
