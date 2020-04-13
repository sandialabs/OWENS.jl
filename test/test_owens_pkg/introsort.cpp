//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: introsort.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "introsort.h"
#include "heapsort.h"
#include "insertionsort.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Type Definitions
struct c_struct_T
{
  int xstart;
  int xend;
  int depth;
};

// Function Definitions

//
// Arguments    : emxArray_int32_T *x
//                int xstart
//                int xend
// Return Type  : void
//
void introsort(emxArray_int32_T *x, int xstart, int xend)
{
  int nsort;
  int pmax;
  int pmin;
  boolean_T exitg1;
  int MAXDEPTH;
  int pivot;
  c_struct_T frame;
  int pow2p;
  int i;
  c_struct_T st_d_data[120];
  int exitg2;
  if (xstart < xend) {
    nsort = (xend - xstart) + 1;
    if (nsort <= 32) {
      insertionsort(x, xstart, xend);
    } else {
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        pivot = (pmin + pmax) >> 1;
        pow2p = 1 << pivot;
        if (pow2p == nsort) {
          pmax = pivot;
          exitg1 = true;
        } else if (pow2p > nsort) {
          pmax = pivot;
        } else {
          pmin = pivot;
        }
      }

      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = xstart;
      frame.xend = xend;
      frame.depth = 0;
      nsort = MAXDEPTH << 1;
      for (i = 0; i < nsort; i++) {
        st_d_data[i] = frame;
      }

      st_d_data[0] = frame;
      pmax = 1;
      while (pmax > 0) {
        frame = st_d_data[pmax - 1];
        pmax--;
        i = frame.xend - frame.xstart;
        if (i + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend);
        } else {
          nsort = (frame.xstart + i / 2) - 1;
          pow2p = frame.xstart - 1;
          if (x->data[nsort] < x->data[pow2p]) {
            pmin = x->data[pow2p];
            x->data[pow2p] = x->data[nsort];
            x->data[nsort] = pmin;
          }

          i = frame.xend - 1;
          if (x->data[i] < x->data[pow2p]) {
            pmin = x->data[pow2p];
            x->data[pow2p] = x->data[i];
            x->data[i] = pmin;
          }

          if (x->data[i] < x->data[nsort]) {
            pmin = x->data[nsort];
            x->data[nsort] = x->data[i];
            x->data[i] = pmin;
          }

          pivot = x->data[nsort];
          i = frame.xend - 2;
          x->data[nsort] = x->data[i];
          x->data[i] = pivot;
          nsort = i;
          do {
            exitg2 = 0;
            for (pow2p++; x->data[pow2p] < pivot; pow2p++) {
            }

            for (nsort--; pivot < x->data[nsort]; nsort--) {
            }

            if (pow2p + 1 >= nsort + 1) {
              exitg2 = 1;
            } else {
              pmin = x->data[pow2p];
              x->data[pow2p] = x->data[nsort];
              x->data[nsort] = pmin;
            }
          } while (exitg2 == 0);

          x->data[i] = x->data[pow2p];
          x->data[pow2p] = pivot;
          if (pow2p + 2 < frame.xend) {
            st_d_data[pmax].xstart = pow2p + 2;
            st_d_data[pmax].xend = frame.xend;
            st_d_data[pmax].depth = frame.depth + 1;
            pmax++;
          }

          if (frame.xstart < pow2p + 1) {
            st_d_data[pmax].xstart = frame.xstart;
            st_d_data[pmax].xend = pow2p + 1;
            st_d_data[pmax].depth = frame.depth + 1;
            pmax++;
          }
        }
      }
    }
  }
}

//
// File trailer for introsort.cpp
//
// [EOF]
//
