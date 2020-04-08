//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: externalForcing.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "externalForcing.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// owens externalForcing function for the OWENS toolkit
//  **********************************************************************
//  *                   Part of the SNL OWENS toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Fexternal, Fdof] = externalForcing(time,aeroLoads)
//
//    This function specifies external forcing for a transient analysis.
//    Fexternal is a vector of loads and Fdof is a corresponding vector of
//    degrees of freedom the concentrated loads in Fexternal correspond to.
//    The input time allows for arbitrary time varying loads
//    The global degree of freedom number corresponding with the local degree
//    of freedom of a node may be calculated by:
//    globalDOFNumber = (nodeNumber-1)*6 + localDOFnumber
//    The localDOFnumber may range from 1 to 6 such that 1 corresponds to a
//    force in "x direction" of the co-rotating hub frame. 2 and 3
//    corresponds to a force in the "y" and "z directions" respectively. 4,
//    5, and 6 correspond to a moment about the "x", "y", and "z" directions
//    respectively.
// Arguments    : double time
//                const double aeroLoads_timeArray_data[]
//                const int aeroLoads_timeArray_size[1]
//                const emxArray_real_T *aeroLoads_ForceValHist
//                const emxArray_real_T *aeroLoads_ForceDof
//                emxArray_real_T *Fexternal
//                emxArray_real_T *Fdof
// Return Type  : void
//
void externalForcing(double time, const double aeroLoads_timeArray_data[], const
                     int aeroLoads_timeArray_size[1], const emxArray_real_T
                     *aeroLoads_ForceValHist, const emxArray_real_T
                     *aeroLoads_ForceDof, emxArray_real_T *Fexternal,
                     emxArray_real_T *Fdof)
{
  int high_i;
  int n;
  double varargin_1_data[2002];
  int i;
  emxArray_real_T *varargin_2;
  int nyrows;
  int vlen;
  int nycols;
  int i1;
  int nx;
  int k;
  int exitg1;
  int j;
  double xtmp;
  int offset;
  int tmp_tmp;
  double b_y1;

  //
  //       input:
  //       time         = simulation time
  //
  //       output:
  //       Fexternal     = vector of external loads (forces/moments)
  //       Fdof          = vector of corresponding DOF numbers to apply loads to 
  //      if(time < 0.2)
  //          Fexternal = 1e6;
  //          Fdof = 20*6+1;
  //      else
  //          Fexternal = [];
  //          Fdof = [];
  //      end
  // temp = load('aeroLoads.mat');
  high_i = aeroLoads_timeArray_size[0];
  n = aeroLoads_timeArray_size[0];
  if (0 <= n - 1) {
    std::memcpy(&varargin_1_data[0], &aeroLoads_timeArray_data[0], n * sizeof
                (double));
  }

  i = Fdof->size[0];
  Fdof->size[0] = aeroLoads_ForceDof->size[0];
  emxEnsureCapacity_real_T(Fdof, i);
  n = aeroLoads_ForceDof->size[0];
  for (i = 0; i < n; i++) {
    Fdof->data[i] = aeroLoads_ForceDof->data[i];
  }

  emxInit_real_T(&varargin_2, 2);
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = aeroLoads_ForceValHist->size[1];
  varargin_2->size[1] = aeroLoads_ForceValHist->size[0];
  emxEnsureCapacity_real_T(varargin_2, i);
  n = aeroLoads_ForceValHist->size[0];
  for (i = 0; i < n; i++) {
    vlen = aeroLoads_ForceValHist->size[1];
    for (i1 = 0; i1 < vlen; i1++) {
      varargin_2->data[i1 + varargin_2->size[0] * i] =
        aeroLoads_ForceValHist->data[i + aeroLoads_ForceValHist->size[0] * i1];
    }
  }

  nyrows = varargin_2->size[0];
  nycols = varargin_2->size[1] - 1;
  nx = aeroLoads_timeArray_size[0] - 1;
  i = Fexternal->size[0] * Fexternal->size[1];
  Fexternal->size[0] = 1;
  Fexternal->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(Fexternal, i);
  n = varargin_2->size[1];
  for (i = 0; i < n; i++) {
    Fexternal->data[i] = rtNaN;
  }

  k = 0;
  do {
    exitg1 = 0;
    if (k <= nx) {
      if (rtIsNaN(aeroLoads_timeArray_data[k])) {
        exitg1 = 1;
      } else {
        k++;
      }
    } else {
      if (aeroLoads_timeArray_data[1] < aeroLoads_timeArray_data[0]) {
        i = (nx + 1) >> 1;
        for (vlen = 0; vlen < i; vlen++) {
          xtmp = varargin_1_data[vlen];
          n = nx - vlen;
          varargin_1_data[vlen] = varargin_1_data[n];
          varargin_1_data[n] = xtmp;
        }

        if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0) &&
            (varargin_2->size[0] > 1)) {
          vlen = varargin_2->size[0];
          n = varargin_2->size[0] - 1;
          nx = varargin_2->size[0] >> 1;
          i = varargin_2->size[1] - 1;
          for (j = 0; j <= i; j++) {
            offset = j * vlen;
            for (k = 0; k < nx; k++) {
              tmp_tmp = offset + k;
              xtmp = varargin_2->data[tmp_tmp];
              i1 = (offset + n) - k;
              varargin_2->data[tmp_tmp] = varargin_2->data[i1];
              varargin_2->data[i1] = xtmp;
            }
          }
        }
      }

      if (rtIsNaN(time)) {
        for (j = 0; j <= nycols; j++) {
          Fexternal->data[j] = rtNaN;
        }
      } else {
        if ((!(time > varargin_1_data[high_i - 1])) && (!(time <
              varargin_1_data[0]))) {
          nx = 1;
          vlen = 2;
          while (high_i > vlen) {
            n = (nx >> 1) + (high_i >> 1);
            if (((nx & 1) == 1) && ((high_i & 1) == 1)) {
              n++;
            }

            if (time >= varargin_1_data[n - 1]) {
              nx = n;
              vlen = n + 1;
            } else {
              high_i = n;
            }
          }

          xtmp = varargin_1_data[nx - 1];
          xtmp = (time - xtmp) / (varargin_1_data[nx] - xtmp);
          if (xtmp == 0.0) {
            for (j = 0; j <= nycols; j++) {
              Fexternal->data[j] = varargin_2->data[(nx + j * nyrows) - 1];
            }
          } else if (xtmp == 1.0) {
            for (j = 0; j <= nycols; j++) {
              Fexternal->data[j] = varargin_2->data[nx + j * nyrows];
            }
          } else {
            for (j = 0; j <= nycols; j++) {
              vlen = nx + j * nyrows;
              b_y1 = varargin_2->data[vlen - 1];
              if (b_y1 == varargin_2->data[vlen]) {
                Fexternal->data[j] = b_y1;
              } else {
                Fexternal->data[j] = (1.0 - xtmp) * b_y1 + xtmp *
                  varargin_2->data[vlen];
              }
            }
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_real_T(&varargin_2);
}

//
// File trailer for externalForcing.cpp
//
// [EOF]
//
