//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: sortLE.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "sortLE.h"
#include "relop.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_rtwutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Declarations
static double rt_atan2d_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(static_cast<double>(b_u0), static_cast<double>(b_u1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

//
// Arguments    : const emxArray_creal_T *v
//                int idx1
//                int idx2
// Return Type  : boolean_T
//
boolean_T sortLE(const emxArray_creal_T *v, int idx1, int idx2)
{
  boolean_T p;
  double ma;
  boolean_T SCALEA;
  double mb;
  boolean_T SCALEB;
  double x;
  double br;
  double absai;
  double absbi;
  double Ma;
  if (rtIsNaN(v->data[idx2 - 1].re) || rtIsNaN(v->data[idx2 - 1].im)) {
    p = true;
  } else if (rtIsNaN(v->data[idx1 - 1].re) || rtIsNaN(v->data[idx1 - 1].im)) {
    p = false;
  } else {
    ma = std::abs(v->data[idx1 - 1].re);
    if ((ma > 8.9884656743115785E+307) || (std::abs(v->data[idx1 - 1].im) >
         8.9884656743115785E+307)) {
      SCALEA = true;
    } else {
      SCALEA = false;
    }

    mb = std::abs(v->data[idx2 - 1].re);
    if ((mb > 8.9884656743115785E+307) || (std::abs(v->data[idx2 - 1].im) >
         8.9884656743115785E+307)) {
      SCALEB = true;
    } else {
      SCALEB = false;
    }

    if (SCALEA || SCALEB) {
      x = rt_hypotd_snf(v->data[idx1 - 1].re / 2.0, v->data[idx1 - 1].im / 2.0);
      br = rt_hypotd_snf(v->data[idx2 - 1].re / 2.0, v->data[idx2 - 1].im / 2.0);
    } else {
      x = rt_hypotd_snf(v->data[idx1 - 1].re, v->data[idx1 - 1].im);
      br = rt_hypotd_snf(v->data[idx2 - 1].re, v->data[idx2 - 1].im);
    }

    if (iseq(x, br)) {
      absai = std::abs(v->data[idx1 - 1].im);
      absbi = std::abs(v->data[idx2 - 1].im);
      if (ma > absai) {
        Ma = ma;
        ma = absai;
      } else {
        Ma = absai;
      }

      if (mb > absbi) {
        absai = mb;
        mb = absbi;
      } else {
        absai = absbi;
      }

      if (Ma > absai) {
        if (ma < mb) {
          x = Ma - absai;
          br = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + absai / 2.0) * (mb - ma);
        } else {
          x = Ma;
          br = absai;
        }
      } else if (Ma < absai) {
        if (ma > mb) {
          br = absai - Ma;
          x = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + absai / 2.0) * (ma - mb);
        } else {
          x = Ma;
          br = absai;
        }
      } else {
        x = ma;
        br = mb;
      }

      if (iseq(x, br)) {
        x = rt_atan2d_snf(v->data[idx1 - 1].im, v->data[idx1 - 1].re);
        br = rt_atan2d_snf(v->data[idx2 - 1].im, v->data[idx2 - 1].re);
        if (iseq(x, br)) {
          absai = v->data[idx1 - 1].re;
          absbi = v->data[idx1 - 1].im;
          br = v->data[idx2 - 1].re;
          Ma = v->data[idx2 - 1].im;
          if (x > 0.78539816339744828) {
            if (x > 2.3561944901923448) {
              absai = -absbi;
              br = -Ma;
            } else {
              absai = -absai;
              br = -br;
            }
          } else if (x > -0.78539816339744828) {
            absai = absbi;
            br = Ma;
          } else {
            if (!(x > -2.3561944901923448)) {
              absai = -absbi;
              br = -Ma;
            }
          }

          x = absai;
          if (iseq(absai, br)) {
            x = 0.0;
            br = 0.0;
          }
        }
      }
    }

    p = (x <= br);
  }

  return p;
}

//
// File trailer for sortLE.cpp
//
// [EOF]
//
