//
// File: calculateTheo.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "calculateTheo.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// calculateTheo Calculates Theodorsen function value for unsteady aero
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Theo] = calculateTheo(k)
//
//    This function accepts a reduced frequency value and calculates the
//    complex Theodorsen function value.
//
//    input:
//    k     = reduced frequency
//
//    output:
//    Theo  = Theodorsen  constant value
// Arguments    : double k
// Return Type  : creal_T
//
creal_T calculateTheo(double k)
{
  creal_T Theo;
  double br;
  double bi;
  double bim;
  double re;
  double im;
  double s;
  br = 1.0 - 0.0455 / k * 0.0;
  bi = 0.0 - 0.0455 / k;
  if (bi == 0.0) {
    re = 0.165 / br;
    im = 0.0;
  } else {
    bim = std::abs(bi);
    if (br > bim) {
      s = bi / br;
      bim = br + s * bi;
      re = (0.165 + s * 0.0) / bim;
      im = (0.0 - s * 0.165) / bim;
    } else if (bim == br) {
      if (br > 0.0) {
        br = 0.5;
      } else {
        br = -0.5;
      }

      if (bi > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      re = 0.165 * br + 0.0 * bim;
      im = 0.0 * br - 0.165 * bim;
    } else {
      s = br / bi;
      bim = bi + s * br;
      re = s * 0.165 / bim;
      im = (s * 0.0 - 0.165) / bim;
    }
  }

  br = 1.0 - 0.3 / k * 0.0;
  bi = 0.0 - 0.3 / k;
  if (bi == 0.0) {
    bi = 0.335 / br;
    bim = 0.0;
  } else {
    bim = std::abs(bi);
    if (br > bim) {
      s = bi / br;
      bim = br + s * bi;
      bi = (0.335 + s * 0.0) / bim;
      bim = (0.0 - s * 0.335) / bim;
    } else if (bim == br) {
      if (br > 0.0) {
        br = 0.5;
      } else {
        br = -0.5;
      }

      if (bi > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      bi = 0.335 * br + 0.0 * bim;
      bim = 0.0 * br - 0.335 * bim;
    } else {
      s = br / bi;
      bim = bi + s * br;
      bi = s * 0.335 / bim;
      bim = (s * 0.0 - 0.335) / bim;
    }
  }

  Theo.re = (1.0 - re) - bi;
  Theo.im = (0.0 - im) - bim;
  if (rtIsNaN(k) || rtIsInf(k)) {
    Theo.re = 0.0;
    Theo.im = 0.0;
  }

  return Theo;
}

//
// File trailer for calculateTheo.cpp
//
// [EOF]
//
