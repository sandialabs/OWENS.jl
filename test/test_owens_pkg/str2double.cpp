//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: str2double.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "str2double.h"
#include "rt_nonfinite.h"
#include "str2double1.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <stdio.h>
#include <string.h>

// Function Definitions

//
// Arguments    : char s
// Return Type  : creal_T
//
creal_T b_str2double(char s)
{
  creal_T x;
  int ntoread;
  int k;
  int i;
  boolean_T isimag1;
  boolean_T b_finite;
  double scanned1;
  int idx;
  char s1[3];
  boolean_T isneg;
  boolean_T exitg1;
  boolean_T success;
  double scanned2;
  boolean_T isfinite2;
  boolean_T foundsign;
  double b_scanned1;
  x.re = rtNaN;
  x.im = 0.0;
  ntoread = 0;
  k = 1;
  i = static_cast<unsigned char>(s) & 127;
  if (bv[i] || (s == '\x00')) {
    k = 2;
  }

  isimag1 = false;
  b_finite = true;
  scanned1 = 0.0;
  idx = 1;
  s1[0] = '\x00';
  s1[1] = '\x00';
  s1[2] = '\x00';
  isneg = false;
  exitg1 = false;
  while ((!exitg1) && (k <= 1)) {
    if (s == '-') {
      isneg = !isneg;
      k = 2;
    } else if ((s == ',') || (s == '+') || bv[i]) {
      k = 2;
    } else {
      exitg1 = true;
    }
  }

  success = (k <= 1);
  if (success && isneg) {
    s1[0] = '-';
    idx = 2;
  }

  if (success) {
    isneg = false;
    if (k <= 1) {
      isneg = ((s == 'j') || (s == 'i'));
    }

    if (isneg) {
      isimag1 = true;
      k++;
      s1[idx - 1] = '1';
      idx++;
    } else {
      b_readNonFinite(s, &k, &b_finite, &scanned1);
      if (b_finite) {
        success = b_copydigits(s1, &idx, s, &k, true);
        if (success) {
          success = b_copyexponent(s1, &idx, s, &k);
        }
      } else {
        if ((idx >= 2) && (s1[0] == '-')) {
          idx = 1;
          s1[0] = ' ';
          scanned1 = -scanned1;
        }
      }

      while ((k <= 1) && (bv[i] || (s == '\x00') || (s == ','))) {
        k = 2;
      }

      if ((k <= 1) && (s == '*')) {
        k = 2;
      }

      if ((k <= 1) && ((s == 'i') || (s == 'j'))) {
        k = 2;
        isimag1 = true;
      }
    }

    while ((k <= 1) && (bv[i] || (s == '\x00') || (s == ','))) {
      k = 2;
    }
  }

  if (b_finite) {
    ntoread = 1;
  }

  if (success && (k <= 1)) {
    s1[idx - 1] = ' ';
    idx++;
    k = 1;
    b_readfloat(s1, &idx, s, &k, true, &isneg, &isfinite2, &scanned2, &foundsign,
                &success);
    if (isfinite2) {
      ntoread++;
    }

    if (success && (k > 1) && ((isimag1 ^ isneg) != 0) && foundsign) {
      success = true;
    } else {
      success = false;
    }
  } else {
    scanned2 = 0.0;
  }

  if (success) {
    s1[idx - 1] = '\x00';
    if (ntoread == 2) {
      ntoread = sscanf(&s1[0], "%lf %lf", &scanned1, &scanned2);
      if (ntoread != 2) {
        scanned1 = rtNaN;
        scanned2 = rtNaN;
      }
    } else {
      if (ntoread == 1) {
        ntoread = sscanf(&s1[0], "%lf", &b_scanned1);
        if (ntoread != 1) {
          b_scanned1 = rtNaN;
        }

        if (b_finite) {
          scanned1 = b_scanned1;
        } else {
          scanned2 = b_scanned1;
        }
      }
    }

    if (isimag1) {
      x.re = scanned2;
      x.im = scanned1;
    } else {
      x.re = scanned1;
      x.im = scanned2;
    }
  }

  return x;
}

//
// Arguments    : const emxArray_char_T *s
// Return Type  : creal_T
//
creal_T c_str2double(const emxArray_char_T *s)
{
  creal_T x;
  emxArray_char_T *s1;
  int ntoread;
  int k;
  boolean_T exitg1;
  int idx;
  char c;
  int loop_ub;
  boolean_T isimag1;
  boolean_T isfinite1;
  double scanned1;
  boolean_T unusedU0;
  boolean_T success;
  double scanned2;
  boolean_T isfinite2;
  boolean_T foundsign;
  double b_scanned1;
  x.re = rtNaN;
  x.im = 0.0;
  if (s->size[1] >= 1) {
    emxInit_char_T(&s1, 2);
    ntoread = 0;
    k = 1;
    exitg1 = false;
    while ((!exitg1) && (k <= s->size[1])) {
      c = s->data[k - 1];
      if (bv[static_cast<unsigned char>(c) & 127] || (c == '\x00')) {
        k++;
      } else {
        exitg1 = true;
      }
    }

    idx = s1->size[0] * s1->size[1];
    s1->size[0] = 1;
    s1->size[1] = s->size[1] + 2;
    emxEnsureCapacity_char_T(s1, idx);
    loop_ub = s->size[1] + 2;
    for (idx = 0; idx < loop_ub; idx++) {
      s1->data[idx] = '\x00';
    }

    idx = 1;
    c_readfloat(s1, &idx, s, &k, s->size[1], true, &isimag1, &isfinite1,
                &scanned1, &unusedU0, &success);
    if (isfinite1) {
      ntoread = 1;
    }

    if (success && (k <= s->size[1])) {
      s1->data[idx - 1] = ' ';
      idx++;
      c_readfloat(s1, &idx, s, &k, s->size[1], true, &unusedU0, &isfinite2,
                  &scanned2, &foundsign, &success);
      if (isfinite2) {
        ntoread++;
      }

      if (success && (k > s->size[1]) && ((isimag1 ^ unusedU0) != 0) &&
          foundsign) {
        success = true;
      } else {
        success = false;
      }
    } else {
      scanned2 = 0.0;
    }

    if (success) {
      s1->data[idx - 1] = '\x00';
      if (ntoread == 2) {
        idx = sscanf(&s1->data[0], "%lf %lf", &scanned1, &scanned2);
        if (idx != 2) {
          scanned1 = rtNaN;
          scanned2 = rtNaN;
        }
      } else {
        if (ntoread == 1) {
          idx = sscanf(&s1->data[0], "%lf", &b_scanned1);
          if (idx != 1) {
            b_scanned1 = rtNaN;
          }

          if (isfinite1) {
            scanned1 = b_scanned1;
          } else {
            scanned2 = b_scanned1;
          }
        }
      }

      if (isimag1) {
        x.re = scanned2;
        x.im = scanned1;
      } else {
        x.re = scanned1;
        x.im = scanned2;
      }
    }

    emxFree_char_T(&s1);
  }

  return x;
}

//
// Arguments    : const emxArray_char_T *s
// Return Type  : creal_T
//
creal_T d_str2double(const emxArray_char_T *s)
{
  creal_T x;
  int n;
  emxArray_char_T *s1;
  int ntoread;
  int k;
  int idx;
  int loop_ub;
  boolean_T isimag1;
  boolean_T isfinite1;
  double scanned1;
  boolean_T unusedU0;
  boolean_T success;
  double scanned2;
  boolean_T isfinite2;
  boolean_T foundsign;
  double b_scanned1;
  x.re = rtNaN;
  x.im = 0.0;
  n = s->size[0] * s->size[1];
  if (n >= 1) {
    emxInit_char_T(&s1, 2);
    ntoread = 0;
    k = 1;
    while ((k <= n) && (bv[static_cast<unsigned char>(s->data[k - 1]) & 127] ||
                        (s->data[k - 1] == '\x00'))) {
      k++;
    }

    idx = s1->size[0] * s1->size[1];
    s1->size[0] = 1;
    s1->size[1] = n + 2;
    emxEnsureCapacity_char_T(s1, idx);
    loop_ub = n + 2;
    for (idx = 0; idx < loop_ub; idx++) {
      s1->data[idx] = '\x00';
    }

    idx = 1;
    d_readfloat(s1, &idx, s, &k, n, true, &isimag1, &isfinite1, &scanned1,
                &unusedU0, &success);
    if (isfinite1) {
      ntoread = 1;
    }

    if (success && (k <= n)) {
      s1->data[idx - 1] = ' ';
      idx++;
      d_readfloat(s1, &idx, s, &k, n, true, &unusedU0, &isfinite2, &scanned2,
                  &foundsign, &success);
      if (isfinite2) {
        ntoread++;
      }

      if (success && (k > n) && ((isimag1 ^ unusedU0) != 0) && foundsign) {
        success = true;
      } else {
        success = false;
      }
    } else {
      scanned2 = 0.0;
    }

    if (success) {
      s1->data[idx - 1] = '\x00';
      if (ntoread == 2) {
        idx = sscanf(&s1->data[0], "%lf %lf", &scanned1, &scanned2);
        if (idx != 2) {
          scanned1 = rtNaN;
          scanned2 = rtNaN;
        }
      } else {
        if (ntoread == 1) {
          idx = sscanf(&s1->data[0], "%lf", &b_scanned1);
          if (idx != 1) {
            b_scanned1 = rtNaN;
          }

          if (isfinite1) {
            scanned1 = b_scanned1;
          } else {
            scanned2 = b_scanned1;
          }
        }
      }

      if (isimag1) {
        x.re = scanned2;
        x.im = scanned1;
      } else {
        x.re = scanned1;
        x.im = scanned2;
      }
    }

    emxFree_char_T(&s1);
  }

  return x;
}

//
// Arguments    : const char s[2]
// Return Type  : creal_T
//
creal_T str2double(const char s[2])
{
  creal_T x;
  int ntoread;
  int k;
  boolean_T isimag1;
  boolean_T b_finite;
  double scanned1;
  int idx;
  char s1[4];
  boolean_T isneg;
  boolean_T exitg1;
  boolean_T success;
  int b_k;
  double scanned2;
  boolean_T unusedU2;
  boolean_T foundsign;
  double b_scanned1;
  char c;
  x.re = rtNaN;
  x.im = 0.0;
  ntoread = 0;
  k = 0;
  while ((k + 1 <= 2) && (bv[static_cast<unsigned char>(s[k]) & 127] || (s[k] ==
           '\x00'))) {
    k++;
  }

  isimag1 = false;
  b_finite = true;
  scanned1 = 0.0;
  idx = 1;
  s1[0] = '\x00';
  s1[1] = '\x00';
  s1[2] = '\x00';
  s1[3] = '\x00';
  isneg = false;
  exitg1 = false;
  while ((!exitg1) && (k + 1 <= 2)) {
    if (s[k] == '-') {
      isneg = !isneg;
      k++;
    } else if ((s[k] == ',') || (s[k] == '+') || bv[static_cast<unsigned char>
               (s[k]) & 127]) {
      k++;
    } else {
      exitg1 = true;
    }
  }

  success = (k + 1 <= 2);
  if (success && isneg) {
    s1[0] = '-';
    idx = 2;
  }

  b_k = k + 1;
  if (success) {
    isneg = false;
    if (k + 1 <= 2) {
      isneg = ((s[k] == 'j') || (s[k] == 'i'));
    }

    if (isneg) {
      isimag1 = true;
      b_k = k + 2;
      while ((b_k <= 2) && (bv[static_cast<unsigned char>(s[1]) & 127] || (s[1] ==
               '\x00') || (s[1] == ','))) {
        b_k = 3;
      }

      if ((b_k <= 2) && (s[1] == '*')) {
        b_k = 3;
        readfloat(s1, &idx, s, &b_k, false, &isneg, &b_finite, &scanned1,
                  &unusedU2, &success);
      } else {
        s1[idx - 1] = '1';
        idx++;
      }
    } else {
      b_k = k + 1;
      readNonFinite(s, &b_k, &b_finite, &scanned1);
      if (b_finite) {
        success = copydigits(s1, &idx, s, &b_k, true);
        if (success) {
          success = copyexponent(s1, &idx, s, &b_k);
        }
      } else {
        if ((idx >= 2) && (s1[0] == '-')) {
          idx = 1;
          s1[0] = ' ';
          scanned1 = -scanned1;
        }
      }

      exitg1 = false;
      while ((!exitg1) && (b_k <= 2)) {
        c = s[b_k - 1];
        if (bv[static_cast<unsigned char>(c) & 127] || (c == '\x00') || (c ==
             ',')) {
          b_k++;
        } else {
          exitg1 = true;
        }
      }

      if ((b_k <= 2) && (s[b_k - 1] == '*')) {
        b_k++;
        exitg1 = false;
        while ((!exitg1) && (b_k <= 2)) {
          c = s[b_k - 1];
          if (bv[static_cast<unsigned char>(c) & 127] || (c == '\x00') || (c ==
               ',')) {
            b_k++;
          } else {
            exitg1 = true;
          }
        }
      }

      if (b_k <= 2) {
        c = s[b_k - 1];
        if ((c == 'i') || (c == 'j')) {
          b_k++;
          isimag1 = true;
        }
      }
    }

    exitg1 = false;
    while ((!exitg1) && (b_k <= 2)) {
      c = s[b_k - 1];
      if (bv[static_cast<unsigned char>(c) & 127] || (c == '\x00') || (c == ','))
      {
        b_k++;
      } else {
        exitg1 = true;
      }
    }
  }

  if (b_finite) {
    ntoread = 1;
  }

  if (success && (b_k <= 2)) {
    s1[idx - 1] = ' ';
    idx++;
    readfloat(s1, &idx, s, &b_k, true, &isneg, &unusedU2, &scanned2, &foundsign,
              &success);
    if (unusedU2) {
      ntoread++;
    }

    if (success && (b_k > 2) && ((isimag1 ^ isneg) != 0) && foundsign) {
      success = true;
    } else {
      success = false;
    }
  } else {
    scanned2 = 0.0;
  }

  if (success) {
    s1[idx - 1] = '\x00';
    if (ntoread == 2) {
      ntoread = sscanf(&s1[0], "%lf %lf", &scanned1, &scanned2);
      if (ntoread != 2) {
        scanned1 = rtNaN;
        scanned2 = rtNaN;
      }
    } else {
      if (ntoread == 1) {
        ntoread = sscanf(&s1[0], "%lf", &b_scanned1);
        if (ntoread != 1) {
          b_scanned1 = rtNaN;
        }

        if (b_finite) {
          scanned1 = b_scanned1;
        } else {
          scanned2 = b_scanned1;
        }
      }
    }

    if (isimag1) {
      x.re = scanned2;
      x.im = scanned1;
    } else {
      x.re = scanned1;
      x.im = scanned2;
    }
  }

  return x;
}

//
// File trailer for str2double.cpp
//
// [EOF]
//
