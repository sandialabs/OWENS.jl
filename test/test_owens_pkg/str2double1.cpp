//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: str2double1.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "str2double1.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include <string.h>

// Function Declarations
static boolean_T c_copydigits(emxArray_char_T *s1, int *idx, const
  emxArray_char_T *s, int *k, int n, boolean_T allowpoint);
static boolean_T d_copydigits(emxArray_char_T *s1, int *idx, const
  emxArray_char_T *s, int *k, int n, boolean_T allowpoint);

// Function Definitions

//
// Arguments    : emxArray_char_T *s1
//                int *idx
//                const emxArray_char_T *s
//                int *k
//                int n
//                boolean_T allowpoint
// Return Type  : boolean_T
//
static boolean_T c_copydigits(emxArray_char_T *s1, int *idx, const
  emxArray_char_T *s, int *k, int n, boolean_T allowpoint)
{
  boolean_T success;
  boolean_T haspoint;
  boolean_T exitg1;
  char c;
  success = (*k <= n);
  haspoint = false;
  exitg1 = false;
  while ((!exitg1) && (success && (*k <= n))) {
    c = s->data[*k - 1];
    if ((c >= '0') && (c <= '9')) {
      s1->data[*idx - 1] = c;
      (*idx)++;
      (*k)++;
    } else if (c == '.') {
      if (allowpoint && (!haspoint)) {
        success = true;
      } else {
        success = false;
      }

      if (success) {
        s1->data[*idx - 1] = '.';
        (*idx)++;
        haspoint = true;
      }

      (*k)++;
    } else if (c == ',') {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return success;
}

//
// Arguments    : emxArray_char_T *s1
//                int *idx
//                const emxArray_char_T *s
//                int *k
//                int n
//                boolean_T allowpoint
// Return Type  : boolean_T
//
static boolean_T d_copydigits(emxArray_char_T *s1, int *idx, const
  emxArray_char_T *s, int *k, int n, boolean_T allowpoint)
{
  boolean_T success;
  boolean_T haspoint;
  boolean_T exitg1;
  success = (*k <= n);
  haspoint = false;
  exitg1 = false;
  while ((!exitg1) && (success && (*k <= n))) {
    if ((s->data[*k - 1] >= '0') && (s->data[*k - 1] <= '9')) {
      s1->data[*idx - 1] = s->data[*k - 1];
      (*idx)++;
      (*k)++;
    } else if (s->data[*k - 1] == '.') {
      if (allowpoint && (!haspoint)) {
        success = true;
      } else {
        success = false;
      }

      if (success) {
        s1->data[*idx - 1] = '.';
        (*idx)++;
        haspoint = true;
      }

      (*k)++;
    } else if (s->data[*k - 1] == ',') {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return success;
}

//
// Arguments    : char s1[3]
//                int *idx
//                char s
//                int *k
//                boolean_T allowpoint
// Return Type  : boolean_T
//
boolean_T b_copydigits(char s1[3], int *idx, char s, int *k, boolean_T
  allowpoint)
{
  boolean_T success;
  boolean_T haspoint;
  boolean_T exitg1;
  success = (*k <= 1);
  haspoint = false;
  exitg1 = false;
  while ((!exitg1) && (success && (*k <= 1))) {
    if ((s >= '0') && (s <= '9')) {
      s1[*idx - 1] = s;
      (*idx)++;
      *k = 2;
    } else if (s == '.') {
      if (allowpoint && (!haspoint)) {
        success = true;
      } else {
        success = false;
      }

      if (success) {
        s1[*idx - 1] = '.';
        (*idx)++;
        haspoint = true;
      }

      *k = 2;
    } else if (s == ',') {
      *k = 2;
    } else {
      exitg1 = true;
    }
  }

  return success;
}

//
// Arguments    : char s1[3]
//                int *idx
//                char s
//                int *k
// Return Type  : boolean_T
//
boolean_T b_copyexponent(char s1[3], int *idx, char s, int *k)
{
  boolean_T success;
  boolean_T b_success;
  success = true;
  if ((*k <= 1) && ((s == 'E') || (s == 'e'))) {
    s1[*idx - 1] = 'e';
    (*idx)++;
    *k = 2;
    b_success = b_copydigits(s1, idx, s, k, false);
    if ((!b_success) || (*k <= 2)) {
      success = false;
    }
  }

  return success;
}

//
// Arguments    : char s
//                int *k
//                boolean_T *b_finite
//                double *fv
// Return Type  : void
//
void b_readNonFinite(char s, int *k, boolean_T *b_finite, double *fv)
{
  int ksaved;
  char c_idx_0;
  char c_idx_1;
  char c_idx_2;
  ksaved = *k;
  c_idx_0 = '\x00';
  while ((*k <= 1) && (s == ',')) {
    *k = 2;
  }

  if (*k <= 1) {
    c_idx_0 = s;
  }

  (*k)++;
  c_idx_1 = '\x00';
  while ((*k <= 1) && (s == ',')) {
    *k = 2;
  }

  if (*k <= 1) {
    c_idx_1 = s;
  }

  (*k)++;
  c_idx_2 = '\x00';
  while ((*k <= 1) && (s == ',')) {
    *k = 2;
  }

  if (*k <= 1) {
    c_idx_2 = s;
  }

  (*k)++;
  if (((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') || (c_idx_1 ==
        'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) {
    *b_finite = false;
    *fv = rtInf;
  } else if (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
              (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n'))) {
    *b_finite = false;
    *fv = rtNaN;
  } else {
    *b_finite = true;
    *fv = 0.0;
    *k = ksaved;
  }
}

//
// Arguments    : char s1[3]
//                int *idx
//                char s
//                int *k
//                boolean_T allowimag
//                boolean_T *isimag
//                boolean_T *b_finite
//                double *nfv
//                boolean_T *foundsign
//                boolean_T *success
// Return Type  : void
//
void b_readfloat(char s1[3], int *idx, char s, int *k, boolean_T allowimag,
                 boolean_T *isimag, boolean_T *b_finite, double *nfv, boolean_T *
                 foundsign, boolean_T *success)
{
  int b_idx;
  boolean_T isneg;
  boolean_T exitg1;
  boolean_T unusedU2;
  *isimag = false;
  *b_finite = true;
  *nfv = 0.0;
  b_idx = *idx - 1;
  isneg = false;
  *foundsign = false;
  exitg1 = false;
  while ((!exitg1) && (*k <= 1)) {
    if (s == '-') {
      isneg = !isneg;
      *foundsign = true;
      *k = 2;
    } else if (s == ',') {
      *k = 2;
    } else if (s == '+') {
      *foundsign = true;
      *k = 2;
    } else if (!bv[static_cast<unsigned char>(s) & 127]) {
      exitg1 = true;
    } else {
      *k = 2;
    }
  }

  *success = (*k <= 1);
  if ((*success) && isneg) {
    if (s1[*idx - 2] == '-') {
      s1[*idx - 2] = ' ';
    } else {
      s1[*idx - 1] = '-';
      b_idx = *idx;
    }
  }

  *idx = b_idx + 1;
  if (*success) {
    isneg = false;
    if (*k <= 1) {
      isneg = ((s == 'j') || (s == 'i'));
    }

    if (isneg) {
      if (allowimag) {
        *isimag = true;
        (*k)++;
        while ((*k <= 1) && (bv[static_cast<unsigned char>(s) & 127] || (s ==
                 '\x00') || (s == ','))) {
          *k = 2;
        }

        if ((*k <= 1) && (s == '*')) {
          *k = 2;
          b_readfloat(s1, idx, '*', k, false, &isneg, b_finite, nfv, &unusedU2,
                      success);
        } else {
          s1[b_idx] = '1';
          *idx = b_idx + 2;
        }
      } else {
        *success = false;
      }
    } else {
      b_readNonFinite(s, k, b_finite, nfv);
      if (*b_finite) {
        *success = b_copydigits(s1, idx, s, k, true);
        if (*success) {
          *success = b_copyexponent(s1, idx, s, k);
        }
      } else {
        if (s1[b_idx - 1] == '-') {
          *idx = b_idx;
          s1[b_idx - 1] = ' ';
          *nfv = -*nfv;
        }
      }

      while ((*k <= 1) && (bv[static_cast<unsigned char>(s) & 127] || (s ==
               '\x00') || (s == ','))) {
        *k = 2;
      }

      if ((*k <= 1) && (s == '*')) {
        *k = 2;
      }

      if ((*k <= 1) && ((s == 'i') || (s == 'j'))) {
        *k = 2;
        *isimag = true;
      }
    }

    while ((*k <= 1) && (bv[static_cast<unsigned char>(s) & 127] || (s == '\x00')
                         || (s == ','))) {
      *k = 2;
    }
  }
}

//
// Arguments    : emxArray_char_T *s1
//                int *idx
//                const emxArray_char_T *s
//                int *k
//                int n
//                boolean_T allowimag
//                boolean_T *isimag
//                boolean_T *b_finite
//                double *nfv
//                boolean_T *foundsign
//                boolean_T *success
// Return Type  : void
//
void c_readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s, int *k,
                 int n, boolean_T allowimag, boolean_T *isimag, boolean_T
                 *b_finite, double *nfv, boolean_T *foundsign, boolean_T
                 *success)
{
  int b_idx;
  boolean_T isneg;
  boolean_T exitg1;
  char c_idx_0;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int b_k;
  int exitg3;
  boolean_T unusedU2;
  char c_idx_1;
  char c_idx_2;
  *isimag = false;
  *b_finite = true;
  *nfv = 0.0;
  b_idx = *idx;
  isneg = false;
  *foundsign = false;
  exitg1 = false;
  while ((!exitg1) && (*k <= n)) {
    c_idx_0 = s->data[*k - 1];
    if (c_idx_0 == '-') {
      isneg = !isneg;
      *foundsign = true;
      (*k)++;
    } else if (c_idx_0 == ',') {
      (*k)++;
    } else if (c_idx_0 == '+') {
      *foundsign = true;
      (*k)++;
    } else if (!bv[static_cast<unsigned char>(c_idx_0) & 127]) {
      exitg1 = true;
    } else {
      (*k)++;
    }
  }

  *success = (*k <= n);
  if ((*success) && isneg) {
    if ((*idx >= 2) && (s1->data[*idx - 2] == '-')) {
      s1->data[*idx - 2] = ' ';
    } else {
      s1->data[*idx - 1] = '-';
      b_idx = *idx + 1;
    }
  }

  *idx = b_idx;
  if (*success) {
    guard1 = false;
    guard2 = false;
    if (*k <= n) {
      c_idx_0 = s->data[*k - 1];
      if (c_idx_0 == 'j') {
        guard2 = true;
      } else if (c_idx_0 == 'i') {
        if (*k >= n - 1) {
          guard2 = true;
        } else {
          b_k = *k;
          c_idx_0 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_0 = s->data[b_k - 1];
          }

          b_k++;
          c_idx_1 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_1 = s->data[b_k - 1];
          }

          b_k++;
          c_idx_2 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_2 = s->data[b_k - 1];
          }

          if ((((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') ||
                (c_idx_1 == 'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) ||
              (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
                (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n')))) {
            guard1 = true;
          } else {
            guard2 = true;
          }
        }
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard2) {
      if (allowimag) {
        *isimag = true;
        (*k)++;
        exitg1 = false;
        while ((!exitg1) && (*k <= n)) {
          c_idx_0 = s->data[*k - 1];
          if (bv[static_cast<unsigned char>(c_idx_0) & 127] || (c_idx_0 ==
               '\x00') || (c_idx_0 == ',')) {
            (*k)++;
          } else {
            exitg1 = true;
          }
        }

        if ((*k <= n) && (s->data[*k - 1] == '*')) {
          (*k)++;
          c_readfloat(s1, idx, s, k, n, false, &isneg, b_finite, nfv, &unusedU2,
                      success);
        } else {
          s1->data[b_idx - 1] = '1';
          *idx = b_idx + 1;
        }
      } else {
        *success = false;
      }
    }

    if (guard1) {
      b_k = *k;
      c_idx_0 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_0 = s->data[b_k - 1];
      }

      b_k++;
      c_idx_1 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_1 = s->data[b_k - 1];
      }

      b_k++;
      c_idx_2 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_2 = s->data[b_k - 1];
      }

      b_k++;
      if (((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') ||
           (c_idx_1 == 'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) {
        *b_finite = false;
        *nfv = rtInf;
      } else if (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
                  (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n'))) {
        *b_finite = false;
        *nfv = rtNaN;
      } else {
        *b_finite = true;
        *nfv = 0.0;
        b_k = *k;
      }

      *k = b_k;
      if (*b_finite) {
        *success = c_copydigits(s1, &b_idx, s, &b_k, n, true);
        *idx = b_idx;
        *k = b_k;
        if (*success) {
          *success = true;
          if (b_k <= n) {
            c_idx_0 = s->data[b_k - 1];
            if ((c_idx_0 == 'E') || (c_idx_0 == 'e')) {
              s1->data[b_idx - 1] = 'e';
              *idx = b_idx + 1;
              b_k++;
              while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
                b_k++;
              }

              if (b_k <= n) {
                c_idx_0 = s->data[b_k - 1];
                if (c_idx_0 == '-') {
                  s1->data[b_idx] = '-';
                  *idx = b_idx + 2;
                  b_k++;
                } else {
                  if (c_idx_0 == '+') {
                    b_k++;
                  }
                }
              }

              *k = b_k;
              isneg = c_copydigits(s1, idx, s, k, n, false);
              if ((!isneg) || (*k <= b_k)) {
                *success = false;
              }
            }
          }
        }
      } else {
        if ((b_idx >= 2) && (s1->data[b_idx - 2] == '-')) {
          *idx = b_idx - 1;
          s1->data[b_idx - 2] = ' ';
          *nfv = -*nfv;
        }
      }

      exitg1 = false;
      while ((!exitg1) && (*k <= n)) {
        if (bv[static_cast<unsigned char>(s->data[*k - 1]) & 127]) {
          (*k)++;
        } else {
          c_idx_0 = s->data[*k - 1];
          if ((c_idx_0 == '\x00') || (c_idx_0 == ',')) {
            (*k)++;
          } else {
            exitg1 = true;
          }
        }
      }

      if ((*k <= n) && (s->data[*k - 1] == '*')) {
        (*k)++;
        while ((*k <= n) && (bv[static_cast<unsigned char>(s->data[*k - 1]) &
                             127] || (s->data[*k - 1] == '\x00') || (s->data[*k
                 - 1] == ','))) {
          (*k)++;
        }
      }

      if (*k <= n) {
        c_idx_0 = s->data[*k - 1];
        if ((c_idx_0 == 'i') || (c_idx_0 == 'j')) {
          (*k)++;
          *isimag = true;
        }
      }
    }

    do {
      exitg3 = 0;
      if (*k <= n) {
        if (bv[static_cast<unsigned char>(s->data[*k - 1]) & 127]) {
          (*k)++;
        } else {
          c_idx_0 = s->data[*k - 1];
          if ((c_idx_0 == '\x00') || (c_idx_0 == ',')) {
            (*k)++;
          } else {
            exitg3 = 1;
          }
        }
      } else {
        exitg3 = 1;
      }
    } while (exitg3 == 0);
  }
}

//
// Arguments    : char s1[4]
//                int *idx
//                const char s[2]
//                int *k
//                boolean_T allowpoint
// Return Type  : boolean_T
//
boolean_T copydigits(char s1[4], int *idx, const char s[2], int *k, boolean_T
                     allowpoint)
{
  boolean_T success;
  boolean_T haspoint;
  boolean_T exitg1;
  char c;
  success = (*k <= 2);
  haspoint = false;
  exitg1 = false;
  while ((!exitg1) && (success && (*k <= 2))) {
    c = s[*k - 1];
    if ((c >= '0') && (c <= '9')) {
      s1[*idx - 1] = c;
      (*idx)++;
      (*k)++;
    } else if (c == '.') {
      if (allowpoint && (!haspoint)) {
        success = true;
      } else {
        success = false;
      }

      if (success) {
        s1[*idx - 1] = '.';
        (*idx)++;
        haspoint = true;
      }

      (*k)++;
    } else if (c == ',') {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return success;
}

//
// Arguments    : char s1[4]
//                int *idx
//                const char s[2]
//                int *k
// Return Type  : boolean_T
//
boolean_T copyexponent(char s1[4], int *idx, const char s[2], int *k)
{
  boolean_T success;
  char c;
  int kexp;
  boolean_T b_success;
  success = true;
  if (*k <= 2) {
    c = s[*k - 1];
    if ((c == 'E') || (c == 'e')) {
      s1[*idx - 1] = 'e';
      (*idx)++;
      (*k)++;
      while ((*k <= 2) && (s[1] == ',')) {
        *k = 3;
      }

      if (*k <= 2) {
        if (s[1] == '-') {
          s1[*idx - 1] = '-';
          (*idx)++;
          *k = 3;
        } else {
          if (s[1] == '+') {
            *k = 3;
          }
        }
      }

      kexp = *k;
      b_success = copydigits(s1, idx, s, k, false);
      if ((!b_success) || (*k <= kexp)) {
        success = false;
      }
    }
  }

  return success;
}

//
// Arguments    : emxArray_char_T *s1
//                int *idx
//                const emxArray_char_T *s
//                int *k
//                int n
//                boolean_T allowimag
//                boolean_T *isimag
//                boolean_T *b_finite
//                double *nfv
//                boolean_T *foundsign
//                boolean_T *success
// Return Type  : void
//
void d_readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s, int *k,
                 int n, boolean_T allowimag, boolean_T *isimag, boolean_T
                 *b_finite, double *nfv, boolean_T *foundsign, boolean_T
                 *success)
{
  int b_idx;
  boolean_T isneg;
  boolean_T exitg1;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int b_k;
  char c_idx_0;
  boolean_T unusedU2;
  char c_idx_1;
  char c_idx_2;
  *isimag = false;
  *b_finite = true;
  *nfv = 0.0;
  b_idx = *idx;
  isneg = false;
  *foundsign = false;
  exitg1 = false;
  while ((!exitg1) && (*k <= n)) {
    if (s->data[*k - 1] == '-') {
      isneg = !isneg;
      *foundsign = true;
      (*k)++;
    } else if (s->data[*k - 1] == ',') {
      (*k)++;
    } else if (s->data[*k - 1] == '+') {
      *foundsign = true;
      (*k)++;
    } else if (!bv[static_cast<unsigned char>(s->data[*k - 1]) & 127]) {
      exitg1 = true;
    } else {
      (*k)++;
    }
  }

  *success = (*k <= n);
  if ((*success) && isneg) {
    if ((*idx >= 2) && (s1->data[*idx - 2] == '-')) {
      s1->data[*idx - 2] = ' ';
    } else {
      s1->data[*idx - 1] = '-';
      b_idx = *idx + 1;
    }
  }

  *idx = b_idx;
  if (*success) {
    guard1 = false;
    guard2 = false;
    if (*k <= n) {
      if (s->data[*k - 1] == 'j') {
        guard2 = true;
      } else if (s->data[*k - 1] == 'i') {
        if (*k >= n - 1) {
          guard2 = true;
        } else {
          b_k = *k;
          c_idx_0 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_0 = s->data[b_k - 1];
          }

          b_k++;
          c_idx_1 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_1 = s->data[b_k - 1];
          }

          b_k++;
          c_idx_2 = '\x00';
          while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
            b_k++;
          }

          if (b_k <= n) {
            c_idx_2 = s->data[b_k - 1];
          }

          if ((((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') ||
                (c_idx_1 == 'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) ||
              (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
                (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n')))) {
            guard1 = true;
          } else {
            guard2 = true;
          }
        }
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard2) {
      if (allowimag) {
        *isimag = true;
        (*k)++;
        while ((*k <= n) && (bv[static_cast<unsigned char>(s->data[*k - 1]) &
                             127] || (s->data[*k - 1] == '\x00') || (s->data[*k
                 - 1] == ','))) {
          (*k)++;
        }

        if ((*k <= n) && (s->data[*k - 1] == '*')) {
          (*k)++;
          d_readfloat(s1, idx, s, k, n, false, &isneg, b_finite, nfv, &unusedU2,
                      success);
        } else {
          s1->data[b_idx - 1] = '1';
          *idx = b_idx + 1;
        }
      } else {
        *success = false;
      }
    }

    if (guard1) {
      b_k = *k;
      c_idx_0 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_0 = s->data[b_k - 1];
      }

      b_k++;
      c_idx_1 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_1 = s->data[b_k - 1];
      }

      b_k++;
      c_idx_2 = '\x00';
      while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
        b_k++;
      }

      if (b_k <= n) {
        c_idx_2 = s->data[b_k - 1];
      }

      b_k++;
      if (((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') ||
           (c_idx_1 == 'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) {
        *b_finite = false;
        *nfv = rtInf;
      } else if (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
                  (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n'))) {
        *b_finite = false;
        *nfv = rtNaN;
      } else {
        *b_finite = true;
        *nfv = 0.0;
        b_k = *k;
      }

      *k = b_k;
      if (*b_finite) {
        *success = d_copydigits(s1, &b_idx, s, &b_k, n, true);
        *idx = b_idx;
        *k = b_k;
        if (*success) {
          *success = true;
          if ((b_k <= n) && ((s->data[b_k - 1] == 'E') || (s->data[b_k - 1] ==
                'e'))) {
            s1->data[b_idx - 1] = 'e';
            *idx = b_idx + 1;
            b_k++;
            while ((b_k <= n) && (s->data[b_k - 1] == ',')) {
              b_k++;
            }

            if (b_k <= n) {
              if (s->data[b_k - 1] == '-') {
                s1->data[b_idx] = '-';
                *idx = b_idx + 2;
                b_k++;
              } else {
                if (s->data[b_k - 1] == '+') {
                  b_k++;
                }
              }
            }

            *k = b_k;
            isneg = d_copydigits(s1, idx, s, k, n, false);
            if ((!isneg) || (*k <= b_k)) {
              *success = false;
            }
          }
        }
      } else {
        if ((b_idx >= 2) && (s1->data[b_idx - 2] == '-')) {
          *idx = b_idx - 1;
          s1->data[b_idx - 2] = ' ';
          *nfv = -*nfv;
        }
      }

      while ((*k <= n) && (bv[static_cast<unsigned char>(s->data[*k - 1]) & 127]
                           || (s->data[*k - 1] == '\x00') || (s->data[*k - 1] ==
               ','))) {
        (*k)++;
      }

      if ((*k <= n) && (s->data[*k - 1] == '*')) {
        (*k)++;
        while ((*k <= n) && (bv[static_cast<unsigned char>(s->data[*k - 1]) &
                             127] || (s->data[*k - 1] == '\x00') || (s->data[*k
                 - 1] == ','))) {
          (*k)++;
        }
      }

      if ((*k <= n) && ((s->data[*k - 1] == 'i') || (s->data[*k - 1] == 'j'))) {
        (*k)++;
        *isimag = true;
      }
    }

    while ((*k <= n) && (bv[static_cast<unsigned char>(s->data[*k - 1]) & 127] ||
                         (s->data[*k - 1] == '\x00') || (s->data[*k - 1] == ',')))
    {
      (*k)++;
    }
  }
}

//
// Arguments    : const char s[2]
//                int *k
//                boolean_T *b_finite
//                double *fv
// Return Type  : void
//
void readNonFinite(const char s[2], int *k, boolean_T *b_finite, double *fv)
{
  int ksaved;
  char c_idx_0;
  char c_idx_1;
  char c_idx_2;
  ksaved = *k;
  c_idx_0 = '\x00';
  while ((*k <= 2) && (s[*k - 1] == ',')) {
    (*k)++;
  }

  if (*k <= 2) {
    c_idx_0 = s[*k - 1];
  }

  (*k)++;
  c_idx_1 = '\x00';
  while ((*k <= 2) && (s[*k - 1] == ',')) {
    (*k)++;
  }

  if (*k <= 2) {
    c_idx_1 = s[*k - 1];
  }

  (*k)++;
  c_idx_2 = '\x00';
  while ((*k <= 2) && (s[*k - 1] == ',')) {
    (*k)++;
  }

  if (*k <= 2) {
    c_idx_2 = s[*k - 1];
  }

  (*k)++;
  if (((c_idx_0 == 'I') || (c_idx_0 == 'i')) && ((c_idx_1 == 'N') || (c_idx_1 ==
        'n')) && ((c_idx_2 == 'F') || (c_idx_2 == 'f'))) {
    *b_finite = false;
    *fv = rtInf;
  } else if (((c_idx_0 == 'N') || (c_idx_0 == 'n')) && ((c_idx_1 == 'A') ||
              (c_idx_1 == 'a')) && ((c_idx_2 == 'N') || (c_idx_2 == 'n'))) {
    *b_finite = false;
    *fv = rtNaN;
  } else {
    *b_finite = true;
    *fv = 0.0;
    *k = ksaved;
  }
}

//
// Arguments    : char s1[4]
//                int *idx
//                const char s[2]
//                int *k
//                boolean_T allowimag
//                boolean_T *isimag
//                boolean_T *b_finite
//                double *nfv
//                boolean_T *foundsign
//                boolean_T *success
// Return Type  : void
//
void readfloat(char s1[4], int *idx, const char s[2], int *k, boolean_T
               allowimag, boolean_T *isimag, boolean_T *b_finite, double *nfv,
               boolean_T *foundsign, boolean_T *success)
{
  int b_idx;
  boolean_T isneg;
  boolean_T exitg1;
  char p_tmp;
  boolean_T unusedU2;
  *isimag = false;
  *b_finite = true;
  *nfv = 0.0;
  b_idx = *idx;
  isneg = false;
  *foundsign = false;
  exitg1 = false;
  while ((!exitg1) && (*k <= 2)) {
    p_tmp = s[*k - 1];
    if (p_tmp == '-') {
      isneg = !isneg;
      *foundsign = true;
      (*k)++;
    } else if (p_tmp == ',') {
      (*k)++;
    } else if (p_tmp == '+') {
      *foundsign = true;
      (*k)++;
    } else if (!bv[static_cast<unsigned char>(p_tmp) & 127]) {
      exitg1 = true;
    } else {
      (*k)++;
    }
  }

  *success = (*k <= 2);
  if ((*success) && isneg) {
    if ((*idx >= 2) && (s1[*idx - 2] == '-')) {
      s1[*idx - 2] = ' ';
    } else {
      s1[*idx - 1] = '-';
      b_idx = *idx + 1;
    }
  }

  *idx = b_idx;
  if (*success) {
    isneg = false;
    if (*k <= 2) {
      p_tmp = s[*k - 1];
      isneg = ((p_tmp == 'j') || (p_tmp == 'i'));
    }

    if (isneg) {
      if (allowimag) {
        *isimag = true;
        (*k)++;
        exitg1 = false;
        while ((!exitg1) && (*k <= 2)) {
          p_tmp = s[*k - 1];
          if (bv[static_cast<unsigned char>(p_tmp) & 127] || (p_tmp == '\x00') ||
              (p_tmp == ',')) {
            (*k)++;
          } else {
            exitg1 = true;
          }
        }

        if ((*k <= 2) && (s[*k - 1] == '*')) {
          (*k)++;
          readfloat(s1, idx, s, k, false, &isneg, b_finite, nfv, &unusedU2,
                    success);
        } else {
          s1[b_idx - 1] = '1';
          *idx = b_idx + 1;
        }
      } else {
        *success = false;
      }
    } else {
      readNonFinite(s, k, b_finite, nfv);
      if (*b_finite) {
        *success = copydigits(s1, idx, s, k, true);
        if (*success) {
          *success = copyexponent(s1, idx, s, k);
        }
      } else {
        if ((b_idx >= 2) && (s1[b_idx - 2] == '-')) {
          *idx = b_idx - 1;
          s1[b_idx - 2] = ' ';
          *nfv = -*nfv;
        }
      }

      exitg1 = false;
      while ((!exitg1) && (*k <= 2)) {
        if (bv[static_cast<unsigned char>(s[*k - 1]) & 127]) {
          (*k)++;
        } else {
          p_tmp = s[*k - 1];
          if ((p_tmp == '\x00') || (p_tmp == ',')) {
            (*k)++;
          } else {
            exitg1 = true;
          }
        }
      }

      if ((*k <= 2) && (s[*k - 1] == '*')) {
        (*k)++;
        while ((*k <= 2) && (bv[static_cast<unsigned char>(s[*k - 1]) & 127] ||
                             (s[*k - 1] == '\x00') || (s[*k - 1] == ','))) {
          (*k)++;
        }
      }

      if (*k <= 2) {
        p_tmp = s[*k - 1];
        if ((p_tmp == 'i') || (p_tmp == 'j')) {
          (*k)++;
          *isimag = true;
        }
      }
    }

    exitg1 = false;
    while ((!exitg1) && (*k <= 2)) {
      p_tmp = s[*k - 1];
      if (bv[static_cast<unsigned char>(p_tmp) & 127] || (p_tmp == '\x00') ||
          (p_tmp == ',')) {
        (*k)++;
      } else {
        exitg1 = true;
      }
    }
  }
}

//
// File trailer for str2double1.cpp
//
// [EOF]
//
