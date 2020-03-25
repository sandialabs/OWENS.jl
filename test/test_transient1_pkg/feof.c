/*
 * File: feof.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "feof.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <string.h>

/* Function Definitions */

/*
 * Arguments    : double fileID
 * Return Type  : double
 */
double b_feof(double fileID)
{
  double st;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T a;
  int b_st;
  b_NULL = NULL;
  getfilestar(fileID, &filestar, &a);
  if (filestar == b_NULL) {
    st = 0.0;
  } else {
    b_st = feof(filestar);
    st = ((int)b_st != 0);
  }

  return st;
}

/*
 * File trailer for feof.c
 *
 * [EOF]
 */
