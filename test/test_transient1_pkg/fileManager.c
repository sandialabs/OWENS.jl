/*
 * File: fileManager.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include "test_transient1_rtwutil.h"
#include <string.h>

/* Variable Definitions */
static FILE * eml_openfiles[20];
static boolean_T eml_autoflush[20];

/* Function Declarations */
static signed char filedata(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : signed char
 */
static signed char filedata(void)
{
  signed char f;
  int k;
  boolean_T exitg1;
  f = 0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (eml_openfiles[k] == NULL) {
      f = (signed char)(k + 1);
      exitg1 = true;
    } else {
      k++;
    }
  }

  return f;
}

/*
 * Arguments    : const emxArray_char_T *cfilename
 *                const char * cpermission
 * Return Type  : signed char
 */
signed char b_cfopen(const emxArray_char_T *cfilename, const char * cpermission)
{
  signed char fileid;
  signed char j;
  emxArray_char_T *ccfilename;
  int input_sizes_idx_1;
  int i;
  FILE * filestar;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    emxInit_char_T(&ccfilename, 2);
    if ((cfilename->size[0] != 0) && (cfilename->size[1] != 0)) {
      input_sizes_idx_1 = cfilename->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }

    i = ccfilename->size[0] * ccfilename->size[1];
    ccfilename->size[0] = 1;
    ccfilename->size[1] = input_sizes_idx_1 + 1;
    emxEnsureCapacity_char_T(ccfilename, i);
    for (i = 0; i < input_sizes_idx_1; i++) {
      ccfilename->data[i] = cfilename->data[i];
    }

    ccfilename->data[input_sizes_idx_1] = '\x00';
    filestar = fopen(&ccfilename->data[0], cpermission);
    emxFree_char_T(&ccfilename);
    if (filestar != NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (i > 127) {
        i = 127;
      }

      fileid = (signed char)i;
    }
  }

  return fileid;
}

/*
 * Arguments    : const emxArray_char_T *cfilename
 *                const char * cpermission
 * Return Type  : signed char
 */
signed char c_cfopen(const emxArray_char_T *cfilename, const char * cpermission)
{
  signed char fileid;
  signed char j;
  emxArray_char_T *ccfilename;
  int i;
  int loop_ub;
  FILE * filestar;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    emxInit_char_T(&ccfilename, 2);
    i = ccfilename->size[0] * ccfilename->size[1];
    ccfilename->size[0] = 1;
    ccfilename->size[1] = cfilename->size[1] + 1;
    emxEnsureCapacity_char_T(ccfilename, i);
    loop_ub = cfilename->size[1];
    for (i = 0; i < loop_ub; i++) {
      ccfilename->data[i] = cfilename->data[i];
    }

    ccfilename->data[cfilename->size[1]] = '\x00';
    filestar = fopen(&ccfilename->data[0], cpermission);
    emxFree_char_T(&ccfilename);
    if (filestar != NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (i > 127) {
        i = 127;
      }

      fileid = (signed char)i;
    }
  }

  return fileid;
}

/*
 * Arguments    : double fid
 * Return Type  : int
 */
int cfclose(double fid)
{
  int st;
  signed char fileid;
  signed char b_fileid;
  FILE * filestar;
  int cst;
  st = -1;
  fileid = (signed char)rt_roundd_snf(fid);
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }

  b_fileid = fileid;
  if (fileid < 0) {
    b_fileid = -1;
  }

  if (b_fileid >= 3) {
    filestar = eml_openfiles[b_fileid - 3];
  } else if (b_fileid == 0) {
    filestar = stdin;
  } else if (b_fileid == 1) {
    filestar = stdout;
  } else if (b_fileid == 2) {
    filestar = stderr;
  } else {
    filestar = NULL;
  }

  if ((filestar != NULL) && (fileid >= 3)) {
    cst = fclose(filestar);
    if (cst == 0) {
      st = 0;
      cst = fileid - 3;
      eml_openfiles[cst] = NULL;
      eml_autoflush[cst] = true;
    }
  }

  return st;
}

/*
 * Arguments    : const char * cfilename
 *                const char * cpermission
 * Return Type  : signed char
 */
signed char cfopen(const char * cfilename, const char * cpermission)
{
  signed char fileid;
  signed char j;
  FILE * filestar;
  int i;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    filestar = fopen(cfilename, cpermission);
    if (filestar != NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (i > 127) {
        i = 127;
      }

      fileid = (signed char)i;
    }
  }

  return fileid;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void filedata_init(void)
{
  FILE * a;
  int i;
  a = NULL;
  for (i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

/*
 * Arguments    : double fid
 *                FILE * *filestar
 *                boolean_T *autoflush
 * Return Type  : void
 */
void getfilestar(double fid, FILE * *filestar, boolean_T *autoflush)
{
  signed char fileid;
  fileid = (signed char)rt_roundd_snf(fid);
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    *filestar = eml_openfiles[fileid - 3];
    *autoflush = eml_autoflush[fileid - 3];
  } else if (fileid == 0) {
    *filestar = stdin;
    *autoflush = true;
  } else if (fileid == 1) {
    *filestar = stdout;
    *autoflush = true;
  } else if (fileid == 2) {
    *filestar = stderr;
    *autoflush = true;
  } else {
    *filestar = NULL;
    *autoflush = true;
  }
}

/*
 * File trailer for fileManager.c
 *
 * [EOF]
 */
