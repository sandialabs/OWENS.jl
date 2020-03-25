/*
 * File: importCactusFile.h
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

#ifndef IMPORTCACTUSFILE_H
#define IMPORTCACTUSFILE_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "test_transient1_types.h"

/* Function Declarations */
extern void b_importCactusFile(const emxArray_char_T *filename, emxArray_real_T *
  data);
extern void importCactusFile(const emxArray_char_T *filename, double data_data[],
  int data_size[2]);

#endif

/*
 * File trailer for importCactusFile.h
 *
 * [EOF]
 */
