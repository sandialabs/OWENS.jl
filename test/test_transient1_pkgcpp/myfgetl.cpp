//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: myfgetl.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "myfgetl.h"
#include "fileManager.h"
#include "fread.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : double fid
// Return Type  : void
//
void b_myfgetl(double fid)
{
  char c_data[1];
  int c_size[1];
  int exitg1;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T a;
  int st;
  int x;

  //  C++ compatable workaround for fget1 from https://www.mathworks.com/matlabcentral/answers/461159-read-text-file-line-by-line-in-deployed-application 
  b_fread(fid, c_data, c_size);
  do {
    exitg1 = 0;
    b_NULL = NULL;
    getfilestar(fid, &filestar, &a);
    if (filestar == b_NULL) {
      x = 0;
    } else {
      st = feof(filestar);
      x = ((int)st != 0);
    }

    if (x == 0) {
      a = false;
      if (c_size[0] == 1) {
        a = !(c_data[0] != '\x0a');
      }

      if (!a) {
        b_fread(fid, c_data, c_size);
      } else {
        exitg1 = 1;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

//
// Arguments    : double fid
//                emxArray_char_T *line
// Return Type  : void
//
void myfgetl(double fid, emxArray_char_T *line)
{
  char c_data[1];
  int c_size[1];
  emxArray_char_T *b_line;
  int exitg1;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T a;
  int st;
  int input_sizes_idx_1;
  signed char b_input_sizes_idx_1;
  int unnamed_idx_1;
  int i;
  int i1;

  //  C++ compatable workaround for fget1 from https://www.mathworks.com/matlabcentral/answers/461159-read-text-file-line-by-line-in-deployed-application 
  line->size[0] = 1;
  line->size[1] = 0;
  b_fread(fid, c_data, c_size);
  emxInit_char_T(&b_line, 2);
  do {
    exitg1 = 0;
    b_NULL = NULL;
    getfilestar(fid, &filestar, &a);
    if (filestar == b_NULL) {
      input_sizes_idx_1 = 0;
    } else {
      st = feof(filestar);
      input_sizes_idx_1 = ((int)st != 0);
    }

    if (input_sizes_idx_1 == 0) {
      a = false;
      if (c_size[0] == 1) {
        a = !(c_data[0] != '\x0a');
      }

      if (!a) {
        a = false;
        if (c_size[0] == 1) {
          a = !(c_data[0] != '\x0d');
        }

        if (!a) {
          if (line->size[1] != 0) {
            input_sizes_idx_1 = line->size[1];
          } else {
            input_sizes_idx_1 = 0;
          }

          b_input_sizes_idx_1 = static_cast<signed char>((c_size[0] != 0));
          if (line->size[1] != 0) {
            unnamed_idx_1 = line->size[1];
          } else {
            unnamed_idx_1 = 0;
          }

          i = b_line->size[0] * b_line->size[1];
          b_line->size[0] = 1;
          b_line->size[1] = input_sizes_idx_1 + b_input_sizes_idx_1;
          emxEnsureCapacity_char_T(b_line, i);
          for (i = 0; i < input_sizes_idx_1; i++) {
            for (i1 = 0; i1 < 1; i1++) {
              b_line->data[b_line->size[0] * i] = line->data[line->size[0] * i];
            }
          }

          input_sizes_idx_1 = b_input_sizes_idx_1;
          for (i = 0; i < input_sizes_idx_1; i++) {
            for (i1 = 0; i1 < 1; i1++) {
              b_line->data[b_line->size[0] * unnamed_idx_1] = c_data[0];
            }
          }

          i = line->size[0] * line->size[1];
          line->size[0] = b_line->size[0];
          line->size[1] = b_line->size[1];
          emxEnsureCapacity_char_T(line, i);
          input_sizes_idx_1 = b_line->size[0] * b_line->size[1];
          for (i = 0; i < input_sizes_idx_1; i++) {
            line->data[i] = b_line->data[i];
          }
        }

        b_fread(fid, c_data, c_size);
      } else {
        exitg1 = 1;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_line);
}

//
// File trailer for myfgetl.cpp
//
// [EOF]
//
