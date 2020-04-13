//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: fwrite.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "fwrite.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : double fileID
//                const emxArray_char_T *x
// Return Type  : void
//
void b_fwrite(double fileID, const emxArray_char_T *x)
{
  FILE * filestar;
  boolean_T autoflush;
  size_t bytesOutSizet;
  getfilestar(fileID, &filestar, &autoflush);
  if (!(fileID != 0.0)) {
    filestar = NULL;
  }

  if (!(filestar == NULL)) {
    bytesOutSizet = fwrite(&x->data[0], sizeof(char), (size_t)x->size[1],
      filestar);
    if (((double)bytesOutSizet > 0.0) && autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : double fileID
// Return Type  : void
//
void c_fwrite(double fileID)
{
  FILE * filestar;
  boolean_T autoflush;
  size_t bytesOutSizet;
  static const char x[8] = { 't', ' ', 'u', 'H', 'i', 's', 't', '\x0a' };

  getfilestar(fileID, &filestar, &autoflush);
  if (!(fileID != 0.0)) {
    filestar = NULL;
  }

  if (!(filestar == NULL)) {
    bytesOutSizet = fwrite(&x[0], sizeof(char), (size_t)8, filestar);
    if (((double)bytesOutSizet > 0.0) && autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : double fileID
// Return Type  : void
//
void d_fwrite(double fileID)
{
  FILE * filestar;
  boolean_T autoflush;
  size_t bytesOutSizet;
  static const char x[11] = { 'n', '_', 't', ' ', 'n', '_', 'e', 'l', 'e', 'm',
    '\x0a' };

  getfilestar(fileID, &filestar, &autoflush);
  if (!(fileID != 0.0)) {
    filestar = NULL;
  }

  if (!(filestar == NULL)) {
    bytesOutSizet = fwrite(&x[0], sizeof(char), (size_t)11, filestar);
    if (((double)bytesOutSizet > 0.0) && autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : double fileID
// Return Type  : void
//
void e_fwrite(double fileID)
{
  FILE * filestar;
  boolean_T autoflush;
  size_t bytesOutSizet;
  static const char x[7] = { 't', ' ', 'e', 'l', 'e', 'm', '\x0a' };

  getfilestar(fileID, &filestar, &autoflush);
  if (!(fileID != 0.0)) {
    filestar = NULL;
  }

  if (!(filestar == NULL)) {
    bytesOutSizet = fwrite(&x[0], sizeof(char), (size_t)7, filestar);
    if (((double)bytesOutSizet > 0.0) && autoflush) {
      fflush(filestar);
    }
  }
}

//
// File trailer for fwrite.cpp
//
// [EOF]
//
