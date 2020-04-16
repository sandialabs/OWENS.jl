//
// File: fread.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "fread.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : double fileID
//                char A_data[]
//                int A_size[1]
// Return Type  : void
//
void b_fread(double fileID, char A_data[], int A_size[1])
{
  FILE * filestar;
  boolean_T a;
  size_t numReadSizeT;
  char b_A_data[1];
  getfilestar(fileID, &filestar, &a);
  if ((!(fileID != 0.0)) || (!(fileID != 1.0)) || (!(fileID != 2.0))) {
    filestar = NULL;
  }

  if (filestar == NULL) {
    A_size[0] = 0;
  } else {
    numReadSizeT = fread(&b_A_data[0], sizeof(char), 1, filestar);
    if ((int)numReadSizeT + 1 <= 1) {
      b_A_data[0] = '\x00';
    }

    A_size[0] = 1;
    A_data[0] = b_A_data[0];
    if ((int)numReadSizeT < 1) {
      A_size[0] = 0;
    }
  }
}

//
// File trailer for fread.cpp
//
// [EOF]
//
