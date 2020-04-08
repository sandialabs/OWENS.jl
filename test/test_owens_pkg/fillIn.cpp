//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: fillIn.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "fillIn.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// Arguments    : coder_internal_sparse *b_this
// Return Type  : void
//
void sparse_fillIn(coder_internal_sparse *b_this)
{
  int idx;
  int i;
  int c;
  int ridx;
  int currRowIdx;
  double val;
  idx = 1;
  i = b_this->colidx->size[0];
  for (c = 0; c <= i - 2; c++) {
    ridx = b_this->colidx->data[c];
    b_this->colidx->data[c] = idx;
    while (ridx < b_this->colidx->data[c + 1]) {
      currRowIdx = b_this->rowidx->data[ridx - 1];
      val = b_this->d->data[ridx - 1];
      ridx++;
      if (val != 0.0) {
        b_this->d->data[idx - 1] = val;
        b_this->rowidx->data[idx - 1] = currRowIdx;
        idx++;
      }
    }
  }

  b_this->colidx->data[b_this->colidx->size[0] - 1] = idx;
}

//
// File trailer for fillIn.cpp
//
// [EOF]
//
