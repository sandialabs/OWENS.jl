//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: processAeroLoadsBLE.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "processAeroLoadsBLE.h"
#include "mapCactusLoadsFile.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// Arguments    : const emxArray_char_T *aeroLoadsFile
//                const emxArray_char_T *OWENSfile
//                double aeroLoads_timeArray_data[]
//                int aeroLoads_timeArray_size[1]
//                emxArray_real_T *aeroLoads_ForceValHist
//                emxArray_real_T *aeroLoads_ForceDof
// Return Type  : void
//
void processAeroLoadsBLE(const emxArray_char_T *aeroLoadsFile, const
  emxArray_char_T *OWENSfile, double aeroLoads_timeArray_data[], int
  aeroLoads_timeArray_size[1], emxArray_real_T *aeroLoads_ForceValHist,
  emxArray_real_T *aeroLoads_ForceDof)
{
  emxArray_char_T *b_aeroLoadsFile;
  int loop_ub;
  int b_loop_ub;
  int i;
  emxArray_char_T *c_aeroLoadsFile;
  static const char b_cv[5] = { '.', 'g', 'e', 'o', 'm' };

  emxArray_char_T *b_OWENSfile;
  static const char cv1[16] = { '_', 'E', 'l', 'e', 'm', 'e', 'n', 't', 'D', 'a',
    't', 'a', '.', 'c', 's', 'v' };

  emxArray_char_T *c_OWENSfile;
  emxArray_char_T *d_OWENSfile;
  emxArray_char_T *e_OWENSfile;
  static const char cv2[5] = { '.', 'm', 'e', 's', 'h' };

  emxInit_char_T(&b_aeroLoadsFile, 2);
  if (1 > aeroLoadsFile->size[1] - 16) {
    loop_ub = 0;
  } else {
    loop_ub = aeroLoadsFile->size[1] - 16;
  }

  // cut off the _ElementData.csv
  if (1 > OWENSfile->size[1] - 6) {
    b_loop_ub = 0;
  } else {
    b_loop_ub = OWENSfile->size[1] - 6;
  }

  // cut off the .owens
  i = b_aeroLoadsFile->size[0] * b_aeroLoadsFile->size[1];
  b_aeroLoadsFile->size[0] = 1;
  b_aeroLoadsFile->size[1] = loop_ub + 5;
  emxEnsureCapacity_char_T(b_aeroLoadsFile, i);
  for (i = 0; i < loop_ub; i++) {
    b_aeroLoadsFile->data[i] = aeroLoadsFile->data[i];
  }

  for (i = 0; i < 5; i++) {
    b_aeroLoadsFile->data[i + loop_ub] = b_cv[i];
  }

  emxInit_char_T(&c_aeroLoadsFile, 2);
  i = c_aeroLoadsFile->size[0] * c_aeroLoadsFile->size[1];
  c_aeroLoadsFile->size[0] = 1;
  c_aeroLoadsFile->size[1] = loop_ub + 16;
  emxEnsureCapacity_char_T(c_aeroLoadsFile, i);
  for (i = 0; i < loop_ub; i++) {
    c_aeroLoadsFile->data[i] = aeroLoadsFile->data[i];
  }

  for (i = 0; i < 16; i++) {
    c_aeroLoadsFile->data[i + loop_ub] = cv1[i];
  }

  emxInit_char_T(&b_OWENSfile, 2);
  i = b_OWENSfile->size[0] * b_OWENSfile->size[1];
  b_OWENSfile->size[0] = 1;
  b_OWENSfile->size[1] = b_loop_ub + 4;
  emxEnsureCapacity_char_T(b_OWENSfile, i);
  for (i = 0; i < b_loop_ub; i++) {
    b_OWENSfile->data[i] = OWENSfile->data[i];
  }

  emxInit_char_T(&c_OWENSfile, 2);
  b_OWENSfile->data[b_loop_ub] = '.';
  b_OWENSfile->data[b_loop_ub + 1] = 'b';
  b_OWENSfile->data[b_loop_ub + 2] = 'l';
  b_OWENSfile->data[b_loop_ub + 3] = 'd';
  i = c_OWENSfile->size[0] * c_OWENSfile->size[1];
  c_OWENSfile->size[0] = 1;
  c_OWENSfile->size[1] = b_loop_ub + 3;
  emxEnsureCapacity_char_T(c_OWENSfile, i);
  for (i = 0; i < b_loop_ub; i++) {
    c_OWENSfile->data[i] = OWENSfile->data[i];
  }

  emxInit_char_T(&d_OWENSfile, 2);
  c_OWENSfile->data[b_loop_ub] = '.';
  c_OWENSfile->data[b_loop_ub + 1] = 'e';
  c_OWENSfile->data[b_loop_ub + 2] = 'l';
  i = d_OWENSfile->size[0] * d_OWENSfile->size[1];
  d_OWENSfile->size[0] = 1;
  d_OWENSfile->size[1] = b_loop_ub + 4;
  emxEnsureCapacity_char_T(d_OWENSfile, i);
  for (i = 0; i < b_loop_ub; i++) {
    d_OWENSfile->data[i] = OWENSfile->data[i];
  }

  emxInit_char_T(&e_OWENSfile, 2);
  d_OWENSfile->data[b_loop_ub] = '.';
  d_OWENSfile->data[b_loop_ub + 1] = 'o';
  d_OWENSfile->data[b_loop_ub + 2] = 'r';
  d_OWENSfile->data[b_loop_ub + 3] = 't';
  i = e_OWENSfile->size[0] * e_OWENSfile->size[1];
  e_OWENSfile->size[0] = 1;
  e_OWENSfile->size[1] = b_loop_ub + 5;
  emxEnsureCapacity_char_T(e_OWENSfile, i);
  for (i = 0; i < b_loop_ub; i++) {
    e_OWENSfile->data[i] = OWENSfile->data[i];
  }

  for (i = 0; i < 5; i++) {
    e_OWENSfile->data[i + b_loop_ub] = cv2[i];
  }

  mapCactusLoadsFile(b_aeroLoadsFile, c_aeroLoadsFile, b_OWENSfile, c_OWENSfile,
                     d_OWENSfile, e_OWENSfile, aeroLoads_timeArray_data,
                     aeroLoads_timeArray_size, aeroLoads_ForceValHist,
                     aeroLoads_ForceDof);

  //  method constrained by not passing the .mat filename
  //  save('aeroLoads.mat','timeArray','ForceValHist','ForceDof');
  //  disp('New aeroLoads.mat file saved.')
  emxFree_char_T(&e_OWENSfile);
  emxFree_char_T(&d_OWENSfile);
  emxFree_char_T(&c_OWENSfile);
  emxFree_char_T(&b_OWENSfile);
  emxFree_char_T(&c_aeroLoadsFile);
  emxFree_char_T(&b_aeroLoadsFile);
}

//
// File trailer for processAeroLoadsBLE.cpp
//
// [EOF]
//
