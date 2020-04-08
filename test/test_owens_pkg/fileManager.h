//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: fileManager.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:47:29
//
#ifndef FILEMANAGER_H
#define FILEMANAGER_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_owens_types.h"

// Function Declarations
extern signed char b_cfopen(const emxArray_char_T *cfilename, const char
  * cpermission);
extern signed char c_cfopen(const emxArray_char_T *cfilename, const char
  * cpermission);
extern int cfclose(double fid);
extern signed char cfopen(const char * cfilename, const char * cpermission);
extern void filedata_init();
extern void getfilestar(double fid, FILE * *filestar, boolean_T *autoflush);

#endif

//
// File trailer for fileManager.h
//
// [EOF]
//
