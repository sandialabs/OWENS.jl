//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readJointData.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//

// Include Files
#include "readJointData.h"
#include "fileManager.h"
#include "myfgetl.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// readJointData reads joint data file
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [joint] = readJointData(inputfile)
//
//    This function reads the joint data file and stores data in the joint
//    object.
//
//       input:
//       inputfile     = string containing joint data filename
// Arguments    : const emxArray_char_T *inputfile
//                emxArray_real_T *joint
// Return Type  : void
//
void readJointData(const emxArray_char_T *inputfile, emxArray_real_T *joint)
{
  boolean_T b_bool;
  int kstr;
  signed char fileid;
  int fid;
  int exitg1;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  double file_length;
  FILE * b_NULL;
  FILE * filestar;
  int st;
  int i;
  int i1;
  int loop_ub;
  unsigned int count;
  emxArray_char_T *line;
  emxArray_uint32_T *delimiter_idx;
  emxArray_real_T *lineinfo;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  emxArray_char_T *b_line;
  int nx;
  int idx;
  boolean_T exitg2;
  unsigned int u;
  double temp[8];
  unsigned int u1;
  creal_T dc;

  //       output:
  //       joint         = array containing joint data
  b_bool = false;
  if ((inputfile->size[0] == 1) && (inputfile->size[1] == 3)) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (inputfile->data[kstr] != b_cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = b_cfopen(inputfile, "rb");
    fid = fileid;
  }

  // open joint file
  file_length = 0.0;
  b_NULL = NULL;
  do {
    exitg1 = 0;
    getfilestar(static_cast<double>(fid), &filestar, &b_bool);
    if (filestar == b_NULL) {
      kstr = 0;
    } else {
      st = feof(filestar);
      kstr = ((int)st != 0);
    }

    if (kstr == 0) {
      b_myfgetl(static_cast<double>(fid));
      file_length++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  cfclose(static_cast<double>(fid));
  i = static_cast<int>(file_length);
  i1 = joint->size[0] * joint->size[1];
  joint->size[0] = i;
  joint->size[1] = 8;
  emxEnsureCapacity_real_T(joint, i1);
  loop_ub = i << 3;
  for (i = 0; i < loop_ub; i++) {
    joint->data[i] = 0.0;
  }

  b_bool = false;
  if ((inputfile->size[0] == 1) && (inputfile->size[1] == 3)) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (inputfile->data[kstr] != b_cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = b_cfopen(inputfile, "rb");
    fid = fileid;
  }

  // open joint file
  count = 1U;
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  emxInit_uint32_T(&delimiter_idx, 2);
  emxInit_real_T(&lineinfo, 1);
  emxInit_boolean_T(&x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_char_T(&b_line, 2);
  do {
    exitg1 = 0;
    getfilestar(static_cast<double>(fid), &filestar, &b_bool);
    if (filestar == b_NULL) {
      kstr = 0;
    } else {
      st = feof(filestar);
      kstr = ((int)st != 0);
    }

    if (kstr == 0) {
      // read in nodal info associated with joint [joint #, master node #, slave node #, joint type] 
      myfgetl(static_cast<double>(fid), line);

      //  Find where all of the delimiters are
      i = x->size[0] * x->size[1];
      x->size[0] = line->size[1];
      x->size[1] = line->size[0];
      emxEnsureCapacity_boolean_T(x, i);
      loop_ub = line->size[0];
      for (i = 0; i < loop_ub; i++) {
        kstr = line->size[1];
        for (i1 = 0; i1 < kstr; i1++) {
          x->data[i1 + x->size[0] * i] = (line->data[i + line->size[0] * i1] ==
            '\x09');
        }
      }

      nx = x->size[0] * x->size[1];
      idx = 0;
      i = ii->size[0];
      ii->size[0] = nx;
      emxEnsureCapacity_int32_T(ii, i);
      kstr = 0;
      exitg2 = false;
      while ((!exitg2) && (kstr <= nx - 1)) {
        if (x->data[kstr]) {
          idx++;
          ii->data[idx - 1] = kstr + 1;
          if (idx >= nx) {
            exitg2 = true;
          } else {
            kstr++;
          }
        } else {
          kstr++;
        }
      }

      if (nx == 1) {
        if (idx == 0) {
          ii->size[0] = 0;
        }
      } else {
        i = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity_int32_T(ii, i);
      }

      if ((line->size[0] == 0) || (line->size[1] == 0)) {
        kstr = 0;
      } else {
        kstr = line->size[1];
      }

      i = delimiter_idx->size[0] * delimiter_idx->size[1];
      delimiter_idx->size[0] = 1;
      delimiter_idx->size[1] = ii->size[0] + 2;
      emxEnsureCapacity_uint32_T(delimiter_idx, i);
      delimiter_idx->data[0] = 0U;
      loop_ub = ii->size[0];
      for (i = 0; i < loop_ub; i++) {
        delimiter_idx->data[i + 1] = static_cast<unsigned int>(ii->data[i]);
      }

      delimiter_idx->data[ii->size[0] + 1] = kstr + 1U;
      i = lineinfo->size[0];
      lineinfo->size[0] = delimiter_idx->size[1] - 1;
      emxEnsureCapacity_real_T(lineinfo, i);
      loop_ub = delimiter_idx->size[1];
      for (i = 0; i <= loop_ub - 2; i++) {
        lineinfo->data[i] = 0.0;
      }

      //  Extract the data from the beginning to the last delimiter
      i = delimiter_idx->size[1];
      for (kstr = 0; kstr <= i - 2; kstr++) {
        u = delimiter_idx->data[kstr];
        u1 = delimiter_idx->data[kstr + 1];
        if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
          i1 = 0;
          nx = 0;
        } else {
          i1 = static_cast<int>(u);
          nx = static_cast<int>((static_cast<double>(u1) - 1.0));
        }

        idx = b_line->size[0] * b_line->size[1];
        b_line->size[0] = 1;
        loop_ub = nx - i1;
        b_line->size[1] = loop_ub;
        emxEnsureCapacity_char_T(b_line, idx);
        for (nx = 0; nx < loop_ub; nx++) {
          b_line->data[nx] = line->data[i1 + nx];
        }

        dc = c_str2double(b_line);
        lineinfo->data[kstr] = dc.re;
      }

      // reads in mass, stiffness, orientation of element attached to master joint 
      temp[0] = lineinfo->data[0];
      temp[1] = lineinfo->data[1];
      temp[2] = lineinfo->data[2];
      temp[3] = lineinfo->data[3];
      temp[4] = lineinfo->data[4];
      temp[5] = lineinfo->data[5];
      temp[6] = lineinfo->data[6];
      temp[7] = lineinfo->data[7];
      for (i = 0; i < 8; i++) {
        joint->data[(static_cast<int>(count) + joint->size[0] * i) - 1] = temp[i];
      }

      // store data in joint array
      count++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_line);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_real_T(&lineinfo);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);
  cfclose(static_cast<double>(fid));

  // close file
}

//
// File trailer for readJointData.cpp
//
// [EOF]
//
