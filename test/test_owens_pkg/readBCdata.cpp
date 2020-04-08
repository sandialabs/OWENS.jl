//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readBCdata.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "readBCdata.h"
#include "fileManager.h"
#include "myfgetl.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Definitions

//
// readBDdata  reads boundary condition file
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [BC] = readBCdata(bcfilename,numNodes,numDofPerNode)
//
//    This function reads the boundray condition file and stores data in the
//    boundary condition object.
//
//       input:
//       bcfilename    = string containing boundary condition filename
//       numNodes      = number of nodes in structural model
//       numDofPerNode = number of degrees of freedom per node
// Arguments    : const emxArray_char_T *bcfilename
//                double numNodes
//                double *BC_numpBC
//                emxArray_real_T *BC_pBC
//                double *BC_numsBC
//                double *BC_nummBC
//                emxArray_real_T *BC_isConstrained
// Return Type  : void
//
void readBCdata(const emxArray_char_T *bcfilename, double numNodes, double
                *BC_numpBC, emxArray_real_T *BC_pBC, double *BC_numsBC, double
                *BC_nummBC, emxArray_real_T *BC_isConstrained)
{
  boolean_T b_bool;
  int kstr;
  signed char fileid;
  int fid;
  int exitg1;
  emxArray_real_T *pBC;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  emxArray_char_T *b_r;
  creal_T dc;
  int i;
  int i1;
  int loop_ub;
  emxArray_char_T *line;
  emxArray_uint32_T *delimiter_idx;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  emxArray_char_T *b_line;
  int b_i;
  int nx;
  int idx;
  emxArray_int8_T *isConstrained;
  boolean_T exitg2;
  emxArray_real_T *constDof;
  double b_index;
  int j;
  double a;
  int k;
  unsigned int u;
  unsigned int u1;
  double absx;
  int exponent;
  creal_T dc1;

  //       output:
  //       BC            = object containing boundary condition data
  b_bool = false;
  if ((bcfilename->size[0] == 1) && (bcfilename->size[1] == 3)) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (bcfilename->data[kstr] != b_cv[kstr]) {
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
    fileid = b_cfopen(bcfilename, "rb");
    fid = fileid;
  }

  emxInit_real_T(&pBC, 2);
  emxInit_char_T(&b_r, 2);

  // open boundary condition file
  myfgetl(static_cast<double>(fid), b_r);
  dc = d_str2double(b_r);

  // read in number of boundary conditions (displacement boundary conditions)
  i = static_cast<int>(dc.re);
  i1 = pBC->size[0] * pBC->size[1];
  pBC->size[0] = i;
  pBC->size[1] = 3;
  emxEnsureCapacity_real_T(pBC, i1);
  loop_ub = i * 3;
  emxFree_char_T(&b_r);
  for (i1 = 0; i1 < loop_ub; i1++) {
    pBC->data[i1] = 0.0;
  }

  // initialize boundary conditions
  emxInit_char_T(&line, 2);
  emxInit_uint32_T(&delimiter_idx, 2);
  emxInit_boolean_T(&x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_char_T(&b_line, 2);
  for (b_i = 0; b_i < i; b_i++) {
    myfgetl(static_cast<double>(fid), line);

    //  Find where all of the delimiters are
    // first two are boundary condition node number and local DOF number
    // third is boundary condition value (typically zero)
    i1 = x->size[0] * x->size[1];
    x->size[0] = line->size[1];
    x->size[1] = line->size[0];
    emxEnsureCapacity_boolean_T(x, i1);
    loop_ub = line->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      kstr = line->size[1];
      for (idx = 0; idx < kstr; idx++) {
        x->data[idx + x->size[0] * i1] = (line->data[i1 + line->size[0] * idx] ==
          ' ');
      }
    }

    nx = x->size[0] * x->size[1];
    idx = 0;
    i1 = ii->size[0];
    ii->size[0] = nx;
    emxEnsureCapacity_int32_T(ii, i1);
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
      i1 = ii->size[0];
      if (1 > idx) {
        ii->size[0] = 0;
      } else {
        ii->size[0] = idx;
      }

      emxEnsureCapacity_int32_T(ii, i1);
    }

    if ((line->size[0] == 0) || (line->size[1] == 0)) {
      kstr = 0;
    } else {
      kstr = line->size[1];
    }

    i1 = delimiter_idx->size[0] * delimiter_idx->size[1];
    delimiter_idx->size[0] = 1;
    delimiter_idx->size[1] = ii->size[0] + 2;
    emxEnsureCapacity_uint32_T(delimiter_idx, i1);
    delimiter_idx->data[0] = 0U;
    loop_ub = ii->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      delimiter_idx->data[i1 + 1] = static_cast<unsigned int>(ii->data[i1]);
    }

    delimiter_idx->data[ii->size[0] + 1] = kstr + 1U;

    //  Extract the data from the beginning to the last delimiter
    i1 = delimiter_idx->size[1];
    for (k = 0; k <= i1 - 2; k++) {
      u = delimiter_idx->data[k];
      u1 = delimiter_idx->data[k + 1];
      if (static_cast<double>(u) + 1.0 > static_cast<double>(u1) - 1.0) {
        idx = 0;
        kstr = 0;
      } else {
        idx = static_cast<int>(u);
        kstr = static_cast<int>((static_cast<double>(u1) - 1.0));
      }

      nx = b_line->size[0] * b_line->size[1];
      b_line->size[0] = 1;
      loop_ub = kstr - idx;
      b_line->size[1] = loop_ub;
      emxEnsureCapacity_char_T(b_line, nx);
      for (kstr = 0; kstr < loop_ub; kstr++) {
        b_line->data[kstr] = line->data[idx + kstr];
      }

      dc1 = c_str2double(b_line);
      pBC->data[b_i + pBC->size[0] * k] = dc1.re;
    }
  }

  emxFree_char_T(&b_line);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_uint32_T(&delimiter_idx);
  emxFree_char_T(&line);

  // store boundary condition data  in boundayr condition object
  i = BC_pBC->size[0] * BC_pBC->size[1];
  BC_pBC->size[0] = pBC->size[0];
  BC_pBC->size[1] = 3;
  emxEnsureCapacity_real_T(BC_pBC, i);
  loop_ub = pBC->size[0] * pBC->size[1];
  for (i = 0; i < loop_ub; i++) {
    BC_pBC->data[i] = pBC->data[i];
  }

  emxInit_int8_T(&isConstrained, 1);
  cfclose(static_cast<double>(fid));

  // create a vector denoting constrained DOFs in the model (0 unconstrained, 1
  // constrained)
  // calculate constrained dof vector
  loop_ub = static_cast<int>((numNodes * 6.0));
  i = isConstrained->size[0];
  isConstrained->size[0] = loop_ub;
  emxEnsureCapacity_int8_T(isConstrained, i);
  for (i = 0; i < loop_ub; i++) {
    isConstrained->data[i] = 0;
  }

  emxInit_real_T(&constDof, 1);
  loop_ub = pBC->size[0];
  i = constDof->size[0];
  constDof->size[0] = pBC->size[0];
  emxEnsureCapacity_real_T(constDof, i);
  for (i = 0; i < loop_ub; i++) {
    constDof->data[i] = (pBC->data[i] - 1.0) * 6.0 + pBC->data[i + pBC->size[0]];
  }

  emxFree_real_T(&pBC);
  b_index = 1.0;
  i = static_cast<int>(numNodes);
  for (b_i = 0; b_i < i; b_i++) {
    for (j = 0; j < 6; j++) {
      a = ((static_cast<double>(b_i) + 1.0) - 1.0) * 6.0 + (static_cast<double>
        (j) + 1.0);
      b_bool = false;
      k = 0;
      exitg2 = false;
      while ((!exitg2) && (k <= constDof->size[0] - 1)) {
        absx = std::abs(constDof->data[k] / 2.0);
        if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
          if (absx <= 2.2250738585072014E-308) {
            absx = 4.94065645841247E-324;
          } else {
            frexp(absx, &exponent);
            absx = std::ldexp(1.0, exponent - 53);
          }
        } else {
          absx = rtNaN;
        }

        if (std::abs(constDof->data[k] - a) < absx) {
          b_bool = true;
          exitg2 = true;
        } else {
          k++;
        }
      }

      if (b_bool) {
        isConstrained->data[static_cast<int>(b_index) - 1] = 1;
      }

      b_index++;
    }
  }

  emxFree_real_T(&constDof);
  i = BC_isConstrained->size[0];
  BC_isConstrained->size[0] = isConstrained->size[0];
  emxEnsureCapacity_real_T(BC_isConstrained, i);
  loop_ub = isConstrained->size[0];
  for (i = 0; i < loop_ub; i++) {
    BC_isConstrained->data[i] = isConstrained->data[i];
  }

  emxFree_int8_T(&isConstrained);
  *BC_numpBC = dc.re;
  *BC_numsBC = 0.0;
  *BC_nummBC = 0.0;
}

//
// File trailer for readBCdata.cpp
//
// [EOF]
//
