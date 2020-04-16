//
// File: test_owens.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//

// Include Files
#include "test_owens.h"
#include "fileManager.h"
#include "owens.h"
#include "rt_nonfinite.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include "test_owens_initialize.h"
#include "tic.h"
#include <cstring>
#include <stdio.h>
#include <string.h>

// Function Definitions

//
// Arguments    : boolean_T test_transient
//                boolean_T test_modal
//                boolean_T test_flutter
// Return Type  : void
//
void test_owens(boolean_T test_transient, boolean_T test_modal, boolean_T
                test_flutter)
{
  double cmkValues[72];
  int idx;
  int jj;
  static const unsigned int uv[6] = { 9808800U, 9781100U, 18914000U, 3635100000U,
    3650900000U, 2436200000U };

  static const double b_dv[6] = { 132900.0, 132900.0, 1.985E+6,
    2.2878204759573773E+8, 2.2889663915476388E+8, 6.1650258756076582E+7 };

  signed char fileid;
  int ii;
  boolean_T exitg1;
  int i;
  boolean_T guard1 = false;
  int i_data[36];
  FILE * b_NULL;
  signed char j_data[36];
  FILE * filestar;
  boolean_T autoflush;
  double v_data[36];
  double b_dv1[2];
  emxArray_real_T *freq;
  double b_freq;
  double damp;
  emxArray_real_T *b_damp;
  double freq_data[4];
  int freq_size[2];
  double damp_data[4];
  int damp_size[2];
  if (isInitialized_test_owens == false) {
    test_owens_initialize();
  }

  printf("%s\n", "Starting");
  fflush(stdout);
  tic();

  //  ****************** COORDINATE SYSTEM DEFINITION *********************** %
  //  1 - x -  surge (pointed downstream)
  //  2 - y -  sway  (right hand rule)
  //  3 - z -  heave (pointed upwards)
  //  4 - Ox - roll
  //  5 - Oy - pitch
  //  6 - Oz - yaw
  //  *********************************************************************** %
  //  use this benchmark file
  //  append this name to the end of the saved files
  //  convert rotational stiffness to N-m/rad
  //  filename root to save the created nodal file
  //  ************************************************************************
  //  perform the transient simulations using OWENS
  //  *************************************************************************
  //  define the filename saving convention
  //  *********************************************************************
  //  perform operations for the nodal file generation
  //  *********************************************************************
  std::memset(&cmkValues[0], 0, 72U * sizeof(double));
  for (idx = 0; idx < 6; idx++) {
    //  set up mass matrix
    jj = idx + 6 * idx;
    cmkValues[jj] = uv[idx];

    //  set up stiffness matrix
    cmkValues[jj + 36] = b_dv[idx];
  }

  //  writeOwensNDL writes a nodal input file for the OWENS Toolkit
  //  **********************************************************************
  //  *                   Part of SNL VAWTGen                              *
  //  * Developed by Sandia National Laboratories Wind Energy Technologies *
  //  *             See license.txt for disclaimer information             *
  //  **********************************************************************
  //    This function writes a boundary condition file for the OWENS Toolkit
  //       input:
  //       fileRoot     = string containing input prefix of file name
  //       output:     (NONE)
  //  **********************************************************************
  //  % % if extension is used in the filename, remove it
  //  % if exist(strfind(fileRoot,'.'),'var')
  //  %     fileRoot = fileRoot(1:find(fileRoot,'.')-1);
  //  % end
  //  open the BC file to save the boundary conditions to
  // construct file name
  fileid = cfopen("./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.ndl",
                  "wb");

  // open boundary condition file
  //  % % first line should list the number of BCs to follow
  //  % totalBCs = 0;
  //  % for nn = 1:length(nodes)
  //  %     totalBCs = totalBCs + numel(find(cmkValues{nn} ~=0));
  //  % end
  //  % fprintf(fid, '%i\n', totalBCs);
  //  write out the boundary conditions into the file
  idx = 0;
  ii = 1;
  jj = 1;
  exitg1 = false;
  while ((!exitg1) && (jj <= 6)) {
    i = (ii + 6 * (jj - 1)) - 1;
    guard1 = false;
    if (cmkValues[i] != 0.0) {
      idx++;
      i_data[idx - 1] = ii;
      j_data[idx - 1] = static_cast<signed char>(jj);
      v_data[idx - 1] = cmkValues[i];
      if (idx >= 36) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
      if (ii > 6) {
        ii = 1;
        jj++;
      }
    }
  }

  if (1 > idx) {
    idx = 0;
  }

  for (ii = 0; ii < idx; ii++) {
    b_NULL = NULL;
    getfilestar(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%.0f %s %.0f %.0f %.2f\n", 1.0, "M6", static_cast<
              double>(i_data[ii]), static_cast<double>(j_data[ii]), v_data[ii]);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  idx = 0;
  ii = 1;
  jj = 1;
  exitg1 = false;
  while ((!exitg1) && (jj <= 6)) {
    i = (ii + 6 * (jj - 1)) - 1;
    guard1 = false;
    if (cmkValues[i] != 0.0) {
      idx++;
      i_data[idx - 1] = ii;
      j_data[idx - 1] = static_cast<signed char>(jj);
      v_data[idx - 1] = cmkValues[i];
      if (idx >= 36) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
      if (ii > 6) {
        ii = 1;
        jj++;
      }
    }
  }

  if (1 > idx) {
    idx = 0;
  }

  for (ii = 0; ii < idx; ii++) {
    b_NULL = NULL;
    getfilestar(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%.0f %s %.0f %.0f %.2f\n", 1.0, "K6", static_cast<
              double>(i_data[ii]), static_cast<double>(j_data[ii]), v_data[ii]);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  //  close boundary condition file
  cfclose(static_cast<double>(fileid));

  //  *********************************************************************
  //  perform operations for the aerodynamic forces file generation
  //  *********************************************************************
  //  *********************************************************************
  //  run a modal analysis of the platform design
  //  *********************************************************************
  //  rpm
  //  number of rpm stations
  //  number of modes to calculate/extract: TODO: Since we can only use eig, push this change through 
  //  [sec]
  //  length of time vector
  if (test_transient) {
    b_dv1[0] = 0.0;
    b_dv1[1] = 1.1;
    owens(b_dv1, &b_freq, &damp);
  }

  if (test_modal) {
    emxInit_real_T(&freq, 2);
    emxInit_real_T(&b_damp, 2);
    b_owens(freq, b_damp);
    emxFree_real_T(&b_damp);
    emxFree_real_T(&freq);
  }

  if (test_flutter) {
    c_owens(freq_data, freq_size, damp_data, damp_size);
  }

  printf("%s\n", "Function Finished");
  fflush(stdout);
}

//
// File trailer for test_owens.cpp
//
// [EOF]
//
