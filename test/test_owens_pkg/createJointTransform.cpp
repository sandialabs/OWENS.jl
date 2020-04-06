//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: createJointTransform.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "createJointTransform.h"
#include "calculateLambda.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <cstring>
#include <math.h>
#include <string.h>

// Variable Definitions
static const signed char iv[36] = { -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
  -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1 };

// Function Declarations
static void b_getNodeMaps(double masterNodeNum, double slaveNodeNum, const
  double slaveDof[3], const double activeDof[6], const double
  slaveActiveDof_data[], const int slaveActiveDof_size[2], double dDOF[3],
  double aDOF_data[], int aDOF_size[1]);
static void c_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2]);
static void c_getNodeMaps(const double Rdd[25], const double Rda[30], double
  masterNodeNum, double slaveNodeNum, const double slaveDof[5], const double
  activeDof[6], const double slaveActiveDof_data[], const int
  slaveActiveDof_size[2], double Tda[30], double dDOF[5], double aDOF_data[],
  int aDOF_size[1]);
static void createTda(double jointType, double slaveNodeNum, double
                      masterNodeNum, double psi, double theta, const double
                      joint[8], double Tda_data[], int Tda_size[2], double
                      dDOF_data[], int dDOF_size[1], double aDOF_data[], int
                      aDOF_size[1]);
static void d_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2]);
static void d_getNodeMaps(const double Rda[36], double masterNodeNum, double
  slaveNodeNum, const double slaveDof[6], const double activeDof[6], const
  double slaveActiveDof_data[], const int slaveActiveDof_size[2], double Tda[36],
  double dDOF[6], double aDOF_data[], int aDOF_size[1]);
static void e_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2]);
static void f_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2]);
static void g_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2]);
static void getNodeMaps(double masterNodeNum, double slaveNodeNum, const double
  slaveDof[6], const double activeDof[6], const double slaveActiveDof_data[],
  const int slaveActiveDof_size[2], double dDOF[6], double aDOF_data[], int
  aDOF_size[1]);

// Function Definitions

//
// Arguments    : double masterNodeNum
//                double slaveNodeNum
//                const double slaveDof[3]
//                const double activeDof[6]
//                const double slaveActiveDof_data[]
//                const int slaveActiveDof_size[2]
//                double dDOF[3]
//                double aDOF_data[]
//                int aDOF_size[1]
// Return Type  : void
//
static void b_getNodeMaps(double masterNodeNum, double slaveNodeNum, const
  double slaveDof[3], const double activeDof[6], const double
  slaveActiveDof_data[], const int slaveActiveDof_size[2], double dDOF[3],
  double aDOF_data[], int aDOF_size[1])
{
  double dDOF_tmp;
  int i;
  double aMap[6];
  double aMap2_data[1];
  int b_i;

  // calculate Tda
  // get number of joint DOFs for this joint
  // get number of active DOFs for this joint
  // initialize arrays
  // get global DOF numbers of slave DOFs for this joint
  dDOF_tmp = (slaveNodeNum - 1.0) * 6.0;
  dDOF[0] = dDOF_tmp + slaveDof[0];

  // get global DOF numbers of slave DOFs for this joint
  dDOF[1] = dDOF_tmp + slaveDof[1];

  // get global DOF numbers of slave DOFs for this joint
  dDOF[2] = dDOF_tmp + slaveDof[2];
  for (i = 0; i < 6; i++) {
    // get global DOF numbers of active DOFs for this joint from master nodes
    aMap[i] = (masterNodeNum - 1.0) * 6.0 + activeDof[i];
  }

  // determine global active DOFs associated with slave node
  if ((slaveActiveDof_size[0] == 0) || (slaveActiveDof_size[1] == 0)) {
    i = 0;
  } else {
    i = 1;
  }

  if (0 <= i - 1) {
    aMap2_data[0] = 0.0;
  }

  if ((slaveActiveDof_size[0] != 0) && (slaveActiveDof_size[1] != 0)) {
    aMap2_data[0] = dDOF_tmp + slaveActiveDof_data[0];
  }

  if (i != 0) {
    // create overall map of active DOFs associated with this joint
    aDOF_size[0] = i + 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }

    if (0 <= i - 1) {
      aDOF_data[6] = aMap2_data[0];
    }
  } else {
    aDOF_size[0] = 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }
  }
}

//
// This function determines the local master DOF associated with a local slave DOF.
//  Get size
// Arguments    : double slaveNodeActiveDof_data[]
//                int slaveNodeActiveDof_size[2]
// Return Type  : void
//
static void c_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2])
{
  double count;
  int i;
  int idx;
  double absx;
  boolean_T tf[6];
  int exponent;
  int ii;
  boolean_T exitg1;
  int b_exponent;
  count = 1.0;
  for (i = 0; i < 6; i++) {
    // loop over number of DOF per node
    for (idx = 0; idx < 6; idx++) {
      tf[idx] = false;
    }

    absx = (static_cast<double>(i) + 1.0) / 2.0;
    for (idx = 0; idx < 6; idx++) {
      frexp(absx, &exponent);
      if (std::abs(static_cast<double>((i - idx))) < std::ldexp(1.0, exponent -
           53)) {
        tf[idx] = true;
      }
    }

    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 6)) {
      if (tf[ii]) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
      //          slaveNodeActiveDof(count) = i;
      count++;
    }
  }

  if (count > 1.0) {
    slaveNodeActiveDof_size[0] = 1;
    slaveNodeActiveDof_size[1] = 1;
    slaveNodeActiveDof_data[0] = 0.0;
    for (i = 0; i < 6; i++) {
      // loop over number of DOF per node
      for (idx = 0; idx < 6; idx++) {
        tf[idx] = false;
      }

      absx = (static_cast<double>(i) + 1.0) / 2.0;
      for (idx = 0; idx < 6; idx++) {
        frexp(absx, &b_exponent);
        if (std::abs(static_cast<double>((i - idx))) < std::ldexp(1.0,
             b_exponent - 53)) {
          tf[idx] = true;
        }
      }

      idx = 0;
      ii = 0;
      exitg1 = false;
      while ((!exitg1) && (ii < 6)) {
        if (tf[ii]) {
          idx = 1;
          exitg1 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
        slaveNodeActiveDof_data[0] = static_cast<double>(i) + 1.0;
      }
    }
  } else {
    slaveNodeActiveDof_size[0] = 0;
    slaveNodeActiveDof_size[1] = 0;
  }
}

//
// Arguments    : const double Rdd[25]
//                const double Rda[30]
//                double masterNodeNum
//                double slaveNodeNum
//                const double slaveDof[5]
//                const double activeDof[6]
//                const double slaveActiveDof_data[]
//                const int slaveActiveDof_size[2]
//                double Tda[30]
//                double dDOF[5]
//                double aDOF_data[]
//                int aDOF_size[1]
// Return Type  : void
//
static void c_getNodeMaps(const double Rdd[25], const double Rda[30], double
  masterNodeNum, double slaveNodeNum, const double slaveDof[5], const double
  activeDof[6], const double slaveActiveDof_data[], const int
  slaveActiveDof_size[2], double Tda[30], double dDOF[5], double aDOF_data[],
  int aDOF_size[1])
{
  int i;
  double y[25];
  double x[25];
  int j;
  signed char ipiv[5];
  int mmj_tmp;
  int b;
  signed char p[5];
  int jj;
  int jy;
  int jp1j;
  int iy;
  int ix;
  double smax;
  int k;
  double s;
  int i1;
  double aMap[6];
  double aMap2_data[1];
  for (i = 0; i < 25; i++) {
    y[i] = 0.0;
    x[i] = Rdd[i];
  }

  for (i = 0; i < 5; i++) {
    ipiv[i] = static_cast<signed char>((i + 1));
  }

  for (j = 0; j < 4; j++) {
    mmj_tmp = 3 - j;
    b = j * 6;
    jj = j * 6;
    jp1j = b + 2;
    iy = 5 - j;
    jy = 0;
    ix = b;
    smax = std::abs(x[jj]);
    for (k = 2; k <= iy; k++) {
      ix++;
      s = std::abs(x[ix]);
      if (s > smax) {
        jy = k - 1;
        smax = s;
      }
    }

    if (x[jj + jy] != 0.0) {
      if (jy != 0) {
        iy = j + jy;
        ipiv[j] = static_cast<signed char>((iy + 1));
        ix = j;
        for (k = 0; k < 5; k++) {
          smax = x[ix];
          x[ix] = x[iy];
          x[iy] = smax;
          ix += 5;
          iy += 5;
        }
      }

      i = (jj - j) + 5;
      for (ix = jp1j; ix <= i; ix++) {
        x[ix - 1] /= x[jj];
      }
    }

    jy = b + 5;
    iy = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = x[jy];
      if (x[jy] != 0.0) {
        ix = jj + 1;
        i = iy + 7;
        i1 = (iy - j) + 10;
        for (jp1j = i; jp1j <= i1; jp1j++) {
          x[jp1j - 1] += x[ix] * -smax;
          ix++;
        }
      }

      jy += 5;
      iy += 5;
    }
  }

  for (i = 0; i < 5; i++) {
    p[i] = static_cast<signed char>((i + 1));
  }

  if (ipiv[0] > 1) {
    jy = ipiv[0] - 1;
    iy = p[jy];
    p[jy] = p[0];
    p[0] = static_cast<signed char>(iy);
  }

  if (ipiv[1] > 2) {
    jy = ipiv[1] - 1;
    iy = p[jy];
    p[jy] = p[1];
    p[1] = static_cast<signed char>(iy);
  }

  if (ipiv[2] > 3) {
    jy = ipiv[2] - 1;
    iy = p[jy];
    p[jy] = p[2];
    p[2] = static_cast<signed char>(iy);
  }

  if (ipiv[3] > 4) {
    jy = ipiv[3] - 1;
    iy = p[jy];
    p[jy] = p[3];
    p[3] = static_cast<signed char>(iy);
  }

  for (k = 0; k < 5; k++) {
    b = 5 * (p[k] - 1);
    y[k + b] = 1.0;
    for (j = k + 1; j < 6; j++) {
      i = (j + b) - 1;
      if (y[i] != 0.0) {
        i1 = j + 1;
        for (ix = i1; ix < 6; ix++) {
          jy = (ix + b) - 1;
          y[jy] -= y[i] * x[(ix + 5 * (j - 1)) - 1];
        }
      }
    }
  }

  for (j = 0; j < 5; j++) {
    jy = 5 * j;
    for (k = 4; k >= 0; k--) {
      iy = 5 * k;
      i = k + jy;
      if (y[i] != 0.0) {
        y[i] /= x[k + iy];
        for (ix = 0; ix < k; ix++) {
          b = ix + jy;
          y[b] -= y[i] * x[ix + iy];
        }
      }
    }
  }

  for (i = 0; i < 25; i++) {
    y[i] = -y[i];
  }

  // calculate Tda
  // get number of joint DOFs for this joint
  // get number of active DOFs for this joint
  // initialize arrays
  for (ix = 0; ix < 5; ix++) {
    for (i = 0; i < 6; i++) {
      smax = 0.0;
      for (i1 = 0; i1 < 5; i1++) {
        smax += y[ix + 5 * i1] * Rda[i1 + 5 * i];
      }

      Tda[ix + 5 * i] = smax;
    }

    // get global DOF numbers of slave DOFs for this joint
    dDOF[ix] = (slaveNodeNum - 1.0) * 6.0 + slaveDof[ix];
  }

  for (ix = 0; ix < 6; ix++) {
    // get global DOF numbers of active DOFs for this joint from master nodes
    aMap[ix] = (masterNodeNum - 1.0) * 6.0 + activeDof[ix];
  }

  // determine global active DOFs associated with slave node
  if ((slaveActiveDof_size[0] == 0) || (slaveActiveDof_size[1] == 0)) {
    iy = 0;
  } else {
    iy = 1;
  }

  if (0 <= iy - 1) {
    aMap2_data[0] = 0.0;
  }

  if ((slaveActiveDof_size[0] != 0) && (slaveActiveDof_size[1] != 0)) {
    aMap2_data[0] = (slaveNodeNum - 1.0) * 6.0 + slaveActiveDof_data[0];
  }

  if (iy != 0) {
    // create overall map of active DOFs associated with this joint
    aDOF_size[0] = iy + 6;
    for (i = 0; i < 6; i++) {
      aDOF_data[i] = aMap[i];
    }

    if (0 <= iy - 1) {
      aDOF_data[6] = aMap2_data[0];
    }
  } else {
    aDOF_size[0] = 6;
    for (i = 0; i < 6; i++) {
      aDOF_data[i] = aMap[i];
    }
  }
}

//
// This function creates a constraint transformation matrix for a single
// joint. Tda is this matrix, dDOF contains a listing of dependent global
// DOFs associated with this joint, and aDOF contains a listing of active
// global DOFs associated with this joint.
// Arguments    : double jointType
//                double slaveNodeNum
//                double masterNodeNum
//                double psi
//                double theta
//                const double joint[8]
//                double Tda_data[]
//                int Tda_size[2]
//                double dDOF_data[]
//                int dDOF_size[1]
//                double aDOF_data[]
//                int aDOF_size[1]
// Return Type  : void
//
static void createTda(double jointType, double slaveNodeNum, double
                      masterNodeNum, double psi, double theta, const double
                      joint[8], double Tda_data[], int Tda_size[2], double
                      dDOF_data[], int dDOF_size[1], double aDOF_data[], int
                      aDOF_size[1])
{
  double d;
  double Lambda[144];
  double slaveActiveDof5_data[1];
  int slaveActiveDof5_size[2];
  int i;
  double activeDof3[6];
  double b_dv[3];
  double b_dv1[6];
  double dv2[6];
  int b_i;
  double slaveDof3[5];
  static const signed char b_iv[5] = { 1, 2, 3, 4, 6 };

  static const signed char b_iv1[5] = { 1, 2, 3, 5, 6 };

  int Rda4_tmp;
  double dDOF[3];
  double Rda5[36];
  signed char b_I[9];
  static const signed char Tda[36] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  double globalConstraintEqMatrix4[60];
  static const signed char a[60] = { -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };

  double Rda4[30];
  static const signed char b_Tda[27] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const signed char b_a[60] = { -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  signed char i1;
  double Rdd4[25];
  double c_Tda[30];
  double b_dDOF[5];
  signed char c_I[9];
  static const signed char c_a[60] = { -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 };

  int globalConstraintEqMatrix4_tmp;
  static const signed char b_iv2[24] = { 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 };

  double d_Tda[36];
  static const signed char b_iv3[24] = { 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 };

  static const signed char iv4[24] = { 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

  if (jointType == 4.0) {
    d = std::abs(theta);
    if ((std::abs(d - 90.0) < 0.001) || (std::abs(d - 270.0) < 0.001)) {
      theta = 0.0;
      jointType = 3.0;
    }
  }

  // calculate transformation matrix from hub frame to joint frame
  calculateLambda(psi * 3.1415926535897931 / 180.0, theta * 3.1415926535897931 /
                  180.0, 0.0, Lambda);

  // Tda is a local mapping of dependent DOFs to active DOFs at a node
  //  u_d = Tda * u_a
  //  such that u_d is a list of local dependent slave DOFs at a jont and
  //  u_a is a list of local dependent slave DOFs at a joint.
  if (jointType == 0.0) {
    // for weld/fixed joint type
    // active DOF list at joint
    // slave DOF list at joint
    // determine local active DOFs associated with slave node
    c_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);

    // from constraint equation for fixed joint
    for (i = 0; i < 6; i++) {
      b_dv1[i] = static_cast<double>(i) + 1.0;
      dv2[i] = static_cast<double>(i) + 1.0;
    }

    getNodeMaps(masterNodeNum, slaveNodeNum, b_dv1, dv2, slaveActiveDof5_data,
                slaveActiveDof5_size, activeDof3, aDOF_data, aDOF_size);
    Tda_size[0] = 6;
    Tda_size[1] = 6;
    for (i = 0; i < 36; i++) {
      Tda_data[i] = Tda[i];
    }

    dDOF_size[0] = 6;
    for (i = 0; i < 6; i++) {
      dDOF_data[i] = activeDof3[i];
    }
  } else if (jointType == 1.0) {
    // for pinned joint type
    // active DOF list at joint
    // slave DOF list at joint
    // determine local active DOFs associated with slave node
    d_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);

    // from constraint equation for pinned joint
    b_dv[0] = 1.0;
    b_dv[1] = 2.0;
    b_dv[2] = 3.0;
    for (i = 0; i < 6; i++) {
      b_dv1[i] = static_cast<double>(i) + 1.0;
    }

    b_getNodeMaps(masterNodeNum, slaveNodeNum, b_dv, b_dv1, slaveActiveDof5_data,
                  slaveActiveDof5_size, dDOF, aDOF_data, aDOF_size);
    Tda_size[0] = 3;
    Tda_size[1] = 9;
    for (i = 0; i < 27; i++) {
      Tda_data[i] = b_Tda[i];
    }

    dDOF_size[0] = 3;
    dDOF_data[0] = dDOF[0];
    dDOF_data[1] = dDOF[1];
    dDOF_data[2] = dDOF[2];
  } else if (jointType == 2.0) {
    // hinge axis along localy "2" frame of joint
    d = std::abs(psi);
    if ((std::abs(d - 90.0) < 0.001) || (std::abs(d - 270.0) < 0.001)) {
      for (i = 0; i < 6; i++) {
        activeDof3[i] = static_cast<double>(i) + 1.0;
      }

      for (i = 0; i < 5; i++) {
        slaveDof3[i] = b_iv1[i];
      }

      // determine local active DOFs associated with slave node
      e_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
      for (i = 0; i < 9; i++) {
        b_I[i] = 0;
      }

      b_I[0] = 1;
      b_I[4] = 1;
      b_I[8] = 1;
      for (i = 0; i < 3; i++) {
        i1 = b_I[3 * i];
        globalConstraintEqMatrix4[5 * i] = -static_cast<double>(i1);
        globalConstraintEqMatrix4_tmp = 5 * (i + 3);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp] = 0.0;
        b_i = 5 * (i + 6);
        globalConstraintEqMatrix4[b_i] = i1;
        Rda4_tmp = 5 * (i + 9);
        globalConstraintEqMatrix4[Rda4_tmp] = 0.0;
        i1 = b_I[3 * i + 1];
        globalConstraintEqMatrix4[5 * i + 1] = -static_cast<double>(i1);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 1] = 0.0;
        globalConstraintEqMatrix4[b_i + 1] = i1;
        globalConstraintEqMatrix4[Rda4_tmp + 1] = 0.0;
        i1 = b_I[3 * i + 2];
        globalConstraintEqMatrix4[5 * i + 2] = -static_cast<double>(i1);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 2] = 0.0;
        globalConstraintEqMatrix4[b_i + 2] = i1;
        globalConstraintEqMatrix4[Rda4_tmp + 2] = 0.0;
      }

      for (i = 0; i < 12; i++) {
        globalConstraintEqMatrix4_tmp = i << 1;
        globalConstraintEqMatrix4[5 * i + 3] =
          b_iv2[globalConstraintEqMatrix4_tmp];
        globalConstraintEqMatrix4[5 * i + 4] =
          b_iv2[globalConstraintEqMatrix4_tmp + 1];
      }
    } else {
      for (i = 0; i < 6; i++) {
        activeDof3[i] = static_cast<double>(i) + 1.0;
      }

      for (i = 0; i < 5; i++) {
        slaveDof3[i] = b_iv[i];
      }

      // determine local active DOFs associated with slave node
      f_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
      for (i = 0; i < 5; i++) {
        for (b_i = 0; b_i < 12; b_i++) {
          d = 0.0;
          for (Rda4_tmp = 0; Rda4_tmp < 12; Rda4_tmp++) {
            d += static_cast<double>(b_a[i + 5 * Rda4_tmp]) * Lambda[Rda4_tmp +
              12 * b_i];
          }

          globalConstraintEqMatrix4[i + 5 * b_i] = d;
        }
      }
    }

    // extract Rda from globalConstraintEqMatrix2
    for (b_i = 0; b_i < 6; b_i++) {
      for (i = 0; i < 5; i++) {
        Rda4[i + 5 * b_i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          (activeDof3[b_i]) - 1)];
      }
    }

    if ((slaveActiveDof5_size[0] != 0) && (slaveActiveDof5_size[1] != 0)) {
      for (i = 0; i < 5; i++) {
        Rda4[i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          ((slaveActiveDof5_data[0] + 6.0)) - 1)];
      }
    }

    // extract Rdd from globalConstraintEqMatrix2
    for (b_i = 0; b_i < 5; b_i++) {
      for (i = 0; i < 5; i++) {
        Rdd4[i + 5 * b_i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          (slaveDof3[b_i]) + 5)];
      }
    }

    c_getNodeMaps(Rdd4, Rda4, masterNodeNum, slaveNodeNum, slaveDof3, activeDof3,
                  slaveActiveDof5_data, slaveActiveDof5_size, c_Tda, b_dDOF,
                  aDOF_data, aDOF_size);
    Tda_size[0] = 5;
    Tda_size[1] = 6;
    std::memcpy(&Tda_data[0], &c_Tda[0], 30U * sizeof(double));
    dDOF_size[0] = 5;
    for (i = 0; i < 5; i++) {
      dDOF_data[i] = b_dDOF[i];
    }
  } else if (jointType == 3.0) {
    // hinge axis along local "1" frame of joint
    d = std::abs(theta);
    if ((std::abs(d - 90.0) < 0.001) || (std::abs(d - 270.0) < 0.001)) {
      for (i = 0; i < 6; i++) {
        activeDof3[i] = static_cast<double>(i) + 1.0;
      }

      for (i = 0; i < 5; i++) {
        slaveDof3[i] = static_cast<double>(i) + 1.0;
      }

      // determine local active DOFs associated with slave node
      g_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
      for (i = 0; i < 9; i++) {
        b_I[i] = 0;
      }

      b_I[0] = 1;
      b_I[4] = 1;
      b_I[8] = 1;
      for (i = 0; i < 9; i++) {
        c_I[i] = 0;
      }

      for (b_i = 0; b_i < 3; b_i++) {
        c_I[b_i + 3 * b_i] = 1;
        globalConstraintEqMatrix4[5 * b_i] = -static_cast<double>(b_I[3 * b_i]);
        globalConstraintEqMatrix4_tmp = 5 * (b_i + 3);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp] = 0.0;
        globalConstraintEqMatrix4[5 * b_i + 1] = -static_cast<double>(b_I[3 *
          b_i + 1]);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 1] = 0.0;
        globalConstraintEqMatrix4[5 * b_i + 2] = -static_cast<double>(b_I[3 *
          b_i + 2]);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 2] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        globalConstraintEqMatrix4_tmp = 5 * (i + 6);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp] = c_I[3 * i];
        b_i = 5 * (i + 9);
        globalConstraintEqMatrix4[b_i] = 0.0;
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 1] = c_I[3 * i
          + 1];
        globalConstraintEqMatrix4[b_i + 1] = 0.0;
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 2] = c_I[3 * i
          + 2];
        globalConstraintEqMatrix4[b_i + 2] = 0.0;
      }

      for (i = 0; i < 12; i++) {
        globalConstraintEqMatrix4_tmp = i << 1;
        globalConstraintEqMatrix4[5 * i + 3] =
          b_iv3[globalConstraintEqMatrix4_tmp];
        globalConstraintEqMatrix4[5 * i + 4] =
          b_iv3[globalConstraintEqMatrix4_tmp + 1];
      }
    } else if ((std::abs(std::abs(psi) - 90.0) < 0.001) || (std::abs(std::abs
                 (psi) - 270.0) < 0.001)) {
      for (i = 0; i < 6; i++) {
        activeDof3[i] = static_cast<double>(i) + 1.0;
      }

      for (i = 0; i < 5; i++) {
        slaveDof3[i] = b_iv[i];
      }

      // determine local active DOFs associated with slave node
      f_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
      for (i = 0; i < 9; i++) {
        b_I[i] = 0;
      }

      b_I[0] = 1;
      b_I[4] = 1;
      b_I[8] = 1;
      for (i = 0; i < 9; i++) {
        c_I[i] = 0;
      }

      for (b_i = 0; b_i < 3; b_i++) {
        c_I[b_i + 3 * b_i] = 1;
        globalConstraintEqMatrix4[5 * b_i] = -static_cast<double>(b_I[3 * b_i]);
        globalConstraintEqMatrix4_tmp = 5 * (b_i + 3);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp] = 0.0;
        globalConstraintEqMatrix4[5 * b_i + 1] = -static_cast<double>(b_I[3 *
          b_i + 1]);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 1] = 0.0;
        globalConstraintEqMatrix4[5 * b_i + 2] = -static_cast<double>(b_I[3 *
          b_i + 2]);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 2] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        globalConstraintEqMatrix4_tmp = 5 * (i + 6);
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp] = c_I[3 * i];
        b_i = 5 * (i + 9);
        globalConstraintEqMatrix4[b_i] = 0.0;
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 1] = c_I[3 * i
          + 1];
        globalConstraintEqMatrix4[b_i + 1] = 0.0;
        globalConstraintEqMatrix4[globalConstraintEqMatrix4_tmp + 2] = c_I[3 * i
          + 2];
        globalConstraintEqMatrix4[b_i + 2] = 0.0;
      }

      for (i = 0; i < 12; i++) {
        globalConstraintEqMatrix4_tmp = i << 1;
        globalConstraintEqMatrix4[5 * i + 3] = iv4[globalConstraintEqMatrix4_tmp];
        globalConstraintEqMatrix4[5 * i + 4] = iv4[globalConstraintEqMatrix4_tmp
          + 1];
      }
    } else {
      for (i = 0; i < 6; i++) {
        activeDof3[i] = static_cast<double>(i) + 1.0;
      }

      for (i = 0; i < 5; i++) {
        slaveDof3[i] = b_iv1[i];
      }

      // determine local active DOFs associated with slave node
      e_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
      for (i = 0; i < 5; i++) {
        for (b_i = 0; b_i < 12; b_i++) {
          d = 0.0;
          for (Rda4_tmp = 0; Rda4_tmp < 12; Rda4_tmp++) {
            d += static_cast<double>(c_a[i + 5 * Rda4_tmp]) * Lambda[Rda4_tmp +
              12 * b_i];
          }

          globalConstraintEqMatrix4[i + 5 * b_i] = d;
        }
      }
    }

    // extract Rda from globalConstraintEqMatrix3
    for (b_i = 0; b_i < 6; b_i++) {
      for (i = 0; i < 5; i++) {
        Rda4[i + 5 * b_i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          (activeDof3[b_i]) - 1)];
      }
    }

    if ((slaveActiveDof5_size[0] != 0) && (slaveActiveDof5_size[1] != 0)) {
      for (i = 0; i < 5; i++) {
        Rda4[i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          ((slaveActiveDof5_data[0] + 6.0)) - 1)];
      }
    }

    // extract Rdd from globalConstraintEqMatrix3
    for (b_i = 0; b_i < 5; b_i++) {
      for (i = 0; i < 5; i++) {
        Rdd4[i + 5 * b_i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          (slaveDof3[b_i]) + 5)];
      }
    }

    c_getNodeMaps(Rdd4, Rda4, masterNodeNum, slaveNodeNum, slaveDof3, activeDof3,
                  slaveActiveDof5_data, slaveActiveDof5_size, c_Tda, b_dDOF,
                  aDOF_data, aDOF_size);
    Tda_size[0] = 5;
    Tda_size[1] = 6;
    std::memcpy(&Tda_data[0], &c_Tda[0], 30U * sizeof(double));
    dDOF_size[0] = 5;
    for (i = 0; i < 5; i++) {
      dDOF_data[i] = b_dDOF[i];
    }
  } else if (jointType == 4.0) {
    // hinge axis along local "3" frame of joint
    // determine local active DOFs associated with slave node
    g_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);
    for (i = 0; i < 5; i++) {
      for (b_i = 0; b_i < 12; b_i++) {
        d = 0.0;
        for (Rda4_tmp = 0; Rda4_tmp < 12; Rda4_tmp++) {
          d += static_cast<double>(a[i + 5 * Rda4_tmp]) * Lambda[Rda4_tmp + 12 *
            b_i];
        }

        globalConstraintEqMatrix4[i + 5 * b_i] = d;
      }
    }

    // extract Rda from globalConstraintEqMatrix4
    for (b_i = 0; b_i < 6; b_i++) {
      for (i = 0; i < 5; i++) {
        Rda4_tmp = i + 5 * b_i;
        Rda4[Rda4_tmp] = globalConstraintEqMatrix4[Rda4_tmp];
      }
    }

    if ((slaveActiveDof5_size[0] != 0) && (slaveActiveDof5_size[1] != 0)) {
      for (i = 0; i < 5; i++) {
        Rda4[i] = globalConstraintEqMatrix4[i + 5 * (static_cast<int>
          ((slaveActiveDof5_data[0] + 6.0)) - 1)];
      }
    }

    // extract Rdd from globalConstraintEqMatrix4
    for (b_i = 0; b_i < 5; b_i++) {
      for (i = 0; i < 5; i++) {
        Rdd4[i + 5 * b_i] = globalConstraintEqMatrix4[i + 5 * (b_i + 6)];
      }

      slaveDof3[b_i] = static_cast<double>(b_i) + 1.0;
    }

    for (i = 0; i < 6; i++) {
      b_dv1[i] = static_cast<double>(i) + 1.0;
    }

    c_getNodeMaps(Rdd4, Rda4, masterNodeNum, slaveNodeNum, slaveDof3, b_dv1,
                  slaveActiveDof5_data, slaveActiveDof5_size, c_Tda, b_dDOF,
                  aDOF_data, aDOF_size);
    Tda_size[0] = 5;
    Tda_size[1] = 6;
    std::memcpy(&Tda_data[0], &c_Tda[0], 30U * sizeof(double));
    dDOF_size[0] = 5;
    for (i = 0; i < 5; i++) {
      dDOF_data[i] = b_dDOF[i];
    }
  } else {
    if (jointType == 5.0) {
      // rigid bar constraint type
      // active DOF list at joint
      // slave DOF list at joint
      // determine local active DOFs associated with slave node
      c_determineActiveDofsFromSlaveN(slaveActiveDof5_data, slaveActiveDof5_size);

      // need to define lx,ly,lz, from mesh level
      for (i = 0; i < 36; i++) {
        Rda5[i] = iv[i];
      }

      Rda5[18] = 0.0;
      Rda5[24] = -joint[6];
      Rda5[30] = joint[5];
      Rda5[19] = joint[6];
      Rda5[25] = 0.0;
      Rda5[31] = -joint[4];
      Rda5[20] = -joint[5];
      Rda5[26] = joint[4];
      Rda5[32] = 0.0;
      for (i = 0; i < 6; i++) {
        b_dv1[i] = static_cast<double>(i) + 1.0;
        dv2[i] = static_cast<double>(i) + 1.0;
      }

      d_getNodeMaps(Rda5, masterNodeNum, slaveNodeNum, b_dv1, dv2,
                    slaveActiveDof5_data, slaveActiveDof5_size, d_Tda,
                    activeDof3, aDOF_data, aDOF_size);
      Tda_size[0] = 6;
      Tda_size[1] = 6;
      std::memcpy(&Tda_data[0], &d_Tda[0], 36U * sizeof(double));
      dDOF_size[0] = 6;
      for (i = 0; i < 6; i++) {
        dDOF_data[i] = activeDof3[i];
      }
    }
  }
}

//
// This function determines the local master DOF associated with a local slave DOF.
//  Get size
// Arguments    : double slaveNodeActiveDof_data[]
//                int slaveNodeActiveDof_size[2]
// Return Type  : void
//
static void d_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2])
{
  double count;
  int i;
  boolean_T tf[3];
  double absx;
  int exponent;
  int b_exponent;
  int idx;
  int ii;
  boolean_T exitg1;
  count = 1.0;
  for (i = 0; i < 6; i++) {
    // loop over number of DOF per node
    tf[0] = false;
    tf[1] = false;
    tf[2] = false;
    absx = (static_cast<double>(i) + 1.0) / 2.0;
    frexp(absx, &exponent);
    if (i < std::ldexp(1.0, exponent - 53)) {
      tf[0] = true;
    }

    frexp(absx, &exponent);
    if (std::abs(static_cast<double>(i) - 1.0) < std::ldexp(1.0, exponent - 53))
    {
      tf[1] = true;
    }

    frexp(absx, &exponent);
    if (std::abs(static_cast<double>(i) - 2.0) < std::ldexp(1.0, exponent - 53))
    {
      tf[2] = true;
    }

    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 3)) {
      if (tf[ii]) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
      //          slaveNodeActiveDof(count) = i;
      count++;
    }
  }

  if (count > 1.0) {
    slaveNodeActiveDof_size[0] = 1;
    slaveNodeActiveDof_size[1] = 1;
    slaveNodeActiveDof_data[0] = 0.0;
    for (i = 0; i < 6; i++) {
      // loop over number of DOF per node
      tf[0] = false;
      tf[1] = false;
      tf[2] = false;
      absx = (static_cast<double>(i) + 1.0) / 2.0;
      frexp(absx, &b_exponent);
      if (i < std::ldexp(1.0, b_exponent - 53)) {
        tf[0] = true;
      }

      frexp(absx, &b_exponent);
      if (std::abs(static_cast<double>(i) - 1.0) < std::ldexp(1.0, b_exponent -
           53)) {
        tf[1] = true;
      }

      frexp(absx, &b_exponent);
      if (std::abs(static_cast<double>(i) - 2.0) < std::ldexp(1.0, b_exponent -
           53)) {
        tf[2] = true;
      }

      idx = 0;
      ii = 0;
      exitg1 = false;
      while ((!exitg1) && (ii < 3)) {
        if (tf[ii]) {
          idx = 1;
          exitg1 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
        slaveNodeActiveDof_data[0] = static_cast<double>(i) + 1.0;
      }
    }
  } else {
    slaveNodeActiveDof_size[0] = 0;
    slaveNodeActiveDof_size[1] = 0;
  }
}

//
// Arguments    : const double Rda[36]
//                double masterNodeNum
//                double slaveNodeNum
//                const double slaveDof[6]
//                const double activeDof[6]
//                const double slaveActiveDof_data[]
//                const int slaveActiveDof_size[2]
//                double Tda[36]
//                double dDOF[6]
//                double aDOF_data[]
//                int aDOF_size[1]
// Return Type  : void
//
static void d_getNodeMaps(const double Rda[36], double masterNodeNum, double
  slaveNodeNum, const double slaveDof[6], const double activeDof[6], const
  double slaveActiveDof_data[], const int slaveActiveDof_size[2], double Tda[36],
  double dDOF[6], double aDOF_data[], int aDOF_size[1])
{
  int i;
  int b_i;
  double d;
  int i1;
  double aMap2_data[1];
  double aMap[6];

  // calculate Tda
  // get number of joint DOFs for this joint
  // get number of active DOFs for this joint
  // initialize arrays
  for (i = 0; i < 6; i++) {
    for (b_i = 0; b_i < 6; b_i++) {
      d = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        d += static_cast<double>(iv[i + 6 * i1]) * Rda[i1 + 6 * b_i];
      }

      Tda[i + 6 * b_i] = d;
    }

    // get global DOF numbers of slave DOFs for this joint
    dDOF[i] = (slaveNodeNum - 1.0) * 6.0 + slaveDof[i];

    // get global DOF numbers of active DOFs for this joint from master nodes
    aMap[i] = (masterNodeNum - 1.0) * 6.0 + activeDof[i];
  }

  // determine global active DOFs associated with slave node
  if ((slaveActiveDof_size[0] == 0) || (slaveActiveDof_size[1] == 0)) {
    i = 0;
  } else {
    i = 1;
  }

  if (0 <= i - 1) {
    aMap2_data[0] = 0.0;
  }

  if ((slaveActiveDof_size[0] != 0) && (slaveActiveDof_size[1] != 0)) {
    aMap2_data[0] = (slaveNodeNum - 1.0) * 6.0 + slaveActiveDof_data[0];
  }

  if (i != 0) {
    // create overall map of active DOFs associated with this joint
    aDOF_size[0] = i + 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }

    if (0 <= i - 1) {
      aDOF_data[6] = aMap2_data[0];
    }
  } else {
    aDOF_size[0] = 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }
  }
}

//
// This function determines the local master DOF associated with a local slave DOF.
//  Get size
// Arguments    : double slaveNodeActiveDof_data[]
//                int slaveNodeActiveDof_size[2]
// Return Type  : void
//
static void e_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2])
{
  double count;
  int i;
  int idx;
  double absx;
  boolean_T tf[5];
  int exponent;
  int ii;
  static const signed char b_iv[5] = { 1, 2, 3, 5, 6 };

  boolean_T exitg1;
  int b_exponent;
  count = 1.0;
  for (i = 0; i < 6; i++) {
    // loop over number of DOF per node
    for (idx = 0; idx < 5; idx++) {
      tf[idx] = false;
    }

    absx = (static_cast<double>(i) + 1.0) / 2.0;
    for (idx = 0; idx < 5; idx++) {
      frexp(absx, &exponent);
      if (std::abs(static_cast<double>((i - b_iv[idx])) + 1.0) < std::ldexp(1.0,
           exponent - 53)) {
        tf[idx] = true;
      }
    }

    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 5)) {
      if (tf[ii]) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
      //          slaveNodeActiveDof(count) = i;
      count++;
    }
  }

  if (count > 1.0) {
    slaveNodeActiveDof_size[0] = 1;
    slaveNodeActiveDof_size[1] = 1;
    slaveNodeActiveDof_data[0] = 0.0;
    for (i = 0; i < 6; i++) {
      // loop over number of DOF per node
      for (idx = 0; idx < 5; idx++) {
        tf[idx] = false;
      }

      absx = (static_cast<double>(i) + 1.0) / 2.0;
      for (idx = 0; idx < 5; idx++) {
        frexp(absx, &b_exponent);
        if (std::abs(static_cast<double>((i - b_iv[idx])) + 1.0) < std::ldexp
            (1.0, b_exponent - 53)) {
          tf[idx] = true;
        }
      }

      idx = 0;
      ii = 0;
      exitg1 = false;
      while ((!exitg1) && (ii < 5)) {
        if (tf[ii]) {
          idx = 1;
          exitg1 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
        slaveNodeActiveDof_data[0] = static_cast<double>(i) + 1.0;
      }
    }
  } else {
    slaveNodeActiveDof_size[0] = 0;
    slaveNodeActiveDof_size[1] = 0;
  }
}

//
// This function determines the local master DOF associated with a local slave DOF.
//  Get size
// Arguments    : double slaveNodeActiveDof_data[]
//                int slaveNodeActiveDof_size[2]
// Return Type  : void
//
static void f_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2])
{
  double count;
  int i;
  int idx;
  double absx;
  boolean_T tf[5];
  int exponent;
  int ii;
  static const signed char b_iv[5] = { 1, 2, 3, 4, 6 };

  boolean_T exitg1;
  int b_exponent;
  count = 1.0;
  for (i = 0; i < 6; i++) {
    // loop over number of DOF per node
    for (idx = 0; idx < 5; idx++) {
      tf[idx] = false;
    }

    absx = (static_cast<double>(i) + 1.0) / 2.0;
    for (idx = 0; idx < 5; idx++) {
      frexp(absx, &exponent);
      if (std::abs(static_cast<double>((i - b_iv[idx])) + 1.0) < std::ldexp(1.0,
           exponent - 53)) {
        tf[idx] = true;
      }
    }

    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 5)) {
      if (tf[ii]) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
      //          slaveNodeActiveDof(count) = i;
      count++;
    }
  }

  if (count > 1.0) {
    slaveNodeActiveDof_size[0] = 1;
    slaveNodeActiveDof_size[1] = 1;
    slaveNodeActiveDof_data[0] = 0.0;
    for (i = 0; i < 6; i++) {
      // loop over number of DOF per node
      for (idx = 0; idx < 5; idx++) {
        tf[idx] = false;
      }

      absx = (static_cast<double>(i) + 1.0) / 2.0;
      for (idx = 0; idx < 5; idx++) {
        frexp(absx, &b_exponent);
        if (std::abs(static_cast<double>((i - b_iv[idx])) + 1.0) < std::ldexp
            (1.0, b_exponent - 53)) {
          tf[idx] = true;
        }
      }

      idx = 0;
      ii = 0;
      exitg1 = false;
      while ((!exitg1) && (ii < 5)) {
        if (tf[ii]) {
          idx = 1;
          exitg1 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
        slaveNodeActiveDof_data[0] = static_cast<double>(i) + 1.0;
      }
    }
  } else {
    slaveNodeActiveDof_size[0] = 0;
    slaveNodeActiveDof_size[1] = 0;
  }
}

//
// This function determines the local master DOF associated with a local slave DOF.
//  Get size
// Arguments    : double slaveNodeActiveDof_data[]
//                int slaveNodeActiveDof_size[2]
// Return Type  : void
//
static void g_determineActiveDofsFromSlaveN(double slaveNodeActiveDof_data[],
  int slaveNodeActiveDof_size[2])
{
  double count;
  int i;
  int idx;
  double absx;
  boolean_T tf[5];
  int exponent;
  int ii;
  boolean_T exitg1;
  int b_exponent;
  count = 1.0;
  for (i = 0; i < 6; i++) {
    // loop over number of DOF per node
    for (idx = 0; idx < 5; idx++) {
      tf[idx] = false;
    }

    absx = (static_cast<double>(i) + 1.0) / 2.0;
    for (idx = 0; idx < 5; idx++) {
      frexp(absx, &exponent);
      if (std::abs(static_cast<double>((i - idx))) < std::ldexp(1.0, exponent -
           53)) {
        tf[idx] = true;
      }
    }

    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 5)) {
      if (tf[ii]) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
      //          slaveNodeActiveDof(count) = i;
      count++;
    }
  }

  if (count > 1.0) {
    slaveNodeActiveDof_size[0] = 1;
    slaveNodeActiveDof_size[1] = 1;
    slaveNodeActiveDof_data[0] = 0.0;
    for (i = 0; i < 6; i++) {
      // loop over number of DOF per node
      for (idx = 0; idx < 5; idx++) {
        tf[idx] = false;
      }

      absx = (static_cast<double>(i) + 1.0) / 2.0;
      for (idx = 0; idx < 5; idx++) {
        frexp(absx, &b_exponent);
        if (std::abs(static_cast<double>((i - idx))) < std::ldexp(1.0,
             b_exponent - 53)) {
          tf[idx] = true;
        }
      }

      idx = 0;
      ii = 0;
      exitg1 = false;
      while ((!exitg1) && (ii < 5)) {
        if (tf[ii]) {
          idx = 1;
          exitg1 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        // if i is not in slaveDof list add it to a list of local active DOFs associated with a slave node 
        slaveNodeActiveDof_data[0] = static_cast<double>(i) + 1.0;
      }
    }
  } else {
    slaveNodeActiveDof_size[0] = 0;
    slaveNodeActiveDof_size[1] = 0;
  }
}

//
// Arguments    : double masterNodeNum
//                double slaveNodeNum
//                const double slaveDof[6]
//                const double activeDof[6]
//                const double slaveActiveDof_data[]
//                const int slaveActiveDof_size[2]
//                double dDOF[6]
//                double aDOF_data[]
//                int aDOF_size[1]
// Return Type  : void
//
static void getNodeMaps(double masterNodeNum, double slaveNodeNum, const double
  slaveDof[6], const double activeDof[6], const double slaveActiveDof_data[],
  const int slaveActiveDof_size[2], double dDOF[6], double aDOF_data[], int
  aDOF_size[1])
{
  int i;
  double aMap[6];
  double aMap2_data[1];
  int b_i;

  // calculate Tda
  // get number of joint DOFs for this joint
  // get number of active DOFs for this joint
  // initialize arrays
  for (i = 0; i < 6; i++) {
    // get global DOF numbers of slave DOFs for this joint
    dDOF[i] = (slaveNodeNum - 1.0) * 6.0 + slaveDof[i];

    // get global DOF numbers of active DOFs for this joint from master nodes
    aMap[i] = (masterNodeNum - 1.0) * 6.0 + activeDof[i];
  }

  // determine global active DOFs associated with slave node
  if ((slaveActiveDof_size[0] == 0) || (slaveActiveDof_size[1] == 0)) {
    i = 0;
  } else {
    i = 1;
  }

  if (0 <= i - 1) {
    aMap2_data[0] = 0.0;
  }

  if ((slaveActiveDof_size[0] != 0) && (slaveActiveDof_size[1] != 0)) {
    aMap2_data[0] = (slaveNodeNum - 1.0) * 6.0 + slaveActiveDof_data[0];
  }

  if (i != 0) {
    // create overall map of active DOFs associated with this joint
    aDOF_size[0] = i + 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }

    if (0 <= i - 1) {
      aDOF_data[6] = aMap2_data[0];
    }
  } else {
    aDOF_size[0] = 6;
    for (b_i = 0; b_i < 6; b_i++) {
      aDOF_data[b_i] = aMap[b_i];
    }
  }
}

//
// createJointTransform   Creates transformation matrix for joint constaints
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [jointTransform,reducedDOF] = createJointTransform(joint,numNodes,
//                                                       numDofPerNode)
//
//    This function calculates the eigenvalues and vectors of a structural
//    dynamic system.
//
//    input:
//    joint         = object containing joint data
//    numModes      = number of nodes in mesh
//    numDofPerNode = number of degrees of freedom per node
//
//    output:
//    jointTransform = joint transformation matrix
//    reducedDOF     = map of original DOF numbering to reduced DOF numbering
// Arguments    : const emxArray_real_T *joint
//                double numNodes
//                emxArray_real_T *output_jointTransform
//                emxArray_real_T *output_reducedDOF
// Return Type  : void
//
void createJointTransform(const emxArray_real_T *joint, double numNodes,
  emxArray_real_T *output_jointTransform, emxArray_real_T *output_reducedDOF)
{
  double adNumDof;
  double dependentCount;
  double count;
  int i;
  int b_i;
  emxArray_real_T *slaveDof;
  double absx;
  int nz;
  int j;
  int i1;
  int i2;
  signed char con[5];
  static const signed char b_iv[5] = { 1, 2, 3, 4, 6 };

  static const signed char b_iv1[5] = { 1, 2, 3, 5, 6 };

  unsigned int unnamed_idx_1;
  emxArray_boolean_T *tf;
  emxArray_uint32_T *reducedDOF;
  int na;
  int exponent;
  int k;
  emxArray_int32_T *b_r;
  emxArray_int32_T *ii;
  int b_exponent;
  double b_joint[8];
  double Tda_data[54];
  int Tda_size[2];
  double dDOF_data[6];
  int dDOF_size[1];
  double aDOF_data[7];
  int aDOF_size[1];
  int i3;
  int c_exponent;
  boolean_T exitg1;

  // get number of joints in model
  // extract number of active DOFs, number of dependent DOFs, slave DOF numbers
  // This function gets the total number of DOFs in the model, active
  // number of DOFs in the model, and a list of slave DOFs that will be
  // eliminated by joint constraints.
  adNumDof = numNodes * 6.0;

  // total number of DOFs (active and dependent)
  // get number of joints
  // Get Count
  dependentCount = 0.0;
  count = 1.0;
  i = joint->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    // loop over number of joints
    absx = joint->data[b_i + joint->size[0] * 3];
    if ((absx == 0.0) || (absx == 5.0)) {
      // for a "fixed/weld" joint type or rigid bar constraint
      for (j = 0; j < 6; j++) {
        count++;
      }
    }

    if (absx == 1.0) {
      // for a "pinned" joint type
      count++;
      count++;
      count++;
    }

    if (absx == 2.0) {
      // for a single axis hinge joint along a local "2" axis of a joint
      for (j = 0; j < 5; j++) {
        count++;
      }
    }

    if (absx == 3.0) {
      // for a single axis hinge =joint along a local "1" axis of a joint
      for (j = 0; j < 5; j++) {
        count++;
      }
    }

    if (absx == 4.0) {
      // for a single axis hinge joint along a local "3" axis of a joint
      for (j = 0; j < 5; j++) {
        count++;
      }
    }
  }

  emxInit_real_T(&slaveDof, 2);
  i = slaveDof->size[0] * slaveDof->size[1];
  slaveDof->size[0] = 1;
  nz = static_cast<int>(count);
  slaveDof->size[1] = nz;
  emxEnsureCapacity_real_T(slaveDof, i);
  for (i = 0; i < nz; i++) {
    slaveDof->data[i] = 0.0;
  }

  count = 1.0;
  i = joint->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    // loop over number of joints
    absx = joint->data[b_i + joint->size[0] * 3];
    if ((absx == 0.0) || (absx == 5.0)) {
      // for a "fixed/weld" joint type or rigid bar constraint
      // all DOFs of a slave node are constrained
      dependentCount += 6.0;

      // increment number of dependent DOFs
      for (j = 0; j < 6; j++) {
        slaveDof->data[static_cast<int>(count) - 1] = 6.0 * (joint->data[b_i +
          joint->size[0] * 2] - 1.0) + (static_cast<double>(j) + 1.0);

        // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
        count++;
      }
    }

    if (joint->data[b_i + joint->size[0] * 3] == 1.0) {
      // for a "pinned" joint type
      // only translational (first 3) DOFs of a slave node are  constrained
      dependentCount += 3.0;

      // increment number of dependent DOFs
      absx = 6.0 * (joint->data[b_i + joint->size[0] * 2] - 1.0);
      slaveDof->data[static_cast<int>(count) - 1] = absx + 1.0;

      // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
      count++;
      slaveDof->data[static_cast<int>(count) - 1] = absx + 2.0;

      // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
      count++;
      slaveDof->data[static_cast<int>(count) - 1] = absx + 3.0;

      // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
      count++;
    }

    if (joint->data[b_i + joint->size[0] * 3] == 2.0) {
      // for a single axis hinge joint along a local "2" axis of a joint
      absx = std::abs(joint->data[b_i + joint->size[0] * 6]);
      if ((std::abs(absx - 90.0) < 0.001) || (std::abs(absx - 270.0) < 0.001)) {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = b_iv1[i1];
        }
      } else {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = b_iv[i1];
        }

        // all but 5th DOF of a  slave node are constrained
      }

      dependentCount += 5.0;

      // increment number of dependent DOFs
      for (j = 0; j < 5; j++) {
        slaveDof->data[static_cast<int>(count) - 1] = 6.0 * (joint->data[b_i +
          joint->size[0] * 2] - 1.0) + static_cast<double>(con[j]);

        // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
        count++;
      }
    }

    if (joint->data[b_i + joint->size[0] * 3] == 3.0) {
      // for a single axis hinge =joint along a local "1" axis of a joint
      absx = std::abs(joint->data[b_i + joint->size[0] * 7]);
      if ((std::abs(absx - 90.0) < 0.001) || (std::abs(absx - 270.0) < 0.001)) {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = static_cast<signed char>((i1 + 1));
        }
      } else if ((std::abs(std::abs(joint->data[b_i + joint->size[0] * 6]) -
                           90.0) < 0.001) || (std::abs(std::abs(joint->data[b_i
                    + joint->size[0] * 6]) - 270.0) < 0.001)) {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = b_iv[i1];
        }
      } else {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = b_iv1[i1];
        }

        // all but the 4th DOF of a slave node are constrained
      }

      dependentCount += 5.0;

      // increment number of dependent DOFs
      for (j = 0; j < 5; j++) {
        slaveDof->data[static_cast<int>(count) - 1] = 6.0 * (joint->data[b_i +
          joint->size[0] * 2] - 1.0) + static_cast<double>(con[j]);

        // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
        count++;
      }
    }

    if (joint->data[b_i + joint->size[0] * 3] == 4.0) {
      // for a single axis hinge joint along a local "3" axis of a joint
      absx = std::abs(joint->data[b_i + joint->size[0] * 7]);
      if ((std::abs(absx - 90.0) < 0.001) || (std::abs(absx - 270.0) < 0.001)) {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = b_iv1[i1];
        }

        if ((std::abs(std::abs(joint->data[b_i + joint->size[0] * 6]) - 90.0) <
             0.001) || (std::abs(std::abs(joint->data[b_i + joint->size[0] * 6])
              - 270.0) < 0.001)) {
          for (i1 = 0; i1 < 5; i1++) {
            con[i1] = b_iv[i1];
          }
        }
      } else {
        for (i1 = 0; i1 < 5; i1++) {
          con[i1] = static_cast<signed char>((i1 + 1));
        }
      }

      dependentCount += 5.0;

      // all but the 6th DOF of a slave node are constrained
      for (j = 0; j < 5; j++) {
        slaveDof->data[static_cast<int>(count) - 1] = 6.0 * (joint->data[b_i +
          joint->size[0] * 2] - 1.0) + static_cast<double>(con[j]);

        // assign slave DOFs (joint(i,3) is the slave node number associated with this joint 
        count++;
      }
    }
  }

  // calculate number of active DOFs in the model
  // initialize joint transformation matrix
  i = static_cast<int>(adNumDof);
  i1 = output_jointTransform->size[0] * output_jointTransform->size[1];
  output_jointTransform->size[0] = i;
  i2 = static_cast<int>((adNumDof - dependentCount));
  output_jointTransform->size[1] = i2;
  emxEnsureCapacity_real_T(output_jointTransform, i1);
  nz = i * i2;
  for (i = 0; i < nz; i++) {
    output_jointTransform->data[i] = 0.0;
  }

  // form reduced DOF vector which maps original DOF numbering to reduced DOF
  // numbering
  // Get Count
  count = 1.0;
  i = static_cast<int>((numNodes * 6.0));
  if (0 <= i - 1) {
    unnamed_idx_1 = static_cast<unsigned int>(slaveDof->size[1]);
  }

  emxInit_boolean_T(&tf, 2);
  for (b_i = 0; b_i < i; b_i++) {
    // loop over total number of DOFs in model
    na = slaveDof->size[1];
    i1 = tf->size[0] * tf->size[1];
    tf->size[0] = 1;
    tf->size[1] = static_cast<int>(unnamed_idx_1);
    emxEnsureCapacity_boolean_T(tf, i1);
    nz = static_cast<int>(unnamed_idx_1);
    for (i1 = 0; i1 < nz; i1++) {
      tf->data[i1] = false;
    }

    for (j = 0; j < na; j++) {
      frexp((static_cast<double>(b_i) + 1.0) / 2.0, &exponent);
      if (std::abs((static_cast<double>(b_i) + 1.0) - slaveDof->data[j]) < std::
          ldexp(1.0, exponent - 53)) {
        tf->data[j] = true;
      }
    }

    j = tf->size[1];
    nz = tf->data[0];
    for (k = 2; k <= j; k++) {
      nz += tf->data[k - 1];
    }

    if (nz == 0) {
      //          reducedDOF(count) = i; %if DOF is NOT a slave DOF include it in reducedDOF 
      count++;
    }
  }

  emxInit_uint32_T(&reducedDOF, 2);
  i1 = reducedDOF->size[0] * reducedDOF->size[1];
  reducedDOF->size[0] = 1;
  nz = static_cast<int>((count - 1.0));
  reducedDOF->size[1] = nz;
  emxEnsureCapacity_uint32_T(reducedDOF, i1);
  for (i1 = 0; i1 < nz; i1++) {
    reducedDOF->data[i1] = 0U;
  }

  count = 1.0;
  if (0 <= i - 1) {
    unnamed_idx_1 = static_cast<unsigned int>(slaveDof->size[1]);
  }

  for (b_i = 0; b_i < i; b_i++) {
    // loop over total number of DOFs in model
    na = slaveDof->size[1];
    i1 = tf->size[0] * tf->size[1];
    tf->size[0] = 1;
    tf->size[1] = static_cast<int>(unnamed_idx_1);
    emxEnsureCapacity_boolean_T(tf, i1);
    nz = static_cast<int>(unnamed_idx_1);
    for (i1 = 0; i1 < nz; i1++) {
      tf->data[i1] = false;
    }

    for (j = 0; j < na; j++) {
      frexp((static_cast<double>(b_i) + 1.0) / 2.0, &b_exponent);
      if (std::abs((static_cast<double>(b_i) + 1.0) - slaveDof->data[j]) < std::
          ldexp(1.0, b_exponent - 53)) {
        tf->data[j] = true;
      }
    }

    j = tf->size[1];
    nz = tf->data[0];
    for (k = 2; k <= j; k++) {
      nz += tf->data[k - 1];
    }

    if (nz == 0) {
      reducedDOF->data[static_cast<int>(count) - 1] = b_i + 1U;

      // if DOF is NOT a slave DOF include it in reducedDOF
      count++;
    }
  }

  emxFree_real_T(&slaveDof);

  // create identity portion of transformation matrix (This is done by Craig,
  // but here the original DOF ordering is retained
  for (b_i = 0; b_i < i2; b_i++) {
    // loop over number of active DOFs
    output_jointTransform->data[(static_cast<int>(reducedDOF->data[b_i]) +
      output_jointTransform->size[0] * b_i) - 1] = 1.0;

    // mapping of active DOFs in full DOF list to reduced DOF list
  }

  // impose Tda portion of identity matrix and map to appropriate locations
  i = joint->size[0];
  emxInit_int32_T(&b_r, 1);
  emxInit_int32_T(&ii, 2);
  for (b_i = 0; b_i < i; b_i++) {
    //  loop of number of joints in the model
    // get joint type
    // get slave node number associated with joint
    // get master node number associated with joint
    // get psi orientation angle associated with joint
    // get theta orientation angle associated with joint
    // Tda is a local transform between dependent and active DOFs for nodes
    // associated with a particular joint, dDOF is a listing of dependent
    // global DOFs associated with this joint, aDOF is a listing of
    // active global DOFs associated with this joint.
    for (i1 = 0; i1 < 8; i1++) {
      b_joint[i1] = joint->data[b_i + joint->size[0] * i1];
    }

    createTda(joint->data[b_i + joint->size[0] * 3], joint->data[b_i +
              joint->size[0] * 2], joint->data[b_i + joint->size[0]],
              joint->data[b_i + joint->size[0] * 6], joint->data[b_i +
              joint->size[0] * 7], b_joint, Tda_data, Tda_size, dDOF_data,
              dDOF_size, aDOF_data, aDOF_size);
    i1 = aDOF_size[0];
    if (0 <= aDOF_size[0] - 1) {
      i3 = dDOF_size[0];
      if (0 <= dDOF_size[0] - 1) {
        unnamed_idx_1 = static_cast<unsigned int>(reducedDOF->size[1]);
      }
    }

    for (exponent = 0; exponent < i1; exponent++) {
      // loop over global active DOFs associated with joint
      for (k = 0; k < i3; k++) {
        // loop over global dependent DOFs associated with joint
        adNumDof = aDOF_data[exponent];
        na = reducedDOF->size[1];
        i2 = tf->size[0] * tf->size[1];
        tf->size[0] = 1;
        tf->size[1] = static_cast<int>(unnamed_idx_1);
        emxEnsureCapacity_boolean_T(tf, i2);
        nz = static_cast<int>(unnamed_idx_1);
        for (i2 = 0; i2 < nz; i2++) {
          tf->data[i2] = false;
        }

        for (j = 0; j < na; j++) {
          absx = std::abs(adNumDof / 2.0);
          if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
            if (absx <= 2.2250738585072014E-308) {
              absx = 4.94065645841247E-324;
            } else {
              frexp(absx, &c_exponent);
              absx = std::ldexp(1.0, c_exponent - 53);
            }
          } else {
            absx = rtNaN;
          }

          if (std::abs(adNumDof - static_cast<double>(reducedDOF->data[j])) <
              absx) {
            tf->data[j] = true;
          }
        }

        i2 = tf->size[1];
        nz = 0;
        j = ii->size[0] * ii->size[1];
        ii->size[0] = 1;
        ii->size[1] = tf->size[1];
        emxEnsureCapacity_int32_T(ii, j);
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= i2 - 1)) {
          if (tf->data[j]) {
            nz++;
            ii->data[nz - 1] = j + 1;
            if (nz >= i2) {
              exitg1 = true;
            } else {
              j++;
            }
          } else {
            j++;
          }
        }

        if (tf->size[1] == 1) {
          if (nz == 0) {
            ii->size[0] = 1;
            ii->size[1] = 0;
          }
        } else {
          i2 = ii->size[0] * ii->size[1];
          if (1 > nz) {
            ii->size[1] = 0;
          } else {
            ii->size[1] = nz;
          }

          emxEnsureCapacity_int32_T(ii, i2);
        }

        // determine reduced DOF associated with active DOF from original DOF listing 
        i2 = b_r->size[0];
        b_r->size[0] = ii->size[1];
        emxEnsureCapacity_int32_T(b_r, i2);
        nz = ii->size[1];
        for (i2 = 0; i2 < nz; i2++) {
          b_r->data[i2] = ii->data[i2];
        }

        nz = b_r->size[0];
        for (i2 = 0; i2 < nz; i2++) {
          output_jointTransform->data[(static_cast<int>(dDOF_data[k]) +
            output_jointTransform->size[0] * (b_r->data[i2] - 1)) - 1] =
            Tda_data[k + Tda_size[0] * exponent];
        }

        // map local joint transformation matrix (Tda) to entries in global transformation matrix (jointTransform) 
      }
    }
  }

  emxFree_int32_T(&ii);
  emxFree_boolean_T(&tf);
  emxFree_int32_T(&b_r);
  i = output_reducedDOF->size[0] * output_reducedDOF->size[1];
  output_reducedDOF->size[0] = 1;
  output_reducedDOF->size[1] = reducedDOF->size[1];
  emxEnsureCapacity_real_T(output_reducedDOF, i);
  nz = reducedDOF->size[0] * reducedDOF->size[1];
  for (i = 0; i < nz; i++) {
    output_reducedDOF->data[i] = reducedDOF->data[i];
  }

  emxFree_uint32_T(&reducedDOF);
}

//
// File trailer for createJointTransform.cpp
//
// [EOF]
//
