//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateLambda.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 15:21:39
//

// Include Files
#include "calculateLambda.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions

//
// calculateLambda Calculates transformation matrix from element to hub frame
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [lambda] = calculateLambda(theta1,theta2,theta3 )
//
//    This function calculates a transformation matrix to transform the
//    element degree of freedom vector (12 DOFs) from the hub frame to
//    the element frame. The transformation matrix is constructed via the
//    direction cosine matrices of a 3-2-1 Euler rotation sequence.
//
//       input:
//       theta1        = angle (rad) of rotation for 1st rotation
//                       of 3-2-1 sequence
//       theta2        = angle (rad) of rotation for 2nd rotation
//                       of 3-2-1 sequence
//       theta3        = angle (rad) of rotation for 3rd rotation
//                       of 3-2-1 sequence
// Arguments    : double theta1
//                double theta2
//                double theta3
//                double lambda[144]
// Return Type  : void
//
void calculateLambda(double theta1, double theta2, double theta3, double lambda
                     [144])
{
  double ct1;
  double st1;
  double ct2;
  double st2;
  double ct3;
  double st3;
  double fac1;
  double fac2;
  double dcm[9];
  int i;
  int lambda_tmp;
  int b_lambda_tmp;
  int c_lambda_tmp;

  //       output:
  //       lambda        = 12 x 12 transformation matrix
  //  dcm that is created is [dcm] = [M1(theta3)][M2(theta2)][M3(theta1)]
  // calculateLambda Calculates transformation matrix from element to hub frame
  //  **********************************************************************
  //  *                   Part of the SNL OWENS Toolkit                    *
  //  * Developed by Sandia National Laboratories Wind Energy Technologies *
  //  *             See license.txt for disclaimer information             *
  //  **********************************************************************
  //    [lambda] = calculateLambda(theta1,theta2,theta3 )
  //
  //    This function calculates a transformation matrix to transform the
  //    element degree of freedom vector (12 DOFs) from the element frame to
  //    the hub frame. The transformation matrix is constructed via the
  //    direction cosine matrices of a 3-2-1 Euler rotation sequence.
  //
  //       input:
  //       theta1        = angle (rad) of rotation for 1st rotation
  //                       of 3-2-1 sequence
  //       theta2        = angle (rad) of rotation for 2nd rotation
  //                       of 3-2-1 sequence
  //       theta3        = angle (rad) of rotation for 3rd rotation
  //                       of 3-2-1 sequence
  //       output:
  //       lambda        = 3 x 3 transformation matrix
  //  dcm that is created is [dcm] = [M1(theta3)][M2(theta2)][M3(theta1)]
  ct1 = std::cos(theta1);
  st1 = std::sin(theta1);
  ct2 = std::cos(theta2);
  st2 = std::sin(theta2);
  ct3 = std::cos(theta3);
  st3 = std::sin(theta3);
  fac1 = st3 * st2;
  fac2 = ct3 * st2;
  dcm[0] = ct2 * ct1;
  dcm[3] = ct2 * st1;
  dcm[6] = -st2;
  dcm[1] = fac1 * ct1 - ct3 * st1;
  dcm[4] = fac1 * st1 + ct3 * ct1;
  dcm[7] = st3 * ct2;
  dcm[2] = fac2 * ct1 + st3 * st1;
  dcm[5] = fac2 * st1 - st3 * ct1;
  dcm[8] = ct3 * ct2;
  std::memset(&lambda[0], 0, 144U * sizeof(double));
  for (i = 0; i < 3; i++) {
    ct1 = dcm[3 * i];
    lambda[12 * i] = ct1;
    lambda_tmp = 12 * (i + 3);
    lambda[lambda_tmp + 3] = ct1;
    b_lambda_tmp = 12 * (i + 6);
    lambda[b_lambda_tmp + 6] = ct1;
    c_lambda_tmp = 12 * (i + 9);
    lambda[c_lambda_tmp + 9] = ct1;
    ct1 = dcm[3 * i + 1];
    lambda[12 * i + 1] = ct1;
    lambda[lambda_tmp + 4] = ct1;
    lambda[b_lambda_tmp + 7] = ct1;
    lambda[c_lambda_tmp + 10] = ct1;
    ct1 = dcm[3 * i + 2];
    lambda[12 * i + 2] = ct1;
    lambda[lambda_tmp + 5] = ct1;
    lambda[b_lambda_tmp + 8] = ct1;
    lambda[c_lambda_tmp + 11] = ct1;
  }
}

//
// File trailer for calculateLambda.cpp
//
// [EOF]
//
