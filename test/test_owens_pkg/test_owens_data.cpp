//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: test_owens_data.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 03-Apr-2020 15:56:19
//

// Include Files
#include "test_owens_data.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Variable Definitions
const boolean_T bv[128] = { false, false, false, false, false, false, false,
  false, false, true, true, true, true, true, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, true, true,
  true, true, true, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false };

const double dv[4] = { 0.33998104358485631, -0.33998104358485631,
  0.86113631159405257, -0.86113631159405257 };

const double dv1[4] = { 0.65214515486254621, 0.65214515486254621,
  0.34785484513745385, 0.34785484513745385 };

const signed char iv1[12] = { 1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12 };

const e_struct_T r1 = { { 0.0, 0.0, 0.0, 0.0 },// K11
  { 0.0, 0.0, 0.0, 0.0 },              // K12
  { 0.0, 0.0, 0.0, 0.0 },              // K13
  { 0.0, 0.0, 0.0, 0.0 },              // K14
  { 0.0, 0.0, 0.0, 0.0 },              // K15
  { 0.0, 0.0, 0.0, 0.0 },              // K16
  { 0.0, 0.0, 0.0, 0.0 },              // K22
  { 0.0, 0.0, 0.0, 0.0 },              // K23
  { 0.0, 0.0, 0.0, 0.0 },              // K24
  { 0.0, 0.0, 0.0, 0.0 },              // K25
  { 0.0, 0.0, 0.0, 0.0 },              // K26
  { 0.0, 0.0, 0.0, 0.0 },              // K33
  { 0.0, 0.0, 0.0, 0.0 },              // K34
  { 0.0, 0.0, 0.0, 0.0 },              // K35
  { 0.0, 0.0, 0.0, 0.0 },              // K36
  { 0.0, 0.0, 0.0, 0.0 },              // K44
  { 0.0, 0.0, 0.0, 0.0 },              // K45
  { 0.0, 0.0, 0.0, 0.0 },              // K46
  { 0.0, 0.0, 0.0, 0.0 },              // K55
  { 0.0, 0.0, 0.0, 0.0 },              // K56
  { 0.0, 0.0, 0.0, 0.0 },              // K66
  { 0.0, 0.0, 0.0, 0.0 },              // M11
  { 0.0, 0.0, 0.0, 0.0 },              // M15
  { 0.0, 0.0, 0.0, 0.0 },              // M16
  { 0.0, 0.0, 0.0, 0.0 },              // M22
  { 0.0, 0.0, 0.0, 0.0 },              // M24
  { 0.0, 0.0, 0.0, 0.0 },              // M33
  { 0.0, 0.0, 0.0, 0.0 },              // M34
  { 0.0, 0.0, 0.0, 0.0 },              // M44
  { 0.0, 0.0, 0.0, 0.0 },              // M55
  { 0.0, 0.0, 0.0, 0.0 },              // M56
  { 0.0, 0.0, 0.0, 0.0 },              // M66
  { 0.0, 0.0, 0.0, 0.0 },              // S11
  { 0.0, 0.0, 0.0, 0.0 },              // S12
  { 0.0, 0.0, 0.0, 0.0 },              // S13
  { 0.0, 0.0, 0.0, 0.0 },              // S15
  { 0.0, 0.0, 0.0, 0.0 },              // S16
  { 0.0, 0.0, 0.0, 0.0 },              // S22
  { 0.0, 0.0, 0.0, 0.0 },              // S23
  { 0.0, 0.0, 0.0, 0.0 },              // S25
  { 0.0, 0.0, 0.0, 0.0 },              // S26
  { 0.0, 0.0, 0.0, 0.0 },              // S33
  { 0.0, 0.0, 0.0, 0.0 },              // S35
  { 0.0, 0.0, 0.0, 0.0 },              // S36
  { 0.0, 0.0, 0.0, 0.0 },              // S55
  { 0.0, 0.0, 0.0, 0.0 },              // S56
  { 0.0, 0.0, 0.0, 0.0 },              // S66
  { 0.0, 0.0, 0.0, 0.0 },              // S14_1
  { 0.0, 0.0, 0.0, 0.0 },              // S14_2
  { 0.0, 0.0, 0.0, 0.0 },              // S24_1
  { 0.0, 0.0, 0.0, 0.0 },              // S24_2
  { 0.0, 0.0, 0.0, 0.0 },              // S34_1
  { 0.0, 0.0, 0.0, 0.0 },              // S34_2
  { 0.0, 0.0, 0.0, 0.0 },              // S45_1
  { 0.0, 0.0, 0.0, 0.0 },              // S45_2
  { 0.0, 0.0, 0.0, 0.0 },              // S46_1
  { 0.0, 0.0, 0.0, 0.0 },              // S46_2
  { 0.0, 0.0, 0.0, 0.0 },              // S44_1
  { 0.0, 0.0, 0.0, 0.0 },              // S44_2
  { 0.0, 0.0, 0.0, 0.0 },              // S44_3
  { 0.0, 0.0, 0.0, 0.0 },              // C12
  { 0.0, 0.0, 0.0, 0.0 },              // C13
  { 0.0, 0.0, 0.0, 0.0 },              // C23
  { 0.0, 0.0, 0.0, 0.0 },              // C24
  { 0.0, 0.0, 0.0, 0.0 },              // C25
  { 0.0, 0.0, 0.0, 0.0 },              // C26
  { 0.0, 0.0, 0.0, 0.0 },              // C34
  { 0.0, 0.0, 0.0, 0.0 },              // C35
  { 0.0, 0.0, 0.0, 0.0 },              // C36
  { 0.0, 0.0, 0.0, 0.0 },              // C14_1
  { 0.0, 0.0, 0.0, 0.0 },              // C14_2
  { 0.0, 0.0, 0.0, 0.0 },              // C45_1
  { 0.0, 0.0, 0.0, 0.0 },              // C45_2
  { 0.0, 0.0, 0.0, 0.0 },              // C46_1
  { 0.0, 0.0, 0.0, 0.0 },              // C46_2
  0.0,                                 // mel
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },// moiel
  { 0.0, 0.0, 0.0 }                    // xmel
};

boolean_T isInitialized_test_owens = false;

//
// File trailer for test_owens_data.cpp
//
// [EOF]
//
