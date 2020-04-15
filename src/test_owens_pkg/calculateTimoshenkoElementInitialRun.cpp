//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: calculateTimoshenkoElementInitialRun.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "calculateTimoshenkoElementInitialRun.h"
#include "calculateLambda.h"
#include "calculateShapeFunctions.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include <cstring>
#include <string.h>

// Function Declarations
static void calculateElementMass(double rhoA, double rhoIyy, double rhoIzz,
  double rhoIyz, double rhoJ, double ycm, double zcm, double x, double y, double
  z, double integrationFactor, double *M, double Itens[9], double xm[3]);

// Function Definitions

//
// This function calculates element mass properties.
// Arguments    : double rhoA
//                double rhoIyy
//                double rhoIzz
//                double rhoIyz
//                double rhoJ
//                double ycm
//                double zcm
//                double x
//                double y
//                double z
//                double integrationFactor
//                double *M
//                double Itens[9]
//                double xm[3]
// Return Type  : void
//
static void calculateElementMass(double rhoA, double rhoIyy, double rhoIzz,
  double rhoIyz, double rhoJ, double ycm, double zcm, double x, double y, double
  z, double integrationFactor, double *M, double Itens[9], double xm[3])
{
  double b_rhoA;
  double rhoA_tmp;
  double b_rhoA_tmp;
  double c_rhoA[9];
  double c_rhoA_tmp;
  double d_rhoA_tmp;
  double b_integrationFactor[9];
  int i;

  //  function [F] = calculateVec1(f,integrationFactor,N,F)
  //  %This function is a general routine to calculate an element vector
  //      len=length(N);
  //      for i=1:len
  //          F(i) = F(i) + f*N(i)*integrationFactor;
  //      end
  //
  //  end
  b_rhoA = rhoA * integrationFactor;
  *M += b_rhoA;
  y += ycm;
  z += zcm;

  // Total MOI = parallel axis theorem + local MOI
  rhoA_tmp = z * z;
  b_rhoA_tmp = y * y;
  c_rhoA[0] = b_rhoA * (b_rhoA_tmp + rhoA_tmp);
  c_rhoA_tmp = b_rhoA * (-x * y);
  c_rhoA[3] = c_rhoA_tmp;
  d_rhoA_tmp = b_rhoA * (-x * z);
  c_rhoA[6] = d_rhoA_tmp;
  c_rhoA[1] = c_rhoA_tmp;
  c_rhoA_tmp = x * x;
  c_rhoA[4] = b_rhoA * (c_rhoA_tmp + rhoA_tmp);
  rhoA_tmp = b_rhoA * (-y * z);
  c_rhoA[7] = rhoA_tmp;
  c_rhoA[2] = d_rhoA_tmp;
  c_rhoA[5] = rhoA_tmp;
  c_rhoA[8] = b_rhoA * (c_rhoA_tmp + b_rhoA_tmp);
  b_integrationFactor[0] = integrationFactor * rhoJ;
  b_integrationFactor[3] = integrationFactor * 0.0;
  b_integrationFactor[6] = integrationFactor * 0.0;
  b_integrationFactor[1] = integrationFactor * 0.0;
  b_integrationFactor[4] = integrationFactor * rhoIyy;
  b_rhoA = integrationFactor * rhoIyz;
  b_integrationFactor[7] = b_rhoA;
  b_integrationFactor[2] = integrationFactor * 0.0;
  b_integrationFactor[5] = b_rhoA;
  b_integrationFactor[8] = integrationFactor * rhoIzz;
  for (i = 0; i < 9; i++) {
    Itens[i] = (Itens[i] + c_rhoA[i]) + b_integrationFactor[i];
  }

  xm[0] += x * rhoA * integrationFactor;
  xm[1] += y * rhoA * integrationFactor;
  xm[2] += z * rhoA * integrationFactor;
}

//
// calculateTimoshenkoElementInitialRun performs initial element calculations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [elStorage] = calculateTimoshenkoElementInitialRun(input)
//
//    This function performs initial element calculations and stores them for
//    later use and efficiency gains.
//
//       input:
//       input      = object containing element input
//
//       output:
//       elStorage  = object containing stored element data
//
// -------- assign input block ----------------
// Arguments    : const double input_xloc[2]
//                const double input_sectionProps_twist[2]
//                const double input_sectionProps_rhoA[2]
//                const double input_sectionProps_EIyy[2]
//                const double input_sectionProps_EIzz[2]
//                const double input_sectionProps_GJ[2]
//                const double input_sectionProps_EA[2]
//                const double input_sectionProps_rhoIyy[2]
//                const double input_sectionProps_rhoIzz[2]
//                const double input_sectionProps_rhoJ[2]
//                const double input_sectionProps_zcm[2]
//                const double input_sectionProps_ycm[2]
//                double input_sweepAngle
//                double input_coneAngle
//                double input_rollAngle
//                const double input_x[2]
//                const double input_y[2]
//                const double input_z[2]
//                boolean_T input_concMassFlag
//                const double input_concMass[8]
//                g_struct_T *elStorage
// Return Type  : void
//
void c_calculateTimoshenkoElementIni(const double input_xloc[2], const double
  input_sectionProps_twist[2], const double input_sectionProps_rhoA[2], const
  double input_sectionProps_EIyy[2], const double input_sectionProps_EIzz[2],
  const double input_sectionProps_GJ[2], const double input_sectionProps_EA[2],
  const double input_sectionProps_rhoIyy[2], const double
  input_sectionProps_rhoIzz[2], const double input_sectionProps_rhoJ[2], const
  double input_sectionProps_zcm[2], const double input_sectionProps_ycm[2],
  double input_sweepAngle, double input_coneAngle, double input_rollAngle, const
  double input_x[2], const double input_y[2], const double input_z[2], boolean_T
  input_concMassFlag, const double input_concMass[8], g_struct_T *elStorage)
{
  double K11_idx_0;
  double K12_idx_0;
  double K13_idx_0;
  double K14_idx_0;
  double K15_idx_0;
  double K16_idx_0;
  double K22_idx_0;
  double K24_idx_0;
  double K33_idx_0;
  double K34_idx_0;
  double K44_idx_0;
  double K45_idx_0;
  double K46_idx_0;
  double K55_idx_0;
  double K56_idx_0;
  double K66_idx_0;
  double S11_idx_0;
  double S12_idx_0;
  double S13_idx_0;
  double S14_1_idx_0;
  double S14_2_idx_0;
  double S15_idx_0;
  double S16_idx_0;
  double S22_idx_0;
  double S23_idx_0;
  double S24_1_idx_0;
  double S24_2_idx_0;
  double S25_idx_0;
  double S26_idx_0;
  double S33_idx_0;
  double S34_1_idx_0;
  double S34_2_idx_0;
  double S35_idx_0;
  double S36_idx_0;
  double S44_1_idx_0;
  double S44_2_idx_0;
  double S44_3_idx_0;
  double S45_1_idx_0;
  double S45_2_idx_0;
  double S46_1_idx_0;
  double S46_2_idx_0;
  double S55_idx_0;
  double S56_idx_0;
  double S66_idx_0;
  double M11_idx_0;
  double M15_idx_0;
  double M16_idx_0;
  double M22_idx_0;
  double M24_idx_0;
  double M33_idx_0;
  double M34_idx_0;
  double M44_idx_0;
  double M55_idx_0;
  double M56_idx_0;
  double M66_idx_0;
  double C12_idx_0;
  double C13_idx_0;
  double C14_1_idx_0;
  double C14_2_idx_0;
  double C23_idx_0;
  double C24_idx_0;
  double C34_idx_0;
  double C25_idx_0;
  double C26_idx_0;
  double C35_idx_0;
  double C36_idx_0;
  double C45_1_idx_0;
  double C45_2_idx_0;
  double C46_1_idx_0;
  double C46_2_idx_0;
  double K11_idx_1;
  double K12_idx_1;
  double K13_idx_1;
  double K14_idx_1;
  double K15_idx_1;
  double K16_idx_1;
  double K22_idx_1;
  double K24_idx_1;
  double K33_idx_1;
  double K34_idx_1;
  double K44_idx_1;
  double K45_idx_1;
  double K46_idx_1;
  double K55_idx_1;
  double K56_idx_1;
  double K66_idx_1;
  double S11_idx_1;
  double S12_idx_1;
  double S13_idx_1;
  double S14_1_idx_1;
  double S14_2_idx_1;
  double S15_idx_1;
  double S16_idx_1;
  double S22_idx_1;
  double S23_idx_1;
  double S24_1_idx_1;
  double S24_2_idx_1;
  double S25_idx_1;
  double S26_idx_1;
  double S33_idx_1;
  double S34_1_idx_1;
  double S34_2_idx_1;
  double S35_idx_1;
  double S36_idx_1;
  double S44_1_idx_1;
  double S44_2_idx_1;
  double S44_3_idx_1;
  double S45_1_idx_1;
  double S45_2_idx_1;
  double S46_1_idx_1;
  double S46_2_idx_1;
  double S55_idx_1;
  double S56_idx_1;
  double S66_idx_1;
  double M11_idx_1;
  double M15_idx_1;
  double M16_idx_1;
  double M22_idx_1;
  double M24_idx_1;
  double M33_idx_1;
  double M34_idx_1;
  double M44_idx_1;
  double M55_idx_1;
  double M56_idx_1;
  double M66_idx_1;
  double C12_idx_1;
  double C13_idx_1;
  double C14_1_idx_1;
  double C14_2_idx_1;
  double C23_idx_1;
  double C24_idx_1;
  double C34_idx_1;
  double C25_idx_1;
  double C26_idx_1;
  double C35_idx_1;
  double C36_idx_1;
  double C45_1_idx_1;
  double C45_2_idx_1;
  double C46_1_idx_1;
  double C46_2_idx_1;
  double K11_idx_2;
  double K12_idx_2;
  double K13_idx_2;
  double K14_idx_2;
  double K15_idx_2;
  double K16_idx_2;
  double K22_idx_2;
  double K24_idx_2;
  double K33_idx_2;
  double K34_idx_2;
  double K44_idx_2;
  double K45_idx_2;
  double K46_idx_2;
  double K55_idx_2;
  double K56_idx_2;
  double K66_idx_2;
  double S11_idx_2;
  double S12_idx_2;
  double S13_idx_2;
  double S14_1_idx_2;
  double S14_2_idx_2;
  double S15_idx_2;
  double S16_idx_2;
  double S22_idx_2;
  double S23_idx_2;
  double S24_1_idx_2;
  double S24_2_idx_2;
  double S25_idx_2;
  double S26_idx_2;
  double S33_idx_2;
  double S34_1_idx_2;
  double S34_2_idx_2;
  double S35_idx_2;
  double S36_idx_2;
  double S44_1_idx_2;
  double S44_2_idx_2;
  double S44_3_idx_2;
  double S45_1_idx_2;
  double S45_2_idx_2;
  double S46_1_idx_2;
  double S46_2_idx_2;
  double S55_idx_2;
  double S56_idx_2;
  double S66_idx_2;
  double M11_idx_2;
  double M15_idx_2;
  double M16_idx_2;
  double M22_idx_2;
  double M24_idx_2;
  double M33_idx_2;
  double M34_idx_2;
  double M44_idx_2;
  double M55_idx_2;
  double M56_idx_2;
  double M66_idx_2;
  double C12_idx_2;
  double C13_idx_2;
  double C14_1_idx_2;
  double C14_2_idx_2;
  double C23_idx_2;
  double C24_idx_2;
  double C34_idx_2;
  double C25_idx_2;
  double C26_idx_2;
  double C35_idx_2;
  double C36_idx_2;
  double C45_1_idx_2;
  double C45_2_idx_2;
  double C46_1_idx_2;
  double C46_2_idx_2;
  double K11_idx_3;
  double K12_idx_3;
  double K13_idx_3;
  double K14_idx_3;
  double K15_idx_3;
  double K16_idx_3;
  double K22_idx_3;
  double K24_idx_3;
  double K33_idx_3;
  double K34_idx_3;
  double K44_idx_3;
  double K45_idx_3;
  double K46_idx_3;
  double K55_idx_3;
  double K56_idx_3;
  double K66_idx_3;
  double S11_idx_3;
  double S12_idx_3;
  double S13_idx_3;
  double S14_1_idx_3;
  double S14_2_idx_3;
  double S15_idx_3;
  double S16_idx_3;
  double S22_idx_3;
  double S23_idx_3;
  double S24_1_idx_3;
  double S24_2_idx_3;
  double S25_idx_3;
  double S26_idx_3;
  double S33_idx_3;
  double S34_1_idx_3;
  double S34_2_idx_3;
  double S35_idx_3;
  double S36_idx_3;
  double S44_1_idx_3;
  double S44_2_idx_3;
  double S44_3_idx_3;
  double S45_1_idx_3;
  double S45_2_idx_3;
  double S46_1_idx_3;
  double S46_2_idx_3;
  double S55_idx_3;
  double S56_idx_3;
  double S66_idx_3;
  double M11_idx_3;
  double M15_idx_3;
  double M16_idx_3;
  double M22_idx_3;
  double M24_idx_3;
  double M33_idx_3;
  double M34_idx_3;
  double M44_idx_3;
  double M55_idx_3;
  double M56_idx_3;
  double M66_idx_3;
  double C12_idx_3;
  double C13_idx_3;
  double C14_1_idx_3;
  double C14_2_idx_3;
  double C23_idx_3;
  double C24_idx_3;
  double C34_idx_3;
  double C25_idx_3;
  double C26_idx_3;
  double C35_idx_3;
  double C36_idx_3;
  double C45_1_idx_3;
  double C45_2_idx_3;
  double C46_1_idx_3;
  double C46_2_idx_3;
  double elementMass;
  double elementItens[9];
  double elxm[3];
  double lambda[144];
  int i;
  double N_data[2];
  int N_size[2];
  double p_N_x_data[2];
  int p_N_x_size[1];
  double Jac;
  double integrationFactor;
  double GA;
  double K14_tmp;
  double K12_tmp_tmp;
  double K11_tmp;
  double K55_tmp;
  double rhoIzz_tmp;
  double EA;
  double rhoA;
  double rhoIyy_tmp;
  double ycm;
  double zcm;
  double rhoIyy;
  double rhoIzz;
  double rhoIyz_tmp;
  double rhoIyz;
  double rhoJ;
  double EA_tmp;
  double K16_tmp;
  double b_EA_tmp;
  double M15_tmp_tmp;
  double M16_tmp_tmp;
  double M24_tmp_tmp;
  double M34_tmp_tmp;
  double M11_tmp;
  double b_M11_tmp;
  double b_M15_tmp_tmp;
  double b_M16_tmp_tmp;
  double b_M24_tmp_tmp;
  double b_M34_tmp_tmp;
  double c_M11_tmp;
  int b_i;
  double posLocal[3];
  double b_EA;
  double c_EA_tmp;
  double c_EA;
  double d_EA_tmp;
  double e_EA_tmp;
  double lamSlimTran[9];
  double d_EA;
  double e_EA;
  double b_lamSlimTran[9];
  double f_EA;
  int elementItens_tmp;
  double g_EA;
  double f_EA_tmp;
  double h_EA;
  double i_EA;
  double j_EA;
  double k_EA;
  double l_EA;
  double m_EA;
  double n_EA;
  double o_EA;
  double p_EA;
  double q_EA;
  double r_EA;
  double C12_tmp;
  double C13_tmp;
  double C14_1_tmp_tmp;
  double C14_2_tmp_tmp;
  double C23_tmp;
  double C24_tmp_tmp;
  double C34_tmp_tmp;
  double C45_1_tmp;
  double C45_2_tmp;
  double C46_1_tmp;
  double C46_2_tmp;
  double S11_tmp_tmp;
  double S14_2_tmp;
  double S15_tmp;
  double S16_tmp;
  double S66_tmp;

  //  CN2H           = eye(3,3);
  // --------------------------------------------
  // calculate quad points
  // Initialize element sub matrices and sub vectors
  //      F1 = zeros(numNodesPerEl,1);
  //      F3 = F1;
  //      F2 = F1;
  //      F4 = F1;
  //      F5 = F1;
  //      F6 = F1;
  K11_idx_0 = 0.0;
  K12_idx_0 = 0.0;
  K13_idx_0 = 0.0;
  K14_idx_0 = 0.0;
  K15_idx_0 = 0.0;
  K16_idx_0 = 0.0;
  K22_idx_0 = 0.0;
  elStorage->K23[0] = 0.0;
  K24_idx_0 = 0.0;
  elStorage->K25[0] = 0.0;
  K33_idx_0 = 0.0;
  K34_idx_0 = 0.0;
  elStorage->K36[0] = 0.0;
  K44_idx_0 = 0.0;
  K45_idx_0 = 0.0;
  K46_idx_0 = 0.0;
  K55_idx_0 = 0.0;
  K56_idx_0 = 0.0;
  K66_idx_0 = 0.0;
  S11_idx_0 = 0.0;
  S12_idx_0 = 0.0;
  S13_idx_0 = 0.0;
  S14_1_idx_0 = 0.0;
  S14_2_idx_0 = 0.0;
  S15_idx_0 = 0.0;
  S16_idx_0 = 0.0;
  S22_idx_0 = 0.0;
  S23_idx_0 = 0.0;
  S24_1_idx_0 = 0.0;
  S24_2_idx_0 = 0.0;
  S25_idx_0 = 0.0;
  S26_idx_0 = 0.0;
  S33_idx_0 = 0.0;
  S34_1_idx_0 = 0.0;
  S34_2_idx_0 = 0.0;
  S35_idx_0 = 0.0;
  S36_idx_0 = 0.0;
  S44_1_idx_0 = 0.0;
  S44_2_idx_0 = 0.0;
  S44_3_idx_0 = 0.0;
  S45_1_idx_0 = 0.0;
  S45_2_idx_0 = 0.0;
  S46_1_idx_0 = 0.0;
  S46_2_idx_0 = 0.0;
  S55_idx_0 = 0.0;
  S56_idx_0 = 0.0;
  S66_idx_0 = 0.0;
  M11_idx_0 = 0.0;
  M15_idx_0 = 0.0;
  M16_idx_0 = 0.0;
  M22_idx_0 = 0.0;
  M24_idx_0 = 0.0;
  M33_idx_0 = 0.0;
  M34_idx_0 = 0.0;
  M44_idx_0 = 0.0;
  M55_idx_0 = 0.0;
  M56_idx_0 = 0.0;
  M66_idx_0 = 0.0;
  C12_idx_0 = 0.0;
  C13_idx_0 = 0.0;
  C14_1_idx_0 = 0.0;
  C14_2_idx_0 = 0.0;
  C23_idx_0 = 0.0;
  C24_idx_0 = 0.0;
  C34_idx_0 = 0.0;
  C25_idx_0 = 0.0;
  C26_idx_0 = 0.0;
  C35_idx_0 = 0.0;
  C36_idx_0 = 0.0;
  C45_1_idx_0 = 0.0;
  C45_2_idx_0 = 0.0;
  C46_1_idx_0 = 0.0;
  C46_2_idx_0 = 0.0;
  K11_idx_1 = 0.0;
  K12_idx_1 = 0.0;
  K13_idx_1 = 0.0;
  K14_idx_1 = 0.0;
  K15_idx_1 = 0.0;
  K16_idx_1 = 0.0;
  K22_idx_1 = 0.0;
  elStorage->K23[1] = 0.0;
  K24_idx_1 = 0.0;
  elStorage->K25[1] = 0.0;
  K33_idx_1 = 0.0;
  K34_idx_1 = 0.0;
  elStorage->K36[1] = 0.0;
  K44_idx_1 = 0.0;
  K45_idx_1 = 0.0;
  K46_idx_1 = 0.0;
  K55_idx_1 = 0.0;
  K56_idx_1 = 0.0;
  K66_idx_1 = 0.0;
  S11_idx_1 = 0.0;
  S12_idx_1 = 0.0;
  S13_idx_1 = 0.0;
  S14_1_idx_1 = 0.0;
  S14_2_idx_1 = 0.0;
  S15_idx_1 = 0.0;
  S16_idx_1 = 0.0;
  S22_idx_1 = 0.0;
  S23_idx_1 = 0.0;
  S24_1_idx_1 = 0.0;
  S24_2_idx_1 = 0.0;
  S25_idx_1 = 0.0;
  S26_idx_1 = 0.0;
  S33_idx_1 = 0.0;
  S34_1_idx_1 = 0.0;
  S34_2_idx_1 = 0.0;
  S35_idx_1 = 0.0;
  S36_idx_1 = 0.0;
  S44_1_idx_1 = 0.0;
  S44_2_idx_1 = 0.0;
  S44_3_idx_1 = 0.0;
  S45_1_idx_1 = 0.0;
  S45_2_idx_1 = 0.0;
  S46_1_idx_1 = 0.0;
  S46_2_idx_1 = 0.0;
  S55_idx_1 = 0.0;
  S56_idx_1 = 0.0;
  S66_idx_1 = 0.0;
  M11_idx_1 = 0.0;
  M15_idx_1 = 0.0;
  M16_idx_1 = 0.0;
  M22_idx_1 = 0.0;
  M24_idx_1 = 0.0;
  M33_idx_1 = 0.0;
  M34_idx_1 = 0.0;
  M44_idx_1 = 0.0;
  M55_idx_1 = 0.0;
  M56_idx_1 = 0.0;
  M66_idx_1 = 0.0;
  C12_idx_1 = 0.0;
  C13_idx_1 = 0.0;
  C14_1_idx_1 = 0.0;
  C14_2_idx_1 = 0.0;
  C23_idx_1 = 0.0;
  C24_idx_1 = 0.0;
  C34_idx_1 = 0.0;
  C25_idx_1 = 0.0;
  C26_idx_1 = 0.0;
  C35_idx_1 = 0.0;
  C36_idx_1 = 0.0;
  C45_1_idx_1 = 0.0;
  C45_2_idx_1 = 0.0;
  C46_1_idx_1 = 0.0;
  C46_2_idx_1 = 0.0;
  K11_idx_2 = 0.0;
  K12_idx_2 = 0.0;
  K13_idx_2 = 0.0;
  K14_idx_2 = 0.0;
  K15_idx_2 = 0.0;
  K16_idx_2 = 0.0;
  K22_idx_2 = 0.0;
  elStorage->K23[2] = 0.0;
  K24_idx_2 = 0.0;
  elStorage->K25[2] = 0.0;
  K33_idx_2 = 0.0;
  K34_idx_2 = 0.0;
  elStorage->K36[2] = 0.0;
  K44_idx_2 = 0.0;
  K45_idx_2 = 0.0;
  K46_idx_2 = 0.0;
  K55_idx_2 = 0.0;
  K56_idx_2 = 0.0;
  K66_idx_2 = 0.0;
  S11_idx_2 = 0.0;
  S12_idx_2 = 0.0;
  S13_idx_2 = 0.0;
  S14_1_idx_2 = 0.0;
  S14_2_idx_2 = 0.0;
  S15_idx_2 = 0.0;
  S16_idx_2 = 0.0;
  S22_idx_2 = 0.0;
  S23_idx_2 = 0.0;
  S24_1_idx_2 = 0.0;
  S24_2_idx_2 = 0.0;
  S25_idx_2 = 0.0;
  S26_idx_2 = 0.0;
  S33_idx_2 = 0.0;
  S34_1_idx_2 = 0.0;
  S34_2_idx_2 = 0.0;
  S35_idx_2 = 0.0;
  S36_idx_2 = 0.0;
  S44_1_idx_2 = 0.0;
  S44_2_idx_2 = 0.0;
  S44_3_idx_2 = 0.0;
  S45_1_idx_2 = 0.0;
  S45_2_idx_2 = 0.0;
  S46_1_idx_2 = 0.0;
  S46_2_idx_2 = 0.0;
  S55_idx_2 = 0.0;
  S56_idx_2 = 0.0;
  S66_idx_2 = 0.0;
  M11_idx_2 = 0.0;
  M15_idx_2 = 0.0;
  M16_idx_2 = 0.0;
  M22_idx_2 = 0.0;
  M24_idx_2 = 0.0;
  M33_idx_2 = 0.0;
  M34_idx_2 = 0.0;
  M44_idx_2 = 0.0;
  M55_idx_2 = 0.0;
  M56_idx_2 = 0.0;
  M66_idx_2 = 0.0;
  C12_idx_2 = 0.0;
  C13_idx_2 = 0.0;
  C14_1_idx_2 = 0.0;
  C14_2_idx_2 = 0.0;
  C23_idx_2 = 0.0;
  C24_idx_2 = 0.0;
  C34_idx_2 = 0.0;
  C25_idx_2 = 0.0;
  C26_idx_2 = 0.0;
  C35_idx_2 = 0.0;
  C36_idx_2 = 0.0;
  C45_1_idx_2 = 0.0;
  C45_2_idx_2 = 0.0;
  C46_1_idx_2 = 0.0;
  C46_2_idx_2 = 0.0;
  K11_idx_3 = 0.0;
  K12_idx_3 = 0.0;
  K13_idx_3 = 0.0;
  K14_idx_3 = 0.0;
  K15_idx_3 = 0.0;
  K16_idx_3 = 0.0;
  K22_idx_3 = 0.0;
  elStorage->K23[3] = 0.0;
  K24_idx_3 = 0.0;
  elStorage->K25[3] = 0.0;
  K33_idx_3 = 0.0;
  K34_idx_3 = 0.0;
  elStorage->K36[3] = 0.0;
  K44_idx_3 = 0.0;
  K45_idx_3 = 0.0;
  K46_idx_3 = 0.0;
  K55_idx_3 = 0.0;
  K56_idx_3 = 0.0;
  K66_idx_3 = 0.0;
  S11_idx_3 = 0.0;
  S12_idx_3 = 0.0;
  S13_idx_3 = 0.0;
  S14_1_idx_3 = 0.0;
  S14_2_idx_3 = 0.0;
  S15_idx_3 = 0.0;
  S16_idx_3 = 0.0;
  S22_idx_3 = 0.0;
  S23_idx_3 = 0.0;
  S24_1_idx_3 = 0.0;
  S24_2_idx_3 = 0.0;
  S25_idx_3 = 0.0;
  S26_idx_3 = 0.0;
  S33_idx_3 = 0.0;
  S34_1_idx_3 = 0.0;
  S34_2_idx_3 = 0.0;
  S35_idx_3 = 0.0;
  S36_idx_3 = 0.0;
  S44_1_idx_3 = 0.0;
  S44_2_idx_3 = 0.0;
  S44_3_idx_3 = 0.0;
  S45_1_idx_3 = 0.0;
  S45_2_idx_3 = 0.0;
  S46_1_idx_3 = 0.0;
  S46_2_idx_3 = 0.0;
  S55_idx_3 = 0.0;
  S56_idx_3 = 0.0;
  S66_idx_3 = 0.0;
  M11_idx_3 = 0.0;
  M15_idx_3 = 0.0;
  M16_idx_3 = 0.0;
  M22_idx_3 = 0.0;
  M24_idx_3 = 0.0;
  M33_idx_3 = 0.0;
  M34_idx_3 = 0.0;
  M44_idx_3 = 0.0;
  M55_idx_3 = 0.0;
  M56_idx_3 = 0.0;
  M66_idx_3 = 0.0;
  C12_idx_3 = 0.0;
  C13_idx_3 = 0.0;
  C14_1_idx_3 = 0.0;
  C14_2_idx_3 = 0.0;
  C23_idx_3 = 0.0;
  C24_idx_3 = 0.0;
  C34_idx_3 = 0.0;
  C25_idx_3 = 0.0;
  C26_idx_3 = 0.0;
  C35_idx_3 = 0.0;
  C36_idx_3 = 0.0;
  C45_1_idx_3 = 0.0;
  C45_2_idx_3 = 0.0;
  C46_1_idx_3 = 0.0;
  C46_2_idx_3 = 0.0;
  elementMass = 0.0;
  std::memset(&elementItens[0], 0, 9U * sizeof(double));
  elxm[0] = 0.0;
  elxm[1] = 0.0;
  elxm[2] = 0.0;

  // Sort displacement vector
  // Written for 2 node element with 6 dof per node
  calculateLambda(input_sweepAngle * 3.1415926535897931 / 180.0, input_coneAngle
                  * 3.1415926535897931 / 180.0, (input_rollAngle + 0.5 *
    (input_sectionProps_twist[0] + input_sectionProps_twist[1])) *
                  3.1415926535897931 / 180.0, lambda);

  // Integration loop
  for (i = 0; i < 4; i++) {
    // Calculate shape functions at quad point i
    calculateShapeFunctions(dv[i], input_xloc, N_data, N_size, p_N_x_data,
      p_N_x_size, &Jac);
    integrationFactor = Jac * dv1[i];

    // ..... interpolate for value at quad point .....
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct stiffness terms
    // int(Ey dA)    v bend - extension
    // int(Ez dA)    w bend - extension
    // int(Gz dA)    w bend - twist
    // int(Gz dA)    v bend - twist
    //  extension twist
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // struct mass terms
    // set to zero to deactivate nonlinearites from initial element calculations 
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // .... end interpolate value at quad points ........
    // adjust moments of inertia for offsets
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    EA = N_data[0] * input_sectionProps_EA[0] + N_data[1] *
      input_sectionProps_EA[1];
    rhoA = N_data[0] * input_sectionProps_rhoA[0] + N_data[1] *
      input_sectionProps_rhoA[1];
    ycm = N_data[0] * input_sectionProps_ycm[0] + N_data[1] *
      input_sectionProps_ycm[1];
    zcm = N_data[0] * input_sectionProps_zcm[0] + N_data[1] *
      input_sectionProps_zcm[1];
    rhoIyy_tmp = zcm * zcm;
    rhoIyy = (N_data[0] * input_sectionProps_rhoIyy[0] + N_data[1] *
              input_sectionProps_rhoIyy[1]) + rhoA * rhoIyy_tmp;

    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    rhoIzz_tmp = ycm * ycm;
    rhoIzz = (N_data[0] * input_sectionProps_rhoIzz[0] + N_data[1] *
              input_sectionProps_rhoIzz[1]) + rhoA * rhoIzz_tmp;
    rhoIyz_tmp = rhoA * ycm;
    rhoIyz = rhoIyz_tmp * zcm;

    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    K11_tmp = EA * p_N_x_data[0];
    K11_idx_0 += K11_tmp * p_N_x_data[0] * integrationFactor;
    K11_idx_2 += K11_tmp * p_N_x_data[1] * integrationFactor;
    K11_tmp = EA * p_N_x_data[1];
    K11_idx_1 += K11_tmp * p_N_x_data[0] * integrationFactor;
    K11_idx_3 += K11_tmp * p_N_x_data[1] * integrationFactor;
    rhoJ = (N_data[0] * input_sectionProps_rhoJ[0] + N_data[1] *
            input_sectionProps_rhoJ[1]) + rhoA * (rhoIzz_tmp + rhoIyy_tmp);

    // Calculate strutural stiffness sub matrices
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    EA_tmp = 0.5 * EA * 0.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    K12_tmp_tmp = EA_tmp * p_N_x_data[0];
    K12_idx_0 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K13_idx_0 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K14_tmp = 0.0 * p_N_x_data[0] * p_N_x_data[0] * integrationFactor;
    K14_idx_0 += K14_tmp;
    K15_idx_0 += K14_tmp;
    K11_tmp = -0.0 * p_N_x_data[0] * p_N_x_data[0] * integrationFactor;
    K16_idx_0 += K11_tmp;
    K22_idx_0 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K24_idx_0 += K11_tmp;
    K33_idx_0 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K34_idx_0 += K14_tmp;
    K12_idx_2 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K13_idx_2 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K14_tmp = 0.0 * p_N_x_data[0] * p_N_x_data[1] * integrationFactor;
    K14_idx_2 += K14_tmp;
    K15_idx_2 += K14_tmp;
    K16_tmp = -0.0 * p_N_x_data[0] * p_N_x_data[1] * integrationFactor;
    K16_idx_2 += K16_tmp;
    K22_idx_2 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K24_idx_2 += K16_tmp;
    K33_idx_2 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K34_idx_2 += K14_tmp;
    K12_tmp_tmp = EA_tmp * p_N_x_data[1];
    K12_idx_1 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K13_idx_1 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K14_tmp = 0.0 * p_N_x_data[1] * p_N_x_data[0] * integrationFactor;
    K14_idx_1 += K14_tmp;
    K15_idx_1 += K14_tmp;
    rhoIyy_tmp = -0.0 * p_N_x_data[1] * p_N_x_data[0] * integrationFactor;
    K16_idx_1 += rhoIyy_tmp;
    K22_idx_1 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K24_idx_1 += rhoIyy_tmp;
    K33_idx_1 += K12_tmp_tmp * p_N_x_data[0] * integrationFactor;
    K34_idx_1 += K14_tmp;
    K12_idx_3 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K13_idx_3 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K14_tmp = 0.0 * p_N_x_data[1] * p_N_x_data[1] * integrationFactor;
    K14_idx_3 += K14_tmp;
    K15_idx_3 += K14_tmp;
    rhoIzz_tmp = -0.0 * p_N_x_data[1] * p_N_x_data[1] * integrationFactor;
    K16_idx_3 += rhoIzz_tmp;
    K22_idx_3 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K24_idx_3 += rhoIzz_tmp;
    K33_idx_3 += K12_tmp_tmp * p_N_x_data[1] * integrationFactor;
    K34_idx_3 += K14_tmp;
    K12_tmp_tmp = N_data[0] * input_sectionProps_GJ[0] + N_data[1] *
      input_sectionProps_GJ[1];

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    Jac = K12_tmp_tmp * p_N_x_data[0];
    K44_idx_0 += Jac * p_N_x_data[0] * integrationFactor;
    GA = 0.0 * p_N_x_data[0] * N_data[0] * integrationFactor;
    K45_idx_0 += GA;
    K46_idx_0 += GA;
    K44_idx_2 += Jac * p_N_x_data[1] * integrationFactor;
    GA = 0.0 * p_N_x_data[0] * N_data[1] * integrationFactor;
    K45_idx_2 += GA;
    K46_idx_2 += GA;
    Jac = K12_tmp_tmp * p_N_x_data[1];
    K44_idx_1 += Jac * p_N_x_data[0] * integrationFactor;
    GA = 0.0 * p_N_x_data[1] * N_data[0] * integrationFactor;
    K45_idx_1 += GA;
    K46_idx_1 += GA;
    K44_idx_3 += Jac * p_N_x_data[1] * integrationFactor;
    GA = 0.0 * p_N_x_data[1] * N_data[1] * integrationFactor;
    K45_idx_3 += GA;
    K46_idx_3 += GA;
    K12_tmp_tmp = N_data[0] * input_sectionProps_EIyy[0] + N_data[1] *
      input_sectionProps_EIyy[1];

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    K55_tmp = K12_tmp_tmp * p_N_x_data[0];
    K55_idx_0 += K55_tmp * p_N_x_data[0] * integrationFactor;
    K56_idx_0 += K11_tmp;
    K55_idx_2 += K55_tmp * p_N_x_data[1] * integrationFactor;
    K56_idx_2 += K16_tmp;
    K55_tmp = K12_tmp_tmp * p_N_x_data[1];
    K55_idx_1 += K55_tmp * p_N_x_data[0] * integrationFactor;
    K56_idx_1 += rhoIyy_tmp;
    K55_idx_3 += K55_tmp * p_N_x_data[1] * integrationFactor;
    K56_idx_3 += rhoIzz_tmp;
    K12_tmp_tmp = N_data[0] * input_sectionProps_EIzz[0] + N_data[1] *
      input_sectionProps_EIzz[1];

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // Calculate structural mass sub matrices
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    EA_tmp = rhoA * zcm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    b_EA_tmp = -rhoA * ycm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    K14_tmp = -rhoA * zcm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // Calculate Centrifugal load vector and gravity load vector
    // eventually incorporate lambda into gp level to account for variable
    // twist
    // these are set to unity to get coefficients for omega components
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    // This function interpolates a value using distinct values at valNode
    // and the corresponding shape function N.
    Jac = K12_tmp_tmp * p_N_x_data[0];
    K66_idx_0 += Jac * p_N_x_data[0] * integrationFactor;
    GA = rhoA * N_data[0];
    K16_tmp = GA * N_data[0] * integrationFactor;
    M11_idx_0 += K16_tmp;
    M15_tmp_tmp = EA_tmp * N_data[0];
    M15_idx_0 += M15_tmp_tmp * N_data[0] * integrationFactor;
    M16_tmp_tmp = b_EA_tmp * N_data[0];
    M16_idx_0 += M16_tmp_tmp * N_data[0] * integrationFactor;
    M22_idx_0 += K16_tmp;
    M24_tmp_tmp = K14_tmp * N_data[0];
    M24_idx_0 += M24_tmp_tmp * N_data[0] * integrationFactor;
    M33_idx_0 += K16_tmp;
    M34_tmp_tmp = rhoIyz_tmp * N_data[0];
    M34_idx_0 += M34_tmp_tmp * N_data[0] * integrationFactor;
    rhoIyy_tmp = rhoJ * N_data[0];
    M44_idx_0 += rhoIyy_tmp * N_data[0] * integrationFactor;
    rhoIzz_tmp = rhoIyy * N_data[0];
    M55_idx_0 += rhoIzz_tmp * N_data[0] * integrationFactor;
    K11_tmp = -rhoIyz * N_data[0];
    M56_idx_0 += K11_tmp * N_data[0] * integrationFactor;
    K55_tmp = rhoIzz * N_data[0];
    M66_idx_0 += K55_tmp * N_data[0] * integrationFactor;
    K66_idx_2 += Jac * p_N_x_data[1] * integrationFactor;
    M11_tmp = GA * N_data[1] * integrationFactor;
    M11_idx_2 += M11_tmp;
    M15_idx_2 += M15_tmp_tmp * N_data[1] * integrationFactor;
    M16_idx_2 += M16_tmp_tmp * N_data[1] * integrationFactor;
    M22_idx_2 += M11_tmp;
    M24_idx_2 += M24_tmp_tmp * N_data[1] * integrationFactor;
    M33_idx_2 += M11_tmp;
    M34_idx_2 += M34_tmp_tmp * N_data[1] * integrationFactor;
    M44_idx_2 += rhoIyy_tmp * N_data[1] * integrationFactor;
    M55_idx_2 += rhoIzz_tmp * N_data[1] * integrationFactor;
    M56_idx_2 += K11_tmp * N_data[1] * integrationFactor;
    M66_idx_2 += K55_tmp * N_data[1] * integrationFactor;
    Jac = K12_tmp_tmp * p_N_x_data[1];
    K66_idx_1 += Jac * p_N_x_data[0] * integrationFactor;
    GA = rhoA * N_data[1];
    b_M11_tmp = GA * N_data[0] * integrationFactor;
    M11_idx_1 += b_M11_tmp;
    b_M15_tmp_tmp = EA_tmp * N_data[1];
    M15_idx_1 += b_M15_tmp_tmp * N_data[0] * integrationFactor;
    b_M16_tmp_tmp = b_EA_tmp * N_data[1];
    M16_idx_1 += b_M16_tmp_tmp * N_data[0] * integrationFactor;
    M22_idx_1 += b_M11_tmp;
    b_M24_tmp_tmp = K14_tmp * N_data[1];
    M24_idx_1 += b_M24_tmp_tmp * N_data[0] * integrationFactor;
    M33_idx_1 += b_M11_tmp;
    b_M34_tmp_tmp = rhoIyz_tmp * N_data[1];
    M34_idx_1 += b_M34_tmp_tmp * N_data[0] * integrationFactor;
    rhoIyy_tmp = rhoJ * N_data[1];
    M44_idx_1 += rhoIyy_tmp * N_data[0] * integrationFactor;
    rhoIzz_tmp = rhoIyy * N_data[1];
    M55_idx_1 += rhoIzz_tmp * N_data[0] * integrationFactor;
    K11_tmp = -rhoIyz * N_data[1];
    M56_idx_1 += K11_tmp * N_data[0] * integrationFactor;
    K55_tmp = rhoIzz * N_data[1];
    M66_idx_1 += K55_tmp * N_data[0] * integrationFactor;
    K66_idx_3 += Jac * p_N_x_data[1] * integrationFactor;
    c_M11_tmp = GA * N_data[1] * integrationFactor;
    M11_idx_3 += c_M11_tmp;
    M15_idx_3 += b_M15_tmp_tmp * N_data[1] * integrationFactor;
    M16_idx_3 += b_M16_tmp_tmp * N_data[1] * integrationFactor;
    M22_idx_3 += c_M11_tmp;
    M24_idx_3 += b_M24_tmp_tmp * N_data[1] * integrationFactor;
    M33_idx_3 += c_M11_tmp;
    M34_idx_3 += b_M34_tmp_tmp * N_data[1] * integrationFactor;
    M44_idx_3 += rhoIyy_tmp * N_data[1] * integrationFactor;
    M55_idx_3 += rhoIzz_tmp * N_data[1] * integrationFactor;
    M56_idx_3 += K11_tmp * N_data[1] * integrationFactor;
    M66_idx_3 += K55_tmp * N_data[1] * integrationFactor;
    K12_tmp_tmp = N_data[0] * input_x[0] + N_data[1] * input_x[1];
    Jac = N_data[0] * input_y[0] + N_data[1] * input_y[1];
    GA = N_data[0] * input_z[0] + N_data[1] * input_z[1];
    for (b_i = 0; b_i < 3; b_i++) {
      posLocal[b_i] = (lambda[b_i] * K12_tmp_tmp + lambda[b_i + 12] * Jac) +
        lambda[b_i + 24] * GA;
    }

    //         g=9.81; %gravitational acceleration [m/s^2]
    //         a_x = 0; %acceleration of body in x and y (hardwired to zero for now) 
    //         a_y = 0;
    //         a_z = -g;
    //         fx = rhoA*a_x; %let these loads be defined in the inertial frame
    //         fy = rhoA*a_y;
    //         fz = rhoA*a_z;
    //         rvec = [ 0; ycm; zcm];
    //
    //         fi_hub = CN2H*[fx;fy;fz];
    //
    //         disLoadgpLocal = lambda(1:3,1:3)*fi_hub;
    //         disMomentgp = cross(rvec,disLoadgpLocal);
    //         f1 = rhoA*((O2^2 + O3^2)*xbarlocal - O1*O2*ybarlocal - O1*O3*zbarlocal) - disLoadgpLocal(1);    %omega dot loading not 
    //         [F1] = calculateVec1(f1,integrationFactor,N1,F1);
    //         f2 = rhoA*((O1^2+O3^2)*ybarlocal - zbarlocal*O2*O3 - xbarlocal*O1*O2) - disLoadgpLocal(2); 
    //         [F2] = calculateVec1(f2,integrationFactor,N2,F2);
    //         f3 = rhoA*((O1^2+O2^2)*zbarlocal - O3*O1*xbarlocal - O2*O3*ybarlocal) - disLoadgpLocal(3); 
    //         [F3] = calculateVec1(f3,integrationFactor,N3,F3);
    //         f4 = rhoA*(xbarlocal*(O1*O2*zcm - ycm*O1*O3)-ybarlocal*(ycm*O2*O3 + zcm*(O1^2+O3^2))... 
    //                    + zbarlocal*(ycm*(O1^2+O2^2)+zcm*O2*O3)) - disMomentgp(1); 
    //         [F4] = calculateVec1(f4,integrationFactor,N4,F4);
    //         f5 = rhoA*zcm*(xbarlocal*(O2^2+O3^2) - ybarlocal*O1*O2 - zbarlocal*O1*O3) - disMomentgp(2); 
    //         [F5] = calculateVec1(f5,integrationFactor,N5,F5);
    //         f6 = rhoA*ycm*((O1*O3*zbarlocal + O1*O2*ybarlocal)-(xbarlocal*(O2^2+O3^2))) - disMomentgp(3); 
    //         [F6] = calculateVec1(f6,integrationFactor,N6,F6);
    //
    // Gyric matrix (Coriolis)
    EA = -2.0 * rhoA;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    b_EA = 2.0 * rhoA;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    b_EA_tmp = 2.0 * rhoA * ycm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    c_EA_tmp = 2.0 * rhoA * zcm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    c_EA = -2.0 * rhoA;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    d_EA_tmp = -2.0 * rhoA * ycm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    e_EA_tmp = -2.0 * rhoA * zcm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    d_EA = -2.0 * rhoIyy;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    e_EA = -2.0 * rhoIyz;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    f_EA = 2.0 * rhoIzz;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    g_EA = 2.0 * rhoIyz;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // Spin softening matrix
    f_EA_tmp = -rhoA * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    h_EA = rhoA * -zcm;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    i_EA = K14_tmp * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    j_EA = rhoIyz_tmp * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    k_EA = EA_tmp * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    l_EA = -rhoA * (ycm * 2.0);

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    m_EA = -(rhoIyy * 2.0);

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    n_EA = -(rhoIzz * 2.0);

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    o_EA = -(2.0 * rhoIyz);

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    p_EA = -rhoIyy * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    q_EA = rhoIyz * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    r_EA = -rhoIzz * 2.0;

    // This function is a general routine to calculate an element matrix
    // Element calculation functions---------------------------------
    C12_tmp = EA * N_data[0];
    C12_idx_0 += C12_tmp * N_data[0] * integrationFactor;
    C13_tmp = b_EA * N_data[0];
    C13_idx_0 += C13_tmp * N_data[0] * integrationFactor;
    C14_1_tmp_tmp = b_EA_tmp * N_data[0];
    C14_1_idx_0 += C14_1_tmp_tmp * N_data[0] * integrationFactor;
    C14_2_tmp_tmp = c_EA_tmp * N_data[0];
    C14_2_idx_0 += C14_2_tmp_tmp * N_data[0] * integrationFactor;
    C23_tmp = c_EA * N_data[0];
    C23_idx_0 += C23_tmp * N_data[0] * integrationFactor;
    C24_tmp_tmp = d_EA_tmp * N_data[0];
    C24_idx_0 += C24_tmp_tmp * N_data[0] * integrationFactor;
    C25_idx_0 += C14_2_tmp_tmp * N_data[0] * integrationFactor;
    C26_idx_0 += C24_tmp_tmp * N_data[0] * integrationFactor;
    C34_tmp_tmp = e_EA_tmp * N_data[0];
    C34_idx_0 += C34_tmp_tmp * N_data[0] * integrationFactor;
    C35_idx_0 += C34_tmp_tmp * N_data[0] * integrationFactor;
    C36_idx_0 += C14_1_tmp_tmp * N_data[0] * integrationFactor;
    C45_1_tmp = d_EA * N_data[0];
    C45_1_idx_0 += C45_1_tmp * N_data[0] * integrationFactor;
    C45_2_tmp = e_EA * N_data[0];
    C45_2_idx_0 += C45_2_tmp * N_data[0] * integrationFactor;
    C46_1_tmp = f_EA * N_data[0];
    C46_1_idx_0 += C46_1_tmp * N_data[0] * integrationFactor;
    C46_2_tmp = g_EA * N_data[0];
    C46_2_idx_0 += C46_2_tmp * N_data[0] * integrationFactor;
    S11_tmp_tmp = f_EA_tmp * N_data[0];
    S11_idx_0 += S11_tmp_tmp * N_data[0] * integrationFactor;
    S12_idx_0 += K16_tmp;
    S13_idx_0 += K16_tmp;
    S14_1_idx_0 += M34_tmp_tmp * N_data[0] * integrationFactor;
    S14_2_tmp = h_EA * N_data[0];
    S14_2_idx_0 += S14_2_tmp * N_data[0] * integrationFactor;
    S15_tmp = i_EA * N_data[0];
    S15_idx_0 += S15_tmp * N_data[0] * integrationFactor;
    S16_tmp = j_EA * N_data[0];
    S16_idx_0 += S16_tmp * N_data[0] * integrationFactor;
    S22_idx_0 += S11_tmp_tmp * N_data[0] * integrationFactor;
    S23_idx_0 += K16_tmp;
    Jac = k_EA * N_data[0];
    S24_1_idx_0 += Jac * N_data[0] * integrationFactor;
    S24_2_idx_0 += M34_tmp_tmp * N_data[0] * integrationFactor;
    S25_idx_0 += M15_tmp_tmp * N_data[0] * integrationFactor;
    S26_idx_0 += M16_tmp_tmp * N_data[0] * integrationFactor;
    S33_idx_0 += S11_tmp_tmp * N_data[0] * integrationFactor;
    GA = l_EA * N_data[0];
    S34_1_idx_0 += GA * N_data[0] * integrationFactor;
    S34_2_idx_0 += M24_tmp_tmp * N_data[0] * integrationFactor;
    S35_idx_0 += M15_tmp_tmp * N_data[0] * integrationFactor;
    S36_idx_0 += M16_tmp_tmp * N_data[0] * integrationFactor;
    rhoIyy_tmp = m_EA * N_data[0];
    S44_1_idx_0 += rhoIyy_tmp * N_data[0] * integrationFactor;
    rhoIzz_tmp = n_EA * N_data[0];
    S44_2_idx_0 += rhoIzz_tmp * N_data[0] * integrationFactor;
    K11_tmp = o_EA * N_data[0];
    S44_3_idx_0 += K11_tmp * N_data[0] * integrationFactor;
    K55_tmp = rhoIyz * N_data[0];
    K12_tmp_tmp = K55_tmp * N_data[0] * integrationFactor;
    S45_1_idx_0 += K12_tmp_tmp;
    K14_tmp = -rhoIyy * N_data[0];
    S45_2_idx_0 += K14_tmp * N_data[0] * integrationFactor;
    S46_1_idx_0 += K12_tmp_tmp;
    K16_tmp = -rhoIzz * N_data[0];
    S46_2_idx_0 += K16_tmp * N_data[0] * integrationFactor;
    EA_tmp = p_EA * N_data[0];
    S55_idx_0 += EA_tmp * N_data[0] * integrationFactor;
    rhoIyz_tmp = q_EA * N_data[0];
    S56_idx_0 += rhoIyz_tmp * N_data[0] * integrationFactor;
    S66_tmp = r_EA * N_data[0];
    S66_idx_0 += S66_tmp * N_data[0] * integrationFactor;
    C12_idx_2 += C12_tmp * N_data[1] * integrationFactor;
    C13_idx_2 += C13_tmp * N_data[1] * integrationFactor;
    C14_1_idx_2 += C14_1_tmp_tmp * N_data[1] * integrationFactor;
    C14_2_idx_2 += C14_2_tmp_tmp * N_data[1] * integrationFactor;
    C23_idx_2 += C23_tmp * N_data[1] * integrationFactor;
    C24_idx_2 += C24_tmp_tmp * N_data[1] * integrationFactor;
    C25_idx_2 += C14_2_tmp_tmp * N_data[1] * integrationFactor;
    C26_idx_2 += C24_tmp_tmp * N_data[1] * integrationFactor;
    C34_idx_2 += C34_tmp_tmp * N_data[1] * integrationFactor;
    C35_idx_2 += C34_tmp_tmp * N_data[1] * integrationFactor;
    C36_idx_2 += C14_1_tmp_tmp * N_data[1] * integrationFactor;
    C45_1_idx_2 += C45_1_tmp * N_data[1] * integrationFactor;
    C45_2_idx_2 += C45_2_tmp * N_data[1] * integrationFactor;
    C46_1_idx_2 += C46_1_tmp * N_data[1] * integrationFactor;
    C46_2_idx_2 += C46_2_tmp * N_data[1] * integrationFactor;
    S11_idx_2 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S12_idx_2 += M11_tmp;
    S13_idx_2 += M11_tmp;
    S14_1_idx_2 += M34_tmp_tmp * N_data[1] * integrationFactor;
    S14_2_idx_2 += S14_2_tmp * N_data[1] * integrationFactor;
    S15_idx_2 += S15_tmp * N_data[1] * integrationFactor;
    S16_idx_2 += S16_tmp * N_data[1] * integrationFactor;
    S22_idx_2 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S23_idx_2 += M11_tmp;
    S24_1_idx_2 += Jac * N_data[1] * integrationFactor;
    S24_2_idx_2 += M34_tmp_tmp * N_data[1] * integrationFactor;
    S25_idx_2 += M15_tmp_tmp * N_data[1] * integrationFactor;
    S26_idx_2 += M16_tmp_tmp * N_data[1] * integrationFactor;
    S33_idx_2 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S34_1_idx_2 += GA * N_data[1] * integrationFactor;
    S34_2_idx_2 += M24_tmp_tmp * N_data[1] * integrationFactor;
    S35_idx_2 += M15_tmp_tmp * N_data[1] * integrationFactor;
    S36_idx_2 += M16_tmp_tmp * N_data[1] * integrationFactor;
    S44_1_idx_2 += rhoIyy_tmp * N_data[1] * integrationFactor;
    S44_2_idx_2 += rhoIzz_tmp * N_data[1] * integrationFactor;
    S44_3_idx_2 += K11_tmp * N_data[1] * integrationFactor;
    K12_tmp_tmp = K55_tmp * N_data[1] * integrationFactor;
    S45_1_idx_2 += K12_tmp_tmp;
    S45_2_idx_2 += K14_tmp * N_data[1] * integrationFactor;
    S46_1_idx_2 += K12_tmp_tmp;
    S46_2_idx_2 += K16_tmp * N_data[1] * integrationFactor;
    S55_idx_2 += EA_tmp * N_data[1] * integrationFactor;
    S56_idx_2 += rhoIyz_tmp * N_data[1] * integrationFactor;
    S66_idx_2 += S66_tmp * N_data[1] * integrationFactor;
    C12_tmp = EA * N_data[1];
    C12_idx_1 += C12_tmp * N_data[0] * integrationFactor;
    C13_tmp = b_EA * N_data[1];
    C13_idx_1 += C13_tmp * N_data[0] * integrationFactor;
    C14_1_tmp_tmp = b_EA_tmp * N_data[1];
    C14_1_idx_1 += C14_1_tmp_tmp * N_data[0] * integrationFactor;
    C14_2_tmp_tmp = c_EA_tmp * N_data[1];
    C14_2_idx_1 += C14_2_tmp_tmp * N_data[0] * integrationFactor;
    C23_tmp = c_EA * N_data[1];
    C23_idx_1 += C23_tmp * N_data[0] * integrationFactor;
    C24_tmp_tmp = d_EA_tmp * N_data[1];
    C24_idx_1 += C24_tmp_tmp * N_data[0] * integrationFactor;
    C25_idx_1 += C14_2_tmp_tmp * N_data[0] * integrationFactor;
    C26_idx_1 += C24_tmp_tmp * N_data[0] * integrationFactor;
    C34_tmp_tmp = e_EA_tmp * N_data[1];
    C34_idx_1 += C34_tmp_tmp * N_data[0] * integrationFactor;
    C35_idx_1 += C34_tmp_tmp * N_data[0] * integrationFactor;
    C36_idx_1 += C14_1_tmp_tmp * N_data[0] * integrationFactor;
    C45_1_tmp = d_EA * N_data[1];
    C45_1_idx_1 += C45_1_tmp * N_data[0] * integrationFactor;
    C45_2_tmp = e_EA * N_data[1];
    C45_2_idx_1 += C45_2_tmp * N_data[0] * integrationFactor;
    C46_1_tmp = f_EA * N_data[1];
    C46_1_idx_1 += C46_1_tmp * N_data[0] * integrationFactor;
    C46_2_tmp = g_EA * N_data[1];
    C46_2_idx_1 += C46_2_tmp * N_data[0] * integrationFactor;
    S11_tmp_tmp = f_EA_tmp * N_data[1];
    S11_idx_1 += S11_tmp_tmp * N_data[0] * integrationFactor;
    S12_idx_1 += b_M11_tmp;
    S13_idx_1 += b_M11_tmp;
    S14_1_idx_1 += b_M34_tmp_tmp * N_data[0] * integrationFactor;
    S14_2_tmp = h_EA * N_data[1];
    S14_2_idx_1 += S14_2_tmp * N_data[0] * integrationFactor;
    S15_tmp = i_EA * N_data[1];
    S15_idx_1 += S15_tmp * N_data[0] * integrationFactor;
    S16_tmp = j_EA * N_data[1];
    S16_idx_1 += S16_tmp * N_data[0] * integrationFactor;
    S22_idx_1 += S11_tmp_tmp * N_data[0] * integrationFactor;
    S23_idx_1 += b_M11_tmp;
    Jac = k_EA * N_data[1];
    S24_1_idx_1 += Jac * N_data[0] * integrationFactor;
    S24_2_idx_1 += b_M34_tmp_tmp * N_data[0] * integrationFactor;
    S25_idx_1 += b_M15_tmp_tmp * N_data[0] * integrationFactor;
    S26_idx_1 += b_M16_tmp_tmp * N_data[0] * integrationFactor;
    S33_idx_1 += S11_tmp_tmp * N_data[0] * integrationFactor;
    GA = l_EA * N_data[1];
    S34_1_idx_1 += GA * N_data[0] * integrationFactor;
    S34_2_idx_1 += b_M24_tmp_tmp * N_data[0] * integrationFactor;
    S35_idx_1 += b_M15_tmp_tmp * N_data[0] * integrationFactor;
    S36_idx_1 += b_M16_tmp_tmp * N_data[0] * integrationFactor;
    rhoIyy_tmp = m_EA * N_data[1];
    S44_1_idx_1 += rhoIyy_tmp * N_data[0] * integrationFactor;
    rhoIzz_tmp = n_EA * N_data[1];
    S44_2_idx_1 += rhoIzz_tmp * N_data[0] * integrationFactor;
    K11_tmp = o_EA * N_data[1];
    S44_3_idx_1 += K11_tmp * N_data[0] * integrationFactor;
    K55_tmp = rhoIyz * N_data[1];
    K12_tmp_tmp = K55_tmp * N_data[0] * integrationFactor;
    S45_1_idx_1 += K12_tmp_tmp;
    K14_tmp = -rhoIyy * N_data[1];
    S45_2_idx_1 += K14_tmp * N_data[0] * integrationFactor;
    S46_1_idx_1 += K12_tmp_tmp;
    K16_tmp = -rhoIzz * N_data[1];
    S46_2_idx_1 += K16_tmp * N_data[0] * integrationFactor;
    EA_tmp = p_EA * N_data[1];
    S55_idx_1 += EA_tmp * N_data[0] * integrationFactor;
    rhoIyz_tmp = q_EA * N_data[1];
    S56_idx_1 += rhoIyz_tmp * N_data[0] * integrationFactor;
    S66_tmp = r_EA * N_data[1];
    S66_idx_1 += S66_tmp * N_data[0] * integrationFactor;
    C12_idx_3 += C12_tmp * N_data[1] * integrationFactor;
    C13_idx_3 += C13_tmp * N_data[1] * integrationFactor;
    C14_1_idx_3 += C14_1_tmp_tmp * N_data[1] * integrationFactor;
    C14_2_idx_3 += C14_2_tmp_tmp * N_data[1] * integrationFactor;
    C23_idx_3 += C23_tmp * N_data[1] * integrationFactor;
    C24_idx_3 += C24_tmp_tmp * N_data[1] * integrationFactor;
    C25_idx_3 += C14_2_tmp_tmp * N_data[1] * integrationFactor;
    C26_idx_3 += C24_tmp_tmp * N_data[1] * integrationFactor;
    C34_idx_3 += C34_tmp_tmp * N_data[1] * integrationFactor;
    C35_idx_3 += C34_tmp_tmp * N_data[1] * integrationFactor;
    C36_idx_3 += C14_1_tmp_tmp * N_data[1] * integrationFactor;
    C45_1_idx_3 += C45_1_tmp * N_data[1] * integrationFactor;
    C45_2_idx_3 += C45_2_tmp * N_data[1] * integrationFactor;
    C46_1_idx_3 += C46_1_tmp * N_data[1] * integrationFactor;
    C46_2_idx_3 += C46_2_tmp * N_data[1] * integrationFactor;
    S11_idx_3 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S12_idx_3 += c_M11_tmp;
    S13_idx_3 += c_M11_tmp;
    S14_1_idx_3 += b_M34_tmp_tmp * N_data[1] * integrationFactor;
    S14_2_idx_3 += S14_2_tmp * N_data[1] * integrationFactor;
    S15_idx_3 += S15_tmp * N_data[1] * integrationFactor;
    S16_idx_3 += S16_tmp * N_data[1] * integrationFactor;
    S22_idx_3 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S23_idx_3 += c_M11_tmp;
    S24_1_idx_3 += Jac * N_data[1] * integrationFactor;
    S24_2_idx_3 += b_M34_tmp_tmp * N_data[1] * integrationFactor;
    S25_idx_3 += b_M15_tmp_tmp * N_data[1] * integrationFactor;
    S26_idx_3 += b_M16_tmp_tmp * N_data[1] * integrationFactor;
    S33_idx_3 += S11_tmp_tmp * N_data[1] * integrationFactor;
    S34_1_idx_3 += GA * N_data[1] * integrationFactor;
    S34_2_idx_3 += b_M24_tmp_tmp * N_data[1] * integrationFactor;
    S35_idx_3 += b_M15_tmp_tmp * N_data[1] * integrationFactor;
    S36_idx_3 += b_M16_tmp_tmp * N_data[1] * integrationFactor;
    S44_1_idx_3 += rhoIyy_tmp * N_data[1] * integrationFactor;
    S44_2_idx_3 += rhoIzz_tmp * N_data[1] * integrationFactor;
    S44_3_idx_3 += K11_tmp * N_data[1] * integrationFactor;
    K12_tmp_tmp = K55_tmp * N_data[1] * integrationFactor;
    S45_1_idx_3 += K12_tmp_tmp;
    S45_2_idx_3 += K14_tmp * N_data[1] * integrationFactor;
    S46_1_idx_3 += K12_tmp_tmp;
    S46_2_idx_3 += K16_tmp * N_data[1] * integrationFactor;
    S55_idx_3 += EA_tmp * N_data[1] * integrationFactor;
    S56_idx_3 += rhoIyz_tmp * N_data[1] * integrationFactor;
    S66_idx_3 += S66_tmp * N_data[1] * integrationFactor;
    calculateElementMass(rhoA, rhoIyy, rhoIzz, rhoIyz, rhoJ, ycm, zcm, posLocal
                         [0], posLocal[1], posLocal[2], integrationFactor,
                         &elementMass, elementItens, elxm);
  }

  // ==========================================================
  // Reduced integration loop
  // Calculate shape functions at quad point i
  b_calculateShapeFunctions(input_xloc, N_data, N_size, p_N_x_data, p_N_x_size,
    &Jac);
  integrationFactor = Jac * 2.0;

  // ..... interpolate for value at quad point .....
  // struct stiffness terms
  // This function interpolates a value using distinct values at valNode
  // and the corresponding shape function N.
  GA = (N_data[0] * input_sectionProps_EA[0] + N_data[1] *
        input_sectionProps_EA[1]) / 2.6 * 5.0 / 6.0;

  // .... end interpolate value at quad points ........
  // Calculate strutural stiffness sub matrices
  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  K14_tmp = GA * p_N_x_data[0];
  Jac = K14_tmp * p_N_x_data[0] * integrationFactor;
  elStorage->K22[0] = K22_idx_0 + Jac;
  K12_tmp_tmp = K14_tmp * p_N_x_data[1] * integrationFactor;
  elStorage->K22[2] = K22_idx_2 + K12_tmp_tmp;
  K11_tmp = GA * p_N_x_data[1];
  K55_tmp = K11_tmp * p_N_x_data[0] * integrationFactor;
  elStorage->K22[1] = K22_idx_1 + K55_tmp;
  rhoIzz_tmp = K11_tmp * p_N_x_data[1] * integrationFactor;
  elStorage->K22[3] = K22_idx_3 + rhoIzz_tmp;

  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  rhoIyy_tmp = -GA * p_N_x_data[0];
  elStorage->K26[0] = rhoIyy_tmp * N_data[0] * integrationFactor;
  elStorage->K26[2] = rhoIyy_tmp * N_data[1] * integrationFactor;
  rhoIyy_tmp = -GA * p_N_x_data[1];
  elStorage->K26[1] = rhoIyy_tmp * N_data[0] * integrationFactor;
  elStorage->K26[3] = rhoIyy_tmp * N_data[1] * integrationFactor;

  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  elStorage->K33[0] = K33_idx_0 + Jac;
  elStorage->K33[2] = K33_idx_2 + K12_tmp_tmp;
  elStorage->K33[1] = K33_idx_1 + K55_tmp;
  elStorage->K33[3] = K33_idx_3 + rhoIzz_tmp;

  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  elStorage->K35[0] = K14_tmp * N_data[0] * integrationFactor;
  elStorage->K35[2] = K14_tmp * N_data[1] * integrationFactor;
  elStorage->K35[1] = K11_tmp * N_data[0] * integrationFactor;
  elStorage->K35[3] = K11_tmp * N_data[1] * integrationFactor;

  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  K14_tmp = GA * N_data[0];
  Jac = K14_tmp * N_data[0] * integrationFactor;
  elStorage->K55[0] = K55_idx_0 + Jac;
  K14_tmp = K14_tmp * N_data[1] * integrationFactor;
  elStorage->K55[2] = K55_idx_2 + K14_tmp;
  K12_tmp_tmp = GA * N_data[1];
  K11_tmp = K12_tmp_tmp * N_data[0] * integrationFactor;
  elStorage->K55[1] = K55_idx_1 + K11_tmp;
  K12_tmp_tmp = K12_tmp_tmp * N_data[1] * integrationFactor;
  elStorage->K55[3] = K55_idx_3 + K12_tmp_tmp;

  // This function is a general routine to calculate an element matrix
  // Element calculation functions---------------------------------
  elStorage->K66[0] = K66_idx_0 + Jac;
  elStorage->K66[2] = K66_idx_2 + K14_tmp;
  elStorage->K66[1] = K66_idx_1 + K11_tmp;
  elStorage->K66[3] = K66_idx_3 + K12_tmp_tmp;

  // Store structural stiffness "K" into elementStorage
  // Store structural stiffness "M" into elementStorage
  // Store spin softening coefficient "S" into element storage
  // Store coriolis coefficient "C" into element sotrage
  elStorage->K11[0] = K11_idx_0;
  elStorage->K12[0] = K12_idx_0;
  elStorage->K13[0] = K13_idx_0;
  elStorage->K14[0] = K14_idx_0;
  elStorage->K15[0] = K15_idx_0;
  elStorage->K16[0] = K16_idx_0;
  elStorage->K24[0] = K24_idx_0;
  elStorage->K34[0] = K34_idx_0;
  elStorage->K44[0] = K44_idx_0;
  elStorage->K45[0] = K45_idx_0;
  elStorage->K46[0] = K46_idx_0;
  elStorage->K56[0] = K56_idx_0;
  elStorage->M11[0] = M11_idx_0;
  elStorage->M15[0] = M15_idx_0;
  elStorage->M16[0] = M16_idx_0;
  elStorage->M22[0] = M22_idx_0;
  elStorage->M24[0] = M24_idx_0;
  elStorage->M33[0] = M33_idx_0;
  elStorage->M34[0] = M34_idx_0;
  elStorage->M44[0] = M44_idx_0;
  elStorage->M55[0] = M55_idx_0;
  elStorage->M56[0] = M56_idx_0;
  elStorage->M66[0] = M66_idx_0;
  elStorage->S11[0] = 0.5 * S11_idx_0;
  elStorage->S12[0] = S12_idx_0;
  elStorage->S13[0] = S13_idx_0;
  elStorage->S15[0] = 0.5 * S15_idx_0;
  elStorage->S16[0] = 0.5 * S16_idx_0;
  elStorage->S22[0] = 0.5 * S22_idx_0;
  elStorage->S23[0] = S23_idx_0;
  elStorage->S25[0] = S25_idx_0;
  elStorage->S26[0] = S26_idx_0;
  elStorage->S33[0] = 0.5 * S33_idx_0;
  elStorage->S35[0] = S35_idx_0;
  elStorage->S36[0] = S36_idx_0;
  elStorage->S55[0] = 0.5 * S55_idx_0;
  elStorage->S56[0] = 0.5 * S56_idx_0;
  elStorage->S66[0] = 0.5 * S66_idx_0;
  elStorage->S14_1[0] = S14_1_idx_0;
  elStorage->S14_2[0] = S14_2_idx_0;
  elStorage->S24_1[0] = S24_1_idx_0;
  elStorage->S24_2[0] = S24_2_idx_0;
  elStorage->S34_1[0] = S34_1_idx_0;
  elStorage->S34_2[0] = S34_2_idx_0;
  elStorage->S45_1[0] = S45_1_idx_0;
  elStorage->S45_2[0] = S45_2_idx_0;
  elStorage->S46_1[0] = S46_1_idx_0;
  elStorage->S46_2[0] = S46_2_idx_0;
  elStorage->S44_1[0] = S44_1_idx_0;
  elStorage->S44_2[0] = S44_2_idx_0;
  elStorage->S44_3[0] = S44_3_idx_0;
  elStorage->C12[0] = C12_idx_0;
  elStorage->C13[0] = C13_idx_0;
  elStorage->C23[0] = C23_idx_0;
  elStorage->C24[0] = C24_idx_0;
  elStorage->C25[0] = C25_idx_0;
  elStorage->C26[0] = C26_idx_0;
  elStorage->C34[0] = C34_idx_0;
  elStorage->C35[0] = C35_idx_0;
  elStorage->C36[0] = C36_idx_0;
  elStorage->C14_1[0] = C14_1_idx_0;
  elStorage->C14_2[0] = C14_2_idx_0;
  elStorage->C45_1[0] = C45_1_idx_0;
  elStorage->C45_2[0] = C45_2_idx_0;
  elStorage->C46_1[0] = C46_1_idx_0;
  elStorage->C46_2[0] = C46_2_idx_0;
  elStorage->K11[1] = K11_idx_1;
  elStorage->K12[1] = K12_idx_1;
  elStorage->K13[1] = K13_idx_1;
  elStorage->K14[1] = K14_idx_1;
  elStorage->K15[1] = K15_idx_1;
  elStorage->K16[1] = K16_idx_1;
  elStorage->K24[1] = K24_idx_1;
  elStorage->K34[1] = K34_idx_1;
  elStorage->K44[1] = K44_idx_1;
  elStorage->K45[1] = K45_idx_1;
  elStorage->K46[1] = K46_idx_1;
  elStorage->K56[1] = K56_idx_1;
  elStorage->M11[1] = M11_idx_1;
  elStorage->M15[1] = M15_idx_1;
  elStorage->M16[1] = M16_idx_1;
  elStorage->M22[1] = M22_idx_1;
  elStorage->M24[1] = M24_idx_1;
  elStorage->M33[1] = M33_idx_1;
  elStorage->M34[1] = M34_idx_1;
  elStorage->M44[1] = M44_idx_1;
  elStorage->M55[1] = M55_idx_1;
  elStorage->M56[1] = M56_idx_1;
  elStorage->M66[1] = M66_idx_1;
  elStorage->S11[1] = 0.5 * S11_idx_1;
  elStorage->S12[1] = S12_idx_1;
  elStorage->S13[1] = S13_idx_1;
  elStorage->S15[1] = 0.5 * S15_idx_1;
  elStorage->S16[1] = 0.5 * S16_idx_1;
  elStorage->S22[1] = 0.5 * S22_idx_1;
  elStorage->S23[1] = S23_idx_1;
  elStorage->S25[1] = S25_idx_1;
  elStorage->S26[1] = S26_idx_1;
  elStorage->S33[1] = 0.5 * S33_idx_1;
  elStorage->S35[1] = S35_idx_1;
  elStorage->S36[1] = S36_idx_1;
  elStorage->S55[1] = 0.5 * S55_idx_1;
  elStorage->S56[1] = 0.5 * S56_idx_1;
  elStorage->S66[1] = 0.5 * S66_idx_1;
  elStorage->S14_1[1] = S14_1_idx_1;
  elStorage->S14_2[1] = S14_2_idx_1;
  elStorage->S24_1[1] = S24_1_idx_1;
  elStorage->S24_2[1] = S24_2_idx_1;
  elStorage->S34_1[1] = S34_1_idx_1;
  elStorage->S34_2[1] = S34_2_idx_1;
  elStorage->S45_1[1] = S45_1_idx_1;
  elStorage->S45_2[1] = S45_2_idx_1;
  elStorage->S46_1[1] = S46_1_idx_1;
  elStorage->S46_2[1] = S46_2_idx_1;
  elStorage->S44_1[1] = S44_1_idx_1;
  elStorage->S44_2[1] = S44_2_idx_1;
  elStorage->S44_3[1] = S44_3_idx_1;
  elStorage->C12[1] = C12_idx_1;
  elStorage->C13[1] = C13_idx_1;
  elStorage->C23[1] = C23_idx_1;
  elStorage->C24[1] = C24_idx_1;
  elStorage->C25[1] = C25_idx_1;
  elStorage->C26[1] = C26_idx_1;
  elStorage->C34[1] = C34_idx_1;
  elStorage->C35[1] = C35_idx_1;
  elStorage->C36[1] = C36_idx_1;
  elStorage->C14_1[1] = C14_1_idx_1;
  elStorage->C14_2[1] = C14_2_idx_1;
  elStorage->C45_1[1] = C45_1_idx_1;
  elStorage->C45_2[1] = C45_2_idx_1;
  elStorage->C46_1[1] = C46_1_idx_1;
  elStorage->C46_2[1] = C46_2_idx_1;
  elStorage->K11[2] = K11_idx_2;
  elStorage->K12[2] = K12_idx_2;
  elStorage->K13[2] = K13_idx_2;
  elStorage->K14[2] = K14_idx_2;
  elStorage->K15[2] = K15_idx_2;
  elStorage->K16[2] = K16_idx_2;
  elStorage->K24[2] = K24_idx_2;
  elStorage->K34[2] = K34_idx_2;
  elStorage->K44[2] = K44_idx_2;
  elStorage->K45[2] = K45_idx_2;
  elStorage->K46[2] = K46_idx_2;
  elStorage->K56[2] = K56_idx_2;
  elStorage->M11[2] = M11_idx_2;
  elStorage->M15[2] = M15_idx_2;
  elStorage->M16[2] = M16_idx_2;
  elStorage->M22[2] = M22_idx_2;
  elStorage->M24[2] = M24_idx_2;
  elStorage->M33[2] = M33_idx_2;
  elStorage->M34[2] = M34_idx_2;
  elStorage->M44[2] = M44_idx_2;
  elStorage->M55[2] = M55_idx_2;
  elStorage->M56[2] = M56_idx_2;
  elStorage->M66[2] = M66_idx_2;
  elStorage->S11[2] = 0.5 * S11_idx_2;
  elStorage->S12[2] = S12_idx_2;
  elStorage->S13[2] = S13_idx_2;
  elStorage->S15[2] = 0.5 * S15_idx_2;
  elStorage->S16[2] = 0.5 * S16_idx_2;
  elStorage->S22[2] = 0.5 * S22_idx_2;
  elStorage->S23[2] = S23_idx_2;
  elStorage->S25[2] = S25_idx_2;
  elStorage->S26[2] = S26_idx_2;
  elStorage->S33[2] = 0.5 * S33_idx_2;
  elStorage->S35[2] = S35_idx_2;
  elStorage->S36[2] = S36_idx_2;
  elStorage->S55[2] = 0.5 * S55_idx_2;
  elStorage->S56[2] = 0.5 * S56_idx_2;
  elStorage->S66[2] = 0.5 * S66_idx_2;
  elStorage->S14_1[2] = S14_1_idx_2;
  elStorage->S14_2[2] = S14_2_idx_2;
  elStorage->S24_1[2] = S24_1_idx_2;
  elStorage->S24_2[2] = S24_2_idx_2;
  elStorage->S34_1[2] = S34_1_idx_2;
  elStorage->S34_2[2] = S34_2_idx_2;
  elStorage->S45_1[2] = S45_1_idx_2;
  elStorage->S45_2[2] = S45_2_idx_2;
  elStorage->S46_1[2] = S46_1_idx_2;
  elStorage->S46_2[2] = S46_2_idx_2;
  elStorage->S44_1[2] = S44_1_idx_2;
  elStorage->S44_2[2] = S44_2_idx_2;
  elStorage->S44_3[2] = S44_3_idx_2;
  elStorage->C12[2] = C12_idx_2;
  elStorage->C13[2] = C13_idx_2;
  elStorage->C23[2] = C23_idx_2;
  elStorage->C24[2] = C24_idx_2;
  elStorage->C25[2] = C25_idx_2;
  elStorage->C26[2] = C26_idx_2;
  elStorage->C34[2] = C34_idx_2;
  elStorage->C35[2] = C35_idx_2;
  elStorage->C36[2] = C36_idx_2;
  elStorage->C14_1[2] = C14_1_idx_2;
  elStorage->C14_2[2] = C14_2_idx_2;
  elStorage->C45_1[2] = C45_1_idx_2;
  elStorage->C45_2[2] = C45_2_idx_2;
  elStorage->C46_1[2] = C46_1_idx_2;
  elStorage->C46_2[2] = C46_2_idx_2;
  elStorage->K11[3] = K11_idx_3;
  elStorage->K12[3] = K12_idx_3;
  elStorage->K13[3] = K13_idx_3;
  elStorage->K14[3] = K14_idx_3;
  elStorage->K15[3] = K15_idx_3;
  elStorage->K16[3] = K16_idx_3;
  elStorage->K24[3] = K24_idx_3;
  elStorage->K34[3] = K34_idx_3;
  elStorage->K44[3] = K44_idx_3;
  elStorage->K45[3] = K45_idx_3;
  elStorage->K46[3] = K46_idx_3;
  elStorage->K56[3] = K56_idx_3;
  elStorage->M11[3] = M11_idx_3;
  elStorage->M15[3] = M15_idx_3;
  elStorage->M16[3] = M16_idx_3;
  elStorage->M22[3] = M22_idx_3;
  elStorage->M24[3] = M24_idx_3;
  elStorage->M33[3] = M33_idx_3;
  elStorage->M34[3] = M34_idx_3;
  elStorage->M44[3] = M44_idx_3;
  elStorage->M55[3] = M55_idx_3;
  elStorage->M56[3] = M56_idx_3;
  elStorage->M66[3] = M66_idx_3;
  elStorage->S11[3] = 0.5 * S11_idx_3;
  elStorage->S12[3] = S12_idx_3;
  elStorage->S13[3] = S13_idx_3;
  elStorage->S15[3] = 0.5 * S15_idx_3;
  elStorage->S16[3] = 0.5 * S16_idx_3;
  elStorage->S22[3] = 0.5 * S22_idx_3;
  elStorage->S23[3] = S23_idx_3;
  elStorage->S25[3] = S25_idx_3;
  elStorage->S26[3] = S26_idx_3;
  elStorage->S33[3] = 0.5 * S33_idx_3;
  elStorage->S35[3] = S35_idx_3;
  elStorage->S36[3] = S36_idx_3;
  elStorage->S55[3] = 0.5 * S55_idx_3;
  elStorage->S56[3] = 0.5 * S56_idx_3;
  elStorage->S66[3] = 0.5 * S66_idx_3;
  elStorage->S14_1[3] = S14_1_idx_3;
  elStorage->S14_2[3] = S14_2_idx_3;
  elStorage->S24_1[3] = S24_1_idx_3;
  elStorage->S24_2[3] = S24_2_idx_3;
  elStorage->S34_1[3] = S34_1_idx_3;
  elStorage->S34_2[3] = S34_2_idx_3;
  elStorage->S45_1[3] = S45_1_idx_3;
  elStorage->S45_2[3] = S45_2_idx_3;
  elStorage->S46_1[3] = S46_1_idx_3;
  elStorage->S46_2[3] = S46_2_idx_3;
  elStorage->S44_1[3] = S44_1_idx_3;
  elStorage->S44_2[3] = S44_2_idx_3;
  elStorage->S44_3[3] = S44_3_idx_3;
  elStorage->C12[3] = C12_idx_3;
  elStorage->C13[3] = C13_idx_3;
  elStorage->C23[3] = C23_idx_3;
  elStorage->C24[3] = C24_idx_3;
  elStorage->C25[3] = C25_idx_3;
  elStorage->C26[3] = C26_idx_3;
  elStorage->C34[3] = C34_idx_3;
  elStorage->C35[3] = C35_idx_3;
  elStorage->C36[3] = C36_idx_3;
  elStorage->C14_1[3] = C14_1_idx_3;
  elStorage->C14_2[3] = C14_2_idx_3;
  elStorage->C45_1[3] = C45_1_idx_3;
  elStorage->C45_2[3] = C45_2_idx_3;
  elStorage->C46_1[3] = C46_1_idx_3;
  elStorage->C46_2[3] = C46_2_idx_3;
  for (b_i = 0; b_i < 3; b_i++) {
    lamSlimTran[3 * b_i] = lambda[b_i];
    lamSlimTran[3 * b_i + 1] = lambda[b_i + 12];
    lamSlimTran[3 * b_i + 2] = lambda[b_i + 24];
  }

  for (b_i = 0; b_i < 3; b_i++) {
    K14_tmp = lamSlimTran[b_i + 3];
    Jac = lamSlimTran[b_i + 6];
    for (i = 0; i < 3; i++) {
      b_lamSlimTran[b_i + 3 * i] = (lamSlimTran[b_i] * elementItens[3 * i] +
        K14_tmp * elementItens[3 * i + 1]) + Jac * elementItens[3 * i + 2];
    }
  }

  for (b_i = 0; b_i < 3; b_i++) {
    K14_tmp = b_lamSlimTran[b_i + 3];
    Jac = b_lamSlimTran[b_i + 6];
    K12_tmp_tmp = 0.0;
    for (i = 0; i < 3; i++) {
      elementItens_tmp = b_i + 3 * i;
      elementItens[elementItens_tmp] = (b_lamSlimTran[b_i] * lambda[12 * i] +
        K14_tmp * lambda[12 * i + 1]) + Jac * lambda[12 * i + 2];
      K12_tmp_tmp += lamSlimTran[elementItens_tmp] * elxm[i];
    }

    posLocal[b_i] = K12_tmp_tmp;
  }

  elxm[0] = posLocal[0];
  elxm[1] = posLocal[1];
  elxm[2] = posLocal[2];

  //
  //
  if (input_concMassFlag) {
    // modify element mass, moi, and xm to account for concentrated terms
    elementMass += input_concMass[0] + input_concMass[4];
    Jac = input_z[0] * input_z[0];
    rhoIyy_tmp = input_z[1] * input_z[1];
    rhoIzz_tmp = input_y[0] * input_y[0];
    K11_tmp = input_y[1] * input_y[1];
    elementItens[0] = (((elementItens[0] + input_concMass[0] * (rhoIzz_tmp + Jac))
                        + input_concMass[4] * (K11_tmp + rhoIyy_tmp)) +
                       input_concMass[1]) + input_concMass[5];
    K55_tmp = input_x[0] * input_x[0];
    K12_tmp_tmp = input_x[1] * input_x[1];
    elementItens[4] = (((elementItens[4] + input_concMass[0] * (K55_tmp + Jac))
                        + input_concMass[4] * (K12_tmp_tmp + rhoIyy_tmp)) +
                       input_concMass[2]) + input_concMass[6];
    elementItens[8] = (((elementItens[8] + input_concMass[0] * (K55_tmp +
      rhoIzz_tmp)) + input_concMass[4] * (K12_tmp_tmp + K11_tmp)) +
                       input_concMass[3]) + input_concMass[7];
    Jac = input_concMass[0] * input_x[0];
    rhoIyy_tmp = input_concMass[4] * input_x[1];
    rhoIzz_tmp = Jac * input_y[0];
    K11_tmp = rhoIyy_tmp * input_y[1];
    elementItens[3] = (elementItens[3] - rhoIzz_tmp) - K11_tmp;
    K55_tmp = Jac * input_z[0];
    K12_tmp_tmp = rhoIyy_tmp * input_z[1];
    elementItens[6] = (elementItens[6] - K55_tmp) - K12_tmp_tmp;
    elementItens[1] = (elementItens[1] - rhoIzz_tmp) - K11_tmp;
    GA = input_concMass[0] * input_y[0];
    rhoIzz_tmp = GA * input_z[0];
    K14_tmp = input_concMass[4] * input_y[1];
    K11_tmp = K14_tmp * input_z[1];
    elementItens[7] = (elementItens[7] - rhoIzz_tmp) - K11_tmp;
    elementItens[2] = (elementItens[2] - K55_tmp) - K12_tmp_tmp;
    elementItens[5] = (elementItens[5] - rhoIzz_tmp) - K11_tmp;
    elxm[0] = (posLocal[0] + Jac) + rhoIyy_tmp;
    elxm[1] = (posLocal[1] + GA) + K14_tmp;
    elxm[2] = (posLocal[2] + input_concMass[0] * input_z[0]) + input_concMass[4]
      * input_z[1];
  }

  // store element mass properties
  elStorage->mel = elementMass;
  std::memcpy(&elStorage->moiel[0], &elementItens[0], 9U * sizeof(double));
  elStorage->xmel[0] = elxm[0];
  elStorage->xmel[1] = elxm[1];
  elStorage->xmel[2] = elxm[2];
}

//
// File trailer for calculateTimoshenkoElementInitialRun.cpp
//
// [EOF]
//
