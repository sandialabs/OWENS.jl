//
// File: test_owens_types.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//
#ifndef TEST_OWENS_TYPES_H
#define TEST_OWENS_TYPES_H

// Include Files
#include "rtwtypes.h"

// Type Definitions
#include <stdio.h>
#include <time.h>

// Type Definitions
struct b_struct_T
{
  double delta_t;
  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
  double a6;
  double a7;
  double a8;
};

struct e_struct_T
{
  double ac[2];
  double twist[2];
  double rhoA[2];
  double EIyy[2];
  double EIzz[2];
  double GJ[2];
  double EA[2];
  double rhoIyy[2];
  double rhoIzz[2];
  double rhoJ[2];
  double zcm[2];
  double ycm[2];
  double a[2];
  double EIyz[2];
  double alpha1[2];
  double alpha2[2];
  double alpha3[2];
  double alpha4[2];
  double alpha5[2];
  double alpha6[2];
  double rhoIyz[2];
  double b[2];
  double a0[2];
  double aeroCenterOffset[2];
};

struct emxArray_real_T_1x12
{
  double data[12];
  int size[2];
};

struct emxArray_real_T_2
{
  double data[2];
  int size[1];
};

struct emxArray_real_T_2x2
{
  double data[4];
  int size[2];
};

struct s5oXUELIgdfYl3uJD5eSioC_tag
{
  char analysisType[3];
  double elementOrder;
  boolean_T modalFlag;
  b_struct_T timeInt;
  double xloc[2];
  e_struct_T sectionProps;
  double sweepAngle;
  double coneAngle;
  double rollAngle;
  double aeroSweepAngle;
  boolean_T firstIteration;
  double concMass[8];
  double concStiff[12];
  double concLoad[12];
  emxArray_real_T_1x12 disp;
  emxArray_real_T_1x12 dispm1;
  emxArray_real_T_1x12 dispdot;
  emxArray_real_T_1x12 dispddot;
  emxArray_real_T_2 x;
  emxArray_real_T_2 y;
  emxArray_real_T_2x2 z;
  double accelVec[3];
  double Omega;
  double OmegaDot;
  double omegaVec[3];
  double omegaDotVec[3];
  emxArray_real_T_1x12 displ_iter;
  boolean_T useDisp;
  boolean_T preStress;
  char iterationType[2];
  double freq;
  boolean_T aeroElasticOn;
  boolean_T aeroForceOn;
  double airDensity;
  boolean_T gravityOn;
  double RayleighAlpha;
  double RayleighBeta;
  double CN2H[9];
};

typedef s5oXUELIgdfYl3uJD5eSioC_tag o_struct_T;
struct emxArray_real_T_12x12
{
  double data[144];
  int size[2];
};

struct soak0T2LBhvoGwKUbBxyXmH_tag
{
  double FhatLessConc[12];
  double Ke[144];
  double Fe[12];
  emxArray_real_T_12x12 Me;
  emxArray_real_T_12x12 Ce;
};

typedef soak0T2LBhvoGwKUbBxyXmH_tag p_struct_T;
struct emxArray_char_T_1x3
{
  char data[3];
  int size[2];
};

struct emxArray_real_T_3x3
{
  double data[9];
  int size[2];
};

struct sOEAe0pL1P7QsnkyxNjo8cH_tag
{
  double elementOrder;
  boolean_T modalFlag;
  b_struct_T timeInt;
  double xloc[2];
  e_struct_T sectionProps;
  double sweepAngle;
  double coneAngle;
  double rollAngle;
  double aeroSweepAngle;
  char iterationType[2];
  boolean_T useDisp;
  boolean_T preStress;
  boolean_T aeroElasticOn;
  boolean_T aeroForceOn;
  double loadStepPrev;
  double loadStep;
  double maxNumLoadSteps;
  double MAXIT;
  double tolerance;
  emxArray_char_T_1x3 analysisType;
  emxArray_real_T_1x12 disp;
  emxArray_real_T_1x12 dispdot;
  emxArray_real_T_1x12 dispddot;
  emxArray_real_T_1x12 displ_iter;
  double concMass[8];
  double concStiff[12];
  double concLoad[12];
  emxArray_real_T_1x12 dispm1;
  emxArray_real_T_2 x;
  emxArray_real_T_2 y;
  emxArray_real_T_2x2 z;
  boolean_T gravityOn;
  double RayleighAlpha;
  double RayleighBeta;
  emxArray_real_T_3x3 accelVec;
  emxArray_real_T_3x3 omegaVec;
  emxArray_real_T_3x3 omegaDotVec;
  double Omega;
  double OmegaDot;
  double CN2H[9];
  double airDensity;
  double freq;
  boolean_T firstIteration;
};

typedef sOEAe0pL1P7QsnkyxNjo8cH_tag q_struct_T;
struct emxArray_real_T_60
{
  double data[60];
  int size[1];
};

struct emxArray_real_T_60x12
{
  double data[720];
  int size[2];
};

struct sjd27qDtc428idbHnJzvpKF_tag
{
  double numBlades;
  emxArray_real_T_60 bladeNum;
  emxArray_real_T_60 h;
  emxArray_real_T_60 nodeNum;
  emxArray_real_T_60 elementNum;
  emxArray_real_T_60x12 remaining;
};

typedef sjd27qDtc428idbHnJzvpKF_tag r_struct_T;
struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct sgIsDuxxBd5cAyrhqZQjvKB_tag
{
  double numpBC;
  emxArray_real_T *pBC;
  double numsBC;
  double nummBC;
  emxArray_real_T *isConstrained;
  emxArray_real_T *map;
  emxArray_real_T *redVectorMap;
};

typedef sgIsDuxxBd5cAyrhqZQjvKB_tag s_struct_T;
struct emxArray_struct_T
{
  e_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct h_struct_T
{
  emxArray_real_T *nodeNum;
  double numEl;
  double numNodes;
  emxArray_real_T *x;
  emxArray_real_T *y;
  emxArray_real_T *z;
  emxArray_real_T *elNum;
  emxArray_real_T *conn;
};

struct emxArray_char_T
{
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct emxArray_boolean_T
{
  boolean_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct i_struct_T
{
  emxArray_struct_T *props;
  emxArray_real_T *elLen;
  emxArray_real_T *psi;
  emxArray_real_T *theta;
  emxArray_real_T *roll;
  emxArray_boolean_T *rotationalEffects;
};

struct emxArray_uint32_T
{
  unsigned int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct emxArray_int8_T
{
  signed char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct f_struct_T
{
  double eps_xx_0[4];
  double eps_xx_z[4];
  double eps_xx_y[4];
  double gam_xz_0[4];
  double gam_xz_y[4];
  double gam_xy_0[4];
  double gam_xy_z[4];
};

struct b_emxArray_struct_T
{
  f_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct g_struct_T
{
  double K11[4];
  double K12[4];
  double K13[4];
  double K14[4];
  double K15[4];
  double K16[4];
  double K22[4];
  double K23[4];
  double K24[4];
  double K25[4];
  double K26[4];
  double K33[4];
  double K34[4];
  double K35[4];
  double K36[4];
  double K44[4];
  double K45[4];
  double K46[4];
  double K55[4];
  double K56[4];
  double K66[4];
  double M11[4];
  double M15[4];
  double M16[4];
  double M22[4];
  double M24[4];
  double M33[4];
  double M34[4];
  double M44[4];
  double M55[4];
  double M56[4];
  double M66[4];
  double S11[4];
  double S12[4];
  double S13[4];
  double S15[4];
  double S16[4];
  double S22[4];
  double S23[4];
  double S25[4];
  double S26[4];
  double S33[4];
  double S35[4];
  double S36[4];
  double S55[4];
  double S56[4];
  double S66[4];
  double S14_1[4];
  double S14_2[4];
  double S24_1[4];
  double S24_2[4];
  double S34_1[4];
  double S34_2[4];
  double S45_1[4];
  double S45_2[4];
  double S46_1[4];
  double S46_2[4];
  double S44_1[4];
  double S44_2[4];
  double S44_3[4];
  double C12[4];
  double C13[4];
  double C23[4];
  double C24[4];
  double C25[4];
  double C26[4];
  double C34[4];
  double C35[4];
  double C36[4];
  double C14_1[4];
  double C14_2[4];
  double C45_1[4];
  double C45_2[4];
  double C46_1[4];
  double C46_2[4];
  double mel;
  double moiel[9];
  double xmel[3];
};

struct c_emxArray_struct_T
{
  g_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct j_struct_T
{
  emxArray_real_T *displ_s;
  emxArray_real_T *displdot_s;
  emxArray_real_T *displddot_s;
};

struct k_struct_T
{
  b_emxArray_struct_T *elStrain;
  emxArray_real_T *displ_sp1;
  emxArray_real_T *displddot_sp1;
  emxArray_real_T *displdot_sp1;
};

struct l_struct_T
{
  double NElem;
  double FlipN;
  emxArray_real_T *QCx;
  emxArray_real_T *QCy;
  emxArray_real_T *QCz;
  emxArray_real_T *tx;
  emxArray_real_T *ty;
  emxArray_real_T *tz;
  emxArray_real_T *CtoR;
  emxArray_real_T *PEx;
  emxArray_real_T *PEy;
  emxArray_real_T *PEz;
  emxArray_real_T *tEx;
  emxArray_real_T *tEy;
  emxArray_real_T *tEz;
  emxArray_real_T *nEx;
  emxArray_real_T *nEy;
  emxArray_real_T *nEz;
  emxArray_real_T *sEx;
  emxArray_real_T *sEy;
  emxArray_real_T *sEz;
  emxArray_real_T *ECtoR;
  emxArray_real_T *EAreaR;
  emxArray_real_T *iSect;
};

struct d_emxArray_struct_T
{
  l_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct n_struct_T
{
  double NElem;
  double TtoC;
  emxArray_real_T *MCx;
  emxArray_real_T *MCy;
  emxArray_real_T *MCz;
  emxArray_real_T *CtoR;
  emxArray_real_T *PEx;
  emxArray_real_T *PEy;
  emxArray_real_T *PEz;
  emxArray_real_T *sEx;
  emxArray_real_T *sEy;
  emxArray_real_T *sEz;
  emxArray_real_T *ECtoR;
  emxArray_real_T *EAreaR;
  double BIndS;
  double EIndS;
  double BIndE;
  double EIndE;
};

struct f_emxArray_struct_T
{
  n_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct m_struct_T
{
  emxArray_real_T *N;
  emxArray_real_T *T;
  emxArray_real_T *M25;
};

struct e_emxArray_struct_T
{
  m_struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct coder_internal_sparse
{
  emxArray_real_T *d;
  emxArray_int32_T *colidx;
  emxArray_int32_T *rowidx;
};

struct emxArray_char_T_1x76
{
  char data[76];
  int size[2];
};

struct d_struct_T
{
  char iterationType[2];
  boolean_T adaptiveLoadSteppingFlag;
  double tolerance;
  double maxIterations;
  double maxNumLoadSteps;
  double minLoadStepDelta;
  double minLoadStep;
  double prescribedLoadStep;
};

struct emxArray_creal_T
{
  creal_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

struct stQEoBnY9NZL2taFdkfL4xB_tag
{
  char analysisType[2];
  double turbineStartup;
  boolean_T aeroElasticOn;
  boolean_T aeroForceOn;
  double airDensity;
  double guessFreq;
  boolean_T gravityOn;
  boolean_T generatorOn;
  double OmegaGenStart;
  boolean_T omegaControl;
  double totalNumDof;
  d_struct_T nlParams;
  boolean_T spinUpOn;
  boolean_T nlOn;
  double numModesToExtract;
  emxArray_char_T *aeroloadfile;
  emxArray_char_T *owensfile;
  double RayleighAlpha;
  double RayleighBeta;
  double elementOrder;
  r_struct_T bladeData;
  s_struct_T BC;
  emxArray_real_T *joint;
  boolean_T hydroOn;
  double c_platformTurbineConnectionNode;
  emxArray_char_T_1x76 outFilename;
  emxArray_real_T *jointTransform;
  emxArray_real_T *reducedDOFList;
};

typedef stQEoBnY9NZL2taFdkfL4xB_tag t_struct_T;

#endif

//
// File trailer for test_owens_types.h
//
// [EOF]
//
