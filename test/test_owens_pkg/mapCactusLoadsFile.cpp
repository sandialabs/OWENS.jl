//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: mapCactusLoadsFile.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 07-Apr-2020 17:21:12
//

// Include Files
#include "mapCactusLoadsFile.h"
#include "calculateLoadVecFromDistForce.h"
#include "fileManager.h"
#include "find.h"
#include "importCactusFile.h"
#include "linear_interp.h"
#include "readCactusGeom.h"
#include "readElementData.h"
#include "readMesh.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include "test_owens_rtwutil.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Declarations
static int div_s32_floor(int numerator, int denominator);
static void readBldFile(const emxArray_char_T *bldFn, double
  *bladeData_numBlades, double bladeData_bladeNum_data[], int
  bladeData_bladeNum_size[1], double bladeData_h_data[], int bladeData_h_size[1],
  double bladeData_nodeNum_data[], int bladeData_nodeNum_size[1], double
  bladeData_elementNum_data[], int bladeData_elementNum_size[1], double
  bladeData_remaining_data[], int bladeData_remaining_size[2], emxArray_real_T
  *structuralSpanLocNorm, emxArray_real_T *structuralNodeNumbers,
  emxArray_real_T *structuralElNumbers);

// Function Definitions

//
// Arguments    : int numerator
//                int denominator
// Return Type  : int
//
static int div_s32_floor(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  boolean_T quotientNeedsNegation;
  unsigned int tempAbsQuotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~static_cast<unsigned int>(numerator) + 1U;
    } else {
      absNumerator = static_cast<unsigned int>(numerator);
    }

    if (denominator < 0) {
      absDenominator = ~static_cast<unsigned int>(denominator) + 1U;
    } else {
      absDenominator = static_cast<unsigned int>(denominator);
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }

      quotient = -static_cast<int>(tempAbsQuotient);
    } else {
      quotient = static_cast<int>(tempAbsQuotient);
    }
  }

  return quotient;
}

//
// READ IN BLD FILE
// Arguments    : const emxArray_char_T *bldFn
//                double *bladeData_numBlades
//                double bladeData_bladeNum_data[]
//                int bladeData_bladeNum_size[1]
//                double bladeData_h_data[]
//                int bladeData_h_size[1]
//                double bladeData_nodeNum_data[]
//                int bladeData_nodeNum_size[1]
//                double bladeData_elementNum_data[]
//                int bladeData_elementNum_size[1]
//                double bladeData_remaining_data[]
//                int bladeData_remaining_size[2]
//                emxArray_real_T *structuralSpanLocNorm
//                emxArray_real_T *structuralNodeNumbers
//                emxArray_real_T *structuralElNumbers
// Return Type  : void
//
static void readBldFile(const emxArray_char_T *bldFn, double
  *bladeData_numBlades, double bladeData_bladeNum_data[], int
  bladeData_bladeNum_size[1], double bladeData_h_data[], int bladeData_h_size[1],
  double bladeData_nodeNum_data[], int bladeData_nodeNum_size[1], double
  bladeData_elementNum_data[], int bladeData_elementNum_size[1], double
  bladeData_remaining_data[], int bladeData_remaining_size[2], emxArray_real_T
  *structuralSpanLocNorm, emxArray_real_T *structuralNodeNumbers,
  emxArray_real_T *structuralElNumbers)
{
  double a_data[960];
  int a_size[2];
  int i;
  int idx;
  double numBlades;
  int k;
  boolean_T exitg1;
  int i1;
  double d;
  int b_i;
  int loop_ub;
  double numNodesPerBlade;
  int i2;
  double d1;
  double a;
  double b_a_data[60];
  importCactusFile(bldFn, a_data, a_size);
  i = a_size[0];
  if (a_size[0] <= 2) {
    if (a_size[0] == 1) {
      numBlades = a_data[0];
    } else if ((a_data[0] < a_data[1]) || (rtIsNaN(a_data[0]) && (!rtIsNaN
                 (a_data[1])))) {
      numBlades = a_data[1];
    } else {
      numBlades = a_data[0];
    }
  } else {
    if (!rtIsNaN(a_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= a_size[0])) {
        if (!rtIsNaN(a_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      numBlades = a_data[0];
    } else {
      numBlades = a_data[idx - 1];
      i1 = idx + 1;
      for (k = i1; k <= i; k++) {
        d = a_data[k - 1];
        if (numBlades < d) {
          numBlades = d;
        }
      }
    }
  }

  //  numStruts = min(bladeNum);
  //  if(numStruts>0)
  //  %     numStruts = 0;
  //  else
  //  %     numStruts = abs(numStruts);
  //  end
  idx = -1;
  b_i = 0;
  exitg1 = false;
  while ((!exitg1) && (b_i <= a_size[0] - 1)) {
    if (rtIsNaN(a_data[b_i + a_size[0] * 15])) {
      idx = b_i;
      exitg1 = true;
    } else {
      b_i++;
    }
  }

  if (idx + 1 == 0) {
    idx = a_size[0];
  } else {
    //      strutDataBlock = a(strutStartIndex:end,:);
    //      [strutEntries, ~] = size(strutDataBlock);
    //      numNodesPerStrut = strutEntries/numStruts;
    //      numElPerStrut = numNodesPerStrut - 1;
  }

  if (1 > idx) {
    loop_ub = -1;
  } else {
    loop_ub = idx - 1;
  }

  numNodesPerBlade = static_cast<double>((loop_ub + 1)) / numBlades;

  //  numElPerBlade = numNodesPerBlade - 1;
  //  [len,~]=size(bladeDataBlock);
  i = static_cast<int>(numBlades);
  i1 = structuralSpanLocNorm->size[0] * structuralSpanLocNorm->size[1];
  structuralSpanLocNorm->size[0] = i;
  i2 = static_cast<int>(numNodesPerBlade);
  structuralSpanLocNorm->size[1] = i2;
  emxEnsureCapacity_real_T(structuralSpanLocNorm, i1);
  idx = i * i2;
  for (i1 = 0; i1 < idx; i1++) {
    structuralSpanLocNorm->data[i1] = 0.0;
  }

  i1 = structuralNodeNumbers->size[0] * structuralNodeNumbers->size[1];
  structuralNodeNumbers->size[0] = i;
  structuralNodeNumbers->size[1] = i2;
  emxEnsureCapacity_real_T(structuralNodeNumbers, i1);
  for (i1 = 0; i1 < idx; i1++) {
    structuralNodeNumbers->data[i1] = 0.0;
  }

  i1 = structuralElNumbers->size[0] * structuralElNumbers->size[1];
  structuralElNumbers->size[0] = i;
  structuralElNumbers->size[1] = i2;
  emxEnsureCapacity_real_T(structuralElNumbers, i1);
  for (i1 = 0; i1 < idx; i1++) {
    structuralElNumbers->data[i1] = 0.0;
  }

  for (b_i = 0; b_i < i; b_i++) {
    d = ((static_cast<double>(b_i) + 1.0) - 1.0) * numNodesPerBlade + 1.0;
    d1 = (static_cast<double>(b_i) + 1.0) * numNodesPerBlade;
    if (d > d1) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(d) - 1;
      i2 = static_cast<int>(d1);
    }

    idx = static_cast<int>(((static_cast<double>(b_i) + 1.0) * numNodesPerBlade));
    a = a_data[(idx + a_size[0]) - 1];
    k = i2 - i1;
    for (i2 = 0; i2 < k; i2++) {
      b_a_data[i2] = a_data[(i1 + i2) + a_size[0]] / a;
    }

    k = structuralSpanLocNorm->size[1];
    for (i1 = 0; i1 < k; i1++) {
      structuralSpanLocNorm->data[b_i + structuralSpanLocNorm->size[0] * i1] =
        b_a_data[i1];
    }

    if (d > d1) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(d) - 1;
      i2 = idx;
    }

    k = i2 - i1;
    for (i2 = 0; i2 < k; i2++) {
      b_a_data[i2] = a_data[(i1 + i2) + a_size[0] * 2];
    }

    k = structuralNodeNumbers->size[1];
    for (i1 = 0; i1 < k; i1++) {
      structuralNodeNumbers->data[b_i + structuralNodeNumbers->size[0] * i1] =
        b_a_data[i1];
    }

    if (d > d1) {
      i1 = 0;
      idx = 0;
    } else {
      i1 = static_cast<int>(d) - 1;
    }

    k = idx - i1;
    for (i2 = 0; i2 < k; i2++) {
      b_a_data[i2] = a_data[(i1 + i2) + a_size[0] * 3];
    }

    k = structuralElNumbers->size[1];
    for (i1 = 0; i1 < k; i1++) {
      structuralElNumbers->data[b_i + structuralElNumbers->size[0] * i1] =
        b_a_data[i1];
    }
  }

  // assign data to bladeData object
  bladeData_bladeNum_size[0] = loop_ub + 1;
  bladeData_h_size[0] = loop_ub + 1;
  bladeData_nodeNum_size[0] = loop_ub + 1;
  bladeData_elementNum_size[0] = loop_ub + 1;
  if (0 <= loop_ub) {
    std::memcpy(&bladeData_bladeNum_data[0], &a_data[0], (loop_ub + 1) * sizeof
                (double));
  }

  for (i = 0; i <= loop_ub; i++) {
    bladeData_h_data[i] = a_data[i + a_size[0]];
    bladeData_nodeNum_data[i] = a_data[i + a_size[0] * 2];
    bladeData_elementNum_data[i] = a_data[i + a_size[0] * 3];
  }

  bladeData_remaining_size[0] = loop_ub + 1;
  bladeData_remaining_size[1] = 12;
  for (i = 0; i < 12; i++) {
    for (i1 = 0; i1 <= loop_ub; i1++) {
      bladeData_remaining_data[i1 + bladeData_remaining_size[0] * i] = a_data[i1
        + a_size[0] * (i + 4)];
    }
  }

  //
  *bladeData_numBlades = numBlades;
}

//
// function [time,aeroDistLoadsArrayTime,aeroDistLoadsNodeMap,aeroDistLoadsElMap] = mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)
// Arguments    : const emxArray_char_T *geomFn
//                const emxArray_char_T *loadsFn
//                const emxArray_char_T *bldFn
//                const emxArray_char_T *elFn
//                const emxArray_char_T *ortFn
//                const emxArray_char_T *meshFn
//                double time_data[]
//                int time_size[1]
//                emxArray_real_T *ForceValHist
//                emxArray_real_T *ForceDof
// Return Type  : void
//
void mapCactusLoadsFile(const emxArray_char_T *geomFn, const emxArray_char_T
  *loadsFn, const emxArray_char_T *bldFn, const emxArray_char_T *elFn, const
  emxArray_char_T *ortFn, const emxArray_char_T *meshFn, double time_data[], int
  time_size[1], emxArray_real_T *ForceValHist, emxArray_real_T *ForceDof)
{
  d_emxArray_struct_T *cactusGeom_blade;
  emxArray_real_T *data;
  emxArray_real_T *RefR;
  emxArray_real_T *expl_temp;
  emxArray_real_T *b_expl_temp;
  emxArray_real_T *c_expl_temp;
  f_emxArray_struct_T *d_expl_temp;
  int cactusGeom_NBlade;
  int m;
  int n;
  int i;
  double numAeroEl;
  int b_i;
  double numAeroTS;
  double node2;
  int i1;
  double uloc_data[2002];
  emxArray_int8_T *init_bladeForce_N;
  double NperSpan_data[2002];
  double TperSpan_data[2002];
  int i2;
  double M25perSpan_data[2002];
  emxArray_int8_T *init_bladeForce_T;
  emxArray_int8_T *init_bladeForce_M25;
  l_struct_T e_expl_temp;
  e_emxArray_struct_T *bladeForce;
  unsigned int b_index;
  emxArray_real_T *spanLocNorm;
  int j;
  int k;
  emxArray_real_T *structuralSpanLocNorm;
  emxArray_real_T *structuralNodeNumbers;
  emxArray_real_T *structuralElNumbers;
  double bladeData_bladeNum_data[60];
  int bladeData_bladeNum_size[1];
  double bladeData_h_data[60];
  int bladeData_h_size[1];
  double bladeData_nodeNum_data[60];
  int bladeData_nodeNum_size[1];
  double bladeData_elementNum_data[60];
  int bladeData_elementNum_size[1];
  double bladeData_remaining_data[720];
  int bladeData_remaining_size[2];
  emxArray_int8_T *init_structuralLoad_N;
  emxArray_int8_T *init_structuralLoad_T;
  emxArray_int8_T *init_structuralLoad_M25;
  l_struct_T f_expl_temp;
  e_emxArray_struct_T *structuralLoad;
  emxArray_real_T *b_spanLocNorm;
  emxArray_real_T *b_bladeForce;
  emxArray_real_T *b_structuralSpanLocNorm;
  emxArray_real_T *b_r;
  emxArray_real_T *mesh_x;
  emxArray_real_T *mesh_y;
  emxArray_real_T *mesh_z;
  emxArray_struct_T *el_props;
  emxArray_real_T *g_expl_temp;
  emxArray_boolean_T *h_expl_temp;
  boolean_T exitg1;
  emxArray_boolean_T *x;
  double b_node2;
  double dofList[12];
  double elInput_xloc[2];
  double elInput_sectionProps_twist[2];
  double elInput_extDistF2Node[2];
  double elInput_extDistF3Node[2];
  double elInput_extDistF4Node[2];
  double output_Fe[12];
  emxInit_struct_T3(&cactusGeom_blade, 2);
  emxInit_real_T(&data, 2);
  emxInit_real_T(&RefR, 1);
  emxInit_real_T(&expl_temp, 1);
  emxInit_real_T(&b_expl_temp, 1);
  emxInit_real_T(&c_expl_temp, 1);
  emxInit_struct_T5(&d_expl_temp, 2);
  readCactusGeom(geomFn, &cactusGeom_NBlade, &m, expl_temp, b_expl_temp,
                 c_expl_temp, RefR, cactusGeom_blade, d_expl_temp);
  b_importCactusFile(loadsFn, data);

  // define these from params file
  //      RefAR = cactusGeom.RefAR*ft2m*ft2m;
  n = RefR->size[0];
  emxFree_struct_T5(&d_expl_temp);
  for (i = 0; i < n; i++) {
    RefR->data[i] *= 0.30478512648582745;
  }

  // m/s
  numAeroEl = 0.0;
  for (b_i = 0; b_i < cactusGeom_NBlade; b_i++) {
    numAeroEl += cactusGeom_blade->data[b_i].NElem;
  }

  numAeroTS = static_cast<double>(data->size[0]) / numAeroEl;

  //  time = zeros(len/numAeroEl,1);
  node2 = rt_roundd_snf(numAeroEl);
  if (node2 < 2.147483648E+9) {
    if (node2 >= -2.147483648E+9) {
      i = static_cast<int>(node2);
    } else {
      i = MIN_int32_T;
    }
  } else if (node2 >= 2.147483648E+9) {
    i = MAX_int32_T;
  } else {
    i = 0;
  }

  if ((i == 0) || (((i > 0) && (1 > data->size[0])) || ((0 > i) && (data->size[0]
         > 1)))) {
    i = 1;
    i1 = -1;
  } else {
    i1 = data->size[0] - 1;
  }

  n = div_s32_floor(i1, i);
  time_size[0] = n + 1;
  for (i1 = 0; i1 <= n; i1++) {
    time_data[i1] = data->data[i * i1] * RefR->data[0] / 25.0;
  }

  n = data->size[0];
  for (i = 0; i < n; i++) {
    uloc_data[i] = data->data[i + data->size[0] * 11] * 25.0;
  }

  //      cl = data(:,13);
  //      cd = data(:,14);
  //
  //      cx = data(:,19);
  //      cy = data(:,20);
  //      cz = data(:,21);
  // calculate element areas
  //      Fx = zeros(len,1);
  //      Fy = Fx;
  //      Fz = Fx;
  i = data->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    numAeroEl = uloc_data[b_i] * uloc_data[b_i];
    m = static_cast<int>(data->data[b_i + data->size[0] * 2]) - 1;
    n = static_cast<int>(data->data[b_i + data->size[0] * 3]) - 1;
    node2 = cactusGeom_blade->data[m].ECtoR->data[n] * RefR->data[0];
    NperSpan_data[b_i] = data->data[b_i + data->size[0] * 16] * 0.5 * 1.225 *
      numAeroEl * node2;
    TperSpan_data[b_i] = data->data[b_i + data->size[0] * 17] * 0.5 * 1.225 *
      numAeroEl * node2;
    M25perSpan_data[b_i] = data->data[b_i + data->size[0] * 14] * 0.5 * 1.225 *
      numAeroEl * node2 * cactusGeom_blade->data[m].ECtoR->data[n] * RefR->data
      [0];
  }

  emxFree_real_T(&data);
  emxInit_int8_T(&init_bladeForce_N, 2);

  //  Initialize bladeForce
  i = static_cast<int>(numAeroTS);
  i1 = init_bladeForce_N->size[0] * init_bladeForce_N->size[1];
  init_bladeForce_N->size[0] = i;
  i2 = static_cast<int>(cactusGeom_blade->data[0].NElem);
  init_bladeForce_N->size[1] = i2;
  emxEnsureCapacity_int8_T(init_bladeForce_N, i1);
  m = i * i2;
  for (i1 = 0; i1 < m; i1++) {
    init_bladeForce_N->data[i1] = 0;
  }

  emxInit_int8_T(&init_bladeForce_T, 2);
  i1 = init_bladeForce_T->size[0] * init_bladeForce_T->size[1];
  init_bladeForce_T->size[0] = i;
  init_bladeForce_T->size[1] = i2;
  emxEnsureCapacity_int8_T(init_bladeForce_T, i1);
  for (i1 = 0; i1 < m; i1++) {
    init_bladeForce_T->data[i1] = 0;
  }

  emxInit_int8_T(&init_bladeForce_M25, 2);
  i1 = init_bladeForce_M25->size[0] * init_bladeForce_M25->size[1];
  init_bladeForce_M25->size[0] = i;
  init_bladeForce_M25->size[1] = i2;
  emxEnsureCapacity_int8_T(init_bladeForce_M25, i1);
  for (i1 = 0; i1 < m; i1++) {
    init_bladeForce_M25->data[i1] = 0;
  }

  emxInitStruct_struct_T4(&e_expl_temp);
  i1 = e_expl_temp.M25->size[0] * e_expl_temp.M25->size[1];
  e_expl_temp.M25->size[0] = init_bladeForce_M25->size[0];
  e_expl_temp.M25->size[1] = init_bladeForce_M25->size[1];
  emxEnsureCapacity_real_T(e_expl_temp.M25, i1);
  n = init_bladeForce_M25->size[0] * init_bladeForce_M25->size[1];
  for (i1 = 0; i1 < n; i1++) {
    e_expl_temp.M25->data[i1] = init_bladeForce_M25->data[i1];
  }

  emxFree_int8_T(&init_bladeForce_M25);
  i1 = e_expl_temp.T->size[0] * e_expl_temp.T->size[1];
  e_expl_temp.T->size[0] = init_bladeForce_T->size[0];
  e_expl_temp.T->size[1] = init_bladeForce_T->size[1];
  emxEnsureCapacity_real_T(e_expl_temp.T, i1);
  n = init_bladeForce_T->size[0] * init_bladeForce_T->size[1];
  for (i1 = 0; i1 < n; i1++) {
    e_expl_temp.T->data[i1] = init_bladeForce_T->data[i1];
  }

  emxFree_int8_T(&init_bladeForce_T);
  i1 = e_expl_temp.N->size[0] * e_expl_temp.N->size[1];
  e_expl_temp.N->size[0] = init_bladeForce_N->size[0];
  e_expl_temp.N->size[1] = init_bladeForce_N->size[1];
  emxEnsureCapacity_real_T(e_expl_temp.N, i1);
  n = init_bladeForce_N->size[0] * init_bladeForce_N->size[1];
  for (i1 = 0; i1 < n; i1++) {
    e_expl_temp.N->data[i1] = init_bladeForce_N->data[i1];
  }

  emxFree_int8_T(&init_bladeForce_N);
  emxInit_struct_T4(&bladeForce, 2);
  c_repmat(e_expl_temp, cactusGeom_NBlade, bladeForce);
  b_index = 1U;
  emxFreeStruct_struct_T5(&e_expl_temp);
  for (b_i = 0; b_i < i; b_i++) {
    for (j = 0; j < cactusGeom_NBlade; j++) {
      i1 = static_cast<int>(cactusGeom_blade->data[j].NElem);
      for (k = 0; k < i1; k++) {
        //             %%
        //                  bladeForce(j).Fx(i,k) = Fx(index);
        //                  bladeForce(j).Fy(i,k) = Fy(index);
        //                  bladeForce(j).Fz(i,k) = Fz(index);
        //                  bladeForce(j).M25(i,k) = M25perSpan(index);
        //                  spanVec = [blade(j).sEx(k);blade(j).sEy(k);blade(j).sEz(k)]; 
        //                  tanVec = [blade(j).tEx(k);blade(j).tEy(k);blade(j).tEz(k)]; 
        //                  normVec = [blade(j).nEx(k);blade(j).nEy(k);blade(j).nEz(k)]; 
        //  %                 dcm = [spanVec, tanVec, normVec];
        //                  dcm = [tanVec, spanVec, -normVec];
        //             %%
        n = static_cast<int>(b_index) - 1;
        bladeForce->data[j].N->data[b_i + bladeForce->data[j].N->size[0] * k] =
          NperSpan_data[n];
        bladeForce->data[j].T->data[b_i + bladeForce->data[j].T->size[0] * k] =
          TperSpan_data[n];
        bladeForce->data[j].M25->data[b_i + bladeForce->data[j].M25->size[0] * k]
          = M25perSpan_data[n];
        b_index = static_cast<unsigned int>((static_cast<int>(b_index) + 1));
      }
    }
  }

  emxInit_real_T(&spanLocNorm, 2);
  i1 = spanLocNorm->size[0] * spanLocNorm->size[1];
  spanLocNorm->size[0] = cactusGeom_NBlade;
  spanLocNorm->size[1] = i2;
  emxEnsureCapacity_real_T(spanLocNorm, i1);
  n = cactusGeom_NBlade * i2;
  for (i1 = 0; i1 < n; i1++) {
    spanLocNorm->data[i1] = 0.0;
  }

  for (b_i = 0; b_i < cactusGeom_NBlade; b_i++) {
    numAeroEl = cactusGeom_blade->data[b_i].QCy->data[static_cast<int>
      ((cactusGeom_blade->data[0].NElem + 1.0)) - 1] * RefR->data[0];
    n = spanLocNorm->size[1];
    for (i1 = 0; i1 < n; i1++) {
      spanLocNorm->data[b_i + spanLocNorm->size[0] * i1] =
        cactusGeom_blade->data[b_i].PEy->data[i1] * RefR->data[0] / numAeroEl;
    }
  }

  emxFree_struct_T3(&cactusGeom_blade);
  emxInit_real_T(&structuralSpanLocNorm, 2);
  emxInit_real_T(&structuralNodeNumbers, 2);
  emxInit_real_T(&structuralElNumbers, 2);
  readBldFile(bldFn, &numAeroEl, bladeData_bladeNum_data,
              bladeData_bladeNum_size, bladeData_h_data, bladeData_h_size,
              bladeData_nodeNum_data, bladeData_nodeNum_size,
              bladeData_elementNum_data, bladeData_elementNum_size,
              bladeData_remaining_data, bladeData_remaining_size,
              structuralSpanLocNorm, structuralNodeNumbers, structuralElNumbers);

  // Initialize structuralLoad
  if ((structuralElNumbers->size[0] == 0) || (structuralElNumbers->size[1] == 0))
  {
    n = 0;
  } else {
    m = structuralElNumbers->size[0];
    n = structuralElNumbers->size[1];
    if (m > n) {
      n = m;
    }
  }

  emxInit_int8_T(&init_structuralLoad_N, 2);
  i1 = init_structuralLoad_N->size[0] * init_structuralLoad_N->size[1];
  init_structuralLoad_N->size[0] = i;
  init_structuralLoad_N->size[1] = n;
  emxEnsureCapacity_int8_T(init_structuralLoad_N, i1);
  n *= i;
  for (i1 = 0; i1 < n; i1++) {
    init_structuralLoad_N->data[i1] = 0;
  }

  if ((structuralElNumbers->size[0] == 0) || (structuralElNumbers->size[1] == 0))
  {
    n = 0;
  } else {
    m = structuralElNumbers->size[0];
    n = structuralElNumbers->size[1];
    if (m > n) {
      n = m;
    }
  }

  emxInit_int8_T(&init_structuralLoad_T, 2);
  i1 = init_structuralLoad_T->size[0] * init_structuralLoad_T->size[1];
  init_structuralLoad_T->size[0] = i;
  init_structuralLoad_T->size[1] = n;
  emxEnsureCapacity_int8_T(init_structuralLoad_T, i1);
  n *= i;
  for (i1 = 0; i1 < n; i1++) {
    init_structuralLoad_T->data[i1] = 0;
  }

  if ((structuralElNumbers->size[0] == 0) || (structuralElNumbers->size[1] == 0))
  {
    n = 0;
  } else {
    m = structuralElNumbers->size[0];
    n = structuralElNumbers->size[1];
    if (m > n) {
      n = m;
    }
  }

  emxInit_int8_T(&init_structuralLoad_M25, 2);
  i1 = init_structuralLoad_M25->size[0] * init_structuralLoad_M25->size[1];
  init_structuralLoad_M25->size[0] = i;
  init_structuralLoad_M25->size[1] = n;
  emxEnsureCapacity_int8_T(init_structuralLoad_M25, i1);
  n *= i;
  for (i1 = 0; i1 < n; i1++) {
    init_structuralLoad_M25->data[i1] = 0;
  }

  emxInitStruct_struct_T4(&f_expl_temp);
  i1 = f_expl_temp.M25->size[0] * f_expl_temp.M25->size[1];
  f_expl_temp.M25->size[0] = init_structuralLoad_M25->size[0];
  f_expl_temp.M25->size[1] = init_structuralLoad_M25->size[1];
  emxEnsureCapacity_real_T(f_expl_temp.M25, i1);
  n = init_structuralLoad_M25->size[0] * init_structuralLoad_M25->size[1];
  for (i1 = 0; i1 < n; i1++) {
    f_expl_temp.M25->data[i1] = init_structuralLoad_M25->data[i1];
  }

  emxFree_int8_T(&init_structuralLoad_M25);
  i1 = f_expl_temp.T->size[0] * f_expl_temp.T->size[1];
  f_expl_temp.T->size[0] = init_structuralLoad_T->size[0];
  f_expl_temp.T->size[1] = init_structuralLoad_T->size[1];
  emxEnsureCapacity_real_T(f_expl_temp.T, i1);
  n = init_structuralLoad_T->size[0] * init_structuralLoad_T->size[1];
  for (i1 = 0; i1 < n; i1++) {
    f_expl_temp.T->data[i1] = init_structuralLoad_T->data[i1];
  }

  emxFree_int8_T(&init_structuralLoad_T);
  i1 = f_expl_temp.N->size[0] * f_expl_temp.N->size[1];
  f_expl_temp.N->size[0] = init_structuralLoad_N->size[0];
  f_expl_temp.N->size[1] = init_structuralLoad_N->size[1];
  emxEnsureCapacity_real_T(f_expl_temp.N, i1);
  n = init_structuralLoad_N->size[0] * init_structuralLoad_N->size[1];
  for (i1 = 0; i1 < n; i1++) {
    f_expl_temp.N->data[i1] = init_structuralLoad_N->data[i1];
  }

  emxFree_int8_T(&init_structuralLoad_N);
  emxInit_struct_T4(&structuralLoad, 2);
  c_repmat(f_expl_temp, cactusGeom_NBlade, structuralLoad);
  emxFreeStruct_struct_T5(&f_expl_temp);
  emxInit_real_T(&b_spanLocNorm, 2);
  emxInit_real_T(&b_bladeForce, 2);
  emxInit_real_T(&b_structuralSpanLocNorm, 2);
  emxInit_real_T(&b_r, 2);
  for (b_i = 0; b_i < cactusGeom_NBlade; b_i++) {
    for (j = 0; j < i; j++) {
      n = spanLocNorm->size[1];
      i1 = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
      b_spanLocNorm->size[0] = 1;
      b_spanLocNorm->size[1] = spanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_spanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_spanLocNorm->data[i1] = spanLocNorm->data[b_i + spanLocNorm->size[0] *
          i1];
      }

      n = bladeForce->data[b_i].N->size[1];
      i1 = b_bladeForce->size[0] * b_bladeForce->size[1];
      b_bladeForce->size[0] = 1;
      b_bladeForce->size[1] = bladeForce->data[b_i].N->size[1];
      emxEnsureCapacity_real_T(b_bladeForce, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_bladeForce->data[i1] = bladeForce->data[b_i].N->data[j +
          bladeForce->data[b_i].N->size[0] * i1];
      }

      n = structuralSpanLocNorm->size[1];
      i1 = b_structuralSpanLocNorm->size[0] * b_structuralSpanLocNorm->size[1];
      b_structuralSpanLocNorm->size[0] = 1;
      b_structuralSpanLocNorm->size[1] = structuralSpanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_structuralSpanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_structuralSpanLocNorm->data[i1] = structuralSpanLocNorm->data[b_i +
          structuralSpanLocNorm->size[0] * i1];
      }

      linear_interp(b_spanLocNorm, b_bladeForce, b_structuralSpanLocNorm, b_r);
      n = b_r->size[1];
      for (i1 = 0; i1 < n; i1++) {
        structuralLoad->data[b_i].N->data[j + structuralLoad->data[b_i].N->size
          [0] * i1] = b_r->data[i1];
      }

      n = spanLocNorm->size[1];
      i1 = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
      b_spanLocNorm->size[0] = 1;
      b_spanLocNorm->size[1] = spanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_spanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_spanLocNorm->data[i1] = spanLocNorm->data[b_i + spanLocNorm->size[0] *
          i1];
      }

      n = bladeForce->data[b_i].T->size[1];
      i1 = b_bladeForce->size[0] * b_bladeForce->size[1];
      b_bladeForce->size[0] = 1;
      b_bladeForce->size[1] = bladeForce->data[b_i].T->size[1];
      emxEnsureCapacity_real_T(b_bladeForce, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_bladeForce->data[i1] = bladeForce->data[b_i].T->data[j +
          bladeForce->data[b_i].T->size[0] * i1];
      }

      n = structuralSpanLocNorm->size[1];
      i1 = b_structuralSpanLocNorm->size[0] * b_structuralSpanLocNorm->size[1];
      b_structuralSpanLocNorm->size[0] = 1;
      b_structuralSpanLocNorm->size[1] = structuralSpanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_structuralSpanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_structuralSpanLocNorm->data[i1] = structuralSpanLocNorm->data[b_i +
          structuralSpanLocNorm->size[0] * i1];
      }

      linear_interp(b_spanLocNorm, b_bladeForce, b_structuralSpanLocNorm, b_r);
      n = b_r->size[1];
      for (i1 = 0; i1 < n; i1++) {
        structuralLoad->data[b_i].T->data[j + structuralLoad->data[b_i].T->size
          [0] * i1] = b_r->data[i1];
      }

      n = spanLocNorm->size[1];
      i1 = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
      b_spanLocNorm->size[0] = 1;
      b_spanLocNorm->size[1] = spanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_spanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_spanLocNorm->data[i1] = spanLocNorm->data[b_i + spanLocNorm->size[0] *
          i1];
      }

      n = bladeForce->data[b_i].M25->size[1];
      i1 = b_bladeForce->size[0] * b_bladeForce->size[1];
      b_bladeForce->size[0] = 1;
      b_bladeForce->size[1] = bladeForce->data[b_i].M25->size[1];
      emxEnsureCapacity_real_T(b_bladeForce, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_bladeForce->data[i1] = bladeForce->data[b_i].M25->data[j +
          bladeForce->data[b_i].M25->size[0] * i1];
      }

      n = structuralSpanLocNorm->size[1];
      i1 = b_structuralSpanLocNorm->size[0] * b_structuralSpanLocNorm->size[1];
      b_structuralSpanLocNorm->size[0] = 1;
      b_structuralSpanLocNorm->size[1] = structuralSpanLocNorm->size[1];
      emxEnsureCapacity_real_T(b_structuralSpanLocNorm, i1);
      for (i1 = 0; i1 < n; i1++) {
        b_structuralSpanLocNorm->data[i1] = structuralSpanLocNorm->data[b_i +
          structuralSpanLocNorm->size[0] * i1];
      }

      linear_interp(b_spanLocNorm, b_bladeForce, b_structuralSpanLocNorm, b_r);
      n = b_r->size[1];
      for (i1 = 0; i1 < n; i1++) {
        structuralLoad->data[b_i].M25->data[j + structuralLoad->data[b_i]
          .M25->size[0] * i1] = b_r->data[i1];
      }
    }
  }

  emxFree_real_T(&b_r);
  emxFree_real_T(&b_structuralSpanLocNorm);
  emxFree_real_T(&b_bladeForce);
  emxFree_real_T(&structuralSpanLocNorm);
  emxFree_struct_T4(&bladeForce);
  emxInit_real_T(&mesh_x, 1);
  emxInit_real_T(&mesh_y, 1);
  emxInit_real_T(&mesh_z, 1);
  emxInit_struct_T(&el_props, 2);
  emxInit_real_T(&g_expl_temp, 2);
  emxInit_boolean_T(&h_expl_temp, 2);

  // integrate over elements
  // read element data in
  b_readMesh(meshFn, expl_temp, &numAeroEl, &node2, mesh_x, mesh_y, mesh_z,
             b_expl_temp, g_expl_temp);
  b_readElementData(numAeroEl, elFn, ortFn, bladeData_nodeNum_data,
                    bladeData_nodeNum_size, bladeData_elementNum_data,
                    bladeData_elementNum_size, bladeData_remaining_data,
                    bladeData_remaining_size, el_props, expl_temp, c_expl_temp,
                    RefR, b_expl_temp, h_expl_temp);

  //      [~,~,timeLen] = size(aeroDistLoadsArrayTime);
  m = structuralNodeNumbers->size[0];
  n = structuralNodeNumbers->size[1];
  i1 = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
  b_spanLocNorm->size[0] = 1;
  b_spanLocNorm->size[1] = structuralNodeNumbers->size[1];
  emxEnsureCapacity_real_T(b_spanLocNorm, i1);
  emxFree_boolean_T(&h_expl_temp);
  emxFree_real_T(&g_expl_temp);
  emxFree_real_T(&expl_temp);
  if (structuralNodeNumbers->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      b_spanLocNorm->data[j] = structuralNodeNumbers->data
        [structuralNodeNumbers->size[0] * j];
      for (b_i = 2; b_i <= m; b_i++) {
        numAeroTS = b_spanLocNorm->data[j];
        numAeroEl = structuralNodeNumbers->data[(b_i +
          structuralNodeNumbers->size[0] * j) - 1];
        if ((!rtIsNaN(numAeroEl)) && (rtIsNaN(numAeroTS) || (numAeroTS <
              numAeroEl))) {
          b_spanLocNorm->data[j] = numAeroEl;
        }
      }
    }
  }

  n = b_spanLocNorm->size[1];
  if (b_spanLocNorm->size[1] <= 2) {
    if (b_spanLocNorm->size[1] == 1) {
      numAeroEl = b_spanLocNorm->data[0];
    } else if ((b_spanLocNorm->data[0] < b_spanLocNorm->data[1]) || (rtIsNaN
                (b_spanLocNorm->data[0]) && (!rtIsNaN(b_spanLocNorm->data[1]))))
    {
      numAeroEl = b_spanLocNorm->data[1];
    } else {
      numAeroEl = b_spanLocNorm->data[0];
    }
  } else {
    if (!rtIsNaN(b_spanLocNorm->data[0])) {
      m = 1;
    } else {
      m = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= b_spanLocNorm->size[1])) {
        if (!rtIsNaN(b_spanLocNorm->data[k - 1])) {
          m = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (m == 0) {
      numAeroEl = b_spanLocNorm->data[0];
    } else {
      numAeroEl = b_spanLocNorm->data[m - 1];
      i1 = m + 1;
      for (k = i1; k <= n; k++) {
        node2 = b_spanLocNorm->data[k - 1];
        if (numAeroEl < node2) {
          numAeroEl = node2;
        }
      }
    }
  }

  i1 = static_cast<int>((numAeroEl * 6.0));
  i2 = spanLocNorm->size[0] * spanLocNorm->size[1];
  spanLocNorm->size[0] = i1;
  spanLocNorm->size[1] = i;
  emxEnsureCapacity_real_T(spanLocNorm, i2);
  n = i1 * i;
  for (i1 = 0; i1 < n; i1++) {
    spanLocNorm->data[i1] = 0.0;
  }

  for (b_i = 0; b_i < i; b_i++) {
    for (j = 0; j < cactusGeom_NBlade; j++) {
      i1 = structuralNodeNumbers->size[1];
      for (k = 0; k <= i1 - 2; k++) {
        // get element data
        //  orientation angle,xloc,sectionProps,element order]
        // get dof map
        numAeroEl = structuralNodeNumbers->data[j + structuralNodeNumbers->size
          [0] * k];
        node2 = structuralNodeNumbers->data[j + structuralNodeNumbers->size[0] *
          (k + 1)];
        numAeroTS = (numAeroEl - 1.0) * 6.0;
        b_node2 = (node2 - 1.0) * 6.0;
        for (i2 = 0; i2 < 6; i2++) {
          dofList[i2] = numAeroTS + (static_cast<double>(i2) + 1.0);
          dofList[i2 + 6] = b_node2 + (static_cast<double>(i2) + 1.0);
        }

        n = static_cast<int>(node2) - 1;
        m = static_cast<int>(numAeroEl) - 1;
        numAeroTS = mesh_x->data[n] - mesh_x->data[m];
        numAeroEl = mesh_y->data[n] - mesh_y->data[m];
        node2 = mesh_z->data[n] - mesh_z->data[m];
        elInput_xloc[0] = 0.0;
        elInput_xloc[1] = std::sqrt((numAeroTS * numAeroTS + numAeroEl *
          numAeroEl) + node2 * node2);
        n = static_cast<int>(structuralElNumbers->data[j +
                             structuralElNumbers->size[0] * k]) - 1;
        elInput_sectionProps_twist[0] = el_props->data[n].twist[0];
        elInput_sectionProps_twist[1] = el_props->data[n].twist[1];
        elInput_extDistF2Node[0] = structuralLoad->data[j].T->data[b_i +
          structuralLoad->data[j].T->size[0] * k];
        elInput_extDistF2Node[1] = structuralLoad->data[j].T->data[b_i +
          structuralLoad->data[j].T->size[0] * (k + 1)];
        elInput_extDistF3Node[0] = -structuralLoad->data[j].N->data[b_i +
          structuralLoad->data[j].N->size[0] * k];
        elInput_extDistF3Node[1] = -structuralLoad->data[j].N->data[b_i +
          structuralLoad->data[j].N->size[0] * (k + 1)];
        elInput_extDistF4Node[0] = -structuralLoad->data[j].M25->data[b_i +
          structuralLoad->data[j].M25->size[0] * k];
        elInput_extDistF4Node[1] = -structuralLoad->data[j].M25->data[b_i +
          structuralLoad->data[j].M25->size[0] * (k + 1)];
        calculateLoadVecFromDistForce(elInput_xloc, elInput_sectionProps_twist,
          c_expl_temp->data[n], RefR->data[n], b_expl_temp->data[n],
          elInput_extDistF2Node, elInput_extDistF3Node, elInput_extDistF4Node,
          output_Fe);

        // asssembly
        for (m = 0; m < 12; m++) {
          n = static_cast<int>(dofList[m]) - 1;
          spanLocNorm->data[n + spanLocNorm->size[0] * b_i] += output_Fe[m];
        }
      }
    }
  }

  emxFree_real_T(&c_expl_temp);
  emxFree_real_T(&b_expl_temp);
  emxFree_real_T(&structuralElNumbers);
  emxFree_struct_T(&el_props);
  emxFree_real_T(&mesh_z);
  emxFree_real_T(&mesh_y);
  emxFree_real_T(&mesh_x);
  emxFree_struct_T4(&structuralLoad);
  emxFree_real_T(&RefR);
  emxInit_boolean_T(&x, 1);

  // reduce Fg to nonzero components
  // assumes any loaded DOF will never be identically zero throughout time
  // history
  n = spanLocNorm->size[0];
  i = x->size[0];
  x->size[0] = spanLocNorm->size[0];
  emxEnsureCapacity_boolean_T(x, i);
  for (i = 0; i < n; i++) {
    x->data[i] = (spanLocNorm->data[i] != 0.0);
  }

  m = x->size[0];
  if (x->size[0] == 0) {
    n = 0;
  } else {
    n = x->data[0];
    for (k = 2; k <= m; k++) {
      n += x->data[k - 1];
    }
  }

  i = ForceValHist->size[0] * ForceValHist->size[1];
  ForceValHist->size[0] = n;
  ForceValHist->size[1] = spanLocNorm->size[1];
  emxEnsureCapacity_real_T(ForceValHist, i);
  n *= spanLocNorm->size[1];
  for (i = 0; i < n; i++) {
    ForceValHist->data[i] = 0.0;
  }

  n = spanLocNorm->size[0];
  i = x->size[0];
  x->size[0] = spanLocNorm->size[0];
  emxEnsureCapacity_boolean_T(x, i);
  for (i = 0; i < n; i++) {
    x->data[i] = (spanLocNorm->data[i] != 0.0);
  }

  m = x->size[0];
  if (x->size[0] == 0) {
    n = 0;
  } else {
    n = x->data[0];
    for (k = 2; k <= m; k++) {
      n += x->data[k - 1];
    }
  }

  emxFree_boolean_T(&x);
  i = ForceDof->size[0];
  ForceDof->size[0] = n;
  emxEnsureCapacity_real_T(ForceDof, i);
  for (i = 0; i < n; i++) {
    ForceDof->data[i] = 0.0;
  }

  b_index = 1U;
  m = structuralNodeNumbers->size[0];
  n = structuralNodeNumbers->size[1];
  i = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
  b_spanLocNorm->size[0] = 1;
  b_spanLocNorm->size[1] = structuralNodeNumbers->size[1];
  emxEnsureCapacity_real_T(b_spanLocNorm, i);
  if (structuralNodeNumbers->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      b_spanLocNorm->data[j] = structuralNodeNumbers->data
        [structuralNodeNumbers->size[0] * j];
      for (b_i = 2; b_i <= m; b_i++) {
        numAeroTS = b_spanLocNorm->data[j];
        numAeroEl = structuralNodeNumbers->data[(b_i +
          structuralNodeNumbers->size[0] * j) - 1];
        if ((!rtIsNaN(numAeroEl)) && (rtIsNaN(numAeroTS) || (numAeroTS <
              numAeroEl))) {
          b_spanLocNorm->data[j] = numAeroEl;
        }
      }
    }
  }

  emxFree_real_T(&structuralNodeNumbers);
  n = b_spanLocNorm->size[1];
  if (b_spanLocNorm->size[1] <= 2) {
    if (b_spanLocNorm->size[1] == 1) {
      numAeroEl = b_spanLocNorm->data[0];
    } else if ((b_spanLocNorm->data[0] < b_spanLocNorm->data[1]) || (rtIsNaN
                (b_spanLocNorm->data[0]) && (!rtIsNaN(b_spanLocNorm->data[1]))))
    {
      numAeroEl = b_spanLocNorm->data[1];
    } else {
      numAeroEl = b_spanLocNorm->data[0];
    }
  } else {
    if (!rtIsNaN(b_spanLocNorm->data[0])) {
      m = 1;
    } else {
      m = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= b_spanLocNorm->size[1])) {
        if (!rtIsNaN(b_spanLocNorm->data[k - 1])) {
          m = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (m == 0) {
      numAeroEl = b_spanLocNorm->data[0];
    } else {
      numAeroEl = b_spanLocNorm->data[m - 1];
      i = m + 1;
      for (k = i; k <= n; k++) {
        node2 = b_spanLocNorm->data[k - 1];
        if (numAeroEl < node2) {
          numAeroEl = node2;
        }
      }
    }
  }

  i = static_cast<int>((numAeroEl * 6.0));
  for (b_i = 0; b_i < i; b_i++) {
    n = spanLocNorm->size[1];
    i1 = b_spanLocNorm->size[0] * b_spanLocNorm->size[1];
    b_spanLocNorm->size[0] = 1;
    b_spanLocNorm->size[1] = spanLocNorm->size[1];
    emxEnsureCapacity_real_T(b_spanLocNorm, i1);
    for (i1 = 0; i1 < n; i1++) {
      b_spanLocNorm->data[i1] = spanLocNorm->data[b_i + spanLocNorm->size[0] *
        i1];
    }

    c_eml_find(b_spanLocNorm, bladeData_bladeNum_size, bladeData_remaining_size);
    if (bladeData_remaining_size[1] != 0) {
      n = spanLocNorm->size[1];
      for (i1 = 0; i1 < n; i1++) {
        ForceValHist->data[(static_cast<int>(b_index) + ForceValHist->size[0] *
                            i1) - 1] = spanLocNorm->data[b_i + spanLocNorm->
          size[0] * i1];
      }

      ForceDof->data[static_cast<int>(b_index) - 1] = static_cast<double>(b_i) +
        1.0;
      b_index++;
    }
  }

  emxFree_real_T(&b_spanLocNorm);
  emxFree_real_T(&spanLocNorm);
}

//
// File trailer for mapCactusLoadsFile.cpp
//
// [EOF]
//
