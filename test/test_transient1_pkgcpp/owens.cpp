//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: owens.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "owens.h"
#include "calculateBCMap.h"
#include "createJointTransform.h"
#include "fileManager.h"
#include "find.h"
#include "getSplitLine.h"
#include "myfgetl.h"
#include "readBCdata.h"
#include "readBladeData.h"
#include "readElementData.h"
#include "readGeneratorProps.h"
#include "readInitCond.h"
#include "readJointData.h"
#include "readMesh.h"
#include "readNLParamsFile.h"
#include "readNodalTerms.h"
#include "rt_nonfinite.h"
#include "str2double.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include "transientExec.h"
#include <cmath>
#include <cstring>
#include <math.h>
#include <string.h>

// Function Declarations
static void constructReducedDispVectorMap(double numNodes, double numReducedDof,
  double BC_numpBC, const emxArray_real_T *BC_pBC, const emxArray_real_T
  *BC_isConstrained, emxArray_real_T *redVectorMap);

// Function Definitions

//
// This function creates a map of unconstrained DOFs between a full
// listing and reduced listing (aftger constraints have been applied)
// Arguments    : double numNodes
//                double numReducedDof
//                double BC_numpBC
//                const emxArray_real_T *BC_pBC
//                const emxArray_real_T *BC_isConstrained
//                emxArray_real_T *redVectorMap
// Return Type  : void
//
static void constructReducedDispVectorMap(double numNodes, double numReducedDof,
  double BC_numpBC, const emxArray_real_T *BC_pBC, const emxArray_real_T
  *BC_isConstrained, emxArray_real_T *redVectorMap)
{
  emxArray_real_T *bcdoflist;
  int loop_ub;
  int i;
  int b_i;
  double b_index;
  emxArray_real_T *dofList;
  int i1;
  double d;
  emxArray_boolean_T *tf;
  emxArray_int32_T *ii;
  boolean_T b_tf;
  boolean_T exitg1;
  int na;
  unsigned int unnamed_idx_1;
  int exponent;
  int b_exponent;
  emxInit_real_T(&bcdoflist, 1);
  loop_ub = static_cast<int>(BC_numpBC);
  i = bcdoflist->size[0];
  bcdoflist->size[0] = loop_ub;
  emxEnsureCapacity_real_T(bcdoflist, i);

  // create a listing of constrained DOFs from boundary condition file
  for (b_i = 0; b_i < loop_ub; b_i++) {
    bcdoflist->data[b_i] = 0.0;
    bcdoflist->data[b_i] = (BC_pBC->data[b_i] - 1.0) * 6.0 + BC_pBC->data[b_i +
      BC_pBC->size[0]];
  }

  // This function searches over all DOFs in a structural model and
  // determines and returns "dofVector" containing only unconstrained DOFs
  // loop over all DOFs in the model checking if constrained by BC or not
  b_index = 1.0;
  i = static_cast<int>(numNodes);
  for (b_i = 0; b_i < i; b_i++) {
    for (loop_ub = 0; loop_ub < 6; loop_ub++) {
      if (!(BC_isConstrained->data[static_cast<int>((((static_cast<double>(b_i)
               + 1.0) - 1.0) * 6.0 + (static_cast<double>(loop_ub) + 1.0))) - 1]
            != 0.0)) {
        //              dofVector(index) = (i-1)*numDofPerNode + j; %DOF vector only contains unconstrained DOFs 
        b_index++;
      }
    }
  }

  emxInit_real_T(&dofList, 2);
  i1 = dofList->size[0] * dofList->size[1];
  dofList->size[0] = 1;
  loop_ub = static_cast<int>(b_index);
  dofList->size[1] = loop_ub;
  emxEnsureCapacity_real_T(dofList, i1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    dofList->data[i1] = 0.0;
  }

  b_index = 1.0;
  for (b_i = 0; b_i < i; b_i++) {
    for (loop_ub = 0; loop_ub < 6; loop_ub++) {
      d = ((static_cast<double>(b_i) + 1.0) - 1.0) * 6.0 + (static_cast<double>
        (loop_ub) + 1.0);
      if (!(BC_isConstrained->data[static_cast<int>(d) - 1] != 0.0)) {
        dofList->data[static_cast<int>(b_index) - 1] = d;

        // DOF vector only contains unconstrained DOFs
        b_index++;
      }
    }
  }

  // calculate a reduced (unconstrained) DOF vector
  i = static_cast<int>(numReducedDof);
  i1 = redVectorMap->size[0];
  redVectorMap->size[0] = i;
  emxEnsureCapacity_real_T(redVectorMap, i1);
  emxInit_boolean_T(&tf, 2);
  emxInit_int32_T(&ii, 2);
  for (b_i = 0; b_i < i; b_i++) {
    b_tf = false;
    loop_ub = 0;
    exitg1 = false;
    while ((!exitg1) && (loop_ub <= bcdoflist->size[0] - 1)) {
      b_index = std::abs(bcdoflist->data[loop_ub] / 2.0);
      if ((!rtIsInf(b_index)) && (!rtIsNaN(b_index))) {
        if (b_index <= 2.2250738585072014E-308) {
          b_index = 4.94065645841247E-324;
        } else {
          frexp(b_index, &exponent);
          b_index = std::ldexp(1.0, exponent - 53);
        }
      } else {
        b_index = rtNaN;
      }

      if (std::abs(bcdoflist->data[loop_ub] - (static_cast<double>(b_i) + 1.0)) <
          b_index) {
        b_tf = true;
        exitg1 = true;
      } else {
        loop_ub++;
      }
    }

    if (b_tf) {
      // creates a map of unconstrained reduced DOFs
      redVectorMap->data[b_i] = -1.0;
    } else {
      na = dofList->size[1];
      unnamed_idx_1 = static_cast<unsigned int>(dofList->size[1]);
      i1 = tf->size[0] * tf->size[1];
      tf->size[0] = 1;
      tf->size[1] = static_cast<int>(unnamed_idx_1);
      emxEnsureCapacity_boolean_T(tf, i1);
      loop_ub = static_cast<int>(unnamed_idx_1);
      for (i1 = 0; i1 < loop_ub; i1++) {
        tf->data[i1] = false;
      }

      for (loop_ub = 0; loop_ub < na; loop_ub++) {
        frexp((static_cast<double>(b_i) + 1.0) / 2.0, &b_exponent);
        if (std::abs((static_cast<double>(b_i) + 1.0) - dofList->data[loop_ub]) <
            std::ldexp(1.0, b_exponent - 53)) {
          tf->data[loop_ub] = true;
        }
      }

      i1 = tf->size[1];
      na = 0;
      loop_ub = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = tf->size[1];
      emxEnsureCapacity_int32_T(ii, loop_ub);
      loop_ub = 0;
      exitg1 = false;
      while ((!exitg1) && (loop_ub <= i1 - 1)) {
        if (tf->data[loop_ub]) {
          na++;
          ii->data[na - 1] = loop_ub + 1;
          if (na >= i1) {
            exitg1 = true;
          } else {
            loop_ub++;
          }
        } else {
          loop_ub++;
        }
      }

      if (tf->size[1] == 1) {
        if (na == 0) {
          ii->size[0] = 1;
          ii->size[1] = 0;
        }
      } else {
        i1 = ii->size[0] * ii->size[1];
        if (1 > na) {
          ii->size[1] = 0;
        } else {
          ii->size[1] = na;
        }

        emxEnsureCapacity_int32_T(ii, i1);
      }

      redVectorMap->data[b_i] = ii->data[0];
    }
  }

  emxFree_int32_T(&ii);
  emxFree_boolean_T(&tf);
  emxFree_real_T(&dofList);
  emxFree_real_T(&bcdoflist);
}

//
// owens Startup function for the OWENS toolkit
//  **********************************************************************
//  *                   Part of the SNL OWENS toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [freq,damp]=owens(varargin)
//
//    This function is a start up function for launching various analysis
//    modes of the OWENS toolkit.
//
//       input:
//       varargin      = input parameter list
//          varargin{1} is the .owens file associated with analysis
//          varargin{2} is a string describing analysis type
//                      'S' = static analysis
//                      'M' = modal analysis
//                      'TNB' = transient analysis with Newmark-Beta time
//                      integration
//                      'ROM' = reduced order model transient analysis
//       output:
//       freq         = array of modal frequencies (when applicable to
//                      analysis type)
//       damp         = array of modal damping (when applicable to analysis
//                      type)
//       displ        = array containing converged solution for static
//                      displacement
// Arguments    : const double varargin_7[2]
// Return Type  : void
//
void owens(const double varargin_7[2])
{
  signed char fileid;
  int last_delimiter_data[73];
  int last_delimiter_size[2];
  int loop_ub;
  char fdirectory_data[73];
  static const char cv[73] = { '.', '/', 'i', 'n', 'p', 'u', 't', '_', 'f', 'i',
    'l', 'e', 's', '_', 't', 'e', 's', 't', '/', '1', '_', 'F', 'o', 'u', 'r',
    'C', 'o', 'l', 'u', 'm', 'n', 'S', 'e', 'm', 'i', '_', '2', 'n', 'd', 'P',
    'a', 's', 's', '_', '1', '5', 'm', 'T', 'o', 'w', 'e', 'r', 'E', 'x', 't',
    '_', 'N', 'O', 'c', 'e', 'n', 't', 'S', 't', 'i', 'f', 'f', '.', 'o', 'w',
    'e', 'n', 's' };

  emxArray_char_T *line;
  emxArray_char_T *varargin_2;
  emxArray_char_T *b_varargin_2;
  emxArray_char_T *c_varargin_2;
  emxArray_char_T *d_varargin_2;
  emxArray_char_T *e_varargin_2;
  emxArray_char_T *f_varargin_2;
  emxArray_char_T *g_varargin_2;
  emxArray_boolean_T *b_line;
  signed char input_sizes_idx_1;
  int sizes_idx_1;
  signed char b_input_sizes_idx_1;
  signed char c_input_sizes_idx_1;
  int b_sizes_idx_1;
  signed char d_input_sizes_idx_1;
  int c_sizes_idx_1;
  signed char e_input_sizes_idx_1;
  int d_sizes_idx_1;
  signed char f_input_sizes_idx_1;
  int e_sizes_idx_1;
  signed char g_input_sizes_idx_1;
  int f_sizes_idx_1;
  int i;
  int b_loop_ub;
  emxArray_real_T *delimiter_idx;
  int g_sizes_idx_1;
  emxArray_int32_T *b_r;
  int i1;
  emxArray_char_T *blddatafilename;
  emxArray_char_T *model_aeroloadfile;
  emxArray_char_T *model_owensfile;
  emxArray_char_T *generatorfilename;
  static const char cv1[6] = { '.', 'o', 'w', 'e', 'n', 's' };

  signed char h_input_sizes_idx_1;
  int unnamed_idx_1;
  emxArray_char_T *fdirectory;
  double model_RayleighAlpha;
  double model_RayleighBeta;
  emxArray_real_T *mesh_x;
  emxArray_real_T *mesh_y;
  emxArray_real_T *mesh_z;
  emxArray_real_T *mesh_conn;
  emxArray_real_T *expl_temp;
  double mesh_numEl;
  double mesh_numNodes;
  double bladeData_numBlades;
  double bladeData_bladeNum_data[60];
  int bladeData_bladeNum_size[1];
  double bladeData_h_data[60];
  int bladeData_h_size[1];
  double bladeData_nodeNum_data[60];
  int bladeData_nodeNum_size[1];
  double bladeData_elementNum_data[60];
  int bladeData_elementNum_size[1];
  double bladeData_remaining_data[720];
  emxArray_real_T *model_BC_pBC;
  emxArray_real_T *BC_pBC;
  double BC_numpBC;
  double b_expl_temp;
  emxArray_struct_T *el_props;
  emxArray_real_T *el_psi;
  emxArray_real_T *el_theta;
  emxArray_real_T *el_roll;
  emxArray_boolean_T *el_rotationalEffects;
  emxArray_real_T *joint;
  creal_T dc;
  emxArray_real_T *jnt_struct_jointTransform;
  emxArray_real_T *jnt_struct_reducedDOF;
  e_struct_T c_expl_temp;
  emxArray_real_T *b_BC_numpBC;
  emxArray_real_T *b_mesh_numNodes;
  char model_nlParams_iterationType[2];
  boolean_T c_model_nlParams_adaptiveLoadSt;
  double model_nlParams_maxNumLoadSteps;
  double model_nlParams_minLoadStepDelta;
  double model_nlParams_minLoadStep;
  double c_model_nlParams_prescribedLoad;
  f_struct_T d_expl_temp;
  static const char model_analysisType[3] = { 'T', 'N', 'B' };

  //  Initialize Model Struct
  //  model = struct('analysisType',char,...
  //  'turbineStartup',false,...
  //  'aeroElasticOn',false,...
  //  'airDensity',zeros,...
  //  'gravityOn',false,...
  //  'nlOn',false,...
  //  'spinUpOn',false,...
  //  'numModesToExtract',zeros,...
  //  'delta_t',zeros,...
  //  'numTS',zeros,...
  //  'OmegaInit',zeros,...
  //  'OmegaGenStart',zeros,...
  //  'usingRotorSpeedFunction',false,...
  //  'tocp',zeros,...
  //  'Omegaocp',zeros,...
  //  'numModesForROM',zeros,...
  //  'guessFreq',zeros,...
  //  'aeroloadfile',char,...
  //  'owensfile',char,...
  //  'RayleighAlpha',zeros,...
  //  'RayleighBeta',zeros,...
  //  'elementOrder',zeros,...
  //  'bladeData',zeros,...
  //  'BC',zeros,...
  //  'joint',zeros,...
  //  'nodalTerms',zeros,...
  //  'initCond',zeros,...
  //  'aeroLoadsOn',zeros,...
  //  'useGeneratorFunction',zeros,...
  //  'generatorProps',zeros,...
  //  'outFilename',zeros,...
  //  'jointTransform',zeros,...
  //  'reducedDOFList',zeros,...
  //  'nlParams',zeros,...
  //  'analysisType',zeros);
  // input file initialization
  // anaysis type intialization
  // initialization of turbine startup,
  //  aeroElastic flags, and air density
  // flag to activate gravity loading in structural dynamics/static simulations
  // Initialize only, gets changed later on
  // Initialize only, gets changed later on
  // Initialize only, gets changed later on
  // Initialize only, gets changed later on
  // TRANSIENT ANALYSIS (TNB = newmark beta time integation, TD =  dean time integration) 
  //  time step size
  //  number of time steps
  //  flag for nonlinear elastic calculation
  // turbine operation flag
  // specified rotor speed profile
  // this option uses a discretely specified rotor speed profile
  // set flag to not use user specified rotor speed function
  // time points for rotor speed provfile
  // rotor speed value at time points (Hz)
  fileid = cfopen("./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.owens",
                  "rb");

  // reads in model file names from .owens file
  eml_find(last_delimiter_data, last_delimiter_size);

  // '
  loop_ub = last_delimiter_data[last_delimiter_size[1] - 1];
  if (1 > loop_ub) {
    loop_ub = 0;
  }

  if (0 <= loop_ub - 1) {
    std::memcpy(&fdirectory_data[0], &cv[0], loop_ub * sizeof(char));
  }

  emxInit_char_T(&line, 2);
  emxInit_char_T(&varargin_2, 2);
  emxInit_char_T(&b_varargin_2, 2);
  emxInit_char_T(&c_varargin_2, 2);
  emxInit_char_T(&d_varargin_2, 2);
  emxInit_char_T(&e_varargin_2, 2);
  emxInit_char_T(&f_varargin_2, 2);
  emxInit_char_T(&g_varargin_2, 2);
  emxInit_boolean_T(&b_line, 2);
  myfgetl(static_cast<double>(fileid), varargin_2);
  if (loop_ub != 0) {
    input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    input_sizes_idx_1 = 0;
  }

  if ((varargin_2->size[0] != 0) && (varargin_2->size[1] != 0)) {
    sizes_idx_1 = varargin_2->size[1];
  } else {
    sizes_idx_1 = 0;
  }

  // mesh file name
  myfgetl(static_cast<double>(fileid), b_varargin_2);
  if (loop_ub != 0) {
    b_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    b_input_sizes_idx_1 = 0;
  }

  // element data file name
  myfgetl(static_cast<double>(fileid), c_varargin_2);
  if (loop_ub != 0) {
    c_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    c_input_sizes_idx_1 = 0;
  }

  if ((c_varargin_2->size[0] != 0) && (c_varargin_2->size[1] != 0)) {
    b_sizes_idx_1 = c_varargin_2->size[1];
  } else {
    b_sizes_idx_1 = 0;
  }

  // element orientation file name
  myfgetl(static_cast<double>(fileid), d_varargin_2);
  if (loop_ub != 0) {
    d_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    d_input_sizes_idx_1 = 0;
  }

  if ((d_varargin_2->size[0] != 0) && (d_varargin_2->size[1] != 0)) {
    c_sizes_idx_1 = d_varargin_2->size[1];
  } else {
    c_sizes_idx_1 = 0;
  }

  // joint data file name
  myfgetl(static_cast<double>(fileid), e_varargin_2);
  if (loop_ub != 0) {
    e_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    e_input_sizes_idx_1 = 0;
  }

  if ((e_varargin_2->size[0] != 0) && (e_varargin_2->size[1] != 0)) {
    d_sizes_idx_1 = e_varargin_2->size[1];
  } else {
    d_sizes_idx_1 = 0;
  }

  // concentrated nodal data file name
  myfgetl(static_cast<double>(fileid), f_varargin_2);
  if (loop_ub != 0) {
    f_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    f_input_sizes_idx_1 = 0;
  }

  if ((f_varargin_2->size[0] != 0) && (f_varargin_2->size[1] != 0)) {
    e_sizes_idx_1 = f_varargin_2->size[1];
  } else {
    e_sizes_idx_1 = 0;
  }

  // boundary condition file name
  myfgetl(static_cast<double>(fileid), line);
  str2double(*(char (*)[2])&line->data[0]);
  myfgetl(static_cast<double>(fileid), g_varargin_2);
  if (loop_ub != 0) {
    g_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    g_input_sizes_idx_1 = 0;
  }

  if ((g_varargin_2->size[0] != 0) && (g_varargin_2->size[1] != 0)) {
    f_sizes_idx_1 = g_varargin_2->size[1];
  } else {
    f_sizes_idx_1 = 0;
  }

  // initial condition filename
  myfgetl(static_cast<double>(fileid), line);
  i = b_line->size[0] * b_line->size[1];
  b_line->size[0] = line->size[1];
  b_line->size[1] = line->size[0];
  emxEnsureCapacity_boolean_T(b_line, i);
  b_loop_ub = line->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    g_sizes_idx_1 = line->size[1];
    for (i1 = 0; i1 < g_sizes_idx_1; i1++) {
      b_line->data[i1 + b_line->size[0] * i] = (line->data[i + line->size[0] *
        i1] == ' ');
    }
  }

  emxInit_real_T(&delimiter_idx, 1);
  emxInit_int32_T(&b_r, 1);
  b_eml_find(b_line, b_r);
  i = delimiter_idx->size[0];
  delimiter_idx->size[0] = b_r->size[0];
  emxEnsureCapacity_real_T(delimiter_idx, i);
  b_loop_ub = b_r->size[0];
  emxFree_boolean_T(&b_line);
  for (i = 0; i < b_loop_ub; i++) {
    delimiter_idx->data[i] = b_r->data[i];
  }

  emxFree_int32_T(&b_r);
  b_str2double(line->data[0]);

  // flag for activating aerodynamic analysis
  if (delimiter_idx->data[0] + 1.0 > delimiter_idx->data[1] - 1.0) {
    i = 0;
    i1 = 0;
  } else {
    i = static_cast<int>((delimiter_idx->data[0] + 1.0)) - 1;
    i1 = static_cast<int>((delimiter_idx->data[1] - 1.0));
  }

  emxInit_char_T(&blddatafilename, 2);
  g_sizes_idx_1 = blddatafilename->size[0] * blddatafilename->size[1];
  blddatafilename->size[0] = 1;
  blddatafilename->size[1] = (loop_ub + i1) - i;
  emxEnsureCapacity_char_T(blddatafilename, g_sizes_idx_1);
  for (g_sizes_idx_1 = 0; g_sizes_idx_1 < loop_ub; g_sizes_idx_1++) {
    blddatafilename->data[g_sizes_idx_1] = fdirectory_data[g_sizes_idx_1];
  }

  b_loop_ub = i1 - i;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    blddatafilename->data[i1 + loop_ub] = line->data[i + i1];
  }

  // blade data file name
  i = line->size[0] * line->size[1];
  if (delimiter_idx->data[1] + 1.0 > i) {
    i1 = 0;
    i = 0;
  } else {
    i1 = static_cast<int>((delimiter_idx->data[1] + 1.0)) - 1;
  }

  emxInit_char_T(&model_aeroloadfile, 2);
  g_sizes_idx_1 = model_aeroloadfile->size[0] * model_aeroloadfile->size[1];
  model_aeroloadfile->size[0] = 1;
  model_aeroloadfile->size[1] = (loop_ub + i) - i1;
  emxEnsureCapacity_char_T(model_aeroloadfile, g_sizes_idx_1);
  for (g_sizes_idx_1 = 0; g_sizes_idx_1 < loop_ub; g_sizes_idx_1++) {
    model_aeroloadfile->data[g_sizes_idx_1] = fdirectory_data[g_sizes_idx_1];
  }

  b_loop_ub = i - i1;
  for (i = 0; i < b_loop_ub; i++) {
    model_aeroloadfile->data[i + loop_ub] = line->data[i1 + i];
  }

  emxInit_char_T(&model_owensfile, 2);

  // .csv file containing CACTUS aerodynamic loads
  if (1 > blddatafilename->size[1] - 4) {
    b_loop_ub = 0;
  } else {
    b_loop_ub = blddatafilename->size[1] - 4;
  }

  i = model_owensfile->size[0] * model_owensfile->size[1];
  model_owensfile->size[0] = 1;
  model_owensfile->size[1] = b_loop_ub + 6;
  emxEnsureCapacity_char_T(model_owensfile, i);
  for (i = 0; i < b_loop_ub; i++) {
    model_owensfile->data[i] = blddatafilename->data[i];
  }

  for (i = 0; i < 6; i++) {
    model_owensfile->data[i + b_loop_ub] = cv1[i];
  }

  emxInit_char_T(&generatorfilename, 2);
  myfgetl(static_cast<double>(fileid), line);

  // flag to include drive shaft effects
  str2double(*(char (*)[2])&line->data[0]);

  // drive shaft file name
  myfgetl(static_cast<double>(fileid), line);
  if (loop_ub != 0) {
    h_input_sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    h_input_sizes_idx_1 = 0;
  }

  if ((line->size[0] != 0) && (line->size[1] != 0)) {
    g_sizes_idx_1 = line->size[1];
  } else {
    g_sizes_idx_1 = 0;
  }

  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = generatorfilename->size[0] * generatorfilename->size[1];
  generatorfilename->size[0] = 1;
  generatorfilename->size[1] = h_input_sizes_idx_1 + g_sizes_idx_1;
  emxEnsureCapacity_char_T(generatorfilename, i);
  b_loop_ub = h_input_sizes_idx_1;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      generatorfilename->data[generatorfilename->size[0] * i] =
        fdirectory_data[i];
    }
  }

  for (i = 0; i < g_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      generatorfilename->data[generatorfilename->size[0] * (i + unnamed_idx_1)] =
        line->data[line->size[0] * i];
    }
  }

  emxInit_char_T(&fdirectory, 2);

  // generator file name
  getSplitLine(static_cast<double>(fileid), delimiter_idx);

  // read in alpha/beta for rayleigh damping
  model_RayleighAlpha = delimiter_idx->data[0];
  model_RayleighBeta = delimiter_idx->data[1];
  cfclose(static_cast<double>(fileid));

  //  close .owens file
  // model definitions
  // linear element order
  // --------------------------------------------
  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  b_loop_ub = input_sizes_idx_1;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        varargin_2->data[varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&varargin_2);
  emxInit_real_T(&mesh_x, 1);
  emxInit_real_T(&mesh_y, 1);
  emxInit_real_T(&mesh_z, 1);
  emxInit_real_T(&mesh_conn, 2);
  emxInit_real_T(&expl_temp, 1);
  readMesh(fdirectory, delimiter_idx, &mesh_numEl, &mesh_numNodes, mesh_x,
           mesh_y, mesh_z, expl_temp, mesh_conn);

  // read mesh file
  readBladeData(blddatafilename, &bladeData_numBlades, bladeData_bladeNum_data,
                bladeData_bladeNum_size, bladeData_h_data, bladeData_h_size,
                bladeData_nodeNum_data, bladeData_nodeNum_size,
                bladeData_elementNum_data, bladeData_elementNum_size,
                bladeData_remaining_data, last_delimiter_size);

  // reads overall blade data file
  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = f_input_sizes_idx_1 + e_sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  b_loop_ub = f_input_sizes_idx_1;
  emxFree_char_T(&blddatafilename);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < e_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        f_varargin_2->data[f_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&f_varargin_2);
  emxInit_real_T(&model_BC_pBC, 2);
  emxInit_real_T(&BC_pBC, 2);
  readBCdata(fdirectory, mesh_numNodes, &BC_numpBC, BC_pBC, &bladeData_numBlades,
             &b_expl_temp, delimiter_idx);

  // read boundary condition file
  i = model_BC_pBC->size[0] * model_BC_pBC->size[1];
  model_BC_pBC->size[0] = BC_pBC->size[0];
  model_BC_pBC->size[1] = BC_pBC->size[1];
  emxEnsureCapacity_real_T(model_BC_pBC, i);
  b_loop_ub = BC_pBC->size[0] * BC_pBC->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    model_BC_pBC->data[i] = BC_pBC->data[i];
  }

  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  if ((b_varargin_2->size[0] != 0) && (b_varargin_2->size[1] != 0)) {
    g_sizes_idx_1 = b_varargin_2->size[1];
  } else {
    g_sizes_idx_1 = 0;
  }

  if (loop_ub != 0) {
    sizes_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    sizes_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = b_input_sizes_idx_1 + g_sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  b_loop_ub = b_input_sizes_idx_1;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < g_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        b_varargin_2->data[b_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&b_varargin_2);
  i = line->size[0] * line->size[1];
  line->size[0] = 1;
  line->size[1] = c_input_sizes_idx_1 + b_sizes_idx_1;
  emxEnsureCapacity_char_T(line, i);
  b_loop_ub = c_input_sizes_idx_1;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      line->data[line->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < b_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      line->data[line->size[0] * (i + sizes_idx_1)] = c_varargin_2->
        data[c_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&c_varargin_2);
  emxInit_struct_T(&el_props, 2);
  emxInit_real_T(&el_psi, 1);
  emxInit_real_T(&el_theta, 1);
  emxInit_real_T(&el_roll, 1);
  emxInit_boolean_T(&el_rotationalEffects, 2);
  readElementData(mesh_numEl, fdirectory, line, bladeData_nodeNum_data,
                  bladeData_nodeNum_size, bladeData_elementNum_data,
                  bladeData_elementNum_size, bladeData_remaining_data,
                  last_delimiter_size, el_props, expl_temp, el_psi, el_theta,
                  el_roll, el_rotationalEffects);

  // read element data file (also reads orientation and blade data file associated with elements) 
  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = d_input_sizes_idx_1 + c_sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  b_loop_ub = d_input_sizes_idx_1;
  emxFree_char_T(&line);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < c_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        d_varargin_2->data[d_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&d_varargin_2);
  emxInit_real_T(&joint, 2);
  readJointData(fdirectory, joint);

  // read joint data file
  //  rbarFileName = [inputfile(1:end-6),'.rbar']; %setrbarfile
  //  [model.joint] = readRBarFile(rbarFileName,model.joint,mesh); %read rbar file name 
  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = e_input_sizes_idx_1 + d_sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  b_loop_ub = e_input_sizes_idx_1;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < d_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        e_varargin_2->data[e_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&e_varargin_2);
  readNodalTerms(fdirectory);

  // read concentrated nodal terms file
  //  [model] = readPlatformFile(model,platformFlag,platfilename);
  // for transient analysis...
  if (loop_ub != 0) {
    unnamed_idx_1 = static_cast<signed char>(loop_ub);
  } else {
    unnamed_idx_1 = 0;
  }

  i = fdirectory->size[0] * fdirectory->size[1];
  fdirectory->size[0] = 1;
  fdirectory->size[1] = g_input_sizes_idx_1 + f_sizes_idx_1;
  emxEnsureCapacity_char_T(fdirectory, i);
  loop_ub = g_input_sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * i] = fdirectory_data[i];
    }
  }

  for (i = 0; i < f_sizes_idx_1; i++) {
    for (i1 = 0; i1 < 1; i1++) {
      fdirectory->data[fdirectory->size[0] * (i + unnamed_idx_1)] =
        g_varargin_2->data[g_varargin_2->size[0] * i];
    }
  }

  emxFree_char_T(&g_varargin_2);
  readInitCond(fdirectory);

  // read initial conditions
  //      [model] = readDriveShaftProps(model,driveShaftFlag,driveshaftfilename); %reads drive shaft properties 
  // set drive shaft unactive
  // set drive shat properties to 0
  // set gear ratio and efficiency to 1
  dc = d_str2double(generatorfilename);
  emxFree_char_T(&fdirectory);
  if (!(dc.re == 1.0)) {
    readGeneratorProps(generatorfilename);

    // reads generator properties
  }

  emxFree_char_T(&generatorfilename);
  emxInit_real_T(&jnt_struct_jointTransform, 2);
  emxInit_real_T(&jnt_struct_reducedDOF, 2);
  emxInitStruct_struct_T(&c_expl_temp);
  emxInit_real_T(&b_BC_numpBC, 1);
  emxInit_real_T(&b_mesh_numNodes, 1);

  // generates an output filename for analysis results %TODO: map to the output location instead of input 
  createJointTransform(joint, mesh_numNodes, jnt_struct_jointTransform,
                       jnt_struct_reducedDOF);

  // creates a joint transform to constrain model degrees of freedom (DOF) consistent with joint constraints 
  calculateBCMap(BC_numpBC, BC_pBC, jnt_struct_reducedDOF, b_BC_numpBC);

  // create boundary condition map from original DOF numbering to reduced/constrained DOF numbering 
  constructReducedDispVectorMap(mesh_numNodes, static_cast<double>
    (jnt_struct_jointTransform->size[1]), BC_numpBC, model_BC_pBC, delimiter_idx,
    b_mesh_numNodes);

  // create a map between reduced and full DOF lists
  // EXECUTE TRANSIENT ANALYSIS
  readNLParamsFile(model_nlParams_iterationType,
                   &c_model_nlParams_adaptiveLoadSt, &bladeData_numBlades,
                   &b_expl_temp, &model_nlParams_maxNumLoadSteps,
                   &model_nlParams_minLoadStepDelta, &model_nlParams_minLoadStep,
                   &c_model_nlParams_prescribedLoad);
  i = c_expl_temp.conn->size[0] * c_expl_temp.conn->size[1];
  c_expl_temp.conn->size[0] = mesh_conn->size[0];
  c_expl_temp.conn->size[1] = 2;
  emxEnsureCapacity_real_T(c_expl_temp.conn, i);
  loop_ub = mesh_conn->size[0] * mesh_conn->size[1];
  emxFree_real_T(&b_mesh_numNodes);
  emxFree_real_T(&b_BC_numpBC);
  emxFree_real_T(&jnt_struct_reducedDOF);
  emxFree_real_T(&BC_pBC);
  emxFree_real_T(&delimiter_idx);
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.conn->data[i] = mesh_conn->data[i];
  }

  emxFree_real_T(&mesh_conn);
  i = c_expl_temp.z->size[0];
  c_expl_temp.z->size[0] = mesh_z->size[0];
  emxEnsureCapacity_real_T(c_expl_temp.z, i);
  loop_ub = mesh_z->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.z->data[i] = mesh_z->data[i];
  }

  emxFree_real_T(&mesh_z);
  i = c_expl_temp.y->size[0];
  c_expl_temp.y->size[0] = mesh_y->size[0];
  emxEnsureCapacity_real_T(c_expl_temp.y, i);
  loop_ub = mesh_y->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.y->data[i] = mesh_y->data[i];
  }

  emxFree_real_T(&mesh_y);
  i = c_expl_temp.x->size[0];
  c_expl_temp.x->size[0] = mesh_x->size[0];
  emxEnsureCapacity_real_T(c_expl_temp.x, i);
  loop_ub = mesh_x->size[0];
  for (i = 0; i < loop_ub; i++) {
    c_expl_temp.x->data[i] = mesh_x->data[i];
  }

  emxFree_real_T(&mesh_x);
  emxInitStruct_struct_T1(&d_expl_temp);
  c_expl_temp.numNodes = mesh_numNodes;
  c_expl_temp.numEl = mesh_numEl;
  i = d_expl_temp.rotationalEffects->size[0] *
    d_expl_temp.rotationalEffects->size[1];
  d_expl_temp.rotationalEffects->size[0] = 1;
  d_expl_temp.rotationalEffects->size[1] = el_rotationalEffects->size[1];
  emxEnsureCapacity_boolean_T(d_expl_temp.rotationalEffects, i);
  loop_ub = el_rotationalEffects->size[0] * el_rotationalEffects->size[1];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.rotationalEffects->data[i] = el_rotationalEffects->data[i];
  }

  emxFree_boolean_T(&el_rotationalEffects);
  i = d_expl_temp.roll->size[0];
  d_expl_temp.roll->size[0] = el_roll->size[0];
  emxEnsureCapacity_real_T(d_expl_temp.roll, i);
  loop_ub = el_roll->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.roll->data[i] = el_roll->data[i];
  }

  emxFree_real_T(&el_roll);
  i = d_expl_temp.theta->size[0];
  d_expl_temp.theta->size[0] = el_theta->size[0];
  emxEnsureCapacity_real_T(d_expl_temp.theta, i);
  loop_ub = el_theta->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.theta->data[i] = el_theta->data[i];
  }

  emxFree_real_T(&el_theta);
  i = d_expl_temp.psi->size[0];
  d_expl_temp.psi->size[0] = el_psi->size[0];
  emxEnsureCapacity_real_T(d_expl_temp.psi, i);
  loop_ub = el_psi->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.psi->data[i] = el_psi->data[i];
  }

  emxFree_real_T(&el_psi);
  i = d_expl_temp.elLen->size[0];
  d_expl_temp.elLen->size[0] = expl_temp->size[0];
  emxEnsureCapacity_real_T(d_expl_temp.elLen, i);
  loop_ub = expl_temp->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.elLen->data[i] = expl_temp->data[i];
  }

  emxFree_real_T(&expl_temp);
  i = d_expl_temp.props->size[0] * d_expl_temp.props->size[1];
  d_expl_temp.props->size[0] = 1;
  d_expl_temp.props->size[1] = el_props->size[1];
  emxEnsureCapacity_struct_T(d_expl_temp.props, i);
  loop_ub = el_props->size[0] * el_props->size[1];
  for (i = 0; i < loop_ub; i++) {
    d_expl_temp.props->data[i] = el_props->data[i];
  }

  emxFree_struct_T(&el_props);
  transientExec(model_analysisType, varargin_7, model_aeroloadfile,
                model_owensfile, model_RayleighAlpha, model_RayleighBeta,
                BC_numpBC, model_BC_pBC, joint, jnt_struct_jointTransform,
                c_expl_temp, d_expl_temp);
  emxFreeStruct_struct_T1(&d_expl_temp);
  emxFreeStruct_struct_T(&c_expl_temp);
  emxFree_real_T(&jnt_struct_jointTransform);
  emxFree_real_T(&joint);
  emxFree_real_T(&model_BC_pBC);
  emxFree_char_T(&model_owensfile);
  emxFree_char_T(&model_aeroloadfile);
}

//
// File trailer for owens.cpp
//
// [EOF]
//
