//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: elementPostProcess.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "elementPostProcess.h"
#include "ConcMassAssociatedWithElement.h"
#include "calculateTimoshenkoElementNL.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include "test_owens.h"
#include <cstring>
#include <string.h>

// Function Definitions

//
// elementPostProcess post processes element for reaction force
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [Fpp] = elementPostProcess(elementNumber,model,mesh,el,elStorage,....
//            timeInt,dispData,displ_iter,rbData,Omega,OmegaDot,CN2H)
//
//    This function calculates the reaction force associated with an element.
//
//    input:
//    elementNumber  = node number joint constraints are desired at
//    model          = object containing model data
//    mesh           = object containing mesh data
//    elStorage      = object containing stored element data
//    el             = object containing element data
//    timeInt        = object containing time integration parameters
//    dispData       = object containing displacement data
//    displ_iter     = converged displacement solution
//    rbData         = vector containing rigid body displacement, velocity,
//                      and acceleration
//    Omega          = rotor speed (Hz)
//    OmegaDot       = rotor acceleratin (Hz)
//    CN2H           = transformation matrix from inertial frame to hub frame
//
//    output:
//    Fpp            = vector containing reaction force vector associated
//                     with element
// Arguments    : double elementNumber
//                const char model_analysisType[3]
//                double model_RayleighAlpha
//                double model_RayleighBeta
//                const emxArray_real_T *model_joint
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                const emxArray_struct_T *el_props
//                const emxArray_real_T *el_elLen
//                const emxArray_real_T *el_psi
//                const emxArray_real_T *el_theta
//                const emxArray_real_T *el_roll
//                const c_emxArray_struct_T *elStorage
//                const b_struct_T *timeInt
//                const i_struct_T dispData
//                const emxArray_real_T *displ_iter
//                const double rbData[9]
//                double Omega
//                double OmegaDot
//                const double CN2H[9]
//                double Fpp[12]
// Return Type  : void
//
void elementPostProcess(double elementNumber, const char model_analysisType[3],
  double model_RayleighAlpha, double model_RayleighBeta, const emxArray_real_T
  *model_joint, const emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y,
  const emxArray_real_T *mesh_z, const emxArray_real_T *mesh_conn, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, const c_emxArray_struct_T *elStorage, const
  b_struct_T *timeInt, const i_struct_T dispData, const emxArray_real_T
  *displ_iter, const double rbData[9], double Omega, double OmegaDot, const
  double CN2H[9], double Fpp[12])
{
  double elInput_z_data[4];
  double eldisp_data[12];
  double elInput_dispdot_data[12];
  double elInput_dispddot_data[12];
  double eldispiter_data[12];
  double b_index;
  int elInput_xloc_idx_1_tmp;
  int elInput_analysisType_size_idx_1;
  char elInput_analysisType_data[3];
  int j;
  double d;
  int elx_tmp;
  double elx_data[2];
  double ely_data[2];
  double b_mesh_conn[2];
  int k;
  double elInput_concMass[8];
  int aoffset;
  int elInput_omegaVec_size_idx_0;
  int elInput_omegaVec_size_idx_1;
  double elInput_accelVec_data[9];
  double elInput_omegaVec_data[9];
  double elInput_omegaDotVec_data[9];
  p_struct_T expl_temp;
  o_struct_T b_expl_temp;
  double elOutput_Ke[144];

  // some initializations
  elInput_z_data[0] = 0.0;
  elInput_z_data[1] = 0.0;
  std::memset(&eldisp_data[0], 0, 12U * sizeof(double));
  std::memset(&elInput_dispdot_data[0], 0, 12U * sizeof(double));
  std::memset(&elInput_dispddot_data[0], 0, 12U * sizeof(double));
  std::memset(&eldispiter_data[0], 0, 12U * sizeof(double));

  //  not always used, but must be declared, specific to TNB and ROM
  // unpack displacement information
  // unpack mesh information
  // construct elInput for elementNumber
  b_index = 1.0;
  elInput_xloc_idx_1_tmp = static_cast<int>(elementNumber) - 1;
  if (c_strcmp(model_analysisType) || b_strcmp(model_analysisType)) {
    elInput_analysisType_size_idx_1 = 3;
    elInput_analysisType_data[0] = 'T';
    elInput_analysisType_data[1] = 'N';
    elInput_analysisType_data[2] = 'B';
  } else {
    // if(strcmp(analysisType,'S'))
    elInput_analysisType_size_idx_1 = 1;
    elInput_analysisType_data[0] = 'M';

    //  else
    //      error('analysisType not supported, choose another')
  }

  // unpack connectivity list, nodal terms, etc.
  for (j = 0; j < 2; j++) {
    d = mesh_conn->data[elInput_xloc_idx_1_tmp + mesh_conn->size[0] * j];
    elx_tmp = static_cast<int>(d) - 1;
    elx_data[j] = mesh_x->data[elx_tmp];
    ely_data[j] = mesh_y->data[elx_tmp];
    elInput_z_data[j] = mesh_z->data[elx_tmp];
    for (k = 0; k < 6; k++) {
      // get element cooridnates
      if (c_strcmp(model_analysisType) || b_strcmp(model_analysisType)) {
        aoffset = static_cast<int>(((d - 1.0) * 6.0 + (static_cast<double>(k) +
          1.0))) - 1;
        elx_tmp = static_cast<int>(b_index) - 1;
        eldisp_data[elx_tmp] = dispData.displ_s->data[aoffset];
        elInput_dispdot_data[elx_tmp] = dispData.displdot_s->data[aoffset];
        elInput_dispddot_data[elx_tmp] = dispData.displddot_s->data[aoffset];
        eldispiter_data[elx_tmp] = displ_iter->data[aoffset];
      }

      b_index++;
    }
  }

  // specific to TNB and ROM
  // specific to TNB and ROM
  if ((!c_strcmp(model_analysisType)) && (!b_strcmp(model_analysisType))) {
    // if(strcmp(analysisType,'S'))
    std::memcpy(&eldispiter_data[0], &eldisp_data[0], 12U * sizeof(double));

    //  else
    //      error('analysisType not supported, choose another')
  }

  // get concentrated terms associated with element
  b_mesh_conn[0] = mesh_conn->data[elInput_xloc_idx_1_tmp];
  b_mesh_conn[1] = mesh_conn->data[elInput_xloc_idx_1_tmp + mesh_conn->size[0]];
  ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);
  if (c_strcmp(model_analysisType) || b_strcmp(model_analysisType)) {
    j = 1;
    k = 3;
    elInput_omegaVec_size_idx_0 = 1;
    elInput_omegaVec_size_idx_1 = 3;
    aoffset = 1;
    elx_tmp = 3;
    elInput_accelVec_data[0] = rbData[0];
    elInput_omegaVec_data[0] = rbData[3];
    elInput_omegaDotVec_data[0] = rbData[6];
    elInput_accelVec_data[1] = rbData[1];
    elInput_omegaVec_data[1] = rbData[4];
    elInput_omegaDotVec_data[1] = rbData[7];
    elInput_accelVec_data[2] = rbData[2];
    elInput_omegaVec_data[2] = rbData[5];
    elInput_omegaDotVec_data[2] = rbData[8];
  } else {
    // if(strcmp(analysisType,'S'))
    j = 3;
    k = 1;
    elInput_omegaVec_size_idx_0 = 3;
    elInput_omegaVec_size_idx_1 = 1;
    aoffset = 3;
    elx_tmp = 1;
    elInput_accelVec_data[0] = 0.0;
    elInput_omegaVec_data[0] = 0.0;
    elInput_omegaDotVec_data[0] = 0.0;
    elInput_accelVec_data[1] = 0.0;
    elInput_omegaVec_data[1] = 0.0;
    elInput_omegaDotVec_data[1] = 0.0;
    elInput_accelVec_data[2] = 0.0;
    elInput_omegaVec_data[2] = 0.0;
    elInput_omegaDotVec_data[2] = 0.0;

    //  else
    //      error('analysisType not supported, choose another')
  }

  //  Specific to TNB and ROM, but must be declared
  // Is not used for this model type, but must be declared.
  // calculate element stiffness matrix and force vector
  // (or effective stiffness matrix and force vector from time integration)
  expl_temp.firstIteration = false;
  expl_temp.freq = 0.0;
  expl_temp.airDensity = 0.0;
  std::memcpy(&expl_temp.CN2H[0], &CN2H[0], 9U * sizeof(double));
  expl_temp.OmegaDot = OmegaDot;
  expl_temp.Omega = Omega;
  expl_temp.omegaDotVec.size[0] = aoffset;
  expl_temp.omegaDotVec.size[1] = elx_tmp;
  aoffset *= elx_tmp;
  if (0 <= aoffset - 1) {
    std::memcpy(&expl_temp.omegaDotVec.data[0], &elInput_omegaDotVec_data[0],
                aoffset * sizeof(double));
  }

  expl_temp.omegaVec.size[0] = elInput_omegaVec_size_idx_0;
  expl_temp.omegaVec.size[1] = elInput_omegaVec_size_idx_1;
  aoffset = elInput_omegaVec_size_idx_0 * elInput_omegaVec_size_idx_1;
  if (0 <= aoffset - 1) {
    std::memcpy(&expl_temp.omegaVec.data[0], &elInput_omegaVec_data[0], aoffset *
                sizeof(double));
  }

  expl_temp.accelVec.size[0] = j;
  expl_temp.accelVec.size[1] = k;
  aoffset = j * k;
  if (0 <= aoffset - 1) {
    std::memcpy(&expl_temp.accelVec.data[0], &elInput_accelVec_data[0], aoffset *
                sizeof(double));
  }

  expl_temp.RayleighBeta = model_RayleighBeta;
  expl_temp.RayleighAlpha = model_RayleighAlpha;
  expl_temp.gravityOn = true;
  expl_temp.z.size[0] = 2;
  expl_temp.z.size[1] = 2;
  expl_temp.z.data[0] = elInput_z_data[0];
  expl_temp.z.data[1] = elInput_z_data[1];
  expl_temp.z.data[2] = 0.0;
  expl_temp.z.data[3] = 0.0;
  expl_temp.y.size[0] = 2;
  expl_temp.x.size[0] = 2;
  expl_temp.y.data[0] = ely_data[0];
  expl_temp.x.data[0] = elx_data[0];
  expl_temp.y.data[1] = ely_data[1];
  expl_temp.x.data[1] = elx_data[1];
  expl_temp.dispm1.size[0] = 1;
  expl_temp.dispm1.size[1] = 12;
  std::memset(&expl_temp.dispm1.data[0], 0, 12U * sizeof(double));
  std::memset(&expl_temp.concLoad[0], 0, 12U * sizeof(double));
  std::memset(&expl_temp.concStiff[0], 0, 12U * sizeof(double));
  std::memcpy(&expl_temp.concMass[0], &elInput_concMass[0], 8U * sizeof(double));
  expl_temp.displ_iter.size[0] = 1;
  expl_temp.displ_iter.size[1] = 12;
  expl_temp.dispddot.size[0] = 1;
  expl_temp.dispddot.size[1] = 12;
  expl_temp.dispdot.size[0] = 1;
  expl_temp.dispdot.size[1] = 12;
  expl_temp.disp.size[0] = 1;
  expl_temp.disp.size[1] = 12;
  std::memcpy(&expl_temp.displ_iter.data[0], &eldispiter_data[0], 12U * sizeof
              (double));
  std::memcpy(&expl_temp.dispddot.data[0], &elInput_dispddot_data[0], 12U *
              sizeof(double));
  std::memcpy(&expl_temp.dispdot.data[0], &elInput_dispdot_data[0], 12U * sizeof
              (double));
  std::memcpy(&expl_temp.disp.data[0], &eldisp_data[0], 12U * sizeof(double));
  expl_temp.analysisType.size[0] = 1;
  expl_temp.analysisType.size[1] = elInput_analysisType_size_idx_1;
  if (0 <= elInput_analysisType_size_idx_1 - 1) {
    std::memcpy(&expl_temp.analysisType.data[0], &elInput_analysisType_data[0],
                elInput_analysisType_size_idx_1 * sizeof(char));
  }

  expl_temp.tolerance = 1.0E-6;
  expl_temp.MAXIT = 50.0;
  expl_temp.maxNumLoadSteps = 20.0;
  expl_temp.loadStep = 1.0;
  expl_temp.loadStepPrev = 0.0;
  expl_temp.aeroForceOn = false;
  expl_temp.aeroElasticOn = false;
  expl_temp.preStress = false;
  expl_temp.useDisp = false;
  expl_temp.aeroSweepAngle = 0.0;
  expl_temp.rollAngle = el_roll->data[elInput_xloc_idx_1_tmp];
  expl_temp.coneAngle = el_theta->data[elInput_xloc_idx_1_tmp];
  expl_temp.sweepAngle = el_psi->data[elInput_xloc_idx_1_tmp];
  expl_temp.sectionProps = el_props->data[elInput_xloc_idx_1_tmp];
  expl_temp.iterationType[0] = 'D';
  expl_temp.xloc[0] = 0.0;
  expl_temp.iterationType[1] = 'I';
  expl_temp.xloc[1] = el_elLen->data[elInput_xloc_idx_1_tmp];
  expl_temp.timeInt = *timeInt;
  expl_temp.modalFlag = true;
  expl_temp.elementOrder = 1.0;
  b_calculateTimoshenkoElementNL(&expl_temp, &elStorage->
    data[elInput_xloc_idx_1_tmp], &b_expl_temp);
  std::memcpy(&eldisp_data[0], &b_expl_temp.FhatLessConc[0], 12U * sizeof(double));
  std::memcpy(&elOutput_Ke[0], &b_expl_temp.Ke[0], 144U * sizeof(double));

  // post process for reaction force
  for (elx_tmp = 0; elx_tmp < 12; elx_tmp++) {
    elInput_dispdot_data[elx_tmp] = eldispiter_data[elx_tmp];
    Fpp[elx_tmp] = 0.0;
  }

  for (k = 0; k < 12; k++) {
    aoffset = k * 12;
    for (elx_tmp = 0; elx_tmp < 12; elx_tmp++) {
      Fpp[elx_tmp] += elInput_dispdot_data[k] * elOutput_Ke[aoffset + elx_tmp];
    }
  }

  // if(strcmp(analysisType,'TNB')||strcmp(analysisType,'ROM')||strcmp(analysisType,'S')) 
  //  else
  //      error('analysisType not supported, choose another')
  for (elx_tmp = 0; elx_tmp < 12; elx_tmp++) {
    Fpp[elx_tmp] -= eldisp_data[elx_tmp];
  }

  // ----------------------------------------------------------------------
}

//
// File trailer for elementPostProcess.cpp
//
// [EOF]
//
