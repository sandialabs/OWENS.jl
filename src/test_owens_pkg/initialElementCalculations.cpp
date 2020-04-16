//
// File: initialElementCalculations.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:31:03
//

// Include Files
#include "initialElementCalculations.h"
#include "ConcMassAssociatedWithElement.h"
#include "calculateTimoshenkoElementInitialRun.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_data.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// initialElementCalculations  performs intitial element calculations
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [elStorage] = initialElementCalculations(model,el,mesh)
//
//    This function performs initial element calculation for use later in
//    analysis for efficiency gains.
//
//       input:
//       model               = object containing model information
//       el                  = object containing element information
//       mesh                = object containing mesh information
//
//       output:
//       elStorage           = object containing stored element data
// Arguments    : const emxArray_real_T *model_joint
//                const emxArray_struct_T *el_props
//                const emxArray_real_T *el_elLen
//                const emxArray_real_T *el_psi
//                const emxArray_real_T *el_theta
//                const emxArray_real_T *el_roll
//                double mesh_numEl
//                const emxArray_real_T *mesh_x
//                const emxArray_real_T *mesh_y
//                const emxArray_real_T *mesh_z
//                const emxArray_real_T *mesh_conn
//                c_emxArray_struct_T *elStorage
// Return Type  : void
//
void initialElementCalculations(const emxArray_real_T *model_joint, const
  emxArray_struct_T *el_props, const emxArray_real_T *el_elLen, const
  emxArray_real_T *el_psi, const emxArray_real_T *el_theta, const
  emxArray_real_T *el_roll, double mesh_numEl, const emxArray_real_T *mesh_x,
  const emxArray_real_T *mesh_y, const emxArray_real_T *mesh_z, const
  emxArray_real_T *mesh_conn, c_emxArray_struct_T *elStorage)
{
  int idx;
  int loop_ub;
  int i;
  double elInput_xloc[2];
  double elInput_sectionProps_twist[2];
  double elInput_sectionProps_rhoA[2];
  double elInput_sectionProps_EIyy[2];
  double elInput_sectionProps_EIzz[2];
  double elInput_sectionProps_GJ[2];
  double elInput_sectionProps_EA[2];
  double elInput_sectionProps_rhoIyy[2];
  double elInput_sectionProps_rhoIzz[2];
  double elInput_sectionProps_rhoJ[2];
  double elInput_sectionProps_zcm[2];
  double elInput_sectionProps_ycm[2];
  double d;
  double elx[2];
  double ely[2];
  double elz[2];
  double b_mesh_conn[2];
  double massConc[8];
  int ii_size_idx_0;
  int ii;
  boolean_T exitg1;

  // initial element calculation
  idx = elStorage->size[0] * elStorage->size[1];
  elStorage->size[0] = 1;
  loop_ub = static_cast<int>(mesh_numEl);
  elStorage->size[1] = loop_ub;
  emxEnsureCapacity_struct_T5(elStorage, idx);
  for (idx = 0; idx < loop_ub; idx++) {
    elStorage->data[idx] = r2;
  }

  for (i = 0; i < loop_ub; i++) {
    // Calculate Ke and Fe for element i
    // assign elInput for element i
    elInput_xloc[0] = 0.0;
    elInput_xloc[1] = el_elLen->data[i];

    // get concentrated terms associated with elemetn
    elInput_sectionProps_twist[0] = el_props->data[i].twist[0];
    elInput_sectionProps_rhoA[0] = el_props->data[i].rhoA[0];
    elInput_sectionProps_EIyy[0] = el_props->data[i].EIyy[0];
    elInput_sectionProps_EIzz[0] = el_props->data[i].EIzz[0];
    elInput_sectionProps_GJ[0] = el_props->data[i].GJ[0];
    elInput_sectionProps_EA[0] = el_props->data[i].EA[0];
    elInput_sectionProps_rhoIyy[0] = el_props->data[i].rhoIyy[0];
    elInput_sectionProps_rhoIzz[0] = el_props->data[i].rhoIzz[0];
    elInput_sectionProps_rhoJ[0] = el_props->data[i].rhoJ[0];
    elInput_sectionProps_zcm[0] = el_props->data[i].zcm[0];
    elInput_sectionProps_ycm[0] = el_props->data[i].ycm[0];

    // get element cooridnates
    d = mesh_conn->data[i];
    idx = static_cast<int>(d) - 1;
    elx[0] = mesh_x->data[idx];
    ely[0] = mesh_y->data[idx];
    elz[0] = mesh_z->data[idx];
    b_mesh_conn[0] = d;
    elInput_sectionProps_twist[1] = el_props->data[i].twist[1];
    elInput_sectionProps_rhoA[1] = el_props->data[i].rhoA[1];
    elInput_sectionProps_EIyy[1] = el_props->data[i].EIyy[1];
    elInput_sectionProps_EIzz[1] = el_props->data[i].EIzz[1];
    elInput_sectionProps_GJ[1] = el_props->data[i].GJ[1];
    elInput_sectionProps_EA[1] = el_props->data[i].EA[1];
    elInput_sectionProps_rhoIyy[1] = el_props->data[i].rhoIyy[1];
    elInput_sectionProps_rhoIzz[1] = el_props->data[i].rhoIzz[1];
    elInput_sectionProps_rhoJ[1] = el_props->data[i].rhoJ[1];
    elInput_sectionProps_zcm[1] = el_props->data[i].zcm[1];
    elInput_sectionProps_ycm[1] = el_props->data[i].ycm[1];

    // get element cooridnates
    d = mesh_conn->data[i + mesh_conn->size[0]];
    idx = static_cast<int>(d) - 1;
    elx[1] = mesh_x->data[idx];
    ely[1] = mesh_y->data[idx];
    elz[1] = mesh_z->data[idx];
    b_mesh_conn[1] = d;
    ConcMassAssociatedWithElement(b_mesh_conn, model_joint, massConc);
    idx = 0;
    ii_size_idx_0 = 1;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 8)) {
      if (massConc[ii] != 0.0) {
        idx = 1;
        exitg1 = true;
      } else {
        ii++;
      }
    }

    if (idx == 0) {
      ii_size_idx_0 = 0;
    }

    // only needed for structure mass props (not used in saved element matrices) 
    c_calculateTimoshenkoElementIni(elInput_xloc, elInput_sectionProps_twist,
      elInput_sectionProps_rhoA, elInput_sectionProps_EIyy,
      elInput_sectionProps_EIzz, elInput_sectionProps_GJ,
      elInput_sectionProps_EA, elInput_sectionProps_rhoIyy,
      elInput_sectionProps_rhoIzz, elInput_sectionProps_rhoJ,
      elInput_sectionProps_zcm, elInput_sectionProps_ycm, el_psi->data[i],
      el_theta->data[i], el_roll->data[i], elx, ely, elz, ii_size_idx_0 != 0,
      massConc, &elStorage->data[i]);

    // initial element calculations for storage
  }
}

//
// File trailer for initialElementCalculations.cpp
//
// [EOF]
//
