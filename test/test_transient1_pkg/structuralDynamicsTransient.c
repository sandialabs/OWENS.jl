/*
 * File: structuralDynamicsTransient.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "structuralDynamicsTransient.h"
#include "ConcMassAssociatedWithElement.h"
#include "assembly.h"
#include "calculateReactionForceAtNode.h"
#include "calculateStrainForElements.h"
#include "calculateTimoshenkoElementNL.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "strcmp.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
#include <string.h>

/* Function Definitions */

/*
 * structuralDynamicsTransient perform transient analysis
 *  **********************************************************************
 *  *                   Part of the SNL OWENS Toolkit                    *
 *  * Developed by Sandia National Laboratories Wind Energy Technologies *
 *  *             See license.txt for disclaimer information             *
 *  **********************************************************************
 *    [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,...
 *                              dispData,Omega,OmegaDot,time,delta_t,...
 *                              elStorage,Fexternal,Fdof,CN2H,rbData)
 *
 *    This function performs transient structural dynamics analysis.
 *
 *    input:
 *    model      = object containing model data
 *    mesh       = object containing mesh data
 *    el         = object containing element data
 *    dispData   = object containing displacement data
 *    Omega      = rotor speed (Hz)
 *    OmegaDot   = rotor acceleratin (Hz)
 *    time       = current simulation time
 *    delta_t    = time step size
 *    elStorage  = object containing stored element data
 *    Fexternal  = vector containing external force values
 *    Fdof       = vector containing global DOF numbering associated with
 *                 external force values
 *    CN2H       = transformation matrix from inertial frame to hub frame
 *    rbData     = vector containing rigid body displacement, velocity, and
 *                 acceleration
 *
 *    output:
 *    dispOut       = object containing displacement data at end of time step
 *    FReaction_sp1 = vector containing reaction force at turbine base at
 *                    end of time step
 * Arguments    : const char model_analysisType[3]
 *                double model_RayleighAlpha
 *                double model_RayleighBeta
 *                double model_BC_numpBC
 *                const emxArray_real_T *model_BC_pBC
 *                const emxArray_real_T *model_joint
 *                const emxArray_real_T *model_jointTransform
 *                double mesh_numEl
 *                const emxArray_real_T *mesh_x
 *                const emxArray_real_T *mesh_y
 *                const emxArray_real_T *mesh_z
 *                const emxArray_real_T *mesh_conn
 *                const f_struct_T el
 *                const g_struct_T dispData
 *                double Omega
 *                double OmegaDot
 *                const c_emxArray_struct_T *elStorage
 *                const emxArray_real_T *Fexternal
 *                const emxArray_real_T *Fdof
 *                const double CN2H[9]
 *                h_struct_T *dispOut
 *                double FReaction_sp1[6]
 * Return Type  : void
 */
void structuralDynamicsTransient(const char model_analysisType[3], double
  model_RayleighAlpha, double model_RayleighBeta, double model_BC_numpBC, const
  emxArray_real_T *model_BC_pBC, const emxArray_real_T *model_joint, const
  emxArray_real_T *model_jointTransform, double mesh_numEl, const
  emxArray_real_T *mesh_x, const emxArray_real_T *mesh_y, const emxArray_real_T *
  mesh_z, const emxArray_real_T *mesh_conn, const f_struct_T el, const
  g_struct_T dispData, double Omega, double OmegaDot, const c_emxArray_struct_T *
  elStorage, const emxArray_real_T *Fexternal, const emxArray_real_T *Fdof,
  const double CN2H[9], h_struct_T *dispOut, double FReaction_sp1[6])
{
  double totalNumDOF;
  double elz_data[4];
  double eldisp_data[12];
  double eldispdot_data[12];
  double eldispddot_data[12];
  double eldispiter_data[12];
  emxArray_real_T *disp_s;
  emxArray_real_T *dispdot_s;
  emxArray_real_T *dispddot_s;
  emxArray_real_T *solution;
  int i;
  emxArray_real_T *Kg;
  emxArray_real_T *Fg;
  emxArray_real_T *y;
  int loop_ub;
  emxArray_real_T *a;
  int m;
  boolean_T elInput_firstIteration;
  int b_i;
  int j;
  double b_mesh_conn[2];
  double elInput_concMass[8];
  double elx_data[2];
  double ely_data[2];
  l_struct_T expl_temp;
  int k;
  int boffset;
  int inner;
  int coffset;
  int aoffset;
  m_struct_T b_expl_temp;
  double elOutput_Ke[144];
  double elOutput_Fe[12];
  double eqNumber;
  struct_T c_expl_temp;

  /* -------- get model information ----------- */
  totalNumDOF = (double)mesh_x->size[0] * 6.0;

  /*  [~,numReducedDOF]=size(model.jointTransform); */
  /* ----------------------------------------- */
  /* initialize displacements, tolerance, uNorm, iteration count for nonlinear */
  /* iteration */
  elz_data[0] = 0.0;
  elz_data[1] = 0.0;
  memset(&eldisp_data[0], 0, 12U * sizeof(double));
  memset(&eldispdot_data[0], 0, 12U * sizeof(double));
  memset(&eldispddot_data[0], 0, 12U * sizeof(double));
  memset(&eldispiter_data[0], 0, 12U * sizeof(double));

  /*  if(model.nlOn) */
  /*       iterationType = model.nlParams.iterationType; */
  /*  else */
  /*       iterationType = 'LINEAR'; */
  /*  end */
  emxInit_real_T(&disp_s, 1);
  emxInit_real_T(&dispdot_s, 1);
  emxInit_real_T(&dispddot_s, 1);
  emxInit_real_T(&solution, 1);
  if (c_strcmp(model_analysisType)) {
    /* ------ newmark integration parameters --------- */
    i = disp_s->size[0];
    disp_s->size[0] = dispData.displ_s->size[0];
    emxEnsureCapacity_real_T(disp_s, i);
    loop_ub = dispData.displ_s->size[0];
    for (i = 0; i < loop_ub; i++) {
      disp_s->data[i] = dispData.displ_s->data[i];
    }

    i = dispdot_s->size[0];
    dispdot_s->size[0] = dispData.displdot_s->size[0];
    emxEnsureCapacity_real_T(dispdot_s, i);
    loop_ub = dispData.displdot_s->size[0];
    for (i = 0; i < loop_ub; i++) {
      dispdot_s->data[i] = dispData.displdot_s->data[i];
    }

    i = dispddot_s->size[0];
    dispddot_s->size[0] = dispData.displddot_s->size[0];
    emxEnsureCapacity_real_T(dispddot_s, i);
    loop_ub = dispData.displddot_s->size[0];
    for (i = 0; i < loop_ub; i++) {
      dispddot_s->data[i] = dispData.displddot_s->data[i];
    }

    i = solution->size[0];
    solution->size[0] = dispData.displ_s->size[0];
    emxEnsureCapacity_real_T(solution, i);
    loop_ub = dispData.displ_s->size[0];
    for (i = 0; i < loop_ub; i++) {
      solution->data[i] = dispData.displ_s->data[i];
    }
  }

  /* ----------------------------------------------- */
  emxInit_real_T(&Kg, 2);
  emxInit_real_T(&Fg, 1);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&a, 2);

  /* iteration loop */
  /* ------- intitialization ----------------- */
  i = Kg->size[0] * Kg->size[1];
  loop_ub = (int)totalNumDOF;
  Kg->size[0] = loop_ub;
  Kg->size[1] = loop_ub;
  emxEnsureCapacity_real_T(Kg, i);
  m = loop_ub * loop_ub;
  for (i = 0; i < m; i++) {
    Kg->data[i] = 0.0;
  }

  /* initialize global stiffness and force vector */
  i = Fg->size[0];
  Fg->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Fg, i);
  for (i = 0; i < loop_ub; i++) {
    Fg->data[i] = 0.0;
  }

  /* ------------------------------------------- */
  /* ---- element  calculation and assembly ---------------------------------- */
  i = (int)mesh_numEl;
  if (0 <= i - 1) {
    elInput_firstIteration = true;
  }

  for (b_i = 0; b_i < i; b_i++) {
    /* Calculate Ke and Fe for element i */
    totalNumDOF = 1.0;

    /* initialize element data */
    /* get concentrated terms associated with elemetn */
    for (j = 0; j < 2; j++) {
      /* get element cooridnates */
      loop_ub = (int)mesh_conn->data[b_i + mesh_conn->size[0] * j] - 1;
      elx_data[j] = mesh_x->data[loop_ub];
      ely_data[j] = mesh_y->data[loop_ub];
      elz_data[j] = mesh_z->data[loop_ub];

      /* get element nodal displacements at s and s-1 time step */
      for (k = 0; k < 6; k++) {
        /*                  if(strcmp(analysisType,'TD')) */
        /*                      eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k); */
        /*                      eldisp_sm1(index) = disp_sm1((conn(i,j)-1)*numDOFPerNode + k); */
        /*                      eldispiter(index) = displ_iter((conn(i,j)-1)*numDOFPerNode + k); */
        /*                  end */
        if (c_strcmp(model_analysisType)) {
          loop_ub = (int)((mesh_conn->data[b_i + mesh_conn->size[0] * j] - 1.0) *
                          6.0 + ((double)k + 1.0)) - 1;
          m = (int)totalNumDOF - 1;
          eldispiter_data[m] = solution->data[loop_ub];
          eldisp_data[m] = disp_s->data[loop_ub];
          eldispdot_data[m] = dispdot_s->data[loop_ub];
          eldispddot_data[m] = dispddot_s->data[loop_ub];
        }

        totalNumDOF++;
      }

      b_mesh_conn[j] = mesh_conn->data[b_i + mesh_conn->size[0] * j];
    }

    ConcMassAssociatedWithElement(b_mesh_conn, model_joint, elInput_concMass);

    /*  specific to 'TD', but must be declared */
    /*  specific to 'TNB' , but must be declared */
    /* Is not used for this model type, but must be declared. */
    expl_temp.freq = 0.0;
    expl_temp.airDensity = 0.0;
    memcpy(&expl_temp.CN2H[0], &CN2H[0], 9U * sizeof(double));
    expl_temp.RayleighBeta = model_RayleighBeta;
    expl_temp.RayleighAlpha = model_RayleighAlpha;
    expl_temp.gravityOn = true;
    expl_temp.aeroForceOn = false;
    expl_temp.aeroElasticOn = false;
    expl_temp.iterationType[0] = 'D';
    expl_temp.iterationType[1] = 'I';
    expl_temp.preStress = false;
    expl_temp.useDisp = false;
    expl_temp.displ_iter.size[0] = 1;
    expl_temp.displ_iter.size[1] = 12;
    expl_temp.OmegaDot = OmegaDot;
    expl_temp.Omega = Omega;
    expl_temp.omegaDotVec[0] = 0.0;
    expl_temp.omegaVec[0] = 0.0;
    expl_temp.accelVec[0] = 0.0;
    expl_temp.omegaDotVec[1] = 0.0;
    expl_temp.omegaVec[1] = 0.0;
    expl_temp.accelVec[1] = 0.0;
    expl_temp.omegaDotVec[2] = 0.0;
    expl_temp.omegaVec[2] = 0.0;
    expl_temp.accelVec[2] = 0.0;
    expl_temp.z.size[0] = 2;
    expl_temp.z.size[1] = 2;
    expl_temp.z.data[0] = elz_data[0];
    expl_temp.z.data[1] = elz_data[1];
    expl_temp.z.data[2] = 0.0;
    expl_temp.z.data[3] = 0.0;
    expl_temp.y.size[0] = 2;
    expl_temp.x.size[0] = 2;
    expl_temp.y.data[0] = ely_data[0];
    expl_temp.x.data[0] = elx_data[0];
    expl_temp.y.data[1] = ely_data[1];
    expl_temp.x.data[1] = elx_data[1];
    expl_temp.dispddot.size[0] = 1;
    expl_temp.dispddot.size[1] = 12;
    expl_temp.dispdot.size[0] = 1;
    expl_temp.dispdot.size[1] = 12;
    expl_temp.dispm1.size[0] = 1;
    expl_temp.dispm1.size[1] = 12;
    expl_temp.disp.size[0] = 1;
    expl_temp.disp.size[1] = 12;
    for (boffset = 0; boffset < 12; boffset++) {
      expl_temp.displ_iter.data[boffset] = eldispiter_data[boffset];
      expl_temp.dispddot.data[boffset] = eldispddot_data[boffset];
      expl_temp.dispdot.data[boffset] = eldispdot_data[boffset];
      expl_temp.dispm1.data[boffset] = 0.0;
      expl_temp.disp.data[boffset] = eldisp_data[boffset];
      expl_temp.concLoad[boffset] = 0.0;
      expl_temp.concStiff[boffset] = 0.0;
    }

    memcpy(&expl_temp.concMass[0], &elInput_concMass[0], 8U * sizeof(double));
    expl_temp.firstIteration = elInput_firstIteration;
    expl_temp.aeroSweepAngle = 0.0;
    expl_temp.rollAngle = el.roll->data[b_i];
    expl_temp.coneAngle = el.theta->data[b_i];
    expl_temp.sweepAngle = el.psi->data[b_i];
    expl_temp.sectionProps = el.props->data[b_i];
    expl_temp.xloc[0] = 0.0;
    expl_temp.xloc[1] = el.elLen->data[b_i];
    expl_temp.timeInt.delta_t = 0.002;
    expl_temp.timeInt.a1 = 0.001;
    expl_temp.timeInt.a2 = 0.001;
    expl_temp.timeInt.a3 = 1.0E+6;
    expl_temp.timeInt.a4 = 2000.0;
    expl_temp.timeInt.a5 = 1.0;
    expl_temp.timeInt.a6 = 1000.0;
    expl_temp.timeInt.a7 = 1.0;
    expl_temp.timeInt.a8 = 0.0;
    expl_temp.modalFlag = true;
    expl_temp.elementOrder = 1.0;
    expl_temp.analysisType[0] = model_analysisType[0];
    expl_temp.analysisType[1] = model_analysisType[1];
    expl_temp.analysisType[2] = model_analysisType[2];
    calculateTimoshenkoElementNL(&expl_temp, &elStorage->data[b_i], &b_expl_temp);
    memcpy(&elOutput_Ke[0], &b_expl_temp.Ke[0], 144U * sizeof(double));
    memcpy(&elOutput_Fe[0], &b_expl_temp.Fe[0], 12U * sizeof(double));

    /* calculate timoshenko element */
    b_mesh_conn[0] = mesh_conn->data[b_i];
    b_mesh_conn[1] = mesh_conn->data[b_i + mesh_conn->size[0]];
    assembly(elOutput_Ke, elOutput_Fe, b_mesh_conn, Kg, Fg);

    /* assemble element stiffness matrix and force vector */
    /*          Erestotal = Erestotal + elOutput.Eres; */
    /* ................................................ */
  }

  /* ------- end element calculation and assembly ------------------ */
  /*     %% */
  /* ---------------------------------------------------------------------- */
  /*     %% */
  /* Apply external loads to structure */
  if ((Fexternal->size[0] == 0) || (Fexternal->size[1] == 0)) {
    loop_ub = 0;
  } else {
    loop_ub = Fexternal->size[1];
  }

  for (b_i = 0; b_i < loop_ub; b_i++) {
    if (c_strcmp(model_analysisType)) {
      Fg->data[(int)Fdof->data[b_i] - 1] += Fexternal->data[b_i];
    }
  }

  /* ------ apply constraints on system ----------------------------------- */
  /* This function transforms a matrix by the transformation matrix to */
  /* enforce joint constraints */
  i = a->size[0] * a->size[1];
  a->size[0] = model_jointTransform->size[1];
  a->size[1] = model_jointTransform->size[0];
  emxEnsureCapacity_real_T(a, i);
  loop_ub = model_jointTransform->size[0];
  for (i = 0; i < loop_ub; i++) {
    m = model_jointTransform->size[1];
    for (boffset = 0; boffset < m; boffset++) {
      a->data[boffset + a->size[0] * i] = model_jointTransform->data[i +
        model_jointTransform->size[0] * boffset];
    }
  }

  if ((a->size[1] == 1) || (Kg->size[0] == 1)) {
    i = y->size[0] * y->size[1];
    y->size[0] = a->size[0];
    y->size[1] = Kg->size[1];
    emxEnsureCapacity_real_T(y, i);
    loop_ub = a->size[0];
    for (i = 0; i < loop_ub; i++) {
      m = Kg->size[1];
      for (boffset = 0; boffset < m; boffset++) {
        y->data[i + y->size[0] * boffset] = 0.0;
        inner = a->size[1];
        for (coffset = 0; coffset < inner; coffset++) {
          y->data[i + y->size[0] * boffset] += a->data[i + a->size[0] * coffset]
            * Kg->data[coffset + Kg->size[0] * boffset];
        }
      }
    }
  } else {
    m = a->size[0];
    inner = a->size[1];
    loop_ub = Kg->size[1];
    i = y->size[0] * y->size[1];
    y->size[0] = a->size[0];
    y->size[1] = Kg->size[1];
    emxEnsureCapacity_real_T(y, i);
    for (j = 0; j < loop_ub; j++) {
      coffset = j * m;
      boffset = j * inner;
      for (b_i = 0; b_i < m; b_i++) {
        y->data[coffset + b_i] = 0.0;
      }

      for (k = 0; k < inner; k++) {
        aoffset = k * m;
        totalNumDOF = Kg->data[boffset + k];
        for (b_i = 0; b_i < m; b_i++) {
          i = coffset + b_i;
          y->data[i] += totalNumDOF * a->data[aoffset + b_i];
        }
      }
    }
  }

  if ((y->size[1] == 1) || (model_jointTransform->size[0] == 1)) {
    i = Kg->size[0] * Kg->size[1];
    Kg->size[0] = y->size[0];
    Kg->size[1] = model_jointTransform->size[1];
    emxEnsureCapacity_real_T(Kg, i);
    loop_ub = y->size[0];
    for (i = 0; i < loop_ub; i++) {
      m = model_jointTransform->size[1];
      for (boffset = 0; boffset < m; boffset++) {
        Kg->data[i + Kg->size[0] * boffset] = 0.0;
        inner = y->size[1];
        for (coffset = 0; coffset < inner; coffset++) {
          Kg->data[i + Kg->size[0] * boffset] += y->data[i + y->size[0] *
            coffset] * model_jointTransform->data[coffset +
            model_jointTransform->size[0] * boffset];
        }
      }
    }
  } else {
    m = y->size[0];
    inner = y->size[1];
    loop_ub = model_jointTransform->size[1];
    i = Kg->size[0] * Kg->size[1];
    Kg->size[0] = y->size[0];
    Kg->size[1] = model_jointTransform->size[1];
    emxEnsureCapacity_real_T(Kg, i);
    for (j = 0; j < loop_ub; j++) {
      coffset = j * m;
      boffset = j * inner;
      for (b_i = 0; b_i < m; b_i++) {
        Kg->data[coffset + b_i] = 0.0;
      }

      for (k = 0; k < inner; k++) {
        aoffset = k * m;
        totalNumDOF = model_jointTransform->data[boffset + k];
        for (b_i = 0; b_i < m; b_i++) {
          i = coffset + b_i;
          Kg->data[i] += totalNumDOF * y->data[aoffset + b_i];
        }
      }
    }
  }

  /* This function transforms a vector by the transformation matrix to */
  /* enforce joint constraints */
  i = a->size[0] * a->size[1];
  a->size[0] = model_jointTransform->size[1];
  a->size[1] = model_jointTransform->size[0];
  emxEnsureCapacity_real_T(a, i);
  loop_ub = model_jointTransform->size[0];
  for (i = 0; i < loop_ub; i++) {
    m = model_jointTransform->size[1];
    for (boffset = 0; boffset < m; boffset++) {
      a->data[boffset + a->size[0] * i] = model_jointTransform->data[i +
        model_jointTransform->size[0] * boffset];
    }
  }

  i = solution->size[0];
  solution->size[0] = Fg->size[0];
  emxEnsureCapacity_real_T(solution, i);
  loop_ub = Fg->size[0];
  for (i = 0; i < loop_ub; i++) {
    solution->data[i] = Fg->data[i];
  }

  if ((a->size[1] == 1) || (Fg->size[0] == 1)) {
    i = solution->size[0];
    solution->size[0] = a->size[0];
    emxEnsureCapacity_real_T(solution, i);
    loop_ub = a->size[0];
    for (i = 0; i < loop_ub; i++) {
      solution->data[i] = 0.0;
      m = a->size[1];
      for (boffset = 0; boffset < m; boffset++) {
        solution->data[i] += a->data[i + a->size[0] * boffset] * Fg->
          data[boffset];
      }
    }

    i = Fg->size[0];
    Fg->size[0] = solution->size[0];
    emxEnsureCapacity_real_T(Fg, i);
    loop_ub = solution->size[0];
    for (i = 0; i < loop_ub; i++) {
      Fg->data[i] = solution->data[i];
    }
  } else {
    m = a->size[0];
    inner = a->size[1];
    i = Fg->size[0];
    Fg->size[0] = a->size[0];
    emxEnsureCapacity_real_T(Fg, i);
    for (b_i = 0; b_i < m; b_i++) {
      Fg->data[b_i] = 0.0;
    }

    for (k = 0; k < inner; k++) {
      aoffset = k * m;
      for (b_i = 0; b_i < m; b_i++) {
        Fg->data[b_i] += solution->data[k] * a->data[aoffset + b_i];
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /*     %% */
  /* Apply BCs to global system */
  /* applyBC Applies boundary conditions to system for static analysis */
  /*  ********************************************************************** */
  /*  *                   Part of the SNL OWENS Toolkit                    * */
  /*  * Developed by Sandia National Laboratories Wind Energy Technologies * */
  /*  *             See license.txt for disclaimer information             * */
  /*  ********************************************************************** */
  /*    [Kg,Fg] = applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode) */
  /*  */
  /*    This function applies boundary conditions to the stiffness matrix and */
  /*    load vector for a static analysis. */
  /*  */
  /*       input: */
  /*       Kg            = assembled global stiffness matrix */
  /*       Fg            = assembled global load vector */
  /*       BC            = struct of boundary condition information */
  /*       u             = global displacement vector */
  /*       iterationType = for nonlinear analysis, not used in BLAST */
  /*       numDofPerNode = number of degrees of freedom per node */
  /*       output: */
  /*       Kg            = global stiffness matrix with boundary conditions */
  /*       Fg            = global load vector with boundary condition */
  m = Kg->size[0] - 1;

  /* APPLY BCs FOR PRIMARY VARIABLE */
  if (model_BC_numpBC > 0.0) {
    i = model_BC_pBC->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      totalNumDOF = model_BC_pBC->data[b_i + model_BC_pBC->size[0] * 2];
      eqNumber = (model_BC_pBC->data[b_i] - 1.0) * 6.0 + model_BC_pBC->data[b_i
        + model_BC_pBC->size[0]];
      for (j = 0; j <= m; j++) {
        loop_ub = (int)eqNumber - 1;
        Kg->data[loop_ub + Kg->size[0] * j] = 0.0;
        Fg->data[j] -= Kg->data[j + Kg->size[0] * loop_ub] * totalNumDOF;
        Kg->data[j + Kg->size[0] * loop_ub] = 0.0;
      }

      loop_ub = (int)eqNumber - 1;
      Fg->data[loop_ub] = totalNumDOF;
      Kg->data[loop_ub + Kg->size[0] * loop_ub] = 1.0;
    }
  }

  /* APPLY BCs FOR SECONDARY VARIABLE */
  /*  if(BC.numsBC > 0) % This does not appear to be used */
  /*      sBC = BC.sBC; */
  /*      [numsBC,~] = size(sBC); */
  /*       */
  /*      for i=1:numsBC */
  /*          nodeNumber = sBC(i,1); */
  /*          dofNumber = sBC(i,2); */
  /*          specVal =  sBC(i,3); */
  /*           */
  /*          eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber; */
  /*           */
  /*          Fg(eqNumber) = Fg(eqNumber) + specVal; */
  /*           */
  /*      end */
  /*  end */
  mldiv(Kg, Fg);

  /* solve for displacements */
  if ((model_jointTransform->size[1] == 1) || (Fg->size[0] == 1)) {
    i = solution->size[0];
    solution->size[0] = model_jointTransform->size[0];
    emxEnsureCapacity_real_T(solution, i);
    loop_ub = model_jointTransform->size[0];
    for (i = 0; i < loop_ub; i++) {
      solution->data[i] = 0.0;
      m = model_jointTransform->size[1];
      for (boffset = 0; boffset < m; boffset++) {
        solution->data[i] += model_jointTransform->data[i +
          model_jointTransform->size[0] * boffset] * Fg->data[boffset];
      }
    }
  } else {
    m = model_jointTransform->size[0];
    inner = model_jointTransform->size[1];
    i = solution->size[0];
    solution->size[0] = model_jointTransform->size[0];
    emxEnsureCapacity_real_T(solution, i);
    for (b_i = 0; b_i < m; b_i++) {
      solution->data[b_i] = 0.0;
    }

    for (k = 0; k < inner; k++) {
      aoffset = k * m;
      for (b_i = 0; b_i < m; b_i++) {
        solution->data[b_i] += Fg->data[k] * model_jointTransform->data[aoffset
          + b_i];
      }
    }
  }

  /* transform to full dof listing */
  emxFree_real_T(&a);
  emxFree_real_T(&y);
  emxFree_real_T(&Fg);
  emxFree_real_T(&Kg);
  emxFree_real_T(&dispddot_s);
  emxFree_real_T(&dispdot_s);
  emxFree_real_T(&disp_s);

  /* Calculate reaction at turbine base (hardwired to node number 1) */
  c_expl_temp.a8 = 0.0;
  c_expl_temp.a7 = 1.0;
  c_expl_temp.a6 = 1000.0;
  c_expl_temp.a5 = 1.0;
  c_expl_temp.a4 = 2000.0;
  c_expl_temp.a3 = 1.0E+6;
  c_expl_temp.a2 = 0.001;
  c_expl_temp.a1 = 0.001;
  c_expl_temp.delta_t = 0.002;
  calculateReactionForceAtNode(model_analysisType, model_RayleighAlpha,
    model_RayleighBeta, model_joint, mesh_x, mesh_y, mesh_z, mesh_conn, el.props,
    el.elLen, el.psi, el.theta, el.roll, elStorage, &c_expl_temp, dispData,
    solution, Omega, OmegaDot, CN2H, FReaction_sp1);

  /* Calculate strain */
  calculateStrainForElements(mesh_numEl, mesh_conn, el.props, el.elLen, el.psi,
    el.theta, el.roll, solution, dispOut->elStrain);
  i = dispOut->displ_sp1->size[0];
  dispOut->displ_sp1->size[0] = solution->size[0];
  emxEnsureCapacity_real_T(dispOut->displ_sp1, i);
  loop_ub = solution->size[0];
  for (i = 0; i < loop_ub; i++) {
    dispOut->displ_sp1->data[i] = solution->data[i];
  }

  /* store displacement vector in dispOut */
  /*  Specific to TNB, but must be declared */
  loop_ub = solution->size[0];
  for (i = 0; i < loop_ub; i++) {
    solution->data[i] = (1.0E+6 * (solution->data[i] - dispData.displ_s->data[i])
                         - 2000.0 * dispData.displdot_s->data[i]) -
      dispData.displddot_s->data[i];
  }

  /* store velocity vector in dispOut */
  i = dispOut->displddot_sp1->size[0];
  dispOut->displddot_sp1->size[0] = solution->size[0];
  emxEnsureCapacity_real_T(dispOut->displddot_sp1, i);
  loop_ub = solution->size[0];
  for (i = 0; i < loop_ub; i++) {
    dispOut->displddot_sp1->data[i] = solution->data[i];
  }

  i = dispOut->displdot_sp1->size[0];
  dispOut->displdot_sp1->size[0] = dispData.displdot_s->size[0];
  emxEnsureCapacity_real_T(dispOut->displdot_sp1, i);
  loop_ub = dispData.displdot_s->size[0];
  for (i = 0; i < loop_ub; i++) {
    dispOut->displdot_sp1->data[i] = (dispData.displdot_s->data[i] + 0.001 *
      dispData.displddot_s->data[i]) + 0.001 * solution->data[i];
  }

  emxFree_real_T(&solution);

  /* store acceleration vector in dispOut */
}

/*
 * File trailer for structuralDynamicsTransient.c
 *
 * [EOF]
 */
