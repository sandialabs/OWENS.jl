//
// File: fminbnd.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 10:41:25
//

// Include Files
#include "fminbnd.h"
#include "linearAnalysisModal.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <string.h>

// Function Definitions

//
// Arguments    : double c_funfcnInput_tunableEnvironmen
//                const char d_funfcnInput_tunableEnvironmen[2]
//                boolean_T e_funfcnInput_tunableEnvironmen
//                double f_funfcnInput_tunableEnvironmen
//                double g_funfcnInput_tunableEnvironmen
//                double h_funfcnInput_tunableEnvironmen
//                const emxArray_real_T *i_funfcnInput_tunableEnvironmen
//                const emxArray_real_T *j_funfcnInput_tunableEnvironmen
//                const emxArray_real_T *k_funfcnInput_tunableEnvironmen
//                const char l_funfcnInput_tunableEnvironmen[]
//                const int m_funfcnInput_tunableEnvironmen[2]
//                const emxArray_real_T *n_funfcnInput_tunableEnvironmen
//                const emxArray_real_T *o_funfcnInput_tunableEnvironmen
//                const h_struct_T p_funfcnInput_tunableEnvironmen
//                const i_struct_T q_funfcnInput_tunableEnvironmen
//                const emxArray_real_T *r_funfcnInput_tunableEnvironmen
//                const c_emxArray_struct_T *s_funfcnInput_tunableEnvironmen
//                double bx
//                double *xf
//                double *fval
//                double *exitflag
//                double *output_iterations
//                double *output_funcCount
// Return Type  : void
//
void fminbnd(double c_funfcnInput_tunableEnvironmen, const char
             d_funfcnInput_tunableEnvironmen[2], boolean_T
             e_funfcnInput_tunableEnvironmen, double
             f_funfcnInput_tunableEnvironmen, double
             g_funfcnInput_tunableEnvironmen, double
             h_funfcnInput_tunableEnvironmen, const emxArray_real_T
             *i_funfcnInput_tunableEnvironmen, const emxArray_real_T
             *j_funfcnInput_tunableEnvironmen, const emxArray_real_T
             *k_funfcnInput_tunableEnvironmen, const char
             l_funfcnInput_tunableEnvironmen[], const int
             m_funfcnInput_tunableEnvironmen[2], const emxArray_real_T
             *n_funfcnInput_tunableEnvironmen, const emxArray_real_T
             *o_funfcnInput_tunableEnvironmen, const h_struct_T
             p_funfcnInput_tunableEnvironmen, const i_struct_T
             q_funfcnInput_tunableEnvironmen, const emxArray_real_T
             *r_funfcnInput_tunableEnvironmen, const c_emxArray_struct_T
             *s_funfcnInput_tunableEnvironmen, double bx, double *xf, double
             *fval, double *exitflag, double *output_iterations, double
             *output_funcCount)
{
  emxArray_real_T *freq;
  emxArray_real_T *unusedU8;
  emxArray_real_T *unusedU9;
  emxArray_real_T *unusedUa;
  emxArray_real_T *unusedUb;
  int iter;
  double a;
  double b;
  double v;
  double w;
  double d;
  double e;
  int fx_tmp;
  int funccount;
  double fv;
  double fw;
  double xm;
  double tol1;
  double tol2;
  boolean_T exitg1;
  boolean_T guard1 = false;
  double p;
  double b_r;
  double x;
  double q;
  emxInit_real_T(&freq, 2);
  emxInit_real_T(&unusedU8, 2);
  emxInit_real_T(&unusedU9, 3);
  emxInit_real_T(&unusedUa, 3);
  emxInit_real_T(&unusedUb, 2);
  iter = 0;
  *exitflag = 1.0;
  a = 0.0;
  b = bx;
  v = 0.3819660112501051 * bx;
  w = v;
  *xf = v;
  d = 0.0;
  e = 0.0;
  b_linearAnalysisModal(d_funfcnInput_tunableEnvironmen,
                        e_funfcnInput_tunableEnvironmen, v,
                        f_funfcnInput_tunableEnvironmen,
                        g_funfcnInput_tunableEnvironmen,
                        h_funfcnInput_tunableEnvironmen,
                        i_funfcnInput_tunableEnvironmen,
                        j_funfcnInput_tunableEnvironmen,
                        k_funfcnInput_tunableEnvironmen,
                        l_funfcnInput_tunableEnvironmen,
                        m_funfcnInput_tunableEnvironmen,
                        n_funfcnInput_tunableEnvironmen,
                        o_funfcnInput_tunableEnvironmen,
                        p_funfcnInput_tunableEnvironmen.numEl,
                        p_funfcnInput_tunableEnvironmen.x,
                        p_funfcnInput_tunableEnvironmen.y,
                        p_funfcnInput_tunableEnvironmen.z,
                        p_funfcnInput_tunableEnvironmen.conn,
                        q_funfcnInput_tunableEnvironmen,
                        r_funfcnInput_tunableEnvironmen,
                        s_funfcnInput_tunableEnvironmen, freq, unusedU8,
                        unusedU9, unusedUa, unusedUb);
  fx_tmp = static_cast<int>(c_funfcnInput_tunableEnvironmen) - 1;
  *fval = std::abs(freq->data[fx_tmp] - v);
  funccount = 1;
  fv = *fval;
  fw = *fval;
  xm = 0.5 * bx;
  tol1 = 1.4901161193847656E-8 * std::abs(v) + 0.00033333333333333332;
  tol2 = 2.0 * tol1;
  exitg1 = false;
  while ((!exitg1) && (std::abs(*xf - xm) > tol2 - 0.5 * (b - a))) {
    guard1 = false;
    if (std::abs(e) > tol1) {
      p = *xf - w;
      b_r = p * (*fval - fv);
      x = *xf - v;
      q = x * (*fval - fw);
      p = x * q - p * b_r;
      q = 2.0 * (q - b_r);
      if (q > 0.0) {
        p = -p;
      }

      q = std::abs(q);
      b_r = e;
      e = d;
      if ((std::abs(p) < std::abs(0.5 * q * b_r)) && (p > q * (a - *xf)) && (p <
           q * (b - *xf))) {
        d = p / q;
        x = *xf + d;
        if ((x - a < tol2) || (b - x < tol2)) {
          p = xm - *xf;
          x = p;
          if (p < 0.0) {
            x = -1.0;
          } else if (p > 0.0) {
            x = 1.0;
          } else {
            if (p == 0.0) {
              x = 0.0;
            }
          }

          d = tol1 * (x + static_cast<double>((p == 0.0)));
        }
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (*xf >= xm) {
        e = a - *xf;
      } else {
        e = b - *xf;
      }

      d = 0.3819660112501051 * e;
    }

    x = d;
    if (d < 0.0) {
      x = -1.0;
    } else if (d > 0.0) {
      x = 1.0;
    } else {
      if (d == 0.0) {
        x = 0.0;
      }
    }

    p = std::abs(d);
    if ((p > tol1) || rtIsNaN(tol1)) {
      tol1 = p;
    }

    x = *xf + (x + static_cast<double>((d == 0.0))) * tol1;
    b_linearAnalysisModal(d_funfcnInput_tunableEnvironmen,
                          e_funfcnInput_tunableEnvironmen, x,
                          f_funfcnInput_tunableEnvironmen,
                          g_funfcnInput_tunableEnvironmen,
                          h_funfcnInput_tunableEnvironmen,
                          i_funfcnInput_tunableEnvironmen,
                          j_funfcnInput_tunableEnvironmen,
                          k_funfcnInput_tunableEnvironmen,
                          l_funfcnInput_tunableEnvironmen,
                          m_funfcnInput_tunableEnvironmen,
                          n_funfcnInput_tunableEnvironmen,
                          o_funfcnInput_tunableEnvironmen,
                          p_funfcnInput_tunableEnvironmen.numEl,
                          p_funfcnInput_tunableEnvironmen.x,
                          p_funfcnInput_tunableEnvironmen.y,
                          p_funfcnInput_tunableEnvironmen.z,
                          p_funfcnInput_tunableEnvironmen.conn,
                          q_funfcnInput_tunableEnvironmen,
                          r_funfcnInput_tunableEnvironmen,
                          s_funfcnInput_tunableEnvironmen, freq, unusedU8,
                          unusedU9, unusedUa, unusedUb);
    p = std::abs(freq->data[fx_tmp] - x);
    funccount++;
    iter++;
    if (p <= *fval) {
      if (x >= *xf) {
        a = *xf;
      } else {
        b = *xf;
      }

      v = w;
      fv = fw;
      w = *xf;
      fw = *fval;
      *xf = x;
      *fval = p;
    } else {
      if (x < *xf) {
        a = x;
      } else {
        b = x;
      }

      if ((p <= fw) || (w == *xf)) {
        v = w;
        fv = fw;
        w = x;
        fw = p;
      } else {
        if ((p <= fv) || (v == *xf) || (v == w)) {
          v = x;
          fv = p;
        }
      }
    }

    xm = 0.5 * (a + b);
    tol1 = 1.4901161193847656E-8 * std::abs(*xf) + 0.00033333333333333332;
    tol2 = 2.0 * tol1;
    if ((funccount >= 500) || (iter >= 500)) {
      *exitflag = 0.0;
      exitg1 = true;
    }
  }

  emxFree_real_T(&unusedUb);
  emxFree_real_T(&unusedUa);
  emxFree_real_T(&unusedU9);
  emxFree_real_T(&unusedU8);
  emxFree_real_T(&freq);
  *output_iterations = iter;
  *output_funcCount = funccount;
}

//
// File trailer for fminbnd.cpp
//
// [EOF]
//
