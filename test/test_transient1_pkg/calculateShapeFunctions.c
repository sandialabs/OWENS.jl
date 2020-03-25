/*
 * File: calculateShapeFunctions.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/* Include Files */
#include "calculateShapeFunctions.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include <string.h>

/* Function Definitions */

/*
 * calculateShapeFunctions Calculates Lagrange shape functions
 *  **********************************************************************
 *  *                   Part of the SNL OWENS toolkit                    *
 *  * Developed by Sandia National Laboratories Wind Energy Technologies *
 *  *             See license.txt for disclaimer information             *
 *  **********************************************************************
 *    [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi,x)
 *
 *    This function calculates the Lagrange shape function, shape
 *    function derivative, and Jacobian to map between the local element
 *    domain and physical length of the element. The shape function
 *    derivative is defined with respect to the physical length domain. The
 *    shape functions may be linear or quadratic in order.
 *
 *       input:
 *       elementOrder = order of element: 1 linear, 2 quadratic
 *       xi           = guass point values to evaluate shape functions at
 *       x            = nodal coordinates in physical length domain
 *
 *       output:
 *       N            = shape function value at specified gauss points
 *       p_N_x        = shape function derivative w.r.t physical length
 *                      domain at specified gauss points
 *       Jac          = Jacobian for mat between local element domain and
 *                      physical length domain.
 * Arguments    : const double x[2]
 *                double N_data[]
 *                int N_size[2]
 *                double p_N_x_data[]
 *                int p_N_x_size[1]
 *                double *Jac
 * Return Type  : void
 */
void b_calculateShapeFunctions(const double x[2], double N_data[], int N_size[2],
  double p_N_x_data[], int p_N_x_size[1], double *Jac)
{
  /*  N shape function */
  /*  p_N_xi partial derivative of shape function w.r.t. xi */
  /* Linear interpolation functions */
  N_size[0] = 1;
  N_size[1] = 2;
  N_data[0] = 0.5;
  N_data[1] = 0.5;

  /* Quadratic interpolation functions */
  *Jac = -0.5 * x[0] + 0.5 * x[1];
  p_N_x_size[0] = 2;
  p_N_x_data[0] = -0.5 / *Jac;
  p_N_x_data[1] = 0.5 / *Jac;
}

/*
 * calculateShapeFunctions Calculates Lagrange shape functions
 *  **********************************************************************
 *  *                   Part of the SNL OWENS toolkit                    *
 *  * Developed by Sandia National Laboratories Wind Energy Technologies *
 *  *             See license.txt for disclaimer information             *
 *  **********************************************************************
 *    [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi,x)
 *
 *    This function calculates the Lagrange shape function, shape
 *    function derivative, and Jacobian to map between the local element
 *    domain and physical length of the element. The shape function
 *    derivative is defined with respect to the physical length domain. The
 *    shape functions may be linear or quadratic in order.
 *
 *       input:
 *       elementOrder = order of element: 1 linear, 2 quadratic
 *       xi           = guass point values to evaluate shape functions at
 *       x            = nodal coordinates in physical length domain
 *
 *       output:
 *       N            = shape function value at specified gauss points
 *       p_N_x        = shape function derivative w.r.t physical length
 *                      domain at specified gauss points
 *       Jac          = Jacobian for mat between local element domain and
 *                      physical length domain.
 * Arguments    : double xi
 *                const double x[2]
 *                double N_data[]
 *                int N_size[2]
 *                double p_N_x_data[]
 *                int p_N_x_size[1]
 *                double *Jac
 * Return Type  : void
 */
void calculateShapeFunctions(double xi, const double x[2], double N_data[], int
  N_size[2], double p_N_x_data[], int p_N_x_size[1], double *Jac)
{
  /*  N shape function */
  /*  p_N_xi partial derivative of shape function w.r.t. xi */
  /* Linear interpolation functions */
  N_size[0] = 1;
  N_size[1] = 2;
  N_data[0] = 0.5 * (1.0 - xi);
  N_data[1] = 0.5 * (xi + 1.0);

  /* Quadratic interpolation functions */
  *Jac = -0.5 * x[0] + 0.5 * x[1];
  p_N_x_size[0] = 2;
  p_N_x_data[0] = -0.5 / *Jac;
  p_N_x_data[1] = 0.5 / *Jac;
}

/*
 * File trailer for calculateShapeFunctions.c
 *
 * [EOF]
 */
