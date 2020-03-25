/*
 * File: rtGetNaN.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 25-Mar-2020 11:46:48
 */

/*
 * Abstract:
 *       MATLAB for code generation function to initialize non-finite, NaN
 */
#include "rtGetNaN.h"

/* Function: rtGetNaN ======================================================================
 * Abstract:
 * Initialize rtNaN needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real_T rtGetNaN(void)
{
  return rtNaN;
}

/* Function: rtGetNaNF =====================================================================
 * Abstract:
 * Initialize rtNaNF needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real32_T rtGetNaNF(void)
{
  return rtNaNF;
}

/*
 * File trailer for rtGetNaN.c
 *
 * [EOF]
 */
