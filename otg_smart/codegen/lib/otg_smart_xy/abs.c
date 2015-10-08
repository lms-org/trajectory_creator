/*
 * File: abs.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "abs.h"

/* Function Definitions */

/*
 * Arguments    : const double x[100]
 *                double y[100]
 * Return Type  : void
 */
void b_abs(const double x[100], double y[100])
{
  int k;
  for (k = 0; k < 100; k++) {
    y[k] = fabs(x[k]);
  }
}

/*
 * File trailer for abs.c
 *
 * [EOF]
 */
