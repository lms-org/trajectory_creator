/*
 * File: power.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "power.h"

/* Function Definitions */

/*
 * Arguments    : const double a[100]
 *                double y[100]
 * Return Type  : void
 */
void power(const double a[100], double y[100])
{
  int k;
  for (k = 0; k < 100; k++) {
    y[k] = a[k] * a[k];
  }
}

/*
 * File trailer for power.c
 *
 * [EOF]
 */
