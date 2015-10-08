/*
 * File: linspace.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 14:13:47
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "linspace.h"

/* Function Definitions */

/*
 * Arguments    : double d2
 *                double y[100]
 * Return Type  : void
 */
void linspace(double d2, double y[100])
{
  double delta1;
  int k;
  y[99] = d2;
  y[0] = 0.0;
  if ((d2 < 0.0) && (fabs(d2) > 8.9884656743115785E+307)) {
    delta1 = d2 / 99.0;
    for (k = 0; k < 98; k++) {
      y[1 + k] = delta1 * (1.0 + (double)k);
    }
  } else {
    delta1 = d2 / 99.0;
    for (k = 0; k < 98; k++) {
      y[1 + k] = (1.0 + (double)k) * delta1;
    }
  }
}

/*
 * File trailer for linspace.c
 *
 * [EOF]
 */
