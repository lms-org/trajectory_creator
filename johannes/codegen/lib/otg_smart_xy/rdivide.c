/*
 * File: rdivide.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "rdivide.h"

/* Function Definitions */

/*
 * Arguments    : const double x[100]
 *                const double y[100]
 *                double z[100]
 * Return Type  : void
 */
void rdivide(const double x[100], const double y[100], double z[100])
{
  int i4;
  for (i4 = 0; i4 < 100; i4++) {
    z[i4] = x[i4] / y[i4];
  }
}

/*
 * File trailer for rdivide.c
 *
 * [EOF]
 */
