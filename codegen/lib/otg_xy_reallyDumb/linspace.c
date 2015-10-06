/*
 * File: linspace.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 11:15:45
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "linspace.h"
#include "otg_xy_reallyDumb_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : double d1
 *                double d2
 *                double n1
 *                emxArray_real_T *y
 * Return Type  : void
 */
void linspace(double d1, double d2, double n1, emxArray_real_T *y)
{
  double delta1;
  int i1;
  double delta2;
  int k;
  if (n1 < 0.0) {
    n1 = 0.0;
  }

  delta1 = floor(n1);
  i1 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)delta1;
  emxEnsureCapacity((emxArray__common *)y, i1, (int)sizeof(double));
  if ((int)delta1 >= 1) {
    y->data[(int)delta1 - 1] = d2;
    if (y->size[1] >= 2) {
      y->data[0] = d1;
      if (y->size[1] >= 3) {
        if (((d1 < 0.0) != (d2 < 0.0)) && ((fabs(d1) > 8.9884656743115785E+307) ||
             (fabs(d2) > 8.9884656743115785E+307))) {
          delta1 = d1 / ((double)y->size[1] - 1.0);
          delta2 = d2 / ((double)y->size[1] - 1.0);
          i1 = y->size[1];
          for (k = 0; k <= i1 - 3; k++) {
            y->data[k + 1] = (d1 + delta2 * (1.0 + (double)k)) - delta1 * (1.0 +
              (double)k);
          }
        } else {
          delta1 = (d2 - d1) / ((double)y->size[1] - 1.0);
          i1 = y->size[1];
          for (k = 0; k <= i1 - 3; k++) {
            y->data[k + 1] = d1 + (1.0 + (double)k) * delta1;
          }
        }
      }
    }
  }
}

/*
 * File trailer for linspace.c
 *
 * [EOF]
 */
