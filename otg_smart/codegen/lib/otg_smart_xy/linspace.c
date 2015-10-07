/*
 * File: linspace.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "linspace.h"
#include "otg_smart_xy_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : double d2
 *                short n1
 *                emxArray_real_T *y
 * Return Type  : void
 */
void linspace(double d2, short n1, emxArray_real_T *y)
{
  int i4;
  double delta1;
  int k;
  i4 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n1;
  emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
  if (n1 >= 1) {
    y->data[n1 - 1] = d2;
    if (y->size[1] >= 2) {
      y->data[0] = 0.0;
      if (y->size[1] >= 3) {
        if ((d2 < 0.0) && (fabs(d2) > 8.9884656743115785E+307)) {
          delta1 = d2 / ((double)y->size[1] - 1.0);
          i4 = y->size[1];
          for (k = 0; k <= i4 - 3; k++) {
            y->data[k + 1] = delta1 * (1.0 + (double)k);
          }
        } else {
          delta1 = d2 / ((double)y->size[1] - 1.0);
          i4 = y->size[1];
          for (k = 0; k <= i4 - 3; k++) {
            y->data[k + 1] = (1.0 + (double)k) * delta1;
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
