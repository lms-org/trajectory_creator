/*
 * File: polyval.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 11:15:45
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "polyval.h"
#include "otg_xy_reallyDumb_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
void polyval(const double p_data[], const int p_size[2], const emxArray_real_T
             *x, emxArray_real_T *y)
{
  unsigned int uv0[2];
  int i2;
  int nc;
  int loop_ub;
  int k;
  for (i2 = 0; i2 < 2; i2++) {
    uv0[i2] = (unsigned int)x->size[i2];
  }

  i2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
  nc = p_size[0] * p_size[1];
  if (!((int)uv0[1] == 0)) {
    i2 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    i2 = y->size[0] * y->size[1];
    y->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    loop_ub = (int)uv0[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      y->data[i2] = p_data[0];
    }

    for (k = 0; k <= nc - 2; k++) {
      i2 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        y->data[i2] = x->data[i2] * y->data[i2] + p_data[k + 1];
      }
    }
  }
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
