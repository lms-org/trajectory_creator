/*
 * File: polyval.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 17:34:35
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
  unsigned int uv1[2];
  int i4;
  int nc;
  int loop_ub;
  int k;
  for (i4 = 0; i4 < 2; i4++) {
    uv1[i4] = (unsigned int)x->size[i4];
  }

  i4 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv1[1];
  emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
  nc = p_size[0] * p_size[1];
  if (!((int)uv1[1] == 0)) {
    i4 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
    i4 = y->size[0] * y->size[1];
    y->size[1] = (int)uv1[1];
    emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
    loop_ub = (int)uv1[1];
    for (i4 = 0; i4 < loop_ub; i4++) {
      y->data[i4] = p_data[0];
    }

    for (k = 0; k <= nc - 2; k++) {
      i4 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i4 = 0; i4 < loop_ub; i4++) {
        y->data[i4] = x->data[i4] * y->data[i4] + p_data[k + 1];
      }
    }
  }
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
