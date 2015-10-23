/*
 * File: polyval.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "polyval.h"
#include "otg_smart_xy_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const double p[6]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
void b_polyval(const double p[6], const double x[100], double y[100])
{
  int i5;
  int k;
  for (i5 = 0; i5 < 100; i5++) {
    y[i5] = p[0];
  }

  for (k = 0; k < 5; k++) {
    for (i5 = 0; i5 < 100; i5++) {
      y[i5] = x[i5] * y[i5] + p[k + 1];
    }
  }
}

/*
 * Arguments    : const double p[5]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
void c_polyval(const double p[5], const double x[100], double y[100])
{
  int i6;
  int k;
  for (i6 = 0; i6 < 100; i6++) {
    y[i6] = p[0];
  }

  for (k = 0; k < 4; k++) {
    for (i6 = 0; i6 < 100; i6++) {
      y[i6] = x[i6] * y[i6] + p[k + 1];
    }
  }
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
void d_polyval(const double p_data[], const int p_size[2], const emxArray_real_T
               *x, emxArray_real_T *y)
{
  unsigned int uv0[2];
  int i7;
  int nc;
  int loop_ub;
  int k;
  for (i7 = 0; i7 < 2; i7++) {
    uv0[i7] = (unsigned int)x->size[i7];
  }

  i7 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
  nc = p_size[0] * p_size[1];
  if (!((int)uv0[1] == 0)) {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
    i7 = y->size[0] * y->size[1];
    y->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
    loop_ub = (int)uv0[1];
    for (i7 = 0; i7 < loop_ub; i7++) {
      y->data[i7] = p_data[0];
    }

    for (k = 0; k <= nc - 2; k++) {
      i7 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i7 = 0; i7 < loop_ub; i7++) {
        y->data[i7] = x->data[i7] * y->data[i7] + p_data[k + 1];
      }
    }
  }
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
void polyval(const double p_data[], const int p_size[2], const double x[100],
             double y[100])
{
  int i3;
  int k;
  if (!(p_size[1] == 0)) {
    for (i3 = 0; i3 < 100; i3++) {
      y[i3] = p_data[0];
    }

    for (k = 0; k <= p_size[1] - 2; k++) {
      for (i3 = 0; i3 < 100; i3++) {
        y[i3] = x[i3] * y[i3] + p_data[k + 1];
      }
    }
  }
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
