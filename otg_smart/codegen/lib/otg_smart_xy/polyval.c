/*
 * File: polyval.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "polyval.h"
#include "otg_smart_xy_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
void b_polyval(const double p_data[], const int p_size[2], const emxArray_real_T
               *x, emxArray_real_T *y)
{
  unsigned short uv0[2];
  int i5;
  int nc;
  int loop_ub;
  int k;
  for (i5 = 0; i5 < 2; i5++) {
    uv0[i5] = (unsigned short)x->size[i5];
  }

  i5 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i5, (int)sizeof(double));
  nc = p_size[0] * p_size[1];
  if (!(uv0[1] == 0)) {
    i5 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i5, (int)sizeof(double));
    i5 = y->size[0] * y->size[1];
    y->size[1] = uv0[1];
    emxEnsureCapacity((emxArray__common *)y, i5, (int)sizeof(double));
    loop_ub = uv0[1];
    for (i5 = 0; i5 < loop_ub; i5++) {
      y->data[i5] = p_data[0];
    }

    for (k = 0; k <= nc - 2; k++) {
      i5 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i5, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i5 = 0; i5 < loop_ub; i5++) {
        y->data[i5] = x->data[i5] * y->data[i5] + p_data[k + 1];
      }
    }
  }
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const creal_T x_data[]
 *                const int x_size[2]
 *                creal_T y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
void polyval(const double p_data[], const int p_size[2], const creal_T x_data[],
             const int x_size[2], creal_T y_data[], int y_size[2])
{
  signed char iv1[2];
  int i3;
  int loop_ub;
  int k;
  double x_data_im;
  for (i3 = 0; i3 < 2; i3++) {
    iv1[i3] = (signed char)x_size[i3];
  }

  y_size[0] = 1;
  y_size[1] = iv1[1];
  if (!(p_size[1] == 0)) {
    y_size[0] = 1;
    y_size[1] = iv1[1];
    loop_ub = iv1[1];
    for (i3 = 0; i3 < loop_ub; i3++) {
      y_data[i3].re = p_data[0];
      y_data[i3].im = 0.0;
    }

    for (k = 0; k <= p_size[1] - 2; k++) {
      y_size[0] = 1;
      y_size[1] = x_size[1];
      loop_ub = x_size[0] * x_size[1];
      for (i3 = 0; i3 < loop_ub; i3++) {
        x_data_im = x_data[i3].re * y_data[i3].im + x_data[i3].im * y_data[i3].
          re;
        y_data[i3].re = (x_data[i3].re * y_data[i3].re - x_data[i3].im *
                         y_data[i3].im) + p_data[k + 1];
        y_data[i3].im = x_data_im;
      }
    }
  }
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
