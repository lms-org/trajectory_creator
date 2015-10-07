/*
 * File: cos.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "cos.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
void b_cos(emxArray_real_T *x)
{
  int i8;
  int k;
  i8 = x->size[1];
  for (k = 0; k < i8; k++) {
    x->data[k] = cos(x->data[k]);
  }
}

/*
 * File trailer for cos.c
 *
 * [EOF]
 */
