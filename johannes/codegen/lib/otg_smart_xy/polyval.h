/*
 * File: polyval.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

#ifndef __POLYVAL_H__
#define __POLYVAL_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void b_polyval(const double p[6], const double x[100], double y[100]);
extern void c_polyval(const double p[5], const double x[100], double y[100]);
extern void d_polyval(const double p_data[], const int p_size[2], const
                      emxArray_real_T *x, emxArray_real_T *y);
extern void polyval(const double p_data[], const int p_size[2], const double x
                    [100], double y[100]);

#endif

/*
 * File trailer for polyval.h
 *
 * [EOF]
 */
