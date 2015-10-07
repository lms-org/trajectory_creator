/*
 * File: polyval.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

#ifndef __POLYVAL_H__
#define __POLYVAL_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void b_polyval(const double p_data[], const int p_size[2], const
                      emxArray_real_T *x, emxArray_real_T *y);
extern void polyval(const double p_data[], const int p_size[2], const creal_T
                    x_data[], const int x_size[2], creal_T y_data[], int y_size
                    [2]);

#endif

/*
 * File trailer for polyval.h
 *
 * [EOF]
 */
