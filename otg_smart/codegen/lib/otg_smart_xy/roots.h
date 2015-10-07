/*
 * File: roots.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

#ifndef __ROOTS_H__
#define __ROOTS_H__

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
extern void b_roots(const double c[5], creal_T r_data[], int r_size[1]);
extern void roots(const double c_data[], const int c_size[2], creal_T r_data[],
                  int r_size[1]);

#endif

/*
 * File trailer for roots.h
 *
 * [EOF]
 */
