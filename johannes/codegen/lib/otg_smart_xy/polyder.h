/*
 * File: polyder.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

#ifndef __POLYDER_H__
#define __POLYDER_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void b_polyder(const double u_data[], const int u_size[2], double a_data[],
                      int a_size[2]);
extern void c_polyder(const double u[5], double a_data[], int a_size[2]);
extern void polyder(const double u[6], double a_data[], int a_size[2]);

#endif

/*
 * File trailer for polyder.h
 *
 * [EOF]
 */
