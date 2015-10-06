/*
 * File: otg_xy_reallyDumb.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 17:34:35
 */

#ifndef __OTG_XY_REALLYDUMB_H__
#define __OTG_XY_REALLYDUMB_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_xy_reallyDumb_types.h"

/* Function Declarations */
extern void otg_xy_reallyDumb(double v1, double d1, double kj, double kT, double
  ks, double kd, double dT, double Tmin, double Tmax, const emxArray_real_T
  *dataVeh, double safetyS, double safetyD, double dt, double m, double kappa,
  double b_y0, double phi, double vx0, double ax0, double w, double *flag,
  emxArray_real_T *x, emxArray_real_T *y, double *T);

#endif

/*
 * File trailer for otg_xy_reallyDumb.h
 *
 * [EOF]
 */
