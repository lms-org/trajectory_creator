/*
 * File: otg_smart_xy.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

#ifndef __OTG_SMART_XY_H__
#define __OTG_SMART_XY_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void otg_smart_xy(double absTOL, int maxIter, double v1, double d1,
  double kj, double kT, double ks, double kd, const double dataVeh[3], double
  safetyS, double safetyD, double kappaMax, double aOrthMax, double m, double
  kappa, double b_y0, double phi, double vx0, double ax0, double w, double
  *flag1, double *flag2, double *flag3, double *flagAll, emxArray_real_T *x,
  emxArray_real_T *y, double T_data[], int T_size[2], double *TOL);

#endif

/*
 * File trailer for otg_smart_xy.h
 *
 * [EOF]
 */
