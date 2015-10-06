/*
 * File: otg_pspdT_reallyDumb.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 16:22:34
 */

#ifndef __OTG_PSPDT_REALLYDUMB_H__
#define __OTG_PSPDT_REALLYDUMB_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_xy_reallyDumb_types.h"

/* Function Declarations */
extern void otg_pspdT_reallyDumb(const double S[3], const double D[4], double kj,
  double kT, double ks, double kd, double dT, double Tmin, double Tmax, const
  emxArray_real_T *dataVeh, double safetyS, double safetyD, double dt, double
  *flag, double ps_data[], int ps_size[2], double pd_data[], int pd_size[2],
  double *T);

#endif

/*
 * File trailer for otg_pspdT_reallyDumb.h
 *
 * [EOF]
 */
