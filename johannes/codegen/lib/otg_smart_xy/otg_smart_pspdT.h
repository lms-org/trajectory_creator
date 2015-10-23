/*
 * File: otg_smart_pspdT.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

#ifndef __OTG_SMART_PSPDT_H__
#define __OTG_SMART_PSPDT_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void otg_smart_pspdT(double absTOL, int maxIter, const double S[3], const
  double D[4], double kj, double kT, double ks, double kd, const double dataVeh
  [3], double safetyS, double safetyD, double kappa, double kappaMax, double
  aOrthMax, double *flag1, double *flag2, double *flag3, double *flagAll, double
  ps_data[], int ps_size[2], double pd_data[], int pd_size[2], double T_data[],
  int T_size[2], double *TOL);

#endif

/*
 * File trailer for otg_smart_pspdT.h
 *
 * [EOF]
 */
