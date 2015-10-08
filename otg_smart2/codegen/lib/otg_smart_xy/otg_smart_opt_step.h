/*
 * File: otg_smart_opt_step.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 14:13:47
 */

#ifndef __OTG_SMART_OPT_STEP_H__
#define __OTG_SMART_OPT_STEP_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "otg_smart_xy_types.h"

/* Function Declarations */
extern void eml_li_find(const boolean_T x[3], int y_data[], int y_size[2]);
extern void otg_smart_opt_step(const double Ts[3], double Cs[3], const double
  notD[3], const double coll[3], const double S[3], const double D[4], double kj,
  double kT, double ks, double kd, const double dataVeh[3], double safetyS,
  double safetyD, double kappa, double kappaMax, double aOrthMax, double *flag1,
  double *flag2, double Ts_new[3], double Cs_new[3], double notD_new[3], double
  coll_new[3]);

#endif

/*
 * File trailer for otg_smart_opt_step.h
 *
 * [EOF]
 */
