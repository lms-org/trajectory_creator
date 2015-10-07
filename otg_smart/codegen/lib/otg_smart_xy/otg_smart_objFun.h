/*
 * File: otg_smart_objFun.h
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

#ifndef __OTG_SMART_OBJFUN_H__
#define __OTG_SMART_OBJFUN_H__

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
extern void b_otg_smart_objFun(double T, const double S[3], const double D[4],
  double kj, double kT, double ks, double kd, const double dataVeh[3], double
  safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double *Ctot, double *notD, double *coll, double *flag);
extern void c_otg_smart_objFun(const double T_data[], const int T_size[2], const
  double S[3], const double D[4], double kj, double kT, double ks, double kd,
  const double dataVeh[3], double safetyS, double safetyD, double kappa, double
  kappaMax, double aOrthMax, double Ctot_data[], int Ctot_size[2], double
  notD_data[], int notD_size[2], double coll_data[], int coll_size[2], double
  *flag);
extern void otg_smart_objFun(const double T[3], const double S[3], const double
  D[4], double kj, double kT, double ks, double kd, const double dataVeh[3],
  double safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double Ctot[3], double notD[3], double coll[3], double *flag);

#endif

/*
 * File trailer for otg_smart_objFun.h
 *
 * [EOF]
 */
