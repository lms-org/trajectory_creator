/*
 * File: otg_pd.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_pd.h"
#include "otg_smart_pspdT.h"
#include "otg_smart_xy_rtwutil.h"

/* Function Definitions */

/*
 * This is the OTG = Optimal Trajectory Generation package. This is based on
 * the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
 * Frenet Frame"(2010) by Werling et. al.
 *
 * --------------------------------------------------------------------------
 *
 * otg_pd returns the coefficients of the unique 5th order Polynomial for d(t)
 *    INPUT:
 *        D = [d0, d0d, d0dd, d1] s.t.
 *            d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1
 *        T = end time
 *
 *    OUTPUT:
 *        pd = coefficients of the polynomial in std. matlab notation i.e.
 *            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
 *            + pd(6)
 *
 *  See also
 * Arguments    : const double D[4]
 *                double T
 *                double pd[6]
 * Return Type  : void
 */
void otg_pd(const double D[4], double T, double pd[6])
{
  /* init */
  /* precomputed formulas */
  pd[0] = -(((D[2] * (T * T) + 6.0 * D[1] * T) + 12.0 * D[0]) - 12.0 * D[3]) /
    (2.0 * rt_powd_snf(T, 5.0));
  pd[1] = (((3.0 * D[2] * (T * T) + 16.0 * D[1] * T) + 30.0 * D[0]) - 30.0 * D[3])
    / (2.0 * rt_powd_snf(T, 4.0));
  pd[2] = -(((3.0 * D[2] * (T * T) + 12.0 * D[1] * T) + 20.0 * D[0]) - 20.0 * D
            [3]) / (2.0 * rt_powd_snf(T, 3.0));
  pd[3] = D[2] / 2.0;
  pd[4] = D[1];
  pd[5] = D[0];
}

/*
 * File trailer for otg_pd.c
 *
 * [EOF]
 */
