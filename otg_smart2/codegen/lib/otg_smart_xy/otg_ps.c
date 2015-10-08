/*
 * File: otg_ps.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 14:13:47
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_ps.h"
#include "otg_smart_pspdT.h"
#include "otg_smart_xy_rtwutil.h"

/* Function Definitions */

/*
 * This is the OTG = Optimal trajectory generation package. This is based on
 * the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
 * Frenet Frame"(2010) by Werling et. al.
 *
 * --------------------------------------------------------------------------
 *
 * otg_ps returns the coefficients of the unique 4th order Polynomial for s(t)
 *    INPUT:
 *        S = [v0, a0, v1] s.t.
 *            s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1
 *        T = end time
 *
 *    OUTPUT:
 *        ps = coefficients of the polynomial in std. matlab notation i.e.
 *            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
 *
 *  See also
 * Arguments    : const double S[3]
 *                double T
 *                double ps[5]
 * Return Type  : void
 */
void otg_ps(const double S[3], double T, double ps[5])
{
  /* init */
  /* formulas precomputed */
  ps[0] = 1.0 / rt_powd_snf(T, 3.0) * ((S[0] * 2.0 - S[2] * 2.0) + T * S[1]) *
    0.25;
  ps[1] = 1.0 / (T * T) * ((S[0] * 3.0 - S[2] * 3.0) + T * S[1] * 2.0) *
    -0.33333333333333331;
  ps[2] = S[1] * 0.5;
  ps[3] = S[0];
  ps[4] = 0.0;
}

/*
 * File trailer for otg_ps.c
 *
 * [EOF]
 */
