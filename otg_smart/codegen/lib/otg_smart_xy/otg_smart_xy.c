/*
 * File: otg_smart_xy.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_smart_xy_emxutil.h"
#include "repmat.h"
#include "cos.h"
#include "sin.h"
#include "polyval.h"
#include "linspace.h"
#include "otg_smart_pspdT.h"

/* Function Definitions */

/*
 * This is the OTG = Optimal trajectory generation package. This is based on
 * the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
 * Frenet Frame"(2010) by Werling et. al.
 *
 * This is the otg_smart subpackage which relies on the assumption that there
 * is only one other vehicle on the road. This improves the optimation
 * drastically in speed and accuracy. If this fails the otg_dumb subpackage
 * should be used
 *
 * --------------------------------------------------------------------------
 *
 * otg_smart_opt_xy gives the coordinates x and y of the optimal trjectory in
 * the car coordinate system
 *    INPUT:
 *
 *        absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed
 *            to lie within a ball with radius absTOL and middle Topt
 *        maxIter = max. Number of Iterations (50 should be suffiecient)
 *
 *
 *        v1 = velocity in s direction at the end of the trajectory
 *        d1 = distance from the center line at the end of the trajectory
 *
 *        kj = weight for the jerk functional
 *        kT = weight for the time
 *        ks = weight of the longitudinal weight function
 *        kd = weight of the lateral weight function
 *
 *        dataVeh = (3x1) data of other vehivles on the road
 *            dataVeh = [s0; v; I] s.t.
 *                s0 = initial distance between obstacle car i and own car
 *                v  = velocity of obstacle car i anlong the road
 *                I  = -1/1: -1: on right lane, +1: on the left lane: this
 *                must be the same lane as the start of the vehicle
 *        safetyS = min. safety distance in s direction
 *        safetyD = min. safety distance in D direction with sign:
 *            if I == -1: d(t) must be bigger  than safetyD
 *            if I ==  1: d(t) must be smaller than safetyD
 *
 *        kappaMax = max. curvature of the road in xy coordinate system
 *        aOrthMax = max. acceleration orthogonal to the trajectory
 *
 *        m = number of points wanted
 *
 *        kappa = curvature of the center line
 *        phi = angle between the x-axis of the car and the tangent on the
 *            center line at the section between center line and y axis
 *
 *    OUTPUT:
 *        x = x in the car coordinate sytsem
 *        y = y in the car coordinate sytsem
 *        T = end time
 *        TOL = limit on the error in T
 *
 *        flag1 = flag for the subroutine otg_smart_objFun
 *                 1: all good
 *                -1: the minimum velocity is negative!
 *                -2: the safety point was not unique
 *                -3: not considered case (by the programmer)
 *        flag2 = flag for the subroutine otg_smart_opt_step
 *                 1: all good in this routine
 *                 0: no drivable trajectory
 *                -1: driveability condition isn't as expected
 *                -2: Ts not monotonically increasing
 *        flag3 = flag for the otg_smart_optT subroutine
 *                 1: all good in this routine
 *                 0: no drivable trajectory
 *                -1: number of iterations not sufficient
 *                    -1.1: number of iterations not sufficient + sol. unique
 *                    -1.2: number of iterations not sufficient + not unique
 *                -2: no drivable solution was found with the given number of iterations
 *                -10: something went wrong in one of the subroutines
 *        flagAll = all together
 *                1: all is good;
 *                -1: something went somwhere wrong
 *
 *
 *  See also
 * Arguments    : double absTOL
 *                short maxIter
 *                double v1
 *                double d1
 *                double kj
 *                double kT
 *                double ks
 *                double kd
 *                const double dataVeh[3]
 *                double safetyS
 *                double safetyD
 *                double kappaMax
 *                double aOrthMax
 *                short m
 *                double kappa
 *                double b_y0
 *                double phi
 *                double vx0
 *                double ax0
 *                double w
 *                double *flag1
 *                double *flag2
 *                double *flag3
 *                double *flagAll
 *                emxArray_real_T *x
 *                emxArray_real_T *y
 *                double *T
 *                double *TOL
 * Return Type  : void
 */
void otg_smart_xy(double absTOL, short maxIter, double v1, double d1, double kj,
                  double kT, double ks, double kd, const double dataVeh[3],
                  double safetyS, double safetyD, double kappaMax, double
                  aOrthMax, short m, double kappa, double b_y0, double phi,
                  double vx0, double ax0, double w, double *flag1, double *flag2,
                  double *flag3, double *flagAll, emxArray_real_T *x,
                  emxArray_real_T *y, double *T, double *TOL)
{
  int ss;
  int loop_ub;
  double cosphi_d;
  double y0_d;
  double D[4];
  double dv0[3];
  int T_size[2];
  double T_data[3];
  int pd_size[2];
  double pd_data[36];
  int ps_size[2];
  double ps_data[25];
  emxArray_real_T *tt;
  emxArray_real_T *b_ss;
  emxArray_real_T *dd;
  emxArray_real_T *XY;
  emxArray_real_T *xy_car;
  emxArray_real_T *r0;
  emxArray_real_T *c_ss;
  int i0;
  emxArray_real_T *b_tt;
  emxArray_real_T *r1;
  double R[4];
  short k;
  double dv1[2];

  /* % init */
  ss = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = m;
  emxEnsureCapacity((emxArray__common *)x, ss, (int)sizeof(double));
  loop_ub = m;
  for (ss = 0; ss < loop_ub; ss++) {
    x->data[ss] = 0.0;
  }

  ss = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = m;
  emxEnsureCapacity((emxArray__common *)y, ss, (int)sizeof(double));
  loop_ub = m;
  for (ss = 0; ss < loop_ub; ss++) {
    y->data[ss] = 0.0;
  }

  /* % converte the state to S, D */
  /* really dirty */
  cosphi_d = -sin(phi) * w;
  y0_d = -sin(phi) * vx0;
  D[0] = -cos(phi) * b_y0;
  D[1] = -(cosphi_d * b_y0 + cos(phi) * y0_d);
  D[2] = -((-(w * w) * cos(phi) * b_y0 + 2.0 * cosphi_d * y0_d) + cos(phi) *
           -(cos(phi) * w * vx0 + sin(phi) * ax0));
  D[3] = d1;

  /* % coefficients */
  dv0[0] = cos(phi) * vx0;
  dv0[1] = -sin(phi) * w * vx0 + cos(phi) * ax0;
  dv0[2] = v1;
  otg_smart_pspdT(absTOL, maxIter, dv0, D, kj, kT, ks, kd, dataVeh, safetyS,
                  safetyD, kappa, kappaMax, aOrthMax, flag1, flag2, flag3,
                  &cosphi_d, ps_data, ps_size, pd_data, pd_size, T_data, T_size,
                  &y0_d);
  *flagAll = cosphi_d;
  *TOL = y0_d;
  *T = T_data[0];
  if (cosphi_d == -1.0) {
  } else {
    emxInit_real_T(&tt, 2);
    emxInit_real_T(&b_ss, 2);
    emxInit_real_T(&dd, 2);

    /* % calculate tt */
    b_linspace(T_data[0], m, tt);

    /* % calc s and d at the times */
    d_polyval(ps_data, ps_size, tt, b_ss);
    d_polyval(pd_data, pd_size, tt, dd);
    emxInit_real_T(&XY, 2);
    emxInit_real_T(&xy_car, 2);
    if (kappa > 0.0) {
      ss = tt->size[0] * tt->size[1];
      tt->size[0] = 1;
      tt->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)tt, ss, (int)sizeof(double));
      loop_ub = b_ss->size[0] * b_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        tt->data[ss] = kappa * b_ss->data[ss];
      }

      emxInit_real_T(&r0, 2);
      b_sin(tt);
      ss = r0->size[0] * r0->size[1];
      r0->size[0] = 1;
      r0->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)r0, ss, (int)sizeof(double));
      loop_ub = b_ss->size[0] * b_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        r0->data[ss] = kappa * b_ss->data[ss];
      }

      emxInit_real_T(&c_ss, 2);
      b_cos(r0);
      ss = c_ss->size[0] * c_ss->size[1];
      c_ss->size[0] = 1;
      c_ss->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)c_ss, ss, (int)sizeof(double));
      loop_ub = b_ss->size[0] * b_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        c_ss->data[ss] = kappa * b_ss->data[ss];
      }

      b_sin(c_ss);
      ss = c_ss->size[0] * c_ss->size[1];
      c_ss->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)c_ss, ss, (int)sizeof(double));
      ss = c_ss->size[0];
      i0 = c_ss->size[1];
      loop_ub = ss * i0;
      for (ss = 0; ss < loop_ub; ss++) {
        c_ss->data[ss] = -c_ss->data[ss];
      }

      ss = b_ss->size[0] * b_ss->size[1];
      b_ss->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)b_ss, ss, (int)sizeof(double));
      ss = b_ss->size[0];
      loop_ub = b_ss->size[1];
      loop_ub *= ss;
      for (ss = 0; ss < loop_ub; ss++) {
        b_ss->data[ss] *= kappa;
      }

      emxInit_real_T(&b_tt, 2);
      b_cos(b_ss);
      repmat(dd, xy_car);
      cosphi_d = 1.0 / kappa;
      ss = b_tt->size[0] * b_tt->size[1];
      b_tt->size[0] = 2;
      b_tt->size[1] = tt->size[1];
      emxEnsureCapacity((emxArray__common *)b_tt, ss, (int)sizeof(double));
      loop_ub = tt->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        b_tt->data[b_tt->size[0] * ss] = tt->data[tt->size[0] * ss];
      }

      loop_ub = r0->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        b_tt->data[1 + b_tt->size[0] * ss] = 1.0 - r0->data[r0->size[0] * ss];
      }

      emxFree_real_T(&r0);
      emxInit_real_T(&r1, 2);
      ss = r1->size[0] * r1->size[1];
      r1->size[0] = 2;
      r1->size[1] = c_ss->size[1];
      emxEnsureCapacity((emxArray__common *)r1, ss, (int)sizeof(double));
      loop_ub = c_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        r1->data[r1->size[0] * ss] = c_ss->data[c_ss->size[0] * ss];
      }

      emxFree_real_T(&c_ss);
      loop_ub = b_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        r1->data[1 + r1->size[0] * ss] = b_ss->data[b_ss->size[0] * ss];
      }

      ss = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = b_tt->size[1];
      emxEnsureCapacity((emxArray__common *)XY, ss, (int)sizeof(double));
      loop_ub = b_tt->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * ss] = cosphi_d * b_tt->data[i0 +
            b_tt->size[0] * ss] + r1->data[i0 + r1->size[0] * ss] * xy_car->
            data[i0 + xy_car->size[0] * ss];
        }
      }

      emxFree_real_T(&r1);
      emxFree_real_T(&b_tt);
    } else {
      emxInit_real_T(&r0, 2);
      emxInit_real_T(&c_ss, 2);
      b_repmat(m, xy_car);
      repmat(dd, r0);
      ss = c_ss->size[0] * c_ss->size[1];
      c_ss->size[0] = 2;
      c_ss->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)c_ss, ss, (int)sizeof(double));
      loop_ub = b_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        c_ss->data[c_ss->size[0] * ss] = b_ss->data[b_ss->size[0] * ss];
      }

      loop_ub = m;
      for (ss = 0; ss < loop_ub; ss++) {
        c_ss->data[1 + c_ss->size[0] * ss] = 0.0;
      }

      ss = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = c_ss->size[1];
      emxEnsureCapacity((emxArray__common *)XY, ss, (int)sizeof(double));
      loop_ub = c_ss->size[1];
      for (ss = 0; ss < loop_ub; ss++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * ss] = c_ss->data[i0 + c_ss->size[0] * ss]
            + xy_car->data[i0 + xy_car->size[0] * ss] * r0->data[i0 + r0->size[0]
            * ss];
        }
      }

      emxFree_real_T(&c_ss);
      emxFree_real_T(&r0);
    }

    emxFree_real_T(&dd);
    emxFree_real_T(&b_ss);
    emxFree_real_T(&tt);
    ss = xy_car->size[0] * xy_car->size[1];
    xy_car->size[0] = 2;
    xy_car->size[1] = m;
    emxEnsureCapacity((emxArray__common *)xy_car, ss, (int)sizeof(double));
    loop_ub = m << 1;
    for (ss = 0; ss < loop_ub; ss++) {
      xy_car->data[ss] = 0.0;
    }

    R[0] = cos(phi);
    R[2] = -sin(phi);
    R[1] = sin(phi);
    R[3] = cos(phi);
    for (k = 1; k <= m; k++) {
      dv1[0] = -sin(phi);
      dv1[1] = cos(phi);
      for (ss = 0; ss < 2; ss++) {
        cosphi_d = 0.0;
        for (i0 = 0; i0 < 2; i0++) {
          cosphi_d += R[ss + (i0 << 1)] * XY->data[i0 + XY->size[0] * (k - 1)];
        }

        xy_car->data[ss + xy_car->size[0] * (k - 1)] = cosphi_d + dv1[ss] * -D[0];
      }
    }

    emxFree_real_T(&XY);
    loop_ub = xy_car->size[1];
    ss = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)x, ss, (int)sizeof(double));
    for (ss = 0; ss < loop_ub; ss++) {
      x->data[x->size[0] * ss] = xy_car->data[xy_car->size[0] * ss];
    }

    loop_ub = xy_car->size[1];
    ss = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)y, ss, (int)sizeof(double));
    for (ss = 0; ss < loop_ub; ss++) {
      y->data[y->size[0] * ss] = xy_car->data[1 + xy_car->size[0] * ss];
    }

    emxFree_real_T(&xy_car);
  }
}

/*
 * File trailer for otg_smart_xy.c
 *
 * [EOF]
 */
