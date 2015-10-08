/*
 * File: otg_smart_pspdT.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_smart_pspdT.h"
#include "otg_smart_xy_emxutil.h"
#include "otg_smart_objFun.h"
#include "mrdivide.h"
#include "otg_smart_opt_step.h"
#include "otg_smart_xy_rtwutil.h"

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
 * otg_smart_opt_pspdT gives ps, pd, T for the optimal trajectory
 * T
 *    INPUT:
 *
 *        absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed
 *            to lie within a ball with radius absTOL and middle Topt
 *        maxIter = max. Number of Iterations (50 should be suffiecient)
 *
 *
 *        S = [v0, a0, v1] s.t.
 *            s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1
 *        D = [d0, d0d, d0dd, d1] s.t.
 *            d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1
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
 *        kappa = curvature of the circle which describes the centerline
 *        kappaMax = max. curvature of the road in xy coordinate system
 *        aOrthMax = max. acceleration orthogonal to the trajectory
 *
 *    OUTPUT:
 *        ps = coefficients of the polynomial in std. matlab notation i.e.
 *            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
 *        pd = coefficients of the polynomial in std. matlab notation i.e.
 *            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
 *            + pd(6)
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
 *                const double S[3]
 *                const double D[4]
 *                double kj
 *                double kT
 *                double ks
 *                double kd
 *                const double dataVeh[3]
 *                double safetyS
 *                double safetyD
 *                double kappa
 *                double kappaMax
 *                double aOrthMax
 *                double *flag1
 *                double *flag2
 *                double *flag3
 *                double *flagAll
 *                double ps_data[]
 *                int ps_size[2]
 *                double pd_data[]
 *                int pd_size[2]
 *                double T_data[]
 *                int T_size[2]
 *                double *TOL
 * Return Type  : void
 */
void otg_smart_pspdT(double absTOL, short maxIter, const double S[3], const
                     double D[4], double kj, double kT, double ks, double kd,
                     const double dataVeh[3], double safetyS, double safetyD,
                     double kappa, double kappaMax, double aOrthMax, double
                     *flag1, double *flag2, double *flag3, double *flagAll,
                     double ps_data[], int ps_size[2], double pd_data[], int
                     pd_size[2], double T_data[], int T_size[2], double *TOL)
{
  int i1;
  double Tstart;
  double Ts[3];
  double coll[3];
  double notD[3];
  double Cs[3];
  short iter;
  int32_T exitg4;
  double b_coll[3];
  double b_notD[3];
  double b_Cs[3];
  double b_Ts[3];
  boolean_T indOk[3];
  int k;
  int ii_size[2];
  int ii_data[3];
  int idx;
  boolean_T exitg5;
  boolean_T guard1 = false;
  int i2;
  double indOk_data[3];
  emxArray_real_T *r2;
  emxArray_real_T *r3;
  emxArray_real_T *r4;
  signed char iv0[2];
  double ps[5];
  double a;
  int indOk_size[2];
  int c_size[2];
  double e;
  boolean_T firstmult;
  int32_T exitg3;
  double ed2;
  double pd[6];
  double b_a;
  int32_T exitg2;
  int32_T exitg1;

  /* % init */
  ps_size[0] = 5;
  ps_size[1] = 1;
  for (i1 = 0; i1 < 5; i1++) {
    ps_data[i1] = 0.0;
  }

  pd_size[0] = 6;
  pd_size[1] = 1;
  for (i1 = 0; i1 < 6; i1++) {
    pd_data[i1] = 0.0;
  }

  *flagAll = -1.0;
  Tstart = dataVeh[0] / (fabs((S[0] + S[2]) / 2.0 - 0.5 * dataVeh[1]) + 0.001);
  if ((Tstart <= 1.0) || (Tstart >= 20.0)) {
    Tstart = 5.0;

    /* just some hard-coded plausibility check */
  }

  /* % get the best T */
  /* This is the OTG = Optimal trajectory generation package. This is based on */
  /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
  /* Frenet Frame"(2010) by Werling et. al. */
  /*  */
  /* This is the otg_smart subpackage which relies on the assumption that there */
  /* is only one other vehicle on the road. This improves the optimation */
  /* drastically in speed and accuracy. If this fails the otg_dumb subpackage */
  /* should be used */
  /*  */
  /* -------------------------------------------------------------------------- */
  /*  */
  /* otg_smart_optT performs the optimization in T */
  /* T */
  /*    INPUT: */
  /*  */
  /*        Tstart = start value for T */
  /*        absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed */
  /*            to lie within a ball with radius absTOL and middle Topt */
  /*        maxIter = max. Number of Iterations (50 should be suffiecient) */
  /*  */
  /*        S = [v0, a0, v1] s.t. */
  /*            s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1 */
  /*        D = [d0, d0d, d0dd, d1] s.t. */
  /*            d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1 */
  /*  */
  /*        kj = weight for the jerk functional */
  /*        kT = weight for the time */
  /*        ks = weight of the longitudinal weight function */
  /*        kd = weight of the lateral weight function */
  /*  */
  /*        dataVeh = (3x1) data of other vehivles on the road */
  /*            dataVeh = [s0; v; I] s.t. */
  /*                s0 = initial distance between obstacle car i and own car */
  /*                v  = velocity of obstacle car i anlong the road */
  /*                I  = -1/1: -1: on right lane, +1: on the left lane: this */
  /*                must be the same lane as the start of the vehicle */
  /*        safetyS = min. safety distance in s direction */
  /*        safetyD = min. safety distance in D direction with sign: */
  /*            if I == -1: d(t) must be bigger  than safetyD */
  /*            if I ==  1: d(t) must be smaller than safetyD */
  /*  */
  /*        kappa = curvature of the circle which describes the centerline */
  /*        kappaMax = max. curvature of the road in xy coordinate system */
  /*        aOrthMax = max. acceleration orthogonal to the trajectory */
  /*  */
  /*    OUTPUT: */
  /*        Topt = optimal time for the trajectory */
  /*        TOL = limit on the error in T */
  /*        flag1 = flag for the subroutine otg_smart_objFun */
  /*                 1: all good */
  /*                -1: the minimum velocity is negative! */
  /*                -2: the safety point was not unique */
  /*                -3: not considered case (by the programmer) */
  /*        flag2 = flag for the subroutine otg_smart_opt_step */
  /*                 1: all good in this routine */
  /*                 0: no drivable trajectory */
  /*                -1: driveability condition isn't as expected */
  /*                -2: Ts not monotonically increasing */
  /*        flag3 = flag for this routine */
  /*                 1: all good in this routine */
  /*                 0: no drivable trajectory */
  /*                -1: number of iterations not sufficient */
  /*                    -1.1: number of iterations not sufficient + sol. unique */
  /*                    -1.2: number of iterations not sufficient + not unique */
  /*                -2: no drivable solution was found with the given number of iterations */
  /*                -10: something went wrong in one of the subroutines */
  /*  */
  /*  */
  /*  See also */
  /* % init */
  *TOL = -100.0;
  T_size[0] = 1;
  T_size[1] = 1;
  T_data[0] = 0.0;
  *flag2 = -100.0;
  for (i1 = 0; i1 < 3; i1++) {
    Ts[i1] = (0.5 + 0.5 * (double)i1) * Tstart;
  }

  otg_smart_objFun(Ts, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa,
                   kappaMax, aOrthMax, Cs, notD, coll);
  *flag1 = 1.0;
  iter = 1;
  do {
    exitg4 = 0;
    if (iter <= maxIter) {
      /*  call the subroutine for one step */
      otg_smart_opt_step(Ts, Cs, notD, coll, S, D, kj, kT, ks, kd, dataVeh,
                         safetyS, safetyD, kappa, kappaMax, aOrthMax, flag1,
                         flag2, b_Ts, b_Cs, b_notD, b_coll);
      for (i1 = 0; i1 < 3; i1++) {
        Ts[i1] = b_Ts[i1];
        Cs[i1] = b_Cs[i1];
        notD[i1] = b_notD[i1];
        coll[i1] = b_coll[i1];
      }

      if (*flag2 <= 0.0) {
        *flag3 = -10.0;
        exitg4 = 1;
      } else if ((fabs(b_Ts[2] - b_Ts[0]) < absTOL) && (b_notD[1] == 0.0) &&
                 (b_coll[1] == 0.0)) {
        *flag3 = 1.0;
        T_size[0] = 1;
        T_size[1] = 1;
        T_data[0] = b_Ts[1];
        *TOL = b_Ts[2] - b_Ts[0];
        exitg4 = 1;
      } else {
        iter++;
      }
    } else {
      /* % number of iterations was not sufficient */
      for (i1 = 0; i1 < 3; i1++) {
        indOk[i1] = ((notD[i1] == 0.0) && (coll[i1] == 0.0));
      }

      Tstart = indOk[0];
      for (k = 0; k < 2; k++) {
        Tstart += (double)indOk[k + 1];
      }

      if (Tstart > 0.0) {
        /* at least one drivable path */
        Tstart = indOk[0];
        for (k = 0; k < 2; k++) {
          Tstart += (double)indOk[k + 1];
        }

        if (Tstart == 1.0) {
          /* the drivable one is unique */
          *flag3 = -1.1;
          eml_li_find(indOk, ii_data, ii_size);
          T_size[0] = 1;
          T_size[1] = ii_size[1];
          k = ii_size[0] * ii_size[1];
          for (i1 = 0; i1 < k; i1++) {
            T_data[i1] = Ts[ii_data[i1] - 1];
          }

          *TOL = Ts[2] - Ts[0];
        } else {
          /*  the drivable path is not unique, just choose the first one */
          *flag3 = -1.1;
          idx = 0;
          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = 1 + (i1 << 1);
          }

          k = 1;
          exitg5 = false;
          while ((!exitg5) && (k < 4)) {
            guard1 = false;
            if (indOk[k - 1]) {
              idx++;
              ii_data[idx - 1] = k;
              if (idx >= 3) {
                exitg5 = true;
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              k++;
            }
          }

          if (1 > idx) {
            i2 = 0;
          } else {
            i2 = idx;
          }

          k = ii_size[0] * i2;
          for (i1 = 0; i1 < k; i1++) {
            indOk_data[i1] = ii_data[i1];
          }

          T_size[0] = 1;
          T_size[1] = 1;
          T_data[0] = Ts[(int)indOk_data[0] - 1];
          *TOL = Ts[2] - Ts[0];
        }
      } else {
        /* no drivable solution was found with the given number of iterations */
        *flag3 = -2.0;
      }

      exitg4 = 1;
    }
  } while (exitg4 == 0);

  emxInit_real_T(&r2, 2);
  emxInit_real_T(&r3, 2);
  emxInit_real_T(&r4, 2);
  if ((*flag2 <= 0.0) || (*flag3 <= 0.0)) {
  } else {
    /*  get the coefficients */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_ps returns the coefficients of the unique 4th order Polynomial for s(t) */
    /*    INPUT: */
    /*        S = [v0, a0, v1] s.t. */
    /*            s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1 */
    /*        T = end time */
    /*  */
    /*    OUTPUT: */
    /*        ps = coefficients of the polynomial in std. matlab notation i.e. */
    /*            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5) */
    /*  */
    /*  See also  */
    /* init */
    /* formulas precomputed */
    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    for (k = 0; k < iv0[1]; k++) {
      indOk_data[k] = rt_powd_snf(T_data[k], 3.0);
    }

    Tstart = S[0] * 2.0 - S[2] * 2.0;
    k = iv0[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = 1.0 / indOk_data[i1] * (Tstart + T_data[i1] * S[1]) *
        0.25;
    }

    ps[0] = indOk_data[0];
    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    for (k = 0; k < iv0[1]; k++) {
      indOk_data[k] = T_data[k] * T_data[k];
    }

    Tstart = S[0] * 3.0 - S[2] * 3.0;
    k = iv0[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = 1.0 / indOk_data[i1] * (Tstart + T_data[i1] * S[1] * 2.0)
        * -0.33333333333333331;
    }

    ps[1] = indOk_data[0];
    ps[2] = S[1] * 0.5;
    ps[3] = S[0];
    ps[4] = 0.0;
    ps_size[0] = 1;
    ps_size[1] = 5;
    for (i1 = 0; i1 < 5; i1++) {
      ps_data[i1] = ps[i1];
    }

    /* This is the OTG = Optimal Trajectory Generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_pd returns the coefficients of the unique 5th order Polynomial for d(t) */
    /*    INPUT: */
    /*        D = [d0, d0d, d0dd, d1] s.t. */
    /*            d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1 */
    /*        T = end time */
    /*  */
    /*    OUTPUT: */
    /*        pd = coefficients of the polynomial in std. matlab notation i.e. */
    /*            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t */
    /*            + pd(6) */
    /*  */
    /*  See also  */
    /* init */
    /* precomputed formulas */
    a = 6.0 * D[1];
    indOk_size[1] = T_size[1];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = T_data[i1];
    }

    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    c_size[0] = 1;
    c_size[1] = iv0[1];
    e = 5.0;
    firstmult = true;
    do {
      exitg3 = 0;
      ed2 = floor(e / 2.0);
      if (2.0 * ed2 != e) {
        if (firstmult) {
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          k = ii_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg3 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        k = indOk_size[1];
        for (i1 = 0; i1 < k; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg3 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 12.0 * D[0];
    e = 12.0 * D[3];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = -(((D[2] * (T_data[0] * T_data[i1]) + a * T_data[i1]) +
                          Tstart) - e);
    }

    i1 = r4->size[0] * r4->size[1];
    r4->size[0] = 1;
    r4->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r4, i1, (int)sizeof(double));
    k = c_size[1];
    for (i1 = 0; i1 < k; i1++) {
      r4->data[i1] = 2.0 * Cs[i1];
    }

    mrdivide(indOk_data, indOk_size, r4->data, r4->size);
    pd[0] = indOk_data[0];
    a = 3.0 * D[2];
    b_a = 16.0 * D[1];
    indOk_size[1] = T_size[1];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = T_data[i1];
    }

    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    c_size[0] = 1;
    c_size[1] = iv0[1];
    e = 4.0;
    firstmult = true;
    do {
      exitg2 = 0;
      ed2 = floor(e / 2.0);
      if (2.0 * ed2 != e) {
        if (firstmult) {
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          k = ii_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg2 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        k = indOk_size[1];
        for (i1 = 0; i1 < k; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg2 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 30.0 * D[0];
    e = 30.0 * D[3];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = ((a * (T_data[0] * T_data[i1]) + b_a * T_data[i1]) +
                        Tstart) - e;
    }

    i1 = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r3, i1, (int)sizeof(double));
    k = c_size[1];
    for (i1 = 0; i1 < k; i1++) {
      r3->data[i1] = 2.0 * Cs[i1];
    }

    mrdivide(indOk_data, indOk_size, r3->data, r3->size);
    pd[1] = indOk_data[0];
    a = 3.0 * D[2];
    b_a = 12.0 * D[1];
    indOk_size[1] = T_size[1];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = T_data[i1];
    }

    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    c_size[0] = 1;
    c_size[1] = iv0[1];
    e = 3.0;
    firstmult = true;
    do {
      exitg1 = 0;
      ed2 = floor(e / 2.0);
      if (2.0 * ed2 != e) {
        if (firstmult) {
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          k = ii_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          k = indOk_size[1];
          for (i1 = 0; i1 < k; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg1 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        k = indOk_size[1];
        for (i1 = 0; i1 < k; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg1 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 20.0 * D[0];
    e = 20.0 * D[3];
    k = T_size[1];
    for (i1 = 0; i1 < k; i1++) {
      indOk_data[i1] = -(((a * (T_data[0] * T_data[i1]) + b_a * T_data[i1]) +
                          Tstart) - e);
    }

    i1 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r2, i1, (int)sizeof(double));
    k = c_size[1];
    for (i1 = 0; i1 < k; i1++) {
      r2->data[i1] = 2.0 * Cs[i1];
    }

    mrdivide(indOk_data, indOk_size, r2->data, r2->size);
    pd[2] = indOk_data[0];
    pd[3] = D[2] / 2.0;
    pd[4] = D[1];
    pd[5] = D[0];
    pd_size[0] = 1;
    pd_size[1] = 6;
    for (i1 = 0; i1 < 6; i1++) {
      pd_data[i1] = pd[i1];
    }

    /*  see if not colliding */
    /* too lazy, should be ok */
    c_otg_smart_objFun(T_data, T_size, S, D, kj, kT, ks, kd, dataVeh, safetyS,
                       safetyD, kappa, kappaMax, aOrthMax, indOk_data,
                       indOk_size, Cs, c_size, Ts, ii_size);
    if ((Cs[0] == 0.0) && (Ts[0] == 0.0)) {
      *flagAll = 1.0;
    }
  }

  emxFree_real_T(&r4);
  emxFree_real_T(&r3);
  emxFree_real_T(&r2);
}

/*
 * File trailer for otg_smart_pspdT.c
 *
 * [EOF]
 */
