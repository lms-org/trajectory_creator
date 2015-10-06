/*
 * File: otg_pspdT_reallyDumb.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 17:34:35
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "otg_pspdT_reallyDumb.h"
#include "otg_xy_reallyDumb_emxutil.h"
#include "linspace.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * This is the OTG = Optimal trajectory generation package. This is based on
 * the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
 * Frenet Frame"(2010) by Werling et. al.
 *
 * --------------------------------------------------------------------------
 *
 * otg_pspdT_reallyDumb returns the coefficients of the polynomial with the
 * best cost function and no collision
 *    INPUT:
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
 *        dT = time interval between two possible end Times T
 *        Tmin = minimal T
 *        Tmax = maximal T
 *
 *        dataVeh = (3xN) data of other vehivles on the road
 *            dataVeh(:,i) = [s0i; vi; Ii] s.t.
 *                s0i = initial distance between obstacle car i and own car
 *                vi  = velocity of obstacle car i anlong the road
 *                Ii  = -1/1: -1: on right lane, +1: on the left lane
 *        safetyS = min. safety distance in s direction
 *        safetyD = min. safety distance in D direction with sign:
 *            if Ii == -1: d(t) must be bigger  than safetyD
 *            if Ii ==  1: d(t) must be smaller than safetyD
 *        dt = sampling interval for the collision detection
 *
 *    OUTPUT:
 *        flag = 1/0/-1
 *            1:  everything ok
 *            0:  solution was not unique. solution with minimal T is choosen
 *            -1: no solution: the vectors are full of zeros
 *        ps = coefficients of the polynomial in std. matlab notation i.e.
 *            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
 *        pd = coefficients of the polynomial in std. matlab notation i.e.
 *            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
 *            + pd(6)
 *        T = end time
 *
 *  See also
 * Arguments    : const double S[3]
 *                const double D[4]
 *                double kj
 *                double kT
 *                double ks
 *                double kd
 *                double dT
 *                double Tmin
 *                double Tmax
 *                const emxArray_real_T *dataVeh
 *                double safetyS
 *                double safetyD
 *                double dt
 *                double *flag
 *                double ps_data[]
 *                int ps_size[2]
 *                double pd_data[]
 *                int pd_size[2]
 *                double *T
 * Return Type  : void
 */
void otg_pspdT_reallyDumb(const double S[3], const double D[4], double kj,
  double kT, double ks, double kd, double dT, double Tmin, double Tmax, const
  emxArray_real_T *dataVeh, double safetyS, double safetyD, double dt, double
  *flag, double ps_data[], int ps_size[2], double pd_data[], int pd_size[2],
  double *T)
{
  int i1;
  emxArray_real_T *TT;
  emxArray_real_T *PS;
  double n;
  int j;
  emxArray_real_T *PD;
  emxArray_real_T *C;
  emxArray_real_T *collisions;
  int k;
  emxArray_real_T *tt;
  emxArray_real_T *ss;
  emxArray_real_T *dd;
  double ps[5];
  double pd[6];
  unsigned int uv0[2];
  int b_n;
  double y;
  boolean_T exitg4;
  int nx;
  int32_T exitg3;
  emxArray_boolean_T *ind_noCol;
  emxArray_int32_T *ii;
  boolean_T exitg2;
  boolean_T exitg1;
  boolean_T guard1 = false;

  /* % init for case of return */
  ps_size[0] = 1;
  ps_size[1] = 5;
  for (i1 = 0; i1 < 5; i1++) {
    ps_data[i1] = 0.0;
  }

  pd_size[0] = 1;
  pd_size[1] = 6;
  for (i1 = 0; i1 < 6; i1++) {
    pd_data[i1] = 0.0;
  }

  emxInit_real_T(&TT, 2);
  emxInit_real_T(&PS, 2);
  *T = 0.0;
  *flag = 1.0;

  /* % init */
  n = ceil(Tmax - Tmin) / dT;

  /* number samples in the T space */
  linspace(Tmin, Tmax, n, TT);

  /* generate all Ts */
  i1 = PS->size[0] * PS->size[1];
  PS->size[0] = 5;
  PS->size[1] = (int)n;
  emxEnsureCapacity((emxArray__common *)PS, i1, (int)sizeof(double));
  j = 5 * (int)n;
  for (i1 = 0; i1 < j; i1++) {
    PS->data[i1] = 0.0;
  }

  emxInit_real_T(&PD, 2);

  /* all the coefficients pd */
  i1 = PD->size[0] * PD->size[1];
  PD->size[0] = 6;
  PD->size[1] = (int)n;
  emxEnsureCapacity((emxArray__common *)PD, i1, (int)sizeof(double));
  j = 6 * (int)n;
  for (i1 = 0; i1 < j; i1++) {
    PD->data[i1] = 0.0;
  }

  emxInit_real_T(&C, 2);

  /* all the coefficients pd */
  i1 = C->size[0] * C->size[1];
  C->size[0] = 1;
  C->size[1] = (int)n;
  emxEnsureCapacity((emxArray__common *)C, i1, (int)sizeof(double));
  j = (int)n;
  for (i1 = 0; i1 < j; i1++) {
    C->data[i1] = 0.0;
  }

  emxInit_real_T(&collisions, 2);

  /* the cost function values */
  i1 = collisions->size[0] * collisions->size[1];
  collisions->size[0] = 1;
  collisions->size[1] = (int)n;
  emxEnsureCapacity((emxArray__common *)collisions, i1, (int)sizeof(double));
  j = (int)n;
  for (i1 = 0; i1 < j; i1++) {
    collisions->data[i1] = 0.0;
  }

  /* collisions (0: no collision, 1: collision) */
  k = 0;
  emxInit_real_T(&tt, 2);
  emxInit_real_T(&ss, 2);
  emxInit_real_T(&dd, 2);
  while (k <= (int)n - 1) {
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
    ps[0] = 1.0 / rt_powd_snf(TT->data[k], 3.0) * ((S[0] * 2.0 - S[2] * 2.0) +
      TT->data[k] * S[1]) * 0.25;
    ps[1] = 1.0 / (TT->data[k] * TT->data[k]) * ((S[0] * 3.0 - S[2] * 3.0) +
      TT->data[k] * S[1] * 2.0) * -0.33333333333333331;
    ps[2] = S[1] * 0.5;
    ps[3] = S[0];
    ps[4] = 0.0;
    for (i1 = 0; i1 < 5; i1++) {
      PS->data[i1 + PS->size[0] * k] = ps[i1];
    }

    /* get coefficients */
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
    pd[0] = -(((D[2] * (TT->data[k] * TT->data[k]) + 6.0 * D[1] * TT->data[k]) +
               12.0 * D[0]) - 12.0 * D[3]) / (2.0 * rt_powd_snf(TT->data[k], 5.0));
    pd[1] = (((3.0 * D[2] * (TT->data[k] * TT->data[k]) + 16.0 * D[1] * TT->
               data[k]) + 30.0 * D[0]) - 30.0 * D[3]) / (2.0 * rt_powd_snf
      (TT->data[k], 4.0));
    pd[2] = -(((3.0 * D[2] * (TT->data[k] * TT->data[k]) + 12.0 * D[1] *
                TT->data[k]) + 20.0 * D[0]) - 20.0 * D[3]) / (2.0 * rt_powd_snf
      (TT->data[k], 3.0));
    pd[3] = D[2] / 2.0;
    pd[4] = D[1];
    pd[5] = D[0];
    for (i1 = 0; i1 < 6; i1++) {
      PD->data[i1 + PD->size[0] * k] = pd[i1];
    }

    /* get coefficients */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_Ctot returns value of the cost total function Ctot = ks*Cs + kd*Cd */
    /*    INPUT: */
    /*        ps = coefficients of the polynomial (otg_ps) s.t. */
    /*            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5) */
    /*        pd = coefficients of the polynomial d (otg_pd) s.t. */
    /*            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t */
    /*            + pd(6) */
    /*        T = end time */
    /*        kj = weight for the jerk functional */
    /*        kT = weight for the time */
    /*        ks = weight of the longitudinal weight function */
    /*        kd = weight of the lateral weight function */
    /*  */
    /*    OUTPUT: */
    /*        Ctot = ks*Cs + kd*Cd */
    /*  */
    /*  See also  */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_Cs returns value of the cost function Cs = kj*Jts + kT*T for s */
    /*    INPUT: */
    /*        ps = coefficients of the polynomial (otg_ps) s.t. */
    /*            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5) */
    /*        T = end time */
    /*        kj = weight for the jerk functional */
    /*        kT = weight for the time */
    /*  */
    /*    OUTPUT: */
    /*        Cs = kj*Jts + kT*T */
    /*  */
    /*  See also  */
    /* precomputed formula */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_Jts returns value of the functional Jt = int_0^T  (s'''(t))^2 dt for s */
    /*    INPUT: */
    /*        ps = coefficients of the polynomial (otg_ps) s.t. */
    /*            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5) */
    /*        T = end time */
    /*  */
    /*    OUTPUT: */
    /*        Jts  = int_0^T  (s'''(t))^2 dt */
    /*  */
    /*  See also  */
    /* precomputed formula */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_Cd returns value of the cost function Cd = kj*Jtd + kT*T for d */
    /*    INPUT: */
    /*        pd = coefficients of the polynomial d (otg_pd) s.t. */
    /*            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t */
    /*            + pd(6) */
    /*        T = end time */
    /*        kj = weight for the jerk functional */
    /*        kT = weight for the time */
    /*  */
    /*    OUTPUT: */
    /*        Cd = kj*Jtd + kT*T */
    /*  */
    /*  See also  */
    /* precomputed formula */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_Jtd returns value of the functional Jt = int_0^T  (d'''(t))^2 dt for d */
    /*    INPUT: */
    /*        pd = coefficients of the polynomial d (otg_pd) s.t. */
    /*            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t */
    /*            + pd(6) */
    /*        T = end time */
    /*  */
    /*    OUTPUT: */
    /*        Jtd = int_0^T  (d'''(t))^2 dt */
    /*  */
    /*  See also  */
    /* precomputed formula */
    C->data[k] = ks * (kj * ((192.0 * rt_powd_snf(TT->data[k], 3.0) * (PS->
      data[PS->size[0] * k] * PS->data[PS->size[0] * k]) + 144.0 * (TT->data[k] *
      TT->data[k]) * PS->data[1 + PS->size[0] * k] * PS->data[PS->size[0] * k])
      + 36.0 * TT->data[k] * (PS->data[1 + PS->size[0] * k] * PS->data[1 +
      PS->size[0] * k])) + kT * TT->data[k]) + kd * (kj * ((((36.0 * TT->data[k]
      * (PD->data[2 + PD->size[0] * k] * PD->data[2 + PD->size[0] * k]) +
      rt_powd_snf(TT->data[k], 3.0) * (192.0 * (PD->data[1 + PD->size[0] * k] *
      PD->data[1 + PD->size[0] * k]) + 240.0 * PD->data[2 + PD->size[0] * k] *
      PD->data[PD->size[0] * k])) + 720.0 * rt_powd_snf(TT->data[k], 5.0) *
      (PD->data[PD->size[0] * k] * PD->data[PD->size[0] * k])) + 144.0 *
      (TT->data[k] * TT->data[k]) * PD->data[2 + PD->size[0] * k] * PD->data[1 +
      PD->size[0] * k]) + 720.0 * rt_powd_snf(TT->data[k], 4.0) * PD->data[1 +
      PD->size[0] * k] * PD->data[PD->size[0] * k]) + kT * TT->data[k]);

    /* compute total cost function */
    /* This is the OTG = Optimal trajectory generation package. This is based on */
    /* the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a */
    /* Frenet Frame"(2010) by Werling et. al. */
    /*  */
    /* -------------------------------------------------------------------------- */
    /*  */
    /* otg_collDet_dumb sees if a collision occurs (collision = 1) */
    /*    INPUT: */
    /*        ps = coefficients of the polynomial (otg_ps) s.t. */
    /*            s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5) */
    /*        pd = coefficients of the polynomial d (otg_pd) s.t. */
    /*            d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t */
    /*            + pd(6) */
    /*        T = end time */
    /*        dataVeh = (3xN) data of other vehivles on the road */
    /*            dataVeh(:,i) = [s0i; vi; Ii] s.t. */
    /*                s0i = initial distance between obstacle car i and own car */
    /*                vi  = velocity of obstacle car i anlong the road */
    /*                Ii  = -1/1: -1: on right lane, +1: on the left lane */
    /*        safetyS = min. safety distance in s direction */
    /*        safetyD = min. safety distance in D direction with sign: */
    /*            if Ii == -1: d(t) must be bigger  than safetyD */
    /*            if Ii ==  1: d(t) must be smaller than safetyD */
    /*        dt = sampling interval for the collision detection */
    /*  */
    /*    OUTPUT: */
    /*        collsion = 0/1: 0: no collsion, 1: collsion */
    /*  */
    /*  See also */
    /* number of obstacle vehicles */
    b_linspace(TT->data[k], ceil(TT->data[k] / dt), tt);

    /* samples */
    for (i1 = 0; i1 < 2; i1++) {
      uv0[i1] = (unsigned int)tt->size[i1];
    }

    i1 = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    ss->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)ss, i1, (int)sizeof(double));
    if (!((int)uv0[1] == 0)) {
      i1 = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)ss, i1, (int)sizeof(double));
      i1 = ss->size[0] * ss->size[1];
      ss->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)ss, i1, (int)sizeof(double));
      j = (int)uv0[1];
      for (i1 = 0; i1 < j; i1++) {
        ss->data[i1] = PS->data[PS->size[0] * k];
      }

      for (b_n = 0; b_n < 4; b_n++) {
        i1 = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        ss->size[1] = tt->size[1];
        emxEnsureCapacity((emxArray__common *)ss, i1, (int)sizeof(double));
        y = PS->data[(b_n + PS->size[0] * k) + 1];
        j = tt->size[0] * tt->size[1];
        for (i1 = 0; i1 < j; i1++) {
          ss->data[i1] = tt->data[i1] * ss->data[i1] + y;
        }
      }
    }

    for (i1 = 0; i1 < 2; i1++) {
      uv0[i1] = (unsigned int)tt->size[i1];
    }

    i1 = dd->size[0] * dd->size[1];
    dd->size[0] = 1;
    dd->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)dd, i1, (int)sizeof(double));
    if (!((int)uv0[1] == 0)) {
      i1 = dd->size[0] * dd->size[1];
      dd->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)dd, i1, (int)sizeof(double));
      i1 = dd->size[0] * dd->size[1];
      dd->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)dd, i1, (int)sizeof(double));
      j = (int)uv0[1];
      for (i1 = 0; i1 < j; i1++) {
        dd->data[i1] = PD->data[PD->size[0] * k];
      }

      for (b_n = 0; b_n < 5; b_n++) {
        i1 = dd->size[0] * dd->size[1];
        dd->size[0] = 1;
        dd->size[1] = tt->size[1];
        emxEnsureCapacity((emxArray__common *)dd, i1, (int)sizeof(double));
        y = PD->data[(b_n + PD->size[0] * k) + 1];
        j = tt->size[0] * tt->size[1];
        for (i1 = 0; i1 < j; i1++) {
          dd->data[i1] = tt->data[i1] * dd->data[i1] + y;
        }
      }
    }

    j = 0;
    b_n = 0;
    exitg4 = false;
    while ((!exitg4) && (b_n <= tt->size[1] - 1)) {
      nx = 0;
      do {
        exitg3 = 0;
        if (nx <= dataVeh->size[1] - 1) {
          /* go over all vehicles */
          /* just for readability */
          if ((dataVeh->data[2 + dataVeh->size[0] * nx] == -1.0) && (dd->
               data[b_n] < safetyD) && ((dataVeh->data[dataVeh->size[0] * nx] +
                dataVeh->data[1 + dataVeh->size[0] * nx] * tt->data[b_n]) -
               ss->data[b_n] < safetyS)) {
            /* anfangsabstand + geschw. Auto * t - position unser Auto */
            j = 1;

            /* more than one collsion doesn't matter */
            exitg3 = 1;
          } else if ((dataVeh->data[2 + dataVeh->size[0] * nx] == 1.0) &&
                     (dd->data[b_n] > safetyD) && ((dataVeh->data[dataVeh->size
                       [0] * nx] + dataVeh->data[1 + dataVeh->size[0] * nx] *
                       tt->data[b_n]) - ss->data[b_n] < safetyS)) {
            /* anfangsabstand + geschw. Auto * t - position unser Auto */
            j = 1;

            /* more than one collsion doesn't matter */
            exitg3 = 1;
          } else {
            nx++;
          }
        } else {
          b_n++;
          exitg3 = 2;
        }
      } while (exitg3 == 0);

      if (exitg3 == 1) {
        exitg4 = true;
      }
    }

    collisions->data[k] = j;

    /* check for collision */
    k++;
  }

  emxFree_real_T(&dd);
  emxFree_real_T(&ss);
  emxFree_real_T(&tt);
  emxInit_boolean_T(&ind_noCol, 2);
  i1 = ind_noCol->size[0] * ind_noCol->size[1];
  ind_noCol->size[0] = 1;
  ind_noCol->size[1] = collisions->size[1];
  emxEnsureCapacity((emxArray__common *)ind_noCol, i1, (int)sizeof(boolean_T));
  j = collisions->size[0] * collisions->size[1];
  for (i1 = 0; i1 < j; i1++) {
    ind_noCol->data[i1] = (collisions->data[i1] == 0.0);
  }

  emxFree_real_T(&collisions);

  /* get all with no collisions */
  if (ind_noCol->size[1] == 0) {
    y = 0.0;
  } else {
    y = ind_noCol->data[0];
    for (k = 2; k <= ind_noCol->size[1]; k++) {
      y += (double)ind_noCol->data[k - 1];
    }
  }

  if (y == 0.0) {
    /* no trajectory with no collision */
    *flag = -1.0;
  } else {
    b_n = ind_noCol->size[1];
    k = 0;
    for (nx = 1; nx <= b_n; nx++) {
      if (ind_noCol->data[nx - 1]) {
        k++;
      }
    }

    emxInit_int32_T(&ii, 2);
    i1 = ii->size[0] * ii->size[1];
    ii->size[0] = 1;
    ii->size[1] = k;
    emxEnsureCapacity((emxArray__common *)ii, i1, (int)sizeof(int));
    j = 0;
    for (nx = 1; nx <= b_n; nx++) {
      if (ind_noCol->data[nx - 1]) {
        ii->data[j] = nx;
        j++;
      }
    }

    j = 1;
    b_n = ii->size[1];
    y = C->data[ii->data[0] - 1];
    if (ii->size[1] > 1) {
      if (rtIsNaN(y)) {
        nx = 2;
        exitg2 = false;
        while ((!exitg2) && (nx <= b_n)) {
          j = nx;
          if (!rtIsNaN(C->data[ii->data[ii->size[0] * (nx - 1)] - 1])) {
            y = C->data[ii->data[ii->size[0] * (nx - 1)] - 1];
            exitg2 = true;
          } else {
            nx++;
          }
        }
      }

      if (j < ii->size[1]) {
        while (j + 1 <= b_n) {
          if (C->data[ii->data[ii->size[0] * j] - 1] < y) {
            y = C->data[ii->data[ii->size[0] * j] - 1];
          }

          j++;
        }
      }
    }

    i1 = ind_noCol->size[0] * ind_noCol->size[1];
    ind_noCol->size[0] = 1;
    ind_noCol->size[1] = C->size[1];
    emxEnsureCapacity((emxArray__common *)ind_noCol, i1, (int)sizeof(boolean_T));
    j = C->size[0] * C->size[1];
    for (i1 = 0; i1 < j; i1++) {
      ind_noCol->data[i1] = ((C->data[i1] == y) && ind_noCol->data[i1]);
    }

    nx = ind_noCol->size[1];
    b_n = 0;
    i1 = ii->size[0] * ii->size[1];
    ii->size[0] = 1;
    ii->size[1] = ind_noCol->size[1];
    emxEnsureCapacity((emxArray__common *)ii, i1, (int)sizeof(int));
    j = 1;
    exitg1 = false;
    while ((!exitg1) && (j <= nx)) {
      guard1 = false;
      if (ind_noCol->data[j - 1]) {
        b_n++;
        ii->data[b_n - 1] = j;
        if (b_n >= nx) {
          exitg1 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        j++;
      }
    }

    if (ind_noCol->size[1] == 1) {
      if (b_n == 0) {
        i1 = ii->size[0] * ii->size[1];
        ii->size[0] = 1;
        ii->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)ii, i1, (int)sizeof(int));
      }
    } else {
      i1 = ii->size[0] * ii->size[1];
      if (1 > b_n) {
        ii->size[1] = 0;
      } else {
        ii->size[1] = b_n;
      }

      emxEnsureCapacity((emxArray__common *)ii, i1, (int)sizeof(int));
    }

    i1 = C->size[0] * C->size[1];
    C->size[0] = 1;
    C->size[1] = ii->size[1];
    emxEnsureCapacity((emxArray__common *)C, i1, (int)sizeof(double));
    j = ii->size[0] * ii->size[1];
    for (i1 = 0; i1 < j; i1++) {
      C->data[i1] = ii->data[i1];
    }

    emxFree_int32_T(&ii);

    /* the second one could be useless but there is the slight chance that a colliding path has the same C value as the noncolliding minimal one */
    /* % give solution */
    j = (int)C->data[0];
    ps_size[0] = 5;
    ps_size[1] = 1;
    for (i1 = 0; i1 < 5; i1++) {
      ps_data[i1] = PS->data[i1 + PS->size[0] * (j - 1)];
    }

    j = (int)C->data[0];
    pd_size[0] = 6;
    pd_size[1] = 1;
    for (i1 = 0; i1 < 6; i1++) {
      pd_data[i1] = PD->data[i1 + PD->size[0] * (j - 1)];
    }

    *T = TT->data[(int)C->data[0] - 1];
  }

  emxFree_boolean_T(&ind_noCol);
  emxFree_real_T(&C);
  emxFree_real_T(&PD);
  emxFree_real_T(&PS);
  emxFree_real_T(&TT);
}

/*
 * File trailer for otg_pspdT_reallyDumb.c
 *
 * [EOF]
 */
