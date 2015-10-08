/*
 * File: otg_smart_objFun.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_smart_objFun.h"
#include "otg_pd.h"
#include "otg_ps.h"
#include "otg_smart_pspdT.h"
#include "power.h"
#include "abs.h"
#include "rdivide.h"
#include "polyval.h"
#include "polyder.h"
#include "linspace.h"
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
 * otg_smart_objFun returns for given T and parameters the value of the cost
 * function C_tot and indicators if the trajectory corresponding to T is
 * drivable
 *    INPUT:
 *
 *        T = [1xn] for n different end times T(i)
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
 *        C_tot = [1xn]: cost function values for the different T(i)
 *        notD = [1xn] = 0/1/2/3 indicator if the trajectory is not drivable: notD(i)
 *            = 0 <==> trajectory i is drivable
 *            = 1 <==> trajectory i not drivable because of curvature
 *            = 2 <==> trajectory i not drivable because of orth. acceleration
 *            = 3 <==> trajectory i not drivable because of both above
 *        coll = [1xn]: indicator if the trajectory is colliding: coll(i)
 *            = 1 <==> trajectory i is colliding
 *        flag
 *            1: all good
 *            -1: the minimum velocity is negative!
 *            -2: the safety point was not unique
 *            -3: not considered case (by the programmer)
 *
 *  See also
 * Arguments    : double T
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
 *                double *Ctot
 *                double *notD
 *                double *coll
 * Return Type  : void
 */
void b_otg_smart_objFun(double T, const double S[3], const double D[4], double
  kj, double kT, double ks, double kd, const double dataVeh[3], double safetyS,
  double safetyD, double kappa, double kappaMax, double aOrthMax, double *Ctot,
  double *notD, double *coll)
{
  double d_d_data[5];
  double PS[5];
  int ixstart;
  double dv4[6];
  double PD[6];
  double tt[100];
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int s_d_size[2];
  double s_d_data[4];
  double d_d_samples[100];
  double d_dd_samples[100];
  double s_d_samples[100];
  double I;
  double dv5[100];
  double kappas_xy[100];
  int ix;
  boolean_T exitg2;
  double mtmp;
  boolean_T exitg1;

  /* % init */
  *notD = 0.0;
  *coll = 0.0;

  /* % Get the values of ps, pd for the Ts */
  otg_ps(S, T, d_d_data);
  for (ixstart = 0; ixstart < 5; ixstart++) {
    PS[ixstart] = d_d_data[ixstart];
  }

  otg_pd(D, T, dv4);
  for (ixstart = 0; ixstart < 6; ixstart++) {
    PD[ixstart] = dv4[ixstart];
  }

  /* % Calculate the costs */
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
  *Ctot = ks * (kj * ((192.0 * rt_powd_snf(T, 3.0) * (PS[0] * PS[0]) + 144.0 *
                       (T * T) * PS[1] * PS[0]) + 36.0 * T * (PS[1] * PS[1])) +
                kT * T) + kd * (kj * ((((36.0 * T * (PD[2] * PD[2]) +
    rt_powd_snf(T, 3.0) * (192.0 * (PD[1] * PD[1]) + 240.0 * PD[2] * PD[0])) +
    720.0 * rt_powd_snf(T, 5.0) * (PD[0] * PD[0])) + 144.0 * (T * T) * PD[2] *
    PD[1]) + 720.0 * rt_powd_snf(T, 4.0) * PD[1] * PD[0]) + kT * T);

  /* % See if drivable: We use a very rough approximation here!! */
  /*  polynomial root finding */
  /* sampling based approach */
  linspace(T, tt);
  polyder(PD, d_d_data, d_d_size);

  /* first derivative w.r.t. t of d(t) */
  b_polyder(d_d_data, d_d_size, d_dd_data, d_dd_size);

  /* second derivative w.r.t. t of d(t) */
  c_polyder(PS, s_d_data, s_d_size);

  /* first derivative w.r.t. t of s(t) = v(t) */
  /* now the samples for the curvature can be calculated */
  polyval(d_d_data, d_d_size, tt, d_d_samples);
  polyval(d_dd_data, d_dd_size, tt, d_dd_samples);
  polyval(s_d_data, s_d_size, tt, s_d_samples);
  I = fabs(kappa);
  rdivide(d_dd_samples, s_d_samples, dv5);
  b_abs(dv5, kappas_xy);
  for (ixstart = 0; ixstart < 100; ixstart++) {
    kappas_xy[ixstart] += I;
  }

  ixstart = 1;
  I = kappas_xy[0];
  if (rtIsNaN(kappas_xy[0])) {
    ix = 2;
    exitg2 = false;
    while ((!exitg2) && (ix < 101)) {
      ixstart = ix;
      if (!rtIsNaN(kappas_xy[ix - 1])) {
        I = kappas_xy[ix - 1];
        exitg2 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 100) {
    while (ixstart + 1 < 101) {
      if (kappas_xy[ixstart] > I) {
        I = kappas_xy[ixstart];
      }

      ixstart++;
    }
  }

  if (I > kappaMax) {
    /* the curvature could be (according to our approximation) to big */
    *notD = 1.0;
  }

  power(s_d_samples, d_dd_samples);
  power(d_d_samples, s_d_samples);
  for (ixstart = 0; ixstart < 100; ixstart++) {
    d_dd_samples[ixstart] += s_d_samples[ixstart];
  }

  ixstart = 1;
  mtmp = d_dd_samples[0];
  if (rtIsNaN(d_dd_samples[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 101)) {
      ixstart = ix;
      if (!rtIsNaN(d_dd_samples[ix - 1])) {
        mtmp = d_dd_samples[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 100) {
    while (ixstart + 1 < 101) {
      if (d_dd_samples[ixstart] > mtmp) {
        mtmp = d_dd_samples[ixstart];
      }

      ixstart++;
    }
  }

  if (mtmp * I > aOrthMax) {
    /* the orthogonal acceleration is too much for the grip */
    *notD += 2.0;

    /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
  }

  /* % check for collision */
  I = dataVeh[2];

  /* just for readability */
  b_polyval(PD, tt, d_dd_samples);
  c_polyval(PS, tt, kappas_xy);
  for (ixstart = 0; ixstart < 100; ixstart++) {
    if ((I == -1.0) && (d_dd_samples[ixstart] < safetyD) && ((dataVeh[0] +
          dataVeh[1] * tt[ixstart]) - kappas_xy[ixstart] < safetyS)) {
      /* anfangsabstand + geschw. Auto * t - position unser Auto */
      *coll = 1.0;

      /* more than one collsion doesn't matter */
    }

    if ((I == 1.0) && (d_dd_samples[ixstart] > safetyD) && ((dataVeh[0] +
          dataVeh[1] * tt[ixstart]) - kappas_xy[ixstart] < safetyS)) {
      /* anfangsabstand + geschw. Auto * t - position unser Auto */
      *coll = 1.0;

      /* more than one collsion doesn't matter */
    }
  }
}

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
 * otg_smart_objFun returns for given T and parameters the value of the cost
 * function C_tot and indicators if the trajectory corresponding to T is
 * drivable
 *    INPUT:
 *
 *        T = [1xn] for n different end times T(i)
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
 *        C_tot = [1xn]: cost function values for the different T(i)
 *        notD = [1xn] = 0/1/2/3 indicator if the trajectory is not drivable: notD(i)
 *            = 0 <==> trajectory i is drivable
 *            = 1 <==> trajectory i not drivable because of curvature
 *            = 2 <==> trajectory i not drivable because of orth. acceleration
 *            = 3 <==> trajectory i not drivable because of both above
 *        coll = [1xn]: indicator if the trajectory is colliding: coll(i)
 *            = 1 <==> trajectory i is colliding
 *        flag
 *            1: all good
 *            -1: the minimum velocity is negative!
 *            -2: the safety point was not unique
 *            -3: not considered case (by the programmer)
 *
 *  See also
 * Arguments    : const double T_data[]
 *                const int T_size[2]
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
 *                double Ctot_data[]
 *                int Ctot_size[2]
 *                double notD_data[]
 *                int notD_size[2]
 *                double coll_data[]
 *                int coll_size[2]
 * Return Type  : void
 */
void c_otg_smart_objFun(const double T_data[], const int T_size[2], const double
  S[3], const double D[4], double kj, double kT, double ks, double kd, const
  double dataVeh[3], double safetyS, double safetyD, double kappa, double
  kappaMax, double aOrthMax, double Ctot_data[], int Ctot_size[2], double
  notD_data[], int notD_size[2], double coll_data[], int coll_size[2])
{
  int n;
  int ixstart;
  int ix;
  double PS_data[15];
  double PD_data[18];
  int i;
  double d_d_data[5];
  double dv6[6];
  double tt[100];
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int s_d_size[2];
  double s_d_data[4];
  double d_d_samples[100];
  double dd[100];
  double ss[100];
  double I;
  double dv7[100];
  boolean_T exitg2;
  double mtmp;
  boolean_T exitg1;

  /* % init */
  n = T_size[1] - 1;
  Ctot_size[0] = 1;
  Ctot_size[1] = T_size[1];
  ixstart = T_size[1];
  for (ix = 0; ix < ixstart; ix++) {
    Ctot_data[ix] = 0.0;
  }

  notD_size[0] = 1;
  notD_size[1] = T_size[1];
  ixstart = T_size[1];
  for (ix = 0; ix < ixstart; ix++) {
    notD_data[ix] = 0.0;
  }

  coll_size[0] = 1;
  coll_size[1] = T_size[1];
  ixstart = T_size[1];
  for (ix = 0; ix < ixstart; ix++) {
    coll_data[ix] = 0.0;
  }

  /* % Get the values of ps, pd for the Ts */
  ixstart = 5 * T_size[1];
  for (ix = 0; ix < ixstart; ix++) {
    PS_data[ix] = 0.0;
  }

  ixstart = 6 * T_size[1];
  for (ix = 0; ix < ixstart; ix++) {
    PD_data[ix] = 0.0;
  }

  for (i = 0; i <= n; i++) {
    otg_ps(S, T_data[i], d_d_data);
    for (ix = 0; ix < 5; ix++) {
      PS_data[ix + 5 * i] = d_d_data[ix];
    }

    otg_pd(D, T_data[i], dv6);
    for (ix = 0; ix < 6; ix++) {
      PD_data[ix + 6 * i] = dv6[ix];
    }
  }

  /* % Calculate the costs */
  for (i = 0; i <= n; i++) {
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
    Ctot_data[i] = ks * (kj * ((192.0 * rt_powd_snf(T_data[i], 3.0) * (PS_data[5
      * i] * PS_data[5 * i]) + 144.0 * (T_data[i] * T_data[i]) * PS_data[1 + 5 *
      i] * PS_data[5 * i]) + 36.0 * T_data[i] * (PS_data[1 + 5 * i] * PS_data[1
      + 5 * i])) + kT * T_data[i]) + kd * (kj * ((((36.0 * T_data[i] * (PD_data
      [2 + 6 * i] * PD_data[2 + 6 * i]) + rt_powd_snf(T_data[i], 3.0) * (192.0 *
      (PD_data[1 + 6 * i] * PD_data[1 + 6 * i]) + 240.0 * PD_data[2 + 6 * i] *
      PD_data[6 * i])) + 720.0 * rt_powd_snf(T_data[i], 5.0) * (PD_data[6 * i] *
      PD_data[6 * i])) + 144.0 * (T_data[i] * T_data[i]) * PD_data[2 + 6 * i] *
      PD_data[1 + 6 * i]) + 720.0 * rt_powd_snf(T_data[i], 4.0) * PD_data[1 + 6 *
      i] * PD_data[6 * i]) + kT * T_data[i]);
  }

  /* % See if drivable: We use a very rough approximation here!! */
  for (i = 0; i <= n; i++) {
    /*  polynomial root finding */
    /* sampling based approach */
    linspace(T_data[i], tt);
    polyder(*(double (*)[6])&PD_data[6 * i], d_d_data, d_d_size);

    /* first derivative w.r.t. t of d(t) */
    b_polyder(d_d_data, d_d_size, d_dd_data, d_dd_size);

    /* second derivative w.r.t. t of d(t) */
    c_polyder(*(double (*)[5])&PS_data[5 * i], s_d_data, s_d_size);

    /* first derivative w.r.t. t of s(t) = v(t) */
    /* now the samples for the curvature can be calculated */
    polyval(d_d_data, d_d_size, tt, d_d_samples);
    polyval(d_dd_data, d_dd_size, tt, dd);
    polyval(s_d_data, s_d_size, tt, ss);
    I = fabs(kappa);
    rdivide(dd, ss, dv7);
    b_abs(dv7, dd);
    for (ix = 0; ix < 100; ix++) {
      dd[ix] += I;
    }

    ixstart = 1;
    I = dd[0];
    if (rtIsNaN(dd[0])) {
      ix = 2;
      exitg2 = false;
      while ((!exitg2) && (ix < 101)) {
        ixstart = ix;
        if (!rtIsNaN(dd[ix - 1])) {
          I = dd[ix - 1];
          exitg2 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 100) {
      while (ixstart + 1 < 101) {
        if (dd[ixstart] > I) {
          I = dd[ixstart];
        }

        ixstart++;
      }
    }

    if (I > kappaMax) {
      /* the curvature could be (according to our approximation) to big */
      notD_data[i] = 1.0;
    }

    power(ss, dd);
    power(d_d_samples, ss);
    for (ix = 0; ix < 100; ix++) {
      dd[ix] += ss[ix];
    }

    ixstart = 1;
    mtmp = dd[0];
    if (rtIsNaN(dd[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 101)) {
        ixstart = ix;
        if (!rtIsNaN(dd[ix - 1])) {
          mtmp = dd[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 100) {
      while (ixstart + 1 < 101) {
        if (dd[ixstart] > mtmp) {
          mtmp = dd[ixstart];
        }

        ixstart++;
      }
    }

    if (mtmp * I > aOrthMax) {
      /* the orthogonal acceleration is too much for the grip */
      notD_data[i] += 2.0;

      /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
    }
  }

  /* % check for collision */
  for (i = 0; i <= n; i++) {
    I = dataVeh[2];

    /* just for readability */
    linspace(T_data[i], tt);
    b_polyval(*(double (*)[6])&PD_data[6 * i], tt, dd);
    c_polyval(*(double (*)[5])&PS_data[5 * i], tt, ss);
    for (ixstart = 0; ixstart < 100; ixstart++) {
      if ((I == -1.0) && (dd[ixstart] < safetyD) && ((dataVeh[0] + dataVeh[1] *
            tt[ixstart]) - ss[ixstart] < safetyS)) {
        /* anfangsabstand + geschw. Auto * t - position unser Auto */
        coll_data[i] = 1.0;

        /* more than one collsion doesn't matter */
      }

      if ((I == 1.0) && (dd[ixstart] > safetyD) && ((dataVeh[0] + dataVeh[1] *
            tt[ixstart]) - ss[ixstart] < safetyS)) {
        /* anfangsabstand + geschw. Auto * t - position unser Auto */
        coll_data[i] = 1.0;

        /* more than one collsion doesn't matter */
      }
    }
  }
}

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
 * otg_smart_objFun returns for given T and parameters the value of the cost
 * function C_tot and indicators if the trajectory corresponding to T is
 * drivable
 *    INPUT:
 *
 *        T = [1xn] for n different end times T(i)
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
 *        C_tot = [1xn]: cost function values for the different T(i)
 *        notD = [1xn] = 0/1/2/3 indicator if the trajectory is not drivable: notD(i)
 *            = 0 <==> trajectory i is drivable
 *            = 1 <==> trajectory i not drivable because of curvature
 *            = 2 <==> trajectory i not drivable because of orth. acceleration
 *            = 3 <==> trajectory i not drivable because of both above
 *        coll = [1xn]: indicator if the trajectory is colliding: coll(i)
 *            = 1 <==> trajectory i is colliding
 *        flag
 *            1: all good
 *            -1: the minimum velocity is negative!
 *            -2: the safety point was not unique
 *            -3: not considered case (by the programmer)
 *
 *  See also
 * Arguments    : const double T[3]
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
 *                double Ctot[3]
 *                double notD[3]
 *                double coll[3]
 * Return Type  : void
 */
void otg_smart_objFun(const double T[3], const double S[3], const double D[4],
                      double kj, double kT, double ks, double kd, const double
                      dataVeh[3], double safetyS, double safetyD, double kappa,
                      double kappaMax, double aOrthMax, double Ctot[3], double
                      notD[3], double coll[3])
{
  int ixstart;
  double PS[15];
  double PD[18];
  int i;
  double d_d_data[5];
  double dv2[6];
  double tt[100];
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int s_d_size[2];
  double s_d_data[4];
  double d_d_samples[100];
  double dd[100];
  double ss[100];
  double I;
  double dv3[100];
  int ix;
  boolean_T exitg2;
  double mtmp;
  boolean_T exitg1;

  /* % init */
  for (ixstart = 0; ixstart < 3; ixstart++) {
    notD[ixstart] = 0.0;
    coll[ixstart] = 0.0;
  }

  /* % Get the values of ps, pd for the Ts */
  for (i = 0; i < 3; i++) {
    otg_ps(S, T[i], d_d_data);
    for (ixstart = 0; ixstart < 5; ixstart++) {
      PS[ixstart + 5 * i] = d_d_data[ixstart];
    }

    otg_pd(D, T[i], dv2);
    for (ixstart = 0; ixstart < 6; ixstart++) {
      PD[ixstart + 6 * i] = dv2[ixstart];
    }
  }

  /* % Calculate the costs */
  for (i = 0; i < 3; i++) {
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
    Ctot[i] = ks * (kj * ((192.0 * rt_powd_snf(T[i], 3.0) * (PS[5 * i] * PS[5 *
      i]) + 144.0 * (T[i] * T[i]) * PS[1 + 5 * i] * PS[5 * i]) + 36.0 * T[i] *
                          (PS[1 + 5 * i] * PS[1 + 5 * i])) + kT * T[i]) + kd *
      (kj * ((((36.0 * T[i] * (PD[2 + 6 * i] * PD[2 + 6 * i]) + rt_powd_snf(T[i],
            3.0) * (192.0 * (PD[1 + 6 * i] * PD[1 + 6 * i]) + 240.0 * PD[2 + 6 *
                    i] * PD[6 * i])) + 720.0 * rt_powd_snf(T[i], 5.0) * (PD[6 *
           i] * PD[6 * i])) + 144.0 * (T[i] * T[i]) * PD[2 + 6 * i] * PD[1 + 6 *
              i]) + 720.0 * rt_powd_snf(T[i], 4.0) * PD[1 + 6 * i] * PD[6 * i])
       + kT * T[i]);
  }

  /* % See if drivable: We use a very rough approximation here!! */
  for (i = 0; i < 3; i++) {
    /*  polynomial root finding */
    /* sampling based approach */
    linspace(T[i], tt);
    polyder(*(double (*)[6])&PD[6 * i], d_d_data, d_d_size);

    /* first derivative w.r.t. t of d(t) */
    b_polyder(d_d_data, d_d_size, d_dd_data, d_dd_size);

    /* second derivative w.r.t. t of d(t) */
    c_polyder(*(double (*)[5])&PS[5 * i], s_d_data, s_d_size);

    /* first derivative w.r.t. t of s(t) = v(t) */
    /* now the samples for the curvature can be calculated */
    polyval(d_d_data, d_d_size, tt, d_d_samples);
    polyval(d_dd_data, d_dd_size, tt, dd);
    polyval(s_d_data, s_d_size, tt, ss);
    I = fabs(kappa);
    rdivide(dd, ss, dv3);
    b_abs(dv3, dd);
    for (ixstart = 0; ixstart < 100; ixstart++) {
      dd[ixstart] += I;
    }

    ixstart = 1;
    I = dd[0];
    if (rtIsNaN(dd[0])) {
      ix = 2;
      exitg2 = false;
      while ((!exitg2) && (ix < 101)) {
        ixstart = ix;
        if (!rtIsNaN(dd[ix - 1])) {
          I = dd[ix - 1];
          exitg2 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 100) {
      while (ixstart + 1 < 101) {
        if (dd[ixstart] > I) {
          I = dd[ixstart];
        }

        ixstart++;
      }
    }

    if (I > kappaMax) {
      /* the curvature could be (according to our approximation) to big */
      notD[i] = 1.0;
    }

    power(ss, dd);
    power(d_d_samples, ss);
    for (ixstart = 0; ixstart < 100; ixstart++) {
      dd[ixstart] += ss[ixstart];
    }

    ixstart = 1;
    mtmp = dd[0];
    if (rtIsNaN(dd[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 101)) {
        ixstart = ix;
        if (!rtIsNaN(dd[ix - 1])) {
          mtmp = dd[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 100) {
      while (ixstart + 1 < 101) {
        if (dd[ixstart] > mtmp) {
          mtmp = dd[ixstart];
        }

        ixstart++;
      }
    }

    if (mtmp * I > aOrthMax) {
      /* the orthogonal acceleration is too much for the grip */
      notD[i] += 2.0;

      /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
    }
  }

  /* % check for collision */
  for (i = 0; i < 3; i++) {
    I = dataVeh[2];

    /* just for readability */
    linspace(T[i], tt);
    b_polyval(*(double (*)[6])&PD[6 * i], tt, dd);
    c_polyval(*(double (*)[5])&PS[5 * i], tt, ss);
    for (ixstart = 0; ixstart < 100; ixstart++) {
      if ((I == -1.0) && (dd[ixstart] < safetyD) && ((dataVeh[0] + dataVeh[1] *
            tt[ixstart]) - ss[ixstart] < safetyS)) {
        /* anfangsabstand + geschw. Auto * t - position unser Auto */
        coll[i] = 1.0;

        /* more than one collsion doesn't matter */
      }

      if ((I == 1.0) && (dd[ixstart] > safetyD) && ((dataVeh[0] + dataVeh[1] *
            tt[ixstart]) - ss[ixstart] < safetyS)) {
        /* anfangsabstand + geschw. Auto * t - position unser Auto */
        coll[i] = 1.0;

        /* more than one collsion doesn't matter */
      }
    }
  }
}

/*
 * File trailer for otg_smart_objFun.c
 *
 * [EOF]
 */
