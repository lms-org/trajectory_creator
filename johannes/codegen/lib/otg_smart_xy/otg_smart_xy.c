/*
 * File: otg_smart_xy.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 14:10:50
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray__common*/

#ifndef typedef_emxArray__common
#define typedef_emxArray__common

typedef struct emxArray__common emxArray__common;

#endif                                 /*typedef_emxArray__common*/

/* Function Declarations */
static void b_abs(const double x[100], double y[100]);
static void b_otg_smart_objFun(double T, const double S[3], const double D[4],
  double kj, double kT, double ks, double kd, const double dataVeh[3], double
  safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double *Ctot, double *notD, double *coll);
static void b_polyder(const double u_data[], const int u_size[2], double a_data[],
                      int a_size[2]);
static void b_polyval(const double p[6], const double x[100], double y[100]);
static void c_otg_smart_objFun(const double T_data[], const int T_size[2], const
  double S[3], const double D[4], double kj, double kT, double ks, double kd,
  const double dataVeh[3], double safetyS, double safetyD, double kappa, double
  kappaMax, double aOrthMax, double Ctot_data[], int Ctot_size[2], double
  notD_data[], int notD_size[2], double coll_data[], int coll_size[2]);
static void c_polyder(const double u[5], double a_data[], int a_size[2]);
static void c_polyval(const double p[5], const double x[100], double y[100]);
static void d_polyval(const double p_data[], const int p_size[2], const
                      emxArray_real_T *x, emxArray_real_T *y);
static void eml_li_find(const boolean_T x[3], int y_data[], int y_size[2]);
static double eml_qrsolve(const double A_data[], const int A_size[1], double
  B_data[]);
static double eml_xnrm2(int n, const double x_data[]);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void linspace(double d2, double y[100]);
static void mrdivide(double A_data[], int A_size[2], const double B_data[],
                     const int B_size[2]);
static void otg_pd(const double D[4], double T, double pd[6]);
static void otg_ps(const double S[3], double T, double ps[5]);
static void otg_smart_objFun(const double T[3], const double S[3], const double
  D[4], double kj, double kT, double ks, double kd, const double dataVeh[3],
  double safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double Ctot[3], double notD[3], double coll[3]);
static void otg_smart_opt_step(const double Ts[3], double Cs[3], const double
  notD[3], const double coll[3], const double S[3], const double D[4], double kj,
  double kT, double ks, double kd, const double dataVeh[3], double safetyS,
  double safetyD, double kappa, double kappaMax, double aOrthMax, double *flag1,
  double *flag2, double Ts_new[3], double Cs_new[3], double notD_new[3], double
  coll_new[3]);
static void otg_smart_pspdT(double absTOL, int maxIter, const double S[3], const
  double D[4], double kj, double kT, double ks, double kd, const double dataVeh
  [3], double safetyS, double safetyD, double kappa, double kappaMax, double
  aOrthMax, double *flag1, double *flag2, double *flag3, double *flagAll, double
  ps_data[], int ps_size[2], double pd_data[], int pd_size[2], double T_data[],
  int T_size[2], double *TOL);
static void polyder(const double u[6], double a_data[], int a_size[2]);
static void polyval(const double p_data[], const int p_size[2], const double x
                    [100], double y[100]);
static void power(const double a[100], double y[100]);
static void rdivide(const double x[100], const double y[100], double z[100]);
static void repmat(const emxArray_real_T *a, emxArray_real_T *b);
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : const double x[100]
 *                double y[100]
 * Return Type  : void
 */
static void b_abs(const double x[100], double y[100])
{
  int k;
  for (k = 0; k < 100; k++) {
    y[k] = fabs(x[k]);
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
static void b_otg_smart_objFun(double T, const double S[3], const double D[4],
  double kj, double kT, double ks, double kd, const double dataVeh[3], double
  safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double *Ctot, double *notD, double *coll)
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
 * Arguments    : const double u_data[]
 *                const int u_size[2]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
static void b_polyder(const double u_data[], const int u_size[2], double a_data[],
                      int a_size[2])
{
  int nymax;
  int nlead0;
  int k;
  if (u_size[1] < 2) {
    nymax = 0;
  } else {
    nymax = u_size[1] - 2;
  }

  a_size[0] = 1;
  a_size[1] = nymax + 1;
  if ((u_size[1] == 0) || (u_size[1] == 1)) {
    a_data[0] = 0.0;
  } else {
    nlead0 = -1;
    k = 1;
    while ((k <= nymax) && (u_data[k - 1] == 0.0)) {
      nlead0++;
      k++;
    }

    nymax -= nlead0;
    a_size[0] = 1;
    a_size[1] = nymax;
    for (k = 1; k <= nymax; k++) {
      a_data[k - 1] = u_data[k + nlead0];
    }
  }

  nymax = a_size[1] - 1;
  for (k = 0; k + 1 <= nymax; k++) {
    a_data[k] *= (double)(nymax - k) + 1.0;
  }

  if ((!(u_size[1] == 0)) && (!((!rtIsInf(u_data[u_size[1] - 1])) && (!rtIsNaN
         (u_data[u_size[1] - 1]))))) {
    a_data[a_size[1] - 1] = rtNaN;
  }
}

/*
 * Arguments    : const double p[6]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
static void b_polyval(const double p[6], const double x[100], double y[100])
{
  int i5;
  int k;
  for (i5 = 0; i5 < 100; i5++) {
    y[i5] = p[0];
  }

  for (k = 0; k < 5; k++) {
    for (i5 = 0; i5 < 100; i5++) {
      y[i5] = x[i5] * y[i5] + p[k + 1];
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
static void c_otg_smart_objFun(const double T_data[], const int T_size[2], const
  double S[3], const double D[4], double kj, double kT, double ks, double kd,
  const double dataVeh[3], double safetyS, double safetyD, double kappa, double
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
 * Arguments    : const double u[5]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
static void c_polyder(const double u[5], double a_data[], int a_size[2])
{
  int nlead0;
  int k;
  nlead0 = 0;
  k = 1;
  while ((k < 4) && (u[k - 1] == 0.0)) {
    nlead0++;
    k++;
  }

  a_size[0] = 1;
  a_size[1] = 4 - nlead0;
  for (k = 0; k + 1 <= 4 - nlead0; k++) {
    a_data[k] = u[k + nlead0];
  }

  for (k = 0; k + 1 <= 3 - nlead0; k++) {
    a_data[k] *= (double)(3 - (nlead0 + k)) + 1.0;
  }

  if (!((!rtIsInf(u[4])) && (!rtIsNaN(u[4])))) {
    a_data[3 - nlead0] = rtNaN;
  }
}

/*
 * Arguments    : const double p[5]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
static void c_polyval(const double p[5], const double x[100], double y[100])
{
  int i6;
  int k;
  for (i6 = 0; i6 < 100; i6++) {
    y[i6] = p[0];
  }

  for (k = 0; k < 4; k++) {
    for (i6 = 0; i6 < 100; i6++) {
      y[i6] = x[i6] * y[i6] + p[k + 1];
    }
  }
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void d_polyval(const double p_data[], const int p_size[2], const
                      emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int uv0[2];
  int i7;
  int nc;
  int loop_ub;
  int k;
  for (i7 = 0; i7 < 2; i7++) {
    uv0[i7] = (unsigned int)x->size[i7];
  }

  i7 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
  nc = p_size[0] * p_size[1];
  if (!((int)uv0[1] == 0)) {
    i7 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
    i7 = y->size[0] * y->size[1];
    y->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
    loop_ub = (int)uv0[1];
    for (i7 = 0; i7 < loop_ub; i7++) {
      y->data[i7] = p_data[0];
    }

    for (k = 0; k <= nc - 2; k++) {
      i7 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i7 = 0; i7 < loop_ub; i7++) {
        y->data[i7] = x->data[i7] * y->data[i7] + p_data[k + 1];
      }
    }
  }
}

/*
 * Arguments    : const boolean_T x[3]
 *                int y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
static void eml_li_find(const boolean_T x[3], int y_data[], int y_size[2])
{
  int k;
  int i;
  k = 0;
  for (i = 0; i < 3; i++) {
    if (x[i]) {
      k++;
    }
  }

  y_size[0] = 1;
  y_size[1] = k;
  k = 0;
  for (i = 0; i < 3; i++) {
    if (x[i]) {
      y_data[k] = i + 1;
      k++;
    }
  }
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[1]
 *                double B_data[]
 * Return Type  : double
 */
static double eml_qrsolve(const double A_data[], const int A_size[1], double
  B_data[])
{
  double Y;
  int m;
  int mn;
  int loop_ub;
  int knt;
  double b_A_data[3];
  int b_m;
  int b_mn;
  double atmp;
  double tau_data_idx_0;
  double xnorm;
  m = A_size[0] - 2;
  mn = (int)fmin(A_size[0], 1.0);
  loop_ub = A_size[0];
  for (knt = 0; knt < loop_ub; knt++) {
    b_A_data[knt] = A_data[knt];
  }

  b_m = A_size[0];
  if (A_size[0] <= 1) {
    b_mn = A_size[0];
  } else {
    b_mn = 1;
  }

  if (A_size[0] == 0) {
  } else {
    loop_ub = 1;
    while (loop_ub <= b_mn) {
      if (1 < b_m) {
        atmp = b_A_data[0];
        tau_data_idx_0 = 0.0;
        xnorm = eml_xnrm2(b_m - 1, b_A_data);
        if (xnorm != 0.0) {
          xnorm = hypot(b_A_data[0], xnorm);
          if (b_A_data[0] >= 0.0) {
            xnorm = -xnorm;
          }

          if (fabs(xnorm) < 1.0020841800044864E-292) {
            knt = 0;
            do {
              knt++;
              for (loop_ub = 1; loop_ub + 1 <= b_m; loop_ub++) {
                b_A_data[loop_ub] *= 9.9792015476736E+291;
              }

              xnorm *= 9.9792015476736E+291;
              atmp *= 9.9792015476736E+291;
            } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

            xnorm = eml_xnrm2(b_m - 1, b_A_data);
            xnorm = hypot(atmp, xnorm);
            if (atmp >= 0.0) {
              xnorm = -xnorm;
            }

            tau_data_idx_0 = (xnorm - atmp) / xnorm;
            atmp = 1.0 / (atmp - xnorm);
            for (loop_ub = 1; loop_ub + 1 <= b_m; loop_ub++) {
              b_A_data[loop_ub] *= atmp;
            }

            for (loop_ub = 1; loop_ub <= knt; loop_ub++) {
              xnorm *= 1.0020841800044864E-292;
            }

            atmp = xnorm;
          } else {
            tau_data_idx_0 = (xnorm - b_A_data[0]) / xnorm;
            atmp = 1.0 / (b_A_data[0] - xnorm);
            for (loop_ub = 1; loop_ub + 1 <= b_m; loop_ub++) {
              b_A_data[loop_ub] *= atmp;
            }

            atmp = xnorm;
          }
        }
      } else {
        atmp = b_A_data[0];
        tau_data_idx_0 = 0.0;
      }

      b_A_data[0] = atmp;
      loop_ub = 2;
    }
  }

  atmp = 0.0;
  if (mn > 0) {
    xnorm = fmax(A_size[0], 1.0) * fabs(b_A_data[0]) * 2.2204460492503131E-16;
    loop_ub = 0;
    while ((loop_ub <= 0) && (!(fabs(b_A_data[0]) <= xnorm))) {
      atmp++;
      loop_ub = 1;
    }
  }

  Y = 0.0;
  loop_ub = 0;
  while (loop_ub <= mn - 1) {
    if (tau_data_idx_0 != 0.0) {
      xnorm = B_data[0];
      for (loop_ub = 1; loop_ub - 1 <= m; loop_ub++) {
        xnorm += b_A_data[loop_ub] * B_data[loop_ub];
      }

      xnorm *= tau_data_idx_0;
      if (xnorm != 0.0) {
        B_data[0] -= xnorm;
        for (loop_ub = 1; loop_ub - 1 <= m; loop_ub++) {
          B_data[loop_ub] -= b_A_data[loop_ub] * xnorm;
        }
      }
    }

    loop_ub = 1;
  }

  for (loop_ub = 0; loop_ub < (int)atmp; loop_ub++) {
    Y = B_data[0];
  }

  for (loop_ub = 0; loop_ub < (int)-(1.0 + (-1.0 - atmp)); loop_ub++) {
    Y /= b_A_data[0];
  }

  return Y;
}

/*
 * Arguments    : int n
 *                const double x_data[]
 * Return Type  : double
 */
static double eml_xnrm2(int n, const double x_data[])
{
  double y;
  double scale;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x_data[1]);
  } else {
    scale = 2.2250738585072014E-308;
    for (k = 2; k <= n + 1; k++) {
      absxk = fabs(x_data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

/*
 * Arguments    : emxArray__common *emxArray
 *                int oldNumel
 *                int elementSize
 * Return Type  : void
 */
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize)
{
  int newNumel;
  int i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      i <<= 1;
    }

    newData = calloc((unsigned int)i, (unsigned int)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : double d2
 *                double y[100]
 * Return Type  : void
 */
static void linspace(double d2, double y[100])
{
  double delta1;
  int k;
  y[99] = d2;
  y[0] = 0.0;
  if ((d2 < 0.0) && (fabs(d2) > 8.9884656743115785E+307)) {
    delta1 = d2 / 99.0;
    for (k = 0; k < 98; k++) {
      y[1 + k] = delta1 * (1.0 + (double)k);
    }
  } else {
    delta1 = d2 / 99.0;
    for (k = 0; k < 98; k++) {
      y[1 + k] = (1.0 + (double)k) * delta1;
    }
  }
}

/*
 * Arguments    : double A_data[]
 *                int A_size[2]
 *                const double B_data[]
 *                const int B_size[2]
 * Return Type  : void
 */
static void mrdivide(double A_data[], int A_size[2], const double B_data[],
                     const int B_size[2])
{
  double b_B_data[3];
  int b_B_size[1];
  int loop_ub;
  int i8;
  double b_A_data[3];
  double d2;
  if ((A_size[1] == 0) || (B_size[1] == 0)) {
    A_size[0] = 1;
    A_size[1] = 1;
    A_data[0] = 0.0;
  } else if (1 == B_size[1]) {
    if (A_size[1] == 0) {
    } else {
      A_data[0] *= 1.0 / B_data[0];
    }
  } else {
    b_B_size[0] = B_size[1];
    loop_ub = B_size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_B_data[i8] = B_data[B_size[0] * i8];
    }

    loop_ub = A_size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_A_data[i8] = A_data[A_size[0] * i8];
    }

    d2 = eml_qrsolve(b_B_data, b_B_size, b_A_data);
    A_size[0] = 1;
    A_size[1] = 1;
    A_data[0] = d2;
  }
}

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
static void otg_pd(const double D[4], double T, double pd[6])
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
static void otg_ps(const double S[3], double T, double ps[5])
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
static void otg_smart_objFun(const double T[3], const double S[3], const double
  D[4], double kj, double kT, double ks, double kd, const double dataVeh[3],
  double safetyS, double safetyD, double kappa, double kappaMax, double aOrthMax,
  double Ctot[3], double notD[3], double coll[3])
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
 * otg_smart_opt_step performs one step of the iterative search for the best
 * T
 *    INPUT:
 *
 *        Ts = [1x3] 3 time points of the current optim step
 *        Cs = [1x3] 3 values of the cost function belonging to the Ts
 *        notD = [1x3] 3 values of drivability cond. belonging to the Ts
 *            = 0 <==> trajectory i is drivable
 *            = 1 <==> trajectory i not drivable because of curvature
 *            = 2 <==> trajectory i not drivable because of orth. acceleration
 *            = 3 <==> trajectory i not drivable because of both above
 *        coll = [1x3] 3 values of collision cond. belonging to the Ts
 *            = 1 <==> trajectory i is colliding
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
 *        Ts_new = [1x3] 3 time points of the next optim step
 *        Cs_new = [1x3] 3 values of the cost function belonging to the Ts_new
 *        notD_new = [1x3] 3 values of drivability cond. belonging to the
 *                Ts_new
 *            = 0 <==> trajectory i is drivable
 *            = 1 <==> trajectory i not drivable because of curvature
 *            = 2 <==> trajectory i not drivable because of orth. acceleration
 *            = 3 <==> trajectory i not drivable because of both above
 *        coll_new = [1x3] 3 values of collision cond. belonging to the
 *                Ts_new
 *            = 1 <==> trajectory i is colliding
 *        flag1 = flag for the subroutine otg_smart_objFun
 *                 1: all good
 *                -1: the minimum velocity is negative!
 *                -2: the safety point was not unique
 *                -3: not considered case (by the programmer)
 *        flag2 = flag for this routine
 *                 1: all good in this routine
 *                 0: no drivable trajectory
 *                -1: driveability condition isn't as expected
 *                -2: Ts not monotonically increasing
 *
 *  See also
 * Arguments    : const double Ts[3]
 *                double Cs[3]
 *                const double notD[3]
 *                const double coll[3]
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
 *                double Ts_new[3]
 *                double Cs_new[3]
 *                double notD_new[3]
 *                double coll_new[3]
 * Return Type  : void
 */
static void otg_smart_opt_step(const double Ts[3], double Cs[3], const double
  notD[3], const double coll[3], const double S[3], const double D[4], double kj,
  double kT, double ks, double kd, const double dataVeh[3], double safetyS,
  double safetyD, double kappa, double kappaMax, double aOrthMax, double *flag1,
  double *flag2, double Ts_new[3], double Cs_new[3], double notD_new[3], double
  coll_new[3])
{
  int ixstart;
  boolean_T b0;
  boolean_T x[3];
  double y;
  double Tb;
  double Ta;
  int ix;
  boolean_T exitg1;
  int tmp_size[2];
  int tmp_data[3];
  double colla;
  double notDa;
  double Ca;
  double collb;
  double notDb;
  double Cb;

  /* % init */
  *flag2 = 1.0;
  for (ixstart = 0; ixstart < 3; ixstart++) {
    Ts_new[ixstart] = 0.0;
    Cs_new[ixstart] = 0.0;
    notD_new[ixstart] = 0.0;
    coll_new[ixstart] = 0.0;
  }

  *flag1 = 0.0;

  /* % check for order */
  if ((Ts[0] < Ts[1]) && (Ts[1] < Ts[2])) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (!b0) {
    *flag2 = -2.0;
  } else {
    /* % check if there is an error with the drivability */
    if ((notD[2] > 0.0) && ((notD[0] == 0.0) || (notD[1] == 0.0))) {
      /* all points should be non drivable */
      /* something went wrong about the drivability */
      *flag2 = -1.0;
    } else {
      /* % check if there is an error with the collisions */
      if ((coll[0] > 0.0) && ((coll[1] == 0.0) || (coll[2] == 0.0))) {
        /* all points should be colliding */
        /* something went wrong about the collsion */
        *flag2 = -3.0;
      } else {
        /* % check if no valid trajectory exists */
        for (ixstart = 0; ixstart < 3; ixstart++) {
          x[ixstart] = ((notD[ixstart] > 0.0) && (coll[ixstart] > 0.0));
        }

        y = x[0];
        for (ixstart = 0; ixstart < 2; ixstart++) {
          y += (double)x[ixstart + 1];
        }

        if (y > 0.0) {
          /* there is at least one traj wich is both not drivable and not colliding */
          /* --> no valid trajectory exists */
          *flag2 = 0.0;
        } else {
          /* % all non drivable */
          if ((notD[0] > 0.0) && (notD[1] > 0.0) && (notD[2] > 0.0)) {
            Ts_new[0] = Ts[1];
            Ts_new[1] = Ts[2];
            Ts_new[2] = Ts[2] + (Ts[2] - Ts[0]);
            b_otg_smart_objFun(Ts_new[2], S, D, kj, kT, ks, kd, dataVeh, safetyS,
                               safetyD, kappa, kappaMax, aOrthMax, &y, &Ta, &Tb);
            *flag1 = 1.0;
            Cs_new[0] = Cs[1];
            Cs_new[1] = Cs[2];
            Cs_new[2] = y;
            notD_new[0] = notD[1];
            notD_new[1] = notD[2];
            notD_new[2] = Ta;
            coll_new[0] = coll[1];
            coll_new[1] = coll[2];
            coll_new[2] = Tb;
          } else {
            /* % all collidiong */
            if ((coll[0] > 0.0) && (coll[1] > 0.0) && (coll[2] > 0.0)) {
              Ts_new[0] = Ts[0] - (Ts[2] - Ts[0]);

              /* fibonacci growth */
              Ts_new[1] = Ts[0];
              Ts_new[2] = Ts[1];
              b_otg_smart_objFun(Ts_new[0], S, D, kj, kT, ks, kd, dataVeh,
                                 safetyS, safetyD, kappa, kappaMax, aOrthMax, &y,
                                 &Ta, &Tb);
              *flag1 = 1.0;
              Cs_new[0] = y;
              Cs_new[1] = Cs[0];
              Cs_new[2] = Cs[1];
              notD_new[0] = Ta;
              notD_new[1] = notD[0];
              notD_new[2] = notD[1];
              coll_new[0] = Tb;
              coll_new[1] = 1.0;
              coll_new[2] = 1.0;
            } else {
              /* % Test all the cases */
              /* workaround for cpp not offering a value inf: if one trajectory is not */
              /* drivable or colliding set its c value to double the max */
              ixstart = 1;
              y = Cs[0];
              if (rtIsNaN(Cs[0])) {
                ix = 2;
                exitg1 = false;
                while ((!exitg1) && (ix < 4)) {
                  ixstart = ix;
                  if (!rtIsNaN(Cs[ix - 1])) {
                    y = Cs[ix - 1];
                    exitg1 = true;
                  } else {
                    ix++;
                  }
                }
              }

              if (ixstart < 3) {
                while (ixstart + 1 < 4) {
                  if (Cs[ixstart] > y) {
                    y = Cs[ixstart];
                  }

                  ixstart++;
                }
              }

              for (ixstart = 0; ixstart < 3; ixstart++) {
                x[ixstart] = ((notD[ixstart] > 0.0) || (coll[ixstart] > 0.0));
              }

              eml_li_find(x, tmp_data, tmp_size);
              ix = tmp_size[0] * tmp_size[1];
              for (ixstart = 0; ixstart < ix; ixstart++) {
                Cs[tmp_data[ixstart] - 1] = y;
              }

              if ((Cs[0] < Cs[1]) && (Cs[1] <= Cs[2])) {
                Ts_new[0] = Ts[0] - (Ts[2] - Ts[0]);

                /* fibonacci growth */
                Ts_new[1] = Ts[0];
                Ts_new[2] = Ts[1];
                b_otg_smart_objFun(Ts_new[0], S, D, kj, kT, ks, kd, dataVeh,
                                   safetyS, safetyD, kappa, kappaMax, aOrthMax,
                                   &y, &Ta, &Tb);
                *flag1 = 1.0;
                Cs_new[0] = y;
                Cs_new[1] = Cs[0];
                Cs_new[2] = Cs[1];
                notD_new[0] = Ta;
                notD_new[1] = notD[0];
                notD_new[2] = notD[1];
                coll_new[0] = Tb;
                coll_new[1] = coll[0];
                coll_new[2] = coll[1];
              } else if ((Cs[0] >= Cs[1]) && (Cs[1] > Cs[2])) {
                Ts_new[0] = Ts[1];
                Ts_new[1] = Ts[2];
                Ts_new[2] = Ts[2] + (Ts[2] - Ts[0]);

                /* fibonacci growth */
                b_otg_smart_objFun(Ts_new[2], S, D, kj, kT, ks, kd, dataVeh,
                                   safetyS, safetyD, kappa, kappaMax, aOrthMax,
                                   &y, &Ta, &Tb);
                *flag1 = 1.0;
                Cs_new[0] = Cs[1];
                Cs_new[1] = Cs[2];
                Cs_new[2] = y;
                notD_new[0] = notD[1];
                notD_new[1] = notD[2];
                notD_new[2] = Ta;
                coll_new[0] = coll[1];
                coll_new[1] = coll[2];
                coll_new[2] = Tb;
              } else {
                if ((Cs[0] >= Cs[1]) && (Cs[1] <= Cs[2])) {
                  Ta = (Ts[0] + Ts[1]) / 2.0;
                  Tb = (Ts[1] + Ts[2]) / 2.0;
                  b_otg_smart_objFun(Ta, S, D, kj, kT, ks, kd, dataVeh, safetyS,
                                     safetyD, kappa, kappaMax, aOrthMax, &Ca,
                                     &notDa, &colla);
                  b_otg_smart_objFun(Tb, S, D, kj, kT, ks, kd, dataVeh, safetyS,
                                     safetyD, kappa, kappaMax, aOrthMax, &Cb,
                                     &notDb, &collb);
                  *flag1 = 1.0;
                  if ((notDa > 0.0) || (colla > 0.0)) {
                    Ca = y;
                  }

                  if ((notDb > 0.0) || (collb > 0.0)) {
                    Cb = y;
                  }

                  if (Cb < Cs[1]) {
                    /* out of the five points the minumum must be somewhere between the */
                    /* three outter right ones */
                    Ts_new[0] = Ts[1];
                    Ts_new[1] = Tb;
                    Ts_new[2] = Ts[2];
                    Cs_new[0] = Cs[1];
                    Cs_new[1] = Cb;
                    Cs_new[2] = Cs[2];
                    notD_new[0] = notD[1];
                    notD_new[1] = notDb;
                    notD_new[2] = notD[2];
                    coll_new[0] = coll[1];
                    coll_new[1] = collb;
                    coll_new[2] = coll[2];
                  } else if (Ca < Cs[1]) {
                    /* out of the five points the minumum must be somewhere between the */
                    /* three outter left ones */
                    Ts_new[0] = Ts[0];
                    Ts_new[1] = Ta;
                    Ts_new[2] = Ts[1];
                    Cs_new[0] = Cs[0];
                    Cs_new[1] = Ca;
                    Cs_new[2] = Cs[1];
                    notD_new[0] = notD[0];
                    notD_new[1] = notDa;
                    notD_new[2] = notD[1];
                    coll_new[0] = coll[0];
                    coll_new[1] = colla;
                    coll_new[2] = coll[1];
                  } else {
                    if ((Ca >= Cs[1]) && (Cb >= Cs[1])) {
                      /* out of the five points the minumum must be somewhere between the */
                      /* three middle ones */
                      Ts_new[0] = Ta;
                      Ts_new[1] = Ts[1];
                      Ts_new[2] = Tb;
                      Cs_new[0] = Ca;
                      Cs_new[1] = Cs[1];
                      Cs_new[2] = Cb;
                      notD_new[0] = notDa;
                      notD_new[1] = notD[1];
                      notD_new[2] = notDb;
                      coll_new[0] = colla;
                      coll_new[1] = coll[1];
                      coll_new[2] = collb;
                    }
                  }
                }
              }
            }
          }
        }
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
 *                int maxIter
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
static void otg_smart_pspdT(double absTOL, int maxIter, const double S[3], const
  double D[4], double kj, double kT, double ks, double kd, const double dataVeh
  [3], double safetyS, double safetyD, double kappa, double kappaMax, double
  aOrthMax, double *flag1, double *flag2, double *flag3, double *flagAll, double
  ps_data[], int ps_size[2], double pd_data[], int pd_size[2], double T_data[],
  int T_size[2], double *TOL)
{
  int i1;
  double Tstart;
  double Ts[3];
  double coll[3];
  double notD[3];
  double Cs[3];
  int iter;
  int32_T exitg4;
  double b_coll[3];
  double b_notD[3];
  double b_Cs[3];
  double b_Ts[3];
  boolean_T indOk[3];
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
      for (iter = 0; iter < 2; iter++) {
        Tstart += (double)indOk[iter + 1];
      }

      if (Tstart > 0.0) {
        /* at least one drivable path */
        Tstart = indOk[0];
        for (iter = 0; iter < 2; iter++) {
          Tstart += (double)indOk[iter + 1];
        }

        if (Tstart == 1.0) {
          /* the drivable one is unique */
          *flag3 = -1.1;
          eml_li_find(indOk, ii_data, ii_size);
          T_size[0] = 1;
          T_size[1] = ii_size[1];
          iter = ii_size[0] * ii_size[1];
          for (i1 = 0; i1 < iter; i1++) {
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

          iter = 1;
          exitg5 = false;
          while ((!exitg5) && (iter < 4)) {
            guard1 = false;
            if (indOk[iter - 1]) {
              idx++;
              ii_data[idx - 1] = iter;
              if (idx >= 3) {
                exitg5 = true;
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              iter++;
            }
          }

          if (1 > idx) {
            i2 = 0;
          } else {
            i2 = idx;
          }

          iter = ii_size[0] * i2;
          for (i1 = 0; i1 < iter; i1++) {
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

    for (iter = 0; iter < iv0[1]; iter++) {
      indOk_data[iter] = rt_powd_snf(T_data[iter], 3.0);
    }

    Tstart = S[0] * 2.0 - S[2] * 2.0;
    iter = iv0[1];
    for (i1 = 0; i1 < iter; i1++) {
      indOk_data[i1] = 1.0 / indOk_data[i1] * (Tstart + T_data[i1] * S[1]) *
        0.25;
    }

    ps[0] = indOk_data[0];
    for (i1 = 0; i1 < 2; i1++) {
      iv0[i1] = (signed char)T_size[i1];
    }

    for (iter = 0; iter < iv0[1]; iter++) {
      indOk_data[iter] = T_data[iter] * T_data[iter];
    }

    Tstart = S[0] * 3.0 - S[2] * 3.0;
    iter = iv0[1];
    for (i1 = 0; i1 < iter; i1++) {
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
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
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
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          iter = ii_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg3 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        iter = indOk_size[1];
        for (i1 = 0; i1 < iter; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg3 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 12.0 * D[0];
    e = 12.0 * D[3];
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
      indOk_data[i1] = -(((D[2] * (T_data[0] * T_data[i1]) + a * T_data[i1]) +
                          Tstart) - e);
    }

    i1 = r4->size[0] * r4->size[1];
    r4->size[0] = 1;
    r4->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r4, i1, (int)sizeof(double));
    iter = c_size[1];
    for (i1 = 0; i1 < iter; i1++) {
      r4->data[i1] = 2.0 * Cs[i1];
    }

    mrdivide(indOk_data, indOk_size, r4->data, r4->size);
    pd[0] = indOk_data[0];
    a = 3.0 * D[2];
    b_a = 16.0 * D[1];
    indOk_size[1] = T_size[1];
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
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
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          iter = ii_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg2 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        iter = indOk_size[1];
        for (i1 = 0; i1 < iter; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg2 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 30.0 * D[0];
    e = 30.0 * D[3];
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
      indOk_data[i1] = ((a * (T_data[0] * T_data[i1]) + b_a * T_data[i1]) +
                        Tstart) - e;
    }

    i1 = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r3, i1, (int)sizeof(double));
    iter = c_size[1];
    for (i1 = 0; i1 < iter; i1++) {
      r3->data[i1] = 2.0 * Cs[i1];
    }

    mrdivide(indOk_data, indOk_size, r3->data, r3->size);
    pd[1] = indOk_data[0];
    a = 3.0 * D[2];
    b_a = 12.0 * D[1];
    indOk_size[1] = T_size[1];
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
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
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Ts[i1] = indOk_data[i1];
          }

          for (i1 = 0; i1 < 2; i1++) {
            ii_size[i1] = c_size[i1];
          }

          c_size[0] = 1;
          c_size[1] = ii_size[1];
          iter = ii_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Ts[ii_size[0] * i1];
          }

          firstmult = false;
        } else {
          Tstart = Cs[0];
          c_size[0] = 1;
          c_size[1] = indOk_size[1];
          iter = indOk_size[1];
          for (i1 = 0; i1 < iter; i1++) {
            Cs[i1] = Tstart * indOk_data[i1];
          }
        }
      }

      if (ed2 == 0.0) {
        exitg1 = 1;
      } else {
        e = ed2;
        Tstart = indOk_data[0];
        iter = indOk_size[1];
        for (i1 = 0; i1 < iter; i1++) {
          indOk_data[i1] *= Tstart;
        }
      }
    } while (exitg1 == 0);

    indOk_size[0] = 1;
    indOk_size[1] = T_size[1];
    Tstart = 20.0 * D[0];
    e = 20.0 * D[3];
    iter = T_size[1];
    for (i1 = 0; i1 < iter; i1++) {
      indOk_data[i1] = -(((a * (T_data[0] * T_data[i1]) + b_a * T_data[i1]) +
                          Tstart) - e);
    }

    i1 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = c_size[1];
    emxEnsureCapacity((emxArray__common *)r2, i1, (int)sizeof(double));
    iter = c_size[1];
    for (i1 = 0; i1 < iter; i1++) {
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
 * Arguments    : const double u[6]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
static void polyder(const double u[6], double a_data[], int a_size[2])
{
  int nlead0;
  int k;
  nlead0 = 0;
  k = 1;
  while ((k < 5) && (u[k - 1] == 0.0)) {
    nlead0++;
    k++;
  }

  a_size[0] = 1;
  a_size[1] = 5 - nlead0;
  for (k = 0; k + 1 <= 5 - nlead0; k++) {
    a_data[k] = u[k + nlead0];
  }

  for (k = 0; k + 1 <= 4 - nlead0; k++) {
    a_data[k] *= (double)(4 - (nlead0 + k)) + 1.0;
  }

  if (!((!rtIsInf(u[5])) && (!rtIsNaN(u[5])))) {
    a_data[4 - nlead0] = rtNaN;
  }
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const double x[100]
 *                double y[100]
 * Return Type  : void
 */
static void polyval(const double p_data[], const int p_size[2], const double x
                    [100], double y[100])
{
  int i3;
  int k;
  if (!(p_size[1] == 0)) {
    for (i3 = 0; i3 < 100; i3++) {
      y[i3] = p_data[0];
    }

    for (k = 0; k <= p_size[1] - 2; k++) {
      for (i3 = 0; i3 < 100; i3++) {
        y[i3] = x[i3] * y[i3] + p_data[k + 1];
      }
    }
  }
}

/*
 * Arguments    : const double a[100]
 *                double y[100]
 * Return Type  : void
 */
static void power(const double a[100], double y[100])
{
  int k;
  for (k = 0; k < 100; k++) {
    y[k] = a[k] * a[k];
  }
}

/*
 * Arguments    : const double x[100]
 *                const double y[100]
 *                double z[100]
 * Return Type  : void
 */
static void rdivide(const double x[100], const double y[100], double z[100])
{
  int i4;
  for (i4 = 0; i4 < 100; i4++) {
    z[i4] = x[i4] / y[i4];
  }
}

/*
 * Arguments    : const emxArray_real_T *a
 *                emxArray_real_T *b
 * Return Type  : void
 */
static void repmat(const emxArray_real_T *a, emxArray_real_T *b)
{
  int unnamed_idx_1;
  int ibmat;
  int itilerow;
  unnamed_idx_1 = b->size[0] * b->size[1];
  b->size[0] = 2;
  b->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)b, unnamed_idx_1, (int)sizeof(double));
  unnamed_idx_1 = a->size[1];
  if (unnamed_idx_1 == 0) {
  } else {
    for (unnamed_idx_1 = 0; unnamed_idx_1 + 1 <= a->size[1]; unnamed_idx_1++) {
      ibmat = unnamed_idx_1 << 1;
      for (itilerow = 0; itilerow < 2; itilerow++) {
        b->data[ibmat + itilerow] = a->data[unnamed_idx_1];
      }
    }
  }
}

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
 * Arguments    : int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : double *data
 *                int numDimensions
 *                int *size
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
  size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : double *data
 *                int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : int rows
 *                int cols
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : emxArray_real_T *emxArray
 * Return Type  : void
 */
void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
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
 *                int maxIter
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
 *                int m
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
 *                double T_data[]
 *                int T_size[2]
 *                double *TOL
 * Return Type  : void
 */
void otg_smart_xy(double absTOL, int maxIter, double v1, double d1, double kj,
                  double kT, double ks, double kd, const double dataVeh[3],
                  double safetyS, double safetyD, double kappaMax, double
                  aOrthMax, int m, double kappa, double b_y0, double phi, double
                  vx0, double ax0, double w, double *flag1, double *flag2,
                  double *flag3, double *flagAll, emxArray_real_T *x,
                  emxArray_real_T *y, double T_data[], int T_size[2], double
                  *TOL)
{
  int ibtile;
  double cosphi_d;
  double y0_d;
  double D[4];
  double dv0[3];
  int pd_size[2];
  double pd_data[36];
  int ps_size[2];
  double ps_data[25];
  emxArray_real_T *tt;
  int k;
  emxArray_real_T *ss;
  emxArray_real_T *dd;
  emxArray_real_T *XY;
  emxArray_real_T *xy_car;
  int jtilecol;
  emxArray_real_T *r0;
  emxArray_real_T *b_ss;
  emxArray_real_T *r1;
  emxArray_real_T *b_x;
  int i0;
  double R[4];
  double dv1[2];

  /* % init */
  ibtile = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = m;
  emxEnsureCapacity((emxArray__common *)x, ibtile, (int)sizeof(double));
  for (ibtile = 0; ibtile < m; ibtile++) {
    x->data[ibtile] = 0.0;
  }

  ibtile = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = m;
  emxEnsureCapacity((emxArray__common *)y, ibtile, (int)sizeof(double));
  for (ibtile = 0; ibtile < m; ibtile++) {
    y->data[ibtile] = 0.0;
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
  if (cosphi_d == -1.0) {
  } else {
    emxInit_real_T(&tt, 2);

    /* % calculate tt */
    ibtile = tt->size[0] * tt->size[1];
    tt->size[0] = 1;
    tt->size[1] = m;
    emxEnsureCapacity((emxArray__common *)tt, ibtile, (int)sizeof(double));
    if (m >= 1) {
      tt->data[m - 1] = T_data[0];
      if (tt->size[1] >= 2) {
        tt->data[0] = 0.0;
        if (tt->size[1] >= 3) {
          if ((T_data[0] < 0.0) && (fabs(T_data[0]) > 8.9884656743115785E+307))
          {
            cosphi_d = T_data[0] / ((double)tt->size[1] - 1.0);
            ibtile = tt->size[1];
            for (k = 0; k <= ibtile - 3; k++) {
              tt->data[k + 1] = cosphi_d * (1.0 + (double)k);
            }
          } else {
            cosphi_d = T_data[0] / ((double)tt->size[1] - 1.0);
            ibtile = tt->size[1];
            for (k = 0; k <= ibtile - 3; k++) {
              tt->data[k + 1] = (1.0 + (double)k) * cosphi_d;
            }
          }
        }
      }
    }

    emxInit_real_T(&ss, 2);
    emxInit_real_T(&dd, 2);

    /* % calc s and d at the times */
    d_polyval(ps_data, ps_size, tt, ss);
    d_polyval(pd_data, pd_size, tt, dd);
    emxInit_real_T(&XY, 2);
    emxInit_real_T(&xy_car, 2);
    if (kappa > 0.0) {
      ibtile = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, ibtile, (int)sizeof(double));
      jtilecol = ss->size[0] * ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        y->data[ibtile] = kappa * ss->data[ibtile];
      }

      emxInit_real_T(&r0, 2);
      ibtile = r0->size[0] * r0->size[1];
      r0->size[0] = 1;
      r0->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)r0, ibtile, (int)sizeof(double));
      jtilecol = y->size[0] * y->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        r0->data[ibtile] = y->data[ibtile];
      }

      for (k = 0; k < y->size[1]; k++) {
        r0->data[k] = sin(r0->data[k]);
      }

      ibtile = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, ibtile, (int)sizeof(double));
      jtilecol = ss->size[0] * ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        y->data[ibtile] = kappa * ss->data[ibtile];
      }

      emxInit_real_T(&b_ss, 2);
      ibtile = b_ss->size[0] * b_ss->size[1];
      b_ss->size[0] = 1;
      b_ss->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)b_ss, ibtile, (int)sizeof(double));
      jtilecol = y->size[0] * y->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        b_ss->data[ibtile] = y->data[ibtile];
      }

      for (k = 0; k < y->size[1]; k++) {
        b_ss->data[k] = cos(b_ss->data[k]);
      }

      ibtile = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, ibtile, (int)sizeof(double));
      jtilecol = ss->size[0] * ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        y->data[ibtile] = kappa * ss->data[ibtile];
      }

      ibtile = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)x, ibtile, (int)sizeof(double));
      jtilecol = y->size[0] * y->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        x->data[ibtile] = y->data[ibtile];
      }

      for (k = 0; k < y->size[1]; k++) {
        x->data[k] = sin(x->data[k]);
      }

      ibtile = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)ss, ibtile, (int)sizeof(double));
      jtilecol = ss->size[0];
      ibtile = ss->size[1];
      jtilecol *= ibtile;
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        ss->data[ibtile] *= kappa;
      }

      ibtile = tt->size[0] * tt->size[1];
      tt->size[0] = 1;
      tt->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)tt, ibtile, (int)sizeof(double));
      jtilecol = ss->size[0] * ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        tt->data[ibtile] = ss->data[ibtile];
      }

      for (k = 0; k < ss->size[1]; k++) {
        tt->data[k] = cos(tt->data[k]);
      }

      emxInit_real_T(&r1, 2);
      repmat(dd, xy_car);
      cosphi_d = 1.0 / kappa;
      ibtile = r1->size[0] * r1->size[1];
      r1->size[0] = 2;
      r1->size[1] = r0->size[1];
      emxEnsureCapacity((emxArray__common *)r1, ibtile, (int)sizeof(double));
      jtilecol = r0->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        r1->data[r1->size[0] * ibtile] = r0->data[r0->size[0] * ibtile];
      }

      emxFree_real_T(&r0);
      jtilecol = b_ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        r1->data[1 + r1->size[0] * ibtile] = 1.0 - b_ss->data[b_ss->size[0] *
          ibtile];
      }

      emxFree_real_T(&b_ss);
      emxInit_real_T(&b_x, 2);
      ibtile = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 2;
      b_x->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)b_x, ibtile, (int)sizeof(double));
      jtilecol = x->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        b_x->data[b_x->size[0] * ibtile] = -x->data[x->size[0] * ibtile];
      }

      jtilecol = tt->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        b_x->data[1 + b_x->size[0] * ibtile] = tt->data[tt->size[0] * ibtile];
      }

      ibtile = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)XY, ibtile, (int)sizeof(double));
      jtilecol = r1->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * ibtile] = cosphi_d * r1->data[i0 +
            r1->size[0] * ibtile] + b_x->data[i0 + b_x->size[0] * ibtile] *
            xy_car->data[i0 + xy_car->size[0] * ibtile];
        }
      }

      emxFree_real_T(&b_x);
      emxFree_real_T(&r1);
    } else {
      ibtile = xy_car->size[0] * xy_car->size[1];
      xy_car->size[0] = 2;
      xy_car->size[1] = m;
      emxEnsureCapacity((emxArray__common *)xy_car, ibtile, (int)sizeof(double));
      if (m == 0) {
      } else {
        for (jtilecol = 1; jtilecol <= m; jtilecol++) {
          ibtile = (jtilecol - 1) << 1;
          for (k = 0; k < 2; k++) {
            xy_car->data[ibtile + k] = ((double)k + 1.0) - 1.0;
          }
        }
      }

      emxInit_real_T(&r0, 2);
      emxInit_real_T(&b_ss, 2);
      repmat(dd, r0);
      ibtile = b_ss->size[0] * b_ss->size[1];
      b_ss->size[0] = 2;
      b_ss->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)b_ss, ibtile, (int)sizeof(double));
      jtilecol = ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        b_ss->data[b_ss->size[0] * ibtile] = ss->data[ss->size[0] * ibtile];
      }

      for (ibtile = 0; ibtile < m; ibtile++) {
        b_ss->data[1 + b_ss->size[0] * ibtile] = 0.0;
      }

      ibtile = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)XY, ibtile, (int)sizeof(double));
      jtilecol = b_ss->size[1];
      for (ibtile = 0; ibtile < jtilecol; ibtile++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * ibtile] = b_ss->data[i0 + b_ss->size[0] *
            ibtile] + xy_car->data[i0 + xy_car->size[0] * ibtile] * r0->data[i0
            + r0->size[0] * ibtile];
        }
      }

      emxFree_real_T(&b_ss);
      emxFree_real_T(&r0);
    }

    emxFree_real_T(&dd);
    emxFree_real_T(&ss);
    emxFree_real_T(&tt);
    ibtile = xy_car->size[0] * xy_car->size[1];
    xy_car->size[0] = 2;
    xy_car->size[1] = m;
    emxEnsureCapacity((emxArray__common *)xy_car, ibtile, (int)sizeof(double));
    jtilecol = m << 1;
    for (ibtile = 0; ibtile < jtilecol; ibtile++) {
      xy_car->data[ibtile] = 0.0;
    }

    R[0] = cos(phi);
    R[2] = -sin(phi);
    R[1] = sin(phi);
    R[3] = cos(phi);
    for (k = 0; k + 1 <= m; k++) {
      dv1[0] = -sin(phi);
      dv1[1] = cos(phi);
      for (ibtile = 0; ibtile < 2; ibtile++) {
        cosphi_d = 0.0;
        for (i0 = 0; i0 < 2; i0++) {
          cosphi_d += R[ibtile + (i0 << 1)] * XY->data[i0 + XY->size[0] * k];
        }

        xy_car->data[ibtile + xy_car->size[0] * k] = cosphi_d + dv1[ibtile] *
          -D[0];
      }
    }

    emxFree_real_T(&XY);
    jtilecol = xy_car->size[1];
    ibtile = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = jtilecol;
    emxEnsureCapacity((emxArray__common *)x, ibtile, (int)sizeof(double));
    for (ibtile = 0; ibtile < jtilecol; ibtile++) {
      x->data[x->size[0] * ibtile] = xy_car->data[xy_car->size[0] * ibtile];
    }

    jtilecol = xy_car->size[1];
    ibtile = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = jtilecol;
    emxEnsureCapacity((emxArray__common *)y, ibtile, (int)sizeof(double));
    for (ibtile = 0; ibtile < jtilecol; ibtile++) {
      y->data[y->size[0] * ibtile] = xy_car->data[1 + xy_car->size[0] * ibtile];
    }

    emxFree_real_T(&xy_car);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void otg_smart_xy_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void otg_smart_xy_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for otg_smart_xy.c
 *
 * [EOF]
 */
