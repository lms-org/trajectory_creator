/*
 * File: otg_smart_objFun.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_smart_objFun.h"
#include "otg_pd.h"
#include "otg_ps.h"
#include "otg_smart_pspdT.h"
#include "sqrt.h"
#include "polyval.h"
#include "roots.h"
#include "polyder.h"
#include "otg_smart_xy_rtwutil.h"

/* Function Declarations */
static void b_eml_li_find(const boolean_T x_data[], const int x_size[1], int
  y_data[], int y_size[1]);
static boolean_T b_eml_relop(const creal_T a, const creal_T b);
static void eml_li_find(const boolean_T x_data[], const int x_size[2], int
  y_data[], int y_size[2]);
static boolean_T eml_relop(const creal_T a, const creal_T b);
static double rt_atan2d_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : const boolean_T x_data[]
 *                const int x_size[1]
 *                int y_data[]
 *                int y_size[1]
 * Return Type  : void
 */
static void b_eml_li_find(const boolean_T x_data[], const int x_size[1], int
  y_data[], int y_size[1])
{
  int n;
  int k;
  int i;
  n = x_size[0];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      y_data[k] = i;
      k++;
    }
  }
}

/*
 * Arguments    : const creal_T a
 *                const creal_T b
 * Return Type  : boolean_T
 */
static boolean_T b_eml_relop(const creal_T a, const creal_T b)
{
  boolean_T p;
  double absbi;
  double y;
  double absxk;
  int exponent;
  double absar;
  double absbr;
  double Ma;
  int b_exponent;
  int c_exponent;
  int d_exponent;
  if ((fabs(a.re) > 8.9884656743115785E+307) || (fabs(a.im) >
       8.9884656743115785E+307) || (fabs(b.re) > 8.9884656743115785E+307) ||
      (fabs(b.im) > 8.9884656743115785E+307)) {
    absbi = rt_hypotd_snf(a.re / 2.0, a.im / 2.0);
    y = rt_hypotd_snf(b.re / 2.0, b.im / 2.0);
  } else {
    absbi = rt_hypotd_snf(a.re, a.im);
    y = rt_hypotd_snf(b.re, b.im);
  }

  absxk = y / 2.0;
  if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      absxk = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      absxk = ldexp(1.0, exponent - 53);
    }
  } else {
    absxk = rtNaN;
  }

  if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
         0.0) == (y > 0.0)))) {
    p = true;
  } else {
    p = false;
  }

  if (p) {
    absar = fabs(a.re);
    absxk = fabs(a.im);
    absbr = fabs(b.re);
    absbi = fabs(b.im);
    if (absar > absxk) {
      Ma = absar;
      absar = absxk;
    } else {
      Ma = absxk;
    }

    if (absbr > absbi) {
      absxk = absbr;
      absbr = absbi;
    } else {
      absxk = absbi;
    }

    if (Ma > absxk) {
      if (absar < absbr) {
        absbi = Ma - absxk;
        y = (absar / 2.0 + absbr / 2.0) / (Ma / 2.0 + absxk / 2.0) * (absbr -
          absar);
      } else {
        absbi = Ma;
        y = absxk;
      }
    } else if (Ma < absxk) {
      if (absar > absbr) {
        y = absxk - Ma;
        absbi = (absar / 2.0 + absbr / 2.0) / (Ma / 2.0 + absxk / 2.0) * (absar
          - absbr);
      } else {
        absbi = Ma;
        y = absxk;
      }
    } else {
      absbi = absar;
      y = absbr;
    }

    absxk = fabs(y / 2.0);
    if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
      if (absxk <= 2.2250738585072014E-308) {
        absxk = 4.94065645841247E-324;
      } else {
        frexp(absxk, &b_exponent);
        absxk = ldexp(1.0, b_exponent - 53);
      }
    } else {
      absxk = rtNaN;
    }

    if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
           0.0) == (y > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      absbi = rt_atan2d_snf(a.im, a.re);
      y = rt_atan2d_snf(b.im, b.re);
      absxk = fabs(y / 2.0);
      if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
        if (absxk <= 2.2250738585072014E-308) {
          absxk = 4.94065645841247E-324;
        } else {
          frexp(absxk, &c_exponent);
          absxk = ldexp(1.0, c_exponent - 53);
        }
      } else {
        absxk = rtNaN;
      }

      if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
             0.0) == (y > 0.0)))) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        if (absbi > 0.78539816339744828) {
          if (absbi > 2.3561944901923448) {
            absbi = -a.im;
            y = -b.im;
          } else {
            absbi = -a.re;
            y = -b.re;
          }
        } else if (absbi > -0.78539816339744828) {
          absbi = a.im;
          y = b.im;
        } else if (absbi > -2.3561944901923448) {
          absbi = a.re;
          y = b.re;
        } else {
          absbi = -a.im;
          y = -b.im;
        }

        absxk = fabs(y / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
          if (absxk <= 2.2250738585072014E-308) {
            absxk = 4.94065645841247E-324;
          } else {
            frexp(absxk, &d_exponent);
            absxk = ldexp(1.0, d_exponent - 53);
          }
        } else {
          absxk = rtNaN;
        }

        if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) &&
             ((absbi > 0.0) == (y > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          absbi = 0.0;
          y = 0.0;
        }
      }
    }
  }

  return absbi < y;
}

/*
 * Arguments    : const boolean_T x_data[]
 *                const int x_size[2]
 *                int y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
static void eml_li_find(const boolean_T x_data[], const int x_size[2], int
  y_data[], int y_size[2])
{
  int n;
  int k;
  int i;
  n = x_size[1];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      k++;
    }
  }

  y_size[0] = 1;
  y_size[1] = k;
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      y_data[k] = i;
      k++;
    }
  }
}

/*
 * Arguments    : const creal_T a
 *                const creal_T b
 * Return Type  : boolean_T
 */
static boolean_T eml_relop(const creal_T a, const creal_T b)
{
  boolean_T p;
  double absbi;
  double y;
  double absxk;
  int exponent;
  double absar;
  double absbr;
  double Ma;
  int b_exponent;
  int c_exponent;
  int d_exponent;
  if ((fabs(a.re) > 8.9884656743115785E+307) || (fabs(a.im) >
       8.9884656743115785E+307) || (fabs(b.re) > 8.9884656743115785E+307) ||
      (fabs(b.im) > 8.9884656743115785E+307)) {
    absbi = rt_hypotd_snf(a.re / 2.0, a.im / 2.0);
    y = rt_hypotd_snf(b.re / 2.0, b.im / 2.0);
  } else {
    absbi = rt_hypotd_snf(a.re, a.im);
    y = rt_hypotd_snf(b.re, b.im);
  }

  absxk = y / 2.0;
  if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      absxk = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      absxk = ldexp(1.0, exponent - 53);
    }
  } else {
    absxk = rtNaN;
  }

  if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
         0.0) == (y > 0.0)))) {
    p = true;
  } else {
    p = false;
  }

  if (p) {
    absar = fabs(a.re);
    absxk = fabs(a.im);
    absbr = fabs(b.re);
    absbi = fabs(b.im);
    if (absar > absxk) {
      Ma = absar;
      absar = absxk;
    } else {
      Ma = absxk;
    }

    if (absbr > absbi) {
      absxk = absbr;
      absbr = absbi;
    } else {
      absxk = absbi;
    }

    if (Ma > absxk) {
      if (absar < absbr) {
        absbi = Ma - absxk;
        y = (absar / 2.0 + absbr / 2.0) / (Ma / 2.0 + absxk / 2.0) * (absbr -
          absar);
      } else {
        absbi = Ma;
        y = absxk;
      }
    } else if (Ma < absxk) {
      if (absar > absbr) {
        y = absxk - Ma;
        absbi = (absar / 2.0 + absbr / 2.0) / (Ma / 2.0 + absxk / 2.0) * (absar
          - absbr);
      } else {
        absbi = Ma;
        y = absxk;
      }
    } else {
      absbi = absar;
      y = absbr;
    }

    absxk = fabs(y / 2.0);
    if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
      if (absxk <= 2.2250738585072014E-308) {
        absxk = 4.94065645841247E-324;
      } else {
        frexp(absxk, &b_exponent);
        absxk = ldexp(1.0, b_exponent - 53);
      }
    } else {
      absxk = rtNaN;
    }

    if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
           0.0) == (y > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      absbi = rt_atan2d_snf(a.im, a.re);
      y = rt_atan2d_snf(b.im, b.re);
      absxk = fabs(y / 2.0);
      if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
        if (absxk <= 2.2250738585072014E-308) {
          absxk = 4.94065645841247E-324;
        } else {
          frexp(absxk, &c_exponent);
          absxk = ldexp(1.0, c_exponent - 53);
        }
      } else {
        absxk = rtNaN;
      }

      if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) && ((absbi >
             0.0) == (y > 0.0)))) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        if (absbi > 0.78539816339744828) {
          if (absbi > 2.3561944901923448) {
            absbi = -a.im;
            y = -b.im;
          } else {
            absbi = -a.re;
            y = -b.re;
          }
        } else if (absbi > -0.78539816339744828) {
          absbi = a.im;
          y = b.im;
        } else if (absbi > -2.3561944901923448) {
          absbi = a.re;
          y = b.re;
        } else {
          absbi = -a.im;
          y = -b.im;
        }

        absxk = fabs(y / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
          if (absxk <= 2.2250738585072014E-308) {
            absxk = 4.94065645841247E-324;
          } else {
            frexp(absxk, &d_exponent);
            absxk = ldexp(1.0, d_exponent - 53);
          }
        } else {
          absxk = rtNaN;
        }

        if ((fabs(y - absbi) < absxk) || (rtIsInf(absbi) && rtIsInf(y) &&
             ((absbi > 0.0) == (y > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          absbi = 0.0;
          y = 0.0;
        }
      }
    }
  }

  return absbi > y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(double)(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
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
 *                double *flag
 * Return Type  : void
 */
void b_otg_smart_objFun(double T, const double S[3], const double D[4], double
  kj, double kT, double ks, double kd, const double dataVeh[3], double safetyS,
  double safetyD, double kappa, double kappaMax, double aOrthMax, double *Ctot,
  double *notD, double *coll, double *flag)
{
  double d_d_data[5];
  double PS[5];
  int ixstart;
  double dv4[6];
  double PD[6];
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int tmp_size[2];
  double tmp_data[4];
  int d_ddd_size[2];
  int n;
  double d_ddd_data[3];
  int b_tmp_size[1];
  creal_T b_tmp_data[3];
  creal_T r_data[3];
  boolean_T b_r_data[2];
  int r_size[2];
  int c_tmp_data[3];
  int b_r_size[2];
  creal_T c_r_data[4];
  int c_tmp_size[2];
  creal_T d_dd_ts_data[5];
  creal_T mtmp;
  int ix;
  boolean_T exitg4;
  double s_d_data[4];
  int c_r_size[2];
  int d_r_size[2];
  int d_tmp_size[2];
  creal_T vx_max;
  boolean_T exitg3;
  creal_T b_mtmp;
  boolean_T exitg2;
  double mtmp_re;
  double d;
  double brm;
  double kappa_xy_max_re;
  boolean_T d_r_data[3];
  int e_r_size[2];
  int f_r_size[2];
  creal_T d_tmp_data[5];
  int e_tmp_size[2];
  boolean_T exitg1;
  double abstandMinusSafety_poly[5];
  boolean_T e_r_data[4];
  int g_r_size[1];
  int e_tmp_data[4];
  creal_T f_r_data[4];
  int h_r_size[1];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;

  /* % init */
  *flag = 1.0;
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
  /* LENKEINSCHLAG */
  polyder(PD, d_d_data, d_d_size);

  /* first derivative w.r.t. t of d(t) */
  b_polyder(d_d_data, d_d_size, d_dd_data, d_dd_size);

  /* second derivative w.r.t. t of d(t) */
  b_polyder(d_dd_data, d_dd_size, tmp_data, tmp_size);
  d_ddd_size[0] = 1;
  d_ddd_size[1] = tmp_size[1];
  n = tmp_size[0] * tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    d_ddd_data[ixstart] = tmp_data[ixstart];
  }

  /* third derivative to find the max of the second derivative */
  roots(d_ddd_data, d_ddd_size, b_tmp_data, b_tmp_size);
  n = b_tmp_size[0];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart].re = b_tmp_data[ixstart].re;
    r_data[ixstart].im = -b_tmp_data[ixstart].im;
  }

  /* get all roots of THIRD derivative */
  /* only real roots */
  r_size[0] = 1;
  r_size[1] = b_tmp_size[0];
  n = b_tmp_size[0];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_r_data[ixstart] = (r_data[ixstart].im == 0.0);
  }

  eml_li_find(b_r_data, r_size, c_tmp_data, tmp_size);
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
  }

  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart] = b_tmp_data[ixstart];
  }

  /* only roots in interval */
  b_r_size[0] = 1;
  b_r_size[1] = tmp_size[1];
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re < T));
  }

  eml_li_find(b_r_data, b_r_size, c_tmp_data, tmp_size);
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
  }

  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart] = b_tmp_data[ixstart];
  }

  /* also max at the border */
  c_tmp_size[0] = 1;
  c_tmp_size[1] = 2 + tmp_size[1];
  c_r_data[0].re = 0.0;
  c_r_data[0].im = 0.0;
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    c_r_data[ixstart + 1] = r_data[ixstart];
  }

  c_r_data[1 + tmp_size[1]].re = T;
  c_r_data[1 + tmp_size[1]].im = 0.0;
  polyval(d_dd_data, d_dd_size, c_r_data, c_tmp_size, d_dd_ts_data, r_size);

  /* HERE ONE NEEDS THE SECOND DERIVATIVE */
  ixstart = 1;
  n = r_size[1];
  mtmp = d_dd_ts_data[0];
  if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
    ix = 1;
    exitg4 = false;
    while ((!exitg4) && (ix + 1 <= n)) {
      ixstart = ix + 1;
      if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
        mtmp = d_dd_ts_data[ix];
        exitg4 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < r_size[1]) {
    while (ixstart + 1 <= n) {
      if (eml_relop(d_dd_ts_data[ixstart], mtmp)) {
        mtmp = d_dd_ts_data[ixstart];
      }

      ixstart++;
    }
  }

  /* get the maximal value of the polynomial on the interval */
  c_polyder(PS, s_d_data, b_r_size);

  /* first derivative w.r.t. t of s(t) = v(t) */
  b_polyder(s_d_data, b_r_size, tmp_data, tmp_size);
  d_ddd_size[0] = 1;
  d_ddd_size[1] = tmp_size[1];
  n = tmp_size[0] * tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    d_ddd_data[ixstart] = tmp_data[ixstart];
  }

  /* second derivative to find the */
  roots(d_ddd_data, d_ddd_size, b_tmp_data, b_tmp_size);
  n = b_tmp_size[0];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart].re = b_tmp_data[ixstart].re;
    r_data[ixstart].im = -b_tmp_data[ixstart].im;
  }

  /* get all roots of first derivative */
  /* only real roots */
  c_r_size[0] = 1;
  c_r_size[1] = b_tmp_size[0];
  n = b_tmp_size[0];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_r_data[ixstart] = (r_data[ixstart].im == 0.0);
  }

  eml_li_find(b_r_data, c_r_size, c_tmp_data, tmp_size);
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
  }

  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart] = b_tmp_data[ixstart];
  }

  /* only roots in interval */
  d_r_size[0] = 1;
  d_r_size[1] = tmp_size[1];
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re < T));
  }

  eml_li_find(b_r_data, d_r_size, c_tmp_data, tmp_size);
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
  }

  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    r_data[ixstart] = b_tmp_data[ixstart];
  }

  /* also min at the border */
  d_tmp_size[0] = 1;
  d_tmp_size[1] = 2 + tmp_size[1];
  c_r_data[0].re = 0.0;
  c_r_data[0].im = 0.0;
  n = tmp_size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    c_r_data[ixstart + 1] = r_data[ixstart];
  }

  c_r_data[1 + tmp_size[1]].re = T;
  c_r_data[1 + tmp_size[1]].im = 0.0;
  polyval(s_d_data, b_r_size, c_r_data, d_tmp_size, d_dd_ts_data, r_size);
  ixstart = 1;
  n = r_size[1];
  vx_max = d_dd_ts_data[0];
  if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
    ix = 1;
    exitg3 = false;
    while ((!exitg3) && (ix + 1 <= n)) {
      ixstart = ix + 1;
      if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
        vx_max = d_dd_ts_data[ix];
        exitg3 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < r_size[1]) {
    while (ixstart + 1 <= n) {
      if (b_eml_relop(d_dd_ts_data[ixstart], vx_max)) {
        vx_max = d_dd_ts_data[ixstart];
      }

      ixstart++;
    }
  }

  /* get the min value of the polynomial on the interval */
  ixstart = 1;
  n = r_size[1];
  b_mtmp = d_dd_ts_data[0];
  if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
    ix = 1;
    exitg2 = false;
    while ((!exitg2) && (ix + 1 <= n)) {
      ixstart = ix + 1;
      if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
        b_mtmp = d_dd_ts_data[ix];
        exitg2 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < r_size[1]) {
    while (ixstart + 1 <= n) {
      if (eml_relop(d_dd_ts_data[ixstart], b_mtmp)) {
        b_mtmp = d_dd_ts_data[ixstart];
      }

      ixstart++;
    }
  }

  /* get the maximal value of the velocity in s direction */
  if (vx_max.re <= 0.0) {
    /* something went terrible wrong. we are driving backwards */
    *flag = -1.0;
  } else {
    /* now the bound for the curvature can be calculated */
    if (vx_max.im == 0.0) {
      if (mtmp.im == 0.0) {
        mtmp_re = mtmp.re / vx_max.re;
        d = 0.0;
      } else if (mtmp.re == 0.0) {
        mtmp_re = 0.0;
        d = mtmp.im / vx_max.re;
      } else {
        mtmp_re = mtmp.re / vx_max.re;
        d = mtmp.im / vx_max.re;
      }
    } else if (vx_max.re == 0.0) {
      if (mtmp.re == 0.0) {
        mtmp_re = mtmp.im / vx_max.im;
        d = 0.0;
      } else if (mtmp.im == 0.0) {
        mtmp_re = 0.0;
        d = -(mtmp.re / vx_max.im);
      } else {
        mtmp_re = mtmp.im / vx_max.im;
        d = -(mtmp.re / vx_max.im);
      }
    } else {
      brm = vx_max.re;
      kappa_xy_max_re = fabs(vx_max.im);
      if (brm > kappa_xy_max_re) {
        kappa_xy_max_re = vx_max.im / vx_max.re;
        d = vx_max.re + kappa_xy_max_re * vx_max.im;
        mtmp_re = (mtmp.re + kappa_xy_max_re * mtmp.im) / d;
        d = (mtmp.im - kappa_xy_max_re * mtmp.re) / d;
      } else if (kappa_xy_max_re == brm) {
        if (vx_max.re > 0.0) {
          kappa_xy_max_re = 0.5;
        } else {
          kappa_xy_max_re = -0.5;
        }

        if (vx_max.im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        mtmp_re = (mtmp.re * kappa_xy_max_re + mtmp.im * d) / brm;
        d = (mtmp.im * kappa_xy_max_re - mtmp.re * d) / brm;
      } else {
        kappa_xy_max_re = vx_max.re / vx_max.im;
        d = vx_max.im + kappa_xy_max_re * vx_max.re;
        mtmp_re = (kappa_xy_max_re * mtmp.re + mtmp.im) / d;
        d = (kappa_xy_max_re * mtmp.im - mtmp.re) / d;
      }
    }

    kappa_xy_max_re = fabs(kappa) + mtmp_re;
    if (kappa_xy_max_re > kappaMax) {
      /* the curvature could be (according to our approximation) to big */
      *notD = 1.0;
    }

    /* MAXIMALE ORTH: BESCHLEUNIGUNG */
    roots(d_dd_data, d_dd_size, b_tmp_data, b_tmp_size);
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      r_data[ixstart].re = b_tmp_data[ixstart].re;
      r_data[ixstart].im = -b_tmp_data[ixstart].im;
    }

    /* get all roots of SECOND derivative */
    /* only real roots */
    e_r_size[0] = 1;
    e_r_size[1] = b_tmp_size[0];
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      d_r_data[ixstart] = (r_data[ixstart].im == 0.0);
    }

    eml_li_find(d_r_data, e_r_size, c_tmp_data, tmp_size);
    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
    }

    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      r_data[ixstart] = b_tmp_data[ixstart];
    }

    /* only roots in interval */
    f_r_size[0] = 1;
    f_r_size[1] = tmp_size[1];
    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      d_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re < T));
    }

    eml_li_find(d_r_data, f_r_size, c_tmp_data, tmp_size);
    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
    }

    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      r_data[ixstart] = b_tmp_data[ixstart];
    }

    /* also max at the border */
    e_tmp_size[0] = 1;
    e_tmp_size[1] = 2 + tmp_size[1];
    d_tmp_data[0].re = 0.0;
    d_tmp_data[0].im = 0.0;
    n = tmp_size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      d_tmp_data[ixstart + 1] = r_data[ixstart];
    }

    d_tmp_data[1 + tmp_size[1]].re = T;
    d_tmp_data[1 + tmp_size[1]].im = 0.0;
    polyval(d_d_data, d_d_size, d_tmp_data, e_tmp_size, d_dd_ts_data, r_size);

    /* HERE ONE NEEDS THE SECOND DERIVATIVE */
    ixstart = 1;
    n = r_size[1];
    mtmp = d_dd_ts_data[0];
    if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
      ix = 1;
      exitg1 = false;
      while ((!exitg1) && (ix + 1 <= n)) {
        ixstart = ix + 1;
        if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
          mtmp = d_dd_ts_data[ix];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < r_size[1]) {
      while (ixstart + 1 <= n) {
        if (eml_relop(d_dd_ts_data[ixstart], mtmp)) {
          mtmp = d_dd_ts_data[ixstart];
        }

        ixstart++;
      }
    }

    /* get the maximal value of the polynomial on the interval */
    vx_max.re = (b_mtmp.re * b_mtmp.re - b_mtmp.im * b_mtmp.im) + (mtmp.re *
      mtmp.re - mtmp.im * mtmp.im);
    vx_max.im = (b_mtmp.re * b_mtmp.im + b_mtmp.im * b_mtmp.re) + (mtmp.re *
      mtmp.im + mtmp.im * mtmp.re);
    b_sqrt(&vx_max);
    if ((vx_max.re * vx_max.re - vx_max.im * vx_max.im) * kappa_xy_max_re -
        (vx_max.re * vx_max.im + vx_max.im * vx_max.re) * d > aOrthMax) {
      /* the orthogonal acceleration is too much for the grip */
      *notD += 2.0;

      /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
    }

    /* % check for collision */
    /* s(t) should be an increasing function. This is guranteed because we */
    /* calculate the min of the 1st derivative above and throw an error if */
    /* this is non positive */
    abstandMinusSafety_poly[0] = 0.0 - PS[0];
    abstandMinusSafety_poly[1] = 0.0 - PS[1];
    abstandMinusSafety_poly[2] = 0.0 - PS[2];
    abstandMinusSafety_poly[3] = dataVeh[1] - PS[3];
    abstandMinusSafety_poly[4] = (dataVeh[0] - safetyS) - PS[4];
    b_roots(abstandMinusSafety_poly, c_r_data, b_tmp_size);

    /* only real roots */
    g_r_size[0] = b_tmp_size[0];
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      e_r_data[ixstart] = (c_r_data[ixstart].im == 0.0);
    }

    b_eml_li_find(e_r_data, g_r_size, e_tmp_data, b_tmp_size);
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      f_r_data[ixstart] = c_r_data[e_tmp_data[ixstart] - 1];
    }

    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      c_r_data[ixstart] = f_r_data[ixstart];
    }

    /* only roots in interval */
    h_r_size[0] = b_tmp_size[0];
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      e_r_data[ixstart] = ((c_r_data[ixstart].re > 0.0) && (c_r_data[ixstart].re
        < T));
    }

    b_eml_li_find(e_r_data, h_r_size, e_tmp_data, b_tmp_size);
    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      f_r_data[ixstart] = c_r_data[e_tmp_data[ixstart] - 1];
    }

    n = b_tmp_size[0];
    for (ixstart = 0; ixstart < n; ixstart++) {
      c_r_data[ixstart] = f_r_data[ixstart];
    }

    guard1 = false;
    guard2 = false;
    guard3 = false;
    guard4 = false;
    if (b_tmp_size[0] == 0) {
      kappa_xy_max_re = T / 2.0;
      d = 0.0 - PS[0];
      for (ixstart = 0; ixstart < 4; ixstart++) {
        d = kappa_xy_max_re * d + abstandMinusSafety_poly[ixstart + 1];
      }

      if (d > 0.0) {
        /* all good no collsion possible */
      } else {
        guard4 = true;
      }
    } else {
      guard4 = true;
    }

    if (guard4) {
      if (b_tmp_size[0] == 0) {
        kappa_xy_max_re = T / 2.0;
        d = 0.0 - PS[0];
        for (ixstart = 0; ixstart < 4; ixstart++) {
          d = kappa_xy_max_re * d + abstandMinusSafety_poly[ixstart + 1];
        }

        if (d <= 0.0) {
          /* fix collsion */
          *coll = 1.0;
        } else {
          guard3 = true;
        }
      } else {
        guard3 = true;
      }
    }

    if (guard3) {
      if (b_tmp_size[0] == 1) {
        if (c_r_data[0].im == 0.0) {
          vx_max.re = c_r_data[0].re / 2.0;
          vx_max.im = 0.0;
        } else if (c_r_data[0].re == 0.0) {
          vx_max.re = 0.0;
          vx_max.im = c_r_data[0].im / 2.0;
        } else {
          vx_max.re = c_r_data[0].re / 2.0;
          vx_max.im = c_r_data[0].im / 2.0;
        }

        mtmp.re = 0.0 - PS[0];
        mtmp.im = 0.0;
        for (ixstart = 0; ixstart < 4; ixstart++) {
          kappa_xy_max_re = vx_max.re * mtmp.im + vx_max.im * mtmp.re;
          mtmp.re = (vx_max.re * mtmp.re - vx_max.im * mtmp.im) +
            abstandMinusSafety_poly[ixstart + 1];
          mtmp.im = kappa_xy_max_re;
        }

        if (mtmp.re > 0.0) {
          /* one unique possible collision possible */
          vx_max.re = PD[0];
          vx_max.im = 0.0;
          for (ixstart = 0; ixstart < 5; ixstart++) {
            kappa_xy_max_re = c_r_data[0].re * vx_max.im + c_r_data[0].im *
              vx_max.re;
            vx_max.re = (c_r_data[0].re * vx_max.re - c_r_data[0].im * vx_max.im)
              + PD[ixstart + 1];
            vx_max.im = kappa_xy_max_re;
          }

          if ((dataVeh[2] == -1.0) && (vx_max.re < safetyD)) {
            *coll = 1.0;
          }

          if ((dataVeh[2] == 1.0) && (vx_max.re > safetyD)) {
            *coll = 1.0;
          }
        } else {
          guard2 = true;
        }
      } else {
        guard2 = true;
      }
    }

    if (guard2) {
      if (b_tmp_size[0] == 1) {
        if (c_r_data[0].im == 0.0) {
          vx_max.re = c_r_data[0].re / 2.0;
          vx_max.im = 0.0;
        } else if (c_r_data[0].re == 0.0) {
          vx_max.re = 0.0;
          vx_max.im = c_r_data[0].im / 2.0;
        } else {
          vx_max.re = c_r_data[0].re / 2.0;
          vx_max.im = c_r_data[0].im / 2.0;
        }

        mtmp.re = 0.0 - PS[0];
        mtmp.im = 0.0;
        for (ixstart = 0; ixstart < 4; ixstart++) {
          kappa_xy_max_re = vx_max.re * mtmp.im + vx_max.im * mtmp.re;
          mtmp.re = (vx_max.re * mtmp.re - vx_max.im * mtmp.im) +
            abstandMinusSafety_poly[ixstart + 1];
          mtmp.im = kappa_xy_max_re;
        }

        if (mtmp.re > 0.0) {
          /* there is a collision directly at the start */
          *coll = -1.0;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
    }

    if (guard1) {
      *flag = -3.0;
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
 *                double *flag
 * Return Type  : void
 */
void c_otg_smart_objFun(const double T_data[], const int T_size[2], const double
  S[3], const double D[4], double kj, double kT, double ks, double kd, const
  double dataVeh[3], double safetyS, double safetyD, double kappa, double
  kappaMax, double aOrthMax, double Ctot_data[], int Ctot_size[2], double
  notD_data[], int notD_size[2], double coll_data[], int coll_size[2], double
  *flag)
{
  int n;
  int ixstart;
  int b_n;
  double PS_data[15];
  double PD_data[18];
  int i;
  double abstandMinusSafety_poly[5];
  double dv5[6];
  boolean_T guard1 = false;
  int32_T exitg2;
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int tmp_size[2];
  double tmp_data[4];
  int d_ddd_size[2];
  double d_ddd_data[3];
  int r_size[1];
  creal_T b_tmp_data[3];
  creal_T r_data[3];
  boolean_T b_r_data[2];
  int b_r_size[2];
  int c_tmp_data[3];
  int c_r_size[2];
  creal_T c_r_data[4];
  int b_tmp_size[2];
  creal_T d_dd_ts_data[5];
  creal_T y;
  int ix;
  boolean_T exitg6;
  int s_d_size[2];
  double s_d_data[4];
  int d_r_size[2];
  int e_r_size[2];
  int c_tmp_size[2];
  creal_T d_r;
  boolean_T exitg5;
  creal_T mtmp;
  boolean_T exitg4;
  double y_re;
  double b_y;
  double brm;
  double kappa_xy_max_re;
  boolean_T d_r_data[3];
  int f_r_size[2];
  int g_r_size[2];
  creal_T d_tmp_data[5];
  int d_tmp_size[2];
  boolean_T exitg3;
  boolean_T exitg1;
  double dv6[5];
  boolean_T e_r_data[4];
  int h_r_size[1];
  int e_tmp_data[4];
  creal_T f_r_data[4];
  int i_r_size[1];
  boolean_T b_guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;

  /* % init */
  n = T_size[1] - 1;
  *flag = 1.0;
  Ctot_size[0] = 1;
  Ctot_size[1] = T_size[1];
  ixstart = T_size[1];
  for (b_n = 0; b_n < ixstart; b_n++) {
    Ctot_data[b_n] = 0.0;
  }

  notD_size[0] = 1;
  notD_size[1] = T_size[1];
  ixstart = T_size[1];
  for (b_n = 0; b_n < ixstart; b_n++) {
    notD_data[b_n] = 0.0;
  }

  coll_size[0] = 1;
  coll_size[1] = T_size[1];
  ixstart = T_size[1];
  for (b_n = 0; b_n < ixstart; b_n++) {
    coll_data[b_n] = 0.0;
  }

  /* % Get the values of ps, pd for the Ts */
  ixstart = 5 * T_size[1];
  for (b_n = 0; b_n < ixstart; b_n++) {
    PS_data[b_n] = 0.0;
  }

  ixstart = 6 * T_size[1];
  for (b_n = 0; b_n < ixstart; b_n++) {
    PD_data[b_n] = 0.0;
  }

  for (i = 0; i <= n; i++) {
    otg_ps(S, T_data[i], abstandMinusSafety_poly);
    for (b_n = 0; b_n < 5; b_n++) {
      PS_data[b_n + 5 * i] = abstandMinusSafety_poly[b_n];
    }

    otg_pd(D, T_data[i], dv5);
    for (b_n = 0; b_n < 6; b_n++) {
      PD_data[b_n + 6 * i] = dv5[b_n];
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
  i = 0;
  guard1 = false;
  do {
    exitg2 = 0;
    if (i <= n) {
      /* LENKEINSCHLAG */
      polyder(*(double (*)[6])&PD_data[6 * i], abstandMinusSafety_poly, d_d_size);

      /* first derivative w.r.t. t of d(t) */
      b_polyder(abstandMinusSafety_poly, d_d_size, d_dd_data, d_dd_size);

      /* second derivative w.r.t. t of d(t) */
      b_polyder(d_dd_data, d_dd_size, tmp_data, tmp_size);
      d_ddd_size[0] = 1;
      d_ddd_size[1] = tmp_size[1];
      ixstart = tmp_size[0] * tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        d_ddd_data[b_n] = tmp_data[b_n];
      }

      /* third derivative to find the max of the second derivative */
      roots(d_ddd_data, d_ddd_size, b_tmp_data, r_size);
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n].re = b_tmp_data[b_n].re;
        r_data[b_n].im = -b_tmp_data[b_n].im;
      }

      /* get all roots of THIRD derivative */
      /* only real roots */
      b_r_size[0] = 1;
      b_r_size[1] = r_size[0];
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_r_data[b_n] = (r_data[b_n].im == 0.0);
      }

      eml_li_find(b_r_data, b_r_size, c_tmp_data, tmp_size);
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
      }

      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n] = b_tmp_data[b_n];
      }

      /* only roots in interval */
      c_r_size[0] = 1;
      c_r_size[1] = tmp_size[1];
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_r_data[b_n] = ((r_data[b_n].re > 0.0) && (r_data[b_n].re < T_data[i]));
      }

      eml_li_find(b_r_data, c_r_size, c_tmp_data, tmp_size);
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
      }

      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n] = b_tmp_data[b_n];
      }

      /* also max at the border */
      b_tmp_size[0] = 1;
      b_tmp_size[1] = 2 + tmp_size[1];
      c_r_data[0].re = 0.0;
      c_r_data[0].im = 0.0;
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        c_r_data[b_n + 1] = r_data[b_n];
      }

      c_r_data[1 + tmp_size[1]].re = T_data[i];
      c_r_data[1 + tmp_size[1]].im = 0.0;
      polyval(d_dd_data, d_dd_size, c_r_data, b_tmp_size, d_dd_ts_data, tmp_size);

      /* HERE ONE NEEDS THE SECOND DERIVATIVE */
      ixstart = 1;
      b_n = tmp_size[1];
      y = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg6 = false;
        while ((!exitg6) && (ix + 1 <= b_n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            y = d_dd_ts_data[ix];
            exitg6 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= b_n) {
          if (eml_relop(d_dd_ts_data[ixstart], y)) {
            y = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the maximal value of the polynomial on the interval */
      c_polyder(*(double (*)[5])&PS_data[5 * i], s_d_data, s_d_size);

      /* first derivative w.r.t. t of s(t) = v(t) */
      b_polyder(s_d_data, s_d_size, tmp_data, tmp_size);
      d_ddd_size[0] = 1;
      d_ddd_size[1] = tmp_size[1];
      ixstart = tmp_size[0] * tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        d_ddd_data[b_n] = tmp_data[b_n];
      }

      /* second derivative to find the */
      roots(d_ddd_data, d_ddd_size, b_tmp_data, r_size);
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n].re = b_tmp_data[b_n].re;
        r_data[b_n].im = -b_tmp_data[b_n].im;
      }

      /* get all roots of first derivative */
      /* only real roots */
      d_r_size[0] = 1;
      d_r_size[1] = r_size[0];
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_r_data[b_n] = (r_data[b_n].im == 0.0);
      }

      eml_li_find(b_r_data, d_r_size, c_tmp_data, tmp_size);
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
      }

      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n] = b_tmp_data[b_n];
      }

      /* only roots in interval */
      e_r_size[0] = 1;
      e_r_size[1] = tmp_size[1];
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_r_data[b_n] = ((r_data[b_n].re > 0.0) && (r_data[b_n].re < T_data[i]));
      }

      eml_li_find(b_r_data, e_r_size, c_tmp_data, tmp_size);
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
      }

      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        r_data[b_n] = b_tmp_data[b_n];
      }

      /* also min at the border */
      c_tmp_size[0] = 1;
      c_tmp_size[1] = 2 + tmp_size[1];
      c_r_data[0].re = 0.0;
      c_r_data[0].im = 0.0;
      ixstart = tmp_size[1];
      for (b_n = 0; b_n < ixstart; b_n++) {
        c_r_data[b_n + 1] = r_data[b_n];
      }

      c_r_data[1 + tmp_size[1]].re = T_data[i];
      c_r_data[1 + tmp_size[1]].im = 0.0;
      polyval(s_d_data, s_d_size, c_r_data, c_tmp_size, d_dd_ts_data, tmp_size);
      ixstart = 1;
      b_n = tmp_size[1];
      d_r = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg5 = false;
        while ((!exitg5) && (ix + 1 <= b_n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            d_r = d_dd_ts_data[ix];
            exitg5 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= b_n) {
          if (b_eml_relop(d_dd_ts_data[ixstart], d_r)) {
            d_r = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the min value of the polynomial on the interval */
      ixstart = 1;
      b_n = tmp_size[1];
      mtmp = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg4 = false;
        while ((!exitg4) && (ix + 1 <= b_n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            mtmp = d_dd_ts_data[ix];
            exitg4 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= b_n) {
          if (eml_relop(d_dd_ts_data[ixstart], mtmp)) {
            mtmp = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the maximal value of the velocity in s direction */
      if (d_r.re <= 0.0) {
        /* something went terrible wrong. we are driving backwards */
        *flag = -1.0;
        exitg2 = 1;
      } else {
        /* now the bound for the curvature can be calculated */
        if (d_r.im == 0.0) {
          if (y.im == 0.0) {
            y_re = y.re / d_r.re;
            b_y = 0.0;
          } else if (y.re == 0.0) {
            y_re = 0.0;
            b_y = y.im / d_r.re;
          } else {
            y_re = y.re / d_r.re;
            b_y = y.im / d_r.re;
          }
        } else if (d_r.re == 0.0) {
          if (y.re == 0.0) {
            y_re = y.im / d_r.im;
            b_y = 0.0;
          } else if (y.im == 0.0) {
            y_re = 0.0;
            b_y = -(y.re / d_r.im);
          } else {
            y_re = y.im / d_r.im;
            b_y = -(y.re / d_r.im);
          }
        } else {
          brm = d_r.re;
          kappa_xy_max_re = fabs(d_r.im);
          if (brm > kappa_xy_max_re) {
            b_y = d_r.im / d_r.re;
            kappa_xy_max_re = d_r.re + b_y * d_r.im;
            y_re = (y.re + b_y * y.im) / kappa_xy_max_re;
            b_y = (y.im - b_y * y.re) / kappa_xy_max_re;
          } else if (kappa_xy_max_re == brm) {
            if (d_r.re > 0.0) {
              b_y = 0.5;
            } else {
              b_y = -0.5;
            }

            if (d_r.im > 0.0) {
              kappa_xy_max_re = 0.5;
            } else {
              kappa_xy_max_re = -0.5;
            }

            y_re = (y.re * b_y + y.im * kappa_xy_max_re) / brm;
            b_y = (y.im * b_y - y.re * kappa_xy_max_re) / brm;
          } else {
            b_y = d_r.re / d_r.im;
            kappa_xy_max_re = d_r.im + b_y * d_r.re;
            y_re = (b_y * y.re + y.im) / kappa_xy_max_re;
            b_y = (b_y * y.im - y.re) / kappa_xy_max_re;
          }
        }

        kappa_xy_max_re = fabs(kappa) + y_re;
        if (kappa_xy_max_re > kappaMax) {
          /* the curvature could be (according to our approximation) to big */
          notD_data[i] = 1.0;
        }

        /* MAXIMALE ORTH: BESCHLEUNIGUNG */
        roots(d_dd_data, d_dd_size, b_tmp_data, r_size);
        ixstart = r_size[0];
        for (b_n = 0; b_n < ixstart; b_n++) {
          r_data[b_n].re = b_tmp_data[b_n].re;
          r_data[b_n].im = -b_tmp_data[b_n].im;
        }

        /* get all roots of SECOND derivative */
        /* only real roots */
        f_r_size[0] = 1;
        f_r_size[1] = r_size[0];
        ixstart = r_size[0];
        for (b_n = 0; b_n < ixstart; b_n++) {
          d_r_data[b_n] = (r_data[b_n].im == 0.0);
        }

        eml_li_find(d_r_data, f_r_size, c_tmp_data, tmp_size);
        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
        }

        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          r_data[b_n] = b_tmp_data[b_n];
        }

        /* only roots in interval */
        g_r_size[0] = 1;
        g_r_size[1] = tmp_size[1];
        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          d_r_data[b_n] = ((r_data[b_n].re > 0.0) && (r_data[b_n].re < T_data[i]));
        }

        eml_li_find(d_r_data, g_r_size, c_tmp_data, tmp_size);
        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          b_tmp_data[b_n] = r_data[c_tmp_data[tmp_size[0] * b_n] - 1];
        }

        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          r_data[b_n] = b_tmp_data[b_n];
        }

        /* also max at the border */
        d_tmp_size[0] = 1;
        d_tmp_size[1] = 2 + tmp_size[1];
        d_tmp_data[0].re = 0.0;
        d_tmp_data[0].im = 0.0;
        ixstart = tmp_size[1];
        for (b_n = 0; b_n < ixstart; b_n++) {
          d_tmp_data[b_n + 1] = r_data[b_n];
        }

        d_tmp_data[1 + tmp_size[1]].re = T_data[i];
        d_tmp_data[1 + tmp_size[1]].im = 0.0;
        polyval(abstandMinusSafety_poly, d_d_size, d_tmp_data, d_tmp_size,
                d_dd_ts_data, tmp_size);

        /* HERE ONE NEEDS THE SECOND DERIVATIVE */
        ixstart = 1;
        b_n = tmp_size[1];
        y = d_dd_ts_data[0];
        if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
          ix = 1;
          exitg3 = false;
          while ((!exitg3) && (ix + 1 <= b_n)) {
            ixstart = ix + 1;
            if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im)))
            {
              y = d_dd_ts_data[ix];
              exitg3 = true;
            } else {
              ix++;
            }
          }
        }

        if (ixstart < tmp_size[1]) {
          while (ixstart + 1 <= b_n) {
            if (eml_relop(d_dd_ts_data[ixstart], y)) {
              y = d_dd_ts_data[ixstart];
            }

            ixstart++;
          }
        }

        /* get the maximal value of the polynomial on the interval */
        d_r.re = (mtmp.re * mtmp.re - mtmp.im * mtmp.im) + (y.re * y.re - y.im *
          y.im);
        d_r.im = (mtmp.re * mtmp.im + mtmp.im * mtmp.re) + (y.re * y.im + y.im *
          y.re);
        b_sqrt(&d_r);
        if ((d_r.re * d_r.re - d_r.im * d_r.im) * kappa_xy_max_re - (d_r.re *
             d_r.im + d_r.im * d_r.re) * b_y > aOrthMax) {
          /* the orthogonal acceleration is too much for the grip */
          notD_data[i] += 2.0;

          /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
        }

        i++;
        guard1 = false;
      }
    } else {
      /* % check for collision */
      i = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    exitg1 = false;
    while ((!exitg1) && (i <= n)) {
      /* s(t) should be an increasing function. This is guranteed because we */
      /* calculate the min of the 1st derivative above and throw an error if */
      /* this is non positive */
      dv6[0] = 0.0;
      dv6[1] = 0.0;
      dv6[2] = 0.0;
      dv6[3] = dataVeh[1];
      dv6[4] = dataVeh[0] - safetyS;
      for (b_n = 0; b_n < 5; b_n++) {
        abstandMinusSafety_poly[b_n] = dv6[b_n] - PS_data[b_n + 5 * i];
      }

      b_roots(abstandMinusSafety_poly, c_r_data, r_size);

      /* only real roots */
      h_r_size[0] = r_size[0];
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        e_r_data[b_n] = (c_r_data[b_n].im == 0.0);
      }

      b_eml_li_find(e_r_data, h_r_size, e_tmp_data, r_size);
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        f_r_data[b_n] = c_r_data[e_tmp_data[b_n] - 1];
      }

      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        c_r_data[b_n] = f_r_data[b_n];
      }

      /* only roots in interval */
      i_r_size[0] = r_size[0];
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        e_r_data[b_n] = ((c_r_data[b_n].re > 0.0) && (c_r_data[b_n].re <
          T_data[i]));
      }

      b_eml_li_find(e_r_data, i_r_size, e_tmp_data, r_size);
      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        f_r_data[b_n] = c_r_data[e_tmp_data[b_n] - 1];
      }

      ixstart = r_size[0];
      for (b_n = 0; b_n < ixstart; b_n++) {
        c_r_data[b_n] = f_r_data[b_n];
      }

      b_guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      if (r_size[0] == 0) {
        kappa_xy_max_re = T_data[i] / 2.0;
        b_y = abstandMinusSafety_poly[0];
        for (ixstart = 0; ixstart < 4; ixstart++) {
          b_y = kappa_xy_max_re * b_y + abstandMinusSafety_poly[ixstart + 1];
        }

        if (b_y > 0.0) {
          /* all good no collsion possible */
          b_guard1 = true;
        } else {
          guard4 = true;
        }
      } else {
        guard4 = true;
      }

      if (guard4) {
        if (r_size[0] == 0) {
          kappa_xy_max_re = T_data[i] / 2.0;
          b_y = abstandMinusSafety_poly[0];
          for (ixstart = 0; ixstart < 4; ixstart++) {
            b_y = kappa_xy_max_re * b_y + abstandMinusSafety_poly[ixstart + 1];
          }

          if (b_y <= 0.0) {
            /* fix collsion */
            coll_data[i] = 1.0;
            b_guard1 = true;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }
      }

      if (guard3) {
        if (r_size[0] == 1) {
          if (c_r_data[0].im == 0.0) {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = 0.0;
          } else if (c_r_data[0].re == 0.0) {
            d_r.re = 0.0;
            d_r.im = c_r_data[0].im / 2.0;
          } else {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = c_r_data[0].im / 2.0;
          }

          y.re = abstandMinusSafety_poly[0];
          y.im = 0.0;
          for (ixstart = 0; ixstart < 4; ixstart++) {
            kappa_xy_max_re = d_r.re * y.im + d_r.im * y.re;
            y.re = (d_r.re * y.re - d_r.im * y.im) +
              abstandMinusSafety_poly[ixstart + 1];
            y.im = kappa_xy_max_re;
          }

          if (y.re > 0.0) {
            /* one unique possible collision possible */
            d_r.re = PD_data[6 * i];
            d_r.im = 0.0;
            for (ixstart = 0; ixstart < 5; ixstart++) {
              kappa_xy_max_re = c_r_data[0].re * d_r.im + c_r_data[0].im *
                d_r.re;
              d_r.re = (c_r_data[0].re * d_r.re - c_r_data[0].im * d_r.im) +
                PD_data[(ixstart + 6 * i) + 1];
              d_r.im = kappa_xy_max_re;
            }

            if ((dataVeh[2] == -1.0) && (d_r.re < safetyD)) {
              coll_data[i] = 1.0;
            }

            if ((dataVeh[2] == 1.0) && (d_r.re > safetyD)) {
              coll_data[i] = 1.0;
            }

            b_guard1 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        if (r_size[0] == 1) {
          if (c_r_data[0].im == 0.0) {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = 0.0;
          } else if (c_r_data[0].re == 0.0) {
            d_r.re = 0.0;
            d_r.im = c_r_data[0].im / 2.0;
          } else {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = c_r_data[0].im / 2.0;
          }

          y.re = abstandMinusSafety_poly[0];
          y.im = 0.0;
          for (ixstart = 0; ixstart < 4; ixstart++) {
            kappa_xy_max_re = d_r.re * y.im + d_r.im * y.re;
            y.re = (d_r.re * y.re - d_r.im * y.im) +
              abstandMinusSafety_poly[ixstart + 1];
            y.im = kappa_xy_max_re;
          }

          if (y.re > 0.0) {
            /* there is a collision directly at the start */
            coll_data[i] = -1.0;
            b_guard1 = true;
          } else {
            guard1 = true;
            exitg1 = true;
          }
        } else {
          guard1 = true;
          exitg1 = true;
        }
      }

      if (b_guard1) {
        i++;
      }
    }
  }

  if (guard1) {
    *flag = -3.0;
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
 *                double *flag
 * Return Type  : void
 */
void otg_smart_objFun(const double T[3], const double S[3], const double D[4],
                      double kj, double kT, double ks, double kd, const double
                      dataVeh[3], double safetyS, double safetyD, double kappa,
                      double kappaMax, double aOrthMax, double Ctot[3], double
                      notD[3], double coll[3], double *flag)
{
  int ixstart;
  double PS[15];
  double PD[18];
  int i;
  double abstandMinusSafety_poly[5];
  double dv2[6];
  boolean_T guard1 = false;
  int32_T exitg2;
  int d_d_size[2];
  int d_dd_size[2];
  double d_dd_data[4];
  int tmp_size[2];
  double tmp_data[4];
  int d_ddd_size[2];
  int n;
  double d_ddd_data[3];
  int r_size[1];
  creal_T b_tmp_data[3];
  creal_T r_data[3];
  boolean_T b_r_data[2];
  int b_r_size[2];
  int c_tmp_data[3];
  int c_r_size[2];
  creal_T c_r_data[4];
  int b_tmp_size[2];
  creal_T d_d_ts_data[5];
  creal_T d_dd_ts_data[4];
  creal_T y;
  int ix;
  boolean_T exitg6;
  int s_d_size[2];
  double s_d_data[4];
  int d_r_size[2];
  int e_r_size[2];
  int c_tmp_size[2];
  creal_T d_r;
  boolean_T exitg5;
  creal_T mtmp;
  boolean_T exitg4;
  double y_re;
  double b_y;
  double brm;
  double kappa_xy_max_re;
  boolean_T d_r_data[3];
  int f_r_size[2];
  int g_r_size[2];
  creal_T d_tmp_data[5];
  int d_tmp_size[2];
  boolean_T exitg3;
  boolean_T exitg1;
  double dv3[5];
  boolean_T e_r_data[4];
  int h_r_size[1];
  int e_tmp_data[4];
  int i_r_size[1];
  boolean_T b_guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;

  /* % init */
  *flag = 1.0;
  for (ixstart = 0; ixstart < 3; ixstart++) {
    notD[ixstart] = 0.0;
    coll[ixstart] = 0.0;
  }

  /* % Get the values of ps, pd for the Ts */
  for (i = 0; i < 3; i++) {
    otg_ps(S, T[i], abstandMinusSafety_poly);
    for (ixstart = 0; ixstart < 5; ixstart++) {
      PS[ixstart + 5 * i] = abstandMinusSafety_poly[ixstart];
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
  i = 0;
  guard1 = false;
  do {
    exitg2 = 0;
    if (i < 3) {
      /* LENKEINSCHLAG */
      polyder(*(double (*)[6])&PD[6 * i], abstandMinusSafety_poly, d_d_size);

      /* first derivative w.r.t. t of d(t) */
      b_polyder(abstandMinusSafety_poly, d_d_size, d_dd_data, d_dd_size);

      /* second derivative w.r.t. t of d(t) */
      b_polyder(d_dd_data, d_dd_size, tmp_data, tmp_size);
      d_ddd_size[0] = 1;
      d_ddd_size[1] = tmp_size[1];
      n = tmp_size[0] * tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_ddd_data[ixstart] = tmp_data[ixstart];
      }

      /* third derivative to find the max of the second derivative */
      roots(d_ddd_data, d_ddd_size, b_tmp_data, r_size);
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart].re = b_tmp_data[ixstart].re;
        r_data[ixstart].im = -b_tmp_data[ixstart].im;
      }

      /* get all roots of THIRD derivative */
      /* only real roots */
      b_r_size[0] = 1;
      b_r_size[1] = r_size[0];
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_r_data[ixstart] = (r_data[ixstart].im == 0.0);
      }

      eml_li_find(b_r_data, b_r_size, c_tmp_data, tmp_size);
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
      }

      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart] = b_tmp_data[ixstart];
      }

      /* only roots in interval */
      c_r_size[0] = 1;
      c_r_size[1] = tmp_size[1];
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re <
          T[i]));
      }

      eml_li_find(b_r_data, c_r_size, c_tmp_data, tmp_size);
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
      }

      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart] = b_tmp_data[ixstart];
      }

      /* also max at the border */
      b_tmp_size[0] = 1;
      b_tmp_size[1] = 2 + tmp_size[1];
      c_r_data[0].re = 0.0;
      c_r_data[0].im = 0.0;
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        c_r_data[ixstart + 1] = r_data[ixstart];
      }

      c_r_data[1 + tmp_size[1]].re = T[i];
      c_r_data[1 + tmp_size[1]].im = 0.0;
      polyval(d_dd_data, d_dd_size, c_r_data, b_tmp_size, d_d_ts_data, tmp_size);
      n = tmp_size[0] * tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_dd_ts_data[ixstart] = d_d_ts_data[ixstart];
      }

      /* HERE ONE NEEDS THE SECOND DERIVATIVE */
      ixstart = 1;
      n = tmp_size[1];
      y = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg6 = false;
        while ((!exitg6) && (ix + 1 <= n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            y = d_dd_ts_data[ix];
            exitg6 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= n) {
          if (eml_relop(d_dd_ts_data[ixstart], y)) {
            y = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the maximal value of the polynomial on the interval */
      c_polyder(*(double (*)[5])&PS[5 * i], s_d_data, s_d_size);

      /* first derivative w.r.t. t of s(t) = v(t) */
      b_polyder(s_d_data, s_d_size, tmp_data, tmp_size);
      d_ddd_size[0] = 1;
      d_ddd_size[1] = tmp_size[1];
      n = tmp_size[0] * tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_ddd_data[ixstart] = tmp_data[ixstart];
      }

      /* second derivative to find the */
      roots(d_ddd_data, d_ddd_size, b_tmp_data, r_size);
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart].re = b_tmp_data[ixstart].re;
        r_data[ixstart].im = -b_tmp_data[ixstart].im;
      }

      /* get all roots of first derivative */
      /* only real roots */
      d_r_size[0] = 1;
      d_r_size[1] = r_size[0];
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_r_data[ixstart] = (r_data[ixstart].im == 0.0);
      }

      eml_li_find(b_r_data, d_r_size, c_tmp_data, tmp_size);
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
      }

      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart] = b_tmp_data[ixstart];
      }

      /* only roots in interval */
      e_r_size[0] = 1;
      e_r_size[1] = tmp_size[1];
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re <
          T[i]));
      }

      eml_li_find(b_r_data, e_r_size, c_tmp_data, tmp_size);
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
      }

      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        r_data[ixstart] = b_tmp_data[ixstart];
      }

      /* also min at the border */
      c_tmp_size[0] = 1;
      c_tmp_size[1] = 2 + tmp_size[1];
      c_r_data[0].re = 0.0;
      c_r_data[0].im = 0.0;
      n = tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        c_r_data[ixstart + 1] = r_data[ixstart];
      }

      c_r_data[1 + tmp_size[1]].re = T[i];
      c_r_data[1 + tmp_size[1]].im = 0.0;
      polyval(s_d_data, s_d_size, c_r_data, c_tmp_size, d_d_ts_data, tmp_size);
      n = tmp_size[0] * tmp_size[1];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_dd_ts_data[ixstart] = d_d_ts_data[ixstart];
      }

      ixstart = 1;
      n = tmp_size[1];
      d_r = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg5 = false;
        while ((!exitg5) && (ix + 1 <= n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            d_r = d_dd_ts_data[ix];
            exitg5 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= n) {
          if (b_eml_relop(d_dd_ts_data[ixstart], d_r)) {
            d_r = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the min value of the polynomial on the interval */
      ixstart = 1;
      n = tmp_size[1];
      mtmp = d_dd_ts_data[0];
      if (rtIsNaN(d_dd_ts_data[0].re) || rtIsNaN(d_dd_ts_data[0].im)) {
        ix = 1;
        exitg4 = false;
        while ((!exitg4) && (ix + 1 <= n)) {
          ixstart = ix + 1;
          if (!(rtIsNaN(d_dd_ts_data[ix].re) || rtIsNaN(d_dd_ts_data[ix].im))) {
            mtmp = d_dd_ts_data[ix];
            exitg4 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < tmp_size[1]) {
        while (ixstart + 1 <= n) {
          if (eml_relop(d_dd_ts_data[ixstart], mtmp)) {
            mtmp = d_dd_ts_data[ixstart];
          }

          ixstart++;
        }
      }

      /* get the maximal value of the velocity in s direction */
      if (d_r.re <= 0.0) {
        /* something went terrible wrong. we are driving backwards */
        *flag = -1.0;
        exitg2 = 1;
      } else {
        /* now the bound for the curvature can be calculated */
        if (d_r.im == 0.0) {
          if (y.im == 0.0) {
            y_re = y.re / d_r.re;
            b_y = 0.0;
          } else if (y.re == 0.0) {
            y_re = 0.0;
            b_y = y.im / d_r.re;
          } else {
            y_re = y.re / d_r.re;
            b_y = y.im / d_r.re;
          }
        } else if (d_r.re == 0.0) {
          if (y.re == 0.0) {
            y_re = y.im / d_r.im;
            b_y = 0.0;
          } else if (y.im == 0.0) {
            y_re = 0.0;
            b_y = -(y.re / d_r.im);
          } else {
            y_re = y.im / d_r.im;
            b_y = -(y.re / d_r.im);
          }
        } else {
          brm = d_r.re;
          kappa_xy_max_re = fabs(d_r.im);
          if (brm > kappa_xy_max_re) {
            b_y = d_r.im / d_r.re;
            kappa_xy_max_re = d_r.re + b_y * d_r.im;
            y_re = (y.re + b_y * y.im) / kappa_xy_max_re;
            b_y = (y.im - b_y * y.re) / kappa_xy_max_re;
          } else if (kappa_xy_max_re == brm) {
            if (d_r.re > 0.0) {
              b_y = 0.5;
            } else {
              b_y = -0.5;
            }

            if (d_r.im > 0.0) {
              kappa_xy_max_re = 0.5;
            } else {
              kappa_xy_max_re = -0.5;
            }

            y_re = (y.re * b_y + y.im * kappa_xy_max_re) / brm;
            b_y = (y.im * b_y - y.re * kappa_xy_max_re) / brm;
          } else {
            b_y = d_r.re / d_r.im;
            kappa_xy_max_re = d_r.im + b_y * d_r.re;
            y_re = (b_y * y.re + y.im) / kappa_xy_max_re;
            b_y = (b_y * y.im - y.re) / kappa_xy_max_re;
          }
        }

        kappa_xy_max_re = fabs(kappa) + y_re;
        if (kappa_xy_max_re > kappaMax) {
          /* the curvature could be (according to our approximation) to big */
          notD[i] = 1.0;
        }

        /* MAXIMALE ORTH: BESCHLEUNIGUNG */
        roots(d_dd_data, d_dd_size, b_tmp_data, r_size);
        n = r_size[0];
        for (ixstart = 0; ixstart < n; ixstart++) {
          r_data[ixstart].re = b_tmp_data[ixstart].re;
          r_data[ixstart].im = -b_tmp_data[ixstart].im;
        }

        /* get all roots of SECOND derivative */
        /* only real roots */
        f_r_size[0] = 1;
        f_r_size[1] = r_size[0];
        n = r_size[0];
        for (ixstart = 0; ixstart < n; ixstart++) {
          d_r_data[ixstart] = (r_data[ixstart].im == 0.0);
        }

        eml_li_find(d_r_data, f_r_size, c_tmp_data, tmp_size);
        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
        }

        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          r_data[ixstart] = b_tmp_data[ixstart];
        }

        /* only roots in interval */
        g_r_size[0] = 1;
        g_r_size[1] = tmp_size[1];
        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          d_r_data[ixstart] = ((r_data[ixstart].re > 0.0) && (r_data[ixstart].re
            < T[i]));
        }

        eml_li_find(d_r_data, g_r_size, c_tmp_data, tmp_size);
        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          b_tmp_data[ixstart] = r_data[c_tmp_data[tmp_size[0] * ixstart] - 1];
        }

        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          r_data[ixstart] = b_tmp_data[ixstart];
        }

        /* also max at the border */
        d_tmp_size[0] = 1;
        d_tmp_size[1] = 2 + tmp_size[1];
        d_tmp_data[0].re = 0.0;
        d_tmp_data[0].im = 0.0;
        n = tmp_size[1];
        for (ixstart = 0; ixstart < n; ixstart++) {
          d_tmp_data[ixstart + 1] = r_data[ixstart];
        }

        d_tmp_data[1 + tmp_size[1]].re = T[i];
        d_tmp_data[1 + tmp_size[1]].im = 0.0;
        polyval(abstandMinusSafety_poly, d_d_size, d_tmp_data, d_tmp_size,
                d_d_ts_data, tmp_size);

        /* HERE ONE NEEDS THE SECOND DERIVATIVE */
        ixstart = 1;
        n = tmp_size[1];
        y = d_d_ts_data[0];
        if (rtIsNaN(d_d_ts_data[0].re) || rtIsNaN(d_d_ts_data[0].im)) {
          ix = 1;
          exitg3 = false;
          while ((!exitg3) && (ix + 1 <= n)) {
            ixstart = ix + 1;
            if (!(rtIsNaN(d_d_ts_data[ix].re) || rtIsNaN(d_d_ts_data[ix].im))) {
              y = d_d_ts_data[ix];
              exitg3 = true;
            } else {
              ix++;
            }
          }
        }

        if (ixstart < tmp_size[1]) {
          while (ixstart + 1 <= n) {
            if (eml_relop(d_d_ts_data[ixstart], y)) {
              y = d_d_ts_data[ixstart];
            }

            ixstart++;
          }
        }

        /* get the maximal value of the polynomial on the interval */
        d_r.re = (mtmp.re * mtmp.re - mtmp.im * mtmp.im) + (y.re * y.re - y.im *
          y.im);
        d_r.im = (mtmp.re * mtmp.im + mtmp.im * mtmp.re) + (y.re * y.im + y.im *
          y.re);
        b_sqrt(&d_r);
        if ((d_r.re * d_r.re - d_r.im * d_r.im) * kappa_xy_max_re - (d_r.re *
             d_r.im + d_r.im * d_r.re) * b_y > aOrthMax) {
          /* the orthogonal acceleration is too much for the grip */
          notD[i] += 2.0;

          /* must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.) */
        }

        i++;
        guard1 = false;
      }
    } else {
      /* % check for collision */
      i = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    exitg1 = false;
    while ((!exitg1) && (i < 3)) {
      /* s(t) should be an increasing function. This is guranteed because we */
      /* calculate the min of the 1st derivative above and throw an error if */
      /* this is non positive */
      dv3[0] = 0.0;
      dv3[1] = 0.0;
      dv3[2] = 0.0;
      dv3[3] = dataVeh[1];
      dv3[4] = dataVeh[0] - safetyS;
      for (ixstart = 0; ixstart < 5; ixstart++) {
        abstandMinusSafety_poly[ixstart] = dv3[ixstart] - PS[ixstart + 5 * i];
      }

      b_roots(abstandMinusSafety_poly, c_r_data, r_size);

      /* only real roots */
      h_r_size[0] = r_size[0];
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        e_r_data[ixstart] = (c_r_data[ixstart].im == 0.0);
      }

      b_eml_li_find(e_r_data, h_r_size, e_tmp_data, r_size);
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_dd_ts_data[ixstart] = c_r_data[e_tmp_data[ixstart] - 1];
      }

      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        c_r_data[ixstart] = d_dd_ts_data[ixstart];
      }

      /* only roots in interval */
      i_r_size[0] = r_size[0];
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        e_r_data[ixstart] = ((c_r_data[ixstart].re > 0.0) && (c_r_data[ixstart].
          re < T[i]));
      }

      b_eml_li_find(e_r_data, i_r_size, e_tmp_data, r_size);
      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        d_dd_ts_data[ixstart] = c_r_data[e_tmp_data[ixstart] - 1];
      }

      n = r_size[0];
      for (ixstart = 0; ixstart < n; ixstart++) {
        c_r_data[ixstart] = d_dd_ts_data[ixstart];
      }

      b_guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      if (r_size[0] == 0) {
        kappa_xy_max_re = T[i] / 2.0;
        b_y = abstandMinusSafety_poly[0];
        for (ixstart = 0; ixstart < 4; ixstart++) {
          b_y = kappa_xy_max_re * b_y + abstandMinusSafety_poly[ixstart + 1];
        }

        if (b_y > 0.0) {
          /* all good no collsion possible */
          b_guard1 = true;
        } else {
          guard4 = true;
        }
      } else {
        guard4 = true;
      }

      if (guard4) {
        if (r_size[0] == 0) {
          kappa_xy_max_re = T[i] / 2.0;
          b_y = abstandMinusSafety_poly[0];
          for (ixstart = 0; ixstart < 4; ixstart++) {
            b_y = kappa_xy_max_re * b_y + abstandMinusSafety_poly[ixstart + 1];
          }

          if (b_y <= 0.0) {
            /* fix collsion */
            coll[i] = 1.0;
            b_guard1 = true;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }
      }

      if (guard3) {
        if (r_size[0] == 1) {
          if (c_r_data[0].im == 0.0) {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = 0.0;
          } else if (c_r_data[0].re == 0.0) {
            d_r.re = 0.0;
            d_r.im = c_r_data[0].im / 2.0;
          } else {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = c_r_data[0].im / 2.0;
          }

          y.re = abstandMinusSafety_poly[0];
          y.im = 0.0;
          for (ixstart = 0; ixstart < 4; ixstart++) {
            kappa_xy_max_re = d_r.re * y.im + d_r.im * y.re;
            y.re = (d_r.re * y.re - d_r.im * y.im) +
              abstandMinusSafety_poly[ixstart + 1];
            y.im = kappa_xy_max_re;
          }

          if (y.re > 0.0) {
            /* one unique possible collision possible */
            d_r.re = PD[6 * i];
            d_r.im = 0.0;
            for (ixstart = 0; ixstart < 5; ixstart++) {
              kappa_xy_max_re = c_r_data[0].re * d_r.im + c_r_data[0].im *
                d_r.re;
              d_r.re = (c_r_data[0].re * d_r.re - c_r_data[0].im * d_r.im) + PD
                [(ixstart + 6 * i) + 1];
              d_r.im = kappa_xy_max_re;
            }

            if ((dataVeh[2] == -1.0) && (d_r.re < safetyD)) {
              coll[i] = 1.0;
            }

            if ((dataVeh[2] == 1.0) && (d_r.re > safetyD)) {
              coll[i] = 1.0;
            }

            b_guard1 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        if (r_size[0] == 1) {
          if (c_r_data[0].im == 0.0) {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = 0.0;
          } else if (c_r_data[0].re == 0.0) {
            d_r.re = 0.0;
            d_r.im = c_r_data[0].im / 2.0;
          } else {
            d_r.re = c_r_data[0].re / 2.0;
            d_r.im = c_r_data[0].im / 2.0;
          }

          y.re = abstandMinusSafety_poly[0];
          y.im = 0.0;
          for (ixstart = 0; ixstart < 4; ixstart++) {
            kappa_xy_max_re = d_r.re * y.im + d_r.im * y.re;
            y.re = (d_r.re * y.re - d_r.im * y.im) +
              abstandMinusSafety_poly[ixstart + 1];
            y.im = kappa_xy_max_re;
          }

          if (y.re > 0.0) {
            /* there is a collision directly at the start */
            coll[i] = -1.0;
            b_guard1 = true;
          } else {
            guard1 = true;
            exitg1 = true;
          }
        } else {
          guard1 = true;
          exitg1 = true;
        }
      }

      if (b_guard1) {
        i++;
      }
    }
  }

  if (guard1) {
    *flag = -3.0;
  }
}

/*
 * File trailer for otg_smart_objFun.c
 *
 * [EOF]
 */
