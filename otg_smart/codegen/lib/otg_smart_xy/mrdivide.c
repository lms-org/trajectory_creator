/*
 * File: mrdivide.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 13:10:03
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "mrdivide.h"

/* Function Declarations */
static double eml_qrsolve(const double A_data[], const int A_size[1], double
  B_data[]);
static double eml_xnrm2(int n, const double x_data[]);
static double rt_hypotd_snf(double u0, double u1);

/* Function Definitions */

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
  int b_A_size;
  m = A_size[0] - 2;
  if (A_size[0] <= 1) {
    mn = A_size[0];
  } else {
    mn = 1;
  }

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
          xnorm = rt_hypotd_snf(b_A_data[0], xnorm);
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
            xnorm = rt_hypotd_snf(atmp, xnorm);
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
    if (A_size[0] >= 1) {
      b_A_size = A_size[0];
    } else {
      b_A_size = 1;
    }

    xnorm = (double)b_A_size * fabs(b_A_data[0]) * 2.2204460492503131E-16;
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
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : double A_data[]
 *                int A_size[2]
 *                const double B_data[]
 *                const int B_size[2]
 * Return Type  : void
 */
void mrdivide(double A_data[], int A_size[2], const double B_data[], const int
              B_size[2])
{
  double b_B_data[3];
  int b_B_size[1];
  int loop_ub;
  int i9;
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
    for (i9 = 0; i9 < loop_ub; i9++) {
      b_B_data[i9] = B_data[B_size[0] * i9];
    }

    loop_ub = A_size[1];
    for (i9 = 0; i9 < loop_ub; i9++) {
      b_A_data[i9] = A_data[A_size[0] * i9];
    }

    d2 = eml_qrsolve(b_B_data, b_B_size, b_A_data);
    A_size[0] = 1;
    A_size[1] = 1;
    A_data[0] = d2;
  }
}

/*
 * File trailer for mrdivide.c
 *
 * [EOF]
 */
