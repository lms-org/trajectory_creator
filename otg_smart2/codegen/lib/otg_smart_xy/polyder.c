/*
 * File: polyder.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 08-Oct-2015 14:13:47
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "polyder.h"

/* Function Definitions */

/*
 * Arguments    : const double u_data[]
 *                const int u_size[2]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
void b_polyder(const double u_data[], const int u_size[2], double a_data[], int
               a_size[2])
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
 * Arguments    : const double u[5]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
void c_polyder(const double u[5], double a_data[], int a_size[2])
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
 * Arguments    : const double u[6]
 *                double a_data[]
 *                int a_size[2]
 * Return Type  : void
 */
void polyder(const double u[6], double a_data[], int a_size[2])
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
 * File trailer for polyder.c
 *
 * [EOF]
 */
