/*
 * File: roots.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "roots.h"
#include "eig.h"

/* Function Definitions */

/*
 * Arguments    : const double c[5]
 *                creal_T r_data[]
 *                int r_size[1]
 * Return Type  : void
 */
void b_roots(const double c[5], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int companDim;
  double ctmp[5];
  boolean_T exitg1;
  int j;
  boolean_T exitg2;
  int a_size[2];
  creal_T a_data[16];
  int eiga_size[1];
  creal_T eiga_data[4];
  memset(&r_data[0], 0, sizeof(creal_T) << 2);
  k1 = 1;
  while ((k1 <= 5) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 5;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }

  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInf(fabs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }

      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 5 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 5 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      k1 = companDim * companDim;
      for (j = 0; j < k1; j++) {
        a_data[j].re = 0.0;
        a_data[j].im = 0.0;
      }

      for (k1 = 0; k1 + 1 < companDim; k1++) {
        a_data[companDim * k1].re = -ctmp[k1];
        a_data[companDim * k1].im = 0.0;
        a_data[(k1 + companDim * k1) + 1].re = 1.0;
        a_data[(k1 + companDim * k1) + 1].im = 0.0;
      }

      a_data[companDim * (companDim - 1)].re = -ctmp[companDim - 1];
      a_data[companDim * (companDim - 1)].im = 0.0;
      for (k1 = 1; k1 <= 5 - k2; k1++) {
        r_data[k1 - 1].re = 0.0;
        r_data[k1 - 1].im = 0.0;
      }

      eig(a_data, a_size, eiga_data, eiga_size);
      for (k1 = 1; k1 <= companDim; k1++) {
        r_data[(k1 - k2) + 4] = eiga_data[k1 - 1];
      }

      r_size[0] = (companDim - k2) + 5;
    }
  } else if (1 > 5 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 5 - k2;
  }
}

/*
 * Arguments    : const double c_data[]
 *                const int c_size[2]
 *                creal_T r_data[]
 *                int r_size[1]
 * Return Type  : void
 */
void roots(const double c_data[], const int c_size[2], creal_T r_data[], int
           r_size[1])
{
  int k2;
  int k1;
  int nTrailingZeros;
  int companDim;
  double ctmp_data[4];
  boolean_T exitg1;
  boolean_T exitg2;
  int a_size[2];
  creal_T a_data[9];
  int tmp_size[1];
  creal_T tmp_data[4];
  creal_T eiga_data[3];
  k2 = c_size[1];
  for (k1 = 0; k1 <= k2 - 2; k1++) {
    r_data[k1].re = 0.0;
    r_data[k1].im = 0.0;
  }

  k1 = 1;
  while ((k1 <= c_size[1]) && (!(c_data[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = c_size[1];
  while ((k2 >= k1) && (!(c_data[k2 - 1] != 0.0))) {
    k2--;
  }

  nTrailingZeros = c_size[1] - k2;
  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      k2 = 0;
      exitg2 = false;
      while ((!exitg2) && (k2 + 1 <= companDim)) {
        ctmp_data[k2] = c_data[k1 + k2] / c_data[k1 - 1];
        if (rtIsInf(fabs(ctmp_data[k2]))) {
          exitg2 = true;
        } else {
          k2++;
        }
      }

      if (k2 + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > nTrailingZeros) {
        r_size[0] = 0;
      } else {
        r_size[0] = nTrailingZeros;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      k2 = companDim * companDim;
      for (k1 = 0; k1 < k2; k1++) {
        a_data[k1].re = 0.0;
        a_data[k1].im = 0.0;
      }

      for (k2 = 0; k2 + 1 < companDim; k2++) {
        a_data[companDim * k2].re = -ctmp_data[k2];
        a_data[companDim * k2].im = 0.0;
        a_data[(k2 + companDim * k2) + 1].re = 1.0;
        a_data[(k2 + companDim * k2) + 1].im = 0.0;
      }

      a_data[companDim * (companDim - 1)].re = -ctmp_data[companDim - 1];
      a_data[companDim * (companDim - 1)].im = 0.0;
      for (k2 = 1; k2 <= nTrailingZeros; k2++) {
        r_data[k2 - 1].re = 0.0;
        r_data[k2 - 1].im = 0.0;
      }

      eig(a_data, a_size, tmp_data, tmp_size);
      k2 = tmp_size[0];
      for (k1 = 0; k1 < k2; k1++) {
        eiga_data[k1] = tmp_data[k1];
      }

      for (k2 = 0; k2 + 1 <= companDim; k2++) {
        r_data[k2 + nTrailingZeros] = eiga_data[k2];
      }

      k2 = nTrailingZeros + companDim;
      if (1 > k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = k2;
      }
    }
  } else if (1 > nTrailingZeros) {
    r_size[0] = 0;
  } else {
    r_size[0] = nTrailingZeros;
  }
}

/*
 * File trailer for roots.c
 *
 * [EOF]
 */
