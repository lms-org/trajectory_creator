/*
 * File: eig.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 11:15:45
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "eig.h"
#include "sqrt.h"
#include "otg_xy_reallyDumb_rtwutil.h"

/* Function Declarations */
static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, double *cs,
  creal_T *sn);
static void eml_matlab_zhgeqz(const creal_T A_data[], const int A_size[2], int
  ilo, int ihi, int *info, creal_T alpha1_data[], int alpha1_size[1], creal_T
  beta1_data[], int beta1_size[1]);
static void eml_matlab_zlartg(const creal_T f, const creal_T g, double *cs,
  creal_T *sn, creal_T *r);

/* Function Definitions */

/*
 * Arguments    : const creal_T f
 *                const creal_T g
 *                double *cs
 *                creal_T *sn
 * Return Type  : void
 */
static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, double *cs,
  creal_T *sn)
{
  double scale;
  double f2s;
  double g2;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  boolean_T guard1 = false;
  double g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        scale = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          scale = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 7.4428285367870146E+137 * f.re;
          g2 = 7.4428285367870146E+137 * f.im;
          scale = rt_hypotd_snf(f2s, g2);
          fs_re = f2s / scale;
          fs_im = g2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      *cs = 1.0 / f2s;
      scale += g2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

/*
 * Arguments    : const creal_T A_data[]
 *                const int A_size[2]
 *                int ilo
 *                int ihi
 *                int *info
 *                creal_T alpha1_data[]
 *                int alpha1_size[1]
 *                creal_T beta1_data[]
 *                int beta1_size[1]
 * Return Type  : void
 */
static void eml_matlab_zhgeqz(const creal_T A_data[], const int A_size[2], int
  ilo, int ihi, int *info, creal_T alpha1_data[], int alpha1_size[1], creal_T
  beta1_data[], int beta1_size[1])
{
  int A_size_idx_0;
  int jp1;
  int jm1;
  creal_T b_A_data[16];
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double rho_re;
  double rho_im;
  double anorm;
  double sumsq;
  boolean_T firstNonZero;
  int j;
  int i;
  double reAij;
  double imAij;
  double temp1;
  double temp2;
  double b_atol;
  boolean_T failed;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T b_guard1 = false;
  creal_T t1;
  creal_T d;
  boolean_T exitg2;
  A_size_idx_0 = A_size[0];
  jp1 = A_size[0] * A_size[1];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    b_A_data[jm1] = A_data[jm1];
  }

  *info = 0;
  if ((A_size[0] == 1) && (A_size[1] == 1)) {
    ihi = 1;
  }

  alpha1_size[0] = A_size[0];
  jp1 = A_size[0];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    alpha1_data[jm1].re = 0.0;
    alpha1_data[jm1].im = 0.0;
  }

  beta1_size[0] = A_size[0];
  jp1 = A_size[0];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    beta1_data[jm1].re = 1.0;
    beta1_data[jm1].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  rho_re = 0.0;
  rho_im = 0.0;
  anorm = 0.0;
  if (ilo > ihi) {
  } else {
    anorm = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jm1 = j + 1;
      if (ihi < j + 1) {
        jm1 = ihi;
      }

      for (i = ilo; i <= jm1; i++) {
        reAij = A_data[(i + A_size[0] * (j - 1)) - 1].re;
        imAij = A_data[(i + A_size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          temp1 = fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            anorm = temp1;
            firstNonZero = false;
          } else if (anorm < temp1) {
            temp2 = anorm / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            anorm = temp1;
          } else {
            temp2 = temp1 / anorm;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          temp1 = fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            anorm = temp1;
            firstNonZero = false;
          } else if (anorm < temp1) {
            temp2 = anorm / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            anorm = temp1;
          } else {
            temp2 = temp1 / anorm;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm *= sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  temp1 = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    temp1 = anorm;
  }

  imAij = 1.0 / temp1;
  failed = true;
  for (j = ihi; j + 1 <= A_size[0]; j++) {
    alpha1_data[j] = A_data[j + A_size[0] * j];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ifrstm = ilo;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = true;
        } else if (fabs(b_A_data[ilast + A_size_idx_0 * ilastm1].re) + fabs
                   (b_A_data[ilast + A_size_idx_0 * ilastm1].im) <= b_atol) {
          b_A_data[ilast + A_size_idx_0 * ilastm1].re = 0.0;
          b_A_data[ilast + A_size_idx_0 * ilastm1].im = 0.0;
          goto60 = true;
        } else {
          j = ilastm1;
          exitg3 = false;
          while ((!exitg3) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = true;
            } else if (fabs(b_A_data[j + A_size_idx_0 * (j - 1)].re) + fabs
                       (b_A_data[j + A_size_idx_0 * (j - 1)].im) <= b_atol) {
              b_A_data[j + A_size_idx_0 * (j - 1)].re = 0.0;
              b_A_data[j + A_size_idx_0 * (j - 1)].im = 0.0;
              firstNonZero = true;
            } else {
              firstNonZero = false;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = true;
              exitg3 = true;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = true;
        } else {
          firstNonZero = false;
        }

        if (!firstNonZero) {
          jp1 = alpha1_size[0];
          for (jm1 = 0; jm1 < jp1; jm1++) {
            alpha1_data[jm1].re = rtNaN;
            alpha1_data[jm1].im = 0.0;
          }

          jp1 = beta1_size[0];
          for (jm1 = 0; jm1 < jp1; jm1++) {
            beta1_data[jm1].re = rtNaN;
            beta1_data[jm1].im = 0.0;
          }

          *info = 1;
          exitg1 = 1;
        } else {
          b_guard1 = false;
          if (goto60) {
            goto60 = false;
            alpha1_data[ilast] = b_A_data[ilast + A_size_idx_0 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              if (ifrstm > ilast + 1) {
                ifrstm = ilo;
              }

              b_guard1 = true;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              ifrstm = ifirst;
              if (iiter - iiter / 10 * 10 != 0) {
                reAij = -(b_A_data[ilast + A_size_idx_0 * ilast].re -
                          b_A_data[ilastm1 + A_size_idx_0 * ilastm1].re);
                temp1 = -(b_A_data[ilast + A_size_idx_0 * ilast].im -
                          b_A_data[ilastm1 + A_size_idx_0 * ilastm1].im);
                if (temp1 == 0.0) {
                  t1.re = reAij / 2.0;
                  t1.im = 0.0;
                } else if (reAij == 0.0) {
                  t1.re = 0.0;
                  t1.im = temp1 / 2.0;
                } else {
                  t1.re = reAij / 2.0;
                  t1.im = temp1 / 2.0;
                }

                d.re = (t1.re * t1.re - t1.im * t1.im) + (b_A_data[ilastm1 +
                  A_size_idx_0 * ilast].re * b_A_data[ilast + A_size_idx_0 *
                  ilastm1].re - b_A_data[ilastm1 + A_size_idx_0 * ilast].im *
                  b_A_data[ilast + A_size_idx_0 * ilastm1].im);
                d.im = (t1.re * t1.im + t1.im * t1.re) + (b_A_data[ilastm1 +
                  A_size_idx_0 * ilast].re * b_A_data[ilast + A_size_idx_0 *
                  ilastm1].im + b_A_data[ilastm1 + A_size_idx_0 * ilast].im *
                  b_A_data[ilast + A_size_idx_0 * ilastm1].re);
                b_sqrt(&d);
                rho_re = b_A_data[ilastm1 + A_size_idx_0 * ilastm1].re - (t1.re
                  - d.re);
                rho_im = b_A_data[ilastm1 + A_size_idx_0 * ilastm1].im - (t1.im
                  - d.im);
                anorm = b_A_data[ilastm1 + A_size_idx_0 * ilastm1].re - (t1.re +
                  d.re);
                reAij = b_A_data[ilastm1 + A_size_idx_0 * ilastm1].im - (t1.im +
                  d.im);
                if (rt_hypotd_snf(rho_re - b_A_data[ilast + A_size_idx_0 * ilast]
                                  .re, rho_im - b_A_data[ilast + A_size_idx_0 *
                                  ilast].im) <= rt_hypotd_snf(anorm -
                     b_A_data[ilast + A_size_idx_0 * ilast].re, reAij -
                     b_A_data[ilast + A_size_idx_0 * ilast].im)) {
                  anorm = rho_re;
                  reAij = rho_im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += b_A_data[ilast + A_size_idx_0 * ilastm1].re;
                eshift_im += b_A_data[ilast + A_size_idx_0 * ilastm1].im;
                anorm = eshift_re;
                reAij = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp.re = b_A_data[j + A_size_idx_0 * j].re - anorm;
                ctemp.im = b_A_data[j + A_size_idx_0 * j].im - reAij;
                temp1 = imAij * (fabs(ctemp.re) + fabs(ctemp.im));
                temp2 = imAij * (fabs(b_A_data[jp1 + A_size_idx_0 * j].re) +
                                 fabs(b_A_data[jp1 + A_size_idx_0 * j].im));
                sumsq = temp1;
                if (temp2 > temp1) {
                  sumsq = temp2;
                }

                if ((sumsq < 1.0) && (sumsq != 0.0)) {
                  temp1 /= sumsq;
                  temp2 /= sumsq;
                }

                if ((fabs(b_A_data[j + A_size_idx_0 * (j - 1)].re) + fabs
                     (b_A_data[j + A_size_idx_0 * (j - 1)].im)) * temp2 <= temp1
                    * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = b_A_data[(ifirst + A_size_idx_0 * (ifirst - 1)) - 1]
                    .re - anorm;
                  ctemp.im = b_A_data[(ifirst + A_size_idx_0 * (ifirst - 1)) - 1]
                    .im - reAij;
                }

                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              b_eml_matlab_zlartg(ctemp, b_A_data[istart + A_size_idx_0 *
                                  (istart - 1)], &anorm, &t1);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  eml_matlab_zlartg(b_A_data[(j + A_size_idx_0 * jm1) - 1],
                                    b_A_data[j + A_size_idx_0 * jm1], &anorm,
                                    &t1, &b_A_data[(j + A_size_idx_0 * jm1) - 1]);
                  b_A_data[j + A_size_idx_0 * jm1].re = 0.0;
                  b_A_data[j + A_size_idx_0 * jm1].im = 0.0;
                }

                for (jp1 = j - 1; jp1 + 1 <= ilastm; jp1++) {
                  d.re = anorm * b_A_data[(j + A_size_idx_0 * jp1) - 1].re +
                    (t1.re * b_A_data[j + A_size_idx_0 * jp1].re - t1.im *
                     b_A_data[j + A_size_idx_0 * jp1].im);
                  d.im = anorm * b_A_data[(j + A_size_idx_0 * jp1) - 1].im +
                    (t1.re * b_A_data[j + A_size_idx_0 * jp1].im + t1.im *
                     b_A_data[j + A_size_idx_0 * jp1].re);
                  reAij = b_A_data[(j + A_size_idx_0 * jp1) - 1].im;
                  temp1 = b_A_data[(j + A_size_idx_0 * jp1) - 1].re;
                  b_A_data[j + A_size_idx_0 * jp1].re = anorm * b_A_data[j +
                    A_size_idx_0 * jp1].re - (t1.re * b_A_data[(j + A_size_idx_0
                    * jp1) - 1].re + t1.im * b_A_data[(j + A_size_idx_0 * jp1) -
                    1].im);
                  b_A_data[j + A_size_idx_0 * jp1].im = anorm * b_A_data[j +
                    A_size_idx_0 * jp1].im - (t1.re * reAij - t1.im * temp1);
                  b_A_data[(j + A_size_idx_0 * jp1) - 1] = d;
                }

                t1.re = -t1.re;
                t1.im = -t1.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = ifrstm - 1; i + 1 <= jp1 + 2; i++) {
                  d.re = anorm * b_A_data[i + A_size_idx_0 * j].re + (t1.re *
                    b_A_data[i + A_size_idx_0 * (j - 1)].re - t1.im * b_A_data[i
                    + A_size_idx_0 * (j - 1)].im);
                  d.im = anorm * b_A_data[i + A_size_idx_0 * j].im + (t1.re *
                    b_A_data[i + A_size_idx_0 * (j - 1)].im + t1.im * b_A_data[i
                    + A_size_idx_0 * (j - 1)].re);
                  reAij = b_A_data[i + A_size_idx_0 * j].im;
                  temp1 = b_A_data[i + A_size_idx_0 * j].re;
                  b_A_data[i + A_size_idx_0 * (j - 1)].re = anorm * b_A_data[i +
                    A_size_idx_0 * (j - 1)].re - (t1.re * b_A_data[i +
                    A_size_idx_0 * j].re + t1.im * b_A_data[i + A_size_idx_0 * j]
                    .im);
                  b_A_data[i + A_size_idx_0 * (j - 1)].im = anorm * b_A_data[i +
                    A_size_idx_0 * (j - 1)].im - (t1.re * reAij - t1.im * temp1);
                  b_A_data[i + A_size_idx_0 * j] = d;
                }

                jm1 = j - 1;
                j++;
              }
            }

            b_guard1 = true;
          }

          if (b_guard1) {
            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1_data[jp1].re = rtNaN;
        alpha1_data[jp1].im = 0.0;
        beta1_data[jp1].re = rtNaN;
        beta1_data[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1_data[j] = b_A_data[j + A_size_idx_0 * j];
    }
  }
}

/*
 * Arguments    : const creal_T f
 *                const creal_T g
 *                double *cs
 *                creal_T *sn
 *                creal_T *r
 * Return Type  : void
 */
static void eml_matlab_zlartg(const creal_T f, const creal_T g, double *cs,
  creal_T *sn, creal_T *r)
{
  double scale;
  double f2s;
  double g2;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  double g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        f2s = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / f2s;
        sn->im = -gs_im / f2s;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          f2s = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / f2s;
          fs_im = f.im / f2s;
        } else {
          g2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          f2s = rt_hypotd_snf(g2, scale);
          fs_re = g2 / f2s;
          fs_im = scale / f2s;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2s = scale + g2;
      g2 = r->re / f2s;
      f2s = r->im / f2s;
      sn->re = g2 * gs_re - f2s * -gs_im;
      sn->im = g2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

/*
 * Arguments    : const creal_T A_data[]
 *                const int A_size[2]
 *                creal_T V_data[]
 *                int V_size[1]
 * Return Type  : void
 */
void eig(const creal_T A_data[], const int A_size[2], creal_T V_data[], int
         V_size[1])
{
  int b_A_size[2];
  int ii;
  int nzcount;
  creal_T b_A_data[16];
  double anrm;
  boolean_T exitg7;
  double absxk;
  creal_T beta1_data[4];
  boolean_T ilascl;
  double anrmto;
  double ctoc;
  boolean_T notdone;
  double stemp_re;
  double cto1;
  double mul;
  int ilo;
  int ihi;
  int32_T exitg2;
  int i;
  int j;
  boolean_T exitg5;
  int jj;
  boolean_T exitg6;
  boolean_T guard2 = false;
  int A_size_idx_1;
  creal_T c_A_data[16];
  creal_T atmp;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = false;
  int beta1_size[1];
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  ii = A_size[0] * A_size[1];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    b_A_data[nzcount].re = A_data[nzcount].re;
    b_A_data[nzcount].im = A_data[nzcount].im;
  }

  anrm = 0.0;
  nzcount = A_size[0] * A_size[1];
  ii = 0;
  exitg7 = false;
  while ((!exitg7) && (ii <= nzcount - 1)) {
    absxk = rt_hypotd_snf(b_A_data[ii].re, b_A_data[ii].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg7 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      ii++;
    }
  }

  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    V_size[0] = A_size[0];
    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      V_data[nzcount].re = rtNaN;
      V_data[nzcount].im = 0.0;
    }

    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      beta1_data[nzcount].re = rtNaN;
      beta1_data[nzcount].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_re = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((stemp_re > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = stemp_re;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = false;
        }

        ii = b_A_size[0] * b_A_size[1];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          b_A_data[nzcount].re *= mul;
          b_A_data[nzcount].im *= mul;
        }
      }
    }

    ilo = 0;
    ihi = A_size[0];
    if (A_size[0] <= 1) {
      ihi = 1;
    } else {
      do {
        exitg2 = 0;
        i = 0;
        j = 0;
        notdone = false;
        ii = ihi;
        exitg5 = false;
        while ((!exitg5) && (ii > 0)) {
          nzcount = 0;
          i = ii;
          j = ihi;
          jj = 1;
          exitg6 = false;
          while ((!exitg6) && (jj <= ihi)) {
            guard2 = false;
            if ((b_A_data[(ii + b_A_size[0] * (jj - 1)) - 1].re != 0.0) ||
                (b_A_data[(ii + b_A_size[0] * (jj - 1)) - 1].im != 0.0) || (ii ==
                 jj)) {
              if (nzcount == 0) {
                j = jj;
                nzcount = 1;
                guard2 = true;
              } else {
                nzcount = 2;
                exitg6 = true;
              }
            } else {
              guard2 = true;
            }

            if (guard2) {
              jj++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg5 = true;
          } else {
            ii--;
          }
        }

        if (!notdone) {
          exitg2 = 2;
        } else {
          jj = b_A_size[0];
          A_size_idx_1 = b_A_size[1];
          ii = b_A_size[0] * b_A_size[1];
          for (nzcount = 0; nzcount < ii; nzcount++) {
            c_A_data[nzcount] = b_A_data[nzcount];
          }

          if (i != ihi) {
            for (ii = 0; ii + 1 <= b_A_size[0]; ii++) {
              atmp = c_A_data[(i + jj * ii) - 1];
              c_A_data[(i + jj * ii) - 1] = c_A_data[(ihi + jj * ii) - 1];
              c_A_data[(ihi + jj * ii) - 1] = atmp;
            }
          }

          if (j != ihi) {
            for (ii = 0; ii + 1 <= ihi; ii++) {
              atmp = c_A_data[ii + jj * (j - 1)];
              c_A_data[ii + jj * (j - 1)] = c_A_data[ii + jj * (ihi - 1)];
              c_A_data[ii + jj * (ihi - 1)] = atmp;
            }
          }

          b_A_size[1] = A_size_idx_1;
          ii = jj * A_size_idx_1;
          for (nzcount = 0; nzcount < ii; nzcount++) {
            b_A_data[nzcount] = c_A_data[nzcount];
          }

          ihi--;
          if (ihi == 1) {
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);

      if (exitg2 == 1) {
      } else {
        do {
          exitg1 = 0;
          i = 0;
          j = 0;
          notdone = false;
          jj = ilo + 1;
          exitg3 = false;
          while ((!exitg3) && (jj <= ihi)) {
            nzcount = 0;
            i = ihi;
            j = jj;
            ii = ilo + 1;
            exitg4 = false;
            while ((!exitg4) && (ii <= ihi)) {
              guard1 = false;
              if ((b_A_data[(ii + b_A_size[0] * (jj - 1)) - 1].re != 0.0) ||
                  (b_A_data[(ii + b_A_size[0] * (jj - 1)) - 1].im != 0.0) || (ii
                   == jj)) {
                if (nzcount == 0) {
                  i = ii;
                  nzcount = 1;
                  guard1 = true;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                ii++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg3 = true;
            } else {
              jj++;
            }
          }

          if (!notdone) {
            exitg1 = 1;
          } else {
            jj = b_A_size[0];
            A_size_idx_1 = b_A_size[1];
            ii = b_A_size[0] * b_A_size[1];
            for (nzcount = 0; nzcount < ii; nzcount++) {
              c_A_data[nzcount] = b_A_data[nzcount];
            }

            if (i != ilo + 1) {
              for (ii = ilo; ii + 1 <= b_A_size[0]; ii++) {
                atmp = c_A_data[(i + jj * ii) - 1];
                c_A_data[(i + jj * ii) - 1] = c_A_data[ilo + jj * ii];
                c_A_data[ilo + jj * ii] = atmp;
              }
            }

            if (j != ilo + 1) {
              for (ii = 0; ii + 1 <= ihi; ii++) {
                atmp = c_A_data[ii + jj * (j - 1)];
                c_A_data[ii + jj * (j - 1)] = c_A_data[ii + jj * ilo];
                c_A_data[ii + jj * ilo] = atmp;
              }
            }

            b_A_size[1] = A_size_idx_1;
            ii = jj * A_size_idx_1;
            for (nzcount = 0; nzcount < ii; nzcount++) {
              b_A_data[nzcount] = c_A_data[nzcount];
            }

            ilo++;
            if (ilo + 1 == ihi) {
              exitg1 = 1;
            }
          }
        } while (exitg1 == 0);
      }
    }

    if ((b_A_size[0] <= 1) || (ihi < ilo + 3)) {
    } else {
      for (ii = ilo; ii + 1 < ihi - 1; ii++) {
        for (nzcount = ihi - 1; nzcount + 1 > ii + 2; nzcount--) {
          eml_matlab_zlartg(b_A_data[(nzcount + b_A_size[0] * ii) - 1],
                            b_A_data[nzcount + b_A_size[0] * ii], &absxk, &atmp,
                            &b_A_data[(nzcount + b_A_size[0] * ii) - 1]);
          b_A_data[nzcount + b_A_size[0] * ii].re = 0.0;
          b_A_data[nzcount + b_A_size[0] * ii].im = 0.0;
          for (j = ii + 1; j + 1 <= ihi; j++) {
            stemp_re = absxk * b_A_data[(nzcount + b_A_size[0] * j) - 1].re +
              (atmp.re * b_A_data[nzcount + b_A_size[0] * j].re - atmp.im *
               b_A_data[nzcount + b_A_size[0] * j].im);
            ctoc = absxk * b_A_data[(nzcount + b_A_size[0] * j) - 1].im +
              (atmp.re * b_A_data[nzcount + b_A_size[0] * j].im + atmp.im *
               b_A_data[nzcount + b_A_size[0] * j].re);
            cto1 = b_A_data[(nzcount + b_A_size[0] * j) - 1].im;
            mul = b_A_data[(nzcount + b_A_size[0] * j) - 1].re;
            b_A_data[nzcount + b_A_size[0] * j].re = absxk * b_A_data[nzcount +
              b_A_size[0] * j].re - (atmp.re * b_A_data[(nzcount + b_A_size[0] *
              j) - 1].re + atmp.im * b_A_data[(nzcount + b_A_size[0] * j) - 1].
              im);
            b_A_data[nzcount + b_A_size[0] * j].im = absxk * b_A_data[nzcount +
              b_A_size[0] * j].im - (atmp.re * cto1 - atmp.im * mul);
            b_A_data[(nzcount + b_A_size[0] * j) - 1].re = stemp_re;
            b_A_data[(nzcount + b_A_size[0] * j) - 1].im = ctoc;
          }

          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = ilo; i + 1 <= ihi; i++) {
            stemp_re = absxk * b_A_data[i + b_A_size[0] * nzcount].re + (atmp.re
              * b_A_data[i + b_A_size[0] * (nzcount - 1)].re - atmp.im *
              b_A_data[i + b_A_size[0] * (nzcount - 1)].im);
            ctoc = absxk * b_A_data[i + b_A_size[0] * nzcount].im + (atmp.re *
              b_A_data[i + b_A_size[0] * (nzcount - 1)].im + atmp.im *
              b_A_data[i + b_A_size[0] * (nzcount - 1)].re);
            cto1 = b_A_data[i + b_A_size[0] * nzcount].im;
            mul = b_A_data[i + b_A_size[0] * nzcount].re;
            b_A_data[i + b_A_size[0] * (nzcount - 1)].re = absxk * b_A_data[i +
              b_A_size[0] * (nzcount - 1)].re - (atmp.re * b_A_data[i +
              b_A_size[0] * nzcount].re + atmp.im * b_A_data[i + b_A_size[0] *
              nzcount].im);
            b_A_data[i + b_A_size[0] * (nzcount - 1)].im = absxk * b_A_data[i +
              b_A_size[0] * (nzcount - 1)].im - (atmp.re * cto1 - atmp.im * mul);
            b_A_data[i + b_A_size[0] * nzcount].re = stemp_re;
            b_A_data[i + b_A_size[0] * nzcount].im = ctoc;
          }
        }
      }
    }

    eml_matlab_zhgeqz(b_A_data, b_A_size, ilo + 1, ihi, &ii, V_data, V_size,
                      beta1_data, beta1_size);
    if ((ii != 0) || (!ilascl)) {
    } else {
      notdone = true;
      while (notdone) {
        stemp_re = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((stemp_re > anrm) && (anrm != 0.0)) {
          mul = 2.0041683600089728E-292;
          anrmto = stemp_re;
        } else if (cto1 > anrmto) {
          mul = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          mul = anrm / anrmto;
          notdone = false;
        }

        ii = V_size[0];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          V_data[nzcount].re *= mul;
          V_data[nzcount].im *= mul;
        }
      }
    }
  }

  ii = V_size[0];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    cto1 = V_data[nzcount].re;
    mul = V_data[nzcount].im;
    if (beta1_data[nzcount].im == 0.0) {
      if (V_data[nzcount].im == 0.0) {
        V_data[nzcount].re /= beta1_data[nzcount].re;
        V_data[nzcount].im = 0.0;
      } else if (V_data[nzcount].re == 0.0) {
        V_data[nzcount].re = 0.0;
        V_data[nzcount].im = mul / beta1_data[nzcount].re;
      } else {
        V_data[nzcount].re /= beta1_data[nzcount].re;
        V_data[nzcount].im = mul / beta1_data[nzcount].re;
      }
    } else if (beta1_data[nzcount].re == 0.0) {
      if (V_data[nzcount].re == 0.0) {
        V_data[nzcount].re = V_data[nzcount].im / beta1_data[nzcount].im;
        V_data[nzcount].im = 0.0;
      } else if (V_data[nzcount].im == 0.0) {
        V_data[nzcount].re = 0.0;
        V_data[nzcount].im = -(cto1 / beta1_data[nzcount].im);
      } else {
        V_data[nzcount].re = V_data[nzcount].im / beta1_data[nzcount].im;
        V_data[nzcount].im = -(cto1 / beta1_data[nzcount].im);
      }
    } else {
      stemp_re = fabs(beta1_data[nzcount].re);
      absxk = fabs(beta1_data[nzcount].im);
      if (stemp_re > absxk) {
        absxk = beta1_data[nzcount].im / beta1_data[nzcount].re;
        ctoc = beta1_data[nzcount].re + absxk * beta1_data[nzcount].im;
        V_data[nzcount].re = (V_data[nzcount].re + absxk * V_data[nzcount].im) /
          ctoc;
        V_data[nzcount].im = (mul - absxk * cto1) / ctoc;
      } else if (absxk == stemp_re) {
        if (beta1_data[nzcount].re > 0.0) {
          absxk = 0.5;
        } else {
          absxk = -0.5;
        }

        if (beta1_data[nzcount].im > 0.0) {
          ctoc = 0.5;
        } else {
          ctoc = -0.5;
        }

        V_data[nzcount].re = (V_data[nzcount].re * absxk + V_data[nzcount].im *
                              ctoc) / stemp_re;
        V_data[nzcount].im = (mul * absxk - cto1 * ctoc) / stemp_re;
      } else {
        absxk = beta1_data[nzcount].re / beta1_data[nzcount].im;
        ctoc = beta1_data[nzcount].im + absxk * beta1_data[nzcount].re;
        V_data[nzcount].re = (absxk * V_data[nzcount].re + V_data[nzcount].im) /
          ctoc;
        V_data[nzcount].im = (absxk * mul - cto1) / ctoc;
      }
    }
  }
}

/*
 * File trailer for eig.c
 *
 * [EOF]
 */
