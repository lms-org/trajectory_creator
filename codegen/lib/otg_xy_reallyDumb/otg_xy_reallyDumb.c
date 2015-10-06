/*
 * File: otg_xy_reallyDumb.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 11:15:45
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "otg_xy_reallyDumb_emxutil.h"
#include "eig.h"
#include "polyval.h"
#include "otg_pspdT_reallyDumb.h"

/* Type Definitions */
#ifndef struct_emxArray_creal_T_4
#define struct_emxArray_creal_T_4

struct emxArray_creal_T_4
{
  creal_T data[4];
  int size[1];
};

#endif                                 /*struct_emxArray_creal_T_4*/

#ifndef typedef_emxArray_creal_T_4
#define typedef_emxArray_creal_T_4

typedef struct emxArray_creal_T_4 emxArray_creal_T_4;

#endif                                 /*typedef_emxArray_creal_T_4*/

/* Function Definitions */

/*
 * This is the OTG = Optimal trajectory generation package. This is based on
 * the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
 * Frenet Frame"(2010) by Werling et. al.
 *
 * --------------------------------------------------------------------------
 *
 * otg_pspdT_reallyDumb returns the points specified by x, y on the
 * trajectory in the car coordinate system
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
 *                Ii  = -1/1: -1: on left lane, +1: on the right lane
 *        safetyS = min. safety distance in s direction
 *        safetyD = min. safety distance in D direction with sign:
 *            if Ii == -1: d(t) must be bigger  than safetyD
 *            if Ii ==  1: d(t) must be smaller than safetyD
 *        dt = sampling interval for the collision detection
 *
 *        ds = distance in s direction between two points
 *        m = number of points wanted
 *
 *        kappa = curvature of the center line
 *        y0 = distance from the center line along the y axis of the cosy of
 *            the car
 *        phi = angle between the x-axis of the car and the tangent on the
 *            center line at the section between center line and y axis
 *
 *    OUTPUT:
 *        flag1 = 1/0/-1
 *            1:  everything ok
 *            0:  solution was not unique. solution with minimal T is choosen
 *            -1: no solution: the vectors are full of zeros
 *        flag2 = 1/0/-1/-2
 *            1:  everything ok
 *            0.5: m too big, smaller one set
 *            0:  solution was not unique. solution with minimal T is choosen
 *            -1: the roots were imaginary
 *            -2: no roots
 *            -3: ds was too big
 *        x = x in the car coordinate sytsem
 *        y = y in the car coordinate sytsem
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
 *                double ds
 *                double m
 *                double kappa
 *                double b_y0
 *                double phi
 *                double *flag1
 *                double *flag2
 *                emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
void otg_xy_reallyDumb(const double S[3], const double D[4], double kj, double
  kT, double ks, double kd, double dT, double Tmin, double Tmax, const
  emxArray_real_T *dataVeh, double safetyS, double safetyD, double dt, double ds,
  double m, double kappa, double b_y0, double phi, double *flag1, double *flag2,
  emxArray_real_T *x, emxArray_real_T *y)
{
  int itilerow;
  int ibmat;
  double T;
  int pd_size[2];
  double pd_data[36];
  int ps_size[2];
  double ps_data[25];
  double b_flag1;
  int nc;
  emxArray_real_T *tt;
  emxArray_real_T *ss;
  emxArray_real_T *dd;
  emxArray_real_T *XY;
  emxArray_real_T *xy_car;
  emxArray_real_T *r0;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *b_x;
  boolean_T guard1 = false;
  double ctmp[5];
  double ps_new[5];
  int k2;
  int companDim;
  boolean_T exitg1;
  boolean_T exitg2;
  int a_size[2];
  creal_T a_data[16];
  emxArray_creal_T_4 b_a_data;
  double R[4];
  double dv0[2];

  /* % init */
  itilerow = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = (int)m;
  emxEnsureCapacity((emxArray__common *)x, itilerow, (int)sizeof(double));
  ibmat = (int)m;
  for (itilerow = 0; itilerow < ibmat; itilerow++) {
    x->data[itilerow] = 0.0;
  }

  itilerow = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)m;
  emxEnsureCapacity((emxArray__common *)y, itilerow, (int)sizeof(double));
  ibmat = (int)m;
  for (itilerow = 0; itilerow < ibmat; itilerow++) {
    y->data[itilerow] = 0.0;
  }

  *flag2 = 1.0;

  /* % coefficients */
  otg_pspdT_reallyDumb(S, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS,
                       safetyD, dt, &b_flag1, ps_data, ps_size, pd_data, pd_size,
                       &T);
  *flag1 = b_flag1;
  if (b_flag1 == -1.0) {
  } else {
    nc = ps_size[0] * ps_size[1];
    b_flag1 = ps_data[0];
    for (ibmat = 0; ibmat <= nc - 2; ibmat++) {
      b_flag1 = T * b_flag1 + ps_data[ibmat + 1];
    }

    emxInit_real_T(&tt, 2);
    emxInit_real_T(&ss, 2);
    emxInit_real_T(&dd, 2);
    emxInit_real_T(&XY, 2);
    emxInit_real_T(&xy_car, 2);
    emxInit_real_T(&r0, 2);
    emxInit_real_T(&r1, 2);
    emxInit_real_T(&r2, 2);
    emxInit_real_T(&b_x, 2);
    guard1 = false;
    if ((m - 1.0) * ds > b_flag1) {
      *flag2 = 0.5;
      nc = ps_size[0] * ps_size[1];
      b_flag1 = ps_data[0];
      for (ibmat = 0; ibmat <= nc - 2; ibmat++) {
        b_flag1 = T * b_flag1 + ps_data[ibmat + 1];
      }

      b_flag1 /= ds;
      m = floor(b_flag1) - 1.0;
      if (floor(b_flag1) - 1.0 < 1.0) {
        *flag2 = -3.0;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      /* % calculate the times where s(t) = i*ds */
      itilerow = tt->size[0] * tt->size[1];
      tt->size[0] = 1;
      tt->size[1] = (int)m;
      emxEnsureCapacity((emxArray__common *)tt, itilerow, (int)sizeof(double));
      ibmat = (int)m;
      for (itilerow = 0; itilerow < ibmat; itilerow++) {
        tt->data[itilerow] = 0.0;
      }

      if (0 <= (int)(m + -1.0) - 1) {
        ctmp[0] = 0.0;
        ctmp[1] = 0.0;
        ctmp[2] = 0.0;
        ctmp[3] = 0.0;
        ctmp[4] = ds;
        for (itilerow = 0; itilerow < 5; itilerow++) {
          ps_new[itilerow] = ps_data[itilerow] - ctmp[itilerow];
        }

        nc = 1;
        while ((nc <= 5) && (!(ps_new[nc - 1] != 0.0))) {
          nc++;
        }

        k2 = 5;
        while ((k2 >= nc) && (!(ps_new[k2 - 1] != 0.0))) {
          k2--;
        }

        if (nc < k2) {
          companDim = k2 - nc;
          exitg1 = false;
          while ((!exitg1) && (companDim > 0)) {
            ibmat = 0;
            exitg2 = false;
            while ((!exitg2) && (ibmat + 1 <= companDim)) {
              ctmp[ibmat] = ps_new[nc + ibmat] / ps_new[nc - 1];
              if (rtIsInf(fabs(ctmp[ibmat]))) {
                exitg2 = true;
              } else {
                ibmat++;
              }
            }

            if (ibmat + 1 > companDim) {
              exitg1 = true;
            } else {
              nc++;
              companDim--;
            }
          }

          if (companDim < 1) {
            if (1 > 5 - k2) {
              nc = 0;
            } else {
              nc = 5 - k2;
            }
          } else {
            a_size[0] = companDim;
            a_size[1] = companDim;
            ibmat = companDim * companDim;
            for (itilerow = 0; itilerow < ibmat; itilerow++) {
              a_data[itilerow].re = 0.0;
              a_data[itilerow].im = 0.0;
            }

            for (ibmat = 0; ibmat + 1 < companDim; ibmat++) {
              a_data[companDim * ibmat].re = -ctmp[ibmat];
              a_data[companDim * ibmat].im = 0.0;
              a_data[(ibmat + companDim * ibmat) + 1].re = 1.0;
              a_data[(ibmat + companDim * ibmat) + 1].im = 0.0;
            }

            a_data[companDim * (companDim - 1)].re = -ctmp[companDim - 1];
            a_data[companDim * (companDim - 1)].im = 0.0;
            eig(a_data, a_size, b_a_data.data, b_a_data.size);
            nc = (companDim - k2) + 5;
          }
        } else if (1 > 5 - k2) {
          nc = 0;
        } else {
          nc = 5 - k2;
        }

        if (nc == 0) {
          *flag2 = -2.0;
        } else {
          *flag2 = -1.0;
        }
      } else {
        /* % calc s and d at the times */
        polyval(ps_data, ps_size, tt, ss);
        polyval(pd_data, pd_size, tt, dd);
        itilerow = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = ss->size[1];
        emxEnsureCapacity((emxArray__common *)y, itilerow, (int)sizeof(double));
        ibmat = ss->size[0] * ss->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          y->data[itilerow] = kappa * ss->data[itilerow];
        }

        itilerow = r0->size[0] * r0->size[1];
        r0->size[0] = 1;
        r0->size[1] = y->size[1];
        emxEnsureCapacity((emxArray__common *)r0, itilerow, (int)sizeof(double));
        ibmat = y->size[0] * y->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          r0->data[itilerow] = y->data[itilerow];
        }

        for (ibmat = 0; ibmat < y->size[1]; ibmat++) {
          r0->data[ibmat] = sin(r0->data[ibmat]);
        }

        itilerow = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = ss->size[1];
        emxEnsureCapacity((emxArray__common *)y, itilerow, (int)sizeof(double));
        ibmat = ss->size[0] * ss->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          y->data[itilerow] = kappa * ss->data[itilerow];
        }

        itilerow = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = y->size[1];
        emxEnsureCapacity((emxArray__common *)r1, itilerow, (int)sizeof(double));
        ibmat = y->size[0] * y->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          r1->data[itilerow] = y->data[itilerow];
        }

        for (ibmat = 0; ibmat < y->size[1]; ibmat++) {
          r1->data[ibmat] = cos(r1->data[ibmat]);
        }

        itilerow = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = ss->size[1];
        emxEnsureCapacity((emxArray__common *)y, itilerow, (int)sizeof(double));
        ibmat = ss->size[0] * ss->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          y->data[itilerow] = kappa * ss->data[itilerow];
        }

        itilerow = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = y->size[1];
        emxEnsureCapacity((emxArray__common *)x, itilerow, (int)sizeof(double));
        ibmat = y->size[0] * y->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          x->data[itilerow] = y->data[itilerow];
        }

        for (ibmat = 0; ibmat < y->size[1]; ibmat++) {
          x->data[ibmat] = sin(x->data[ibmat]);
        }

        itilerow = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)ss, itilerow, (int)sizeof(double));
        nc = ss->size[0];
        ibmat = ss->size[1];
        ibmat *= nc;
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          ss->data[itilerow] *= kappa;
        }

        itilerow = tt->size[0] * tt->size[1];
        tt->size[0] = 1;
        tt->size[1] = ss->size[1];
        emxEnsureCapacity((emxArray__common *)tt, itilerow, (int)sizeof(double));
        ibmat = ss->size[0] * ss->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          tt->data[itilerow] = ss->data[itilerow];
        }

        for (ibmat = 0; ibmat < ss->size[1]; ibmat++) {
          tt->data[ibmat] = cos(tt->data[ibmat]);
        }

        itilerow = xy_car->size[0] * xy_car->size[1];
        xy_car->size[0] = 2;
        xy_car->size[1] = dd->size[1];
        emxEnsureCapacity((emxArray__common *)xy_car, itilerow, (int)sizeof
                          (double));
        nc = dd->size[1];
        if (nc == 0) {
        } else {
          for (nc = 0; nc + 1 <= dd->size[1]; nc++) {
            ibmat = nc << 1;
            for (itilerow = 0; itilerow < 2; itilerow++) {
              xy_car->data[ibmat + itilerow] = dd->data[nc];
            }
          }
        }

        b_flag1 = 1.0 / kappa;
        itilerow = r2->size[0] * r2->size[1];
        r2->size[0] = 2;
        r2->size[1] = r0->size[1];
        emxEnsureCapacity((emxArray__common *)r2, itilerow, (int)sizeof(double));
        ibmat = r0->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          r2->data[r2->size[0] * itilerow] = r0->data[r0->size[0] * itilerow];
        }

        ibmat = r1->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          r2->data[1 + r2->size[0] * itilerow] = 1.0 - r1->data[r1->size[0] *
            itilerow];
        }

        itilerow = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 2;
        b_x->size[1] = x->size[1];
        emxEnsureCapacity((emxArray__common *)b_x, itilerow, (int)sizeof(double));
        ibmat = x->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          b_x->data[b_x->size[0] * itilerow] = -x->data[x->size[0] * itilerow];
        }

        ibmat = tt->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          b_x->data[1 + b_x->size[0] * itilerow] = tt->data[tt->size[0] *
            itilerow];
        }

        itilerow = XY->size[0] * XY->size[1];
        XY->size[0] = 2;
        XY->size[1] = r2->size[1];
        emxEnsureCapacity((emxArray__common *)XY, itilerow, (int)sizeof(double));
        ibmat = r2->size[1];
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          for (nc = 0; nc < 2; nc++) {
            XY->data[nc + XY->size[0] * itilerow] = b_flag1 * r2->data[nc +
              r2->size[0] * itilerow] + b_x->data[nc + b_x->size[0] * itilerow] *
              xy_car->data[nc + xy_car->size[0] * itilerow];
          }
        }

        itilerow = xy_car->size[0] * xy_car->size[1];
        xy_car->size[0] = 2;
        xy_car->size[1] = (int)m;
        emxEnsureCapacity((emxArray__common *)xy_car, itilerow, (int)sizeof
                          (double));
        ibmat = (int)m << 1;
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          xy_car->data[itilerow] = 0.0;
        }

        R[0] = cos(phi);
        R[2] = -sin(phi);
        R[1] = sin(phi);
        R[3] = cos(phi);
        for (ibmat = 0; ibmat < (int)m; ibmat++) {
          b_flag1 = cos(phi) * b_y0;
          dv0[0] = -sin(phi);
          dv0[1] = cos(phi);
          for (itilerow = 0; itilerow < 2; itilerow++) {
            T = 0.0;
            for (nc = 0; nc < 2; nc++) {
              T += R[itilerow + (nc << 1)] * XY->data[nc + XY->size[0] * ibmat];
            }

            xy_car->data[itilerow + xy_car->size[0] * ibmat] = T + dv0[itilerow]
              * b_flag1;
          }
        }

        ibmat = xy_car->size[1];
        itilerow = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = ibmat;
        emxEnsureCapacity((emxArray__common *)x, itilerow, (int)sizeof(double));
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          x->data[x->size[0] * itilerow] = xy_car->data[xy_car->size[0] *
            itilerow];
        }

        ibmat = xy_car->size[1];
        itilerow = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = ibmat;
        emxEnsureCapacity((emxArray__common *)y, itilerow, (int)sizeof(double));
        for (itilerow = 0; itilerow < ibmat; itilerow++) {
          y->data[y->size[0] * itilerow] = xy_car->data[1 + xy_car->size[0] *
            itilerow];
        }
      }
    }

    emxFree_real_T(&b_x);
    emxFree_real_T(&r2);
    emxFree_real_T(&r1);
    emxFree_real_T(&r0);
    emxFree_real_T(&xy_car);
    emxFree_real_T(&XY);
    emxFree_real_T(&dd);
    emxFree_real_T(&ss);
    emxFree_real_T(&tt);
  }
}

/*
 * File trailer for otg_xy_reallyDumb.c
 *
 * [EOF]
 */
