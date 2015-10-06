/*
 * File: otg_xy_reallyDumb.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 16:22:34
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_xy_reallyDumb.h"
#include "otg_xy_reallyDumb_emxutil.h"
#include "repmat.h"
#include "polyval.h"
#include "linspace.h"
#include "otg_pspdT_reallyDumb.h"

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
 *        v1 = velocity in s direction at the end of the trajectory
 *        d1 = distance from the center line at the end of the trajectory
 *
 *        kj = weight for the jerk functional
 *        kT = weight for the time
 *        ks = weight of the longitudinal weight function
 *        kd = weight of the lateral weight function
 *
 *        dT = time interval between two possible end Times T
 *        Tmin = minimal T > 0
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
 *        ds = distance in s direction between two points
 *        m = number of points wanted
 *
 *        kappa = curvature of the center line
 *        phi = angle between the x-axis of the car and the tangent on the
 *            center line at the section between center line and y axis
 *
 *    OUTPUT:
 *        flag = 1/0/-1
 *            1:  everything ok
 *            0:  solution was not unique. solution with minimal T is choosen
 *            -1: no solution: the vectors are full of zeros
 *        x = x in the car coordinate sytsem
 *        y = y in the car coordinate sytsem
 *
 *  See also
 * Arguments    : double v1
 *                double d1
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
 *                double m
 *                double kappa
 *                double b_y0
 *                double phi
 *                double vx0
 *                double ax0
 *                double w
 *                double *flag
 *                emxArray_real_T *x
 *                emxArray_real_T *y
 *                double *T
 * Return Type  : void
 */
void otg_xy_reallyDumb(double v1, double d1, double kj, double kT, double ks,
  double kd, double dT, double Tmin, double Tmax, const emxArray_real_T *dataVeh,
  double safetyS, double safetyD, double dt, double m, double kappa, double b_y0,
  double phi, double vx0, double ax0, double w, double *flag, emxArray_real_T *x,
  emxArray_real_T *y, double *T)
{
  int jtilecol;
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
  emxArray_real_T *ss;
  emxArray_real_T *dd;
  emxArray_real_T *XY;
  emxArray_real_T *xy_car;
  emxArray_real_T *r0;
  int k;
  emxArray_real_T *b_ss;
  emxArray_real_T *r1;
  emxArray_real_T *b_x;
  int i0;
  double R[4];
  double dv1[2];

  /* % init */
  jtilecol = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = (int)m;
  emxEnsureCapacity((emxArray__common *)x, jtilecol, (int)sizeof(double));
  ibtile = (int)m;
  for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
    x->data[jtilecol] = 0.0;
  }

  jtilecol = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)m;
  emxEnsureCapacity((emxArray__common *)y, jtilecol, (int)sizeof(double));
  ibtile = (int)m;
  for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
    y->data[jtilecol] = 0.0;
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
  otg_pspdT_reallyDumb(dv0, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS,
                       safetyD, dt, &cosphi_d, ps_data, ps_size, pd_data,
                       pd_size, &y0_d);
  *flag = cosphi_d;
  *T = y0_d;
  if (cosphi_d == -1.0) {
  } else {
    emxInit_real_T(&tt, 2);
    emxInit_real_T(&ss, 2);
    emxInit_real_T(&dd, 2);

    /* % calculate tt */
    b_linspace(y0_d, m, tt);

    /* % calc s and d at the times */
    polyval(ps_data, ps_size, tt, ss);
    polyval(pd_data, pd_size, tt, dd);
    emxInit_real_T(&XY, 2);
    emxInit_real_T(&xy_car, 2);
    if (kappa > 0.0) {
      jtilecol = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, jtilecol, (int)sizeof(double));
      ibtile = ss->size[0] * ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        y->data[jtilecol] = kappa * ss->data[jtilecol];
      }

      emxInit_real_T(&r0, 2);
      jtilecol = r0->size[0] * r0->size[1];
      r0->size[0] = 1;
      r0->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)r0, jtilecol, (int)sizeof(double));
      ibtile = y->size[0] * y->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        r0->data[jtilecol] = y->data[jtilecol];
      }

      for (k = 0; k < y->size[1]; k++) {
        r0->data[k] = sin(r0->data[k]);
      }

      jtilecol = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, jtilecol, (int)sizeof(double));
      ibtile = ss->size[0] * ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        y->data[jtilecol] = kappa * ss->data[jtilecol];
      }

      emxInit_real_T(&b_ss, 2);
      jtilecol = b_ss->size[0] * b_ss->size[1];
      b_ss->size[0] = 1;
      b_ss->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)b_ss, jtilecol, (int)sizeof(double));
      ibtile = y->size[0] * y->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        b_ss->data[jtilecol] = y->data[jtilecol];
      }

      for (k = 0; k < y->size[1]; k++) {
        b_ss->data[k] = cos(b_ss->data[k]);
      }

      jtilecol = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)y, jtilecol, (int)sizeof(double));
      ibtile = ss->size[0] * ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        y->data[jtilecol] = kappa * ss->data[jtilecol];
      }

      jtilecol = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)x, jtilecol, (int)sizeof(double));
      ibtile = y->size[0] * y->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        x->data[jtilecol] = y->data[jtilecol];
      }

      for (k = 0; k < y->size[1]; k++) {
        x->data[k] = sin(x->data[k]);
      }

      jtilecol = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)ss, jtilecol, (int)sizeof(double));
      ibtile = ss->size[0];
      jtilecol = ss->size[1];
      ibtile *= jtilecol;
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        ss->data[jtilecol] *= kappa;
      }

      jtilecol = tt->size[0] * tt->size[1];
      tt->size[0] = 1;
      tt->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)tt, jtilecol, (int)sizeof(double));
      ibtile = ss->size[0] * ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        tt->data[jtilecol] = ss->data[jtilecol];
      }

      for (k = 0; k < ss->size[1]; k++) {
        tt->data[k] = cos(tt->data[k]);
      }

      emxInit_real_T(&r1, 2);
      repmat(dd, xy_car);
      cosphi_d = 1.0 / kappa;
      jtilecol = r1->size[0] * r1->size[1];
      r1->size[0] = 2;
      r1->size[1] = r0->size[1];
      emxEnsureCapacity((emxArray__common *)r1, jtilecol, (int)sizeof(double));
      ibtile = r0->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        r1->data[r1->size[0] * jtilecol] = r0->data[r0->size[0] * jtilecol];
      }

      emxFree_real_T(&r0);
      ibtile = b_ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        r1->data[1 + r1->size[0] * jtilecol] = 1.0 - b_ss->data[b_ss->size[0] *
          jtilecol];
      }

      emxFree_real_T(&b_ss);
      emxInit_real_T(&b_x, 2);
      jtilecol = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 2;
      b_x->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)b_x, jtilecol, (int)sizeof(double));
      ibtile = x->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        b_x->data[b_x->size[0] * jtilecol] = -x->data[x->size[0] * jtilecol];
      }

      ibtile = tt->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        b_x->data[1 + b_x->size[0] * jtilecol] = tt->data[tt->size[0] * jtilecol];
      }

      jtilecol = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)XY, jtilecol, (int)sizeof(double));
      ibtile = r1->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * jtilecol] = cosphi_d * r1->data[i0 +
            r1->size[0] * jtilecol] + b_x->data[i0 + b_x->size[0] * jtilecol] *
            xy_car->data[i0 + xy_car->size[0] * jtilecol];
        }
      }

      emxFree_real_T(&b_x);
      emxFree_real_T(&r1);
    } else {
      jtilecol = xy_car->size[0] * xy_car->size[1];
      xy_car->size[0] = 2;
      xy_car->size[1] = (int)m;
      emxEnsureCapacity((emxArray__common *)xy_car, jtilecol, (int)sizeof(double));
      if ((int)m == 0) {
      } else {
        for (jtilecol = 1; jtilecol <= (int)m; jtilecol++) {
          ibtile = (jtilecol - 1) << 1;
          for (k = 0; k < 2; k++) {
            xy_car->data[ibtile + k] = ((double)k + 1.0) - 1.0;
          }
        }
      }

      emxInit_real_T(&r0, 2);
      emxInit_real_T(&b_ss, 2);
      repmat(dd, r0);
      jtilecol = b_ss->size[0] * b_ss->size[1];
      b_ss->size[0] = 2;
      b_ss->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)b_ss, jtilecol, (int)sizeof(double));
      ibtile = ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        b_ss->data[b_ss->size[0] * jtilecol] = ss->data[ss->size[0] * jtilecol];
      }

      ibtile = (int)m;
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        b_ss->data[1 + b_ss->size[0] * jtilecol] = 0.0;
      }

      jtilecol = XY->size[0] * XY->size[1];
      XY->size[0] = 2;
      XY->size[1] = b_ss->size[1];
      emxEnsureCapacity((emxArray__common *)XY, jtilecol, (int)sizeof(double));
      ibtile = b_ss->size[1];
      for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
        for (i0 = 0; i0 < 2; i0++) {
          XY->data[i0 + XY->size[0] * jtilecol] = b_ss->data[i0 + b_ss->size[0] *
            jtilecol] + xy_car->data[i0 + xy_car->size[0] * jtilecol] * r0->
            data[i0 + r0->size[0] * jtilecol];
        }
      }

      emxFree_real_T(&b_ss);
      emxFree_real_T(&r0);
    }

    emxFree_real_T(&dd);
    emxFree_real_T(&ss);
    emxFree_real_T(&tt);
    jtilecol = xy_car->size[0] * xy_car->size[1];
    xy_car->size[0] = 2;
    xy_car->size[1] = (int)m;
    emxEnsureCapacity((emxArray__common *)xy_car, jtilecol, (int)sizeof(double));
    ibtile = (int)m << 1;
    for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
      xy_car->data[jtilecol] = 0.0;
    }

    R[0] = cos(phi);
    R[2] = -sin(phi);
    R[1] = sin(phi);
    R[3] = cos(phi);
    for (k = 0; k < (int)m; k++) {
      dv1[0] = -sin(phi);
      dv1[1] = cos(phi);
      for (jtilecol = 0; jtilecol < 2; jtilecol++) {
        cosphi_d = 0.0;
        for (i0 = 0; i0 < 2; i0++) {
          cosphi_d += R[jtilecol + (i0 << 1)] * XY->data[i0 + XY->size[0] * k];
        }

        xy_car->data[jtilecol + xy_car->size[0] * k] = cosphi_d + dv1[jtilecol] *
          D[0];
      }
    }

    emxFree_real_T(&XY);
    ibtile = xy_car->size[1];
    jtilecol = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = ibtile;
    emxEnsureCapacity((emxArray__common *)x, jtilecol, (int)sizeof(double));
    for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
      x->data[x->size[0] * jtilecol] = xy_car->data[xy_car->size[0] * jtilecol];
    }

    ibtile = xy_car->size[1];
    jtilecol = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = ibtile;
    emxEnsureCapacity((emxArray__common *)y, jtilecol, (int)sizeof(double));
    for (jtilecol = 0; jtilecol < ibtile; jtilecol++) {
      y->data[y->size[0] * jtilecol] = xy_car->data[1 + xy_car->size[0] *
        jtilecol];
    }

    emxFree_real_T(&xy_car);
  }
}

/*
 * File trailer for otg_xy_reallyDumb.c
 *
 * [EOF]
 */
