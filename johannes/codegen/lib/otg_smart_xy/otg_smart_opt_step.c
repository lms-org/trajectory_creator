/*
 * File: otg_smart_opt_step.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "otg_smart_opt_step.h"
#include "otg_smart_objFun.h"

/* Function Definitions */

/*
 * Arguments    : const boolean_T x[3]
 *                int y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
void eml_li_find(const boolean_T x[3], int y_data[], int y_size[2])
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
void otg_smart_opt_step(const double Ts[3], double Cs[3], const double notD[3],
  const double coll[3], const double S[3], const double D[4], double kj, double
  kT, double ks, double kd, const double dataVeh[3], double safetyS, double
  safetyD, double kappa, double kappaMax, double aOrthMax, double *flag1, double
  *flag2, double Ts_new[3], double Cs_new[3], double notD_new[3], double
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
 * File trailer for otg_smart_opt_step.c
 *
 * [EOF]
 */
