//
// Created by Lukas Koestler on 23.10.15.
//

#include <types.h>
#include "trajectory.h"


bool Trajectory::doesCollide() {

    if (collisionDetected == true)
    {
        return collision;

    }else {

        collisionDetected = true;

        T dt = tend / numberOfsamples;
        T t = 0;

        T s_local;
        T d_local;

        for (int i = 0; i <= numberOfsamples; i++) {
            s_local = mPtr_s->evalAtPoint(t);
            d_local = mPtr_d->evalAtPoint(t);


            for (auto const &obs:obstacles) {
                if (s_local >= -safetyS + obs.s0 + obs.v0 * t) {
                    //std::cout << "v0 of obstacle  " << obs.v0 << std::endl;
                    if (obs.leftLane) {
                        if (d_local >= -safetyD) {
                            collision = true;
                            //std::cout << "There is something colliding" << std::endl;
                            return true;
                        }
                    } else {
                        //std::cout << "coll. det. right lane" << std::endl;
                        if (d_local <= safetyD) {
                            collision = true;
                            //std::cout << "There is something colliding" << std::endl;
                            return true;
                        }
                    }
                }
            }

            t = t + dt;

        }

        collision = false;
        return false;
    }
}

bool Trajectory::isDrivable() {
    if (drivabilityDetected == true)
    {
        return drivable;
    } else
    {
        drivabilityDetected = true;

        T dt = tend / numberOfsamples;
        T t = 0;

        // differentiate the polynomials
        Poly<4> poly_d_d = mPtr_d->differentiate();
        Poly<3> poly_d_dd = poly_d_d.differentiate();

        Poly<3> poly_s_d = mPtr_s->differentiate();

        for (int i = 0; i <= numberOfsamples; i++) {
            // evaluate the polynomials for the local time
            T d_d_local  = poly_d_d.evalAtPoint(t);
            T d_dd_local = poly_d_dd.evalAtPoint(t);

            T s_d_local = poly_s_d.evalAtPoint(t);

            // update time
            t = t + dt;

            T kappa_xy_local;


            if (s_d_local == 0)
            {
                if (d_dd_local == 0)
                {
                    // 0 divided by 0 should be 0 here: so ok
                    kappa_xy_local = 0;
                }else
                {
                    // to not divide by 0
                    drivable = false;
                    return drivable;
                }

            } else
            {
                // normal calculation for kappa_xy_local
                kappa_xy_local = kappa + d_dd_local/s_d_local;

                // Check if the local curvature is bigger than the max. curvature drivable (speed independent)
                if (abs(kappa_xy_local) > kappa_max)
                {
                    drivable = false;
                    return drivable;
                }

                // check if the orthogonal acceleration is ok (speed dependent)
                T v_local_squarred = pow(d_d_local,2) + pow(s_d_local,2);

                if (abs(kappa_xy_local * v_local_squarred) > aOrthMax)
                {
                    drivable = false;
                    return drivable;
                }
            }

        }
        drivable = true;
        return drivable;
    }
}

T Trajectory::ctot() {
    if (ctot_alreadyCalc) {
        return ctot_value;
    }
    else {
        ctot_alreadyCalc = true;
        ctot_value = coeffCtot1.ks * (coeffCtot1.kj * this->Jts() + coeffCtot1.kT * this->tend) +
                     coeffCtot1.kd * (coeffCtot1.kj * this->Jtd() + coeffCtot1.kT * this->tend);
        return ctot_value;
    }
}


T Trajectory::Jtd() {
    const Vector<6> &coeff_d = mPtr_d->getCoeff();
    return 36 * tend * pow(coeff_d(3), 2) + pow(tend, 3) * (192 * pow(coeff_d(4), 2) + 240 * coeff_d(3) * coeff_d(5)) +
           720 * pow(tend, 5) * pow(coeff_d(5), 2) + 144 * pow(tend, 2) * coeff_d(3) * coeff_d(4) +
           720 * pow(tend, 4) * coeff_d(4) * coeff_d(5);
}

T Trajectory::Jts() {
    const Vector<5> &coeff_s = mPtr_s->getCoeff();
    return 192 * pow(tend, 3) * pow(coeff_s(4), 2) + 144 * pow(tend, 2) * coeff_s(3) * coeff_s(4) +
           36 * tend * pow(coeff_s(3), 2);
}


Trajectory::Trajectory(const T &_v1, const T &_d1, const T & _safetyS, const T& _safetyD, const roadData &_roadData1, const std::vector<obstacleData> &_obstacles, T _tend, const coeffCtot& _coeffCtot) : roadData1(_roadData1),
                                                                                                                                                                                                           obstacles(_obstacles),
                                                                                                                                                                                                           tend(_tend),
                                                                                                                                                                                                           coeffCtot1(_coeffCtot),
                                                                                                                                                                                                           safetyS(_safetyS),
                                                                                                                                                                                                           safetyD(_safetyD)

{
    //init all variables
    ctot_alreadyCalc = false;

    // convert roadData to initial conditions in Frenet space
    S_initialCond S;
    S.v0 = cos(roadData1.phi) * roadData1.vx0;
    S.a0 = -sin(roadData1.phi) * roadData1.w * roadData1.vx0 + cos(roadData1.phi) * roadData1.ax0;
    S.v1 = _v1;

    D_initialCond D;
    double cosphi_d = -sin(roadData1.phi) * roadData1.w;
    double cosphi_dd = -pow(roadData1.w, 2) * cos(roadData1.phi);
    double y0_d = -sin(roadData1.phi) * roadData1.vx0;
    double y0_dd = -(cos(roadData1.phi) * roadData1.w * roadData1.vx0 + sin(roadData1.phi) * roadData1.ax0);

    D.d0 = -cos(roadData1.phi) * roadData1.y0;
    D.d0d = -(cosphi_d * roadData1.y0 + cos(roadData1.phi) * y0_d);
    D.d0dd = -(cosphi_dd * roadData1.y0 + 2 * cosphi_d * y0_d + cos(roadData1.phi) * y0_dd);
    D.d1 = _d1;

    // init polynomial s
    Vector<5> coeff_s;
    // precomputed formulas
    coeff_s << 0, S.v0, S.a0 / 2, -pow(1 / tend, 2) * (3 * S.v0 - 3 * S.v1 + 2 * tend * S.a0) / 3, pow(1 / tend, 3) *
                                                                                                   (2 * S.v0 -
                                                                                                    2 * S.v1 +
                                                                                                    tend * S.a0) / 4;

    Vector<6> coeff_d;
    // precomputed formulas
    coeff_d << D.d0, D.d0d, D.d0dd / 2, -pow(1 / tend, 3) *
                                        (3 * D.d0dd * pow(tend, 2) + 12 * D.d0d * tend + 20 * D.d0 - 20 * D.d1) / 2,
            pow(1 /
                tend, 4) * (3 * D.d0dd * pow(
                    tend, 2) + 16 * D.d0d * tend + 30 * D.d0 - 30 * D.d1) / 2, -pow(1 / tend, 5) *
                                                                               (D.d0dd * pow(tend, 2) +
                                                                                6 * D.d0d * tend + 12 * D.d0 -
                                                                                12 * D.d1) / 2;

    //Init polynomials
    /*mPtr_s = std::unique_ptr<Poly<4>>(new Poly<4>(coeff_s));
    mPtr_d = std::unique_ptr<Poly<5>>(new Poly<5>(coeff_d));*/

    mPtr_s = new Poly<4>(coeff_s);
    mPtr_d = new Poly<5>(coeff_d);
}


template<size_t m>
points2d<m> Trajectory::sampleXY() {
    return points2d<m>();
}
