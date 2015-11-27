//
// Created by Lukas Koestler on 23.10.15.
//

#ifndef PROJECT_TRAJECTORY_H
#define PROJECT_TRAJECTORY_H

#include "types.h"
#include "poly.h"
#include <vector>
#include <types.h>
#include <iomanip>

class Trajectory {
public:
    bool doesCollide();
    bool isDrivable();

    T kappa_max = 1.5;
    T aOrthMax = 0.5*9.81;

    T ctot();

    template<size_t m>
    points2d<m> sampleXY();

    Trajectory()
    {

    };

    Trajectory(const T& _v1, const T& _d1, const T & _safetyS, const T& _safetyD, const roadData& _roadData1, const std::vector<obstacleData>& _obstacles, T _tend, const coeffCtot& _coeffCtot);


    friend std::ostream& operator<<(std::ostream& os,  Trajectory& trajectory)
    {
        os << "Trajectory with the following propteries:" << std::endl
        << "time for Trajectory: " << trajectory.tend << std::endl
        << "Polynomials in Frenet space:" <<std::endl
        << "s: \t" << *trajectory.mPtr_s << std::endl
        << "d: \t" << *trajectory.mPtr_d << std::endl
        << "total value of the cost function: " << std::setprecision(30) << trajectory.ctot() <<std::endl
        << "this trajectory is colliding:" << trajectory.doesCollide() << std::endl
        << "this traj. is physically drivable: " << trajectory.drivable << std::endl;
        return os;

    }

private:
    /*std::unique_ptr<Poly<5>> mPtr_d;
    std::unique_ptr<Poly<4>> mPtr_s;*/


    Poly<5>* mPtr_d;
    Poly<4>* mPtr_s;

    coeffCtot coeffCtot1;

    T tend;

    T kappa;

    T ctot_value;
    bool ctot_alreadyCalc;

    roadData roadData1;

    std::vector<obstacleData> obstacles;

    bool collision;
    bool collisionDetected;

    bool drivable;
    bool drivabilityDetected;

    int numberOfsamples = 300;


    T safetyS;
    T safetyD;

    T Jtd();

    T Jts();


};


#endif //PROJECT_TRAJECTORY_H
