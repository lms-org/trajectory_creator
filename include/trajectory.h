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
    bool doesCollide(); //does the trajectory collide
    bool isDrivable(); //is the traj. drivable

    T kappa_max = 1.5; //max. drivable curvature of the traj.
    T aOrthMax = 0.5*9.81; //max. acc. orthogonal to the traj.

    T ctot(); //function: gives back the value of the total cost function

    template<size_t m>
    points2d<m> sampleXY()
    {

        points2d<m> points;

        T dt = tend/(m-1); //time steps
        Vector<m> tt;
        for (int i = 0; i < m; i++)
        {
            tt(i) = dt*i;
        }

        Vector<m> ss;
        Vector<m> dd;

        ss = mPtr_s->eval<m>(tt);
        dd = mPtr_d->eval<m>(tt);

        Vector<2> XY;
        Vector<2> vec_help;

        Matrix<2, 2> R; //rot. matrix
        R(0,0) =  cos(this->roadData1.phi);
        R(0,1) = -sin(this->roadData1.phi);
        R(1,0) =  sin(this->roadData1.phi);
        R(1,1) =  cos(this->roadData1.phi);

        // generate the XY points
        if (this->kappa > 0)
        {
            for (int i = 0; i <m; i++)
            {
                // Position of the center line in some other XY coordinate system
                XY(0) = sin(this->kappa*ss(i))/this->kappa  - sin(this->kappa*ss(i))*dd(i);
                XY(1) = (1-cos(this->kappa*ss(i)))/this->kappa +  cos(this->kappa*ss(i))*dd(i);

                vec_help = R*XY;

                points.x(i) = vec_help(0) + sin(this->roadData1.phi)*this->D.d0;
                points.y(i) = vec_help(1) - cos(this->roadData1.phi)*this->D.d0;

            }

        } else{
            for (int i = 0; i <m; i++)
            {
                // don't know exactly what this is
                XY(0) = ss(i);
                XY(1) = dd(i);

                vec_help = R*XY;

                //stays the same
                points.x(i) = vec_help(0) + sin(this->roadData1.phi)*this->D.d0;
                points.y(i) = vec_help(1) - cos(this->roadData1.phi)*this->D.d0;

            }


        }

        return points;
    }

    /**
     * empty constructor
     */
    Trajectory()
    {

    };

    Trajectory(const T& _v1, const T& _d1, const T & _safetyS, const T& _safetyD, const roadData& _roadData1, const std::vector<obstacleData>& _obstacles, T _tend, const coeffCtot& _coeffCtot);


    /**
     * Output to iostream
     */
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


    Poly<5>* mPtr_d; //pointer to the poly. object for the d (lateral) direction in Frenet space
    Poly<4>* mPtr_s; //pointer to the poly. object for the s (longitudinal) direction in Frenet space

    coeffCtot coeffCtot1; //coefficients of the cost function (see also types.h)

    T tend; //time the traj. needs (as tstart = 0 by definition)

    T kappa; //curvature of the road (approx. as a circle)

    T ctot_value; //value of the total cost function
    bool ctot_alreadyCalc; //was the value already calculated

    roadData roadData1; //data of the road (see also types.h)

    std::vector<obstacleData> obstacles; //vector: each obstacleData is the data corresp. to one obstacle (see also types.h)

    bool collision = true; //is there a collision
    bool collisionDetected = false; //was collision already checked

    bool drivable = false; //is the traj. physically (no obstacles considered)
    bool drivabilityDetected = false; //was drivability already calculated

    int nSamplesCollisionAndDrivabilityDetection = 300;


    T safetyS; //safety distance in s (longitudinal) direction: IMPORTANT: Both cars are assumed as points so this must at least include half their lengths
    T safetyD; //safety distance in d (lateral) direction: IMPORTANT: this is measured from the center line

    T Jtd(); //function that returns the smoothness-functional for the d direction

    T Jts(); //function that returns the smoothness-functional for the s direction

    S_initialCond S;
    D_initialCond D;

};


#endif //PROJECT_TRAJECTORY_H
