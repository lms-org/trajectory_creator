//
// Created by Lukas Koestler on 23.10.15.
//

#ifndef PROJECT_TRAJECTORY_H
#define PROJECT_TRAJECTORY_H

#include "types.h"
#include "poly.h"
#include "BezierCurve.h"
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
        if (fabs(this->kappa) > 0)
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
     * Projects the trajectory on the unique Bezier Curve that is defined by the k points in points (points2d<k> struct points)
     * m = number of points to output
     * l = diastance between two successive points
     */
    template<size_t m, size_t k>
    points2d<m> projectOntoBezierCurve(const points2d<k> pointsIn, const T l)
    {
        T s_begin = 0;
        T s_end = l*(k-1); //as there are (k-1) segemnts of length l between k points

        BezierCurve<k-1> centerLine = BezierCurve<k-1>(pointsIn, s_begin, s_end);

        //generate the vector with the times
        points2d<m> pointsOut;

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

        for(size_t i = 0; i < m; i++)
        {
            //different cases
            if (ss(i) < 0)
            {
                // this should not be: throw error
                pointsOut.x(i) = 0;
                pointsOut.y(i) = 0;


            }else if(ss(i) >= 0 && ss(i) <= s_end)
            {
                //point on the center line
                Vector<2> centerLinePoint = centerLine.evalAtPoint(ss(i));

                //normal
                Vector<2> centerLineNormal = centerLine.normalAtPoint(ss(i));

                Vector<2> trajPoint = centerLinePoint + dd(i)*centerLineNormal; //the trajectoray point is the centerLinePoint plus d times the normal (should be oriented the right way and normalized)

                pointsOut.x(i) = trajPoint(0);
                pointsOut.y(i) = trajPoint(1);


            } else if (ss(i) > s_end)
            {
                //the point is behind the trajectory --> linear extrapolation
                //point on the center line
                Vector<2> centerLinePoint = centerLine.evalAtPoint(s_end);

                //normal
                Vector<2> centerLineNormal = centerLine.normalAtPoint(s_end);

                //normal
                Vector<2> centerLineTangent = centerLine.tangentAtPoint(s_end);


                Vector<2> trajPoint = centerLinePoint + (ss(i)-s_end)*centerLineTangent + dd(i)*centerLineNormal;

                pointsOut.x(i) = trajPoint(0);
                pointsOut.y(i) = trajPoint(1);

            } else
            {
                // case not considered: use zeros
                // this should not be: throw error
                pointsOut.x(i) = 0;
                pointsOut.y(i) = 0;
            }
        }


        return pointsOut;

    }




    /**
     * empty constructor
     */
    Trajectory()
    {

    }

    Trajectory(const T& _v1, const T& _d1, const T & _safetyS, const T& _safetyD, const RoadData& _roadData1, const std::vector<ObstacleData>& _obstacles, T _tend, const CoeffCtot& _coeffCtot);


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

    CoeffCtot coeffCtot1; //coefficients of the cost function (see also types.h)

    T tend; //time the traj. needs (as tstart = 0 by definition)

    T kappa; //curvature of the road (approx. as a circle)

    T ctot_value; //value of the total cost function
    bool ctot_alreadyCalc; //was the value already calculated

    RoadData roadData1; //data of the road (see also types.h)

    std::vector<ObstacleData> obstacles; //vector: each obstacleData is the data corresp. to one obstacle (see also types.h)

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
