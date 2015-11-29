#include <types.h>
#include <lms/math/vertex.h>
#include <BezierPolynomial.h>
#include "trajectory_generator.h"
#include "lms/math/math.h"
/**
 * creates a trajectory using the sampling method both for different trajectories in end time and to check for collision/drivability
 */
bool TrajectoryGenerator::createTrajectorySample(Trajectory &trajectory,T v1, T d1, T safetyS, T safetyD, T tmin, T tmax, int nSamplesTraj, const RoadData& roadDataIn, std::vector<ObstacleData>& obstacleDataIn, const CoeffCtot& coeffCtotIn) {


    T dt = (tmax - tmin) / (nSamplesTraj - 1); //time increment

    //find the minimum cost function while still drivable and collision free
    bool flagFound = false;
    T bestCtot = 0;

    for (int i = 0; i < nSamplesTraj; i++) {
        //create trajectory
        Trajectory sampleTrajectory = Trajectory(v1, d1, safetyS, safetyD, roadDataIn, obstacleDataIn, tmin + i * dt, coeffCtotIn);
        if (sampleTrajectory.isDrivable() && !sampleTrajectory.doesCollide()) {
            // drivable and not colliding
            if (!flagFound) {
                // the first one found
                flagFound = true;
                bestCtot = sampleTrajectory.ctot();
                trajectory = sampleTrajectory;
            } else {
                // not the first --> compare ctot
                if (sampleTrajectory.ctot() <= bestCtot) {
                    // better --> change
                    bestCtot = sampleTrajectory.ctot();
                    trajectory = sampleTrajectory;
                }
            }

        }
    }

    if (flagFound) {
        // trajectory was already set
        return true;
    }
    else {
        return false;
    }

}

TrajectoryGenerator::TrajectoryGenerator(lms::logging::Logger& _logger) : logger(_logger){
    // Test of the Bezier Polynomials

    const size_t n = 2;
    Vector<3> controlPointsIn;
    controlPointsIn(0) = 1; //some random numbers
    controlPointsIn(1) = 3;
    controlPointsIn(2) = -1;


    T t_begin = 0;
    T t_end = 1;

    const size_t m = 100;

    Vector<m> tt;

    T dt = (t_end-t_begin)/(m-1);

    for (size_t i = 0; i < m; i++)
    {
        tt(i) = t_begin + i*dt;
    }

    BezierPolynomial<2> bezPoly = BezierPolynomial<2>(controlPointsIn, t_begin, t_end);
    BezierPolynomial<1> DerBezPoly = bezPoly.differentiate();

    T valueAt0p5 = bezPoly.evalAtPoint(0.5);

    auto y = bezPoly.eval<m>(tt);
    auto dy = DerBezPoly.eval<m>(tt);


    std::cout << "Time points" << tt << std::endl;
    std::cout << "y: " << y << std::endl;
    std::cout << "dy: " << dy << std::endl;
}

float TrajectoryGenerator::circleCurvature(lms::math::vertex2f p1, lms::math::vertex2f p2, lms::math::vertex2f p3)
{
    // look at Arndt Brunner for explanation: http://www.arndt-bruenner.de/mathe/scripts/kreis3p.htm

    float kappa_est = 0;

    //set up A matrix
    Matrix<3,3> A;
    A << 1, -p1.x, -p1.y,
        1, -p2.x, -p2.y,
        1, -p3.x, -p3.y;

    Vector<3> b;
    b <<    -(pow(p1.x, 2) + pow(p1.y,2)),
            -(pow(p2.x, 2) + pow(p2.y,2)),
            -(pow(p3.x, 2) + pow(p3.y,2));

    Vector<3> sol;

    sol = A.colPivHouseholderQr().solve(b);

    T xm = sol(1)/2;
    T ym = sol(2)/2;
    T r = sqrt(pow(xm,2) + pow(ym,2) - sol(0));

    if (r == 0)
    {
        //throw some kind of error here
        kappa_est = 0;
    }else
    {
        kappa_est = lms::math::sgn<float>(ym)*(float)1/r;


    }

    return kappa_est;

}


