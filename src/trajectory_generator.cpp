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

    /*const size_t n = 2;
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
    std::cout << "dy: " << dy << std::endl;*/

    //Test of the whole new algo with Bezier Curves

    points2d<10> pointsCenter;
    //pointsCenter.x << 0, 0.1995, 0.3966, 0.5797, 0.7582, 0.9121, 1.0686, 1.2260, 1.3857, 1.5638;
    //pointsCenter.x << 0    ,0.1999,    0.3988,    0.5936,    0.7814,    0.9583,    1.1364,    1.2948,    1.4168,    1.5274;
    pointsCenter.x << 0,    0.1960,    0.3802,    0.5453,    0.6846,    0.7927,    0.8652,    0.8992,    0.8933,    0.8479;
    //pointsCenter.y << 0.2500, 0.2637, 0.2977, 0.3783, 0.4684, 0.5961, 0.7206, 0.8441, 0.9644, 1.0555;
    //pointsCenter.y << 0.2500,    0.2577,    0.2778,    0.3232,    0.3921,    0.4853,    0.5763,    0.6985,    0.8570,    1.0236;
    pointsCenter.y << 0.2500,    0.2897,    0.3676,    0.4805,    0.6240,    0.7923,    0.9787,    1.1758,    1.3757,    1.5705;

    // Generate Data
    T v1 = 2;
    T d1 = 0.2;

    T safetyS = 0.55;
    T safetyD = 0.15;

    T tmin = 1;
    T tmax = 5;

    int nSamplesTraj = 100;

    RoadData roadDataIn;
    roadDataIn.vx0 = 1;
    roadDataIn.ax0 = 0;
    roadDataIn.phi = 0;
    roadDataIn.w = 0;
    roadDataIn.y0 = 0.25;
    roadDataIn.kappa = 0.05;

    RoadData roadDataCenter;
    roadDataCenter.vx0 = 1;
    roadDataCenter.ax0 = 0;
    roadDataCenter.phi = 0;
    roadDataCenter.w = 0;
    roadDataCenter.y0 = 0;
    roadDataCenter.kappa = 0.05;

    ObstacleData obstacle1;
    obstacle1.s0 = 2;
    obstacle1.v0 = 0;
    obstacle1.leftLane = false;

    std::vector<ObstacleData> obstaclesIn;
    obstaclesIn.push_back(obstacle1);

    std::vector<ObstacleData> noObs;

    CoeffCtot coeffCtotIn;
    coeffCtotIn.kj = 5;
    coeffCtotIn.kT = 2;
    coeffCtotIn.ks = 20;
    coeffCtotIn.kd = 100;

    Trajectory trajectory;

    T l = 0.2;

    Trajectory trajCenterLine = Trajectory(2, 0, 0, 0, roadDataCenter, noObs, 4, coeffCtotIn);

    BezierCurve<9> centerLine = BezierCurve<9>(pointsCenter, 0, 9*l);

    if(TrajectoryGenerator::createTrajectorySample(trajectory, v1, d1, safetyS, safetyD, tmin, tmax, nSamplesTraj, roadDataIn, obstaclesIn, coeffCtotIn))
    {
        std::cout << "---------- SUCCESS ----------" << std::endl;
        // do rest here

        points2d<100> pointsOut = trajectory.projectOntoBezierCurve<100>(pointsCenter, l);
        std::cout << "---x: " << std::endl;
        std::cout << pointsOut.x << std::endl << std::endl;
        std::cout << "---y: " << std::endl;
        std::cout << pointsOut.y << std::endl;

        std::cout  << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

        Vector<100> tt;

        for(int i = 0; i <100; i++)
        {
            tt(i) = i*(9*0.2/99);
        }

        points2d<100> pointsCenter = centerLine.eval<100>(tt);

        Vector<100> curvature = centerLine.curvature<100>(tt);

        std::cout << "curvature: " << std::endl << curvature << std::endl;

        std::cout << "x center: " << std::endl << pointsCenter.x << std::endl;
        std::cout << "y center: " << std::endl << pointsCenter.y << std::endl;

    }else
    {
        std::cout << "---------- FAILURE ----------" << std::endl;
    }

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


