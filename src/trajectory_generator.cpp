#include <types.h>
#include <curses.h>
#include "trajectory_generator.h"
/**
 * creates a trajectory using the sampling method both for different trajectories in end time and to check for collision/drivability
 */
bool trajectory_generator::create_trajectory_sample(Trajectory &trajectory,T v1, T d1, T safetyS, T safetyD, T tmin, T tmax, int nSamplesTraj, const roadData& roadDataIn, std::vector<obstacleData>& obstacleDataIn, const coeffCtot& coeffCtotIn) {


    T dt = (tmax - tmin) / (nSamplesTraj - 1); //time increment

    //vector containing the trajectories
    std::vector<Trajectory> sampleTrajectories;

    /* //not good in performance
    for (int i = 0; i < nSamplesTraj; i++) {
        sampleTrajectories.push_back(
                Trajectory(v1, d1, safetyS, safetyD, roadDataIn, obstacleDataIn, tmin + i * dt, coeffCtotIn));
    }*/

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

trajectory_generator::trajectory_generator(lms::logging::Logger& _logger) : logger(_logger)
{
    /*
    Vector<8> x;
    x << 0, 0, 0, 0, 1, 2, 500, 100;

    Vector<3> coeff;
    coeff << 1, 2, 3;

    Vector<8> y;

    Poly<2>* poly1 = new Poly<2>(coeff);

    y = poly1->eval<8>(x);

    Poly<3> poly2 = poly1->integrate();
    Poly<1> poly3 = poly1->differentiate();

    std::cout << y << "\n" << "integrated:  " << poly2.getCoeff().transpose() << "\n" << "differentiated:  " << poly3.getCoeff().transpose() << std::endl;


    // test of the trajectory constructor

    T v1 = 2;
    T d1 = 0.2;

    T tend = 10;

    roadData roadDataTest;
    roadDataTest.ax0 = 0;
    roadDataTest.kappa = 0.01;
    roadDataTest.phi = 0;
    roadDataTest.vx0 = 1;
    roadDataTest.w = 0;
    roadDataTest.y0 = 0.25;

    coeffCtot coeffCtotTest;
    coeffCtotTest.kj = 5;
    coeffCtotTest.kT = 2;
    coeffCtotTest.ks = 20;
    coeffCtotTest.kd = 100;

    std::vector<obstacleData> obstaclesEmpty;

    Trajectory* ptr_testTraj = new Trajectory(v1, d1, roadDataTest, obstaclesEmpty, tend, coeffCtotTest);

    std::cout << poly2 << std::endl;

    std::cout << std::endl;

    std::cout << "---------- test fuer trajectory   ----------" <<std::endl << *ptr_testTraj << std::endl;
    std::cout << "---------- test fuer trajectory 2 ----------" <<std::endl << *ptr_testTraj << std::endl;
     */

    // -----------------------------------------------------
    // load DATA from confiq. This is emulated at the moment
    // -----------------------------------------------------

    T v1 = 2;
    T d1 = 0.2;

    T safetyS = 0.55;
    T safetyD = 0.15;

    roadData roadDataIn;
    roadDataIn.ax0 = 0;
    roadDataIn.kappa = 0.25;
    roadDataIn.phi = 0.1;
    roadDataIn.vx0 = 1;
    roadDataIn.w = 0;
    roadDataIn.y0 = 0.25;

    coeffCtot coeffCtotIn;
    coeffCtotIn.kj = 5;
    coeffCtotIn.kT = 2;
    coeffCtotIn.ks = 20;
    coeffCtotIn.kd = 100;

    obstacleData obstacleData1;
    obstacleData1.leftLane = FALSE;
    obstacleData1.s0 = 2;
    obstacleData1.v0 = 0;

    std::vector<obstacleData> obstacleDataIn;

    obstacleDataIn.push_back(obstacleData1);

    Trajectory bestTrajectory1;
    if(create_trajectory_sample(bestTrajectory1,v1, d1, safetyS, safetyD, 0.1, 4, 390, roadDataIn, obstacleDataIn, coeffCtotIn))
    {
        std::cout << "---------- test fuer trajectory   ----------" <<std::endl << bestTrajectory1 << std::endl;

        std::cout << "aOrthMAx:     " << bestTrajectory1.aOrthMax << std::endl;
        std::cout << "kappaMax:     " << bestTrajectory1.kappa_max << std::endl;

        std::cout << "does the best traj. collide:   " << bestTrajectory1.doesCollide() << std::endl;

        points2d<10> points = bestTrajectory1.sampleXY<10>();

        std::cout << "POINTS: x coordinates:  " << points.x << std::endl;
        std::cout << "POINTS: y coordinates:  " << points.y << std::endl;
    } else{
        std::cout<<"No valid trajectory found"<<std::endl;
    }
}


