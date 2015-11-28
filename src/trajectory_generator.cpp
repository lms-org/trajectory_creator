#include <types.h>
#include "trajectory_generator.h"
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
}


