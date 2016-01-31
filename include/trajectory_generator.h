#ifndef TRAJECTORYGENERATOR_H
#define TRAJECTORYGENERATOR_H

#include <lms/logging/logger.h>

#include "types.h"
#include "poly.h"
#include "trajectory.h"

class TrajectoryGenerator {
public:

    static float circleCurvature(lms::math::vertex2f p1, lms::math::vertex2f p2, lms::math::vertex2f p3);

    /**
     * @brief Constructor
     */
    TrajectoryGenerator(lms::logging::Logger& _logger);
    /**
     *creates a trajectory using the sampling method both for different trajectories in end time and to check for collision/drivability
     * TODO beschreibung der parameter
     * @brief createTrajectorySample
     * @param trajectory
     * @param v1
     * @param d1
     * @param safetyS
     * @param safetyD
     * @param tmin
     * @param tmax
     * @param nSamplesTraj
     * @param roadDataIn
     * @param obstacleDataIn
     * @param coeffCtotIn
     * @return
     */
    bool createTrajectorySample(Trajectory &trajectory,T v1, T d1, T safetyS, T safetyD, T tmin, T tmax, int nSamplesTraj, const RoadData& roadDataIn, std::vector<ObstacleData>& obstacleDataIn, const CoeffCtot& coeffCtotIn);


protected:
    lms::logging::Logger logger;
};


#endif //TRAJECTORYGENERATOR_H
