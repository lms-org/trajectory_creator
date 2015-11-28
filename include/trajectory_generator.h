#ifndef TRAJECTORYGENERATOR_H
#define TRAJECTORYGENERATOR_H

#include <lms/logging/logger.h>

#include "types.h"
#include "poly.h"
#include "trajectory.h"

class TrajectoryGenerator {
public:

    /**
     * @brief Constructor
     */
    TrajectoryGenerator(lms::logging::Logger& _logger);
    bool createTrajectorySample(Trajectory &trajectory,T v1, T d1, T safetyS, T safetyD, T tmin, T tmax, int nSamplesTraj, const RoadData& roadDataIn, std::vector<ObstacleData>& obstacleDataIn, const CoeffCtot& coeffCtotIn);


protected:
    lms::logging::Logger logger;
};


#endif //TRAJECTORYGENERATOR_H
