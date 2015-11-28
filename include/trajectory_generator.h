#ifndef TRAJECTORYGENERATOR_H
#define TRAJECTORYGENERATOR_H

#include <lms/logging/logger.h>

#include "types.h"
#include "poly.h"
#include "trajectory.h"

class trajectory_generator {
public:

    /**
     * @brief Constructor
     */
    trajectory_generator(lms::logging::Logger& _logger);

protected:
    lms::logging::Logger logger;

    bool create_trajectory_sample(Trajectory &trajectory,T v1, T d1, T safetyS, T safetyD, T tmin, T tmax, int nSamplesTraj, const roadData& roadDataIn, std::vector<obstacleData>& obstacleDataIn, const coeffCtot& coeffCtotIn);
};


#endif //TRAJECTORYGENERATOR_H
