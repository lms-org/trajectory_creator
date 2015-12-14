#ifndef IMAGE_HINT_TRANSFORMER_H
#define IMAGE_HINT_TRANSFORMER_H

#include "lms/module.h"
#include "lms/math/polyline.h"
#include "street_environment/road.h"
#include "sensor_utils/car.h"
#include "trajectory_generator.h"
#include "street_environment/trajectory.h"

class TrajectoryLineCreator : public lms::Module {
public:
    bool initialize() override;
    bool deinitialize() override;
    bool cycle() override;
private:
    street_environment::Trajectory simpleTrajectory(float trajectoryMaxLength,const int obstacleTrustThreshold);
    bool advancedTrajectory(lms::math::polyLine2f &trajectory);
    lms::ReadDataChannel<street_environment::EnvironmentObjects> envObstacles;
    lms::ReadDataChannel<street_environment::RoadLane> road;
    lms::ReadDataChannel<sensor_utils::Car> car;
    lms::WriteDataChannel<street_environment::Trajectory> trajectory;

    TrajectoryGenerator* generator;

};

#endif /* IMAGE_HINT_TRANSFORMER_H */
