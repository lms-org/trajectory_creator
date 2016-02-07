#include "trajectory_line_creator.h"
#include "lms/math/math.h"
#include "street_environment/obstacle.h"
#include "street_environment/crossing.h"
#include "lms/math/mathEigen.h"
bool TrajectoryLineCreator::initialize() {
    envObstacles = readChannel<street_environment::EnvironmentObjects>("ENVIRONMENT_OBSTACLE");
    roadStates= readChannel<street_environment::RoadStates>("ROAD_STATES");
    road = readChannel<street_environment::RoadLane>("ROAD");
    trajectory = writeChannel<street_environment::Trajectory>("LINE");
    debug_trajectory = writeChannel<lms::math::polyLine2f>("DEBUG_TRAJECTORY");
    debug_trajectory2 = writeChannel<lms::math::polyLine2f>("DEBUG_TRAJECTORY_2");
    car = readChannel<sensor_utils::Car>("CAR");

    generator = new TrajectoryGenerator(logger);

    return true;
}

bool TrajectoryLineCreator::deinitialize() {
    delete generator;
    return true;
}

float TrajectoryLineCreator::targetVelocity(){
    const float obstacleTrustThreshold = config().get<float>("obstacleTrustThreshold",0.5);
    float velocity = 0;
    bool obstacleInSight = false;
    float distanceToObstacle = 0;
    for(street_environment::EnvironmentObjectPtr obj:envObstacles->objects){
        if(obj->getType() != street_environment::Obstacle::TYPE)
            continue;
        street_environment::ObstaclePtr obst = std::static_pointer_cast<street_environment::Obstacle>(obj);
        //Only looking for obstacles on the right side
        if(obst->distanceOrth() < 0.1 && obst->trust() > obstacleTrustThreshold){
            obstacleInSight = true;
            if(obst->position().x > -0.1 && (distanceToObstacle > obst->distanceTang())){
                distanceToObstacle = obst->distanceTang();
            }
        }
    }

    Eigen::Vector3f stateVelocities;
    stateVelocities(0) = config().get<float>("velocity_straight", 6);

    float aOrthMax = config().get<float>("aOrthMax", 9.81*0.5);
    float curve_minVelocity = config().get<float>("curve_minVelocity", 1.8);
    float curve_maxVelocity = config().get<float>("curve_maxVelocity", 5);
    float curveStraight_minVelocity = config().get<float>("curveStraight_minVelocity", 1.5);
    float curveStraight_maxVelocity = config().get<float>("curveStraight_maxVelocity", 2.0);

    float curveVelocity = sqrt(aOrthMax/fabs(roadStates->states[2].curvature));
    float straightCurveVelocity = sqrt(aOrthMax/fabs(roadStates->states[1].curvature));

    if (roadStates->states[2].curvature == 0)
    {
        curveVelocity = curve_maxVelocity;
    }
    if (roadStates->states[1].curvature == 0)
    {
        straightCurveVelocity = straightCurveVelocity;
    }


    if (curveVelocity < curve_minVelocity)
    {
        curveVelocity = curve_minVelocity;
    }
    if (curveVelocity > curve_maxVelocity)
    {
        curveVelocity = curve_maxVelocity;
    }
    if (straightCurveVelocity < curveStraight_minVelocity)
    {
        straightCurveVelocity = curveStraight_minVelocity;
    }
    if (straightCurveVelocity > curveStraight_maxVelocity)
    {
        straightCurveVelocity = curveStraight_maxVelocity;
    }

    stateVelocities(1) = straightCurveVelocity;
    stateVelocities(2) = curveVelocity;

    Eigen::Vector3f stateProbabilities;
    stateProbabilities(0) = roadStates->states[0].probability;
    stateProbabilities(1) = roadStates->states[1].probability;
    stateProbabilities(2) = roadStates->states[2].probability;

    velocity = (stateProbabilities.cwiseProduct(stateVelocities)).sum() / stateProbabilities.sum();

    if (velocity > stateVelocities(0))
    {
        velocity = stateVelocities(0);
    }

    if (velocity < curve_minVelocity)
    {
        velocity = curve_minVelocity;
    }

    if(obstacleInSight){
        velocity = velocity * config().get<float>("obstacleVelocitySafetyFactor", 0.65);
    }else{
        //do nothing
    }

    return velocity;
}

bool TrajectoryLineCreator::cycle() {
    //clear old trajectory
    trajectory->clear();
    //calculate data for creating the trajectory    //calculate the speed without obstacles
    float velocity = targetVelocity();
    logger.debug("set velocity: ") << velocity;


    bool advancedTraj = false;
    street_environment::Trajectory traj;
    if(config().get<bool>("simpleTraj",true)){
        traj= simpleTrajectory(config().get<float>("distanceObstacleBeforeChangeLine",0),velocity);
    }else{
        traj= simpleTrajectory(0,velocity);
        //detect if we have to go left or right
        if(traj[traj.size()-1].velocity == 0){
            //Stop it, we won't go for an advancedTrajectory :)
            logger.warn("STOP");
        }else{
            bool initRightSide (traj[0].distanceToMiddleLane <= 0);
            for(int i = 1; i < (int)traj.size(); i++){
                //try to find change of the line
                bool rightSide =(traj[i].distanceToMiddleLane <= 0);
                if(initRightSide != rightSide){
                    //we change the line
                    //TODO we don't cover obstacles close to each other
                    initRightSide = rightSide;
                    break;
                }
            }
            float a = a;
            if(a == 0){

            }
            //What happens if the current car velocity and the end-velocity is close to zero?
            float endVelocity = traj[traj.size()-1].velocity;
            float averageVelocity = 0;
            float maxVelocity = 0;
            for(int i = 0; i < (int)traj.size(); i++){
                averageVelocity += traj[i].velocity;
                if(maxVelocity < traj[i].velocity){
                    maxVelocity = traj[i].velocity;
                }
            }
            if(maxVelocity < 0.1){ //HACK
                maxVelocity = 1;
            }
            averageVelocity /= traj.size();
            float minTime = 0.05; //TODO #IMPORTANT
            float maxTime = 1.5;
            if(endVelocity < 1){ //TODO HACK
                endVelocity = 1;
            }

            traj.clear();
            if(!advancedTrajectory(traj,initRightSide,endVelocity,minTime,maxTime)){
                logger.warn("advancedTrajectory")<<"FAILED";
                traj = simpleTrajectory(config().get<float>("distanceObstacleBeforeChangeLine",0),velocity);
            }else{
                advancedTraj = true;
            }
        }
    }
    *trajectory = traj;
    debug_trajectory->points().clear();
    debug_trajectory2->points().clear();


    for(street_environment::TrajectoryPoint &v:traj){
        if (advancedTraj)
        {
            debug_trajectory->points().push_back(v.position);
        }else{
            debug_trajectory2->points().push_back(v.position);
        }

    }

    return true;
}

bool TrajectoryLineCreator::advancedTrajectory(street_environment::Trajectory &trajectory, bool rightSide, float endVelocity, float tMin,float tMax){
    if(road->polarDarstellung.size() < 8){
        logger.error("invalid middle")<<road->polarDarstellung.size();
        return false;
    }


    //INPUT
    //hindernisse: vector Abstand-Straße,geschwindigkeit-Straße(absolut), -1 (links) +1 (rechts) (alle hindernisse hintereinander)
    //eigenes auto, vx,vy, dw -winkelgeschwindigkeit dw (zunächst mal 0)
    //vector mit x-koordinaten
    //vector mit y-koordinaten

    int nSamplesTraj = config().get<int>("nSamplesTraj",50);
    double d1; //Abstand zur Mittellinie
    if(rightSide)
        d1= -0.2;
    else
        d1= 0.2;

    //aktuelle daten der fahrspur und des autos relativ dazu
    RoadData dataRoad;
    dataRoad.ax0 = 0; //beschl. am anfang
    dataRoad.kappa = lms::math::circleCurvature(road->points()[2],road->points()[5],road->points()[7]);
    dataRoad.phi = road->polarDarstellung[1];
    float velocity = car->velocity();//Sollte nicht 0 sein, wegen smoothem start
    if(velocity < 1) //TODO HACK
        velocity = 1;
    dataRoad.vx0 = velocity;
    dataRoad.w = 0; //aktuelle winkelgeschwindigkeit
    dataRoad.y0 = road->polarDarstellung[0];
    logger.debug("kappa")<<dataRoad.kappa<< " distanceToMiddle: "<<d1;
    logger.debug("vx0 ") << dataRoad.vx0 << ",  v1 " << endVelocity;
    logger.debug("y0 ") << dataRoad.y0 << ",  phi " << dataRoad.phi;
    logger.debug("tmin ") << tMin << ",  tMAx " << tMax;


    float tangLength = 0;

    float s0_closest;
    bool leftLane_closest;
    bool obstacleInTheWay = false;

    std::vector<ObstacleData> dataObstacle;
    float obstacleTrustThreshold = config().get<float>("obstacleTrustThreshold",0.5);
    for(const street_environment::EnvironmentObjectPtr objPtr:envObstacles->objects){
        if(objPtr->trust() < obstacleTrustThreshold)
            continue;
        if(objPtr->getType() == street_environment::Obstacle::TYPE){

            street_environment::ObstaclePtr obstPtr = std::static_pointer_cast<street_environment::Obstacle>(objPtr);
            if (obstPtr->distanceTang()-tangLength < 0.15)
            {
                continue;
            }
            ObstacleData toAdd;
            toAdd.s0 = obstPtr->distanceTang();
            toAdd.v0 = 0;
            toAdd.leftLane = obstPtr->distanceOrth()> 0;
            dataObstacle.push_back(toAdd);
            logger.debug("obstacle: s0: ") << toAdd.s0 << ",  v0: " << toAdd.v0 << ",  leftLane: " << toAdd.leftLane;


            if ((!obstacleInTheWay) || toAdd.s0 < s0_closest)
            {
                s0_closest = toAdd.s0;
                leftLane_closest = toAdd.leftLane;
            }
            obstacleInTheWay = true;
        }
    }
    CoeffCtot tot;

    // generate d1
    if (!obstacleInTheWay)
    {
        d1 = -0.2;
    }else
    {
        logger.debug("obstacle in the way: s0 = ") << s0_closest << ",  leftLane: " << leftLane_closest;
        if (leftLane_closest)
        {
            d1 = -0.2;
        }else{
            d1 = 0.2;
        }
    }

    tot.kj = config().get<double>("kj",1.0);
    tot.kT = config().get<double>("kT",10.0);
    tot.ks = config().get<double>("ks",0.1);
    tot.kd = config().get<double>("kd",10.0);


    double safetyS = config().get<double>("safetyS",0.1); //Sicherheitsabstand tangential zur Straße
    double safetyD = config().get<double>("safetyD",0.1); //Sicherheitsabstand orthogonal zur Straße

    Trajectory result;
    if(generator->createTrajectorySample(result,endVelocity,d1,safetyS,safetyD,tMin,tMax,nSamplesTraj,dataRoad,dataObstacle,tot)){
        //convert data
        //TODO anzahl der punkte
        //points2d<20> points = result.sampleXY<20>();
        //float middleLength = road->length();
        float middleStepLength = road->length()/10;
        lms::math::polyLine2f myroad= road->getWithDistanceBetweenPoints(middleStepLength);
        points2d<10> pointsMiddle;

        for (int i = 0; i < 10; i++)
        {
            pointsMiddle.x(i) = myroad.points()[i].x;
            pointsMiddle.y(i) = myroad.points()[i].y;
        }
        //trajectory = result.projectOntoBezierCurvePlusVelocity<20>(pointsMiddle, 0.1);
        lms::math::polyLine2f lineMiddle = *road;
        trajectory = result.projectOntoLineSegments(lineMiddle, 0.05);
        return true;
    }else{
        return false;
    }
    return true;
}



street_environment::Trajectory TrajectoryLineCreator::simpleTrajectory(float distanceObstacleBeforeChangeLine, float endVelocity){
    //Mindestabstand zwischen zwei Hindernissen 1m
    //Maximalabstand von der Kreuzung: 15cm
    //An der Kreuzung warten: 2s
    const float obstacleTrustThreshold = config().get<float>("obstacleTrustThreshold",0.5);
    const float crossingTrustThreshold = config().get<float>("crossingTrustThreshold",0.5);
    bool useFixedSpeed = false;
    float fixedSpeed= 0;
    bool firstAdd = true;

    street_environment::Trajectory tempTrajectory;
    // translate the middle lane to the right with a quarter of the street width
    const float translation = config().get<float>("street_width", 0.8)/4.0f;
    const float obstacleLength = config().get<float>("obstacleLength",0.5);

    using lms::math::vertex2f;
    if(road->points().size() == 0){
        logger.error("cycle") << "no valid environment given";
        return tempTrajectory;
    }
    //we start in the car - NOT TODAY MY FRIEND
    const float trajectoryStartDistance = config().get<float>("trajectoryStartDistance",0.3);
    const float distanceBetweenTrajectoryPoints = config().get<float>("obstacleResolution",0.05);
    const lms::math::polyLine2f middle = road->getWithDistanceBetweenPoints(distanceBetweenTrajectoryPoints);
    float tangLength = 0;
    //ignore first point
    for(int i = 1; i <(int) middle.points().size(); i++) {
        const vertex2f p1 = middle.points()[i - 1];
        const vertex2f p2 = middle.points()[i];
        if(p1 == p2){ //should never happen
            logger.error("simpleTrajectory")<<"p1 == p2";
            continue;
        }

        const vertex2f along = p2 - p1;
        //check if the trajectory is long enough
        //TODO, get endpoint
        tangLength +=along.length();
        if(tangLength < trajectoryStartDistance /* tangLength > trajectoryMaxLength*/){//TODO trajectoryMaxLength
            continue;
        }
        const vertex2f normAlong = along / along.length();
        const vertex2f orthogonal(normAlong.y, -normAlong.x);
        const vertex2f orthogonalTrans = orthogonal*translation;

        bool rightSide = true;
        //check all obstacles
        LaneState rightState = getLaneState(tangLength,true);
        LaneState leftState = getLaneState(tangLength,false);
        if(rightState > leftState){
            rightSide = false;
        }

        if(rightSide && rightState == LaneState::BLOCKED){
            useFixedSpeed = true;
            fixedSpeed = 0;
        }

        if(useFixedSpeed){
            endVelocity = fixedSpeed;
        }
        vertex2f result;
        if(rightSide){
            if(firstAdd){
                firstAdd = false;
                tempTrajectory.push_back(street_environment::TrajectoryPoint(lms::math::vertex2f(0,0),normAlong,endVelocity,-0.2)); //TODO
            }
            result= p1 + orthogonalTrans;
            tempTrajectory.push_back(street_environment::TrajectoryPoint(result,normAlong,endVelocity,-0.2)); //TODO
        }else{
            if(firstAdd){
                firstAdd = false;
                tempTrajectory.push_back(street_environment::TrajectoryPoint(lms::math::vertex2f(0,0),normAlong,endVelocity,0.2)); //TODO
            }
            result= p1 - orthogonalTrans;
            tempTrajectory.push_back(street_environment::TrajectoryPoint(result,normAlong,endVelocity,0.2)); //TODO
        }
    }

    return tempTrajectory;

}


LaneState TrajectoryLineCreator::getLaneState(float tangDistance, bool rightSide){
    LaneState result = LaneState::CLEAR;
    const float obstacleTrustThreshold = config().get<float>("obstacleTrustThreshold",0.5);
    const float obstacleLength = config().get<float>("obstacleLength",0.5);
    const float obstacleSavetyDistance = config().get<float>("distanceObstacleBeforeChangeLine",0);
    const float crossingTrustThreshold = config().get<float>("crossingTrustThreshold",0.5);
    for(street_environment::EnvironmentObjectPtr objPtr: envObstacles->objects){
        if(objPtr->getType() == street_environment::Obstacle::TYPE){
            if(objPtr->trust() < obstacleTrustThreshold){
                continue;
            }
            street_environment::ObstaclePtr obstPtr = std::static_pointer_cast<street_environment::Obstacle>(objPtr);
            if(!((rightSide && obstPtr->distanceOrth()< 0)||(rightSide && obstPtr->distanceOrth() > 0))){
                continue;
            }
            float distanceToObstacle = tangDistance-obstPtr->distanceTang();//abstand zum Punkt p2
            if(obstPtr->trust() < obstacleTrustThreshold){
                continue;
            }
            if(distanceToObstacle < obstacleLength){
                if(distanceToObstacle >= 0 && distanceToObstacle){
                    result = LaneState::BLOCKED;
                    break;
                }else if(distanceToObstacle >= -obstacleSavetyDistance){
                    result = LaneState::DANGEROUS;
                }
            }
        }else if(objPtr->getType() == street_environment::Crossing::TYPE){
            //Check the trust
            if(objPtr->trust() < crossingTrustThreshold){
                continue;
            }
            const street_environment::CrossingPtr crossing = std::static_pointer_cast<street_environment::Crossing>(objPtr);
            if(crossing->position().x < config().get<float>("crossingMinDistance",0.3)){ //we already missed the trajectory!
                continue;
            }
            //TODO doesn't work that good if(crossing->foundOppositeStopLine || !config().get<bool>("crossingUseOppositeLine",false)){
                if(car->velocity() < 0.1){//TODO HACK but may work
                    if(const_cast<street_environment::Crossing*>(crossing.get())->startStop()){//TODO HACK
                        logger.info("start waiting in front of crossing");
                    }
                }
                logger.info("simpleTrajectory")<<"crossing: stop "<< crossing->hasToStop() << " blocked:"<<crossing->blocked()<< " waiting for:"<<crossing->stopTime().since().toFloat();

                //check if we have to stop or if crossing is blocked
                if(crossing->hasToStop() || crossing->blocked()){
                    //Check if we are waiting for to long
                    if(!crossing->hasToStop() && crossing->stopTime().since().toFloat()>config().get<float>("maxStopTimeAtCrossing",10)){
                        logger.warn("ignoring crossing")<<"I was waiting for "<<crossing->stopTime().since()<<"s";
                    }else{
                        //check if the Crossing is close enough
                        //As there won't be an obstacle in front of the crossing we can go on the right
                        //TODO we won't indicate if we change line
                        if(crossing->distanceTang()-tangDistance < config().get<float>("minDistanceToCrossing",0.1)){
                            result = LaneState::BLOCKED;
                            continue;
                        }
                    }
                }
        }
    }
    return result;
}



