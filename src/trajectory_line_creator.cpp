#include "trajectory_line_creator.h"
#include "lms/math/math.h"
#include "street_environment/obstacle.h"
#include "street_environment/crossing.h"
bool TrajectoryLineCreator::initialize() {
    envObstacles = readChannel<street_environment::EnvironmentObjects>("ENVIRONMENT_OBSTACLE");
    road = readChannel<street_environment::RoadLane>("ROAD");
    trajectory = writeChannel<street_environment::Trajectory>("LINE");
    car = readChannel<sensor_utils::Car>("CAR");

    generator = new TrajectoryGenerator(logger);

    return true;
}

bool TrajectoryLineCreator::deinitialize() {
    delete generator;
    return true;
}
bool TrajectoryLineCreator::cycle() {
    //clear old trajectory
    trajectory->points().clear();
    //calculate data for creating the trajectory
    float trajectoryMaxLength = config().get<float>("trajectoryMaxLength",2);
    float obstacleTrustThreshold = config().get<float>("obstacleTrustThreshold",0.5);
    //TODO not smart
    street_environment::Trajectory traj;
    if(config().get<bool>("simpleTraj",true)){
        traj= simpleTrajectory(trajectoryMaxLength,obstacleTrustThreshold);
    }else{
        if(!advancedTrajectory(traj)){
            traj = simpleTrajectory(trajectoryMaxLength,obstacleTrustThreshold);
        }
    }
    *trajectory = traj;
    /*trajectory->points().clear();
    for(lms::math::vertex2f &v:traj.points()){
        trajectory->points().push_back(v);
    }
    */
    return true;
}

bool TrajectoryLineCreator::advancedTrajectory(lms::math::polyLine2f &trajectory){
    if(road->polarDarstellung.size() < 8){
        return false;
    }
    //TODO Blinker setzen

    //INPUT
    //hindernisse: vector Abstand-Straße,geschwindigkeit-Straße(absolut), -1 (links) +1 (rechts) (alle hindernisse hintereinander)
    //eigenes auto, vx,vy, dw -winkelgeschwindigkeit dw (zunächst mal 0)
    //vector mit x-koordinaten
    //vector mit y-koordinaten

    int nSamplesTraj = config().get<int>("nSamplesTraj",50);
    double v1 = 2;//endgeschwindigkeit
    double d1 = -0.2;//TODO berechnenAbweichung, die man am Ende haben will
    //aktuelle daten der fahrspur und des autos relativ dazu
    RoadData dataRoad;
    dataRoad.ax0 = 0; //beschl. am anfang

    dataRoad.kappa = generator->circleCurvature(road->points()[1],road->points()[3],road->points()[5]);//->polarDarstellung[4]+road->polarDarstellung[7]+road->polarDarstellung[9];//TODO
    //dataRoad.kappa /= 3.0;
    logger.debug("kappa")<<dataRoad.kappa;
    dataRoad.phi = road->polarDarstellung[1];
    float velocity = 0.001;//Sollte nicht 0 sein, wegen smoothem start
    if(car->velocity() != 0)
        velocity = car->velocity();
    dataRoad.vx0 = velocity;
    dataRoad.w = 0; //aktuelle winkelgeschwindigkeit
    dataRoad.y0 = road->polarDarstellung[0]; //TODO
    std::vector<ObstacleData> dataObstacle;
    CoeffCtot tot;

    tot.kj = config().get<double>("kj",1.0);
    tot.kT = config().get<double>("kT",1.0);
    tot.ks = config().get<double>("ks",1.0);
    tot.kd = config().get<double>("kd",1.0);

    double tMin = 0.1;//Minimal benötigte Zeit
    double tMax = 2; //Maximal benötigte Zeit

    double safetyS = config().get<double>("safetyS",0.1); //Sicherheitsabstand tangential zur Straße
    double safetyD = config().get<double>("safetyD",0.1); //Sicherheitsabstand orthogonal zur Straße

    Trajectory result;
    if(generator->createTrajectorySample(result,v1,d1,safetyS,safetyD,tMin,tMax,nSamplesTraj,dataRoad,dataObstacle,tot)){
        //convert data
        //TODO anzahl der punkte
        points2d<20> points = result.sampleXY<20>();
        for(int i = 0; i < 20; i++){
            trajectory.points().push_back(lms::math::vertex2f(points.x(i),points.y(i)));
        }
        return true;
    }else{
        return false;
    }
    return true;
}

street_environment::Trajectory TrajectoryLineCreator::simpleTrajectory(float trajectoryMaxLength,const float obstacleTrustThreshold){
    //Mindestabstand zwischen zwei Hindernissen 1m
    //Maximalabstand von der Kreuzung: 15cm
    //An der Kreuzung warten: 2s

    street_environment::Trajectory tempTrajectory;
    // translate the middle lane to the right with a quarter of the street width
    const float translation = config().get<float>("street_width", 0.8)/4.0f;
    const float obstacleLength = config().get<float>("obstacleLength",0.5);
    //Könnte von der aktuellen geschwindigkeit abhängen
    float distanceObstacleBeforeChangeLine = config().get<float>("distanceObstacleBeforeChangeLine",0.6);

    using lms::math::vertex2f;
    if(road->points().size() == 0){
        logger.error("cycle") << "no valid environment given";
        return tempTrajectory;
    }
    //we start in the car
    tempTrajectory.points().push_back(lms::math::vertex2f(0,0));
    tempTrajectory.viewDirs.points().push_back(lms::math::vertex2f(1,0));

    const street_environment::RoadLane &middle = *road;
    //length along the road
    float tangLength = 0;
    bool lastWasLeft = false;
    for(size_t i = 1; i < middle.points().size(); i++) {
        const vertex2f p1 = middle.points()[i - 1];
        const vertex2f p2 = middle.points()[i];
        if(p1 == p2) //should never happen
            continue;

        vertex2f along = p2 - p1;
        //check if the trajectory is long enough
        //TODO, get endpoint
        tangLength +=along.length();
        if(tangLength > trajectoryMaxLength){
            break;
        }
        vertex2f mid((p1.x + p2.x) / 2., (p1.y + p2.y) / 2.);
        vertex2f normAlong = along / along.length();
        vertex2f orthogonal(normAlong.y, -normAlong.x);

        bool left = false;
        //man geht davon aus, dass die Abstand, in dem man ausweicht deutlich größer ist als das hinderniss lang!
        float distanceToObstacle = 0;
        //check all obstacles
        //TODO not smart at all, won't work in all cases
        for(const std::shared_ptr<street_environment::EnvironmentObject> obj : envObstacles->objects){
            if(obj->getType() == street_environment::Obstacle::TYPE){
               const street_environment::ObstaclePtr obst = std::static_pointer_cast<street_environment::Obstacle>(obj);
                //check if the obstacle is trusted
                if(obst->trust() < obstacleTrustThreshold){
                    continue;
                }

                distanceToObstacle = obst->distanceTang()-tangLength;//abstand zum Punkt p2

                if((distanceToObstacle >-distanceObstacleBeforeChangeLine) && (distanceToObstacle < obstacleLength)){
                    left = true;
                    break;
                }
            }else if(obj->getType() == street_environment::Crossing::TYPE){
                const street_environment::CrossingPtr crossing = std::static_pointer_cast<street_environment::Crossing>(obj);
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
                        float x = crossing->position().x-config().get<float>("minDistanceToCrossing",0.1);//Wir gehen davon aus, dass crossing.distanceTang() == crossing.position.x ist
                        float y= crossing->position().y;
                        //As there won't be an obstacle in front of the crossing we can go on the right
                        //TODO we won't indicate if we change line
                        if(pow(x*x+y*y,0.5)-mid.length() < distanceObstacleBeforeChangeLine ){
                            orthogonal = orthogonal * translation;
                            vertex2f result = mid + orthogonal;
                            tempTrajectory.points().push_back(result);
                            tempTrajectory.viewDirs.points().push_back(normAlong);
                            //add endPoint
                            //TODO wir gehen davon aus, dass die Kreuzung in der Mitte der rechten Linie ihre Position hat!
                            tempTrajectory.points().push_back(lms::math::vertex2f(x,y));
                            tempTrajectory.viewDirs.points().push_back(normAlong);
                            return tempTrajectory;
                        }
                    }
                }
                //add endPoint
            }else{
                //I don't care about astartLine/whatever
                continue;
            }
        }

        orthogonal = orthogonal * translation;
        if(left){
            orthogonal *= -1;
        }
        bool addPoint = true;
        //set line indicator
        if(i != 1){
            if(lastWasLeft != left){
                tempTrajectory.addChange(tempTrajectory.points().size(),left);
                if(left){
                    //wir fahren von rechts nach links
                    float distanceToObstacleChange = distanceToObstacle-distanceObstacleBeforeChangeLine;
                    //if(fabs(distanceToObstacleChange)<along.length()){
                        lms::math::vertex2f temp = p2+normAlong*distanceToObstacleChange;
                        tempTrajectory.points().push_back(temp-orthogonal);//rechter eckpunkt
                        tempTrajectory.points().push_back(temp +orthogonal);//linker eckpunkt
                        tempTrajectory.viewDirs.points().push_back(normAlong);
                        tempTrajectory.viewDirs.points().push_back(normAlong);
                        addPoint = false; //we don't need the other point
                        //TODO left to right lane change
                    //}
                }else{
                    //wir fahren von links nach rechts
                    float distanceToObstacleChange = distanceToObstacle+obstacleLength;
                    lms::math::vertex2f temp = p2+normAlong*distanceToObstacleChange;
                    tempTrajectory.points().push_back(temp+orthogonal);//linker eckpunkt
                    tempTrajectory.points().push_back(temp -orthogonal);//rechter eckpunkt
                    tempTrajectory.viewDirs.points().push_back(normAlong);
                    tempTrajectory.viewDirs.points().push_back(normAlong);
                }
            }
        }
        lastWasLeft=left;
        if(addPoint){
            vertex2f result = p1 + orthogonal;
            tempTrajectory.points().push_back(result);
            tempTrajectory.viewDirs.points().push_back(normAlong);
        }
    }

    /*Doesn't work iwth viewDirs
    //remove invalid points
    tempTrajectory.reduce([](const lms::math::vertex2f& p1){
        return p1.x < 0;
    });
    */
    float maxAngle = config().get<float>("maxAngleBetweenTrajectoryPoints",M_PI/6);
    //float distanceBetweenTrajectoryPoints = config().get<float>("distanceBetweenTrajectoryPoints",0.3);
    int maxPointsRemoved = config().get<int>("maxPointsRemoved",2);
    //TODO doesn't work with viewDirstempTrajectory.points() = tempTrajectory.getWithDistanceBetweenPoints(distanceBetweenTrajectoryPoints).points();
    int currentPointsRemoved = 0;
    //Smooth the trajectory (in some super professional way)
    for(int i = 2; i <(int) tempTrajectory.points().size();){
        lms::math::vertex2f v1 = tempTrajectory.points()[i-1]-tempTrajectory.points()[i-2];
        lms::math::vertex2f v2 = tempTrajectory.points()[i]-tempTrajectory.points()[i-1];
        float angle = v2.angleBetween(v1);
        if(angle > maxAngle && currentPointsRemoved < maxPointsRemoved){
            //remove the middle
            //TODO handle the change!
            tempTrajectory.points().erase(tempTrajectory.points().begin()+i-1);
            tempTrajectory.viewDirs.points().erase(tempTrajectory.viewDirs.points().begin()+i-1);
            currentPointsRemoved++;
        }else{
            currentPointsRemoved = 0;
            i++;
        }
    }
    if(tempTrajectory.points().size()!=tempTrajectory.viewDirs.points().size()){
        logger.error("INVALID viewDir count")<<tempTrajectory.points().size()<<" "<<tempTrajectory.viewDirs.points().size();
    }
    return tempTrajectory;

}

