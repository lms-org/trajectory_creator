#include "trajectory_line_creator.h"
#include "lms/datamanager.h"
#include "lms/math/math.h"
#include "street_environment/obstacle.h"
#include "street_environment/crossing.h"

bool TrajectoryLineCreator::initialize() {
    envObstacles = datamanager()->readChannel<street_environment::EnvironmentObjects>(this,"ENVIRONMENT_OBSTACLE");
    road = datamanager()->readChannel<street_environment::RoadLane>(this,"ROAD");
    trajectory = datamanager()->writeChannel<lms::math::polyLine2f>(this,"LINE");
    config = getConfig();
    kappa_old = 0;

    //   generator = new trajectory_generator(logger);

    return true;
}

bool TrajectoryLineCreator::deinitialize() {

    //    delete generator;

    return true;
}
bool TrajectoryLineCreator::cycle() {
    //clear old trajectory
    trajectory->points().clear();
    //calculate data for creating the trajectory
    float trajectoryMaxLength = 1;
    float endX;
    float endY;
    float endVx = 0;
    float endVy = 0;
    simpleTrajectory(trajectoryMaxLength,endVx,endVy);
    return true;
}

bool TrajectoryLineCreator::advancedTrajectory(){

    //INPUT
    //hindernisse: vector Abstand-Straße,geschwindigkeit-Straße(absolut), -1 (links) +1 (rechts) (alle hindernisse hintereinander)
    //eigenes auto, vx,vy, dw -winkelgeschwindigkeit dw (zunächst mal 0)
    //vector mit x-koordinaten
    //vector mit y-koordinaten
    double v1 = 2;//endgeschwindigkeit
    double d1 = 0.2;//Abweichung, die man am Ende haben will
    double vx0 = 1; //Sollte nicht 0 sein, wegen smoothem start
    double ax0 = 0; //beschl. am anfang
    double w = 0; //aktuelle winkelgeschwindigkeit

    double kj = config->get<double>("kj",1.0);
    double kT = config->get<double>("kT",1.0);
    double ks = config->get<double>("ks",1.0);
    double kd = config->get<double>("kd",1.0);

    double dT = config->get<double>("dT",0.1); //Intervall zwischen den Endzeiten
    //TODO berechnen
    double tMin = 0.1;//Minimal benötigte Zeit
    double tMax = 10; //Maximal benötigte Zeit

    double safetyS = config->get<double>("safetyS",0.1); //Sicherheitsabstand tangential zur Straße
    double safetyD = config->get<double>("safetyD",0.1); //Sicherheitsabstand orthogonal zur Straße

    double dt = config->get<double>("dt",0.01); //Zeitintervall zwischen den Kollisionsabfragen

    double m = config->get<double>("m",20); //Anzahl der Punkte im Streckenzug

    double y0 = road->polarDarstellung[0];
    double phi = road->polarDarstellung[1];


    int obstacle_count = envObstacles->objects.size();

    /*
    emxArray_real_T *dataVeh = emxCreate_real_T(3,obstacle_count);
    for(int i = 0; i < obstacle_count; i++){
        const std::shared_ptr<street_environment::EnvironmentObject> &obj = envObstacles->objects[i];
        if(obj->getType() != 1){
            logger.warn("cycle")<<"invalid obstacle-type given: "<<obj->name();
            continue;
        }
        const street_environment::Obstacle &obst = obj->getAsReference<const street_environment::Obstacle>();
        float x = obst.position().x;
        float y= obst.position().y;

        dataVeh->data[i*3] = x;//abstand zum hindernis
        logger.debug("Abstand zum hinderniss: ")<<x;
        dataVeh->data[i*3 +1] = 0; //geschwindigkeit vom hindernis
        //TODO rechte oder linke spur
        dataVeh->data[i*3+ 2] = -1;  //TODO
    }

    //get kappa from circle
    float kappa = 0;
    //Wie stark der alte radius ins gewicht fallen soll
    float kappa_ratio = config->get<float>("kappa_ratio",0);

    int kappaCount = 3;
    //std::cout << "trajec-creator: kappa-values: ";
    for(int i = 0; i < kappaCount; i++){
        kappa += road->polarDarstellung[5+i];
        //std::cout<< std::to_string(road->polarDarstellung[2+i])<< " , ";
    }
    //std::cout<<std::endl;
    kappa /= kappaCount;
    logger.debug("kappa")<<kappa_old <<" , "<< kappa << " ratio: "<<kappa_ratio;
    kappa = kappa_ratio*kappa_old+(1-kappa_ratio)*kappa;
    kappa_old = kappa;

    logger.debug("advancedTrajectory")<<"kappa: "<<kappa;
    //output
    double flag;
    double T = -1; //gesamtzeit

    //punkte
    emxArray_real_T *x = emxCreate_real_T(1,m);
    emxArray_real_T *y = emxCreate_real_T(1,m);


    otg_xy_reallyDumb(v1,d1, kj,
                      kT,  ks,  kd,  dT,  tMin,  tMax,
                      dataVeh,  safetyS, safetyD, dt, m,  kappa,  y0,  phi,vx0,ax0,w,  &flag,x, y,&T);


    logger.debug("advancedTrajectory")<<"flag: "<<flag;
    if(flag < 0)
        return false;

    for(int i = 0; i < m; i++){
        lms::math::vertex2f result(x->data[i],y->data[i]);
        trajectory->points().push_back(result);
    }
     */
    return true;

    //OUTPUT
    //gibt x-y koodinaten zurück
}

lms::math::polyLine2f TrajectoryLineCreator::simpleTrajectory(float trajectoryMaxLength,float &endVx,float &endVy){
    //TODO use trajectoryMaxLength
    lms::math::polyLine2f tempTrajectory;
    // translate the middle lane to the right with a quarter of the street width
    const float translation = config->get<float>("street.width", 0.8)/4.0f;
    //TODO das sollte von der aktuellen geschwindigkeit abhängen!
    float distanceObstacleBeforeChangeLine = 0.4;

    using lms::math::vertex2f;
    if(road->points().size() == 0){
        logger.error("cycle") << "no valid environment given";
        return tempTrajectory;
    }

    const street_environment::RoadLane &middle = *road;
    float currentTrajectoryLength = 0;
    for(size_t i = 1; i < middle.points().size(); i++) {
        vertex2f p1 = middle.points()[i - 1];
        vertex2f p2 = middle.points()[i];
        if(p1 == p2)
            continue;

        vertex2f along = p2 - p1;
        //check if the trajectory is long enough
        //TODO, get endpoint
        currentTrajectoryLength +=along.length();
        if(currentTrajectoryLength > trajectoryMaxLength){
            break;
        }
        vertex2f mid((p1.x + p2.x) / 2., (p1.y + p2.y) / 2.);
        vertex2f normAlong = along / along.length();
        vertex2f orthogonal(normAlong.y, -normAlong.x);

        bool left = false;
        //man geht davon aus, dass die Abstand, in dem man ausweicht deutlich größer ist als das hinderniss lang!
        float obstacleLength = 0.3;
        //check all obstacles
        for(const std::shared_ptr<street_environment::EnvironmentObject> obj : envObstacles->objects){
            if(obj->getType() == street_environment::Obstacle::TYPE){
                const street_environment::Obstacle &obst = obj->getAsReference<const street_environment::Obstacle>();
                float x = obst.position().x;
                float y= obst.position().y;
                if(x < 0){
                    x += obstacleLength;
                }
                if(pow(x*x+y*y,0.5)-mid.length() < distanceObstacleBeforeChangeLine ){
                    left = true;
                    break;
                }


                if(left){
                    orthogonal *= -1;
                }
                orthogonal = orthogonal * translation;
                vertex2f result = mid + orthogonal;
                tempTrajectory.points().push_back(result);
            }else if(obj->getType() == street_environment::Crossing::TYPE){
                const street_environment::Crossing &crossing = obj->getAsReference<const street_environment::Crossing>();
                //check if the Crossing is close enough
                //TODO
                if(crossing.getStreetDistanceTangential() < trajectoryMaxLength){

                    endVx = 0;
                    endVy = 0;
                }
                //add endPoint
            }else{
                logger.warn("cycle")<<"invalid obstacle-type given: "<<obj->name();
                continue;
            }
        }

        tempTrajectory.reduce([](const lms::math::vertex2f& p1){
            return p1.x < 0;
        });
    }
    return tempTrajectory;

}

