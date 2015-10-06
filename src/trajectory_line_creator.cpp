#include "trajectory_line_creator.h"
#include "lms/datamanager.h"
#include "lms/math/math.h"
#include "street_environment/obstacle.h"

extern "C"{
#include "lib/otg_xy_reallyDumb/otg_xy_reallyDumb.h"
#include "lib/otg_xy_reallyDumb/otg_xy_reallyDumb_emxAPI.h"
}

bool TrajectoryLineCreator::initialize() {
    envObstacles = datamanager()->readChannel<street_environment::EnvironmentObjects>(this,"ENVIRONMENT_OBSTACLE");
    road = datamanager()->readChannel<street_environment::RoadLane>(this,"ROAD");
    trajectory = datamanager()->writeChannel<lms::math::polyLine2f>(this,"LINE");
    //TODO für was hat der das car?
    car = datamanager()->writeChannel<sensor_utils::Car>(this,"CAR");
    config = getConfig();

    return true;
}

bool TrajectoryLineCreator::deinitialize() {
    return true;
}
bool TrajectoryLineCreator::cycle() {

    //clear old trajectory
    trajectory->points().clear();
    if(!advancedTrajectory()){
        simpleTrajectory();
    }
    return true;
}

bool TrajectoryLineCreator::advancedTrajectory(){

    //INPUT
    //hindernisse: vector Abstand-Straße,geschwindigkeit-Straße(absolut), -1 (links) +1 (rechts) (alle hindernisse hintereinander)
    //eigenes auto, vx,vy, dw -winkelgeschwindigkeit dw (zunächst mal 0)
    //vector mit x-koordinaten
    //vector mit y-koordinaten
    double S[3];
    S[0] = 0;//Anfangsgeschwindigkeit
    S[1] = 0;//Anfangsbeschleunigung
    S[2] = 1;//Endgeschwindigkeit

    //Winkelgeschwindigkeiten
    double D[4];
    D[0] = -0.2;//
    D[1] = 0;//
    D[2] = 0;//
    D[3] = -0.2;//

    double kj = config->get<double>("kj",1.0);
    double kT = config->get<double>("kT",1.0);
    double ks = config->get<double>("ks",1.0);
    double kd = config->get<double>("kd",1.0);

    double dT = config->get<double>("dT",0.1); //Intervall zwischen den Endzeiten
    //TODO berechnen
    double tMin = 0.1;//Minimal benötigte Zeit
    double tMax = 10; //Maximal benötigte Zeit

    double safetyS = config->get<double>("safetyS",0.2); //Sicherheitsabstand tangential zur Straße
    double safetyD = config->get<double>("safetyS",0.2); //Sicherheitsabstand orthogonal zur Straße

    double dt = config->get<double>("dt",0.01); //Zeitintervall zwischen den Kollisionsabfragen

    double ds = config->get<double>("ds",0.005); //Abstand des Streckenzugs
    double m = config->get<double>("m",20); //Anzahl der Punkte im Streckenzug

    double y0 = road->polarDarstellung[0];
    double phi = road->polarDarstellung[1];


    int obstacle_count = envObstacles->objects.size();

    emxArray_real_T *dataVeh =
            dataVeh = emxCreate_real_T(3,obstacle_count);
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
        dataVeh->data[i*3+ 2] = 1; //ob linke (-1) oder rechte spur (1)
    }

    //get kappa from circle
    static float kappa_old =0;
    float kappa = 0;
    //Wie stark der alte radius ins gewicht fallen soll
    float kappa_ratio = 0.5;

    int kappaCount = 3;
    for(int i = 0; i < kappaCount; i++){
        kappa += road->polarDarstellung[2+i];
    }
    kappa /= kappaCount;
    if(kappa_old == 0){
        kappa_old = kappa;
    }
    kappa = kappa_ratio*kappa_old+(1-kappa_ratio)*kappa;
    kappa_old = kappa;

    logger.debug("advancedTrajectory")<<"kappa: "<<kappa;
    //output
    double flag1;
    double flag2;

    //punkte
    emxArray_real_T *x = emxCreate_real_T(1,m);
    emxArray_real_T *y = emxCreate_real_T(1,m);


    otg_xy_reallyDumb(S,D, kj,
                      kT,  ks,  kd,  dT,  tMin,  tMax,
                      dataVeh,  safetyS, safetyD, dt, ds,
                      m,  kappa,  y0,  phi,  &flag1,  &flag2,
                      x, y);


    logger.debug("advancedTrajectory")<<"flag1: "<<flag1 << " flag2: "<<flag2;
    if(flag1 < 0 || flag2 < 0)
        return false;

    for(int i = 0; i < m; i++){
        lms::math::vertex2f result(x->data[i],y->data[i]);
        trajectory->points().push_back(result);
    }
    return true;

    //OUTPUT
    //gibt x-y koodinaten zurück
}

void TrajectoryLineCreator::simpleTrajectory(){


    // translate the middle lane to the right with a quarter of the street width
    const float translation = config->get<float>("street.width", 0.8)/4.0f;
    //TODO das sollte von der aktuellen geschwindigkeit abhängen!
    float distanceObstacleBeforeChangeLine = 0.4;

    using lms::math::vertex2f;
    if(road->points().size() == 0){
        logger.error("cycle") << "no valid environment given";
        return;
    }

    const street_environment::RoadLane &middle = *road;
    logger.debug("simpleTrajectory")<<"number of obstacles: "<<envObstacles->objects.size();
    for(size_t i = 1; i < middle.points().size(); i++) {
        vertex2f p1 = middle.points()[i - 1];
        vertex2f p2 = middle.points()[i];
        if(p1 == p2)
            continue;

        vertex2f along = p2 - p1;
        vertex2f mid((p1.x + p2.x) / 2., (p1.y + p2.y) / 2.);
        vertex2f normAlong = along / along.length();
        vertex2f orthogonal(normAlong.y, -normAlong.x);

        bool left = false;
        //man geht davon aus, dass die Abstand, in dem man ausweicht deutlich größer ist als das hinderniss lang!
        float obstacleLength = 0.3;
        //check all obstacles
        for(const std::shared_ptr<street_environment::EnvironmentObject> obj : envObstacles->objects){
            if(obj->name().find("OBSTACLE") == std::string::npos){
                logger.warn("cycle")<<"invalid obstacle-type given: "<<obj->name();
                continue;
            }
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
        }

        if(left){
            orthogonal *= -1;
        }
        orthogonal = orthogonal * translation;
        vertex2f result = mid + orthogonal;
        trajectory->points().push_back(result);
    }

    trajectory->reduce([](const lms::math::vertex2f& p1){
        return p1.x < 0;
    });


}

