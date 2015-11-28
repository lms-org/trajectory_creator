//
// Created by Lukas Koestler on 23.10.15.
//

#ifndef PROJECT_TYPES_H
#define PROJECT_TYPES_H

#include <Eigen/Dense>

// TYPEDEF
typedef float T;

template<int rows, int cols>
using Matrix = Eigen::Matrix<T, rows, cols>;

template<int rows>
using Vector = Matrix<rows, 1>;

typedef struct
{
    // captures the current state of the road and the car on it
    T y0; // intercept of the center line and the y-axis of the car coordinates (different sign than street coordinates)
    T phi; //angle of the centerline relative to the car x-axis
    T w; // angular velocity of the car (z axis positive rotation)
    T vx0; //initial velocity of the car in x (car) direction
    T ax0; //initial acceleration of the car in x (car) direction
    T kappa; //curvature of the road (approximated as a circle)

} roadData;

typedef struct
{
    // caprutres the data of one obstacle on the road
    T s0; // initial distance along the road of the obstacle to the car
    T v0; // initial velocity of the obstacle (NOT REALTIVE TO THE CAR) along the road
    bool leftLane; //lane of the obstacle
} obstacleData;

typedef struct
{
    //coefficients for the cost function
    T kj; //cost of the non-smoothness (kj high --> gives smooth traj.)
    T kT; //cost of the time the traj. needs (kT high --> gives fast (non-smooth) tral.)
    T kd; //cost. of the lateral direction
    T ks; //cost. of the longitudinal direction
} coeffCtot;

typedef struct
{
    // captures the initial conditions for the polynomial in s (longitudinal) Frenet space
    T v0; //initial velocity in s-direction (!= vx0)
    T a0; //initial acceleration in s-direction (!= ax0)
    T v1; //velocity in s direction at the end of the traj.
} S_initialCond;

typedef struct
{
    // captures the initial conditions for the polynomial in d (lateral) Frenet space
    T d0; //initial offset in d-direction (!= y0, and sign is different from y0)
    T d0d; //initial velocity in d-direction
    T d0dd; //initial acceleration in d-direction
    T d1; //offset in d-direction at the end of the trj.
} D_initialCond;

template<size_t N>
struct  points2d
{
    Vector<N> x;
    Vector<N> y;
};

#endif //PROJECT_TYPES_H
