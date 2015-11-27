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
    T y0;
    T phi;
    T w;
    T vx0;
    T ax0;
    T kappa;
} roadData;

typedef struct
{
    T s0;
    T v0;
    bool leftLane;
} obstacleData;

typedef struct
{
    T kj;
    T kT;
    T kd;
    T ks;
} coeffCtot;

typedef struct
{
    T v0;
    T a0;
    T v1;
} S_initialCond;

typedef struct
{
    T d0;
    T d0d;
    T d0dd;
    T d1;
} D_initialCond;

template<size_t N>
struct  points2d
{
    Vector<N> x;
    Vector<N> y;
};

#endif //PROJECT_TYPES_H
