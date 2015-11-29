//
// Created by Lukas Koestler on 29.11.15.
//

#ifndef PROJECT_BEZIERCURVE_H
#define PROJECT_BEZIERCURVE_H

#include <vector>

#include "BezierPolynomial.h"

template<size_t n>
class BezierCurve {
public:
    BezierCurve(const points2d<n+1> controlPointsIn, const T t_beginIn, const T t_endIn) : t_begin(initBegin(t_beginIn, t_endIn)),
                                                                                                  t_end(initEnd(t_beginIn, t_endIn))
    {

        //init the two polynomials
        poly_x = BezierPolynomial<n>(controlPointsIn.x, t_begin, t_end);
        poly_y = BezierPolynomial<n>(controlPointsIn.y, t_begin, t_end);

        poly_dx = poly_x.differentiate();
        poly_dy = poly_y.differentiate();
    }

    Vector<2> evalAtPoint(const T t)
    {
        Vector<2> result;
        result(0) = 0;
        result(1) = 1;

        if (t < t_begin || t > t_end)
        {
            // throw error
            return result;
        }

        result(0) = poly_x.evalAtPoint(t);
        result(1) = poly_y.evalAtPoint(t);

        return result;

    }

    template <size_t m>
    points2d<m> eval(const Vector<m> tt)
    {
        points2d<m> points;
        points.x = poly_x.eval<m>(tt);
        points.y = poly_y.eval<m>(tt);

        return points;
    }

    Vector<2> normalAtPoint(const T t) {
        Vector<2> result;
        result(0) = 0;
        result(1) = 1;

        if (t < t_begin || t > t_end) {
            // throw error
            return result;
        }

        T dy = poly_dy.evalAtPoint(t);
        T dx = poly_dx.evalAtPoint(t);

        T norm = sqrt(pow(dy, 2) + pow(dx, 2));

        if (norm == 0)
        {
            result(0) = -dy;
            result(1) = dx;
        }else
        {
            result(0) = -dy/norm;
            result(1) = dx/norm;
        }
        return result;
    }

    Vector<2> tangentAtPoint(const T t) {
        Vector<2> result;
        result(0) = 0;
        result(1) = 1;

        if (t < t_begin || t > t_end) {
            // throw error
            return result;
        }

        T dy = poly_dy.evalAtPoint(t);
        T dx = poly_dx.evalAtPoint(t);

        T norm = sqrt(pow(dy, 2) + pow(dx, 2));

        if (norm == 0)
        {
            result(0) = dx;
            result(1) = dy;
        }else
        {
            result(0) = dx/norm;
            result(1) = dy/norm;
        }
        return result;
    }

private:
    BezierPolynomial<n> poly_x;
    BezierPolynomial<n> poly_y;

    BezierPolynomial<n-1> poly_dx;
    BezierPolynomial<n-1> poly_dy;

    const T t_begin;
    const T t_end;


    T initBegin(const T t_beginIn, const T t_endIn)
    {
        if (t_beginIn >= t_endIn)
        {
            return 0;
        }else
        {
            return t_beginIn;
        }
    }

    T initEnd(const T t_beginIn, const T t_endIn)
    {
        if (t_beginIn >= t_endIn)
        {
            return 1;
        }else
        {
            return t_endIn;
        }
    }
};


#endif //PROJECT_BEZIERCURVE_H
