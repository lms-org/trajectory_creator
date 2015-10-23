#ifndef TRAJECTORYGENERATOR_H
#define TRAJECTORYGENERATOR_H

#include <lms/logging/logger.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

class trajectory_generator {
public:

    typedef float T;

    template<int rows, int cols>
    using Matrix = Eigen::Matrix<T, rows, cols>;

    template<int rows>
    using Vector = Matrix<rows, 1>;

public:

    /**
     * @brief Constructor
     */
    trajectory_generator(lms::logging::Logger& _logger);

protected:
    lms::logging::Logger logger;

};


#endif //TRAJECTORYGENERATOR_H
