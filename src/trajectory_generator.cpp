#include "trajectory_generator.h"

trajectory_generator::trajectory_generator(lms::logging::Logger& _logger) : logger(_logger)
{
    Vector<3> a;
    a << 1, 2, 3;

    //logger.warn("construct") << a;
    std::cout << a << "/ " << a.transpose() << std::endl;

}
