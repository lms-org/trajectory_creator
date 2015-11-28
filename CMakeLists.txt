set ( SOURCES
        "src/trajectory_line_creator.cpp"
        "src/trajectory_generator.cpp"
        "src/interface.cpp"
)

set (HEADERS
        "include/trajectory_line_creator.h"
         include/trajectory_generator.h include/poly.h include/types.h src/trajectory.cpp include/trajectory.h)



include_directories("include" /usr/local/opt/eigen/include/eigen3)

add_library (trajectory_creator MODULE ${SOURCES} ${HEADERS})
target_link_libraries(trajectory_creator PRIVATE lmscore math_lib sensor_utils street_environment)
