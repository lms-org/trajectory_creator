set ( SOURCES
        "src/trajectory_line_creator.cpp"
	"src/interface.cpp"
)

set (HEADERS
        "include/trajectory_line_creator.h"
)
include_directories("include" "codegen")


file(GLOB matlab "codegen/lib/otg_xy_reallyDumb/*.c" "codegen/lib/otg_xy_reallyDumb/*.h")
add_library (trajectory_creator MODULE ${SOURCES} ${HEADERS} ${matlab})
target_link_libraries(trajectory_creator PRIVATE lmscore math_lib sensor_utils street_environment)
