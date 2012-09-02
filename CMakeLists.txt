cmake_minimum_required(VERSION 2.8)
project(qfilt)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(DEFAULT_LIBRARIES m)

set(CMAKE_CONFIGURATION_TYPES Release)

include_directories(
    src/
)

add_executable(
    qfilt
    src/main.cpp
    src/argparse.cpp
    src/seq.cpp
)

target_link_libraries (qfilt ${DEFAULT_LIBRARIES})

# do not remove -Wall and -Werror -- fix the errors instead of being lazy
set(DEFAULT_COMPILE_FLAGS "-O0 -g -Wall -Werror")
set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")

set_target_properties(
    qfilt
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}"
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}"
)

install(
    TARGETS qfilt
	RUNTIME DESTINATION bin
	OPTIONAL
)
