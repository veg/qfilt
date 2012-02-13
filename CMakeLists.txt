cmake_minimum_required(VERSION 2.8)
project(454filter)

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
    454filter
    src/454_filter.c
)

target_link_libraries (454filter ${DEFAULT_LIBRARIES})


set(DEFAULT_COMPILE_FLAGS "-O3")
set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")


install(
    TARGETS 454filter
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(454-filter 
	DEPENDS 454filter
)



