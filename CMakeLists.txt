cmake_minimum_required (VERSION 3.1)
project (inertial_tools)

if (EXISTS "$ENV{WORKSPACE}/sysconfig/cmake/physycom_config.cmake")
   message(STATUS "Enabling physycom settings")
   include("$ENV{WORKSPACE}/sysconfig/cmake/physycom_config.cmake")
else()
   message(STATUS "Unable to find physycom settings file")
endif()


find_package(Boost REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(LIBBBMUTILS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src)
set(JSONCONS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src/jsoncons/include)

include_directories(${JSONCONS_INCLUDE_DIR} ${LIBBBMUTILS_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
add_library(math_func ${CMAKE_SOURCE_DIR}/src/math_func.cpp)
add_library(bbmutils ${CMAKE_SOURCE_DIR}/src/libbbmutils/src/bbmutils.c)

add_executable(data_calibrator ${CMAKE_SOURCE_DIR}/src/data_calibrator.cpp)
add_executable(data_filter ${CMAKE_SOURCE_DIR}/src/data_filter.cpp)
add_executable(data_gen ${CMAKE_SOURCE_DIR}/src/data_gen.cpp)
add_executable(data_purge ${CMAKE_SOURCE_DIR}/src/data_purge.cpp)
add_executable(data_reco ${CMAKE_SOURCE_DIR}/src/data_reco.cpp)
add_executable(data_rotator ${CMAKE_SOURCE_DIR}/src/data_rotator.cpp)
add_executable(inertial_stats ${CMAKE_SOURCE_DIR}/src/inertial_stats.cpp)
add_executable(mutual_entropy ${CMAKE_SOURCE_DIR}/src/mutual_entropy.cpp)
add_executable(speed_compare ${CMAKE_SOURCE_DIR}/src/speed_compare.cpp)

target_link_libraries(data_calibrator math_func bbmutils)
target_link_libraries(data_filter math_func bbmutils)
target_link_libraries(data_gen math_func bbmutils)
target_link_libraries(data_purge math_func bbmutils)
target_link_libraries(data_reco math_func bbmutils)
target_link_libraries(data_rotator math_func bbmutils)
target_link_libraries(inertial_stats math_func bbmutils)
target_link_libraries(mutual_entropy math_func bbmutils)
target_link_libraries(speed_compare math_func bbmutils)

install (TARGETS
    data_calibrator
    data_filter
    data_gen
    data_purge
    data_reco
    data_rotator
    inertial_stats
    mutual_entropy
    speed_compare
    DESTINATION
    ${CMAKE_SOURCE_DIR}/bin/
)
