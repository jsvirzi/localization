cmake_minimum_required(VERSION 3.5)
project(demo)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
include_directories("../include")
include_directories("../lidar")
include_directories("../dataAcquisition")
include_directories("../utils")

set(SOURCE_FILES main.cpp ../lidar/lidar.cpp ../dataAcquisition/dataAcquisition.cpp ../utils/utils.cpp)
add_executable(demo ${SOURCE_FILES})
