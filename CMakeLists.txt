cmake_minimum_required(VERSION 3.1)
project(MPC)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_BUILD_TYPE Debug)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#############################
#Dependencies
#############################

find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(./include)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(mpc src/MPC.cpp)
target_link_libraries(mpc ${OpenCV_LIBS} ipopt)
