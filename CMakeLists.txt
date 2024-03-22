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


set(Python3_EXECUTABLE /usr/bin/python3.10)

# Find Python3 interpreter and development components
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

add_executable(mpc src/MPC.cpp)
target_link_libraries(mpc ${OpenCV_LIBS} ipopt)
# Link the Python3 interpreter and development components to the target
target_link_libraries(mpc 
    Python3::Python
    Python3::Module
)

target_link_libraries(mpc 
        Python3::NumPy
    )
