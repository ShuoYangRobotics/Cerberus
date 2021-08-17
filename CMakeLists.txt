cmake_minimum_required(VERSION 3.0.0)
project(vileom)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14")
#-DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    geometry_msgs
    nav_msgs
    tf
    cv_bridge
    camera_models
    image_transport

    # must install unitree_ros
    unitree_legged_msgs
)

find_package(OpenCV REQUIRED)
# message(WARNING "OpenCV_VERSION: ${OpenCV_VERSION}")

find_package(Ceres REQUIRED)
#message(WARNING "CERES_LIBRARIES: ${CERES_LIBRARIES}")

include_directories(${catkin_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3)
include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

catkin_package()

add_library(vileom_lib
    src/estimator/estimator.cpp
    src/featureTracker/feature_tracker.cpp
    src/featureTracker/feature_manager.cpp

    src/utils/parameters.cpp
    src/utils/utility.cpp
    src/utils/visualization.cpp
    src/utils/CameraPoseVisualization.cpp

    src/factor/pose_local_parameterization.cpp
    src/factor/projectionTwoFrameOneCamFactor.cpp
    src/factor/projectionTwoFrameTwoCamFactor.cpp
    src/factor/projectionOneFrameTwoCamFactor.cpp
    src/factor/marginalization_factor.cpp

    src/legKinematics/A1Kinematics.cpp

    src/initial/solve_5pts.cpp
    src/initial/initial_aligment.cpp
    src/initial/initial_sfm.cpp
    src/initial/initial_ex_rotation.cpp
)
target_link_libraries(vileom_lib ${catkin_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES})


add_executable(vileom_feature_track_test src/feature_track_test.cpp )
target_link_libraries(vileom_feature_track_test vileom_lib)


add_executable(ceres_test src/test/ceres_test.cpp src/legKinematics/A1Kinematics.cpp)
target_link_libraries(ceres_test ${CERES_LIBRARIES})

