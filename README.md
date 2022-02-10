This ROS package contains a visual-inertial-leg odometry (VILO) for Unitree A1 robot. The package name is currently called vileom, which is an obsolete name, which will be revised in later versions. 

Our goal is to achieve low-cost minimum sensing on legged robots (sensing suites only has one IMU, one stereo camera, leg joint sensors. Total cost <$1000).

This odometry is based on the VINS-Fusion, one of the most popular visual inertial odometry.

Here are two videos comparing the performance of VILO and VINS
[![video1](https://img.youtube.com/vi/jq9ijL9z3RI/0.jpg)](https://www.youtube.com/watch?v=jq9ijL9z3RI)

[![video2](https://img.youtube.com/vi/oNwwQ0l-O4U/0.jpg)](https://www.youtube.com/watch?v=oNwwQ0l-O4U)


The package is only tested in Ubuntu 18.04 and ROS melodic. 

## Installation
Please follow VINS-Fusion's installation guide. 


This package must also working together with 

https://github.com/paulyang1990/A1_visualizer

and

https://github.com/paulyang1990/a1_launch

The A1_visualizer has its own dependencies. A1_launch is just a collection of launch files. Download these two packages and compile them together with the VILO (vileom).

The A1_launch contains scripts to run the VILO. And the VILO uses A1_visualizer to visualize the estimation results. 

## Dataset
A Google drive folder https://drive.google.com/drive/folders/13GsFDaBkDrslOl9BfE4AJnOn3ECDXVnc

contains several dataset to test the VILO. 

1. outdoor_snow.bag. The bag contains sensor data collected during the snow walking run shown in the second video above. 

2. indoor_with_ground_truth_1.bag. The robot moves forward and back quickly. Groundtruh data is 

## Test On Dataset

```
roslaunch a1_launch hardware_a1_vilo_together.launch
# wait for a few seconds
rosbag play --clock -r 0.5 outdoor_snow.bag
```
Notice the rosbag play should be slow down for slow computers, otherwise the VILO cannot finish computation in time. 
## Sensor Setup
The VILO can only works properly when sensor topics are received correctly and all sensor transformations are set properly. 

We use the "outdoor_snow.bag" dataset to overview the sensor data frequency and format. Command "rosbag info outdoor_snow.bag" shows the following information
```
version:     2.0
duration:    1:41s (101s)
start:       Jan 18 2022 15:27:10.62 (1642537630.62)
end:         Jan 18 2022 15:28:52.55 (1642537732.55)
size:        949.0 MB
messages:    120363
compression: none [1020/1020 chunks]
types:       sensor_msgs/Image      [060021388200f6f0f447d0fcd9c64743]
             sensor_msgs/Imu        [6a62c6daae103f4ff57a132d6f95cec2]
             sensor_msgs/JointState [3066dcd76a6cfaef579bd0f34173e9fd]
topics:      /camera_forward/infra1/image_rect_raw    1529 msgs    : sensor_msgs/Image     
             /camera_forward/infra2/image_rect_raw    1529 msgs    : sensor_msgs/Image     
             /hardware_a1/imu                        48466 msgs    : sensor_msgs/Imu       
             /hardware_a1/joint_foot                 48452 msgs    : sensor_msgs/JointState
             /imu/data                               20387 msgs    : sensor_msgs/Imu

```
For dataset with ground truth position, an additonal topic "/mocap_node/mocap/pose" is in the bag. The data type is "geometry_msgs/PoseStamped".
