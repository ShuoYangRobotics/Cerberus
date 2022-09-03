# Introduction
This ROS package contains a visual-inertial-leg odometry (VILO) for Unitree A1 and Go1 robot. Our goal is to provide a compact and low-cost long term position sensing suite for legged robots (A sensing solution only has one IMU, one stereo camera, and leg sensors. Total cost <$1000).

The focus of this work is adding body velocity calculated from leg joint sensors and calibrate potential kinematic parameter errors to improve its accuracy when used in VILO. This odometry uses the optimization framework from the VINS-Fusion, one of the most popular visual inertial odometry. Additional to VINS's image and IMU measurement model, we add a special contact preintegration term. It achieves lower than 1\% position estimation drift on various datasets. More details of the theoratical contribution can be found in our recent paper:

Here are two videos comparing the performance of the VILO and the VINS-Fusion:
[![video1](https://img.youtube.com/vi/jq9ijL9z3RI/0.jpg)](https://www.youtube.com/watch?v=jq9ijL9z3RI)

[![video2](https://img.youtube.com/vi/oNwwQ0l-O4U/0.jpg)](https://www.youtube.com/watch?v=oNwwQ0l-O4U)



## Installation
use Docker and VSCode ''Remote - Containers''. A dockerfile that configures an individual development environment is shown in .devcontainer/Dockerfile.


## Dataset
A Google drive folder https://drive.google.com/drive/folders/13GsFDaBkDrslOl9BfE4AJnOn3ECDXVnc

contains several dataset to test the VILO. 

1. outdoor_snow.bag. The bag contains sensor data collected during the snow walking run shown in the second video above. 

2. indoor_with_ground_truth_1.bag. The robot moves forward and back quickly. Groundtruh data is 

## Test On Dataset
Since the Remote container 

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

From the info list, the most important topics are 

1. **/camera_forward/infra1/image_rect_raw**   -  sensor_msgs/Image - 15 Hz

   **/camera_forward/infra2/image_rect_raw**   -  sensor_msgs/Image - 15 Hz

   They comes from a realsense camera. It is very hard to get the image right because infra cameras will not be useful if infra light emitter is not turned off. In A1_launch there is a launch file for reading realsense camera and handle the emitter setting (a1_launch/launch/front_camera_read.launch).

2. **/hardware_a1/imu**  -  sensor_msgs/Imu - 500Hz

   This is the output from A1 robot. The topic comes from Unitree's API. In Unitree API a topic with type "LowState" contains all the robot data. But we separate IMU data and leg data to make it generic. Because other robots may not be archtected in this way 

3. **/imu/data**   -  sensor_msgs/Imu - 200Hz

   This IMU data comes from an expensive IMU called XSens Mti-200, which is not directly used in our algorithm. It is just used for reference to benchmark the orientation measurement of A1 IMU. It turns out the A1 IMU performs very well comparing to XSens Mti-200, which is a $3000 sensor. So you can be convinced that A1's performance is good

4. **/hardware_a1/joint_foot**  - sensor_msgs/JointState - 500Hz
 

   This is another part of the "LowState". The structure of this JointState can be seen from its name list: There are 16 JointStates in this topic. The first 12 contains motor joint information, and the last four contains foot contact sensor reading. There is not a standard ROS topic to represent all of them so I created this topic. 

   One caviet when users create the same topic is, the joint velocities directly generated by A1 robot is very noisy. My controller differentiates joint angles to generate joint velocities. 


