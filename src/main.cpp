// standard C
#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>

// opencv and eigen, eigen must go first
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

// ros related
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/exact_time.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
#include <sensor_msgs/JointState.h>

// project related
#include "estimator/estimator.h"
#include "utils/parameters.h"
#include "utils/visualization.h"

#include "kalmanFilter/A1KFCombineLOWithFoot.h"

// global variables
Estimator estimator;

// kf for estimating contact 
// sensor process to filter imu and leg data
A1SensorData data;
A1KFCombineLOWithFoot kf;  // Kalman filter Baseline 3 with foot
double curr_t;

// debug print filtered data
ros::Publisher filterd_imu_pub;
ros::Publisher filterd_joint_pub;
ros::Publisher filterd_pos_pub;

bool first_sensor_received = false;

queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::ImageConstPtr> img0_buf;
queue<sensor_msgs::ImageConstPtr> img1_buf;
std::mutex m_buf;


void img0_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();
    img0_buf.push(img_msg);
    m_buf.unlock();
}

void img1_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();
    img1_buf.push(img_msg);
    m_buf.unlock();
}


cv::Mat getImageFromMsg(const sensor_msgs::ImageConstPtr &img_msg)
{
    cv_bridge::CvImageConstPtr ptr;
    if (img_msg->encoding == "8UC1")
    {
        sensor_msgs::Image img;
        img.header = img_msg->header;
        img.height = img_msg->height;
        img.width = img_msg->width;
        img.is_bigendian = img_msg->is_bigendian;
        img.step = img_msg->step;
        img.data = img_msg->data;
        img.encoding = "mono8";
        ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
    }
    else
        ptr = cv_bridge::toCvCopy(img_msg, sensor_msgs::image_encodings::MONO8);

    cv::Mat img = ptr->image.clone();
    return img;
}

// extract images with same timestamp from two topics
int record_counter = 0;
int record_counter_interval = 50;
void sync_process()
{
    while(1)
    {
        if(STEREO)
        {
            cv::Mat image0, image1;
            std_msgs::Header header;
            double time = 0;
            m_buf.lock();
            if (!img0_buf.empty() && !img1_buf.empty())
            {
                double time0 = img0_buf.front()->header.stamp.toSec();
                double time1 = img1_buf.front()->header.stamp.toSec();
                // 0.003s sync tolerance
                if(time0 < time1 - 0.003)
                {
                    img0_buf.pop();
                    printf("throw img0\n");
                }
                else if(time0 > time1 + 0.003)
                {
                    img1_buf.pop();
                    printf("throw img1\n");
                }
                else
                {
                    time = img0_buf.front()->header.stamp.toSec();
                    header = img0_buf.front()->header;
                    image0 = getImageFromMsg(img0_buf.front());
                    img0_buf.pop();
                    image1 = getImageFromMsg(img1_buf.front());
                    img1_buf.pop();
                    //printf("find img0 and img1\n");
                }
            }
            m_buf.unlock();
            if(!image0.empty())
                estimator.inputImage(time, image0, image1);
        }
        else
        {
            cv::Mat image;
            std_msgs::Header header;
            double time = 0;
            m_buf.lock();
            if(!img0_buf.empty())
            {
                time = img0_buf.front()->header.stamp.toSec();
                header = img0_buf.front()->header;
                image = getImageFromMsg(img0_buf.front());
                img0_buf.pop();
            }
            m_buf.unlock();
            if(!image.empty())
                estimator.inputImage(time, image);
        }
        record_counter++;
        if (record_counter % record_counter_interval == 0) 
        {
            // write result to file
            ofstream foutC(VILO_RESULT_PATH, ios::app);
            foutC.setf(ios::fixed, ios::floatfield);
            foutC.precision(0);
            double time = ros::Time::now().toSec();
            foutC << time * 1e9 << ",";           // 1
            foutC.precision(5);

            // vilo: convert IMU position to robot body position
            Quaterniond tmp_Q;
            tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
            Eigen::Vector3d p_wb(estimator.Ps[WINDOW_SIZE].x(), estimator.Ps[WINDOW_SIZE].y(), estimator.Ps[WINDOW_SIZE].z());
            Eigen::Vector3d v_wb(estimator.Vs[WINDOW_SIZE].x(), estimator.Vs[WINDOW_SIZE].y(), estimator.Vs[WINDOW_SIZE].z());
            Eigen::Vector3d omega = estimator.latest_gyr_0 - estimator.latest_Bg;
            Eigen::Vector3d p_wr = p_wb + tmp_Q.toRotationMatrix()*estimator.R_br*estimator.p_br;

            Eigen::Vector3d v_wr = v_wb + tmp_Q.toRotationMatrix()*Utility::skewSymmetric(omega)*estimator.R_br*estimator.p_br;

            // kf 
            Eigen::Matrix<double, EKF_STATE_SIZE,1> kf_state = kf.get_state();
            foutC
                << p_wr.x() << ","                                // 2
                << p_wr.y() << ","                                // 3
                << p_wr.z() << ","                                // 4
                << v_wr.x() << ","                                // 5
                << v_wr.y() << ","                                // 6
                << v_wr.z() << ","                                // 7
                << kf_state[0] << ","                                // 8
                << kf_state[1] << ","                                // 9
                << kf_state[2] << ","                                // 10
                << kf_state[3] << ","                                // 11
                << kf_state[4] << ","                                // 12
                << kf_state[5] << ","                                // 13
                << data.opti_pos[0] << ","                                // 14
                << data.opti_pos[1] << ","                                // 15
                << data.opti_pos[2] << ","                                // 16
                << endl;
            foutC.close();
        }


        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }
}

void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    for (unsigned int i = 0; i < feature_msg->points.size(); i++)
    {
        int feature_id = feature_msg->channels[0].values[i];
        int camera_id = feature_msg->channels[1].values[i];
        double x = feature_msg->points[i].x;
        double y = feature_msg->points[i].y;
        double z = feature_msg->points[i].z;
        double p_u = feature_msg->channels[2].values[i];
        double p_v = feature_msg->channels[3].values[i];
        double velocity_x = feature_msg->channels[4].values[i];
        double velocity_y = feature_msg->channels[5].values[i];
        if(feature_msg->channels.size() > 5)
        {
            double gx = feature_msg->channels[6].values[i];
            double gy = feature_msg->channels[7].values[i];
            double gz = feature_msg->channels[8].values[i];
            pts_gt[feature_id] = Eigen::Vector3d(gx, gy, gz);
            //printf("receive pts gt %d %f %f %f\n", feature_id, gx, gy, gz);
        }
        ROS_ASSERT(z == 1);
        Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
        xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
        featureFrame[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
    }
    double t = feature_msg->header.stamp.toSec();
    estimator.inputFeature(t, featureFrame);
    return;
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while(!feature_buf.empty())
            feature_buf.pop();
        while(!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        estimator.clearState();
        estimator.setParameter();
    }
    return;
}

// we assume IMU and leg have the same timestamp 
void sensor_callback(const sensor_msgs::Imu::ConstPtr& imu_msg, const sensor_msgs::JointState::ConstPtr& joint_msg) {

    // std::cout<<"sensor_callback"<<std::endl;
    double t = imu_msg->header.stamp.toSec();

    // assemble sensor data
    Eigen::Vector3d acc = Eigen::Vector3d(imu_msg->linear_acceleration.x, imu_msg->linear_acceleration.y, imu_msg->linear_acceleration.z);
    Eigen::Vector3d ang_vel = Eigen::Vector3d(imu_msg->angular_velocity.x, imu_msg->angular_velocity.y, imu_msg->angular_velocity.z);

    Eigen::Matrix<double, NUM_DOF,1> joint_pos;
    Eigen::Matrix<double, NUM_DOF,1> joint_vel;
    Eigen::Matrix<double, NUM_LEG,1> plan_contacts;
    Eigen::Matrix<double, NUM_LEG,1> foot_force_sensor_readings;
    for (int i = 0; i < NUM_DOF; ++i) {
        joint_pos[i] = joint_msg->position[i];
        joint_vel[i] = joint_msg->velocity[i];
    }
    for (int i = 0; i < NUM_LEG; ++i) {
        plan_contacts[i] = joint_msg->velocity[NUM_DOF + i];
        foot_force_sensor_readings[i] = joint_msg->effort[NUM_DOF + i];
    }

    double dt;
    data.input_imu(acc, ang_vel);
    data.input_leg(joint_pos, joint_vel, plan_contacts);

    if ( !kf.is_inited() && first_sensor_received == false ) {
        // the callback is called the first time, filter may not be inited
        dt = 0;
        curr_t = t;
        data.input_dt(dt);
        // init the filter 
        kf.init_filter(data);
    } else if ( !kf.is_inited()) {
        // filter may not be inited even after the callback is called multiple times
        dt = t- curr_t;
        data.input_dt(dt);
        curr_t = t;
    } else {
        dt = t- curr_t;
        
        data.input_dt(dt);
        kf.update_filter(data);
        curr_t = t;
    }

    // get filtered data from data and kf to estimator
    estimator.inputIMU(t, data.acc, data.ang_vel);

    if (CONTACT_SENSOR_TYPE == 0) {
        estimator.inputLeg(t, data.joint_pos, data.joint_vel, kf.get_contacts());
    } else if (CONTACT_SENSOR_TYPE == 1) {
        estimator.inputLeg(t, data.joint_pos, data.joint_vel, data.plan_contacts);
    } else if (CONTACT_SENSOR_TYPE == 2) {
        estimator.inputLeg(t, data.joint_pos, data.joint_vel, foot_force_sensor_readings);
    }
    


    // debug print filtered data
    sensor_msgs::Imu filterd_imu_msg;
    sensor_msgs::JointState filterd_joint_msg;
    filterd_imu_msg.header.stamp = ros::Time::now();
    filterd_imu_msg.linear_acceleration.x = data.acc[0];
    filterd_imu_msg.linear_acceleration.y = data.acc[1];
    filterd_imu_msg.linear_acceleration.z = data.acc[2];

    filterd_imu_msg.angular_velocity.x = data.ang_vel[0];
    filterd_imu_msg.angular_velocity.y = data.ang_vel[1];
    filterd_imu_msg.angular_velocity.z = data.ang_vel[2];

    filterd_joint_msg.header.stamp = ros::Time::now();

    filterd_joint_msg.name = {"FL0", "FL1", "FL2",
                           "FR0", "FR1", "FR2",
                           "RL0", "RL1", "RL2",
                           "RR0", "RR1", "RR2",
                           "FL_foot", "FR_foot", "RL_foot", "RR_foot"};
    filterd_joint_msg.position.resize(NUM_DOF + NUM_LEG);
    filterd_joint_msg.velocity.resize(NUM_DOF + NUM_LEG);
    filterd_joint_msg.effort.resize(NUM_DOF + NUM_LEG);
    for (int i = 0; i < NUM_DOF; ++i) {
        filterd_joint_msg.position[i] = data.joint_pos[i];
        filterd_joint_msg.velocity[i] = data.joint_vel[i];
    }
    Eigen::Vector4d estimated_contact = kf.get_contacts();
    for (int i = 0; i < NUM_LEG; ++i) {
        filterd_joint_msg.velocity[NUM_DOF+i] = estimated_contact[i];
    }
    filterd_imu_pub.publish(filterd_imu_msg);
    filterd_joint_pub.publish(filterd_joint_msg);

    Eigen::Matrix<double, EKF_STATE_SIZE,1> kf_state = kf.get_state();
    nav_msgs::Odometry filterd_pos_msg;
    filterd_pos_msg.header.stamp = ros::Time::now();
    filterd_pos_msg.pose.pose.position.x = kf_state[0];
    filterd_pos_msg.pose.pose.position.y = kf_state[1];
    filterd_pos_msg.pose.pose.position.z = kf_state[2];
    filterd_pos_msg.twist.twist.linear.x = kf_state[3];
    filterd_pos_msg.twist.twist.linear.y = kf_state[4];
    filterd_pos_msg.twist.twist.linear.z = kf_state[5];

    filterd_pos_pub.publish(filterd_pos_msg);

    first_sensor_received = true;
    return;

}


// if optitrack data is available, use it as ground truth
// notice this is sort of an asynchronous callback
double opti_dt = 0;
double opti_curr_t = 0;
ros::Publisher filterd_opti_vel_pub;
bool opti_callback_first_received = false;
void opti_callback(const geometry_msgs::PoseStamped::ConstPtr& opti_msg) {
    std::cout<<"opti_callback"<<std::endl;
    double opti_t = opti_msg->header.stamp.toSec();

    Eigen::Matrix<double, 3, 1> opti_pos; 
    opti_pos << opti_msg->pose.position.x, opti_msg->pose.position.y, opti_msg->pose.position.z;

    // update sensor data
    if (opti_callback_first_received == false) {
        opti_curr_t = opti_t;
        opti_dt = 0;
        data.input_opti_dt(opti_dt);
        data.input_opti_pos(opti_pos);
    } else {
        // only send data to KF if it is initialized and optitrack generates reliable vel data

        opti_dt = opti_t - opti_curr_t;
        data.input_opti_dt(opti_dt);
        data.input_opti_pos(opti_pos);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "vilo");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);

    if(argc != 2)
    {
        printf("please intput: rosrun vilo vilo_feature_track_test [config file] \n"
               "for example: rosrun vilo vilo_feature_track_test "
               "~/catkin_ws/src/VINS-Fusion/config/euroc/euroc_stereo_imu_config.yaml \n");
        return 1;
    }

    string config_file = argv[1];
    printf("config_file: %s\n", argv[1]);

    readParameters(config_file);
    estimator.setParameter();

#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    ROS_WARN("waiting for image and imu...");

    /* subscriber */
    ros::Subscriber sub_feature = n.subscribe("/feature_tracker/feature", 2000, feature_callback);
    ros::Subscriber sub_img0 = n.subscribe(IMAGE0_TOPIC, 100, img0_callback);
    ros::Subscriber sub_img1 = n.subscribe(IMAGE1_TOPIC, 100, img1_callback);
    ros::Subscriber sub_restart = n.subscribe("/vins_restart", 100, restart_callback);
    // optitrack as ground truth
    ros::Subscriber opti_sub = n.subscribe("/mocap_node/Robot_1/pose", 30, opti_callback);

    // sync IMU and leg, we assume that IMU and leg, although come as two separate topics, are actually has the same time stamp
    message_filters::Subscriber<sensor_msgs::Imu> imu_sub;
    message_filters::Subscriber<sensor_msgs::JointState> joint_state_sub;
    imu_sub.subscribe(n, IMU_TOPIC, 200);
    joint_state_sub.subscribe(n, LEG_TOPIC, 200);
    typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Imu, sensor_msgs::JointState> MySyncPolicy;
    // ExactTime takes a queue size as its constructor argument, hence MySyncPolicy(10)
    message_filters::Synchronizer<MySyncPolicy> sync(MySyncPolicy(30), imu_sub, joint_state_sub);

    sync.registerCallback(boost::bind(&sensor_callback, _1, _2));
    
    /* publishers */
    filterd_imu_pub = n.advertise<sensor_msgs::Imu>("/a1_filterd_imu", 30);
    filterd_joint_pub = n.advertise<sensor_msgs::JointState>("/a1_filterd_joint", 30);
    filterd_pos_pub = n.advertise<nav_msgs::Odometry>("/a1_filterd_pos", 30);
    registerPub(n);


    std::thread sync_thread{sync_process};
    ros::spin();

    return 0;
}
