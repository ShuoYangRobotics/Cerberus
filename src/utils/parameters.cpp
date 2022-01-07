/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#include "parameters.h"

double INIT_DEPTH;
double MIN_PARALLAX;
double ACC_N, ACC_W;
double GYR_N, GYR_W;

std::vector<Eigen::Matrix3d> RIC;
std::vector<Eigen::Vector3d> TIC;

Eigen::Vector3d G{0.0, 0.0, 9.8};

double BIAS_ACC_THRESHOLD;
double BIAS_GYR_THRESHOLD;
double SOLVER_TIME;
int NUM_ITERATIONS;
int ESTIMATE_EXTRINSIC;
int ESTIMATE_TD;
int ROLLING_SHUTTER;
std::string EX_CALIB_RESULT_PATH;
std::string VINS_RESULT_PATH;
std::string OUTPUT_FOLDER;
std::string IMU_TOPIC;
int ROW, COL;
double TD;
int NUM_OF_CAM;
int STEREO;
int USE_IMU;
int USE_LEG;
int MULTIPLE_THREAD;
map<int, Eigen::Vector3d> pts_gt;
std::string IMAGE0_TOPIC, IMAGE1_TOPIC;
std::string FISHEYE_MASK;
std::vector<std::string> CAM_NAMES;
std::string LEG_TOPIC;
int NUM_OF_LEG;
int OPTIMIZE_LEG_BIAS;
// temporarily write some parameters here
double PHI_N;
double DPHI_N;
double RHO_XY_N;
double RHO_Z_N;

double V_N_FORCE_THRES_RATIO;
double V_N_MIN;
double V_N_MAX;
double V_N_W1;
double V_N_W2;
double V_N_W3;
double V_N_TERM1_STEEP;
double V_N_TERM2_BOUND_FORCE;
double V_N_TERM3_VEL_DIFF_XY;
double V_N_TERM3_VEL_DIFF_Z;
double V_N_TERM2_VAR_RESCALE;
double V_N_TERM3_DISTANCE_RESCALE;
double V_N_FINAL_RATIO;

double LOWER_LEG_LENGTH;

int MAX_CNT;
int MIN_DIST;
double F_THRESHOLD;
int SHOW_TRACK;
int FLOW_BACK;


template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
    T ans;
    if (n.getParam(name, ans))
    {
        ROS_INFO_STREAM("Loaded " << name << ": " << ans);
    }
    else
    {
        ROS_ERROR_STREAM("Failed to load " << name);
        n.shutdown();
    }
    return ans;
}

void readParameters(std::string config_file)
{
    FILE *fh = fopen(config_file.c_str(),"r");
    if(fh == NULL){
        ROS_WARN("config_file dosen't exist; wrong config_file path");
        ROS_BREAK();
        return;          
    }
    fclose(fh);

    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        std::cerr << "ERROR: Wrong path to settings" << std::endl;
    }

    fsSettings["image0_topic"] >> IMAGE0_TOPIC;
    fsSettings["image1_topic"] >> IMAGE1_TOPIC;
    MAX_CNT = fsSettings["max_cnt"];
    MIN_DIST = fsSettings["min_dist"];
    F_THRESHOLD = fsSettings["F_threshold"];
    SHOW_TRACK = fsSettings["show_track"];
    FLOW_BACK = fsSettings["flow_back"];

    MULTIPLE_THREAD = fsSettings["multiple_thread"];

    USE_IMU = fsSettings["imu"];
    printf("USE_IMU: %d\n", USE_IMU);
    if(USE_IMU)
    {
        fsSettings["imu_topic"] >> IMU_TOPIC;
        printf("IMU_TOPIC: %s\n", IMU_TOPIC.c_str());
        ACC_N = fsSettings["acc_n"];
        ACC_W = fsSettings["acc_w"];
        GYR_N = fsSettings["gyr_n"];
        GYR_W = fsSettings["gyr_w"];
        G.z() = fsSettings["g_norm"];
    }

    USE_LEG = fsSettings["use_leg_odom"];
    printf("USE_LEG: %d\n", USE_LEG);
    if(USE_LEG)
    {
        NUM_OF_LEG = fsSettings["num_of_leg"];
        printf("leg number %d\n", NUM_OF_LEG);
        fsSettings["leg_topic"] >> LEG_TOPIC;
        printf("LEG_TOPIC: %s\n", LEG_TOPIC.c_str());

        OPTIMIZE_LEG_BIAS = fsSettings["optimize_leg_bias"];
        printf("leg optimize_leg_bias %d\n", OPTIMIZE_LEG_BIAS);

        PHI_N = fsSettings["joint_angle_n"];
        DPHI_N = fsSettings["joint_velocity_n"];
        RHO_XY_N = fsSettings["leg_bias_xy_n"];
        RHO_Z_N = fsSettings["leg_bias_z_n"];

        V_N_FORCE_THRES_RATIO = fsSettings["v_n_force_thres_ratio"];
        V_N_MIN = fsSettings["v_n_min"];
        V_N_MAX = fsSettings["v_n_max"];
        V_N_W1 = fsSettings["v_n_w1"];
        V_N_W2 = fsSettings["v_n_w2"];
        V_N_W3 = fsSettings["v_n_w3"];
        V_N_TERM1_STEEP = fsSettings["v_n_term1_steep"];
        V_N_TERM2_BOUND_FORCE = fsSettings["v_n_term2_bound_force"];
        V_N_TERM2_VAR_RESCALE = fsSettings["v_n_term2_var_rescale"];
        V_N_TERM3_VEL_DIFF_XY = fsSettings["v_n_term3_vel_diff_xy"];
        V_N_TERM3_VEL_DIFF_Z = fsSettings["v_n_term3_vel_diff_z"];
        V_N_TERM3_DISTANCE_RESCALE = fsSettings["v_n_term3_distance_rescale"];
        V_N_FINAL_RATIO = fsSettings["v_n_final_ratio"];

        LOWER_LEG_LENGTH = fsSettings["lower_leg_length"];
    }

    SOLVER_TIME = fsSettings["max_solver_time"];
    NUM_ITERATIONS = fsSettings["max_num_iterations"];
    MIN_PARALLAX = fsSettings["keyframe_parallax"];
    MIN_PARALLAX = MIN_PARALLAX / FOCAL_LENGTH;

    fsSettings["output_path"] >> OUTPUT_FOLDER;
    if (OPTIMIZE_LEG_BIAS) {
        VINS_RESULT_PATH = OUTPUT_FOLDER + "/vilo_wb"+ Utility::GetCurrentTimeForFileName() + "-lc-" + to_string(LOWER_LEG_LENGTH) + ".csv";
    } else {
        VINS_RESULT_PATH = OUTPUT_FOLDER + "/vilo_wob"+ Utility::GetCurrentTimeForFileName() + ".csv";
    }
    std::cout << "result path " << VINS_RESULT_PATH << std::endl;
    std::ofstream fout(VINS_RESULT_PATH, std::ios::out);
    fout.close();

    ESTIMATE_EXTRINSIC = fsSettings["estimate_extrinsic"];
    if (ESTIMATE_EXTRINSIC == 2)
    {
        ROS_WARN("have no prior about extrinsic param, calibrate extrinsic param");
        RIC.push_back(Eigen::Matrix3d::Identity());
        TIC.push_back(Eigen::Vector3d::Zero());
        EX_CALIB_RESULT_PATH = OUTPUT_FOLDER + "/extrinsic_parameter.csv";
    }
    else 
    {
        if ( ESTIMATE_EXTRINSIC == 1)
        {
            ROS_WARN(" Optimize extrinsic param around initial guess!");
            EX_CALIB_RESULT_PATH = OUTPUT_FOLDER + "/extrinsic_parameter.csv";
        }
        if (ESTIMATE_EXTRINSIC == 0)
            ROS_WARN(" fix extrinsic param ");

        cv::Mat cv_T;
        fsSettings["body_T_cam0"] >> cv_T;
        Eigen::Matrix4d T;
        cv::cv2eigen(cv_T, T);
        RIC.push_back(T.block<3, 3>(0, 0));
        TIC.push_back(T.block<3, 1>(0, 3));
    } 
    
    NUM_OF_CAM = fsSettings["num_of_cam"];
    printf("camera number %d\n", NUM_OF_CAM);

    if(NUM_OF_CAM != 1 && NUM_OF_CAM != 2)
    {
        printf("num_of_cam should be 1 or 2\n");
        assert(0);
    }


    int pn = config_file.find_last_of('/');
    std::string configPath = config_file.substr(0, pn);
    
    std::string cam0Calib;
    fsSettings["cam0_calib"] >> cam0Calib;
    std::string cam0Path = configPath + "/" + cam0Calib;
    CAM_NAMES.push_back(cam0Path);

    if(NUM_OF_CAM == 2)
    {
        STEREO = 1;
        std::string cam1Calib;
        fsSettings["cam1_calib"] >> cam1Calib;
        std::string cam1Path = configPath + "/" + cam1Calib; 
        //printf("%s cam1 path\n", cam1Path.c_str() );
        CAM_NAMES.push_back(cam1Path);
        
        cv::Mat cv_T;
        fsSettings["body_T_cam1"] >> cv_T;
        Eigen::Matrix4d T;
        cv::cv2eigen(cv_T, T);
        RIC.push_back(T.block<3, 3>(0, 0));
        TIC.push_back(T.block<3, 1>(0, 3));
    }

    INIT_DEPTH = 5.0;
    BIAS_ACC_THRESHOLD = 0.1;
    BIAS_GYR_THRESHOLD = 0.1;

    TD = fsSettings["td"];
    ESTIMATE_TD = fsSettings["estimate_td"];
    if (ESTIMATE_TD)
        ROS_INFO_STREAM("Unsynchronized sensors, online estimate time offset, initial td: " << TD);
    else
        ROS_INFO_STREAM("Synchronized sensors, fix time offset: " << TD);

    ROW = fsSettings["image_height"];
    COL = fsSettings["image_width"];
    ROS_INFO("ROW: %d COL: %d ", ROW, COL);

    if(!USE_IMU)
    {
        ESTIMATE_EXTRINSIC = 0;
        ESTIMATE_TD = 0;
        printf("no imu, fix extrinsic param; no time offset calibration\n");
    }

    // VILEMO requires that we must use IMU and stereo camera
    assert(USE_IMU == 1 && STEREO == 1);

    fsSettings.release();
}
