/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#pragma once

#include <ros/ros.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <fstream>
#include <map>

using namespace std;

const double FOCAL_LENGTH = 460.0;
const int WINDOW_SIZE = 10;
const int NUM_OF_F = 1000;
//#define UNIT_SPHERE_ERROR

extern double INIT_DEPTH;
extern double MIN_PARALLAX;
extern int ESTIMATE_EXTRINSIC;

extern double ACC_N, ACC_N_Z, ACC_W;
extern double GYR_N, GYR_W;

extern std::vector<Eigen::Matrix3d> RIC;     // num of cam, imu to camera rotation
extern std::vector<Eigen::Vector3d> TIC;     // num of cam, imu to camera position
extern Eigen::Vector3d G;

extern double BIAS_ACC_THRESHOLD;
extern double BIAS_GYR_THRESHOLD;
extern double SOLVER_TIME;
extern int NUM_ITERATIONS;
extern std::string EX_CALIB_RESULT_PATH;
extern std::string VINS_RESULT_PATH;
extern std::string OUTPUT_FOLDER;
extern std::string IMU_TOPIC;
extern double TD;
extern int ESTIMATE_TD;
extern int ROLLING_SHUTTER;
extern int ROW, COL;
extern int NUM_OF_CAM;
extern int STEREO;
extern int USE_IMU;
extern std::string LEG_TOPIC;
extern int NUM_OF_LEG;
extern int USE_LEG;
extern int OPTIMIZE_LEG_BIAS;
// temporarily write some parameters here
extern double PHI_N;
extern double DPHI_N;
extern double RHO_XY_N;
extern double RHO_Z_N;
extern double V_N_FORCE_THRES_RATIO;
extern double V_N_MIN_XY;
extern double V_N_MIN_Z;
extern double V_N_MIN;
extern double V_N_MAX;
extern double V_N_W1;
extern double V_N_W2;
extern double V_N_W3;
extern double V_N_TERM1_STEEP;
extern double V_N_TERM2_BOUND_FORCE;
extern double V_N_TERM3_VEL_DIFF_XY;
extern double V_N_TERM3_VEL_DIFF_Z;
extern double V_N_TERM2_VAR_RESCALE;
extern double V_N_TERM3_DISTANCE_RESCALE;
extern double V_N_FINAL_RATIO;

extern double VINS_LOWER_LEG_LENGTH;

extern int MULTIPLE_THREAD;
// pts_gt for debug purpose;
extern map<int, Eigen::Vector3d> pts_gt;

extern std::string IMAGE0_TOPIC, IMAGE1_TOPIC;
extern std::string FISHEYE_MASK;
extern std::vector<std::string> CAM_NAMES;
extern int MAX_CNT;
extern int MIN_DIST;
extern double F_THRESHOLD;
extern int SHOW_TRACK;
extern int FLOW_BACK;

void readParameters(std::string config_file);

#define NUM_OF_LEG 4
#define NUM_OF_DOF 12
// parameters in the leg kinematics  and imu_leg_integration_base
#define RHO_OPT_SIZE  1
#define TOTAL_RHO_OPT_SIZE  4   //4xRHO_OPT_SIZE
#define RHO_FIX_SIZE  4
#define D_FK_DRHO_SIZE  3   // 3xRHO_OPT_SIZE
#define D_J_DRHO_SIZE  9    // 9xRHO_OPT_SIZE
#define RESIDUAL_STATE_SIZE 31  // 3*9 + 4xRHO_OPT_SIZE
#define NOISE_SIZE 46           // 3*14 + 4xRHO_OPT_SIZE

typedef Eigen::Matrix<double, 12, 1> Vector12d; //4xRHO_OPT_SIZE
typedef Eigen::Matrix<double, 4, 1> Vector_rho; //4xRHO_OPT_SIZE

typedef Eigen::Matrix<double, NUM_OF_LEG, 1> Vector_leg;
typedef Eigen::Matrix<double, NUM_OF_DOF, 1> Vector_dof;

enum SIZE_PARAMETERIZATION
{
    SIZE_POSE = 7,       // p3, q4
    SIZE_SPEEDBIAS = 9,
    SIZE_LEG_BIAS = 4,  // 4 x RHO_OPT_SIZE
    SIZE_FEATURE = 1
};

enum StateOrder
{
    O_P = 0,
    O_R = 3,
    O_V = 6,
    O_BA = 9,
    O_BG = 12
};

enum NoiseOrder
{
    O_AN = 0,
    O_GN = 3,
    O_AW = 6,
    O_GW = 9
};

enum ILStateOrder // error state, total is RESIDUAL_STATE_SIZE
{
    ILO_P = 0,
    ILO_R = 3,
    ILO_V = 6,
    ILO_EPS1 = 9,
    ILO_EPS2 = 12,
    ILO_EPS3 = 15,
    ILO_EPS4 = 18,
    ILO_BA = 21,
    ILO_BG = 24,
    ILO_RHO1 = 27,     // change according to RHO_OPT_SIZE
    ILO_RHO2 = 28,     // change according to RHO_OPT_SIZE
    ILO_RHO3 = 29,     // change according to RHO_OPT_SIZE
    ILO_RHO4 = 30,     // change according to RHO_OPT_SIZE
};

enum ILNoiseStateOrder // noise state, total is NOISE_SIZE
{
    ILNO_Ai = 0,
    ILNO_Gi = 3,
    ILNO_Ai1 = 6,
    ILNO_Gi1 = 9,
    ILNO_BA = 12,
    ILNO_BG = 15,
    ILNO_PHIi = 18,
    ILNO_PHIi1 = 21,
    ILNO_DPHIi = 24,
    ILNO_DPHIi1 = 27,
    ILNO_V1 = 30,
    ILNO_V2 = 33,
    ILNO_V3 = 36,
    ILNO_V4 = 39,
    ILNO_NRHO1 = 42,
    ILNO_NRHO2 = 43,
    ILNO_NRHO3 = 44,
    ILNO_NRHO4 = 45,
};
