//
// Created by shuoy on 7/19/21.
//

#ifndef VILENS_ESTIMATOR_H
#define VILENS_ESTIMATOR_H
#pragma once
#include <thread>
#include <iomanip>
#include <mutex>
#include <algorithm>
#include <std_msgs/Header.h>
#include <std_msgs/Float32.h>
#include <ceres/ceres.h>
#include <unordered_map>
#include <vector>
#include <deque>
#include <queue>
#include <opencv2/core/eigen.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include "../utils/parameters.h"
#include "../featureTracker/feature_manager.h"
#include "../featureTracker/feature_tracker.h"
#include "../factor/imu_factor.h"
#include "../factor/imu_leg_factor.h"
#include "../factor/pose_local_parameterization.h"
#include "../factor/marginalization_factor.h"
#include "../factor/projectionTwoFrameOneCamFactor.h"
#include "../factor/projectionTwoFrameTwoCamFactor.h"
#include "../factor/projectionOneFrameTwoCamFactor.h"

#include "../initial/solve_5pts.h"
#include "../initial/initial_sfm.h"
#include "../initial/initial_alignment.h"
#include "../initial/initial_ex_rotation.h"

#include "../legKinematics/A1Kinematics.h"

class Estimator {
public:
    Estimator();
    ~Estimator();
    void setParameter();

    // interface
    void initFirstPose(Eigen::Vector3d p, Eigen::Matrix3d r);
    void inputIMU(double t, const Vector3d &linearAcceleration, const Vector3d &angularVelocity);
    void inputLeg(double t, const Eigen::Ref<const Vector_dof>& jointAngles, const Eigen::Ref<const Vector_dof>& jointVels, const Eigen::Ref<const Vector_leg>& footForces);
    void inputFeature(double t, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &featureFrame);
    void inputImage(double t, const cv::Mat &_img, const cv::Mat &_img1 = cv::Mat());
    void processIMU(double t, double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
    void processIMULeg(double t, double dt,
                       const Vector3d &linear_acceleration, const Vector3d &angular_velocity,
                       const Ref<const Vector_dof> &joint_angle, const Ref<const Vector_dof> &joint_velocity,
                       const Ref<const Vector_leg> &foot_contact);
    // // most important function, trigger optimization
    void processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const double header);
    void processMeasurements();
    void changeSensorType(int use_imu, int use_stereo);


    // internal
    void clearState();
//    bool initialStructure();
//    bool visualInitialAlign();
//    bool relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l);
    void slideWindow();
    void slideWindowNew();
    void slideWindowOld();
    void optimization();
    void vector2double();
    void double2vector();
    bool failureDetection();
    bool getIMUInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &accVector,
                                              vector<pair<double, Eigen::Vector3d>> &gyrVector);
    bool getIMUAndLegInterval(double t0, double t1, 
                        vector<pair<double, Eigen::Vector3d>> &accVector,
                        vector<pair<double, Eigen::Vector3d>> &gyrVector,
                        vector<pair<double, Vector_dof>> &jointAngVector,
                        vector<pair<double, Vector_dof>> &jointVelVector,
                        vector<pair<double, Vector_leg>> &contactFlagVector);
//    bool getLegInterval(double t0, double t1,
//                                   vector<pair<double, Eigen::VectorXd>> &jointAngVector,
//                                   vector<pair<double, Eigen::VectorXd>> &jointVelVector,
//                                   vector<pair<double, Eigen::VectorXd>> &contactFlagVector);
    void getPoseInWorldFrame(Eigen::Matrix4d &T);
    void getPoseInWorldFrame(int index, Eigen::Matrix4d &T);
    void predictPtsInNextFrame();
    void outliersRejection(set<int> &removeIndex);
    double reprojectionError(Matrix3d &Ri, Vector3d &Pi, Matrix3d &rici, Vector3d &tici,
                                     Matrix3d &Rj, Vector3d &Pj, Matrix3d &ricj, Vector3d &ticj,
                                     double depth, Vector3d &uvi, Vector3d &uvj);
    void updateLatestStates();
    void fastPredictIMU(double t, Eigen::Vector3d linear_acceleration, Eigen::Vector3d angular_velocity);
    bool IMUAvailable(double t);
    void initFirstIMUPose(vector<pair<double, Eigen::Vector3d>> &accVector);

    // get ground truth data
    void receiveGroundTruthData(Vector3d &P, Quaterniond &Q, Vector3d &V);

    enum SolverFlag
    {
        INITIAL,
        NON_LINEAR
    };

    enum MarginalizationFlag
    {
        MARGIN_OLD = 0,
        MARGIN_SECOND_NEW = 1
    };

    std::mutex mProcess;
    std::mutex mBuf;
    std::mutex mPropagate;
    queue<pair<double, Eigen::Vector3d>> accBuf;
    queue<pair<double, Eigen::Vector3d>> gyrBuf;

    queue<pair<double, Vector_dof>> legAngBufList;
    queue<pair<double, Vector_dof>> legAngVelBufList;
    queue<pair<double, Vector_leg>> contactFlagBufList;

    queue<pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > > featureBuf;
    double prevTime, curTime;
    bool openExEstimation;

    std::thread trackThread;
    std::thread processThread;

    FeatureTracker featureTracker;

    SolverFlag solver_flag;
    MarginalizationFlag  marginalization_flag;
    Vector3d g;


    // the actual solved results, Ps Vs Rs are the pose of the imu link
    Matrix3d ric[2];
    Vector3d tic[2];
    Vector3d        Ps[(WINDOW_SIZE + 1)];
    Vector3d        Vs[(WINDOW_SIZE + 1)];
    Matrix3d        Rs[(WINDOW_SIZE + 1)];
    Vector3d        Bas[(WINDOW_SIZE + 1)];
    Vector3d        Bgs[(WINDOW_SIZE + 1)];
    Eigen::Matrix<double, RHO_OPT_SIZE, 1>        Rho1[(WINDOW_SIZE + 1)];
    Eigen::Matrix<double, RHO_OPT_SIZE, 1>        Rho2[(WINDOW_SIZE + 1)];
    Eigen::Matrix<double, RHO_OPT_SIZE, 1>        Rho3[(WINDOW_SIZE + 1)];
    Eigen::Matrix<double, RHO_OPT_SIZE, 1>        Rho4[(WINDOW_SIZE + 1)];
    double td;

    Matrix3d back_R0, last_R, last_R0;
    Vector3d back_P0, last_P, last_P0;
    double Headers[(WINDOW_SIZE + 1)];

    IntegrationBase *pre_integrations[(WINDOW_SIZE + 1)];
    IMULegIntegrationBase * il_pre_integrations[(WINDOW_SIZE + 1)];
    Vector3d acc_0, gyr_0;
    Vector_dof phi_0, dphi_0;
    Vector_leg c_0;

    vector<double> dt_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> linear_acceleration_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> angular_velocity_buf[(WINDOW_SIZE + 1)];
    vector<Vector_dof> joint_angle_buf[(WINDOW_SIZE + 1)];
    vector<Vector_dof> joint_velocity_buf[(WINDOW_SIZE + 1)];
    vector<Vector_leg> foot_contact_buf[(WINDOW_SIZE + 1)];

    int leg_msg_counter;


    int frame_count;
    int sum_of_outlier, sum_of_back, sum_of_front, sum_of_invalid;
    int inputImageCnt;

    FeatureManager f_manager;
    MotionEstimator m_estimator;
    InitialEXRotation initial_ex_rotation;

    bool first_imu;
    bool is_valid, is_key;
    bool failure_occur;

    vector<Vector3d> point_cloud;
    vector<Vector3d> margin_cloud;
    vector<Vector3d> key_poses;
    double initial_timestamp;


    // variable of optimization
    double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS]; // velocity and IMU bias
    double para_LegBias[WINDOW_SIZE + 1][SIZE_LEG_BIAS];    // leg bias
    double para_Feature[NUM_OF_F][SIZE_FEATURE];
    double para_Ex_Pose[2][SIZE_POSE];
    double para_Retrive_Pose[SIZE_POSE];
    double para_Td[1][1];
    double para_Tr[1][1];

    int loop_window_index;

    MarginalizationInfo *last_marginalization_info;
    vector<double *> last_marginalization_parameter_blocks;

    map<double, ImageFrame> all_image_frame;
    IntegrationBase *tmp_pre_integration;
    IMULegIntegrationBase *tmp_il_pre_integration;

    Eigen::Vector3d initP;
    Eigen::Matrix3d initR;

    double latest_time;  // used in fastPredictIMU to get dt from last function call
    // used in fastPredictIMU and updateLatestStates
    Eigen::Vector3d latest_P, latest_V, latest_Ba, latest_Bg, latest_acc_0, latest_gyr_0;
    Eigen::Quaterniond latest_Q;

    bool initFirstPoseFlag;
    bool initThreadFlag;

    // add leg kinematics
    // the leg kinematics is relative to body frame, which is the center of the robot
    // following are some parameters that defines the transformation between IMU frame(b) and robot body frame(r)
    Eigen::Vector3d p_br;
    Eigen::Matrix3d R_br;
    // for each leg, there is an offset between the body frame and the hip motor (fx, fy)
    double leg_offset_x[4] = {};
    double leg_offset_y[4] = {};
    // for each leg, there is an offset between the body frame and the hip motor (fx, fy)
    double motor_offset[4] = {};
    double upper_leg_length[4] = {};
    double lower_leg_length[4] = {};
    std::vector<Eigen::VectorXd> rho_fix_list;
    A1Kinematics a1_kin;

    // add ground truth pose and orientation
    Vector3d gt_position;
    Quaterniond gt_orientation;
    Vector3d gt_velocity;
    Vector3d lo_velocity;
    Vector12d lo_velocity_each_leg;
    Vector3d lo_velocity_with_bias;
    Vector12d lo_velocity_with_bias_each_leg;
    Vector4d foot_contact_flag;
};


#endif //VILENS_ESTIMATOR_H
