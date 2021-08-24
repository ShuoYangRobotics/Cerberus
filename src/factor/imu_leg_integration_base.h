//
// Created by shuoy on 8/23/21.
//

#ifndef VILEOM_IMU_LEG_INTEGRATION_BASE_H
#define VILEOM_IMU_LEG_INTEGRATION_BASE_H
#include<algorithm>
#include<cmath>
#include "../utils/utility.h"
#include "../utils/parameters.h"
#include "../legKinematics/A1Kinematics.h"

#include <ceres/ceres.h>
using namespace Eigen;

class IMULegIntegrationBase {
public:
    IMULegIntegrationBase() = delete;
    IMULegIntegrationBase(const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                          const Ref<const VectorXd>& _phi_0, const Ref<const VectorXd>& _dphi_0, const Ref<const VectorXd>& _c_0,
                          const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg, const Ref<const VectorXd>& _linearized_rho,
                          std::vector<Eigen::VectorXd> _rho_fix_list,  const Eigen::Vector3d &_p_br,  const Eigen::Matrix3d &_R_br);

    void push_back(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr,
                   const Ref<const VectorXd>& phi, const Ref<const VectorXd>& dphi, const Ref<const VectorXd>& c);

    void propagate(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
                   const Ref<const VectorXd>& _phi_1, const Ref<const VectorXd>& _dphi_1, const Ref<const VectorXd>& _c_1);

    void midPointIntegration(double _dt,
                             const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                             const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
                             const Ref<const VectorXd> &_phi_0, const Ref<const VectorXd> &_dphi_0, const Ref<const VectorXd> &_c_0,
                             const Ref<const VectorXd> &_phi_1, const Ref<const VectorXd> &_dphi_1, const Ref<const VectorXd> &_c_1,
                             const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
                             const std::vector<Eigen::Vector3d>& delta_epsilon,
                             const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg, const Ref<const VectorXd> &linearized_rho,
                             Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
                             std::vector<Eigen::Vector3d> &result_delta_epsilon,
                             Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, Ref<Eigen::VectorXd> &result_linearized_rho,
                             bool update_jacobian);

private:
    double dt;
    Eigen::Vector3d acc_0, gyr_0;
    Eigen::Vector3d acc_1, gyr_1;

    Eigen::VectorXd phi_0, dphi_0, c_0;  // joint angle, joint velocity, contact force, each has dim 12
    Eigen::VectorXd phi_1, dphi_1, c_1;  // joint angle, joint velocity, contact force, each has dim 12

    // hold the very first acc_0, gyr_0, use in repropagate
    const Eigen::Vector3d linearized_acc, linearized_gyr;
    Eigen::VectorXd linearized_phi, linearized_dphi, linearized_c;

    // the output of optimization,
    Eigen::Vector3d linearized_ba, linearized_bg;
    Eigen::VectorXd linearized_rho;

    // state size
    const static int STATE_SIZE = 39;
    const static int NOISE_SIZE = 45;

    Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> jacobian, covariance;
    Eigen::Matrix<double, NOISE_SIZE, NOISE_SIZE> noise;

    double sum_dt;
    Eigen::Vector3d delta_p;  // alpha
    Eigen::Quaterniond delta_q;  // gamma
    Eigen::Vector3d delta_v;     // beta
    std::vector<Eigen::Vector3d> delta_epsilon;     // epsilon

    std::vector<double> dt_buf;
    std::vector<Eigen::Vector3d> acc_buf;
    std::vector<Eigen::Vector3d> gyr_buf;
    std::vector<Eigen::VectorXd> phi_buf;
    std::vector<Eigen::VectorXd> dphi_buf;
    std::vector<Eigen::VectorXd> c_buf;


    A1Kinematics a1_kin;
    std::vector<Eigen::VectorXd> rho_fix_list;
    // following are some parameters that defines the transformation between IMU frame(b) and robot body frame(r)
    Eigen::Vector3d p_br;
    Eigen::Matrix3d R_br;
};


#endif //VILEOM_IMU_LEG_INTEGRATION_BASE_H
