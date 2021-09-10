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

#define RESIDUAL_STATE_SIZE 39
#define NOISE_SIZE 54

class IMULegIntegrationBase {
public:
    IMULegIntegrationBase() = delete;
    IMULegIntegrationBase(const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                          const Ref<const Vector12d>& _phi_0, const Ref<const Vector12d>& _dphi_0, const Ref<const Vector12d>& _c_0,
                          const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg, const Ref<const Vector12d>& _linearized_rho,
                          std::vector<Eigen::VectorXd> _rho_fix_list,  const Eigen::Vector3d &_p_br,  const Eigen::Matrix3d &_R_br);

    void push_back(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr,
                   const Ref<const Vector12d>& phi, const Ref<const Vector12d>& dphi, const Ref<const Vector12d>& c);
    void repropagate(const Eigen::Vector3d &_linearized_ba,
                     const Eigen::Vector3d &_linearized_bg,
                     const Ref<const Vector12d> &_linearized_rho);
    void propagate(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
                   const Ref<const Vector12d>& _phi_1, const Ref<const Vector12d>& _dphi_1, const Ref<const Vector12d>& _c_1);

    Eigen::Matrix<double, 39, 1> evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi,
                                          const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                          const Vector12d &rhoi,
                                          const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj,
                                          const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj,
                                          const Vector12d &rhoj);


    void midPointIntegration(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                             const Vector3d &_acc_1, const Vector3d &_gyr_1,
                             const Ref<const Vector12d> &_phi_0, const Ref<const Vector12d> &_dphi_0,
                             const Ref<const Vector12d> &_c_0, const Ref<const Vector12d> &_phi_1,
                             const Ref<const Vector12d> &_dphi_1, const Ref<const Vector12d> &_c_1,
                             const Vector3d &delta_p, const Quaterniond &delta_q,
                             const Vector3d &delta_v, const vector<Eigen::Vector3d> &delta_epsilon,
                             const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                             const Ref<const Vector12d> &linearized_rho, Vector3d &result_delta_p,
                             Quaterniond &result_delta_q, Vector3d &result_delta_v,
                             vector<Eigen::Vector3d> &result_delta_epsilon,
                             Vector3d &result_linearized_ba, Vector3d &result_linearized_bg,
                             Vector12d &result_linearized_rho, bool update_jacobian);

    void checkJacobian(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                  const Vector3d &_acc_1, const Vector3d &_gyr_1,
                  const Ref<const Vector12d> &_phi_0, const Ref<const Vector12d> &_dphi_0,
                  const Ref<const Vector12d> &_c_0, const Ref<const Vector12d> &_phi_1,
                  const Ref<const Vector12d> &_dphi_1, const Ref<const Vector12d> &_c_1,
                  const Vector3d &delta_p, const Quaterniond &delta_q,
                  const Vector3d &delta_v, const vector<Eigen::Vector3d> &delta_epsilon,
                  const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                  const Ref<const Vector12d> &linearized_rho);

    // state size

    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> jacobian, covariance;
    double sum_dt;
    Eigen::Vector3d delta_p;  // alpha
    Eigen::Quaterniond delta_q;  // gamma
    Eigen::Vector3d delta_v;     // beta
    std::vector<Eigen::Vector3d> delta_epsilon;     // epsilon

    // the output of optimization,
    Eigen::Vector3d linearized_ba, linearized_bg;
    Vector12d linearized_rho;

private:
    double dt;
    Eigen::Vector3d acc_0, gyr_0;
    Eigen::Vector3d acc_1, gyr_1;

    Vector12d phi_0, dphi_0, c_0;  // joint angle, joint velocity, contact force, each has dim 12
    Vector12d phi_1, dphi_1, c_1;  // joint angle, joint velocity, contact force, each has dim 12

    // hold the very first acc_0, gyr_0, use in repropagate
    const Eigen::Vector3d linearized_acc, linearized_gyr;
    Vector12d linearized_phi, linearized_dphi, linearized_c;



    Eigen::Matrix<double, NOISE_SIZE, NOISE_SIZE> noise;
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> step_jacobian;
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> step_V;


    std::vector<double> dt_buf;
    std::vector<Eigen::Vector3d> acc_buf;
    std::vector<Eigen::Vector3d> gyr_buf;
    std::vector<Vector12d> phi_buf;
    std::vector<Vector12d> dphi_buf;
    std::vector<Vector12d> c_buf;


    A1Kinematics a1_kin;
    std::vector<Eigen::VectorXd> rho_fix_list;
    // following are some parameters that defines the transformation between IMU frame(b) and robot body frame(r)
    Eigen::Vector3d p_br;
    Eigen::Matrix3d R_br;
};


#endif //VILEOM_IMU_LEG_INTEGRATION_BASE_H
