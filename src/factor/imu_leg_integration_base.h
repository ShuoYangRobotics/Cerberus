//
// Created by shuoy on 8/23/21.
//

#ifndef VILO_IMU_LEG_INTEGRATION_BASE_H
#define VILO_IMU_LEG_INTEGRATION_BASE_H
#include <algorithm>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/SVD>

#include "../utils/utility.h"
#include "../utils/parameters.h"
#include "../legKinematics/A1Kinematics.h"

#include <ceres/ceres.h>
using namespace Eigen;

#define FOOT_VAR_WINDOW_SIZE 5
typedef Eigen::Triplet<int> Trip;

class IMULegIntegrationBase
{
public:
    IMULegIntegrationBase() = delete;
    IMULegIntegrationBase(const Eigen::Vector3d &_base_v, const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                          const Ref<const Vector_dof> &_phi_0, const Ref<const Vector_dof> &_dphi_0, const Ref<const Vector_leg> &_c_0,
                          const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg, const Ref<const Vector_rho> &_linearized_rho,
                          std::vector<Eigen::VectorXd> _rho_fix_list, const Eigen::Vector3d &_p_br, const Eigen::Matrix3d &_R_br);

    void push_back(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr,
                   const Ref<const Vector_dof> &phi, const Ref<const Vector_dof> &dphi, const Ref<const Vector_leg> &c);
    void repropagate(const Eigen::Vector3d &_linearized_ba,
                     const Eigen::Vector3d &_linearized_bg,
                     const Ref<const Vector_rho> &_linearized_rho);
    void propagate(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
                   const Ref<const Vector_dof> &_phi_1, const Ref<const Vector_dof> &_dphi_1, const Ref<const Vector_leg> &_c_1);

    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, 1> evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi,
                                                           const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                                           const Vector_rho &rhoi,
                                                           const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj,
                                                           const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj,
                                                           const Vector_rho &rhoj);

    void midPointIntegration(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                             const Vector3d &_acc_1, const Vector3d &_gyr_1,
                             const Ref<const Vector_dof> &_phi_0, const Ref<const Vector_dof> &_dphi_0,
                             const Ref<const Vector_leg> &_c_0, const Ref<const Vector_dof> &_phi_1,
                             const Ref<const Vector_dof> &_dphi_1, const Ref<const Vector_leg> &_c_1,
                             const Vector3d &delta_p, const Quaterniond &delta_q,
                             const Vector3d &delta_v, const vector<Eigen::Vector3d> &delta_epsilon, Vector3d &sum_delta_epsilon,
                             const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                             const Ref<const Vector_rho> &linearized_rho, Vector3d &result_delta_p,
                             Quaterniond &result_delta_q, Vector3d &result_delta_v,
                             vector<Eigen::Vector3d> &result_delta_epsilon, Vector3d &result_sum_delta_epsilon,
                             Vector3d &result_linearized_ba, Vector3d &result_linearized_bg,
                             Vector_rho &result_linearized_rho, bool update_jacobian);

    void checkJacobian(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                       const Vector3d &_acc_1, const Vector3d &_gyr_1,
                       const Ref<const Vector_dof> &_phi_0, const Ref<const Vector_dof> &_dphi_0,
                       const Ref<const Vector_leg> &_c_0, const Ref<const Vector_dof> &_phi_1,
                       const Ref<const Vector_dof> &_dphi_1, const Ref<const Vector_leg> &_c_1,
                       const Vector3d &delta_p, const Quaterniond &delta_q,
                       const Vector3d &delta_v, const vector<Eigen::Vector3d> &delta_epsilon,
                       const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                       const Ref<const Vector_rho> &linearized_rho);

    // state size

    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> jacobian, covariance;
    double sum_dt;
    Eigen::Vector3d delta_p;                    // alpha
    Eigen::Quaterniond delta_q;                 // gamma
    Eigen::Vector3d delta_v;                    // beta
    std::vector<Eigen::Vector3d> delta_epsilon; // epsilon, displacement calculated from each leg
    Vector3d sum_delta_epsilon;                 // sum of epsilon, consider foot contact

    // the output of optimization,
    Eigen::Vector3d linearized_ba, linearized_bg;
    Vector_rho linearized_rho;
    Vector4i foot_contact_flag;

private:
    double dt;
    Eigen::Vector3d acc_0, gyr_0;
    Eigen::Vector3d acc_1, gyr_1;

    Vector_dof phi_0, dphi_0; // joint angle, joint velocity each has dim 12
    Vector_leg c_0;

    Vector_dof phi_1, dphi_1; // joint angle, joint velocity each has dim 12
    Vector_leg c_1;           //  contact flag, dim 4

    // hold the very first acc_0, gyr_0, use in repropagate
    const Eigen::Vector3d linearized_acc, linearized_gyr;
    Vector_dof linearized_phi, linearized_dphi;
    Vector_leg linearized_c;

    // variables to filter the contact force to get the contact flag
    Vector_leg foot_force;
    Vector_leg foot_force_min;
    Vector_leg foot_force_max;
    Vector_leg foot_force_contact_threshold;
    // variables to calculate variance of foot force
    Eigen::Matrix<double, 4, FOOT_VAR_WINDOW_SIZE> foot_force_window;
    Vector4i foot_force_window_idx;
    Vector_leg foot_force_var;

    // keep track of whether the foot is not in contact through out the integration
    std::vector<bool> integration_contact_flag;

    //    Eigen::Matrix<double, NOISE_SIZE, NOISE_SIZE> noise;
    Eigen::DiagonalMatrix<double, NOISE_SIZE> noise_diag;
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> step_jacobian;
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> step_V;

    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> F;
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> V;

    std::vector<double> dt_buf;
    std::vector<Eigen::Vector3d> acc_buf;
    std::vector<Eigen::Vector3d> gyr_buf;
    std::vector<Vector_dof> phi_buf;
    std::vector<Vector_dof> dphi_buf;
    std::vector<Vector_leg> c_buf;

    A1Kinematics a1_kin;
    std::vector<Eigen::VectorXd> rho_fix_list;
    // following are some parameters that defines the transformation between IMU frame(b) and robot body frame(r)
    Eigen::Vector3d p_br;
    Eigen::Matrix3d R_br;
};

#endif // VILO_IMU_LEG_INTEGRATION_BASE_H
