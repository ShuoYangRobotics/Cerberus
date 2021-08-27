//
// Created by shuoy on 8/23/21.
//

#include "imu_leg_factor.h"

bool IMULegFactor::Evaluate(const double *const *parameters, double *residuals, double **jacobians) const {
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);
//        Vector12d  rhoi;rhoi.setZero();
    Vector12d  rhoi; rhoi << parameters[1][9], parameters[1][10], parameters[1][11],parameters[1][12],
            parameters[1][13], parameters[1][14], parameters[1][15],parameters[1][16],
            parameters[1][17], parameters[1][18], parameters[1][19],parameters[1][20];

    Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
    Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
    Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);
//        Vector12d  rhoj;rhoj.setZero();
    Vector12d  rhoj; rhoj << parameters[3][9], parameters[3][10], parameters[3][11],parameters[3][12],
            parameters[3][13], parameters[3][14], parameters[3][15],parameters[3][16],
            parameters[3][17], parameters[3][18], parameters[3][19],parameters[3][20];

    Eigen::Map<Eigen::Matrix<double, 39, 1>> residual(residuals);
    residual = il_pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi, rhoi,
                                            Pj, Qj, Vj, Baj, Bgj, rhoj);
    Eigen::Matrix<double, 39, 39> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 39, 39>>(il_pre_integration->covariance.inverse()).matrixL().transpose();
    //sqrt_info.setIdentity();
    residual = sqrt_info * residual;

    if (jacobians)
    {
        double sum_dt = il_pre_integration->sum_dt;
        Eigen::Matrix3d dp_dba = il_pre_integration->jacobian.template block<3, 3>(ILO_P, ILO_BA);
        Eigen::Matrix3d dp_dbg = il_pre_integration->jacobian.template block<3, 3>(ILO_P, ILO_BG);

        Eigen::Matrix3d dq_dbg = il_pre_integration->jacobian.template block<3, 3>(ILO_R, ILO_BG);

        Eigen::Matrix3d dv_dba = il_pre_integration->jacobian.template block<3, 3>(ILO_V, ILO_BA);
        Eigen::Matrix3d dv_dbg = il_pre_integration->jacobian.template block<3, 3>(ILO_V, ILO_BG);

        Eigen::Matrix3d dep1_dbg = il_pre_integration->jacobian.block<3, 3>(ILO_EPS1, ILO_BG);
        Eigen::Matrix3d dep1_drho1 = il_pre_integration->jacobian.block<3, 3>(ILO_EPS1, ILO_RHO1);
        Eigen::Matrix3d dep2_dbg = il_pre_integration->jacobian.block<3, 3>(ILO_EPS2, ILO_BG);
        Eigen::Matrix3d dep2_drho2 = il_pre_integration->jacobian.block<3, 3>(ILO_EPS2, ILO_RHO2);
        Eigen::Matrix3d dep3_dbg = il_pre_integration->jacobian.block<3, 3>(ILO_EPS3, ILO_BG);
        Eigen::Matrix3d dep3_drho3 = il_pre_integration->jacobian.block<3, 3>(ILO_EPS3, ILO_RHO3);
        Eigen::Matrix3d dep4_dbg = il_pre_integration->jacobian.block<3, 3>(ILO_EPS4, ILO_BG);
        Eigen::Matrix3d dep4_drho4 = il_pre_integration->jacobian.block<3, 3>(ILO_EPS4, ILO_RHO4);

        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 39, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            jacobian_pose_i.setZero();

            jacobian_pose_i.block<3, 3>(ILO_P, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(ILO_P, 3) = Utility::skewSymmetric(
                    Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

            Eigen::Quaterniond corrected_delta_q = il_pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - il_pre_integration->linearized_bg));
            jacobian_pose_i.block<3, 3>(ILO_R, 3) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();
            jacobian_pose_i.block<3, 3>(ILO_V, 3) = Utility::skewSymmetric(Qi.inverse() * (G * sum_dt + Vj - Vi));

            // related to epsilon
            jacobian_pose_i.block<3, 3>(ILO_EPS1, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(ILO_EPS2, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(ILO_EPS3, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(ILO_EPS4, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(ILO_EPS1, 3) = Utility::skewSymmetric(Qi.inverse() * (Pj - Pi));
            jacobian_pose_i.block<3, 3>(ILO_EPS2, 3) = Utility::skewSymmetric(Qi.inverse() * (Pj - Pi));
            jacobian_pose_i.block<3, 3>(ILO_EPS3, 3) = Utility::skewSymmetric(Qi.inverse() * (Pj - Pi));
            jacobian_pose_i.block<3, 3>(ILO_EPS4, 3) = Utility::skewSymmetric(Qi.inverse() * (Pj - Pi));


            jacobian_pose_i = sqrt_info * jacobian_pose_i;

            if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
            {
                ROS_WARN("numerical unstable in preintegration");
                //std::cout << sqrt_info << std::endl;
                //ROS_BREAK();
            }
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 39, 21, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
            jacobian_speedbias_i.setZero();
            jacobian_speedbias_i.block<3, 3>(ILO_P, 0) = -Qi.inverse().toRotationMatrix() * sum_dt;
            jacobian_speedbias_i.block<3, 3>(ILO_P, 3) = -dp_dba;
            jacobian_speedbias_i.block<3, 3>(ILO_P, 6) = -dp_dbg;
            //
            jacobian_speedbias_i.block<3, 3>(ILO_R, 6) = -Utility::Qleft(Qj.inverse() * Qi * il_pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg;

            jacobian_speedbias_i.block<3, 3>(ILO_V, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(ILO_V, 3) = -dv_dba;
            jacobian_speedbias_i.block<3, 3>(ILO_V, 6) = -dv_dbg;

            jacobian_speedbias_i.block<3, 3>(ILO_EPS1, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(ILO_EPS1, 6) = -dep1_dbg;
            jacobian_speedbias_i.block<3, 3>(ILO_EPS1, 9) = -dep1_drho1;

            jacobian_speedbias_i.block<3, 3>(ILO_EPS2, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(ILO_EPS2, 6) = -dep2_dbg;
            jacobian_speedbias_i.block<3, 3>(ILO_EPS2, 12) = -dep2_drho2;

            jacobian_speedbias_i.block<3, 3>(ILO_EPS3, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(ILO_EPS3, 6) = -dep3_dbg;
            jacobian_speedbias_i.block<3, 3>(ILO_EPS3, 15) = -dep3_drho3;

            jacobian_speedbias_i.block<3, 3>(ILO_EPS4, 0) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(ILO_EPS4, 6) = -dep4_dbg;
            jacobian_speedbias_i.block<3, 3>(ILO_EPS4, 18) = -dep4_drho4;

            jacobian_speedbias_i.block<3, 3>(ILO_BA, 3) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i.block<3, 3>(ILO_BG, 6) = -Eigen::Matrix3d::Identity();

            jacobian_speedbias_i.block<3, 3>(ILO_RHO1, 9) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i.block<3, 3>(ILO_RHO2, 12) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i.block<3, 3>(ILO_RHO3, 15) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i.block<3, 3>(ILO_RHO4, 18) = -Eigen::Matrix3d::Identity();

            jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;
        }

        if (jacobians[2]) {
            Eigen::Map<Eigen::Matrix<double, 39, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
            jacobian_pose_j.setZero();

            jacobian_pose_j.block<3, 3>(ILO_P, 0) = Qi.inverse().toRotationMatrix();
            Eigen::Quaterniond corrected_delta_q = il_pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - il_pre_integration->linearized_bg));
            jacobian_pose_j.block<3, 3>(ILO_R, 3) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

            jacobian_pose_j = sqrt_info * jacobian_pose_j;
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 39, 21, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
            jacobian_speedbias_j.setZero();

            jacobian_speedbias_j.block<3, 3>(ILO_V, 0) = Qi.inverse().toRotationMatrix();

            jacobian_speedbias_j.block<3, 3>(ILO_BA, 3) = Eigen::Matrix3d::Identity();

            jacobian_speedbias_j.block<3, 3>(ILO_BG, 6) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j.block<3, 3>(ILO_RHO1, 9) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j.block<3, 3>(ILO_RHO2, 12) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j.block<3, 3>(ILO_RHO3, 15) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j.block<3, 3>(ILO_RHO4, 18) = Eigen::Matrix3d::Identity();

            jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

            //ROS_ASSERT(fabs(jacobian_speedbias_j.maxCoeff()) < 1e8);
            //ROS_ASSERT(fabs(jacobian_speedbias_j.minCoeff()) < 1e8);
        }
    }

    return true;
}
