//
// Created by shuoy on 8/23/21.
//

#include "imu_leg_factor.h"


void IMULegFactor::checkJacobian(const double *const *parameters) {
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);
//        Vector12d  rhoi;rhoi.setZero();
    Vector12d rhoi;
    rhoi << parameters[2][0], parameters[2][1], parameters[2][2], parameters[2][3],
            parameters[2][4], parameters[2][5], parameters[2][6], parameters[2][7],
            parameters[2][8], parameters[2][9], parameters[2][10], parameters[2][11];

    Eigen::Vector3d Pj(parameters[3][0], parameters[3][1], parameters[3][2]);
    Eigen::Quaterniond Qj(parameters[3][6], parameters[3][3], parameters[3][4], parameters[3][5]);

    Eigen::Vector3d Vj(parameters[4][0], parameters[4][1], parameters[4][2]);
    Eigen::Vector3d Baj(parameters[4][3], parameters[4][4], parameters[4][5]);
    Eigen::Vector3d Bgj(parameters[4][6], parameters[4][7], parameters[4][8]);
//        Vector12d  rhoj;rhoj.setZero();
    Vector12d rhoj;
    rhoj << parameters[5][0], parameters[5][1], parameters[5][2], parameters[5][3],
            parameters[5][4], parameters[5][5], parameters[5][6], parameters[5][7],
            parameters[5][8], parameters[5][9], parameters[5][10], parameters[5][11];

    std::vector<int> block_sizes{7,9,12,7,9,12};
    Eigen::Matrix<double, 39, 1> residual;
    double **raw_jacobians = new double *[block_sizes.size()];
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobians;
    jacobians.resize(block_sizes.size());
    for (int xx = 0; xx < static_cast<int>(block_sizes.size()); xx++)
    {
        jacobians[xx].resize(39, block_sizes[xx]);
        raw_jacobians[xx] = jacobians[xx].data();
        //dim += block_sizes[i] == 7 ? 6 : block_sizes[i];
    }

    Evaluate(parameters, residual.data(), raw_jacobians);

    // perturb Pi
    Vector3d turb(0.0001, -0.003, 0.003);
    Eigen::Matrix<double, 7, 1> turb_vec; turb_vec.setZero(); turb_vec.segment<3>(0) = turb;

    Eigen::Matrix<double, 39, 1> turb_residual =
            il_pre_integration->evaluate(Pi+turb, Qi, Vi, Bai, Bgi, rhoi, Pj, Qj, Vj, Baj, Bgj, rhoj);
    Eigen::Matrix<double, 39, 39> sqrt_info2 = Eigen::LLT<Eigen::Matrix<double, 39, 39>>(
            il_pre_integration->covariance.inverse()).matrixL().transpose();
//    sqrt_info2.setIdentity();
    turb_residual = sqrt_info2 * turb_residual;

    Eigen::VectorXd tmp = turb_residual - (jacobians[0]* turb_vec + residual);
    std::cout << "perturb Pi" << std::endl;
    std::cout << tmp.transpose() << std::endl;

    // perturb rhoi
    Vector12d rhoi_rand = 0.003*Eigen::Matrix<double,12,1>::Random();
    turb_residual =
            il_pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi, rhoi+rhoi_rand, Pj, Qj, Vj, Baj, Bgj, rhoj);
    sqrt_info2 = Eigen::LLT<Eigen::Matrix<double, 39, 39>>(
            il_pre_integration->covariance.inverse()).matrixL().transpose();
//    sqrt_info2.setIdentity();
    turb_residual = sqrt_info2 * turb_residual;

    tmp = turb_residual - (jacobians[2]* rhoi_rand + residual);
    std::cout << "perturb rhoi" << std::endl;
    std::cout << tmp.transpose() << std::endl;

    return;
}
