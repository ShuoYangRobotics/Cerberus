//
// Created by shuoy on 7/19/21.
//

#include <ceres/ceres.h>

#include "../legKinematics/A1Kinematics.h"

int main(int argc, char **argv) {

    A1Kinematics a1_kin;
    Eigen::Vector3d q(0.1,0.1,0.3);
    Eigen::VectorXd rho_fix(5); rho_fix << 0.3,0.15,0.08,0.2,0.2;
    Eigen::VectorXd rho_opt(3); rho_opt << 0.0,0.0,0.0;
    Eigen::Vector3d p_bf = a1_kin.fk(q, rho_opt,rho_fix);
    std::cout << "forward kinematics \n" << p_bf << std::endl;
    Eigen::MatrixXd Jac = a1_kin.jac(q, rho_opt,rho_fix);
    std::cout << "Jacobian \n" <<Jac << std::endl;

    Eigen::MatrixXd dfk_drho = a1_kin.dfk_drho(q, rho_opt,rho_fix);
    std::cout << "dfk_drho \n" <<dfk_drho << std::endl;

    Eigen::MatrixXd dJ_dq = a1_kin.dJ_dq(q, rho_opt,rho_fix);
    std::cout << "dJ_dq \n" <<dJ_dq << std::endl;
    Eigen::MatrixXd dJ_drho = a1_kin.dJ_drho(q, rho_opt,rho_fix);
    std::cout << "dJ_drho \n" <<dJ_drho << std::endl;
    return 0;
}