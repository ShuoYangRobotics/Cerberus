//
// Created by shuoy on 8/10/21.
//

#include "A1Kinematics.h"

Eigen::Vector3d A1Kinematics::fk(Eigen::Vector3d q, Eigen::VectorXd rho_opt, Eigen::VectorXd rho_fix) {
    Eigen::Vector3d out;
    autoFunc_fk_derive(q.data(), rho_opt.data(), rho_fix.data(), out.data());
    return out;
}

Eigen::MatrixXd A1Kinematics::jac(Eigen::Vector3d q, Eigen::VectorXd rho_opt, Eigen::VectorXd rho_fix) {
    Eigen::MatrixXd mtx(3,3);
    autoFunc_d_fk_dq(q.data(), rho_opt.data(), rho_fix.data(), mtx.data());
    return mtx;
}

Eigen::MatrixXd A1Kinematics::dfk_drho(Eigen::Vector3d q, Eigen::VectorXd rho_opt, Eigen::VectorXd rho_fix) {
    Eigen::MatrixXd mtx(3,3);
    autoFunc_d_fk_dc(q.data(), rho_opt.data(), rho_fix.data(), mtx.data());
    return mtx;
}

Eigen::MatrixXd A1Kinematics::dJ_dq(Eigen::Vector3d q, Eigen::VectorXd rho_opt, Eigen::VectorXd rho_fix) {
    Eigen::MatrixXd mtx(9,3);
    autoFunc_dJ_dq(q.data(), rho_opt.data(), rho_fix.data(), mtx.data());
    return mtx;
}

Eigen::MatrixXd A1Kinematics::dJ_drho(Eigen::Vector3d q, Eigen::VectorXd rho_opt, Eigen::VectorXd rho_fix) {
    Eigen::MatrixXd mtx(9,3);
    autoFunc_dJ_dpho(q.data(), rho_opt.data(), rho_fix.data(), mtx.data());
    return mtx;
}


// functions generated by matlab
void A1Kinematics::autoFunc_fk_derive(const double in1[3], const double in2[3], const double
in3[5], double p_bf[3])
{
    double p_bf_tmp;
    double t2;
    double t3;
    double t4;
    double t5;
    double t6;
    double t7;
    double t8;
    double t9;

    //     This function was generated by the Symbolic Math Toolbox version 8.6.
    //     10-Aug-2021 14:48:21
    t2 = std::cos(in1[0]);
    t3 = std::cos(in1[1]);
    t4 = std::cos(in1[2]);
    t5 = std::sin(in1[0]);
    t6 = std::sin(in1[1]);
    t7 = std::sin(in1[2]);
    t8 = in1[1] + in1[2];
    t9 = std::sin(t8);
    p_bf[0] = (((in3[0] + in2[2] * t9) - in3[4] * t9) - t6 * in3[3]) + in2[0] *
                                                                       std::cos(t8);
    p_bf[1] = ((((((((in3[1] + in2[1] * t2) + in3[2] * t2) + t3 * t5 * in3[3]) +
                   in2[0] * t3 * t5 * t7) + in2[0] * t4 * t5 * t6) - in2[2] * t3 *
                                                                     t4 * t5) + in2[2] * t5 * t6 * t7) + in3[4] * t3 * t4 * t5) - in3
                                                                                                                                  [4] * t5 * t6 * t7;
    t8 = in2[0] * t2;
    t9 = in2[2] * t2;
    p_bf_tmp = in3[4] * t2;
    p_bf[2] = (((((((in2[1] * t5 + in3[2] * t5) - t2 * t3 * in3[3]) - t8 * t3 * t7)
                  - t8 * t4 * t6) + t9 * t3 * t4) - t9 * t6 * t7) - p_bf_tmp * t3 *
                                                                    t4) + p_bf_tmp * t6 * t7;
}

void A1Kinematics::autoFunc_d_fk_dq(const double in1[3], const double in2[3], const double
in3[5], double jacobian[9])
{
    double b_jacobian_tmp;
    double c_jacobian_tmp;
    double jacobian_tmp;
    double t12;
    double t16;
    double t17;
    double t2;
    double t22;
    double t3;
    double t4;
    double t5;
    double t6;
    double t7;
    double t8;
    double t9;

    //     This function was generated by the Symbolic Math Toolbox version 8.6.
    //     10-Aug-2021 14:48:21
    t2 = std::cos(in1[0]);
    t3 = std::cos(in1[1]);
    t4 = std::cos(in1[2]);
    t5 = std::sin(in1[0]);
    t6 = std::sin(in1[1]);
    t7 = std::sin(in1[2]);
    t8 = in1[1] + in1[2];
    t9 = std::cos(t8);
    t8 = std::sin(t8);
    t12 = in2[0] * t9;
    t16 = in2[2] * t8;
    t17 = in3[4] * t8;
    t22 = (t12 + t16) + -t17;
    jacobian[0] = 0.0;
    jacobian_tmp = in2[0] * t2;
    b_jacobian_tmp = in2[2] * t2;
    c_jacobian_tmp = in3[4] * t2;
    jacobian[1] = (((((((-in2[1] * t5 - in3[2] * t5) + t2 * t3 * in3[3]) +
                       jacobian_tmp * t3 * t7) + jacobian_tmp * t4 * t6) -
                     b_jacobian_tmp * t3 * t4) + b_jacobian_tmp * t6 * t7) +
                   c_jacobian_tmp * t3 * t4) - c_jacobian_tmp * t6 * t7;
    jacobian[2] = (((((((in2[1] * t2 + in3[2] * t2) + t3 * t5 * in3[3]) + in2[0] *
                                                                          t3 * t5 * t7) + in2[0] * t4 * t5 * t6) - in2[2] * t3 * t4 *
                                                                                                                   t5) + in2[2] * t5 * t6 * t7) + in3[4] * t3 * t4 * t5) - in3[4]
                                                                                                                                                                           * t5 * t6 * t7;
    jacobian_tmp = (in2[2] * t9 + -(in3[4] * t9)) + -(in2[0] * t8);
    jacobian[3] = jacobian_tmp - t3 * in3[3];
    b_jacobian_tmp = ((t6 * in3[3] - t12) - t16) + t17;
    jacobian[4] = -t5 * b_jacobian_tmp;
    jacobian[5] = t2 * b_jacobian_tmp;
    jacobian[6] = jacobian_tmp;
    jacobian[7] = t5 * t22;
    jacobian[8] = -t2 * t22;
}

void A1Kinematics::autoFunc_d_fk_dc(const double in1[3], const double [3], const double [5],
                      double d_fk_dc[9])
{
    double t2;
    double t3;
    double t4;
    double t5;

    //     This function was generated by the Symbolic Math Toolbox version 8.6.
    //     10-Aug-2021 14:48:21
    t2 = std::cos(in1[0]);
    t3 = std::sin(in1[0]);
    t4 = in1[1] + in1[2];
    t5 = std::cos(t4);
    t4 = std::sin(t4);
    d_fk_dc[0] = t5;
    d_fk_dc[1] = t3 * t4;
    d_fk_dc[2] = -t2 * t4;
    d_fk_dc[3] = 0.0;
    d_fk_dc[4] = t2;
    d_fk_dc[5] = t3;
    d_fk_dc[6] = t4;
    d_fk_dc[7] = -t3 * t5;
    d_fk_dc[8] = t2 * t5;
}

void A1Kinematics::autoFunc_dJ_dq(const double in1[3], const double in2[3], const double in3[5],
                    double dJ_dq[27])
{
    double t10;
    double t13;
    double t17;
    double t18;
    double t2;
    double t26;
    double t27;
    double t28;
    double t3;
    double t31;
    double t34;
    double t34_tmp;
    double t35;
    double t4;
    double t5;
    double t6;
    double t7;
    double t8;
    double t9;

    //     This function was generated by the Symbolic Math Toolbox version 8.6.
    //     10-Aug-2021 14:48:21
    t2 = std::cos(in1[0]);
    t3 = std::cos(in1[1]);
    t4 = std::cos(in1[2]);
    t5 = std::sin(in1[0]);
    t6 = std::sin(in1[1]);
    t7 = std::sin(in1[2]);
    t8 = in1[1] + in1[2];
    t9 = std::cos(t8);
    t10 = t3 * in3[3];
    t8 = std::sin(t8);
    t13 = in2[0] * t9;
    t17 = in2[2] * t8;
    t18 = in3[4] * t8;
    t9 = (in3[4] * t9 + in2[0] * t8) + -(in2[2] * t9);
    t8 = (t13 + t17) + -t18;
    t26 = (t18 + -t13) + -t17;
    t27 = t5 * t8;
    t28 = t2 * t9;
    t18 = t2 * t8;
    t31 = t10 + t9;
    t34_tmp = t6 * in3[3] + t26;
    t34 = -t2 * t34_tmp;
    t35 = -t5 * t34_tmp;
    t8 = -(t5 * t9);
    dJ_dq[0] = 0.0;
    dJ_dq[1] = (((((((-in2[1] * t2 - in3[2] * t2) - t5 * t10) - in2[0] * t3 * t5 *
                                                                t7) - in2[0] * t4 * t5 * t6) + in2[2] * t3 * t4 * t5) - in2[2]
                                                                                                                        * t5 * t6 * t7) - in3[4] * t3 * t4 * t5) + in3[4] * t5 * t6 * t7;
    t9 = in2[0] * t2;
    t13 = in2[2] * t2;
    t17 = in3[4] * t2;
    dJ_dq[2] = (((((((-in2[1] * t5 - in3[2] * t5) + t2 * t10) + t9 * t3 * t7) + t9
                                                                                * t4 * t6) - t13 * t3 * t4) + t13 * t6 * t7) + t17 * t3 * t4) -
               t17 * t6 * t7;
    dJ_dq[3] = 0.0;
    dJ_dq[4] = t34;
    dJ_dq[5] = t35;
    dJ_dq[6] = 0.0;
    dJ_dq[7] = t18;
    dJ_dq[8] = t27;
    dJ_dq[9] = 0.0;
    dJ_dq[10] = t34;
    dJ_dq[11] = t35;
    dJ_dq[12] = t34_tmp;
    dJ_dq[13] = -t5 * t31;
    dJ_dq[14] = t2 * t31;
    dJ_dq[15] = t26;
    dJ_dq[16] = t8;
    dJ_dq[17] = t28;
    dJ_dq[18] = 0.0;
    dJ_dq[19] = t18;
    dJ_dq[20] = t27;
    dJ_dq[21] = t26;
    dJ_dq[22] = t8;
    dJ_dq[23] = t28;
    dJ_dq[24] = t26;
    dJ_dq[25] = t8;
    dJ_dq[26] = t28;
}

void A1Kinematics::autoFunc_dJ_dpho(const double in1[3], const double [3], const double [5],
                      double dJ_dpho[27])
{
    double t11;
    double t12;
    double t13;
    double t2;
    double t3;
    double t4;
    double t5;
    double t9;

    //     This function was generated by the Symbolic Math Toolbox version 8.6.
    //     10-Aug-2021 14:48:22
    t2 = std::cos(in1[0]);
    t3 = std::sin(in1[0]);
    t4 = in1[1] + in1[2];
    t5 = std::cos(t4);
    t4 = std::sin(t4);
    t9 = t3 * t5;
    t11 = t3 * t4;
    t12 = -(t2 * t5);
    t13 = t2 * -t4;
    dJ_dpho[0] = 0.0;
    dJ_dpho[1] = t2 * t4;
    dJ_dpho[2] = t11;
    dJ_dpho[3] = -t4;
    dJ_dpho[4] = t9;
    dJ_dpho[5] = t12;
    dJ_dpho[6] = -t4;
    dJ_dpho[7] = t9;
    dJ_dpho[8] = t12;
    dJ_dpho[9] = 0.0;
    dJ_dpho[10] = -t3;
    dJ_dpho[11] = t2;
    dJ_dpho[12] = 0.0;
    dJ_dpho[13] = 0.0;
    dJ_dpho[14] = 0.0;
    dJ_dpho[15] = 0.0;
    dJ_dpho[16] = 0.0;
    dJ_dpho[17] = 0.0;
    dJ_dpho[18] = 0.0;
    dJ_dpho[19] = t12;
    dJ_dpho[20] = -t9;
    dJ_dpho[21] = t5;
    dJ_dpho[22] = t11;
    dJ_dpho[23] = t13;
    dJ_dpho[24] = t5;
    dJ_dpho[25] = t11;
    dJ_dpho[26] = t13;
}