//
// Created by shuoy on 8/23/21.
//

#include "imu_leg_integration_base.h"

IMULegIntegrationBase::IMULegIntegrationBase(const Vector3d &_base_v, const Vector3d &_acc_0, const Vector3d &_gyr_0, const Ref<const Vector12d>& _phi_0,
                                             const Ref<const Vector12d>& _dphi_0, const Ref<const Vector12d>& _c_0,
                                             const Vector3d &_linearized_ba, const Vector3d &_linearized_bg, const Vector3d &_linearized_bv, 
                                             const Ref<const Vector_rho>& _linearized_rho,
                                             std::vector<Eigen::VectorXd> _rho_fix_list, const Eigen::Vector3d &_p_br,  const Eigen::Matrix3d &_R_br)
        : acc_0{_acc_0}, gyr_0{_gyr_0}, linearized_acc{_acc_0}, linearized_gyr{_gyr_0},
          linearized_ba{_linearized_ba}, linearized_bg{_linearized_bg}, linearized_bv{_linearized_bv},
          sum_dt{0.0}, delta_p{Eigen::Vector3d::Zero()}, delta_q{Eigen::Quaterniond::Identity()}, delta_v{Eigen::Vector3d::Zero()}
{
    jacobian.setIdentity();
    covariance.setZero();

    base_v = _base_v;

    phi_0 = _phi_0;
    dphi_0 = _dphi_0;
    c_0 = _c_0;

    linearized_phi = _phi_0;
    linearized_dphi = _dphi_0;
    linearized_c = _c_0;

    linearized_rho = _linearized_rho;

    foot_force_min.setZero();
    foot_force_max.setZero();

    for (int j = 0; j < NUM_OF_LEG; j++) {
        integration_contact_flag.push_back(true);
    }

    foot_force_window.setZero();
    foot_force_window_idx.setZero();
    foot_force_var.setZero();


    // the fixed kinematics parameter
    rho_fix_list = _rho_fix_list;
    p_br = _p_br;
    R_br = _R_br;


    noise_diag.diagonal() <<
            (ACC_N * ACC_N), (ACC_N * ACC_N), (ACC_N * ACC_N),
            (GYR_N * GYR_N), (GYR_N * GYR_N), (GYR_N * GYR_N),
            (ACC_N * ACC_N), (ACC_N * ACC_N), (ACC_N * ACC_N),
            (GYR_N * GYR_N), (GYR_N * GYR_N), (GYR_N * GYR_N),
            (ACC_W * ACC_W), (ACC_W * ACC_W), (ACC_W * ACC_W),
            (GYR_W * GYR_W), (GYR_W * GYR_W), (GYR_W * GYR_W),
            (LBV_W * LBV_W), (LBV_W * LBV_W), (LBV_W * LBV_W),
            (PHI_N * PHI_N), (PHI_N * PHI_N), (PHI_N * PHI_N),
            (PHI_N * PHI_N), (PHI_N * PHI_N), (PHI_N * PHI_N),
            (DPHI_N * DPHI_N), (DPHI_N * DPHI_N), (DPHI_N * DPHI_N),
            (DPHI_N * DPHI_N), (DPHI_N * DPHI_N), (DPHI_N * DPHI_N),
            (LBV_N * LBV_N), (LBV_N * LBV_N), (LBV_N * LBV_N),
            RHO_N, RHO_N, RHO_N, RHO_N;

}

void IMULegIntegrationBase::push_back(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr,
                                      const Ref<const Vector12d>& phi, const Ref<const Vector12d>& dphi, const Ref<const Vector12d>& c) {
    dt_buf.push_back(dt);
    acc_buf.push_back(acc);
    gyr_buf.push_back(gyr);
    phi_buf.push_back(phi);
    dphi_buf.push_back(dphi);
    c_buf.push_back(c);
    propagate(dt, acc, gyr, phi, dphi, c);
}

// repropagate uses new bias
void IMULegIntegrationBase::repropagate(
        const Eigen::Vector3d &_linearized_ba, 
        const Eigen::Vector3d &_linearized_bg,  
        const Eigen::Vector3d &_linearized_bv,
        const Ref<const Vector_rho> &_linearized_rho)
{
    sum_dt = 0.0;
    acc_0 = linearized_acc;
    gyr_0 = linearized_gyr;
    phi_0 = linearized_phi;
    dphi_0 = linearized_dphi;
    c_0 = linearized_c;

    delta_p.setZero();
    delta_q.setIdentity();
    delta_v.setZero();
    delta_epsilon.setZero();

    linearized_ba = _linearized_ba;
    linearized_bg = _linearized_bg;
    linearized_bv = _linearized_bv;
    linearized_rho = _linearized_rho;
    jacobian.setIdentity();
    covariance.setZero();
    for (int i = 0; i < static_cast<int>(dt_buf.size()); i++)
        propagate(dt_buf[i], acc_buf[i], gyr_buf[i], phi_buf[i], dphi_buf[i], c_buf[i]);
}


void IMULegIntegrationBase::propagate(double _dt, const Vector3d &_acc_1, const Vector3d &_gyr_1,
                                      const Ref<const Vector12d>& _phi_1, const Ref<const Vector12d>& _dphi_1, const Ref<const Vector12d>& _c_1) {
    dt = _dt;
    acc_1 = _acc_1;
    gyr_1 = _gyr_1;
    phi_1 = _phi_1;
    dphi_1 = _dphi_1;
    c_1 = _c_1;
    Vector3d result_delta_p;
    Quaterniond result_delta_q;
    Vector3d result_delta_v;
    Vector3d result_delta_epsilon;
    Vector3d result_sum_delta_epsilon; 
    Vector3d result_linearized_ba;
    Vector3d result_linearized_bg;
    Vector3d result_linearized_bv;
    Vector_rho result_linearized_rho;
    // midPointIntegration
    midPointIntegration(_dt, acc_0, gyr_0, acc_1, gyr_1,
                        phi_0, dphi_0, c_0, phi_1, dphi_1, c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        result_delta_p, result_delta_q, result_delta_v, result_delta_epsilon, 
                        result_linearized_ba, result_linearized_bg, result_linearized_bv, result_linearized_rho, 1);
    // checkJacobian
//    checkJacobian(_dt, acc_0, gyr_0, acc_1, gyr_1,
//                        phi_0, dphi_0, c_0, phi_1, dphi_1, c_1,
//                        delta_p, delta_q, delta_v, delta_epsilon,
//                        linearized_ba, linearized_bg, linearized_bv, linearized_rho);


    // std::cout << "inside propagate" << std::endl;
    // std::cout << delta_epsilon.transpose() << std::endl;
    // std::cout << result_delta_epsilon.transpose() << std::endl;

    delta_p = result_delta_p;
    delta_q = result_delta_q;
    delta_v = result_delta_v;
    delta_epsilon = result_delta_epsilon;

    linearized_ba = result_linearized_ba;
    linearized_bg = result_linearized_bg;
    linearized_bv = result_linearized_bv;
    linearized_rho = result_linearized_rho;
    delta_q.normalize();
    sum_dt += dt;
    acc_0 = acc_1;
    gyr_0 = gyr_1;
    phi_0 = phi_1;
    dphi_0 = dphi_1;
    c_0 = c_1;

}



void IMULegIntegrationBase::midPointIntegration(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                             const Vector3d &_acc_1, const Vector3d &_gyr_1,
                             const Ref<const Vector12d> &_phi_0, const Ref<const Vector12d> &_dphi_0,
                             const Ref<const Vector12d> &_c_0, const Ref<const Vector12d> &_phi_1,
                             const Ref<const Vector12d> &_dphi_1, const Ref<const Vector12d> &_c_1,
                             const Vector3d &delta_p, const Quaterniond &delta_q,
                             const Vector3d &delta_v, const Vector3d delta_epsilon, 
                             const Vector3d &linearized_ba, const Vector3d &linearized_bg, const Vector3d &linearized_bv,
                             const Ref<const Vector_rho> &linearized_rho, Vector3d &result_delta_p,
                             Quaterniond &result_delta_q, Vector3d &result_delta_v,
                             Vector3d &result_delta_epsilon, 
                             Vector3d &result_linearized_ba, Vector3d &result_linearized_bg, Vector3d &result_linearized_bv,
                             Vector_rho &result_linearized_rho, bool update_jacobian) {
    Vector3d un_acc_0 = delta_q * (_acc_0 - linearized_ba);
    Vector3d un_gyr = 0.5 * (_gyr_0 + _gyr_1) - linearized_bg;
    result_delta_q = delta_q * Quaterniond(1, un_gyr(0) * _dt / 2, un_gyr(1) * _dt / 2, un_gyr(2) * _dt / 2);
    Vector3d un_acc_1 = result_delta_q * (_acc_1 - linearized_ba);
    Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    result_delta_p = delta_p + delta_v * _dt + 0.5 * un_acc * _dt * _dt;
    result_delta_v = delta_v + un_acc * _dt;

    // some helper structs
    Vector3d w_0_x = _gyr_0 - linearized_bg;
    Vector3d w_1_x = _gyr_1 - linearized_bg;
    Matrix3d R_w_0_x,R_w_1_x;

    R_w_0_x<<0, -w_0_x(2), w_0_x(1),
            w_0_x(2), 0, -w_0_x(0),
            -w_0_x(1), w_0_x(0), 0;
    R_w_1_x<<0, -w_1_x(2), w_1_x(1),
            w_1_x(2), 0, -w_1_x(0),
            -w_1_x(1), w_1_x(0), 0;

    std::vector<Eigen::Vector3d> gj(NUM_OF_LEG), gjp1(NUM_OF_LEG);
    std::vector<Eigen::Matrix3d> dgdphi_j(NUM_OF_LEG), dgdphi_jp1(NUM_OF_LEG);
    std::vector<Eigen::Matrix<double, 3, RHO_OPT_SIZE>> dgdrho_j(NUM_OF_LEG), dgdrho_jp1(NUM_OF_LEG);
    std::vector<Eigen::Vector3d> gbj(NUM_OF_LEG), gbjp1(NUM_OF_LEG);


    // velocity from individual leg
    std::vector<Eigen::Vector3d> vmj(NUM_OF_LEG), vmjp1(NUM_OF_LEG);
    // total velocity
    Eigen::Vector3d vm, vmp1;
    vm.setZero();
    vmp1.setZero();

    // from foot contact force infer a contact flag
    // calculate variance
    for (int j = 0; j < NUM_OF_LEG; j++) {
        // get z directional contact force ( contact foot sensor reading)
        double force_mag = 0.5 * (_c_0(3*j+2) + _c_1(3*j+2));

       force_mag = std::max(std::min(force_mag, 1000.0),-300.0); // limit the range of the force mag
        if (force_mag < foot_force_min[j]) {
            foot_force_min[j] = 0.9*foot_force_min[j] + 0.1*force_mag;
        }
        if (force_mag > foot_force_max[j]) {
            foot_force_max[j] = 0.9*foot_force_max[j] + 0.1*force_mag;
        }
        // exponential decay, max force decays faster
        foot_force_min[j] *= 0.9991;
        foot_force_max[j] *= 0.997;
        foot_force_contact_threshold[j] = foot_force_min[j] + V_N_FORCE_THRES_RATIO*(foot_force_max[j]-foot_force_min[j]);


        foot_contact_flag[j] = 1.0/(1+exp(-V_N_TERM1_STEEP*(force_mag-foot_force_contact_threshold[j])));

        // // get z force variance
        // foot_force_window_idx[j] ++;
        // foot_force_window_idx[j] %= FOOT_VAR_WINDOW_SIZE;
        // foot_force_window(j, foot_force_window_idx[j]) = force_mag;
        // Eigen::Matrix<double, 1, FOOT_VAR_WINDOW_SIZE> ys = foot_force_window.row(j);
        // foot_force_var[j] = (ys.array() - ys.mean()).square().sum() / (ys.size() - 1);

        // if (foot_contact_flag[j] < 0.5) {
        //     integration_contact_flag[j] = false;
        // }
    }
//    std::cout << "foot force process" << std::endl;
//    std::cout << "foot_force_max " << foot_force_max.transpose() << std::endl;
//    std::cout << "foot_force_contact_threshold " << foot_force_contact_threshold.transpose() << std::endl;
//    std::cout << "foot_force_min " << foot_force_min.transpose() << std::endl;
//    std::cout << "foot_contact_flag " << foot_contact_flag.transpose() <<std::endl;
//    std::cout << "foot_force_var " << foot_force_var.transpose() <<std::endl;

    // caution: never use leg RL because on my robot the leg is wrong
//    foot_contact_flag[2] = 0;
//    integration_contact_flag[2] = false;

//    std::cout << foot_force_var << std::endl;

    // get velocity measurement
    for (int j = 0; j < NUM_OF_LEG; j++) {
        // calculate fk of each leg
        // gj[j] = a1_kin.fk(_phi_0.segment<3>(3 * j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE * j), rho_fix_list[j]);
        // gjp1[j] = a1_kin.fk(_phi_1.segment<3>(3 * j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE * j), rho_fix_list[j]);
        // // calculate jacobian of each leg
        // dgdphi_j[j] = a1_kin.jac(_phi_0.segment<3>(3 * j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE * j), rho_fix_list[j]);
        // dgdphi_jp1[j] = a1_kin.jac(_phi_1.segment<3>(3 * j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE * j), rho_fix_list[j]);
        VectorXd test(1); test << 0.21;
        gj[j] = a1_kin.fk(_phi_0.segment<3>(3 * j), test, rho_fix_list[j]);
        gjp1[j] = a1_kin.fk(_phi_1.segment<3>(3 * j), test, rho_fix_list[j]);
        // calculate jacobian of each leg
        dgdphi_j[j] = a1_kin.jac(_phi_0.segment<3>(3 * j), test, rho_fix_list[j]);
        dgdphi_jp1[j] = a1_kin.jac(_phi_1.segment<3>(3 * j), test, rho_fix_list[j]);

        gbj[j] = p_br + R_br * gj[j];
        gbjp1[j] = p_br + R_br * gjp1[j];
        // calculate vm
        // vmj[j]   = -R_br * dgdphi_j[j]   * _dphi_0.segment<3>(3 * j) - R_w_0_x * gbj[j];
        // vmjp1[j] = -R_br * dgdphi_jp1[j] * _dphi_1.segment<3>(3 * j) - R_w_1_x * gbjp1[j];
        // Vector3d tmp = ( _phi_1.segment<3>(3 * j) - _phi_0.segment<3>(3 * j) )/_dt;
        vmj[j]   = -R_br * dgdphi_j[j]   * _dphi_0.segment<3>(3 * j) - R_w_0_x * gbj[j];
        vmjp1[j] = -R_br * dgdphi_jp1[j] * _dphi_1.segment<3>(3 * j) - R_w_1_x * gbjp1[j];

    }
    // std::cout << "debug" << std::endl;
    // std::cout << _phi_0.transpose() << std::endl;
    // std::cout << _dphi_0.transpose() << std::endl;
    // std::cout << _phi_1.transpose() << std::endl;
    // std::cout << _dphi_1.transpose() << std::endl;

    // combine velocities
    double total_weight = 0;
    double weight_j[NUM_OF_LEG];
    for (int j = 0; j < NUM_OF_LEG; j++) {
        weight_j[j] = foot_contact_flag[j];
        total_weight += weight_j[j];
        vm += weight_j[j]*vmj[j];
        vmp1 += weight_j[j]*vmjp1[j];
    }

    if (total_weight < 1.0) {
        // no foot is in contact with ground, do something to nv
        noise_diag.diagonal()[33] = 99999;
        noise_diag.diagonal()[34] = 99999;
        noise_diag.diagonal()[35] = 99999;
        vm.setZero(); 
        vmp1.setZero();  
    } else {
        vm /= total_weight; 
        vmp1 /= total_weight; 
        noise_diag.diagonal()[33] = LBV_N;
        noise_diag.diagonal()[34] = LBV_N;
        noise_diag.diagonal()[35] = LBV_N;
    }

    // propagate epsilon
    // std::cout << "vm inside mid integration" << std::endl;
    // std::cout << vmj[0]  << std::endl;
    // std::cout << weight_j[0] << weight_j[1] << weight_j[2]  << weight_j[3]  << std::endl;
    // std::cout << "total_weight " << total_weight  << std::endl;
    // std::cout << vm << std::endl;
    // std::cout << vmp1 << std::endl;
    // std::cout << delta_q.toRotationMatrix() << std::endl;
    // std::cout << result_delta_q.toRotationMatrix() << std::endl;
    Vector3d un_v_0 = delta_q * (vm - linearized_bv);
    Vector3d un_v_1 = result_delta_q * (vmp1 - linearized_bv);
    Vector3d un_v = 0.5 * (un_v_0 + un_v_1);
    result_delta_epsilon = delta_epsilon + un_v * _dt;
    // problem: delta_epsilon is very large

    result_linearized_ba = linearized_ba;
    result_linearized_bg = linearized_bg;
    result_linearized_bv = linearized_bv;
    result_linearized_rho = linearized_rho;


//     // design a new uncertainty function

//     // record all four lo velocities, examine their difference to average
//     // only choose the most accurate two
//     Matrix<double, 3, NUM_OF_LEG> lo_veocities; lo_veocities.setZero();
//     for (int j = 0; j < NUM_OF_LEG; j++) {
//         Eigen::Vector3d lo_v = 0.5 * (delta_q * vi[j] + result_delta_q * vip1[j]);
//         lo_veocities.col(j) = lo_v;
//         // base_v is the current velocity estimation in body frame
//     }

//     Vector12d uncertainties;
//     for (int j = 0; j < NUM_OF_LEG; j++) {
//         double n1 = V_N_MAX*(1-foot_contact_flag[j])+V_N_MIN;
//         double n2 = V_N_TERM2_VAR_RESCALE*foot_force_var[j];
//         Eigen::Vector3d n3; n3.setZero();
//         Eigen::Vector3d tmp = lo_veocities.col(j) - delta_v;
//         for (int k = 0; k < 3; k++) {
// //            if (fabs(tmp(k)) < 0.2) {
//                 n3(k) = V_N_TERM3_DISTANCE_RESCALE*std::pow(tmp(k),2);
// //            } else {
// //                n3(k) = 10e10;
// //            }

//         }
//         Eigen::Vector3d n = n1*Eigen::Vector3d::Ones() + n2*Eigen::Vector3d::Ones();
//         n = n + n3;
//         // we only believe
//         uncertainties.segment<3>(3*j) = n;
//     }
// //    std::cout << uncertainties.transpose() << std::endl;

//     Vector4d rho_uncertainty;
//     for (int j = 0; j < NUM_OF_LEG; j++) {
//         rho_uncertainty[j] = 5 * foot_contact_flag[j] + 0.001;
//     }

//     // use uncertainty to combine LO velocity
//     Vector3d average_delta_epsilon; average_delta_epsilon.setZero();
//     Vector3d average_count; average_count.setZero();
//     Vector12d weight_list; weight_list.setZero();



//     for (int j = 0; j < NUM_OF_LEG; j++) {
//         // large uncertainty, small weight
//         Vector3d weight = (V_N_MAX + V_N_TERM2_VAR_RESCALE + V_N_TERM3_DISTANCE_RESCALE) /  uncertainties.segment<3>(3*j).array();
//         for (int k = 0; k < 3; k++) {
//             if (weight(k) < 0.001) weight(k) = 0.001;
//         }
//         average_delta_epsilon += weight.cwiseProduct(lo_veocities.col(j)) * _dt;
//         average_count += weight;
//         weight_list.segment<3>(3*j) = weight;
//     }
// //    std::cout << weight_list.transpose() << std::endl;
// //    std::cout << "showed lists" << std::endl;
// //    std::cout << _c_0.transpose() <<std::endl;
// //    std::cout << _c_1.transpose() <<std::endl;
// //    std::cout << weight_list.transpose() <<std::endl;
// //    std::cout << uncertainties.transpose() <<std::endl;

//     for (int k = 0; k < 3; k++) {
//         average_delta_epsilon(k) /= average_count(k);
//     }


    // // abnormal case: all four feet are not on ground, in this case the residual must be all 0, we give them small uncertainty to prevent
    // if (foot_contact_flag.sum()<1e-6) {
    //     rho_uncertainty.setConstant(0.00001);
    //     uncertainties.setConstant(10e10);
    // }

    if(update_jacobian)
    {

        Eigen::Matrix3d kappa1_j[NUM_OF_LEG], kappa1_jp1[NUM_OF_LEG];
        Eigen::Matrix<double, 3, RHO_OPT_SIZE> kappa2_j[NUM_OF_LEG], kappa2_jp1[NUM_OF_LEG];
        Eigen::Matrix<double, 3, RHO_OPT_SIZE> kappa3_j[NUM_OF_LEG], kappa3_jp1[NUM_OF_LEG];
        Eigen::Matrix3d kappa4_j[NUM_OF_LEG], kappa4_jp1[NUM_OF_LEG];

        for (int j = 0; j < NUM_OF_LEG; j++) {            
            Eigen::Matrix<double, 3, 9> kron_dphi0; kron_dphi0.setZero();
            kron_dphi0(0,0) = kron_dphi0(1,1) = kron_dphi0(2,2) = _dphi_0(0+3*j);
            kron_dphi0(0,3) = kron_dphi0(1,4) = kron_dphi0(2,5) = _dphi_0(1+3*j);
            kron_dphi0(0,6) = kron_dphi0(1,7) = kron_dphi0(2,8) = _dphi_0(2+3*j);
            // calculate derivative of fk wrt rho
            dgdrho_j[j] = a1_kin.dfk_drho(_phi_0.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);
            
            Eigen::Matrix<double, 9, 3> dJdphi0 = a1_kin.dJ_dq(_phi_0.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);
            kappa1_j[j] = kron_dphi0*dJdphi0;
            
            Eigen::Matrix<double, 9, RHO_OPT_SIZE> dJdrho0 = a1_kin.dJ_drho(_phi_0.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);
            kappa2_j[j] = kron_dphi0*dJdrho0;
            kappa3_j[j] = R_br*kappa2_j[j] + R_w_0_x*R_br*dgdrho_j[j];
 
            Eigen::Matrix<double, 3, 9> kron_dphi1; kron_dphi1.setZero();
            kron_dphi1(0,0) = kron_dphi1(1,1) = kron_dphi1(2,2) = _dphi_1(0+3*j);
            kron_dphi1(0,3) = kron_dphi1(1,4) = kron_dphi1(2,5) = _dphi_1(1+3*j);
            kron_dphi1(0,6) = kron_dphi1(1,7) = kron_dphi1(2,8) = _dphi_1(2+3*j);
            dgdrho_jp1[j] = a1_kin.dfk_drho(_phi_1.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);     
            Eigen::Matrix<double, 9, 3> dJdphi1 = a1_kin.dJ_dq(_phi_1.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);
            
            kappa1_jp1[j] = kron_dphi1*dJdphi1;
            Eigen::Matrix<double, 9, RHO_OPT_SIZE> dJdrho1 = a1_kin.dJ_drho(_phi_1.segment<3>(3*j), linearized_rho.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j), rho_fix_list[j]);


            kappa2_jp1[j] = kron_dphi1*dJdrho1;
            kappa3_jp1[j] = R_br*kappa2_jp1[j] + R_w_1_x*R_br*dgdrho_jp1[j];
        }

        // a very complicated term
        Eigen::Matrix3d kappa_10, kappa_11;
        kappa_11 = 0.25 * _dt * _dt * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vmp1 - linearized_bv);
        kappa_10 = 2.0*kappa_11;
        if (total_weight > 0) {
            for (int j = 0; j < NUM_OF_LEG; j++) {   
                kappa_10 = kappa_10 - 0.5 * _dt *(
                    delta_q.toRotationMatrix()*weight_j[j]/total_weight * Utility::skewSymmetric(gbj[j]) +
                    result_delta_q.toRotationMatrix()*weight_j[j]/total_weight * Utility::skewSymmetric(gbjp1[j])
                );
            }
        }


        Vector3d w_x = 0.5 * (_gyr_0 + _gyr_1) - linearized_bg;
        Vector3d a_0_x = _acc_0 - linearized_ba;
        Vector3d a_1_x = _acc_1 - linearized_ba;
        Matrix3d R_w_x, R_a_0_x, R_a_1_x;

        R_w_x<<0, -w_x(2), w_x(1),
                w_x(2), 0, -w_x(0),
                -w_x(1), w_x(0), 0;
        R_a_0_x<<0, -a_0_x(2), a_0_x(1),
                a_0_x(2), 0, -a_0_x(0),
                -a_0_x(1), a_0_x(0), 0;
        R_a_1_x<<0, -a_1_x(2), a_1_x(1),
                a_1_x(2), 0, -a_1_x(0),
                -a_1_x(1), a_1_x(0), 0;
        Eigen::Matrix3d kappa_9 = (Matrix3d::Identity() - R_w_x * _dt);
        Eigen::Matrix3d kappa_5 = -0.5 * _dt * delta_q.toRotationMatrix() * R_a_0_x +
                                  -0.5 * _dt * result_delta_q.toRotationMatrix() * R_a_1_x * kappa_9;
        Eigen::Matrix3d kappa_7 = 0.25 * _dt * _dt * result_delta_q.toRotationMatrix() * R_a_1_x;
        Eigen::Matrix3d kappa_6 = 2.0*kappa_7;

        Eigen::Matrix3d kappa_8 = -0.5 * _dt * delta_q.toRotationMatrix() * Utility::skewSymmetric(vm - linearized_bv)
                                  -0.5 * _dt * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vmp1 - linearized_bv)*kappa_9;

        Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> F; F.setZero();

        // F row 1
        F.block<3, 3>(ILO_P, ILO_P) = Matrix3d::Identity();
        F.block<3, 3>(ILO_P, ILO_R) = 0.5 * _dt * kappa_5;
        F.block<3, 3>(ILO_P, ILO_V) = Matrix3d::Identity() * _dt;
        F.block<3, 3>(ILO_P, ILO_BA) = -0.25 * _dt * _dt * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix());

        // F row 2
        F.block<3, 3>(ILO_R, ILO_R) = kappa_9;
        F.block<3, 3>(ILO_R, ILO_BG) = -1.0 * Matrix3d::Identity() * _dt;

        // F row 3
        F.block<3, 3>(ILO_V, ILO_R) = kappa_5;
        F.block<3, 3>(ILO_V, ILO_V) = Matrix3d::Identity();
        F.block<3, 3>(ILO_V, ILO_BA) = -0.5 * _dt * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix());
        F.block<3, 3>(ILO_V, ILO_BG) = kappa_6;

        // F row 4 
        F.block<3, 3>(ILO_EPS, ILO_R)   = kappa_8;
        F.block<3, 3>(ILO_EPS, ILO_EPS) = Matrix3d::Identity();
        F.block<3, 3>(ILO_EPS, ILO_BG)  = kappa_10;
        F.block<3, 3>(ILO_EPS, ILO_BV)  = -0.5 * _dt * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix());
        if (total_weight > 0) {
            for (int j = 0; j < NUM_OF_LEG; j++) {   
                F.block<3, RHO_OPT_SIZE>(ILO_EPS, ILO_RHO1+j) = -0.5 * _dt * (
                    delta_q.toRotationMatrix() * weight_j[j]/total_weight*kappa3_j[j] + 
                    result_delta_q.toRotationMatrix() * weight_j[j]/total_weight*kappa3_jp1[j]
                );
            }
        }

        // F rest rows
        F.block<3, 3>(ILO_BA, ILO_BA) = Matrix3d::Identity();
        F.block<3, 3>(ILO_BG, ILO_BG) = Matrix3d::Identity();
        F.block<3, 3>(ILO_BV, ILO_BV) = Matrix3d::Identity();
        F.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO1, ILO_RHO1) = Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity();
        F.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO2, ILO_RHO2) = Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity();
        F.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO3, ILO_RHO3) = Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity();
        F.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO4, ILO_RHO4) = Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity();

        // get V

        Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> V; V.setZero();

        V.block<3, 3>(ILO_P, ILNO_Ai) =  - 0.25 * _dt * _dt * delta_q.toRotationMatrix();
        V.block<3, 3>(ILO_P, ILNO_Gi)  =  0.5 * _dt * kappa_7;
        V.block<3, 3>(ILO_P, ILNO_Ai1) =  - 0.25 * _dt * _dt * result_delta_q.toRotationMatrix();
        V.block<3, 3>(ILO_P, ILNO_Gi1) =  V.block<3, 3>(ILO_P, ILNO_Gi);

        V.block<3, 3>(ILO_R, ILNO_Gi)  = -0.5 * Matrix3d::Identity() * _dt;
        V.block<3, 3>(ILO_R, ILNO_Gi1) = -0.5 * Matrix3d::Identity() * _dt;

        V.block<3, 3>(ILO_V, ILNO_Ai)  =  -0.5 * _dt * delta_q.toRotationMatrix();
        V.block<3, 3>(ILO_V, ILNO_Gi)  =  kappa_7;
        V.block<3, 3>(ILO_V, ILNO_Ai1) =  -0.5 * _dt * result_delta_q.toRotationMatrix();
        V.block<3, 3>(ILO_V, ILNO_Gi1) =  V.block<3, 3>(ILO_V, ILNO_Gi);

        // V row 4 
        V.block<3, 3>(ILO_EPS, ILNO_Gi) = kappa_11;
        V.block<3, 3>(ILO_EPS, ILNO_Gi1) = kappa_11;
        V.block<3, 3>(ILO_EPS, ILNO_PHIi).setZero();
        V.block<3, 3>(ILO_EPS, ILNO_PHIi1).setZero();
        V.block<3, 3>(ILO_EPS, ILNO_DPHIi).setZero();
        V.block<3, 3>(ILO_EPS, ILNO_DPHIi1).setZero();
        if (total_weight > 0) {
            for (int j = 0; j < NUM_OF_LEG; j++) {
                V.block<3, 3>(ILO_EPS, ILNO_Gi)  += - 0.5*_dt*delta_q.toRotationMatrix()*weight_j[j]/total_weight* Utility::skewSymmetric(gbj[j]);
                V.block<3, 3>(ILO_EPS, ILNO_Gi1) += - 0.5*_dt*result_delta_q.toRotationMatrix()*weight_j[j]/total_weight* Utility::skewSymmetric(gbjp1[j]);
                V.block<3, 3>(ILO_EPS, ILNO_PHIi) += 0.5*_dt*delta_q.toRotationMatrix()*weight_j[j]/total_weight*kappa4_j[j];
                V.block<3, 3>(ILO_EPS, ILNO_PHIi1) += 0.5*_dt*result_delta_q.toRotationMatrix()*weight_j[j]/total_weight*kappa4_jp1[j];
                V.block<3, 3>(ILO_EPS, ILNO_DPHIi) += 0.5*_dt*delta_q.toRotationMatrix()*R_br* weight_j[j]/total_weight*dgdphi_j[j];
                V.block<3, 3>(ILO_EPS, ILNO_DPHIi1) += 0.5*_dt*result_delta_q.toRotationMatrix()*R_br* weight_j[j]/total_weight*dgdphi_jp1[j];
            }
        }
        V.block<3, 3>(ILO_EPS, ILNO_V) = - Matrix3d::Identity() * _dt;

        // V rest rows
        V.block<3, 3>(ILO_BA, ILNO_BA) = Matrix3d::Identity() * _dt;
        V.block<3, 3>(ILO_BG, ILNO_BG) = Matrix3d::Identity() * _dt;
        V.block<3, 3>(ILO_BV, ILNO_BV) = Matrix3d::Identity() * _dt;

        V.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO1, ILNO_NRHO1) = -Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity() * _dt;
        V.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO2, ILNO_NRHO2) = -Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity() * _dt;
        V.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO3, ILNO_NRHO3) = -Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity() * _dt;
        V.block<RHO_OPT_SIZE, RHO_OPT_SIZE>(ILO_RHO4, ILNO_NRHO4) = -Eigen::Matrix<double, RHO_OPT_SIZE, RHO_OPT_SIZE>::Identity() * _dt;


        jacobian = F * jacobian;
        covariance = F * covariance * F.transpose() + V * noise_diag * V.transpose();
        step_jacobian = F;
        step_V = V;
    }
//    std::cout << "foot_force" << foot_force.transpose() << std::endl;
//    std::cout << "foot_force_min" << foot_force_min.transpose() << std::endl;
//    std::cout << "foot_force_max" << foot_force_max.transpose() << std::endl;
//    std::cout << "foot_contact_flag" << foot_contact_flag.transpose() << std::endl;
//    std::cout << "foot_force_contact_threshold" << foot_force_contact_threshold.transpose() << std::endl;
//    std::cout << "noise_diag" << noise_diag.diagonal().segment<12>(30).transpose() << std::endl;
}

void IMULegIntegrationBase::checkJacobian(double _dt, const Vector3d &_acc_0, const Vector3d &_gyr_0,
                  const Vector3d &_acc_1, const Vector3d &_gyr_1,
                  const Ref<const Vector12d> &_phi_0, const Ref<const Vector12d> &_dphi_0,
                  const Ref<const Vector12d> &_c_0, const Ref<const Vector12d> &_phi_1,
                  const Ref<const Vector12d> &_dphi_1, const Ref<const Vector12d> &_c_1,
                  const Vector3d &delta_p, const Quaterniond &delta_q,
                  const Vector3d &delta_v, const Vector3d delta_epsilon, 
                  const Vector3d &linearized_ba, const Vector3d &linearized_bg, const Vector3d &linearized_bv,
                  const Ref<const Vector_rho> &linearized_rho) {
    Vector3d result_delta_p;
    Quaterniond result_delta_q;
    Vector3d result_delta_v;
    Vector3d result_delta_epsilon;
    Vector3d result_linearized_ba;
    Vector3d result_linearized_bg;
    Vector3d result_linearized_bv;
    Vector_rho result_linearized_rho;
    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        result_delta_p, result_delta_q, result_delta_v, result_delta_epsilon,
                        result_linearized_ba, result_linearized_bg,result_linearized_bv, result_linearized_rho, 1);

    Vector3d turb_delta_p;
    Quaterniond turb_delta_q;
    Vector3d turb_delta_v;
    Vector3d turb_delta_epsilon;
    Vector3d turb_linearized_ba;
    Vector3d turb_linearized_bg;
    Vector3d turb_linearized_bv;
    Vector_rho turb_linearized_rho;

    Vector3d turb(0.0001, -0.003, 0.003);

    cout << "------------------- check jacobian --------------------" << endl;
    cout << "------------------- check jacobian --------------------" << endl;
    cout << "------------------- check jacobian --------------------" << endl;
    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p + turb, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb p-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_P) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_P) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_P) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_P) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_P) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_P) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_P) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_P) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q * Quaterniond(1, turb(0) / 2, turb(1) / 2, turb(2) / 2), delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb q-------F col 2--------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_R) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_R) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_R) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_R) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_R) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_R) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_R) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_R) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v + turb, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb v-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_V) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_V) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_V) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_V) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_V) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_V) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_V) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_V) * turb).transpose() << endl;


    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon + turb, 
                        linearized_ba, linearized_bg, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb eps-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_EPS) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_EPS) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_EPS) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_EPS) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_EPS) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_EPS) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_EPS) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_EPS) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba + turb, linearized_bg, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb ba-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BA) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BA) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BA) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BA) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BA) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BA) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BA) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_BA) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg + turb, linearized_bv, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb bg-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BG) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BG) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BG) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BG) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BG) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BG) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BG) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_BG) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv + turb, linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb bv-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BV) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BV) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BV) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BV) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BV) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BV) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BV) * turb).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, 3>(ILO_RHO1, ILO_BV) * turb).transpose() << endl;


    Vector_rho input_turb_linearized_rho = linearized_rho;
    Eigen::Matrix<double, TOTAL_RHO_OPT_SIZE, 1> turb_rho;
    turb_rho << 0.003, 0.001, -0.0008, 0.001;
    input_turb_linearized_rho = linearized_rho + turb_rho;
    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, delta_epsilon, 
                        linearized_ba, linearized_bg, linearized_bv, input_turb_linearized_rho,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg,turb_linearized_bv, turb_linearized_rho, 0);
    cout << "turb rho-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_P, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_R, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_V, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "eps diff       " << (turb_delta_epsilon - result_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_EPS, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_BA, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_BG, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, TOTAL_RHO_OPT_SIZE>(ILO_BV, ILO_RHO1) * turb_rho).transpose() << endl;
    cout << "rho diff " << (turb_linearized_rho - result_linearized_rho).transpose() << endl;
    cout << "rho jacob diff " << (step_jacobian.block<TOTAL_RHO_OPT_SIZE, TOTAL_RHO_OPT_SIZE>(ILO_RHO1, ILO_RHO1) * turb_rho).transpose() << endl;
}

Eigen::Matrix<double, RESIDUAL_STATE_SIZE, 1>
IMULegIntegrationBase::evaluate(const Vector3d &Pi, const Quaterniond &Qi, const Vector3d &Vi, 
                                const Vector3d &Bai, const Vector3d &Bgi, const Vector3d &Bvi, const Vector_rho &rhoi, 
                                const Vector3d &Pj, const Quaterniond &Qj, const Vector3d &Vj, 
                                const Vector3d &Baj, const Vector3d &Bgj, const Vector3d &Bvj, const Vector_rho &rhoj) {
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, 1> residuals;


    Eigen::Matrix3d dp_dba = jacobian.block<3, 3>(ILO_P, ILO_BA);
    Eigen::Matrix3d dp_dbg = jacobian.block<3, 3>(ILO_P, ILO_BG);

    Eigen::Matrix3d dq_dbg = jacobian.block<3, 3>(ILO_R, ILO_BG);

    Eigen::Matrix3d dv_dba = jacobian.block<3, 3>(ILO_V, ILO_BA);
    Eigen::Matrix3d dv_dbg = jacobian.block<3, 3>(ILO_V, ILO_BG);

    Eigen::Matrix3d deps_dbg = jacobian.block<3, 3>(ILO_EPS, ILO_BG);
    Eigen::Matrix3d deps_dbv = jacobian.block<3, 3>(ILO_EPS, ILO_BV);
    Eigen::Matrix<double, 3, RHO_OPT_SIZE> deps_drho1 = jacobian.block<3, RHO_OPT_SIZE>(ILO_EPS, ILO_RHO1);
    Eigen::Matrix<double, 3, RHO_OPT_SIZE> deps_drho2 = jacobian.block<3, RHO_OPT_SIZE>(ILO_EPS, ILO_RHO2);
    Eigen::Matrix<double, 3, RHO_OPT_SIZE> deps_drho3 = jacobian.block<3, RHO_OPT_SIZE>(ILO_EPS, ILO_RHO3);
    Eigen::Matrix<double, 3, RHO_OPT_SIZE> deps_drho4 = jacobian.block<3, RHO_OPT_SIZE>(ILO_EPS, ILO_RHO4);

    Eigen::Vector3d dba = Bai - linearized_ba; // Bai is the new optimization result
    Eigen::Vector3d dbg = Bgi - linearized_bg;
    Eigen::Vector3d dbv = Bvi - linearized_bv;

    Vector_rho drho = rhoi  - linearized_rho;

    Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * dbg);
    Eigen::Vector3d corrected_delta_v = delta_v + dv_dba * dba + dv_dbg * dbg;
    Eigen::Vector3d corrected_delta_p = delta_p + dp_dba * dba + dp_dbg * dbg;

    Eigen::Vector3d corrected_delta_epsilon = delta_epsilon + deps_dbg * dbg + deps_dbv * dbv; 
    corrected_delta_epsilon += deps_drho1 * drho.segment<RHO_OPT_SIZE>(0*RHO_OPT_SIZE);
    corrected_delta_epsilon += deps_drho2 * drho.segment<RHO_OPT_SIZE>(1*RHO_OPT_SIZE);
    corrected_delta_epsilon += deps_drho3 * drho.segment<RHO_OPT_SIZE>(2*RHO_OPT_SIZE);
    corrected_delta_epsilon += deps_drho4 * drho.segment<RHO_OPT_SIZE>(3*RHO_OPT_SIZE);

    residuals.block<3, 1>(ILO_P, 0) = Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt) - corrected_delta_p;
    residuals.block<3, 1>(ILO_R, 0) = 2 * (corrected_delta_q.inverse() * (Qi.inverse() * Qj)).vec();
    residuals.block<3, 1>(ILO_V, 0) = Qi.inverse() * (G * sum_dt + Vj - Vi) - corrected_delta_v;
    
    residuals.block<3, 1>(ILO_EPS, 0) = Qi.inverse() * (Pj - Pi) - corrected_delta_epsilon;

    for (int j = 0; j < NUM_OF_LEG; j++) {
            residuals.block<RHO_OPT_SIZE, 1>(ILO_RHO1+RHO_OPT_SIZE*j, 0) = rhoj.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j) - rhoi.segment<RHO_OPT_SIZE>(RHO_OPT_SIZE*j);
    }
    residuals.block<3, 1>(ILO_BA, 0) = Baj - Bai;
    residuals.block<3, 1>(ILO_BG, 0) = Bgj - Bgi;
    residuals.block<3, 1>(ILO_BV, 0) = Bvj - Bvi;

    return residuals;
}
