//
// Created by shuoy on 8/23/21.
//

#include "imu_leg_integration_base.h"

IMULegIntegrationBase::IMULegIntegrationBase(const Vector3d &_base_v, const Vector3d &_acc_0, const Vector3d &_gyr_0, const Ref<const Vector12d>& _phi_0,
                                             const Ref<const Vector12d>& _dphi_0, const Ref<const Vector12d>& _c_0,
                                             const Vector3d &_linearized_ba, const Vector3d &_linearized_bg,
                                             const Vector3d &_linearized_bv,
                                             std::vector<Eigen::VectorXd> _rho_fix_list, const Eigen::Vector3d &_p_br,  const Eigen::Matrix3d &_R_br)
        : acc_0{_acc_0}, gyr_0{_gyr_0}, linearized_acc{_acc_0}, linearized_gyr{_gyr_0},
          linearized_ba{_linearized_ba}, linearized_bg{_linearized_bg}, linearized_bv{_linearized_bv},
          sum_dt{0.0}, delta_p{Eigen::Vector3d::Zero()}, delta_q{Eigen::Quaterniond::Identity()}, delta_v{Eigen::Vector3d::Zero()},
          sum_delta_epsilon{Eigen::Vector3d::Zero()}
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

    foot_force_min.setZero();
    foot_force_max.setZero();

    // the fixed kinematics parameter
    rho_fix_list = _rho_fix_list;
    p_br = _p_br;
    R_br = _R_br;
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
void IMULegIntegrationBase::repropagate(const Eigen::Vector3d &_linearized_ba,
                                        const Eigen::Vector3d &_linearized_bg,
                                        const Eigen::Vector3d &_linearized_bv)
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
    sum_delta_epsilon.setZero();
    linearized_ba = _linearized_ba;
    linearized_bg = _linearized_bg;
    linearized_bv = _linearized_bv;
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
    Vector3d result_sum_delta_epsilon;
    Vector3d result_linearized_ba;
    Vector3d result_linearized_bg;
    Vector3d result_linearized_bv;
    // midPointIntegration
    midPointIntegration(_dt, acc_0, gyr_0, acc_1, gyr_1,
                        phi_0, dphi_0, c_0, phi_1, dphi_1, c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg, linearized_bv,
                        result_delta_p, result_delta_q, result_delta_v, result_sum_delta_epsilon,
                        result_linearized_ba, result_linearized_bg, result_linearized_bv, 1);
    // checkJacobian
//    checkJacobian(_dt, acc_0, gyr_0, acc_1, gyr_1,
//                        phi_0, dphi_0, c_0, phi_1, dphi_1, c_1,
//                        delta_p, delta_q, delta_v, sum_delta_epsilon,
//                        linearized_ba, linearized_bg, linearized_bv);

    delta_p = result_delta_p;
    delta_q = result_delta_q;
    delta_v = result_delta_v;
    sum_delta_epsilon = result_sum_delta_epsilon;
    linearized_ba = result_linearized_ba;
    linearized_bg = result_linearized_bg;
    linearized_bv = result_linearized_bv;
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
                                                const Vector3d &delta_v, const Vector3d &sum_delta_epsilon,
                                                const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                                                const Vector3d &linearized_bv, Vector3d &result_delta_p,
                                                Quaterniond &result_delta_q, Vector3d &result_delta_v, Vector3d &result_sum_delta_epsilon,
                                                Vector3d &result_linearized_ba, Vector3d &result_linearized_bg,
                                                Vector3d &result_linearized_bv, bool update_jacobian) {
    Vector3d un_acc_0 = delta_q * (_acc_0 - linearized_ba);
    Vector3d un_gyr = 0.5 * (_gyr_0 + _gyr_1) - linearized_bg;
    result_delta_q = delta_q * Quaterniond(1, un_gyr(0) * _dt / 2, un_gyr(1) * _dt / 2, un_gyr(2) * _dt / 2);
    Vector3d un_acc_1 = result_delta_q * (_acc_1 - linearized_ba);
    Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    result_delta_p = delta_p + delta_v * _dt + 0.5 * un_acc * _dt * _dt;
    result_delta_v = delta_v + un_acc * _dt;
    result_linearized_ba = linearized_ba;
    result_linearized_bg = linearized_bg;
    result_linearized_bv = linearized_bv;

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

    std::vector<Eigen::Vector3d> fi, fip1;
    std::vector<Eigen::Matrix3d> Ji, Jip1;
    std::vector<Eigen::Matrix3d> hi, hip1;
    std::vector<Eigen::Vector3d> vi, vip1;   //velocity of each leg
    Vector3d vmi, vmi1;                      // weighted sum of all legs

    // from foot contact force infer a contact flag
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
        if ( force_mag > foot_force_contact_threshold[j]) {
            foot_contact_flag[j] = 1;
        } else {
            foot_contact_flag[j] = 0;
        }
    }
    // get velocity measurement
    for (int j = 0; j < NUM_OF_LEG; j++) {
        // calculate fk of each leg
        fi.push_back(a1_kin.fk(_phi_0.segment<3>(3 * j), Eigen::Vector3d::Zero(), rho_fix_list[j]));
        fip1.push_back(a1_kin.fk(_phi_1.segment<3>(3 * j), Eigen::Vector3d::Zero(), rho_fix_list[j]));
        // calculate jacobian of each leg
        Ji.push_back(a1_kin.jac(_phi_0.segment<3>(3 * j), Eigen::Vector3d::Zero(), rho_fix_list[j]));
        Jip1.push_back(a1_kin.jac(_phi_1.segment<3>(3 * j), Eigen::Vector3d::Zero(), rho_fix_list[j]));

        // calculate vm
        vi.push_back(-R_br * Ji[j] * _dphi_0.segment<3>(3 * j) - R_w_0_x * (p_br + R_br * fi[j]));
        vip1.push_back(-R_br * Jip1[j] * _dphi_1.segment<3>(3 * j) - R_w_1_x * (p_br + R_br * fip1[j]));

    }

    // design a new uncertainty function
    Vector12d uncertainties;
    for (int j = 0; j < NUM_OF_LEG; j++) {
        double force_mag = 0.5 * (_c_0(3*j+2) + _c_1(3*j+2));
        double diff_c = (_c_1(3*j+2) - _c_0(3*j+2)) / _dt;
        double diff_c_mag = diff_c;
        // term 1
        double n1 = V_N_MAX*(1-1/(1+exp(-V_N_TERM1_STEEP*(force_mag-foot_force_contact_threshold[j]))))+V_N_MIN;

        // term 2
        if (force_mag>foot_force_contact_threshold[j]) {
            diff_c_mag /= V_N_TERM2_BOUND_FORCE/_dt;
            if (diff_c_mag > 0)
                diff_c_mag = V_N_TERM2_VAR_RESCALE*diff_c_mag;
        } else {
            diff_c_mag = 0;
        }
        double n2 = std::fmax(0, sqrt(fabs(diff_c_mag))*V_N_MAX - 80)+V_N_MIN; //relu

        // term 3
        Eigen::Vector3d n3 = V_N_MIN*Eigen::Vector3d::Ones();
        if (force_mag>foot_force_contact_threshold[j]) {
            // velocity of each individual leg
            Eigen::Vector3d lo_v = 0.5 * (delta_q * (vi[j]-linearized_bv) + result_delta_q * (vip1[j]-linearized_bv));
            Eigen::Vector3d diff;
            diff.setZero();
            // prevent abnormal base_v
            if (base_v.norm() / 3 < 5) {
                diff = lo_v - base_v;
                n3 = diff.array().square().sqrt();
                n3 -= Eigen::Vector3d(V_N_TERM3_VEL_DIFF_XY * fabs(base_v(0)),
                                      V_N_TERM3_VEL_DIFF_XY * fabs(base_v(1)),
                                      V_N_TERM3_VEL_DIFF_Z * fabs(base_v(2)));
                n3 = n3.cwiseMax(0);//relu
                n3 = V_N_TERM3_DISTANCE_RESCALE * n3 + V_N_MIN * Eigen::Vector3d::Ones();
            }
        }

        Eigen::Vector3d n = V_N_FINAL_RATIO*(n1*Eigen::Vector3d::Ones() +
                n2*Eigen::Vector3d::Ones() + n3
                );
        uncertainties.segment<3>(3*j) = n;
    }
//    std::cout << uncertainties.transpose() << std::endl;

    // use uncertainty to combine LO velocity
    Vector4d w;  // weight of each leg
    double total_weight = 0; // total weight of all legs
    vmi.setZero();
    vmi1.setZero();
    for (int j = 0; j < NUM_OF_LEG; j++) {
        w(j) = 1.0/uncertainties.segment<3>(3*j).norm();
        total_weight += w(j);
        vmi += w(j) * vi[j];
        vmi1 += w(j) * vip1[j];
    }

    vmi /= total_weight;
    vmi1 /= total_weight;

    Eigen::Vector3d average_delta_epsilon = 0.5 * (delta_q.toRotationMatrix() * (vmi-linearized_bv) + result_delta_q.toRotationMatrix() * (vmi1-linearized_bv))* _dt;
    result_sum_delta_epsilon = sum_delta_epsilon + average_delta_epsilon;



    noise_diag.diagonal() <<
            (ACC_N * ACC_N), (ACC_N * ACC_N), (ACC_N * ACC_N),
            (GYR_N * GYR_N), (GYR_N * GYR_N), (GYR_N * GYR_N),
            (ACC_N * ACC_N), (ACC_N * ACC_N), (ACC_N * ACC_N),
            (GYR_N * GYR_N), (GYR_N * GYR_N), (GYR_N * GYR_N),
            (ACC_W * ACC_W), (ACC_W * ACC_W), (ACC_W * ACC_W),
            (GYR_W * GYR_W), (GYR_W * GYR_W), (GYR_W * GYR_W),
            (PHI_N * PHI_N), (PHI_N * PHI_N), (PHI_N * PHI_N),
            (PHI_N * PHI_N), (PHI_N * PHI_N), (PHI_N * PHI_N),
            (DPHI_N * DPHI_N), (DPHI_N * DPHI_N), (DPHI_N * DPHI_N),
            (DPHI_N * DPHI_N), (DPHI_N * DPHI_N), (DPHI_N * DPHI_N),
            (RHO_XY_N * RHO_XY_N), (RHO_XY_N * RHO_XY_N), (RHO_XY_N * RHO_XY_N),
            (RHO_XY_N * RHO_XY_N), (RHO_XY_N * RHO_XY_N), (RHO_XY_N * RHO_XY_N),
            (RHO_Z_N * RHO_Z_N), (RHO_Z_N * RHO_Z_N), (RHO_Z_N * RHO_Z_N);
    if(update_jacobian)
    {
        for (int j = 0; j < NUM_OF_LEG; j++) {

            Eigen::Matrix<double, 3, 9> kron_dphi0; kron_dphi0.setZero();
            kron_dphi0(0,0) = kron_dphi0(1,1) = kron_dphi0(2,2) = _dphi_0(0+3*j);
            kron_dphi0(0,3) = kron_dphi0(1,4) = kron_dphi0(2,5) = _dphi_0(1+3*j);
            kron_dphi0(0,6) = kron_dphi0(1,7) = kron_dphi0(2,8) = _dphi_0(2+3*j);

            Eigen::Matrix<double, 3, 9> kron_dphi1; kron_dphi1.setZero();
            kron_dphi1(0,0) = kron_dphi1(1,1) = kron_dphi1(2,2) = _dphi_1(0+3*j);
            kron_dphi1(0,3) = kron_dphi1(1,4) = kron_dphi1(2,5) = _dphi_1(1+3*j);
            kron_dphi1(0,6) = kron_dphi1(1,7) = kron_dphi1(2,8) = _dphi_1(2+3*j);

            // calculate h
            Eigen::Matrix<double, 9, 3> dJdphi0 = a1_kin.dJ_dq(_phi_0.segment<3>(3*j), Eigen::Vector3d::Zero(), rho_fix_list[j]);
            hi.push_back( (R_br*kron_dphi0*dJdphi0 + R_w_0_x*R_br*Ji[j]) );
            Eigen::Matrix<double, 9, 3> dJdphi1 = a1_kin.dJ_dq(_phi_1.segment<3>(3*j), Eigen::Vector3d::Zero(), rho_fix_list[j]);
            hip1.push_back((R_br*kron_dphi1*dJdphi1 + R_w_1_x*R_br*Jip1[j]) );
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
        Eigen::Matrix3d kappa_7 = (Matrix3d::Identity() - R_w_x * _dt);
        Eigen::Matrix3d kappa_1 = -0.5 * delta_q.toRotationMatrix() * R_a_0_x * _dt +
                                  -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * kappa_7 * _dt;
        Eigen::Matrix3d kappa_8 = -0.25 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt * _dt;
        Eigen::Matrix3d kappa_2 = 0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * _dt * _dt;
        Eigen::Matrix3d kappa_3 = 0.5 * kappa_2;
        Eigen::Matrix3d kappa_9 = -0.5 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt;

        Eigen::Matrix3d kappa_4 = -0.5* delta_q.toRotationMatrix() * Utility::skewSymmetric(vmi-linearized_bv) * _dt +
                                  -0.5 * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vmi1-linearized_bv) * kappa_7 * _dt;

        Eigen::Matrix3d w_sum_skew_fb;    w_sum_skew_fb.setZero();
        Eigen::Matrix3d w_sum_skew_fbp1;  w_sum_skew_fbp1.setZero();
        Eigen::Matrix3d w_sum_hi;         w_sum_hi.setZero();
        Eigen::Matrix3d w_sum_hip1;       w_sum_hip1.setZero();
        Eigen::Matrix3d w_sum_Ji;         w_sum_Ji.setZero();
        Eigen::Matrix3d w_sum_Jip1;       w_sum_Jip1.setZero();
        for (int j = 0; j < NUM_OF_LEG; j++) {
            w_sum_skew_fb   += w[j] * Utility::skewSymmetric(p_br + R_br*fi[j]);
            w_sum_skew_fbp1 += w[j] * Utility::skewSymmetric(p_br + R_br*fip1[j]);
            w_sum_hi        += w[j] * hi[j];
            w_sum_hip1      += w[j] * hip1[j];
            w_sum_Ji        += w[j] * Ji[j];
            w_sum_Jip1      += w[j] * Jip1[j];
        }
        w_sum_skew_fb /= total_weight;
        w_sum_skew_fbp1 /= total_weight;
        w_sum_hi /= total_weight;
        w_sum_hip1 /= total_weight;
        w_sum_Ji /= total_weight;
        w_sum_Jip1 /= total_weight;
        Eigen::Matrix3d kappa_5 = 0.5 * _dt * _dt * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vmi1-linearized_bv)
                                  - 0.5 * _dt * (delta_q.toRotationMatrix()*w_sum_skew_fb + result_delta_q.toRotationMatrix()*w_sum_skew_fbp1);


        Eigen::Matrix3d kappa_6 = 0.25* result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vmi-linearized_bv) * _dt * _dt;

        // change to sparse matrix later otherwise they are too large
        Eigen::Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE> F; F.setZero();
//        Eigen::SparseMatrix<double> F(RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE);
//        std::vector<Trip> trp;
//        Eigen::Matrix3d tmp33;
        // F row 1
        F.block<3, 3>(ILO_P, ILO_P) = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_P+s_i,ILO_P+s_i,1));
//        tmp33 = 0.5 * _dt * kappa_1;
        F.block<3, 3>(ILO_P, ILO_R) = 0.5 * _dt * kappa_1;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_P+s_i,ILO_R+s_j,tmp33(s_i,s_j)));
        F.block<3, 3>(ILO_P, ILO_V) = Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_P+s_i,ILO_V+s_i,_dt));
        F.block<3, 3>(ILO_P, ILO_BA) = kappa_8;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_P+s_i,ILO_BA+s_j,kappa_8(s_i,s_j)));
        F.block<3, 3>(ILO_P, ILO_BG) = 0.5 * _dt * kappa_2;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_P+s_i,ILO_BG+s_j,0.5*_dt*kappa_2(s_i,s_j)));

//        // F row 2
        F.block<3, 3>(ILO_R, ILO_R) = kappa_7;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_R+s_i,ILO_R+s_j,kappa_7(s_i,s_j)));
        F.block<3, 3>(ILO_R, ILO_BG) = -1.0 * Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_R+s_i,ILO_BG+s_i,-_dt));
//        // F row 3
        F.block<3, 3>(ILO_V, ILO_R) = kappa_1;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_V+s_i,ILO_R+s_j,kappa_1(s_i,s_j)));
        F.block<3, 3>(ILO_V, ILO_V) = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_V+s_i,ILO_V+s_i,1));
        F.block<3, 3>(ILO_V, ILO_BA) = kappa_9;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_V+s_i,ILO_BA+s_j,kappa_9(s_i,s_j)));

        F.block<3, 3>(ILO_V, ILO_BG) = kappa_2;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_V+s_i,ILO_BG+s_j,kappa_2(s_i,s_j)));
//
//        // F row 4
        F.block<3, 3>(ILO_EPS, ILO_R) = kappa_4;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_EPS+s_i,ILO_R+s_j,kappa_4(s_i,s_j)));

            F.block<3, 3>(ILO_EPS, ILO_EPS)  = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_EPS+s_i,ILO_EPS+s_i,1));

        F.block<3, 3>(ILO_EPS, ILO_BG) = kappa_5;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_EPS+s_i,ILO_BG+s_j,kappa_5(s_i,s_j)));

        F.block<3, 3>(ILO_EPS, ILO_BV) = kappa_9;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(ILO_EPS+s_i,ILO_BV+s_j,kappa_9(s_i,s_j)));
//        // F row 5 6 7

        F.block<3, 3>(ILO_BA, ILO_BA) = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_BA+s_i,ILO_BA+s_i,1));
        F.block<3, 3>(ILO_BG, ILO_BG) = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_BG+s_i,ILO_BG+s_i,1));
        F.block<3, 3>(ILO_BV, ILO_BV) = Matrix3d::Identity();
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(ILO_BV+s_i,ILO_BV+s_i,1));

//        F.setFromTriplets(trp.begin(),trp.end());

        // get V
        Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> V; V.setZero();
//        Eigen::SparseMatrix<double> V(RESIDUAL_STATE_SIZE, NOISE_SIZE);
//        trp.clear();
//        trp.resize(0);
        // V row 1
        V.block<3, 3>(ILO_P, ILNOISE_AI) =  -0.25 * delta_q.toRotationMatrix() * _dt * _dt;
//        tmp33 =  0.25 * delta_q.toRotationMatrix() * _dt * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(0+s_i,0+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_P, ILNOISE_GI) =  kappa_3 * 0.5 * _dt;
//        tmp33 =  0.25 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * _dt * 0.5 * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(0+s_i,3+s_j,tmp33(s_i,s_j)));
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(0+s_i,9+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_P, ILNOISE_AI1) =  -0.25 * result_delta_q.toRotationMatrix() * _dt * _dt;
//        tmp33 =  0.25 * result_delta_q.toRotationMatrix() * _dt * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(0+s_i,6+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_P, ILNOISE_GI1) =  kappa_3 * 0.5 * _dt;

        // V row 2
        V.block<3, 3>(ILO_R, ILNOISE_GI) =  -0.5 * Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(3+s_i,3+s_i,0.5*_dt));
        V.block<3, 3>(ILO_R, ILNOISE_GI1) =  -0.5 * Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(3+s_i,9+s_i,0.5*_dt));

        // V row 3
        V.block<3, 3>(ILO_V, ILNOISE_AI) =  -0.5 * delta_q.toRotationMatrix() * _dt;
//        tmp33 = 0.5 * delta_q.toRotationMatrix() * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(6+s_i,0+s_j,tmp33(s_i,s_j)));

        V.block<3, 3>(ILO_V, ILNOISE_GI) =  kappa_3;
//        tmp33 = 0.5 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * 0.5 * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(6+s_i,3+s_j,tmp33(s_i,s_j)));

//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(6+s_i,9+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_V, ILNOISE_AI1) =  -0.5 * result_delta_q.toRotationMatrix() * _dt;
//        tmp33 = 0.5 * result_delta_q.toRotationMatrix() * _dt;
//        for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(6+s_i,6+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_V, ILNOISE_GI1) =  kappa_3;

        // V row 4
        V.block<3, 3>(ILO_EPS, ILNOISE_GI) = kappa_6 - 0.5 * _dt * delta_q.toRotationMatrix() * w_sum_skew_fb;
//            tmp33 = - 0.25 * _dt * _dt * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vip1[j])
//                    + 0.5 * _dt * delta_q.toRotationMatrix()* Utility::skewSymmetric(p_br + R_br*fi[j]);
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,3+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_EPS, ILNOISE_GI1) = kappa_6 - 0.5 * _dt * result_delta_q.toRotationMatrix() * w_sum_skew_fbp1;
//            tmp33 = - 0.25 * _dt * _dt * result_delta_q.toRotationMatrix() * Utility::skewSymmetric(vip1[j])
//                    + 0.5 * _dt * result_delta_q.toRotationMatrix()* Utility::skewSymmetric(p_br + R_br*fip1[j]);
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,9+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_EPS, ILNOISE_PHII) = 0.5 * _dt * delta_q.toRotationMatrix() * w_sum_hi;
//            tmp33 =  - 0.5 * _dt * hi[j];
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,18+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_EPS, ILNOISE_PHII1) = 0.5 * _dt * result_delta_q.toRotationMatrix() * w_sum_hip1;
//            tmp33 =  - 0.5 * _dt * hip1[j];
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,21+s_j,tmp33(s_i,s_j)));

        V.block<3, 3>(ILO_EPS, ILNOISE_DPHII) = 0.5 * _dt * delta_q.toRotationMatrix() * w_sum_Ji;
//            tmp33 =  - 0.5 * _dt * delta_q.toRotationMatrix() * R_br * Ji[j];
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,24+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_EPS, ILNOISE_DPHII1) = 0.5 * _dt * result_delta_q.toRotationMatrix() * w_sum_Jip1;
//            tmp33 =  - 0.5 * _dt * result_delta_q.toRotationMatrix() * R_br * Jip1[j];
//            for (int s_i=0; s_i<3;s_i++) for (int s_j=0; s_j<3;s_j++) trp.push_back(Trip(9+3*j+s_i,27+s_j,tmp33(s_i,s_j)));
        V.block<3, 3>(ILO_EPS, ILNOISE_VI) = - 0.5 * _dt * delta_q.toRotationMatrix();
//            for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(9+3*j+s_i,30+3*j+s_i,-_dt));
        V.block<3, 3>(ILO_EPS, ILNOISE_VI1) = - 0.5 * _dt * result_delta_q.toRotationMatrix();


        V.block<3, 3>(ILO_BA, ILNOISE_BA) = -Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(21+s_i,12+s_i,-1*_dt));
        V.block<3, 3>(ILO_BG, ILNOISE_BG) = -Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(24+s_i,15+s_i,-1*_dt));

        V.block<3, 3>(ILO_BV, ILNOISE_BV) = -Matrix3d::Identity() * _dt;
//        for (int s_i=0; s_i<3;s_i++)  trp.push_back(Trip(27+s_i,42+s_i,-1*_dt));

//        Eigen::SparseMatrix<double> V2(RESIDUAL_STATE_SIZE, NOISE_SIZE);
//        V2.setFromTriplets(trp.begin(),trp.end());
//        Eigen::Matrix<double, RESIDUAL_STATE_SIZE, NOISE_SIZE> ss = V2-V;
//        double sum1 = ss.block<3,3>(6,3).sum();
//        std::cout <<"----1-----" << sum1 << std::endl;
//        double sum2 = ss.block<3,3>(6,6).sum();
//        std::cout <<"-----2----" << sum2 << std::endl;
//        double sum3 = ss.block<3,3>(6,9).sum();
//        std::cout <<"-----3----" << sum3 << std::endl;

        jacobian = F * jacobian;

        noise_diag.diagonal()(2) = (ACC_N * ACC_N);
        noise_diag.diagonal()(8) = (ACC_N * ACC_N);
        // scale IMU noise using the diff of the foot contact
        for (int j = 0; j < NUM_OF_LEG; j++) {
            // get z directional contact force ( contact foot sensor reading)
//            double force_mag = 0.5 * (_c_0(3*j+2) + _c_1(3*j+2));
            double diff_c = (_c_1(3*j+2) - _c_0(3*j+2)) / _dt;
            double diff_c_mag = diff_c;
            noise_diag.diagonal()(0) *= (1.0f + 0.0002*abs(diff_c_mag));
            noise_diag.diagonal()(1) *= (1.0f + 0.0002*abs(diff_c_mag));
            noise_diag.diagonal()(2) *= (1.0f + 0.001*abs(diff_c_mag));
            noise_diag.diagonal()(6) *= (1.0f + 0.0002*abs(diff_c_mag));
            noise_diag.diagonal()(7) *= (1.0f + 0.0002*abs(diff_c_mag));
            noise_diag.diagonal()(8) *= (1.0f + 0.001*abs(diff_c_mag));
//            noise.block<3, 3>(30+3*j, 30+3*j) = (uncertainty * uncertainty) * coeff;
            // calculate a sum of delta_epsilon
        }
//        std::cout << "The noise is  " << noise.diagonal().transpose() << std::endl;
//        auto tmp = V * noise * V.transpose();
//        covariance = F * covariance * F.transpose() + tmp;
        covariance = F * covariance * F.transpose() + V * noise_diag * V.transpose();
//        SelfAdjointEigenSolver<Matrix<double, RESIDUAL_STATE_SIZE, RESIDUAL_STATE_SIZE>> eigensolver(tmp);
//        std::cout << "The determinant of V * noise * V.transpose() is " << tmp.determinant() << std::endl;
//        std::cout << eigensolver.eigenvalues().transpose() << std::endl;
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
                                          const Vector3d &delta_v, const Vector3d &sum_delta_epsilon,
                                          const Vector3d &linearized_ba, const Vector3d &linearized_bg,
                                          const Vector3d &linearized_bv) {
    Vector3d result_delta_p;
    Quaterniond result_delta_q;
    Vector3d result_delta_v;
    Vector3d result_sum_delta_epsilon;
    Vector3d result_linearized_ba;
    Vector3d result_linearized_bg;
    Vector3d result_linearized_bv;
    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg, linearized_bv,
                        result_delta_p, result_delta_q, result_delta_v, result_sum_delta_epsilon,
                        result_linearized_ba, result_linearized_bg, result_linearized_bv, 0);

    Vector3d turb_delta_p;
    Quaterniond turb_delta_q;
    Vector3d turb_delta_v;
    Vector3d turb_sum_delta_epsilon;
    Vector3d turb_linearized_ba;
    Vector3d turb_linearized_bg;
    Vector3d turb_linearized_bv;

    Vector3d turb(0.001, -0.03, 0.03);

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p + turb, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg, linearized_bv,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "------------------- check jacobian --------------------" << endl;
    cout << "------------------- check jacobian --------------------" << endl;
    cout << "------------------- check jacobian --------------------" << endl;
    cout << "turb p-----F col 1-------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_P) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_P) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_P) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_P) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_P) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_P) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_P) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q * Quaterniond(1, turb(0) / 2, turb(1) / 2, turb(2) / 2), delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg, linearized_bv,
                        turb_delta_p, turb_delta_q, turb_delta_v, turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "turb q-------F col 2--------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_R) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_R) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_R) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_R) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_R) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_R) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_R) * turb).transpose() << endl;


    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon + turb,
                        linearized_ba, linearized_bg, linearized_bv,
                        turb_delta_p, turb_delta_q, turb_delta_v,  turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "turb epsilon-----------------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_EPS) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_EPS) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_EPS) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_EPS) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_EPS) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_EPS) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_EPS) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba + turb, linearized_bg, linearized_bv,
                        turb_delta_p, turb_delta_q, turb_delta_v,  turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "turb ba-----------------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BA) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BA) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BA) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BA) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BA) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BA) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BA) * turb).transpose() << endl;

    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg + turb, linearized_bv,
                        turb_delta_p, turb_delta_q, turb_delta_v,  turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "turb bg-----------------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BG) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BG) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BG) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BG) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BG) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BG) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BG) * turb).transpose() << endl;


    midPointIntegration(_dt, _acc_0, _gyr_0, _acc_1, _gyr_1,
                        _phi_0, _dphi_0, _c_0, _phi_1, _dphi_1, _c_1,
                        delta_p, delta_q, delta_v, sum_delta_epsilon,
                        linearized_ba, linearized_bg, linearized_bv + turb,
                        turb_delta_p, turb_delta_q, turb_delta_v,  turb_sum_delta_epsilon,
                        turb_linearized_ba, turb_linearized_bg, turb_linearized_bv, 0);
    cout << "turb bv-----------------------------------       " << endl;
    cout << "p diff       " << (turb_delta_p - result_delta_p).transpose() << endl;
    cout << "p jacob diff " << (step_jacobian.block<3, 3>(ILO_P, ILO_BV) * turb).transpose() << endl;
    cout << "q diff       " << ((result_delta_q.inverse() * turb_delta_q).vec() * 2).transpose() << endl;
    cout << "q jacob diff " << (step_jacobian.block<3, 3>(ILO_R, ILO_BV) * turb).transpose() << endl;
    cout << "v diff       " << (turb_delta_v - result_delta_v).transpose() << endl;
    cout << "v jacob diff " << (step_jacobian.block<3, 3>(ILO_V, ILO_BV) * turb).transpose() << endl;
    cout << "eps diff       " << (turb_sum_delta_epsilon - result_sum_delta_epsilon).transpose() << endl;
    cout << "eps jacob diff " << (step_jacobian.block<3, 3>(ILO_EPS, ILO_BV) * turb).transpose() << endl;
    cout << "ba diff      " << (turb_linearized_ba - result_linearized_ba).transpose() << endl;
    cout << "ba jacob diff" << (step_jacobian.block<3, 3>(ILO_BA, ILO_BV) * turb).transpose() << endl;
    cout << "bg diff " << (turb_linearized_bg - result_linearized_bg).transpose() << endl;
    cout << "bg jacob diff " << (step_jacobian.block<3, 3>(ILO_BG, ILO_BV) * turb).transpose() << endl;
    cout << "bv diff " << (turb_linearized_bv - result_linearized_bv).transpose() << endl;
    cout << "bv jacob diff " << (step_jacobian.block<3, 3>(ILO_BV, ILO_BV) * turb).transpose() << endl;

}

Eigen::Matrix<double, RESIDUAL_STATE_SIZE, 1>
IMULegIntegrationBase::evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi,
                                const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
                                const Eigen::Vector3d &Bvi,
                                const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj,
                                const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj,
                                const Eigen::Vector3d &Bvj) {
    Eigen::Matrix<double, RESIDUAL_STATE_SIZE, 1> residuals;


    Eigen::Matrix3d dp_dba = jacobian.block<3, 3>(ILO_P, ILO_BA);
    Eigen::Matrix3d dp_dbg = jacobian.block<3, 3>(ILO_P, ILO_BG);

    Eigen::Matrix3d dq_dbg = jacobian.block<3, 3>(ILO_R, ILO_BG);

    Eigen::Matrix3d dv_dba = jacobian.block<3, 3>(ILO_V, ILO_BA);
    Eigen::Matrix3d dv_dbg = jacobian.block<3, 3>(ILO_V, ILO_BG);

    Eigen::Matrix3d dep1_dbg = jacobian.block<3, 3>(ILO_EPS, ILO_BG);
    Eigen::Matrix3d dep1_dbv = jacobian.block<3, 3>(ILO_EPS, ILO_BV);

    Eigen::Vector3d dba = Bai - linearized_ba; // Bai is the new optimization result
    Eigen::Vector3d dbg = Bgi - linearized_bg;
    Eigen::Vector3d dbv = Bvi - linearized_bv;


    Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * dbg);
    Eigen::Vector3d corrected_delta_v = delta_v + dv_dba * dba + dv_dbg * dbg;
    Eigen::Vector3d corrected_delta_p = delta_p + dp_dba * dba + dp_dbg * dbg;

    Eigen::Vector3d corrected_sum_delta_epsilon;
    corrected_sum_delta_epsilon = sum_delta_epsilon + dep1_dbg * dbg + dep1_dbv * dbv;

    // Test: compare them with repropogarated result
//    repropagate()
//    compare delta_q, delta_v,delta_p, delta_epsilon[j] now with
//    repropagate(Bai, Bgi, linearized_rho + rhoi_rand);
//    std::cout << corrected_delta_q.coeffs() << std::endl;
//    std::cout << delta_q.coeffs() << std::endl;
//    std::cout << corrected_delta_epsilon[0] << std::endl;
//    std::cout << delta_epsilon[0] << std::endl;
//    std::cout << corrected_delta_epsilon[1] << std::endl;
//    std::cout << delta_epsilon[1] << std::endl;
//    std::cout << corrected_delta_epsilon[2] << std::endl;
//    std::cout << delta_epsilon[2] << std::endl;

    residuals.block<3, 1>(ILO_P, 0) = Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt) - corrected_delta_p;
    residuals.block<3, 1>(ILO_R, 0) = 2 * (corrected_delta_q.inverse() * (Qi.inverse() * Qj)).vec();
    residuals.block<3, 1>(ILO_V, 0) = Qi.inverse() * (G * sum_dt + Vj - Vi) - corrected_delta_v;
    residuals.block<3, 1>(ILO_EPS, 0) = Qi.inverse() * (Pj - Pi) - corrected_sum_delta_epsilon;
    residuals.block<3, 1>(ILO_BA, 0) = Baj - Bai;
    residuals.block<3, 1>(ILO_BG, 0) = Bgj - Bgi;
    residuals.block<3, 1>(ILO_BV, 0) = Bvj - Bvi;

    return residuals;
}
