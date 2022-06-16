//
// Created by shuoy on 8/23/21.
//

#ifndef VILO_IMU_LEG_FACTOR_H
#define VILO_IMU_LEG_FACTOR_H
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "../utils/utility.h"
#include "../utils/parameters.h"
#include "imu_leg_integration_base.h"

class IMULegFactor : public ceres::SizedCostFunction<ERROR_STATE_SIZE, 7, 9, TOTAL_BIAS_SIZE, 7, 9, TOTAL_BIAS_SIZE> {
public:
    IMULegFactor() = delete;

    IMULegFactor(IMULegIntegrationBase *_il_pre_integration) : il_pre_integration(_il_pre_integration) {
    }

    // para_Pose[i], para_SpeedBias[i], para_LegBias[i], para_Pose[j], para_SpeedBias[j], para_LegBias[j]
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void checkJacobian(double const *const *parameters);
    IMULegIntegrationBase *il_pre_integration;
};
#endif //VILO_IMU_LEG_FACTOR_H
