//
// Created by shuoy on 8/23/21.
//

#ifndef VILEOM_IMU_LEG_FACTOR_H
#define VILEOM_IMU_LEG_FACTOR_H
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "../utils/utility.h"
#include "../utils/parameters.h"
#include "imu_leg_integration_base.h"

class IMULegFactor : public ceres::SizedCostFunction<39, 7, 21, 7, 21>{
public:
    IMULegFactor() = delete;
    IMULegFactor(IMULegIntegrationBase* _il_pre_integration):il_pre_integration(_il_pre_integration)
            {
            }

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    IMULegIntegrationBase* il_pre_integration;
};


#endif //VILEOM_IMU_LEG_FACTOR_H
