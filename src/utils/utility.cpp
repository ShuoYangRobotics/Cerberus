/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#include "utility.h"

Eigen::Matrix3d Utility::g2R(const Eigen::Vector3d &g)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d ng1 = g.normalized();
    Eigen::Vector3d ng2{0, 0, 1.0};
    R0 = Eigen::Quaterniond::FromTwoVectors(ng1, ng2).toRotationMatrix();
    double yaw = Utility::R2ypr(R0).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    // R0 = Utility::ypr2R(Eigen::Vector3d{-90, 0, 0}) * R0;
    return R0;
}

Eigen::Vector3d Utility::lerpGyro(double t, std::vector<std::pair<double, Eigen::Vector3d>> gyroVector)
{
    int idx1, idx2;
    double t1, t2;
    if (t < gyroVector.front().first)
    {
        return gyroVector.front().second;
    }
    else if (t >= gyroVector.back().first)
    {
        return gyroVector.back().second;
    }
    else
    {
        // do the interpolation
        idx1 = 0; idx2 = 1;
        while (t >= gyroVector[idx2].first)
        {
            idx1++; idx2++;
        }
        t1 = gyroVector[idx1].first;
        t2 = gyroVector[idx2].first;
        Eigen::Vector3d vec1 = gyroVector[idx1].second;
        Eigen::Vector3d vec2 = gyroVector[idx2].second;
        return vec1 + (t-t1)*(vec2-vec1)/(t2-t);
    }
}