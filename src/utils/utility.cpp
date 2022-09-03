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
        idx1 = 0;
        idx2 = 1;
        while (t >= gyroVector[idx2].first)
        {
            idx1++;
            idx2++;
        }
        t1 = gyroVector[idx1].first;
        t2 = gyroVector[idx2].first;
        Eigen::Vector3d vec1 = gyroVector[idx1].second;
        Eigen::Vector3d vec2 = gyroVector[idx2].second;
        return vec1 + (t - t1) * (vec2 - vec1) / (t2 - t);
    }
}

Eigen::Matrix<double, 12, 3> Utility::lerpLegSensors(double t, int &starting_idx,
                                                     std::deque<std::pair<double, Vector12d>> jointAngVector,
                                                     std::deque<std::pair<double, Vector12d>> jointAngVelVector,
                                                     std::deque<std::pair<double, Vector12d>> footForceVector)
{
    int idx1, idx2;
    double t1, t2;
    Eigen::Matrix<double, 12, 3> out;
    out.setZero();
    if (t < jointAngVector.front().first)
    {
        out.col(0) = jointAngVector.front().second;
        out.col(1) = jointAngVelVector.front().second;
        out.col(2) = footForceVector.front().second;
        starting_idx = 0;
        return out;
    }
    else if (t >= jointAngVector.back().first)
    {
        out.col(0) = jointAngVector.back().second;
        out.col(1) = jointAngVelVector.back().second;
        out.col(2) = footForceVector.back().second;
        starting_idx = jointAngVector.size();
        return out;
    }
    else
    {
        // do the interpolation
        idx1 = starting_idx;
        idx2 = starting_idx + 1;
        while (t >= jointAngVector[idx2].first)
        {
            idx1++;
            idx2++;
        }
        t1 = jointAngVector[idx1].first;
        t2 = jointAngVector[idx2].first;
        Eigen::Matrix<double, 12, 1> vec1 = jointAngVector[idx1].second;
        Eigen::Matrix<double, 12, 1> vec2 = jointAngVector[idx2].second;
        out.col(0) = vec1 + (t - t1) * (vec2 - vec1) / (t2 - t1);
        vec1 = jointAngVelVector[idx1].second;
        vec2 = jointAngVelVector[idx2].second;
        out.col(1) = vec1 + (t - t1) * (vec2 - vec1) / (t2 - t1);
        vec1 = footForceVector[idx1].second;
        vec2 = footForceVector[idx2].second;
        out.col(2) = vec1 + (t - t1) * (vec2 - vec1) / (t2 - t1);

        starting_idx = idx1;
        return out;
    }
}