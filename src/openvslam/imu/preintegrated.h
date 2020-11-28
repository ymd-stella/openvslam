#ifndef OPENVSLAM_IMU_PREINTEGRATED_H
#define OPENVSLAM_IMU_PREINTEGRATED_H

#include <vector>
#include <memory>

#include "openvslam/type.h"
#include "openvslam/imu/bias.h"
#include <nlohmann/json_fwd.hpp>

namespace openvslam {
namespace imu {

class preintegrated {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    preintegrated(const bias& b);
    preintegrated(
        double dt,
        const MatRC_t<15, 15>& covariance,
        const bias& b,
        const Mat33_t& delta_rotation,
        const Vec3_t& delta_velocity,
        const Vec3_t& delta_position,
        const Mat33_t& jacob_rotation_gyr,
        const Mat33_t& jacob_velocity_gyr,
        const Mat33_t& jacob_velocity_acc,
        const Mat33_t& jacob_position_gyr,
        const Mat33_t& jacob_position_acc);
    explicit preintegrated(const nlohmann::json& json_preintegrated);
    void initialize();

    void integrate(const Vec3_t& acc, const Vec3_t& gyr, const double dt, const Mat66_t& initial_covariance, const Mat66_t& bias_covariance);

    Mat33_t get_delta_rotation_on_bias(const imu::bias& b);
    Vec3_t get_delta_velocity_on_bias(const imu::bias& b);
    Vec3_t get_delta_position_on_bias(const imu::bias& b);

    MatRC_t<15, 15> get_information();

    nlohmann::json to_json() const;

    double dt_;
    MatRC_t<15, 15> covariance_;
    bias b_;
    Mat33_t delta_rotation_;
    Vec3_t delta_velocity_;
    Vec3_t delta_position_;
    Mat33_t jacob_rotation_gyr_;
    Mat33_t jacob_velocity_gyr_;
    Mat33_t jacob_velocity_acc_;
    Mat33_t jacob_position_gyr_;
    Mat33_t jacob_position_acc_;
};

} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_PREINTEGRATED_H
