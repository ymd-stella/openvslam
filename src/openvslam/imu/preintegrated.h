#ifndef OPENVSLAM_IMU_PREINTEGRATED_H
#define OPENVSLAM_IMU_PREINTEGRATED_H

#include <vector>
#include <memory>

#include "openvslam/type.h"
#include "openvslam/imu/bias.h"
#include "openvslam/imu/config.h"

namespace openvslam {
namespace imu {

class measurement {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor for scalar inputs
    measurement(const double acc_x, const double acc_y, const double acc_z,
                const double gyr_x, const double gyr_y, const double gyr_z,
                const double dt);

    //! Constructor for vector inputs
    measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt);

    //! acceleration [m/s^2]
    const Vec3_t acc_;
    //! gyroscope [rad/s]
    const Vec3_t gyr_;
    //! dt [s]
    const double dt_;
};

class preintegrated {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    preintegrated(const bias& b, const std::shared_ptr<config>& conf);
    preintegrated(
        double dt,
        const MatRC_t<15, 15>& covariance,
        const Mat66_t& initial_covariance,
        const Mat66_t& bias_covariance,
        const bias& b,
        const Mat33_t& delta_rotation,
        const Vec3_t& delta_velocity,
        const Vec3_t& delta_position,
        const Mat33_t& jacob_rotation_gyr,
        const Mat33_t& jacob_velocity_gyr,
        const Mat33_t& jacob_velocity_acc,
        const Mat33_t& jacob_position_gyr,
        const Mat33_t& jacob_position_acc,
        const std::vector<measurement>& measurements);
    void initialize();
    void reintegrate();
    void merge_previous(const preintegrated& prev);
    void integrate_new_measurement(const measurement& m);
    void integrate_new_measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt);

    Mat33_t get_delta_rotation_on_bias(const imu::bias& b);
    Vec3_t get_delta_velocity_on_bias(const imu::bias& b);
    Vec3_t get_delta_position_on_bias(const imu::bias& b);

    double dt_;
    MatRC_t<15, 15> covariance_;
    MatRC_t<15, 15> information_;
    Mat66_t initial_covariance_;
    Mat66_t bias_covariance_;
    bias b_;
    Mat33_t delta_rotation_;
    Vec3_t delta_velocity_;
    Vec3_t delta_position_;
    Mat33_t jacob_rotation_gyr_;
    Mat33_t jacob_velocity_gyr_;
    Mat33_t jacob_velocity_acc_;
    Mat33_t jacob_position_gyr_;
    Mat33_t jacob_position_acc_;
    std::vector<measurement> measurements_;

private:
    void integrate(const Vec3_t& acc, const Vec3_t& gyr, const double dt);
};

} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_PREINTEGRATED_H
