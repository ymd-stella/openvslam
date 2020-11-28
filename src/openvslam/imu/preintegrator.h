#ifndef OPENVSLAM_IMU_PREINTEGRATOR_H
#define OPENVSLAM_IMU_PREINTEGRATOR_H

#include <vector>
#include <memory>

#include "openvslam/type.h"
#include "openvslam/imu/bias.h"
#include "openvslam/imu/config.h"

namespace openvslam {
namespace imu {

class preintegrated;

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

class preintegrator {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    preintegrator(const bias& b, const std::shared_ptr<config>& conf);
    void reintegrate();
    void merge_previous(const preintegrator& prev);
    void integrate_new_measurement(const measurement& m);
    void integrate_new_measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt);

    Mat66_t initial_covariance_;
    Mat66_t bias_covariance_;
    std::shared_ptr<preintegrated> preintegrated_;
    eigen_alloc_vector<measurement> measurements_;
};

} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_PREINTEGRATOR_H
