#include "openvslam/imu/preintegrator.h"
#include "openvslam/imu/preintegrated.h"

namespace openvslam {
namespace imu {

measurement::measurement(const double acc_x, const double acc_y, const double acc_z,
                         const double gyr_x, const double gyr_y, const double gyr_z,
                         const double dt)
    : acc_(acc_x, acc_y, acc_z), gyr_(gyr_x, gyr_y, gyr_z), dt_(dt) {}

measurement::measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt)
    : acc_(acc), gyr_(gyr), dt_(dt) {}

preintegrator::preintegrator(const bias& b, const std::shared_ptr<config>& cfg) {
    initial_covariance_ << cfg->get_gyr_covariance(), Mat33_t::Zero(), Mat33_t::Zero(), cfg->get_acc_covariance();
    bias_covariance_ << cfg->get_gyr_bias_covariance(), Mat33_t::Zero(), Mat33_t::Zero(), cfg->get_acc_bias_covariance();
    preintegrated_ = eigen_alloc_shared<preintegrated>(b);
    preintegrated_->initialize();
}

void preintegrator::reintegrate() {
    const auto tmp = measurements_;
    measurements_.clear();
    preintegrated_->initialize();
    for (const auto& m : tmp) {
        integrate_new_measurement(m);
    }
}

void preintegrator::merge_previous(const preintegrator& prev) {
    const auto tmp = measurements_;
    measurements_.clear();
    preintegrated_->initialize();
    for (const auto& m : prev.measurements_) {
        integrate_new_measurement(m);
    }
    for (const auto& m : tmp) {
        integrate_new_measurement(m);
    }
}

void preintegrator::integrate_new_measurement(const measurement& m) {
    measurements_.push_back(m);
    preintegrated_->integrate(m.acc_, m.gyr_, m.dt_, initial_covariance_, bias_covariance_);
}

void preintegrator::integrate_new_measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt) {
    measurements_.emplace_back(acc, gyr, dt);
    preintegrated_->integrate(acc, gyr, dt, initial_covariance_, bias_covariance_);
}

} // namespace imu
} // namespace openvslam
