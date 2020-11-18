#include "openvslam/imu/preintegrated.h"
#include "openvslam/util/converter.h"

namespace openvslam {
namespace imu {

measurement::measurement(const double acc_x, const double acc_y, const double acc_z,
                         const double gyr_x, const double gyr_y, const double gyr_z,
                         const double dt)
    : acc_(acc_x, acc_y, acc_z), gyr_(gyr_x, gyr_y, gyr_z), dt_(dt) {}

measurement::measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt)
    : acc_(acc), gyr_(gyr), dt_(dt) {}

preintegrated::preintegrated(const bias& b, const std::shared_ptr<config>& cfg)
    : b_(b) {
    Mat33_t acc_cov = cfg->get_acc_covariance();
    Mat33_t gyr_cov = cfg->get_gyr_covariance();
    initial_covariance_ << gyr_cov, Mat33_t::Zero(), Mat33_t::Zero(), acc_cov;
    bias_covariance_ << cfg->get_gyr_bias_covariance(), Mat33_t::Zero(), Mat33_t::Zero(), cfg->get_acc_bias_covariance();
    initialize();
}

preintegrated::preintegrated(
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
    const std::vector<measurement>& measurements)
    : dt_(dt), covariance_(covariance), initial_covariance_(initial_covariance), bias_covariance_(bias_covariance),
      b_(b), delta_rotation_(delta_rotation), delta_velocity_(delta_velocity), delta_position_(delta_position),
      jacob_rotation_gyr_(jacob_rotation_gyr), jacob_velocity_gyr_(jacob_velocity_gyr), jacob_velocity_acc_(jacob_velocity_acc),
      jacob_position_gyr_(jacob_position_gyr), jacob_position_acc_(jacob_position_acc), measurements_(measurements) {}

void preintegrated::initialize() {
    delta_rotation_.setIdentity();
    delta_velocity_.setZero();
    delta_position_.setZero();
    jacob_rotation_gyr_.setZero();
    jacob_velocity_gyr_.setZero();
    jacob_velocity_acc_.setZero();
    jacob_position_gyr_.setZero();
    jacob_position_acc_.setZero();
    covariance_.setZero();
    dt_ = 0.0f;
    measurements_.clear();
}

void preintegrated::reintegrate() {
    const auto tmp = measurements_;
    initialize();
    for (const auto& m : tmp) {
        integrate_new_measurement(m);
    }
}

void preintegrated::merge_previous(const preintegrated& prev) {
    const auto tmp = measurements_;
    initialize();
    for (const auto& m : prev.measurements_) {
        integrate_new_measurement(m);
    }
    for (const auto& m : tmp) {
        integrate_new_measurement(m);
    }
}

void preintegrated::integrate_new_measurement(const measurement& m) {
    measurements_.push_back(m);
    integrate(m.acc_, m.gyr_, m.dt_);
}

void preintegrated::integrate_new_measurement(const Vec3_t& acc, const Vec3_t& gyr, const double dt) {
    measurements_.emplace_back(acc, gyr, dt);
    integrate(acc, gyr, dt);
}

void preintegrated::integrate(const Vec3_t& acc, const Vec3_t& gyr, const double dt) {
    // The reference is "On-Manifold Preintegration for Real-Time Visual-Inertial Odometry"
    Vec3_t unbiased_acc = acc - b_.acc_;
    Vec3_t unbiased_gyr = gyr - b_.gyr_;

    // (33)
    delta_position_ = delta_position_ + delta_velocity_ * dt + 0.5 * delta_rotation_ * unbiased_acc * dt * dt;
    delta_velocity_ = delta_velocity_ + delta_rotation_ * unbiased_acc * dt;

    // Iterative Noise Propagation (Appendix A)
    Mat33_t w_acc;
    w_acc << 0, -unbiased_acc(2), unbiased_acc(1),
        unbiased_acc(2), 0, -unbiased_acc(0),
        -unbiased_acc(1), unbiased_acc(0), 0;
    MatRC_t<9, 9> A = MatRC_t<9, 9>::Identity();
    MatRC_t<9, 6> B = MatRC_t<9, 6>::Zero();
    A.block<3, 3>(3, 0) = -delta_rotation_ * dt * w_acc;
    A.block<3, 3>(6, 0) = -0.5 * delta_rotation_ * dt * dt * w_acc;
    A.block<3, 3>(6, 3) = Mat33_t::Identity() * dt;
    B.block<3, 3>(3, 3) = delta_rotation_ * dt;
    B.block<3, 3>(6, 3) = 0.5 * delta_rotation_ * dt * dt;

    // Bias Correction via First-Order Updates (Appendix B)
    jacob_position_acc_ = jacob_position_acc_ + jacob_velocity_acc_ * dt - 0.5 * delta_rotation_ * dt * dt;
    jacob_position_gyr_ = jacob_position_gyr_ + jacob_velocity_gyr_ * dt
                          - 0.5 * delta_rotation_ * dt * dt * w_acc * jacob_rotation_gyr_;
    jacob_velocity_acc_ = jacob_velocity_acc_ - delta_rotation_ * dt;
    jacob_velocity_gyr_ = jacob_velocity_gyr_ - delta_rotation_ * dt * w_acc * jacob_rotation_gyr_;

    // (33)
    const Vec3_t v = (unbiased_gyr - b_.gyr_) * dt;
    Mat33_t delta_rotation = util::converter::exp_so3(v);
    Mat33_t right_jacobian = util::converter::right_jacobian_so3(v);

    delta_rotation_ = util::converter::normalize_rotation(delta_rotation_ * delta_rotation);

    // Compute rotation parts of matrices A and B (Appendix A)
    A.block<3, 3>(0, 0) = delta_rotation.transpose();
    B.block<3, 3>(0, 0) = right_jacobian * dt;

    // Update covariance matrix (63)
    covariance_.block<9, 9>(0, 0) = A * covariance_.block<9, 9>(0, 0) * A.transpose()
                                    + B * initial_covariance_ * B.transpose();
    // Update noise density parts of covariance matrix
    covariance_.block<6, 6>(9, 9) = covariance_.block<6, 6>(9, 9) + bias_covariance_;

    // Bias Correction via First-Order Updates (Appendix B)
    jacob_rotation_gyr_ = A.block<3, 3>(0, 0) * jacob_rotation_gyr_ - B.block<3, 3>(0, 0);

    // Update total integrated time
    dt_ += dt;
}

Mat33_t preintegrated::get_delta_rotation_on_bias(const imu::bias& b) {
    Vec3_t dbg = b.gyr_ - b_.gyr_;
    return util::converter::normalize_rotation(delta_rotation_ * util::converter::exp_so3(jacob_rotation_gyr_ * dbg));
}

Vec3_t preintegrated::get_delta_velocity_on_bias(const imu::bias& b) {
    Vec3_t dbg = b.gyr_ - b_.gyr_;
    Vec3_t dba = b.acc_ - b_.acc_;
    return delta_velocity_ + jacob_velocity_gyr_ * dbg + jacob_velocity_acc_ * dba;
}

Vec3_t preintegrated::get_delta_position_on_bias(const imu::bias& b) {
    Vec3_t dbg = b.gyr_ - b_.gyr_;
    Vec3_t dba = b.acc_ - b_.acc_;
    return delta_position_ + jacob_position_gyr_ * dbg + jacob_position_acc_ * dba;
}

MatRC_t<15, 15> preintegrated::get_information() {
    MatRC_t<15, 15> information;
    information.block<9, 9>(0, 0) = covariance_.block<9, 9>(0, 0).ldlt().solve(MatRC_t<9, 9>::Identity());
    information.block<6, 6>(9, 9) = (covariance_.block<6, 6>(9, 9).diagonal().cwiseInverse()).asDiagonal();
    return information;
}

} // namespace imu
} // namespace openvslam
