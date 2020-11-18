#ifndef OPENVSLAM_IMU_INTERNAL_GRAVITY_SCALE_EDGE_H
#define OPENVSLAM_IMU_INTERNAL_GRAVITY_SCALE_EDGE_H

#include "openvslam/type.h"
#include "openvslam/util/converter.h"
#include "openvslam/imu/bias.h"
#include "openvslam/imu/constant.h"
#include "openvslam/imu/config.h"
#include "openvslam/imu/preintegrated.h"
#include "openvslam/optimize/internal/se3/shot_vertex.h"
#include "openvslam/imu/internal/velocity_vertex.h"
#include "openvslam/imu/internal/bias_vertex.h"
#include "openvslam/imu/internal/gravity_dir_vertex.h"
#include "openvslam/imu/internal/scale_vertex.h"

#include <g2o/core/base_multi_edge.h>

namespace openvslam {
namespace imu {
namespace internal {

class inertial_gravity_scale_edge final : public g2o::BaseMultiEdge<9, VecR_t<9>> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inertial_gravity_scale_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override;

    void linearizeOplus(std::vector<MatX_t>& jacobianOplus);
    void linearizeOplus() override;

    const Vec3_t jacob_rotation_gyr, jacob_velocity_gyr, jacob_position_gyr;
    const Vec3_t jacob_velocity_acc, jacob_position_acc;
    std::shared_ptr<preintegrated> imu_preintegrated_ = nullptr;
    Vec3_t gravity;
};

inline inertial_gravity_scale_edge::inertial_gravity_scale_edge()
    : g2o::BaseMultiEdge<9, VecR_t<9>>() {
    resize(8);
    gravity << 0, 0, -imu::constant::gravity();
}

inline bool inertial_gravity_scale_edge::read(std::istream& is) {
    (void)is;
    return false;
}

inline bool inertial_gravity_scale_edge::write(std::ostream& os) const {
    (void)os;
    return false;
}

inline void inertial_gravity_scale_edge::computeError() {
    const auto keyfrm_vtx1 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[0]);
    const auto velocity_vtx1 = static_cast<const velocity_vertex*>(_vertices[1]);
    const auto gyr_bias_vtx = static_cast<const bias_vertex*>(_vertices[2]);
    const auto acc_bias_vtx = static_cast<const bias_vertex*>(_vertices[3]);
    const auto keyfrm_vtx2 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[4]);
    const auto velocity_vtx2 = static_cast<const velocity_vertex*>(_vertices[5]);
    const auto gravity_dir_vtx = static_cast<const gravity_dir_vertex*>(_vertices[6]);
    const auto scale_vtx = static_cast<const scale_vertex*>(_vertices[7]);

    const bias b(acc_bias_vtx->estimate(), gyr_bias_vtx->estimate());
    const Mat33_t delta_rotation = imu_preintegrated_->get_delta_rotation_on_bias(b);
    const Vec3_t delta_velocity = imu_preintegrated_->get_delta_velocity_on_bias(b);
    const Vec3_t delta_position = imu_preintegrated_->get_delta_position_on_bias(b);
    const double dt = imu_preintegrated_->dt_;

    const Mat33_t Riw1 = keyfrm_vtx1->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi1 = Riw1.transpose();
    const Vec3_t twi1 = -Rwi1 * keyfrm_vtx1->estimate().translation();
    const Mat33_t Riw2 = keyfrm_vtx2->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi2 = Riw2.transpose();
    const Vec3_t twi2 = -Rwi2 * keyfrm_vtx2->estimate().translation();
    const Vec3_t v1 = velocity_vtx1->estimate();
    const Vec3_t v2 = velocity_vtx2->estimate();

    const Vec3_t g = gravity_dir_vtx->estimate() * gravity;
    const double s = scale_vtx->estimate();

    // The reference is "Inertial-Only Optimization for Visual-Inertial Initialization"
    const Vec3_t error_rotation = util::converter::log_so3(delta_rotation.transpose() * Riw1 * Rwi2);        // (7)
    const Vec3_t error_velocity = Riw1 * (s * (v2 - v1) - g * dt) - delta_velocity;                          // (8)
    const Vec3_t error_position = Riw1 * (s * (twi2 - twi1 - v1 * dt) - 0.5 * g * dt * dt) - delta_position; // (9)

    _error << error_rotation, error_velocity, error_position;
}

inline void inertial_gravity_scale_edge::linearizeOplus(std::vector<MatX_t>& jacobianOplus) {
    const auto keyfrm_vtx1 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[0]);
    const auto velocity_vtx1 = static_cast<const velocity_vertex*>(_vertices[1]);
    const auto gyr_bias_vtx = static_cast<const bias_vertex*>(_vertices[2]);
    const auto acc_bias_vtx = static_cast<const bias_vertex*>(_vertices[3]);
    const auto keyfrm_vtx2 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[4]);
    const auto velocity_vtx2 = static_cast<const velocity_vertex*>(_vertices[5]);
    const auto gravity_dir_vtx = static_cast<const gravity_dir_vertex*>(_vertices[6]);
    const auto scale_vtx = static_cast<const scale_vertex*>(_vertices[7]);

    const imu::bias b(acc_bias_vtx->estimate(), gyr_bias_vtx->estimate());
    const imu::bias& b0 = imu_preintegrated_->b_;
    const Vec3_t delta_bias_gyr = b.gyr_ - b0.gyr_;

    const Mat33_t jacob_rotation_gyr = imu_preintegrated_->jacob_rotation_gyr_;
    const Mat33_t jacob_velocity_gyr = imu_preintegrated_->jacob_velocity_gyr_;
    const Mat33_t jacob_position_gyr = imu_preintegrated_->jacob_position_gyr_;
    const Mat33_t jacob_velocity_acc = imu_preintegrated_->jacob_velocity_acc_;
    const Mat33_t jacob_position_acc = imu_preintegrated_->jacob_position_acc_;

    const Mat33_t Riw1 = keyfrm_vtx1->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi1 = Riw1.transpose();
    const Vec3_t twi1 = -Rwi1 * keyfrm_vtx1->estimate().translation();
    const Mat33_t Riw2 = keyfrm_vtx2->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi2 = Riw2.transpose();
    const Vec3_t twi2 = -Rwi2 * keyfrm_vtx2->estimate().translation();

    const Mat33_t Rwg = gravity_dir_vtx->estimate();
    MatRC_t<3, 2> Gm;
    Gm << 0.0, -imu::constant::gravity(),
        imu::constant::gravity(), 0.0,
        0.0, 0.0;
    const MatRC_t<3, 2> dGdTheta = Rwg * Gm;

    const Vec3_t g = Rwg * gravity;
    const double s = scale_vtx->estimate();

    const Vec3_t v1 = velocity_vtx1->estimate();
    const Vec3_t v2 = velocity_vtx2->estimate();

    const Mat33_t delta_rotation = imu_preintegrated_->get_delta_rotation_on_bias(b);
    const Mat33_t error_rotation = delta_rotation.transpose() * Riw1 * Rwi2;
    const Mat33_t inv_right_jacobian = util::converter::inverse_right_jacobian_so3(util::converter::log_so3(error_rotation));
    const double dt = imu_preintegrated_->dt_;

    // The reference is "On-Manifold Preintegration for Real-Time Visual-Inertial Odometry", Appendix C
    // Jacobians wrt Pose 1
    jacobianOplus[0].setZero();
    // rotation
    jacobianOplus[0].block<3, 3>(0, 0) = -inv_right_jacobian * Riw2 * Rwi1;
    jacobianOplus[0].block<3, 3>(3, 0) = util::converter::to_skew_symmetric_mat(Riw1 * (s * (v2 - v1) - g * dt));
    jacobianOplus[0].block<3, 3>(6, 0) = util::converter::to_skew_symmetric_mat(Riw1 * (s * (twi2 - twi1 - v1 * dt) - 0.5 * g * dt * dt));
    // translation
    jacobianOplus[0].block<3, 3>(6, 3) = -s * Eigen::Matrix3d::Identity();

    // Jacobians wrt Velocity 1
    jacobianOplus[1].setZero();
    jacobianOplus[1].block<3, 3>(3, 0) = -s * Riw1;
    jacobianOplus[1].block<3, 3>(6, 0) = -s * Riw1 * dt;

    // Jacobians wrt Gyro 1
    jacobianOplus[2].setZero();
    jacobianOplus[2].block<3, 3>(0, 0) = -inv_right_jacobian * error_rotation.transpose()
                                         * util::converter::right_jacobian_so3(jacob_rotation_gyr * delta_bias_gyr) * jacob_rotation_gyr;
    jacobianOplus[2].block<3, 3>(3, 0) = -jacob_velocity_gyr;
    jacobianOplus[2].block<3, 3>(6, 0) = -jacob_position_gyr;

    // Jacobians wrt Accelerometer 1
    jacobianOplus[3].setZero();
    jacobianOplus[3].block<3, 3>(3, 0) = -jacob_velocity_acc;
    jacobianOplus[3].block<3, 3>(6, 0) = -jacob_position_acc;

    // Jacobians wrt Pose 2
    jacobianOplus[4].setZero();
    // rotation
    jacobianOplus[4].block<3, 3>(0, 0) = inv_right_jacobian;
    // translation
    jacobianOplus[4].block<3, 3>(6, 3) = s * Riw1 * Rwi2;

    // Jacobians wrt Velocity 2
    jacobianOplus[5].setZero();
    jacobianOplus[5].block<3, 3>(3, 0) = s * Riw1;

    // The reference is "Inertial-Only Optimization for Visual-Inertial Initialization"
    // Jacobians wrt Gravity direction
    jacobianOplus[6].setZero();
    jacobianOplus[6].block<3, 2>(3, 0) = -Riw1 * dGdTheta * dt;            // (17)
    jacobianOplus[6].block<3, 2>(6, 0) = -0.5 * Riw1 * dGdTheta * dt * dt; // (18)

    // Jacobians wrt scale factor
    jacobianOplus[7].setZero();
    jacobianOplus[7].block<3, 1>(3, 0) = Riw1 * (v2 - v1);               // (14)
    jacobianOplus[7].block<3, 1>(6, 0) = Riw1 * (twi2 - twi1 - v1 * dt); // (15)
}

inline void inertial_gravity_scale_edge::linearizeOplus() {
    const auto keyfrm_vtx1 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[0]);
    const auto velocity_vtx1 = static_cast<const velocity_vertex*>(_vertices[1]);
    const auto gyr_bias_vtx = static_cast<const bias_vertex*>(_vertices[2]);
    const auto acc_bias_vtx = static_cast<const bias_vertex*>(_vertices[3]);
    const auto keyfrm_vtx2 = static_cast<const optimize::internal::se3::shot_vertex*>(_vertices[4]);
    const auto velocity_vtx2 = static_cast<const velocity_vertex*>(_vertices[5]);
    const auto gravity_dir_vtx = static_cast<const gravity_dir_vertex*>(_vertices[6]);
    const auto scale_vtx = static_cast<const scale_vertex*>(_vertices[7]);

    const imu::bias b(acc_bias_vtx->estimate(), gyr_bias_vtx->estimate());
    const imu::bias& b0 = imu_preintegrated_->b_;
    const Vec3_t delta_bias_gyr = b.gyr_ - b0.gyr_;

    const Mat33_t jacob_rotation_gyr = imu_preintegrated_->jacob_rotation_gyr_;
    const Mat33_t jacob_velocity_gyr = imu_preintegrated_->jacob_velocity_gyr_;
    const Mat33_t jacob_position_gyr = imu_preintegrated_->jacob_position_gyr_;
    const Mat33_t jacob_velocity_acc = imu_preintegrated_->jacob_velocity_acc_;
    const Mat33_t jacob_position_acc = imu_preintegrated_->jacob_position_acc_;

    const Mat33_t Riw1 = keyfrm_vtx1->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi1 = Riw1.transpose();
    const Vec3_t twi1 = -Rwi1 * keyfrm_vtx1->estimate().translation();
    const Mat33_t Riw2 = keyfrm_vtx2->estimate().rotation().toRotationMatrix();
    const Mat33_t Rwi2 = Riw2.transpose();
    const Vec3_t twi2 = -Rwi2 * keyfrm_vtx2->estimate().translation();

    const Mat33_t Rwg = gravity_dir_vtx->estimate();
    MatRC_t<3, 2> Gm;
    Gm << 0.0, -imu::constant::gravity(),
        imu::constant::gravity(), 0.0,
        0.0, 0.0;
    const MatRC_t<3, 2> dGdTheta = Rwg * Gm;

    const Vec3_t g = Rwg * gravity;
    const double s = scale_vtx->estimate();

    const Vec3_t v1 = velocity_vtx1->estimate();
    const Vec3_t v2 = velocity_vtx2->estimate();

    const Mat33_t delta_rotation = imu_preintegrated_->get_delta_rotation_on_bias(b);
    const Mat33_t error_rotation = delta_rotation.transpose() * Riw1 * Rwi2;
    const Mat33_t inv_right_jacobian = util::converter::inverse_right_jacobian_so3(util::converter::log_so3(error_rotation));
    const double dt = imu_preintegrated_->dt_;

    // Jacobians wrt Pose 1
    _jacobianOplus[0].setZero();
    // rotation
    _jacobianOplus[0].block<3, 3>(0, 0) = -inv_right_jacobian * Riw2 * Rwi1;
    _jacobianOplus[0].block<3, 3>(3, 0) = util::converter::to_skew_symmetric_mat(Riw1 * (s * (v2 - v1) - g * dt));
    _jacobianOplus[0].block<3, 3>(6, 0) = util::converter::to_skew_symmetric_mat(Riw1 * (s * (twi2 - twi1 - v1 * dt) - 0.5 * g * dt * dt));
    // translation
    _jacobianOplus[0].block<3, 3>(6, 3) = -s * Eigen::Matrix3d::Identity();

    // Jacobians wrt Velocity 1
    _jacobianOplus[1].setZero();
    _jacobianOplus[1].block<3, 3>(3, 0) = -s * Riw1;
    _jacobianOplus[1].block<3, 3>(6, 0) = -s * Riw1 * dt;

    // Jacobians wrt Gyro 1
    _jacobianOplus[2].setZero();
    _jacobianOplus[2].block<3, 3>(0, 0) = -inv_right_jacobian * error_rotation.transpose()
                                          * util::converter::right_jacobian_so3(jacob_rotation_gyr * delta_bias_gyr) * jacob_rotation_gyr;
    _jacobianOplus[2].block<3, 3>(3, 0) = -jacob_velocity_gyr;
    _jacobianOplus[2].block<3, 3>(6, 0) = -jacob_position_gyr;

    // Jacobians wrt Accelerometer 1
    _jacobianOplus[3].setZero();
    _jacobianOplus[3].block<3, 3>(3, 0) = -jacob_velocity_acc;
    _jacobianOplus[3].block<3, 3>(6, 0) = -jacob_position_acc;

    // Jacobians wrt Pose 2
    _jacobianOplus[4].setZero();
    // rotation
    _jacobianOplus[4].block<3, 3>(0, 0) = inv_right_jacobian;
    // translation
    _jacobianOplus[4].block<3, 3>(6, 3) = s * Riw1 * Rwi2;

    // Jacobians wrt Velocity 2
    _jacobianOplus[5].setZero();
    _jacobianOplus[5].block<3, 3>(3, 0) = s * Riw1;

    // Jacobians wrt Gravity direction
    _jacobianOplus[6].setZero();
    _jacobianOplus[6].block<3, 2>(3, 0) = -Riw1 * dGdTheta * dt;
    _jacobianOplus[6].block<3, 2>(6, 0) = -0.5 * Riw1 * dGdTheta * dt * dt;

    // Jacobians wrt scale factor
    _jacobianOplus[7].setZero();
    _jacobianOplus[7].block<3, 1>(3, 0) = Riw1 * (v2 - v1);
    _jacobianOplus[7].block<3, 1>(6, 0) = Riw1 * (twi2 - twi1 - v1 * dt);
}

} // namespace internal
} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_INTERNAL_GRAVITY_SCALE_EDGE_H
