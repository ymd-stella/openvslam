#include "openvslam/imu/imu_util.h"
#include "openvslam/imu/data.h"
#include "openvslam/data/keyframe.h"
#include "openvslam/util/converter.h"
#include "openvslam/imu/preintegrator.h"
#include "openvslam/imu/preintegrated.h"

namespace openvslam {
namespace imu {

void imu_util::preprocess_imu(const imu::data& imu1, const imu::data& imu2,
                              bool interpolate_last, bool interpolate_curr,
                              double last_stamp, double curr_stamp,
                              Vec3_t& acc, Vec3_t& gyr, double& dt) {
    // Take the average of two consecutive values for numerical integration with the trapezoidal rule.
    Vec3_t acc1, gyr1, acc2, gyr2;
    if (interpolate_last && std::abs(last_stamp - imu1.ts_) > 1e-6) {
        double dt1f = last_stamp - imu1.ts_;
        double dt12 = imu2.ts_ - imu1.ts_;
        acc1 = imu2.acc_;
        gyr1 = imu2.gyr_;
        // interpolation
        acc2 = imu1.acc_ + (imu2.acc_ - imu1.acc_) * (dt1f / dt12);
        gyr2 = imu1.gyr_ + (imu2.gyr_ - imu1.gyr_) * (dt1f / dt12);
        dt = dt1f;
    }
    else if (interpolate_curr && std::abs(curr_stamp - imu2.ts_) > 1e-6) {
        double dt1f = imu2.ts_ - curr_stamp;
        double dt12 = imu2.ts_ - imu1.ts_;
        acc1 = imu1.acc_;
        gyr1 = imu1.gyr_;
        // interpolation
        acc2 = imu1.acc_ + (imu2.acc_ - imu1.acc_) * (dt1f / dt12);
        gyr2 = imu1.gyr_ + (imu2.gyr_ - imu1.gyr_) * (dt1f / dt12);
        dt = dt1f;
    }
    else {
        acc1 = imu1.acc_;
        gyr1 = imu1.gyr_;
        acc2 = imu2.acc_;
        gyr2 = imu2.gyr_;
        dt = imu2.ts_ - imu1.ts_;
    }
    acc = (acc1 + acc2) * 0.5;
    gyr = (gyr1 + gyr2) * 0.5;
}

std::vector<openvslam::data::keyframe*> imu_util::gather_intertial_ref_keyframes(openvslam::data::keyframe* keyfrm) {
    std::vector<openvslam::data::keyframe*> keyfrms;
    keyfrms.push_back(keyfrm);
    while (keyfrms.back()->inertial_ref_keyfrm_) {
        keyfrms.push_back(keyfrms.back()->inertial_ref_keyfrm_);
    }
    return keyfrms;
}

void imu_util::compute_velocity(std::vector<openvslam::data::keyframe*>& keyfrms) {
    for (const auto& keyfrm : keyfrms) {
        assert(keyfrm->imu_preintegrator_from_inertial_ref_keyfrm_);
        assert(keyfrm->inertial_ref_keyfrm_);

        const Mat44_t pose1_wi = keyfrm->inertial_ref_keyfrm_->get_imu_pose_inv();
        const Vec3_t twi1 = pose1_wi.block<3, 1>(0, 3);
        const Mat44_t pose2_wi = keyfrm->get_imu_pose_inv();
        const Vec3_t twi2 = pose2_wi.block<3, 1>(0, 3);
        const Vec3_t velocity = (twi2 - twi1) / keyfrm->imu_preintegrator_from_inertial_ref_keyfrm_->preintegrated_->dt_;
        keyfrm->inertial_ref_keyfrm_->velocity_ = velocity;
    }
}

Mat33_t imu_util::compute_gravity_dir(std::vector<openvslam::data::keyframe*>& keyfrms) {
    Vec3_t integrated_gravity = Vec3_t::Zero();
    for (const auto keyfrm : keyfrms) {
        assert(keyfrm->imu_preintegrator_from_inertial_ref_keyfrm_);
        assert(keyfrm->inertial_ref_keyfrm_);

        const Mat44_t pose1_wi = keyfrm->inertial_ref_keyfrm_->get_imu_pose_inv();
        const Mat33_t Rwi1 = pose1_wi.block<3, 3>(0, 0);
        integrated_gravity -= Rwi1 * keyfrm->imu_preintegrator_from_inertial_ref_keyfrm_->preintegrated_->delta_velocity_;
    }

    integrated_gravity = integrated_gravity / integrated_gravity.norm();
    const Vec3_t ez = (Vec3_t() << 0.0, 0.0, -1.0).finished();
    const Vec3_t normal = ez.cross(integrated_gravity);
    const Vec3_t dir_normal = normal / normal.norm();
    const double cos_gravity = ez.dot(integrated_gravity);
    const Vec3_t rot_wg = dir_normal * acos(cos_gravity);
    Mat33_t Rwg = util::converter::exp_so3(rot_wg);
    return Rwg;
}

} // namespace imu
} // namespace openvslam
