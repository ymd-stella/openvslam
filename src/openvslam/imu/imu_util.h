#ifndef OPENVSLAM_IMU_UTIL_H
#define OPENVSLAM_IMU_UTIL_H

#include <vector>

#include "openvslam/type.h"

namespace openvslam {
namespace data {
class keyframe;
} // namespace data

namespace imu {
class data;

class imu_util {
public:
    static void preprocess_imu(const imu::data& imu1, const imu::data& imu2,
                               bool interpolate_last, bool interpolate_curr,
                               double last_stamp, double curr_stamp,
                               Vec3_t& acc, Vec3_t& gyr, double& dt);
    static std::vector<openvslam::data::keyframe*> gather_intertial_ref_keyframes(openvslam::data::keyframe* keyfrm);
    static void compute_velocity(std::vector<openvslam::data::keyframe*>& keyfrms);
    static Mat33_t compute_gravity_dir(std::vector<openvslam::data::keyframe*>& keyfrms);
};

} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_UTIL_H
