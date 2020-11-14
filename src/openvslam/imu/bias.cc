#include "openvslam/imu/bias.h"

namespace openvslam {
namespace imu {

bias::bias(const float bax, const float bay, const float baz,
           const float bwx, const float bwy, const float bwz)
    : acc_((Vec3_t() << bax, bay, baz).finished()), gyr_((Vec3_t() << bwx, bwy, bwz).finished()) {}

} // namespace imu
} // namespace openvslam
