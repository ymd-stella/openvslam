#ifndef OPENVSLAM_IMU_BIAS_H
#define OPENVSLAM_IMU_BIAS_H

#include "openvslam/type.h"

#include <nlohmann/json_fwd.hpp>

namespace openvslam {
namespace imu {

class bias {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    bias()
        : bias(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}
    bias(const float bax, const float bay, const float baz,
         const float bwx, const float bwy, const float bwz);
    bias(const Vec3_t& acc, const Vec3_t& gyr)
        : acc_(acc), gyr_(gyr) {}
    explicit bias(const nlohmann::json& json_bias);

    nlohmann::json to_json() const;

    Vec3_t acc_;
    Vec3_t gyr_;
};

} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_BIAS_H
