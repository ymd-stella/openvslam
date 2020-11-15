#ifndef OPENVSLAM_IMU_INTERNAL_PRIOR_BIAS_EDGE_H
#define OPENVSLAM_IMU_INTERNAL_PRIOR_BIAS_EDGE_H

#include <g2o/core/base_unary_edge.h>

namespace openvslam {
namespace imu {
namespace internal {

class prior_bias_edge final : public g2o::BaseUnaryEdge<3, Vec3_t, bias_vertex> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    prior_bias_edge();

    bool read(std::istream& is) override;

    bool write(std::ostream& os) const override;

    void computeError() override;

    void linearizeOplus() override;

    Vec3_t prior_bias_;
};

inline prior_bias_edge::prior_bias_edge()
    : g2o::BaseUnaryEdge<3, Vec3_t, bias_vertex>() {
}

inline bool prior_bias_edge::read(std::istream& is) {
    (void)is;
    return false;
}

inline bool prior_bias_edge::write(std::ostream& os) const {
    (void)os;
    return false;
}

inline void prior_bias_edge::computeError() {
    const auto bias_vtx = static_cast<const bias_vertex*>(_vertices[0]);
    _error = prior_bias_ - bias_vtx->estimate();
}

inline void prior_bias_edge::linearizeOplus() {
    _jacobianOplusXi.block<3, 3>(0, 0) = Mat33_t::Identity();
}

} // namespace internal
} // namespace imu
} // namespace openvslam

#endif // OPENVSLAM_IMU_INTERNAL_PRIOR_BIAS_EDGE_H
