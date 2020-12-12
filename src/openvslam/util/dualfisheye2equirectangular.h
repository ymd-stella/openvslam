#ifndef DUALFISHEYE_TO_EQUIRECTANGULER_H
#define DUALFISHEYE_TO_EQUIRECTANGULER_H

#include <openvslam/config.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <memory>

namespace openvslam {
namespace util {

class dualfisheye2equirectangular {
public:
    dualfisheye2equirectangular(const std::shared_ptr<config>& cfg);

    void convert(const cv::Mat& in_img, cv::Mat& out_img) const;

    cv::Mat map_x_, map_y_;
};

} // namespace util
} // namespace openvslam

#endif // DUALFISHEYE_TO_EQUIRECTANGULER_H
