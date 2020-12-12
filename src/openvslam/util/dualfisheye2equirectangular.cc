#include "openvslam/util/dualfisheye2equirectangular.h"

namespace openvslam {
namespace util {

dualfisheye2equirectangular::dualfisheye2equirectangular(const std::shared_ptr<config>& cfg) {
    const float radius = cfg->yaml_node_["Camera.radius"].as<double>();
    const int cols = cfg->camera_->cols_;
    const int rows = cfg->camera_->rows_;
    map_x_ = cv::Mat(cv::Size(cols, rows), CV_32FC1);
    map_y_ = cv::Mat(cv::Size(cols, rows), CV_32FC1);
    const int dst_rows = cols / 2;
    const float src_cx = dst_rows / 2.0;
    const float src_cy = dst_rows / 2.0;
    const float src_cx2 = src_cx * 3.0;

    // make mapping table
    for(int y = 0; y < dst_rows; y++) {
        for(int x = 0; x < cols; x++) {
            const float ux = x / (float)cols - 0.5;
            const float uy = - (y / (float)dst_rows - 0.5);
            const float longitude = 2 * M_PI * ux;
            const float latitude = M_PI * uy;

            // convert to bearing vector
            const float px = cos(latitude) * sin(longitude);
            const float py = -sin(latitude);
            const float pz = cos(latitude) * cos(longitude);

            // xy plane to bearing vector
            const float phi = atan2(sqrt(px * px + py * py), pz);
            float phi2;
            if (phi > M_PI / 2) {
                phi2 = M_PI - phi;
            }
            else {
                phi2 = phi;
            }
            const float r = 2.0 * phi2 / M_PI * radius;
            const float th = atan2(py, px);

            float mx, my;
            if (0 > pz) {
                mx = r * cos(-th + M_PI / 2) + src_cx;
                my = r * sin(-th + M_PI / 2) + src_cy;
            }
            else {
                mx = -r * cos(th - M_PI / 2) + src_cx2;
                my = -r * sin(th - M_PI / 2) + src_cy;
            }
            map_x_.at<float>(y, x) = mx;
            map_y_.at<float>(y, x) = my;
        }
    }
}

void dualfisheye2equirectangular::convert(const cv::Mat& in_img, cv::Mat& out_img) const {
    remap(in_img, out_img, map_x_, map_y_, cv::INTER_LINEAR);
}

} // namespace util
} // namespace openvslam
