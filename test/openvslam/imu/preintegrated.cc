#include "openvslam/imu/preintegrated.h"

#include <gtest/gtest.h>

using namespace openvslam;

TEST(data, reintegrate) {
    // basic information
    const std::string name = "IMU";
    const double rate_hz = 200.0;

    // create the relative pose "from IMU to camera" (_ic) and "from camera to IMU" (_ci)
    Mat44_t rel_pose_ic;
    rel_pose_ic << 0.0148655429818, -0.999880929698, 0.00414029679422, -0.0216401454975,
        0.999557249008, 0.0149672133247, 0.025715529948, -0.064676986768,
        -0.0257744366974, 0.00375618835797, 0.999660727178, 0.00981073058949,
        0.0, 0.0, 0.0, 1.0;

    // acc noise [(m/s^2)/sqrt(Hz)]
    const double ns_acc = 2.0000e-3;
    // gyr noise [(rad/s]/sqrt(Hz)]
    const double ns_gyr = 1.7e-4;
    // acc bias random walk [(m/s^3)/sqrt(Hz)]
    const double rw_acc_bias = 3.0000e-03;
    // gyr bias random walk [(rad/s^2]/sqrt(Hz)]
    const double rw_gyr_bias = 1.9393e-05;

    const auto cfg = eigen_alloc_shared<imu::config>(name, rate_hz, rel_pose_ic, ns_acc, ns_gyr, rw_acc_bias, rw_gyr_bias);

    {
        imu::preintegrator p(imu::bias(), cfg);
        double dT = 0.25;
        MatRC_t<15, 15> C;
        C << 7.225e-09, -1.27462e-16, 6.0359e-17, -7.31754e-11, 4.33882e-09, -1.67116e-10, -6.12624e-12, 3.48724e-10, -9.93679e-12, 0, 0, 0, 0, 0, 0, -1.27394e-16, 7.225e-09, -5.36869e-16, -3.69204e-09, -2.7057e-11, -9.38122e-09, -2.96494e-10, -2.15074e-12, -7.5784e-10, 0, 0, 0, 0, 0, 0, 5.89788e-17, -5.36918e-16, 7.225e-09, 7.40036e-12, 9.10212e-09, 5.82739e-11, -2.97597e-12, 7.35347e-10, 4.45787e-12, 0, 0, 0, 0, 0, 0, -7.31754e-11, -3.69204e-09, 7.40037e-12, 1.00251e-06, -4.61832e-11, 6.38393e-09, 1.25226e-07, -4.83166e-12, 5.78846e-10, 0, 0, 0, 0, 0, 0, 4.33882e-09, -2.70571e-11, 9.10212e-09, -4.61832e-11, 1.01875e-06, 1.81194e-11, -9.6421e-12, 1.267e-07, 3.78773e-12, 0, 0, 0, 0, 0, 0, -1.67116e-10, -9.38121e-09, 5.82739e-11, 6.38393e-09, 1.81194e-11, 1.01624e-06, 5.76379e-10, 1.88945e-12, 1.26474e-07, 0, 0, 0, 0, 0, 0, -6.12624e-12, -2.96494e-10, -2.97597e-12, 1.25226e-07, -9.64209e-12, 5.76379e-10, 2.08531e-08, -9.14081e-13, 5.58616e-11, 0, 0, 0, 0, 0, 0, 3.48724e-10, -2.15074e-12, 7.35346e-10, -4.83166e-12, 1.267e-07, 1.88946e-12, -9.14081e-13, 2.09959e-08, 3.57514e-13, 0, 0, 0, 0, 0, 0, -9.93678e-12, -7.5784e-10, 4.45787e-12, 5.78846e-10, 3.78774e-12, 1.26474e-07, 5.58616e-11, 3.57514e-13, 2.09741e-08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.76089e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.76089e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.76089e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09;
        MatRC_t<6, 6> Nga;
        Nga << 5.78e-06, 0, 0, 0, 0, 0, 0, 5.78e-06, 0, 0, 0, 0, 0, 0, 5.78e-06, 0, 0, 0, 0, 0, 0, 0.0008, 0, 0, 0, 0, 0, 0, 0.0008, 0, 0, 0, 0, 0, 0, 0.0008;
        MatRC_t<6, 6> NgaWalk;
        NgaWalk << 7.52177e-08, 0, 0, 0, 0, 0, 0, 7.52177e-08, 0, 0, 0, 0, 0, 0, 7.52177e-08, 0, 0, 0, 0, 0, 0, 0.0018, 0, 0, 0, 0, 0, 0, 0.0018, 0, 0, 0, 0, 0, 0, 0.0018;
        Mat33_t dR;
        dR << 0.997363, -0.0199074, -0.0697836, 0.0195129, 0.99979, -0.00632964, 0.0698949, 0.00495126, 0.997542;
        Vec3_t dV;
        dV << 2.60269, 0.0085997, -1.02248;
        Vec3_t dP;
        dP << 0.319498, 0.00157434, -0.125288;
        Mat33_t JRg;
        JRg << -0.249855, -0.000426387, -0.00647477, 0.000439231, -0.249993, -0.000361433, 0.0064727, 0.000405747, -0.249861;
        Mat33_t JVg;
        JVg << 0.00188871, 0.12776, 4.47947e-05, -0.138147, 0.000815634, -0.320308, 0.00423786, 0.324622, -0.00129467;
        Mat33_t JVa;
        JVa << -0.249664, 0.00447337, 0.0108113, -0.00442846, -0.249952, 0.00107803, -0.0108299, -0.000850437, -0.249708;
        Mat33_t JPg;
        JPg << 0.000129729, 0.0102601, 0.000131011, -0.0109212, 4.29987e-05, -0.0259519, 0.000152358, 0.0262243, -9.10077e-05;
        Mat33_t JPa;
        JPa << -0.0312242, 0.000453093, 0.000978035, -0.000448955, -0.0312455, 0.000121724, -0.00097996, -0.000102707, -0.0312282;
        {
            Vec3_t a;
            a << 10.0273, 0.040861, -4.05342;
            Vec3_t w;
            w << 0.0959931, -0.440521, 0.242252;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.95375, 0.077636, -3.97987;
            Vec3_t w;
            w << 0.0994838, -0.444012, 0.233874;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.97827, 0.0245166, -3.95535;
            Vec3_t w;
            w << 0.0949459, -0.449946, 0.219213;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0069, -0.0163444, -3.95943;
            Vec3_t w;
            w << 0.0921534, -0.455182, 0.211883;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0069, -0.0817221, -3.94718;
            Vec3_t w;
            w << 0.0907571, -0.453786, 0.211185;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0314, -0.155272, -3.94309;
            Vec3_t w;
            w << 0.0858702, -0.450993, 0.211534;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0436, -0.143014, -4.01255;
            Vec3_t w;
            w << 0.0816814, -0.445059, 0.214676;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.011, -0.106239, -4.13105;
            Vec3_t w;
            w << 0.0771435, -0.432493, 0.223751;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.9987, -0.0898943, -4.22912;
            Vec3_t w;
            w << 0.0694641, -0.41853, 0.235619;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.96192, -0.0572055, -4.31493;
            Vec3_t w;
            w << 0.0586431, -0.410152, 0.244695;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.9088, -0.0490333, -4.38848;
            Vec3_t w;
            w << 0.0460767, -0.401077, 0.252026;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.9088, -0.0898943, -4.42934;
            Vec3_t w;
            w << 0.0335103, -0.395841, 0.25412;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.97009, -0.110325, -4.47837;
            Vec3_t w;
            w << 0.0226893, -0.400728, 0.241903;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0559, -0.118497, -4.47428;
            Vec3_t w;
            w << 0.0129154, -0.402822, 0.219213;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.154, -0.16753, -4.38439;
            Vec3_t w;
            w << 0.00453786, -0.398633, 0.196524;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.2725, -0.208391, -4.31493;
            Vec3_t w;
            w << 0, -0.392001, 0.175231;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.2847, -0.175702, -4.2332;
            Vec3_t w;
            w << -0.00244346, -0.386416, 0.15324;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.2766, -0.143014, -4.15557;
            Vec3_t w;
            w << -0.00383972, -0.376293, 0.131947;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.2929, -0.130755, -4.13105;
            Vec3_t w;
            w << -0.00244346, -0.361632, 0.114145;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.1948, -0.130755, -4.09836;
            Vec3_t w;
            w << 0.000698132, -0.347321, 0.101578;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.0559, -0.224736, -4.08202;
            Vec3_t w;
            w << 0.00314159, -0.333707, 0.0945968;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.93332, -0.310544, -4.1474;
            Vec3_t w;
            w << 0.00698132, -0.322537, 0.0928515;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.81074, -0.306458, -4.23729;
            Vec3_t w;
            w << 0.00907571, -0.314159, 0.095644;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.74536, -0.290113, -4.32718;
            Vec3_t w;
            w << 0.00837758, -0.301244, 0.101927;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.77396, -0.245166, -4.45385;
            Vec3_t w;
            w << 0.00698132, -0.28414, 0.102974;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.83117, -0.118497, -4.56826;
            Vec3_t w;
            w << 0.0010472, -0.270875, 0.0977384;
            double t = 0.00500035;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.84342, 0.040861, -4.63773;
            Vec3_t w;
            w << -0.00733038, -0.260752, 0.090059;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.85568, 0.122583, -4.69902;
            Vec3_t w;
            w << -0.0181514, -0.252026, 0.0788889;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 9.96601, 0.143014, -4.76031;
            Vec3_t w;
            w << -0.0303687, -0.245393, 0.0649262;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.1008, 0.114411, -4.78074;
            Vec3_t w;
            w << -0.0411898, -0.239459, 0.047473;
            double t = 0.00500035;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.2153, 0.0653777, -4.76848;
            Vec3_t w;
            w << -0.0502655, -0.231082, 0.0244346;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.3542, 0.0367749, -4.78483;
            Vec3_t w;
            w << -0.0548033, -0.223053, -0.00279253;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.4523, 0.0367749, -4.78483;
            Vec3_t w;
            w << -0.0551524, -0.214326, -0.0289725;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5054, 0, -4.72762;
            Vec3_t w;
            w << -0.0544543, -0.203505, -0.0499164;
            double t = 0.00500035;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5503, -0.126669, -4.71128;
            Vec3_t w;
            w << -0.0527089, -0.19059, -0.0659734;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5667, -0.236994, -4.7358;
            Vec3_t w;
            w << -0.0457276, -0.177675, -0.0747001;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5544, -0.384094, -4.75214;
            Vec3_t w;
            w << -0.0338594, -0.166853, -0.0747001;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5503, -0.53528, -4.79709;
            Vec3_t w;
            w << -0.0171042, -0.155683, -0.0684169;
            double t = 0.00500035;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5421, -0.539366, -4.84203;
            Vec3_t w;
            w << 0.00139626, -0.149051, -0.058294;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5176, -0.478074, -4.87472;
            Vec3_t w;
            w << 0.0188496, -0.141372, -0.0488692;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.4768, -0.416783, -4.95236;
            Vec3_t w;
            w << 0.0328122, -0.131947, -0.0415388;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.4563, -0.343233, -5.03408;
            Vec3_t w;
            w << 0.042237, -0.12706, -0.03735;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.489, -0.306458, -5.09129;
            Vec3_t w;
            w << 0.047473, -0.119381, -0.037001;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5299, -0.302372, -5.14849;
            Vec3_t w;
            w << 0.0481711, -0.109607, -0.0408407;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5789, -0.277855, -5.17301;
            Vec3_t w;
            w << 0.0485202, -0.104022, -0.0485202;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.5953, -0.216564, -5.15258;
            Vec3_t w;
            w << 0.0492183, -0.0994838, -0.0568977;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.583, -0.24108, -5.13623;
            Vec3_t w;
            w << 0.0495673, -0.0938987, -0.0649263;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.6443, -0.339147, -5.14032;
            Vec3_t w;
            w << 0.0527089, -0.0893609, -0.074002;
            double t = 0.00499988;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.7015, -0.400438, -5.1158;
            Vec3_t w;
            w << 0.0593412, -0.0837758, -0.0767945;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        {
            Vec3_t a;
            a << 10.7546, -0.457644, -5.0586;
            Vec3_t w;
            w << 0.0659734, -0.0771435, -0.069115;
            double t = 0.00500011;
            p.measurements_.emplace_back(a, w, t);
        }
        p.reintegrate();
        EXPECT_LT((p.preintegrated_->jacob_rotation_gyr_ - JRg).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->jacob_velocity_gyr_ - JVg).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->jacob_velocity_acc_ - JVa).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->jacob_position_gyr_ - JPg).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->jacob_position_acc_ - JPa).cwiseAbs().sum(), 1e-5);
        EXPECT_LT(std::abs(p.preintegrated_->dt_ - dT), 1e-5);
        EXPECT_LT((p.preintegrated_->delta_rotation_ - dR).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->delta_velocity_ - dV).cwiseAbs().sum(), 1e-5);
        EXPECT_LT((p.preintegrated_->delta_position_ - dP).cwiseAbs().sum(), 1e-5);
    }
}
