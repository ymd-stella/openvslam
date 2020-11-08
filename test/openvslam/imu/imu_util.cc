#include "openvslam/imu/imu_util.h"
#include "openvslam/imu/data.h"
#include "openvslam/imu/preintegrated.h"

#include <gtest/gtest.h>

using namespace openvslam;

class keyframe_t1 {
public:
    keyframe_t1(unsigned int id, keyframe_t1* inertial_ref_keyfrm)
        : id_(id), inertial_ref_keyfrm_(inertial_ref_keyfrm) {}
    unsigned int id_;
    keyframe_t1* inertial_ref_keyfrm_;
};

TEST(data, preprocess_imu) {
    {
        // int size = 11;
        {
            Vec3_t a1;
            a1 << 9.3898677825927734, 0.040861040353775024, -3.0482337474822998;
            Vec3_t w1;
            w1 << 0.15358898043632507, -0.063529983162879944, 0.17174039781093597;
            double t1 = 1403636581.8135555;
            Vec3_t a2;
            a2 << 9.4225559234619141, 0.1634441614151001, -2.9011340141296387;
            Vec3_t w2;
            w2 << 0.13892820477485657, -0.050963614135980606, 0.15917402505874634;
            double t2 = 1403636581.8185554;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.4062118530273438, 0.10215260088443756, -2.9746837615966797;
            Vec3_t gyr;
            gyr << 0.14625859260559082, -0.057246796786785126, 0.16545721888542175;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          true, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.4225559234619141, 0.1634441614151001, -2.9011340141296387;
            Vec3_t w1;
            w1 << 0.13892820477485657, -0.050963614135980606, 0.15917402505874634;
            double t1 = 1403636581.8185554;
            Vec3_t a2;
            a2 << 9.4552450180053711, 0.17978858947753906, -2.8357563018798828;
            Vec3_t w2;
            w2 << 0.11100293695926666, -0.050265483558177948, 0.14032447338104248;
            double t2 = 1403636581.8235555;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.4388999938964844, 0.17161637544631958, -2.8684451580047607;
            Vec3_t gyr;
            gyr << 0.12496557086706161, -0.050614550709724426, 0.14974924921989441;
            double dt = 0.0050001144409179688;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.4552450180053711, 0.17978858947753906, -2.8357563018798828;
            Vec3_t w1;
            w1 << 0.11100293695926666, -0.050265483558177948, 0.14032447338104248;
            double t1 = 1403636581.8235555;
            Vec3_t a2;
            a2 << 9.3898677825927734, 0.24516625702381134, -2.7867231369018555;
            Vec3_t w2;
            w2 << 0.084473937749862671, -0.056548666208982468, 0.11938051879405975;
            double t2 = 1403636581.8285556;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.4225559234619141, 0.2124774158000946, -2.8112397193908691;
            Vec3_t gyr;
            gyr << 0.097738437354564667, -0.053407073020935059, 0.12985250353813171;
            double dt = 0.0050001144409179688;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.3898677825927734, 0.24516625702381134, -2.7867231369018555;
            Vec3_t w1;
            w1 << 0.084473937749862671, -0.056548666208982468, 0.11938051879405975;
            double t1 = 1403636581.8285556;
            Vec3_t a2;
            a2 << 9.2345952987670898, 0.2124774158000946, -2.876617431640625;
            Vec3_t w2;
            w2 << 0.05864306166768074, -0.060737457126379013, 0.098436571657657623;
            double t2 = 1403636581.8335555;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.3122310638427734, 0.22882184386253357, -2.8316702842712402;
            Vec3_t gyr;
            gyr << 0.071558497846126556, -0.05864306166768074, 0.10890854895114899;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.2345952987670898, 0.2124774158000946, -2.876617431640625;
            Vec3_t w1;
            w1 << 0.05864306166768074, -0.060737457126379013, 0.098436571657657623;
            double t1 = 1403636581.8335555;
            Vec3_t a2;
            a2 << 9.2509393692016602, 0.1634441614151001, -3.0972669124603271;
            Vec3_t w2;
            w2 << 0.02932153083384037, -0.060039326548576355, 0.085870198905467987;
            double t2 = 1403636581.8385553;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.242767333984375, 0.18796078860759735, -2.9869422912597656;
            Vec3_t gyr;
            gyr << 0.04398229718208313, -0.060388393700122833, 0.092153385281562805;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.2509393692016602, 0.1634441614151001, -3.0972669124603271;
            Vec3_t w1;
            w1 << 0.02932153083384037, -0.060039326548576355, 0.085870198905467987;
            double t1 = 1403636581.8385553;
            Vec3_t a2;
            a2 << 9.2019062042236328, 0.13075533509254456, -3.2607111930847168;
            Vec3_t w2;
            w2 << 0.00069813168374821544, -0.067020639777183533, 0.079587012529373169;
            double t2 = 1403636581.8435557;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.2264232635498047, 0.14709974825382233, -3.1789889335632324;
            Vec3_t gyr;
            gyr << 0.015009831637144089, -0.063529983162879944, 0.082728609442710876;
            double dt = 0.0050003528594970703;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.2019062042236328, 0.13075533509254456, -3.2607111930847168;
            Vec3_t w1;
            w1 << 0.00069813168374821544, -0.067020639777183533, 0.079587012529373169;
            double t1 = 1403636581.8435557;
            Vec3_t a2;
            a2 << 9.1692180633544922, -0.098066501319408417, -3.3832943439483643;
            Vec3_t w2;
            w2 << -0.019547687843441963, -0.074700094759464264, 0.078190751373767853;
            double t2 = 1403636581.8485556;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.1855621337890625, 0.016344416886568069, -3.3220028877258301;
            Vec3_t gyr;
            gyr << -0.0094247777014970779, -0.070860370993614197, 0.078888878226280212;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.1692180633544922, -0.098066501319408417, -3.3832943439483643;
            Vec3_t w1;
            w1 << -0.019547687843441963, -0.074700094759464264, 0.078190751373767853;
            double t1 = 1403636581.8485556;
            Vec3_t a2;
            a2 << 9.12835693359375, -0.39226600527763367, -3.4159832000732422;
            Vec3_t w2;
            w2 << -0.028623400256037712, -0.083775803446769714, 0.082379542291164398;
            double t2 = 1403636581.8535554;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.1487874984741211, -0.24516625702381134, -3.3996386528015137;
            Vec3_t gyr;
            gyr << -0.024085544049739838, -0.079237952828407288, 0.080285146832466125;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.12835693359375, -0.39226600527763367, -3.4159832000732422;
            Vec3_t w1;
            w1 << -0.028623400256037712, -0.083775803446769714, 0.082379542291164398;
            double t1 = 1403636581.8535554;
            Vec3_t a2;
            a2 << 9.1692180633544922, -0.7518431544303894, -3.4813606739044189;
            Vec3_t w2;
            w2 << -0.023038346320390701, -0.10402162373065948, 0.082379542291164398;
            double t2 = 1403636581.8585553;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.1487874984741211, -0.57205456495285034, -3.448671817779541;
            Vec3_t gyr;
            gyr << -0.025830872356891632, -0.0938987135887146, 0.082379542291164398;
            double dt = 0.0049998760223388672;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, false, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
        {
            Vec3_t a1;
            a1 << 9.1692180633544922, -0.7518431544303894, -3.4813606739044189;
            Vec3_t w1;
            w1 << -0.023038346320390701, -0.10402162373065948, 0.082379542291164398;
            double t1 = 1403636581.8585553;
            Vec3_t a2;
            a2 << 9.2100791931152344, -0.95614838600158691, -3.5385661125183105;
            Vec3_t w2;
            w2 << -0.0027925267349928617, -0.11938051879405975, 0.094945907592773438;
            double t2 = 1403636581.8635557;
            double t_prev = 1403636581.8135555;
            double t_cur = 1403636581.8635557;
            Vec3_t acc;
            acc << 9.1896486282348633, -0.85399580001831055, -3.5099635124206543;
            Vec3_t gyr;
            gyr << -0.012915436178445816, -0.11170107126235962, 0.088662728667259216;
            double dt = 0.0050003528594970703;
            Vec3_t acc2, gyr2;
            double dt2;
            imu::imu_util::preprocess_imu(imu::data(a1, w1, t1), imu::data(a2, w2, t2),
                                          false, true, t_prev, t_cur,
                                          acc2, gyr2, dt2);
            EXPECT_LT((acc2 - acc).cwiseAbs().sum(), 1e-6);
            EXPECT_LT((gyr2 - gyr).cwiseAbs().sum(), 1e-6);
            EXPECT_LT(std::abs(dt2 - dt), 1e-6);
        }
    }
}
