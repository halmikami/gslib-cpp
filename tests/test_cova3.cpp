#include <gtest/gtest.h>
#include "cova3.h"
#include "setrot.h"
#include "common.h"
#include <cmath>

using namespace gslib;

class Cova3Test : public ::testing::Test {
protected:
    int maxrot = 2;
    std::vector<double> rotmat;

    void SetUp() override {
        rotmat.resize(maxrot * 9, 0.0);
        setrot(0.0, 0.0, 0.0, 1.0, 1.0, 0, maxrot, rotmat);
    }
};

TEST_F(Cova3Test, SphericalZeroDistance) {
    std::vector<double> c0 = {0.1};  // nugget
    std::vector<int> it = {1};       // spherical
    std::vector<double> cc = {0.9};  // sill contribution
    std::vector<double> aa = {100.0}; // range

    auto r = cova3(0, 0, 0, 0, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    EXPECT_NEAR(r.cmax, 1.0, 1e-10);
    EXPECT_NEAR(r.cova, 1.0, 1e-10); // cmax at zero distance
}

TEST_F(Cova3Test, SphericalBeyondRange) {
    std::vector<double> c0 = {0.0};
    std::vector<int> it = {1};
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {50.0};

    auto r = cova3(0, 0, 0, 100, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    EXPECT_NEAR(r.cova, 0.0, 1e-10);
}

TEST_F(Cova3Test, SphericalHalfRange) {
    std::vector<double> c0 = {0.0};
    std::vector<int> it = {1};
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    // At h=50, hr=0.5: C = 1*(1 - 0.5*(1.5 - 0.5*0.25)) = 1*(1-0.5*1.375) = 0.3125
    auto r = cova3(0, 0, 0, 50, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    double expected = 1.0 * (1.0 - 0.5 * (1.5 - 0.5 * 0.25));
    EXPECT_NEAR(r.cova, expected, 1e-6);
}

TEST_F(Cova3Test, ExponentialModel) {
    std::vector<double> c0 = {0.0};
    std::vector<int> it = {2};  // exponential
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    auto r = cova3(0, 0, 0, 50, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    double expected = std::exp(-3.0 * 50.0 / 100.0);
    EXPECT_NEAR(r.cova, expected, 1e-6);
}

TEST_F(Cova3Test, GaussianModel) {
    std::vector<double> c0 = {0.0};
    std::vector<int> it = {3};  // gaussian
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    auto r = cova3(0, 0, 0, 50, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    double ratio = 3.0 * 50.0 / 100.0;
    double expected = std::exp(-ratio * ratio);
    EXPECT_NEAR(r.cova, expected, 1e-6);
}

TEST_F(Cova3Test, HoleEffectModel) {
    std::vector<double> c0 = {0.0};
    std::vector<int> it = {5};  // hole effect
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    auto r = cova3(0, 0, 0, 50, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
    double expected = std::cos(50.0 / 100.0 * gslib::PI);
    EXPECT_NEAR(r.cova, expected, 1e-6);
}
