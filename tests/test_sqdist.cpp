#include <gtest/gtest.h>
#include "setrot.h"
#include "sqdist.h"
#include <cmath>

using namespace gslib;

TEST(SetrotSqdist, IdentityRotation) {
    // No anisotropy, azimuth=0, dip=0 => should give identity-like rotation
    int maxrot = 1;
    std::vector<double> rotmat(maxrot * 9, 0.0);
    setrot(0.0, 0.0, 0.0, 1.0, 1.0, 0, maxrot, rotmat);

    // Distance between (0,0,0) and (1,0,0) should be 1.0
    double d = sqdist(0, 0, 0, 1, 0, 0, 0, maxrot, rotmat);
    EXPECT_NEAR(d, 1.0, 1e-10);

    // Distance between (0,0,0) and (1,1,1) should be 3.0
    d = sqdist(0, 0, 0, 1, 1, 1, 0, maxrot, rotmat);
    EXPECT_NEAR(d, 3.0, 1e-10);
}

TEST(SetrotSqdist, AnisotropyScaling) {
    // anis1=0.5 means second axis is scaled by 1/0.5=2
    // anis2=0.5 means third axis is scaled by 1/0.5=2
    // azm=0 => major axis is N-S, so x maps to row1 (cos(90°-0)=cos90=0 for x, sin90=1 for y)
    // With azm=0: major axis along Y (North), second along X, third along Z
    int maxrot = 1;
    std::vector<double> rotmat(maxrot * 9, 0.0);
    setrot(0.0, 0.0, 0.0, 0.5, 0.5, 0, maxrot, rotmat);

    // Along major axis (Y/North), distance should be Euclidean
    double d_major = sqdist(0, 0, 0, 0, 1, 0, 0, maxrot, rotmat);
    // Along minor axis (X), distance scaled by 1/anis1=2, so sqd = 4
    double d_minor = sqdist(0, 0, 0, 1, 0, 0, 0, maxrot, rotmat);
    // Along Z, distance scaled by 1/anis2=2, so sqd = 4
    double d_z = sqdist(0, 0, 0, 0, 0, 1, 0, maxrot, rotmat);

    EXPECT_NEAR(d_major, 1.0, 1e-10);
    EXPECT_NEAR(d_minor, 4.0, 1e-10);
    EXPECT_NEAR(d_z, 4.0, 1e-10);
}

TEST(SetrotSqdist, ZeroDistance) {
    int maxrot = 1;
    std::vector<double> rotmat(maxrot * 9, 0.0);
    setrot(45.0, 0.0, 0.0, 1.0, 1.0, 0, maxrot, rotmat);

    double d = sqdist(5, 5, 5, 5, 5, 5, 0, maxrot, rotmat);
    EXPECT_NEAR(d, 0.0, 1e-15);
}
