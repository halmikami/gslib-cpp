#include <gtest/gtest.h>
#include "kt3d.h"
#include <cmath>

using namespace gslib;

TEST(Kt3d, ConstantField) {
    // If all values are the same, kriging should return that value
    int n = 25;
    std::vector<double> x(n), y(n), z(n, 0.0), values(n, 10.0);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i % 5) * 10.0;
        y[i] = static_cast<double>(i / 5) * 10.0;
    }

    std::vector<double> xout = {25.0};
    std::vector<double> yout = {25.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0,
                        {1}, {1.0}, {100.0});

    EXPECT_NEAR(result.estimates[0], 10.0, 1e-4);
}

TEST(Kt3d, SimpleKriging) {
    int n = 9;
    std::vector<double> x = {0, 10, 20, 0, 10, 20, 0, 10, 20};
    std::vector<double> y = {0, 0, 0, 10, 10, 10, 20, 20, 20};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {1, 2, 3, 2, 3, 4, 3, 4, 5};

    std::vector<double> xout = {10.0};
    std::vector<double> yout = {10.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0,
                        {1}, {1.0}, {50.0},
                        {}, {}, {}, {}, {},
                        0,   // SK
                        3.0  // skmean
                        );

    EXPECT_NE(result.estimates[0], -999.0);
    EXPECT_GT(result.estimates[0], 0.0);
}

TEST(Kt3d, OrdinaryKriging) {
    int n = 9;
    std::vector<double> x = {0, 10, 20, 0, 10, 20, 0, 10, 20};
    std::vector<double> y = {0, 0, 0, 10, 10, 10, 20, 20, 20};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {1, 2, 3, 2, 3, 4, 3, 4, 5};

    std::vector<double> xout = {10.0};
    std::vector<double> yout = {10.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0,
                        {1}, {1.0}, {50.0},
                        {}, {}, {}, {}, {},
                        1 // OK
                        );

    // At (10,10), surrounded by known data, estimate should be close to 3.0
    EXPECT_NEAR(result.estimates[0], 3.0, 0.5);
    EXPECT_GE(result.variances[0], 0.0);
}

TEST(Kt3d, InsufficientData) {
    // Point far from data with small search radius
    std::vector<double> x = {0};
    std::vector<double> y = {0};
    std::vector<double> z = {0};
    std::vector<double> values = {5.0};

    std::vector<double> xout = {1000.0};
    std::vector<double> yout = {1000.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        10.0, -1, -1, 0, 0, 0, 16, 4, 0,
                        0.0, {1}, {1.0}, {50.0});

    // Should be unestimated
    EXPECT_EQ(result.estimates[0], -999.0);
}

TEST(Kt3d, ExactInterpolation) {
    // OK at a known data location should return exact value
    int n = 4;
    std::vector<double> x = {0, 10, 0, 10};
    std::vector<double> y = {0, 0, 10, 10};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {1, 2, 3, 4};

    std::vector<double> xout = {0.0};
    std::vector<double> yout = {0.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0, {1}, {1.0}, {50.0});

    EXPECT_NEAR(result.estimates[0], 1.0, 0.1);
}
