#include <gtest/gtest.h>
#include "gamv.h"
#include <cmath>
#include <random>

using namespace gslib;

TEST(Gamv, ConstantValues) {
    // Constant values => semivariance should be 0
    int n = 100;
    std::vector<double> x(n), y(n), z(n, 0.0), values(n, 5.0);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i);
        y[i] = 0.0;
    }

    auto result = gamv(x, y, z, values, 10, 5.0, -1.0, 0.0, 90.0);

    for (int i = 0; i < 10; i++) {
        if (result.pair_counts[i] > 0) {
            EXPECT_NEAR(result.semivariance[i], 0.0, 1e-10);
        }
    }
}

TEST(Gamv, LinearTrend) {
    // Linear values along x: values = x
    // Semivariogram should increase with lag
    int n = 100;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), values(n);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i);
        values[i] = static_cast<double>(i);
    }

    auto result = gamv(x, y, z, values, 10, 5.0, -1.0, 0.0, 90.0);

    // Semivariance should be monotonically increasing
    for (int i = 1; i < 10; i++) {
        if (result.pair_counts[i] > 0 && result.pair_counts[i - 1] > 0) {
            EXPECT_GE(result.semivariance[i], result.semivariance[i - 1] - 1e-10);
        }
    }
}

TEST(Gamv, PairCounts) {
    // 1D regular grid: n points, lag=1
    int n = 50;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), values(n, 1.0);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i);
    }

    // Use xlag=1.0 and xltol=0.5 so pairs at distance 1 fall in first lag bin
    auto result = gamv(x, y, z, values, 5, 1.0, 0.5, 0.0, 90.0);

    // Some lag bins should have pairs
    int total_pairs = 0;
    for (int i = 0; i < 5; i++) total_pairs += result.pair_counts[i];
    EXPECT_GT(total_pairs, 0);
}

TEST(Gamv, EmptyInput) {
    std::vector<double> x, y, z, values;
    auto result = gamv(x, y, z, values);
    EXPECT_TRUE(result.lags.empty());
}

TEST(Gamv, DownholeConstraint) {
    // Two drillholes, pairs should only be within same hole
    int n = 20;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), values(n, 1.0);
    std::vector<int> bhid(n);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i);
        bhid[i] = (i < 10) ? 1 : 2;
    }

    auto result_no_bhid = gamv(x, y, z, values, 5, 1.0, 0.5, 0.0, 90.0);
    auto result_bhid = gamv(x, y, z, values, 5, 1.0, 0.5, 0.0, 90.0, 1e10, 0.0, 90.0, 1e10, bhid);

    // With bhid constraint, should have fewer pairs at larger lags
    // (cross-hole pairs are excluded)
    EXPECT_LE(result_bhid.pair_counts[4], result_no_bhid.pair_counts[4]);
}
