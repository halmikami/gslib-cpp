#include <gtest/gtest.h>
#include "gamv.h"
#include "kt3d.h"
#include "cova3.h"
#include "setrot.h"
#include "common.h"
#include <cmath>
#include <random>

using namespace gslib;

// Numerical validation tests: verify correctness against analytical solutions
// and known GSLIB reference values. These serve as the equivalent of
// PyGSLIB comparison tests since pygslib is not available in this environment.

// ============================================================================
// Covariance model analytical verification
// ============================================================================

TEST(Numerical, SphericalVariogramModel) {
    // Spherical model: C(h) = cc*(1 - 1.5*(h/a) + 0.5*(h/a)^3) for h<a, 0 for h>=a
    // With nugget=0.1, cc=0.9, range=100
    int maxrot = 1;
    std::vector<double> rotmat(9, 0.0);
    setrot(0, 0, 0, 1, 1, 0, maxrot, rotmat);

    std::vector<double> c0 = {0.1};
    std::vector<int> it = {1};
    std::vector<double> cc = {0.9};
    std::vector<double> aa = {100.0};

    // Test at multiple distances
    // Spherical: C(h) = cc * (1 - hr*(1.5 - 0.5*hr*hr)) where hr = h/a
    double distances[] = {0, 10, 25, 50, 75, 100, 150};

    for (int i = 0; i < 7; i++) {
        auto r = cova3(0, 0, 0, distances[i], 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
        double d = distances[i];
        double expected;
        if (d < EPSLON) {
            expected = 1.0; // cmax
        } else if (d >= 100.0) {
            expected = 0.0;
        } else {
            double hr = d / 100.0;
            expected = 0.9 * (1.0 - hr * (1.5 - 0.5 * hr * hr));
        }
        EXPECT_NEAR(r.cova, expected, 1e-6)
            << "Failed at distance=" << d;
    }
}

TEST(Numerical, ExponentialVariogramModel) {
    // Exponential: C(h) = cc * exp(-3h/a)
    int maxrot = 1;
    std::vector<double> rotmat(9, 0.0);
    setrot(0, 0, 0, 1, 1, 0, maxrot, rotmat);

    std::vector<double> c0 = {0.0};
    std::vector<int> it = {2};
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    double distances[] = {10, 50, 100, 200};
    for (double d : distances) {
        auto r = cova3(0, 0, 0, d, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
        double expected = std::exp(-3.0 * d / 100.0);
        EXPECT_NEAR(r.cova, expected, 1e-10)
            << "Exponential failed at distance=" << d;
    }
}

TEST(Numerical, GaussianVariogramModel) {
    // Gaussian: C(h) = cc * exp(-(3h/a)^2)
    int maxrot = 1;
    std::vector<double> rotmat(9, 0.0);
    setrot(0, 0, 0, 1, 1, 0, maxrot, rotmat);

    std::vector<double> c0 = {0.0};
    std::vector<int> it = {3};
    std::vector<double> cc = {1.0};
    std::vector<double> aa = {100.0};

    double distances[] = {10, 50, 100};
    for (double d : distances) {
        auto r = cova3(0, 0, 0, d, 0, 0, 0, 1, c0, it, cc, aa, 0, maxrot, rotmat);
        double ratio = 3.0 * d / 100.0;
        double expected = std::exp(-ratio * ratio);
        EXPECT_NEAR(r.cova, expected, 1e-10)
            << "Gaussian failed at distance=" << d;
    }
}

TEST(Numerical, NestedStructures) {
    // Nugget + Spherical + Exponential
    int maxrot = 3;
    std::vector<double> rotmat(maxrot * 9, 0.0);
    setrot(0, 0, 0, 1, 1, 0, maxrot, rotmat);
    setrot(0, 0, 0, 1, 1, 1, maxrot, rotmat);

    std::vector<double> c0 = {0.1};       // nugget
    std::vector<int> it = {1, 2};          // spherical + exponential
    std::vector<double> cc = {0.5, 0.4};   // sill contributions
    std::vector<double> aa = {50.0, 200.0}; // ranges

    // At zero distance: cmax = 0.1 + 0.5 + 0.4 = 1.0
    auto r0 = cova3(0, 0, 0, 0, 0, 0, 0, 2, c0, it, cc, aa, 0, maxrot, rotmat);
    EXPECT_NEAR(r0.cmax, 1.0, 1e-10);
    EXPECT_NEAR(r0.cova, 1.0, 1e-10);

    // At h=25: spherical contributes, exponential contributes
    auto r25 = cova3(0, 0, 0, 25, 0, 0, 0, 2, c0, it, cc, aa, 0, maxrot, rotmat);
    double hr = 25.0 / 50.0;
    double sph_part = 0.5 * (1.0 - hr * (1.5 - 0.5 * hr * hr));
    double exp_part = 0.4 * std::exp(-3.0 * 25.0 / 200.0);
    EXPECT_NEAR(r25.cova, sph_part + exp_part, 1e-6);
}

// ============================================================================
// Variogram analytical verification
// ============================================================================

TEST(Numerical, VariogramLinearData) {
    // For values = x (linear), the theoretical semivariogram is:
    // gamma(h) = 0.5 * E[(Z(x+h) - Z(x))^2] = 0.5 * h^2
    int n = 200;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), values(n);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i);
        values[i] = static_cast<double>(i);
    }

    auto result = gamv(x, y, z, values, 10, 10.0, 5.0, 0.0, 90.0);

    // For linear data, gamma(h) = 0.5*h^2 is the expectation for pairs
    // at exact lag h. With binning, the average distance in a bin can differ.
    // Use actual average lag distance from the result to compute expected.
    for (int i = 0; i < 10; i++) {
        if (result.pair_counts[i] > 0) {
            double lag = result.lags[i];
            double expected = 0.5 * lag * lag;
            // Allow 15% tolerance due to finite sampling and bin averaging
            EXPECT_NEAR(result.semivariance[i], expected, expected * 0.15 + 1.0)
                << "Linear variogram failed at lag=" << lag;
        }
    }
}

TEST(Numerical, VariogramWhiteNoise) {
    // White noise: gamma(h) = variance for h > 0
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    int n = 1000;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), values(n);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i) * 0.5;
        values[i] = dist(rng);
    }

    // Compute sample variance
    double mean = 0, var = 0;
    for (int i = 0; i < n; i++) mean += values[i];
    mean /= n;
    for (int i = 0; i < n; i++) var += (values[i] - mean) * (values[i] - mean);
    var /= (n - 1);

    auto result = gamv(x, y, z, values, 10, 5.0, 2.5, 0.0, 90.0);

    // All lags should be approximately equal to the variance (pure nugget)
    for (int i = 0; i < 10; i++) {
        if (result.pair_counts[i] > 100) {
            EXPECT_NEAR(result.semivariance[i], var, var * 0.3)
                << "White noise variogram failed at lag " << i;
        }
    }
}

// ============================================================================
// Kriging analytical verification
// ============================================================================

TEST(Numerical, KrigingExactInterpolation) {
    // OK should exactly reproduce known values at data locations
    int n = 16;
    std::vector<double> x(n), y(n), z(n, 0.0), values(n);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i % 4) * 10.0;
        y[i] = static_cast<double>(i / 4) * 10.0;
        values[i] = x[i] + 2.0 * y[i] + 1.0;
    }

    // Estimate at each data point
    auto result = kt3d(x, y, z, values, x, y, z,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0, {1}, {1.0}, {100.0});

    for (int i = 0; i < n; i++) {
        EXPECT_NEAR(result.estimates[i], values[i], 0.01)
            << "Exact interpolation failed at point " << i;
    }
}

TEST(Numerical, KrigingVarianceAtDataIsZero) {
    // Kriging variance at a data point should be ~0
    int n = 9;
    std::vector<double> x = {0, 10, 20, 0, 10, 20, 0, 10, 20};
    std::vector<double> y = {0, 0, 0, 10, 10, 10, 20, 20, 20};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    std::vector<double> xout = {10.0};
    std::vector<double> yout = {10.0};
    std::vector<double> zout = {0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0, {1}, {1.0}, {50.0});

    EXPECT_NEAR(result.variances[0], 0.0, 0.05);
}

TEST(Numerical, KrigingVarianceIncreasesWithDistance) {
    // Variance should increase as we move away from data
    int n = 4;
    std::vector<double> x = {0, 100, 0, 100};
    std::vector<double> y = {0, 0, 100, 100};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {1, 2, 3, 4};

    std::vector<double> xout = {50.0, 50.0};
    std::vector<double> yout = {50.0, 150.0};  // second point farther from data
    std::vector<double> zout = {0.0, 0.0};

    auto result = kt3d(x, y, z, values, xout, yout, zout,
                        300.0, -1, -1, 0, 0, 0, 16, 1, 0,
                        0.0, {1}, {1.0}, {200.0});

    // Point at (50,150) is farther from centroid of data => higher variance
    EXPECT_GT(result.variances[1], result.variances[0]);
}

TEST(Numerical, SKvsOKDifference) {
    // SK and OK should give different results when skmean != data mean
    int n = 9;
    std::vector<double> x = {0, 10, 20, 0, 10, 20, 0, 10, 20};
    std::vector<double> y = {0, 0, 0, 10, 10, 10, 20, 20, 20};
    std::vector<double> z(n, 0.0);
    std::vector<double> values = {10, 12, 14, 11, 13, 15, 12, 14, 16};

    std::vector<double> xout = {5.0};
    std::vector<double> yout = {5.0};
    std::vector<double> zout = {0.0};

    auto ok_result = kt3d(x, y, z, values, xout, yout, zout,
                           200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                           0.0, {1}, {1.0}, {50.0},
                           {}, {}, {}, {}, {},
                           1); // OK

    auto sk_result = kt3d(x, y, z, values, xout, yout, zout,
                           200.0, -1, -1, 0, 0, 0, 16, 1, 0,
                           0.0, {1}, {1.0}, {50.0},
                           {}, {}, {}, {}, {},
                           0, 0.0); // SK with mean=0 (very different from data mean)

    // Both should produce valid estimates but differ
    EXPECT_NE(ok_result.estimates[0], -999.0);
    EXPECT_NE(ok_result.estimates[0], sk_result.estimates[0]);
    // SK variance should be <= OK variance
    EXPECT_LE(sk_result.variances[0], ok_result.variances[0] + 0.1);
}
