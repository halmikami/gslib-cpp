#include <gtest/gtest.h>
#include "declus.h"
#include <cmath>
#include <numeric>

using namespace gslib;

TEST(Declus, UniformGrid) {
    // Regular grid: all weights should be approximately 1.0
    int n = 25;
    std::vector<double> x(n), y(n), z(n, 0.0), vr(n, 10.0);
    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i % 5) * 10.0;
        y[i] = static_cast<double>(i / 5) * 10.0;
    }

    auto result = declus(x, y, z, vr, 1.0, 1.0, 0, 10, 5.0, 50.0, 5);

    // Weights should average to 1.0 (sum = n)
    double sum = std::accumulate(result.weights.begin(), result.weights.end(), 0.0);
    EXPECT_NEAR(sum, static_cast<double>(n), 1e-6);
}

TEST(Declus, ClusteredData) {
    // Create clustered data: many points near (0,0), few far away
    std::vector<double> x, y, z, vr;

    // Cluster of 10 points near origin with value 1
    for (int i = 0; i < 10; i++) {
        x.push_back(static_cast<double>(i) * 0.1);
        y.push_back(0.0);
        z.push_back(0.0);
        vr.push_back(1.0);
    }
    // 3 isolated points with value 10
    x.push_back(100.0); y.push_back(0.0); z.push_back(0.0); vr.push_back(10.0);
    x.push_back(200.0); y.push_back(0.0); z.push_back(0.0); vr.push_back(10.0);
    x.push_back(300.0); y.push_back(0.0); z.push_back(0.0); vr.push_back(10.0);

    auto result = declus(x, y, z, vr, 1.0, 1.0, 1, 10, 10.0, 200.0, 5);

    // Isolated points should get higher weights
    double avg_cluster_wt = 0;
    for (int i = 0; i < 10; i++) avg_cluster_wt += result.weights[i];
    avg_cluster_wt /= 10.0;

    double avg_isolated_wt = (result.weights[10] + result.weights[11] + result.weights[12]) / 3.0;

    EXPECT_GT(avg_isolated_wt, avg_cluster_wt);
}

TEST(Declus, SingleCellSize) {
    int n = 10;
    std::vector<double> x(n), y(n, 0.0), z(n, 0.0), vr(n, 5.0);
    for (int i = 0; i < n; i++) x[i] = static_cast<double>(i) * 10.0;

    auto result = declus(x, y, z, vr, 1.0, 1.0, 0, 1, 20.0, 20.0, 5);

    EXPECT_EQ(result.cell_sizes.size(), 2u);
    EXPECT_NEAR(result.declustered_mean, 5.0, 1e-6);
}
