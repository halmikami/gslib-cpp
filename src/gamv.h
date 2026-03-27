#pragma once

#include <vector>
#include <cstdint>

namespace gslib {

struct VariogramResult {
    std::vector<double> lags;
    std::vector<double> semivariance;
    std::vector<int> pair_counts;
    std::vector<double> head_mean;
    std::vector<double> tail_mean;
    std::vector<double> head_var;
    std::vector<double> tail_var;
};

// Experimental variogram for irregularly spaced 3D data.
// Supports semivariogram, cross-semivariogram, covariance, correlogram,
// general relative, pairwise relative, log variogram, madogram,
// and indicator variograms.
//
// ivtype: 1=semivariogram, 2=cross-semi, 3=covariance, 4=correlogram,
//         5=general relative, 6=pairwise relative, 7=log variogram,
//         8=madogram, 9=indicator continuous, 10=indicator categorical
VariogramResult gamv(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& values,
    int n_lags = 15,
    double xlag = 10.0,
    double xltol = -1.0,       // <0 means auto = xlag/2
    double azm = 0.0,
    double atol = 90.0,
    double bandwh = 1.0e10,
    double dip = 0.0,
    double dtol = 90.0,
    double bandwd = 1.0e10,
    const std::vector<int>& bhid = {},
    int ivtype = 1
);

} // namespace gslib
