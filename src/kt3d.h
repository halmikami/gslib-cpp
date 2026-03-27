#pragma once

#include <vector>
#include <cstdint>

namespace gslib {

struct KrigingResult {
    std::vector<double> estimates;
    std::vector<double> variances;
};

// Ordinary/Simple Kriging estimation in 3D.
// ktype: 0=Simple Kriging, 1=Ordinary Kriging
KrigingResult kt3d(
    // Input data
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& values,
    // Output locations
    const std::vector<double>& xout,
    const std::vector<double>& yout,
    const std::vector<double>& zout,
    // Search parameters
    double search_radius = 200.0,
    double search_radius2 = -1.0,   // <0 means same as search_radius
    double search_radius3 = -1.0,   // <0 means same as search_radius
    double sang1 = 0.0,
    double sang2 = 0.0,
    double sang3 = 0.0,
    int ndmax = 16,
    int ndmin = 4,
    int noct = 0,
    // Variogram model
    double nugget = 0.0,
    const std::vector<int>& model_types = {1},
    const std::vector<double>& model_cc = {1.0},
    const std::vector<double>& model_aa = {100.0},
    const std::vector<double>& model_aa1 = {},  // empty = same as aa
    const std::vector<double>& model_aa2 = {},  // empty = same as aa
    const std::vector<double>& model_ang1 = {},
    const std::vector<double>& model_ang2 = {},
    const std::vector<double>& model_ang3 = {},
    // Kriging type
    int ktype = 1,
    double skmean = 0.0,
    // Unestimated value
    double unest = -999.0
);

} // namespace gslib
