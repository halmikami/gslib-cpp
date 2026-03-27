#pragma once

#include <vector>

namespace gslib {

struct DeclusResult {
    std::vector<double> weights;       // Optimal declustering weights
    double declustered_mean;           // Optimal declustered mean
    double weight_min;
    double weight_max;
    std::vector<double> cell_sizes;    // Cell sizes tested (x)
    std::vector<double> cell_means;    // Declustered mean for each cell size
};

// Cell declustering (Deutsch, 1989).
// minmax: 0=minimize mean, 1=maximize mean
// ncell: number of cell sizes to test
// cmin, cmax: min/max cell size
// noff: number of origin offsets
// anisy, anisz: Y and Z cell anisotropy ratios
// maxcel: maximum number of cells (0 = no limit)
DeclusResult declus(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& vr,
    double anisy = 1.0,
    double anisz = 1.0,
    int minmax = 0,
    int ncell = 24,
    double cmin = 1.0,
    double cmax = 100.0,
    int noff = 8,
    int maxcel = 0
);

} // namespace gslib
