#pragma once

#include <vector>

namespace gslib {

// Covariance between two points given a variogram model.
// Returns {cmax, cova}.
// Parameters match Fortran dcova3 (double precision version):
//   nst: number of nested structures
//   ivarg: variogram number (1-based in Fortran; here 0-based)
//   c0: nugget array (size >= ivarg+1)
//   it: structure types, cc: sill contributions, aa: ranges
//       (all size nst * (ivarg+1), indexed from ivarg*nst)
//   irot: rotation matrix index for first structure (0-based)
//   rotmat: rotation matrices (flat, maxrot*9)
struct CovResult {
    double cmax;
    double cova;
};

CovResult cova3(double x1, double y1, double z1,
                double x2, double y2, double z2,
                int ivarg, int nst,
                const std::vector<double>& c0,
                const std::vector<int>& it,
                const std::vector<double>& cc,
                const std::vector<double>& aa,
                int irot, int maxrot,
                const std::vector<double>& rotmat);

} // namespace gslib
