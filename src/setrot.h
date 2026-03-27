#pragma once

#include <vector>

namespace gslib {

// Rotation matrix: indexed as rotmat[ind][i][j] where ind is structure index,
// i is row (0-2), j is col (0-2). Stored as flat vector of size maxrot*3*3.
// Access: rotmat[ind*9 + i*3 + j]

void setrot(double ang1, double ang2, double ang3,
            double anis1, double anis2,
            int ind, int maxrot, std::vector<double>& rotmat);

} // namespace gslib
