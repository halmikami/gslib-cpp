#pragma once

#include <vector>

namespace gslib {

double sqdist(double x1, double y1, double z1,
              double x2, double y2, double z2,
              int ind, int maxrot, const std::vector<double>& rotmat);

} // namespace gslib
