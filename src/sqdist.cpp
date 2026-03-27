#include "sqdist.h"
#include "common.h"

namespace gslib {

double sqdist(double x1, double y1, double z1,
              double x2, double y2, double z2,
              int ind, int maxrot, const std::vector<double>& rotmat) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    double result = 0.0;
    int base = ind * 9;
    for (int i = 0; i < 3; i++) {
        double cont = rotmat[base + i * 3 + 0] * dx
                    + rotmat[base + i * 3 + 1] * dy
                    + rotmat[base + i * 3 + 2] * dz;
        result += cont * cont;
    }
    return result;
}

} // namespace gslib
