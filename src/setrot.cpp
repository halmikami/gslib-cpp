#include "setrot.h"
#include "common.h"

namespace gslib {

void setrot(double ang1, double ang2, double ang3,
            double anis1, double anis2,
            int ind, int maxrot, std::vector<double>& rotmat) {
    // Convert GSLIB angles to mathematical angles:
    //   alpha: angle between major axis and E-W (counter-clockwise positive)
    //   beta:  dip of ellipsoid (positive down)
    //   theta: rotation of minor axis about major axis
    double alpha;
    if (ang1 >= 0.0 && ang1 < 270.0) {
        alpha = (90.0 - ang1) * DEG2RAD;
    } else {
        alpha = (450.0 - ang1) * DEG2RAD;
    }
    double beta = -1.0 * ang2 * DEG2RAD;
    double theta = ang3 * DEG2RAD;

    double sina = std::sin(alpha);
    double sinb = std::sin(beta);
    double sint = std::sin(theta);
    double cosa = std::cos(alpha);
    double cosb = std::cos(beta);
    double cost = std::cos(theta);

    double afac1 = 1.0 / std::max(anis1, 1.0e-20);
    double afac2 = 1.0 / std::max(anis2, 1.0e-20);

    int base = ind * 9;
    rotmat[base + 0] = cosb * cosa;
    rotmat[base + 1] = cosb * sina;
    rotmat[base + 2] = -sinb;
    rotmat[base + 3] = afac1 * (-cost * sina + sint * sinb * cosa);
    rotmat[base + 4] = afac1 * (cost * cosa + sint * sinb * sina);
    rotmat[base + 5] = afac1 * (sint * cosb);
    rotmat[base + 6] = afac2 * (sint * sina + cost * sinb * cosa);
    rotmat[base + 7] = afac2 * (-sint * cosa + cost * sinb * sina);
    rotmat[base + 8] = afac2 * (cost * cosb);
}

} // namespace gslib
