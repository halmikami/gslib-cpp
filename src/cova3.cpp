#include "cova3.h"
#include "sqdist.h"
#include "common.h"
#include <cmath>

namespace gslib {

CovResult cova3(double x1, double y1, double z1,
                double x2, double y2, double z2,
                int ivarg, int nst,
                const std::vector<double>& c0,
                const std::vector<int>& it,
                const std::vector<double>& cc,
                const std::vector<double>& aa,
                int irot, int maxrot,
                const std::vector<double>& rotmat) {
    // Fortran uses 1-based ivarg; C++ uses 0-based
    int istart = ivarg * nst;

    // Maximum covariance
    double cmax = c0[ivarg];
    for (int is = 0; is < nst; is++) {
        int ist = istart + is;
        if (it[ist] == 4) {
            cmax += PMX;
        } else {
            cmax += cc[ist];
        }
    }

    // Check for zero distance
    double hsqd = sqdist(x1, y1, z1, x2, y2, z2, irot, maxrot, rotmat);
    if (hsqd < EPSLON) {
        return {cmax, cmax};
    }

    // Loop over structures
    double cova = 0.0;
    for (int is = 0; is < nst; is++) {
        int ist = istart + is;

        // Compute distance with appropriate rotation matrix
        if (ist != 0) {
            int ir = std::min(irot + is, maxrot - 1);
            hsqd = sqdist(x1, y1, z1, x2, y2, z2, ir, maxrot, rotmat);
        }
        double h = std::sqrt(hsqd);

        switch (it[ist]) {
            case 1: { // Spherical
                double hr = h / aa[ist];
                if (hr < 1.0) {
                    cova += cc[ist] * (1.0 - hr * (1.5 - 0.5 * hr * hr));
                }
                break;
            }
            case 2: { // Exponential
                cova += cc[ist] * std::exp(-3.0 * h / aa[ist]);
                break;
            }
            case 3: { // Gaussian
                double ratio = 3.0 * h / aa[ist];
                cova += cc[ist] * std::exp(-ratio * ratio);
                break;
            }
            case 4: { // Power
                cova += cmax - cc[ist] * std::pow(h, aa[ist]);
                break;
            }
            case 5: { // Hole effect
                cova += cc[ist] * std::cos(h / aa[ist] * PI);
                break;
            }
        }
    }

    return {cmax, cova};
}

} // namespace gslib
