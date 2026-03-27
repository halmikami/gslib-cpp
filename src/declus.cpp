#include "declus.h"
#include "common.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace gslib {

DeclusResult declus(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& vr,
    double anisy,
    double anisz,
    int minmax,
    int ncell,
    double cmin,
    double cmax,
    int noff,
    int maxcel)
{
    int nd = static_cast<int>(x.size());
    if (nd == 0) {
        throw std::invalid_argument("declus: no data points");
    }

    bool minimize = (minmax == 0);
    double roff = static_cast<double>(noff);

    // Compute data extents and average
    double xmin = 1.0e21, ymin = 1.0e21, zmin = 1.0e21;
    double xmax = -1.0e21, ymax = -1.0e21, zmax = -1.0e21;
    double vrav = 0.0;

    std::vector<double> wtopt(nd, 1.0);
    std::vector<double> wt(nd, 0.0);
    std::vector<int> index(nd);

    for (int i = 0; i < nd; i++) {
        vrav += vr[i];
        xmin = std::min(xmin, x[i]);
        xmax = std::max(xmax, x[i]);
        ymin = std::min(ymin, y[i]);
        ymax = std::max(ymax, y[i]);
        zmin = std::min(zmin, z[i]);
        zmax = std::max(zmax, z[i]);
    }
    vrav /= static_cast<double>(nd);

    double vrop = vrav;
    double best = 0.0;

    double xo1 = xmin - 0.01;
    double yo1 = ymin - 0.01;
    double zo1 = zmin - 0.01;

    double tcmax = (ncell == 1) ? cmin : cmax;

    double xinc = (tcmax - cmin) / static_cast<double>(ncell);
    double yinc = anisy * xinc;
    double zinc = anisz * xinc;

    double xcs = cmin - xinc;
    double ycs = cmin * anisy - yinc;
    double zcs = cmin * anisz - zinc;

    DeclusResult result;
    result.cell_sizes.resize(ncell + 1);
    result.cell_means.resize(ncell + 1);

    for (int lp = 0; lp <= ncell; lp++) {
        xcs += xinc;
        ycs += yinc;
        zcs += zinc;

        std::fill(wt.begin(), wt.end(), 0.0);

        int ncellx = static_cast<int>((xmax - (xo1 - xcs)) / xcs) + 1;
        int ncelly = static_cast<int>((ymax - (yo1 - ycs)) / ycs) + 1;
        int ncellz = static_cast<int>((zmax - (zo1 - zcs)) / zcs) + 1;
        int ncellt = ncellx * ncelly * ncellz;

        if (maxcel > 0 && ncellt > maxcel) {
            throw std::runtime_error("declus: ncellt > maxcel, increase cmin or maxcel");
        }

        double xfac = std::min(xcs / roff, 0.5 * (xmax - xmin));
        double yfac = std::min(ycs / roff, 0.5 * (ymax - ymin));
        double zfac = std::min(zcs / roff, 0.5 * (zmax - zmin));

        for (int kp = 0; kp < noff; kp++) {
            double xo = xo1 - static_cast<double>(kp) * xfac;
            double yo = yo1 - static_cast<double>(kp) * yfac;
            double zo = zo1 - static_cast<double>(kp) * zfac;

            std::vector<double> cellwt(ncellt, 0.0);

            for (int i = 0; i < nd; i++) {
                int icellx = static_cast<int>((x[i] - xo) / xcs);
                int icelly = static_cast<int>((y[i] - yo) / ycs);
                int icellz = static_cast<int>((z[i] - zo) / zcs);
                int icell = icellx + icelly * ncellx + icellz * ncelly * ncellx;
                icell = std::max(0, std::min(icell, ncellt - 1));
                index[i] = icell;
                cellwt[icell] += 1.0;
            }

            double sumw = 0.0;
            for (int i = 0; i < nd; i++) {
                sumw += 1.0 / cellwt[index[i]];
            }
            sumw = 1.0 / sumw;

            for (int i = 0; i < nd; i++) {
                wt[i] += (1.0 / cellwt[index[i]]) * sumw;
            }
        }

        // Compute weighted average
        double sumw = 0.0, sumwg = 0.0;
        for (int i = 0; i < nd; i++) {
            sumw += wt[i];
            sumwg += wt[i] * vr[i];
        }
        double vrcr = sumwg / sumw;

        result.cell_sizes[lp] = xcs;
        result.cell_means[lp] = vrcr;

        if ((minimize && vrcr < vrop) || (!minimize && vrcr > vrop) || ncell == 1) {
            best = xcs;
            vrop = vrcr;
            wtopt = wt;
        }
    }

    // Normalize weights to sum to nd
    double sumw = 0.0;
    for (int i = 0; i < nd; i++) sumw += wtopt[i];
    double facto = static_cast<double>(nd) / sumw;
    double wtmin = 1.0e30, wtmax = -1.0e30;
    for (int i = 0; i < nd; i++) {
        wtopt[i] *= facto;
        wtmin = std::min(wtmin, wtopt[i]);
        wtmax = std::max(wtmax, wtopt[i]);
    }

    result.weights = std::move(wtopt);
    result.declustered_mean = vrop;
    result.weight_min = wtmin;
    result.weight_max = wtmax;

    return result;
}

} // namespace gslib
