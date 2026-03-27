#include "kt3d.h"
#include "common.h"
#include "setrot.h"
#include "sqdist.h"
#include "cova3.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

namespace gslib {

// Simple Gauss elimination solver for symmetric positive-definite kriging system
static bool ktsol(int neq, std::vector<double>& a, std::vector<double>& r) {
    // Forward elimination
    for (int k = 0; k < neq - 1; k++) {
        double akk = a[k * neq + k];
        if (std::abs(akk) < 1.0e-20) return false;
        for (int i = k + 1; i < neq; i++) {
            double factor = a[i * neq + k] / akk;
            for (int j = k + 1; j < neq; j++) {
                a[i * neq + j] -= factor * a[k * neq + j];
            }
            r[i] -= factor * r[k];
        }
    }
    // Back substitution
    for (int i = neq - 1; i >= 0; i--) {
        double aii = a[i * neq + i];
        if (std::abs(aii) < 1.0e-20) return false;
        for (int j = i + 1; j < neq; j++) {
            r[i] -= a[i * neq + j] * r[j];
        }
        r[i] /= aii;
    }
    return true;
}

KrigingResult kt3d(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& values,
    const std::vector<double>& xout,
    const std::vector<double>& yout,
    const std::vector<double>& zout,
    double search_radius,
    double search_radius2,
    double search_radius3,
    double sang1,
    double sang2,
    double sang3,
    int ndmax,
    int ndmin,
    int noct,
    double nugget,
    const std::vector<int>& model_types,
    const std::vector<double>& model_cc,
    const std::vector<double>& model_aa,
    const std::vector<double>& model_aa1,
    const std::vector<double>& model_aa2,
    const std::vector<double>& model_ang1,
    const std::vector<double>& model_ang2,
    const std::vector<double>& model_ang3,
    int ktype,
    double skmean,
    double unest)
{
    int nd = static_cast<int>(x.size());
    int nout = static_cast<int>(xout.size());
    int nst = static_cast<int>(model_types.size());

    if (search_radius2 < 0.0) search_radius2 = search_radius;
    if (search_radius3 < 0.0) search_radius3 = search_radius;

    // Compute anisotropy ratios for search
    double sanis1 = search_radius2 / std::max(search_radius, EPSLON);
    double sanis2 = search_radius3 / std::max(search_radius, EPSLON);

    // Build variogram parameter arrays (single variogram, ivarg=0)
    std::vector<double> c0_arr = {nugget};
    std::vector<int> it_arr = model_types;
    std::vector<double> cc_arr = model_cc;
    std::vector<double> aa_arr = model_aa;

    // Build anisotropy ratios for variogram structures
    std::vector<double> anis1_arr(nst), anis2_arr(nst);
    std::vector<double> vang1(nst, 0.0), vang2(nst, 0.0), vang3(nst, 0.0);

    for (int i = 0; i < nst; i++) {
        double a1 = (i < static_cast<int>(model_aa1.size())) ? model_aa1[i] : model_aa[i];
        double a2 = (i < static_cast<int>(model_aa2.size())) ? model_aa2[i] : model_aa[i];
        anis1_arr[i] = a1 / std::max(model_aa[i], EPSLON);
        anis2_arr[i] = a2 / std::max(model_aa[i], EPSLON);
        if (i < static_cast<int>(model_ang1.size())) vang1[i] = model_ang1[i];
        if (i < static_cast<int>(model_ang2.size())) vang2[i] = model_ang2[i];
        if (i < static_cast<int>(model_ang3.size())) vang3[i] = model_ang3[i];
    }

    // Set up rotation matrices: nst structures + 1 for search
    int maxrot = nst + 1;
    std::vector<double> rotmat(maxrot * 9, 0.0);

    double covmax = nugget;
    for (int is = 0; is < nst; is++) {
        setrot(vang1[is], vang2[is], vang3[is],
               anis1_arr[is], anis2_arr[is], is, maxrot, rotmat);
        if (it_arr[is] == 4) {
            covmax += PMX;
        } else {
            covmax += cc_arr[is];
        }
    }
    // Search rotation matrix
    int isrot = nst;
    setrot(sang1, sang2, sang3, sanis1, sanis2, isrot, maxrot, rotmat);

    double radsqd = search_radius * search_radius;
    double unbias = covmax;  // For unbiasedness constraint scaling

    // Point kriging: cbb = point covariance at zero distance
    auto cov_result = cova3(0, 0, 0, 0, 0, 0,
                            0, nst, c0_arr, it_arr, cc_arr, aa_arr,
                            0, maxrot, rotmat);
    double cbb = cov_result.cova;

    KrigingResult result;
    result.estimates.resize(nout, unest);
    result.variances.resize(nout, unest);

    // Main loop over output locations
    for (int idx = 0; idx < nout; idx++) {
        double xloc = xout[idx];
        double yloc = yout[idx];
        double zloc = zout[idx];

        // Find nearby samples using search ellipsoid
        struct NearSample {
            int index;
            double dist;
        };
        std::vector<NearSample> nearby;
        nearby.reserve(nd);

        for (int i = 0; i < nd; i++) {
            double hsqd = sqdist(x[i], y[i], z[i], xloc, yloc, zloc,
                                 isrot, maxrot, rotmat);
            if (hsqd <= radsqd) {
                nearby.push_back({i, hsqd});
            }
        }

        // Sort by distance
        std::sort(nearby.begin(), nearby.end(),
                  [](const NearSample& a, const NearSample& b) {
                      return a.dist < b.dist;
                  });

        // Octant search if requested
        if (noct > 0 && static_cast<int>(nearby.size()) > ndmax) {
            std::vector<int> oct_count(8, 0);
            std::vector<NearSample> filtered;
            filtered.reserve(ndmax);
            for (const auto& ns : nearby) {
                double dx = x[ns.index] - xloc;
                double dy = y[ns.index] - yloc;
                double dz = z[ns.index] - zloc;
                int ioct = (dx >= 0 ? 0 : 1) + (dy >= 0 ? 0 : 2) + (dz >= 0 ? 0 : 4);
                if (oct_count[ioct] < noct) {
                    oct_count[ioct]++;
                    filtered.push_back(ns);
                    if (static_cast<int>(filtered.size()) >= ndmax) break;
                }
            }
            nearby = std::move(filtered);
        }

        // Limit to ndmax
        int na = std::min(static_cast<int>(nearby.size()), ndmax);
        if (na < ndmin) continue;

        // Build local coordinate arrays (relative to estimation point)
        std::vector<double> xa(na), ya(na), za(na), vra(na);
        for (int i = 0; i < na; i++) {
            int ind = nearby[i].index;
            xa[i] = x[ind];
            ya[i] = y[ind];
            za[i] = z[ind];
            vra[i] = values[ind];
        }

        // Handle single sample case
        if (na == 1) {
            auto cr1 = cova3(xa[0], ya[0], za[0], xa[0], ya[0], za[0],
                             0, nst, c0_arr, it_arr, cc_arr, aa_arr,
                             0, maxrot, rotmat);
            double cb1 = cr1.cova;
            auto cr2 = cova3(xa[0], ya[0], za[0], xloc, yloc, zloc,
                             0, nst, c0_arr, it_arr, cc_arr, aa_arr,
                             0, maxrot, rotmat);
            double cb = cr2.cova;
            result.estimates[idx] = vra[0];
            result.variances[idx] = cbb - 2.0 * cb + cb1;
            continue;
        }

        // Set up kriging system
        int neq = (ktype == 0) ? na : na + 1;  // OK adds unbiasedness constraint

        std::vector<double> a_mat(neq * neq, 0.0);
        std::vector<double> r_vec(neq, 0.0);
        std::vector<double> rr(neq, 0.0);

        // Fill covariance matrix
        for (int i = 0; i < na; i++) {
            for (int j = i; j < na; j++) {
                auto cr = cova3(xa[i], ya[i], za[i], xa[j], ya[j], za[j],
                                0, nst, c0_arr, it_arr, cc_arr, aa_arr,
                                0, maxrot, rotmat);
                a_mat[i * neq + j] = cr.cova;
                a_mat[j * neq + i] = cr.cova;
            }
        }

        // Unbiasedness constraints for OK
        if (ktype == 1) {
            for (int i = 0; i < na; i++) {
                a_mat[i * neq + na] = unbias;
                a_mat[na * neq + i] = unbias;
            }
        }

        // Right-hand side: covariance between data and estimation point
        for (int i = 0; i < na; i++) {
            auto cr = cova3(xa[i], ya[i], za[i], xloc, yloc, zloc,
                            0, nst, c0_arr, it_arr, cc_arr, aa_arr,
                            0, maxrot, rotmat);
            r_vec[i] = cr.cova;
        }
        if (ktype == 1) {
            r_vec[na] = unbias;
        }

        // Save RHS for variance computation
        rr = r_vec;

        // Solve kriging system
        if (!ktsol(neq, a_mat, r_vec)) {
            // Singular matrix
            continue;
        }

        // Compute estimate and variance
        double est = 0.0;
        double estv = cbb;
        for (int j = 0; j < neq; j++) {
            estv -= r_vec[j] * rr[j];
            if (j < na) {
                if (ktype == 0) {
                    est += r_vec[j] * (vra[j] - skmean);
                } else {
                    est += r_vec[j] * vra[j];
                }
            }
        }
        if (ktype == 0) est += skmean;

        result.estimates[idx] = est;
        result.variances[idx] = estv;
    }

    return result;
}

} // namespace gslib
