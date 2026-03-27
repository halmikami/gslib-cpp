#include "gamv.h"
#include "common.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace gslib {

VariogramResult gamv(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& values,
    int n_lags,
    double xlag,
    double xltol,
    double azm,
    double atol,
    double bandwh,
    double dip,
    double dtol,
    double bandwd,
    const std::vector<int>& bhid,
    int ivtype)
{
    int nd = static_cast<int>(x.size());
    if (nd == 0) {
        return {};
    }

    // Default tolerance
    double xltoll = xltol;
    if (xltoll <= 0.0) xltoll = 0.5 * xlag;

    // Single direction, single variogram setup
    // Compute direction vectors
    double azmuth = (90.0 - azm) * PI / 180.0;
    double uvxazm = std::cos(azmuth);
    double uvyazm = std::sin(azmuth);
    double csatol = (atol <= 0.0) ? std::cos(45.0 * PI / 180.0)
                                  : std::cos(atol * PI / 180.0);

    double declin = (90.0 - dip) * PI / 180.0;
    double uvzdec = std::cos(declin);
    double uvhdec = std::sin(declin);
    double csdtol = (dtol <= 0.0) ? std::cos(45.0 * PI / 180.0)
                                  : std::cos(dtol * PI / 180.0);

    // Allocate output arrays: n_lags + 2 bins (bin 0 = zero distance, bins 1..n_lags+1)
    int nbins = n_lags + 2;
    std::vector<double> np_arr(nbins, 0.0);
    std::vector<double> dis(nbins, 0.0);
    std::vector<double> gam(nbins, 0.0);
    std::vector<double> hm(nbins, 0.0);
    std::vector<double> tm(nbins, 0.0);
    std::vector<double> hv(nbins, 0.0);
    std::vector<double> tv(nbins, 0.0);

    double dismxs = std::pow((static_cast<double>(n_lags) + 0.5 - 1.0e-20) * xlag, 2);

    bool use_bhid = !bhid.empty();
    bool omni = (atol >= 90.0);

    // Main loop over all pairs
    for (int i = 0; i < nd; i++) {
        for (int j = i + 1; j < nd; j++) {
            // Downhole constraint
            if (use_bhid && bhid[j] != bhid[i]) continue;

            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];
            double dxs = dx * dx;
            double dys = dy * dy;
            double dzs = dz * dz;
            double hs = dxs + dys + dzs;
            if (hs > dismxs) continue;
            if (hs < 0.0) hs = 0.0;
            double h = std::sqrt(hs);

            // Determine lag bin
            int lagbeg, lagend;
            if (h <= 1.0e-20) {
                lagbeg = 0;
                lagend = 0;
            } else {
                lagbeg = -1;
                lagend = -1;
                for (int ilag = 1; ilag < nbins; ilag++) {
                    double lag_dist = xlag * static_cast<double>(ilag - 1);
                    if (h >= (lag_dist - xltoll) && h <= (lag_dist + xltoll)) {
                        if (lagbeg < 0) lagbeg = ilag;
                        lagend = ilag;
                    }
                }
                if (lagend < 0) continue;
            }

            // Check azimuth
            double dxy = std::sqrt(std::max(dxs + dys, 0.0));
            double dcazm;
            if (dxy < 1.0e-20) {
                dcazm = 1.0;
            } else {
                dcazm = (dx * uvxazm + dy * uvyazm) / dxy;
            }
            if (std::abs(dcazm) < csatol) continue;

            // Check horizontal bandwidth
            double band = uvxazm * dy - uvyazm * dx;
            if (std::abs(band) > bandwh) continue;

            // Check dip
            if (dcazm < 0.0) dxy = -dxy;
            double dcdec;
            if (lagbeg == 0) {
                dcdec = 0.0;
            } else {
                dcdec = (dxy * uvhdec + dz * uvzdec) / h;
                if (std::abs(dcdec) < csdtol) continue;
            }

            // Check vertical bandwidth
            band = uvhdec * dz - uvzdec * dxy;
            if (std::abs(band) > bandwd) continue;

            // Determine head/tail values
            double vrh, vrt, vrhpr, vrtpr;
            if (dcazm >= 0.0 && dcdec >= 0.0) {
                vrh = values[i];
                vrt = values[j];
                if (omni || ivtype == 2) {
                    vrtpr = values[i];
                    vrhpr = values[j];
                }
            } else {
                vrh = values[j];
                vrt = values[i];
                if (omni || ivtype == 2) {
                    vrtpr = values[j];
                    vrhpr = values[i];
                }
            }

            // Accumulate by variogram type
            int it = ivtype;

            if (it == 1 || it == 5 || it >= 9) {
                // Semivariogram (and indicator, general relative)
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis[il] += h;
                    tm[il] += vrt;
                    hm[il] += vrh;
                    gam[il] += (vrh - vrt) * (vrh - vrt);
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrtpr;
                        hm[il] += vrhpr;
                        gam[il] += (vrhpr - vrtpr) * (vrhpr - vrtpr);
                    }
                }
            } else if (it == 2) {
                // Cross-semivariogram
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis[il] += h;
                    tm[il] += 0.5 * (vrt + vrtpr);
                    hm[il] += 0.5 * (vrh + vrhpr);
                    gam[il] += (vrhpr - vrh) * (vrt - vrtpr);
                }
            } else if (it == 3) {
                // Covariance
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis[il] += h;
                    tm[il] += vrt;
                    hm[il] += vrh;
                    gam[il] += vrh * vrt;
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrtpr;
                        hm[il] += vrhpr;
                        gam[il] += vrhpr * vrtpr;
                    }
                }
            } else if (it == 4) {
                // Correlogram
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis[il] += h;
                    tm[il] += vrt;
                    hm[il] += vrh;
                    hv[il] += vrh * vrh;
                    tv[il] += vrt * vrt;
                    gam[il] += vrh * vrt;
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrtpr;
                        hm[il] += vrhpr;
                        hv[il] += vrhpr * vrhpr;
                        tv[il] += vrtpr * vrtpr;
                        gam[il] += vrhpr * vrtpr;
                    }
                }
            } else if (it == 6) {
                // Pairwise relative
                for (int il = lagbeg; il <= lagend; il++) {
                    if (std::abs(vrt + vrh) > 1.0e-20) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrt;
                        hm[il] += vrh;
                        double gamma = 2.0 * (vrt - vrh) / (vrt + vrh);
                        gam[il] += gamma * gamma;
                    }
                }
            } else if (it == 7) {
                // Log variogram
                for (int il = lagbeg; il <= lagend; il++) {
                    if (vrt > 1.0e-20 && vrh > 1.0e-20) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrt;
                        hm[il] += vrh;
                        double gamma = std::log(vrt) - std::log(vrh);
                        gam[il] += gamma * gamma;
                    }
                }
            } else if (it == 8) {
                // Madogram
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis[il] += h;
                    tm[il] += vrt;
                    hm[il] += vrh;
                    gam[il] += std::abs(vrh - vrt);
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis[il] += h;
                        tm[il] += vrtpr;
                        hm[il] += vrhpr;
                        gam[il] += std::abs(vrhpr - vrtpr);
                    }
                }
            }
        }
    }

    // Compute averages and correct variogram measures
    for (int il = 0; il < nbins; il++) {
        if (np_arr[il] <= 0.0) continue;
        double rnum = np_arr[il];
        dis[il] /= rnum;
        gam[il] /= rnum;
        hm[il] /= rnum;
        tm[il] /= rnum;
        hv[il] /= rnum;
        tv[il] /= rnum;

        if (ivtype == 1 || ivtype == 2) {
            gam[il] = 0.5 * gam[il];
        } else if (ivtype == 3) {
            gam[il] = gam[il] - hm[il] * tm[il];
        } else if (ivtype == 4) {
            hv[il] = hv[il] - hm[il] * hm[il];
            if (hv[il] < 0.0) hv[il] = 0.0;
            hv[il] = std::sqrt(hv[il]);
            tv[il] = tv[il] - tm[il] * tm[il];
            if (tv[il] < 0.0) tv[il] = 0.0;
            tv[il] = std::sqrt(tv[il]);
            if (hv[il] * tv[il] < 1.0e-20) {
                gam[il] = 0.0;
            } else {
                gam[il] = (gam[il] - hm[il] * tm[il]) / (hv[il] * tv[il]);
            }
            hv[il] = hv[il] * hv[il];
            tv[il] = tv[il] * tv[il];
        } else if (ivtype == 5) {
            double htave = 0.5 * (hm[il] + tm[il]);
            htave = htave * htave;
            if (htave < 1.0e-20) {
                gam[il] = 0.0;
            } else {
                gam[il] = gam[il] / htave;
            }
        } else if (ivtype >= 6) {
            gam[il] = 0.5 * gam[il];
        }
    }

    // Build result: skip bin 0 (zero-distance), return bins 1..n_lags+1
    VariogramResult result;
    result.lags.resize(n_lags);
    result.semivariance.resize(n_lags);
    result.pair_counts.resize(n_lags);
    result.head_mean.resize(n_lags);
    result.tail_mean.resize(n_lags);
    result.head_var.resize(n_lags);
    result.tail_var.resize(n_lags);

    for (int il = 0; il < n_lags; il++) {
        int bin = il + 1;  // bins 1..n_lags
        result.pair_counts[il] = static_cast<int>(np_arr[bin]);
        if (np_arr[bin] > 0) {
            result.lags[il] = dis[bin];
        } else {
            result.lags[il] = xlag * static_cast<double>(il);
        }
        result.semivariance[il] = gam[bin];
        result.head_mean[il] = hm[bin];
        result.tail_mean[il] = tm[bin];
        result.head_var[il] = hv[bin];
        result.tail_var[il] = tv[bin];
    }

    return result;
}

} // namespace gslib
