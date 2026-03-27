#include "gamv.h"
#include "common.h"
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

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

    int nbins = n_lags + 2;
    double dismxs = std::pow((static_cast<double>(n_lags) + 0.5 - 1.0e-20) * xlag, 2);

    bool use_bhid = !bhid.empty();
    bool omni = (atol >= 90.0);

    // Determine number of threads for OpenMP reduction
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif

    // Per-thread accumulators
    std::vector<std::vector<double>> t_np(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_dis(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_gam(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_hm(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_tm(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_hv(nthreads, std::vector<double>(nbins, 0.0));
    std::vector<std::vector<double>> t_tv(nthreads, std::vector<double>(nbins, 0.0));

    // Main loop over all pairs — OpenMP parallelized over outer loop
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < nd; i++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        auto& np_arr = t_np[tid];
        auto& dis_l  = t_dis[tid];
        auto& gam_l  = t_gam[tid];
        auto& hm_l   = t_hm[tid];
        auto& tm_l   = t_tm[tid];
        auto& hv_l   = t_hv[tid];
        auto& tv_l   = t_tv[tid];

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
            double vrh, vrt, vrhpr = 0.0, vrtpr = 0.0;
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
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis_l[il] += h;
                    tm_l[il] += vrt;
                    hm_l[il] += vrh;
                    gam_l[il] += (vrh - vrt) * (vrh - vrt);
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrtpr;
                        hm_l[il] += vrhpr;
                        gam_l[il] += (vrhpr - vrtpr) * (vrhpr - vrtpr);
                    }
                }
            } else if (it == 2) {
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis_l[il] += h;
                    tm_l[il] += 0.5 * (vrt + vrtpr);
                    hm_l[il] += 0.5 * (vrh + vrhpr);
                    gam_l[il] += (vrhpr - vrh) * (vrt - vrtpr);
                }
            } else if (it == 3) {
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis_l[il] += h;
                    tm_l[il] += vrt;
                    hm_l[il] += vrh;
                    gam_l[il] += vrh * vrt;
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrtpr;
                        hm_l[il] += vrhpr;
                        gam_l[il] += vrhpr * vrtpr;
                    }
                }
            } else if (it == 4) {
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis_l[il] += h;
                    tm_l[il] += vrt;
                    hm_l[il] += vrh;
                    hv_l[il] += vrh * vrh;
                    tv_l[il] += vrt * vrt;
                    gam_l[il] += vrh * vrt;
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrtpr;
                        hm_l[il] += vrhpr;
                        hv_l[il] += vrhpr * vrhpr;
                        tv_l[il] += vrtpr * vrtpr;
                        gam_l[il] += vrhpr * vrtpr;
                    }
                }
            } else if (it == 6) {
                for (int il = lagbeg; il <= lagend; il++) {
                    if (std::abs(vrt + vrh) > 1.0e-20) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrt;
                        hm_l[il] += vrh;
                        double gamma = 2.0 * (vrt - vrh) / (vrt + vrh);
                        gam_l[il] += gamma * gamma;
                    }
                }
            } else if (it == 7) {
                for (int il = lagbeg; il <= lagend; il++) {
                    if (vrt > 1.0e-20 && vrh > 1.0e-20) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrt;
                        hm_l[il] += vrh;
                        double gamma = std::log(vrt) - std::log(vrh);
                        gam_l[il] += gamma * gamma;
                    }
                }
            } else if (it == 8) {
                for (int il = lagbeg; il <= lagend; il++) {
                    np_arr[il] += 1.0;
                    dis_l[il] += h;
                    tm_l[il] += vrt;
                    hm_l[il] += vrh;
                    gam_l[il] += std::abs(vrh - vrt);
                    if (omni) {
                        np_arr[il] += 1.0;
                        dis_l[il] += h;
                        tm_l[il] += vrtpr;
                        hm_l[il] += vrhpr;
                        gam_l[il] += std::abs(vrhpr - vrtpr);
                    }
                }
            }
        }
    }

    // Merge thread-local accumulators
    std::vector<double> np_arr(nbins, 0.0), dis(nbins, 0.0), gam(nbins, 0.0);
    std::vector<double> hm(nbins, 0.0), tm(nbins, 0.0), hv(nbins, 0.0), tv(nbins, 0.0);
    for (int t = 0; t < nthreads; t++) {
        for (int il = 0; il < nbins; il++) {
            np_arr[il] += t_np[t][il];
            dis[il]    += t_dis[t][il];
            gam[il]    += t_gam[t][il];
            hm[il]     += t_hm[t][il];
            tm[il]     += t_tm[t][il];
            hv[il]     += t_hv[t][il];
            tv[il]     += t_tv[t][il];
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

    // Build result
    VariogramResult result;
    result.lags.resize(n_lags);
    result.semivariance.resize(n_lags);
    result.pair_counts.resize(n_lags);
    result.head_mean.resize(n_lags);
    result.tail_mean.resize(n_lags);
    result.head_var.resize(n_lags);
    result.tail_var.resize(n_lags);

    for (int il = 0; il < n_lags; il++) {
        int bin = il + 1;
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
