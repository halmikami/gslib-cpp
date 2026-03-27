#include "sort.h"
#include <algorithm>
#include <numeric>

namespace gslib {

std::vector<int> sortem(
    std::vector<double>& primary,
    std::vector<std::vector<double>*> secondary)
{
    int n = static_cast<int>(primary.size());
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&primary](int a, int b) { return primary[a] < primary[b]; });

    // Apply permutation to primary
    std::vector<double> tmp(n);
    for (int i = 0; i < n; i++) tmp[i] = primary[indices[i]];
    primary = tmp;

    // Apply same permutation to secondary arrays
    for (auto* sec : secondary) {
        if (sec && static_cast<int>(sec->size()) == n) {
            for (int i = 0; i < n; i++) tmp[i] = (*sec)[indices[i]];
            *sec = tmp;
        }
    }

    return indices;
}

} // namespace gslib
