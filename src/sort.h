#pragma once

#include <vector>
#include <cstdint>

namespace gslib {

// Sort primary array and rearrange secondary arrays in the same order.
// Returns sorted indices.
std::vector<int> sortem(
    std::vector<double>& primary,
    std::vector<std::vector<double>*> secondary = {}
);

} // namespace gslib
