# gslib-cpp

Modern C++17 implementation of GSLIB (Geostatistical Software Library) with Python bindings via pybind11.

Ported from the Fortran 77/90 code in [PyGSLIB](https://github.com/opengeostat/pygslib).

## Features

- **gamv** — Experimental variogram for irregularly spaced 3D data (10 variogram types)
- **kt3d** — Ordinary/Simple Kriging estimation
- **cova3** — Covariance models (spherical, exponential, Gaussian, power, hole effect)
- **setrot** — Anisotropic rotation matrices
- **sqdist** — Squared anisotropic distance
- **declus** — Cell declustering (Deutsch, 1989)
- **sortem** — Sort utility

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
ctest --output-on-failure
```

### With Python bindings

```bash
cmake .. -DGSLIB_BUILD_PYTHON=ON
make -j$(nproc)
```

## Quick Start (C++)

```cpp
#include "gamv.h"
#include "kt3d.h"

// Experimental variogram
auto vario = gslib::gamv(x, y, z, values, /*n_lags=*/15, /*xlag=*/10.0);

// Ordinary Kriging
auto krig = gslib::kt3d(x, y, z, values, xout, yout, zout,
                         /*search_radius=*/200.0);
```

## Quick Start (Python)

```python
import gslib_cpp

result = gslib_cpp.gamv(x, y, z, values, n_lags=15, xlag=20.0)
krig = gslib_cpp.kt3d(x, y, z, values, xout, yout, zout,
                       search_radius=200.0,
                       model_types=[1], model_cc=[1.0], model_aa=[100.0])
```

## References

- Deutsch, C.V. & Journel, A.G. (1997). *GSLIB: Geostatistical Software Library and User's Guide*. Oxford University Press.
- Journel, A.G. & Huijbregts, C.J. (1978). *Mining Geostatistics*. Academic Press.

## License

MIT
