"""Basic tests for gslib_cpp Python bindings."""
import numpy as np


def test_gamv_constant():
    import gslib_cpp

    n = 50
    x = list(range(n))
    y = [0.0] * n
    z = [0.0] * n
    values = [5.0] * n

    result = gslib_cpp.gamv(x, y, z, values, n_lags=10, xlag=2.0)
    for i, sv in enumerate(result.semivariance):
        if result.pair_counts[i] > 0:
            assert abs(sv) < 1e-10, f"Lag {i}: semivariance should be 0 for constant"


def test_kt3d_constant():
    import gslib_cpp

    x = [0, 10, 20, 0, 10, 20, 0, 10, 20]
    y = [0, 0, 0, 10, 10, 10, 20, 20, 20]
    z = [0.0] * 9
    values = [10.0] * 9

    result = gslib_cpp.kt3d(
        x, y, z, values,
        xout=[10.0], yout=[10.0], zout=[0.0],
        search_radius=200.0,
        ndmin=1,
        nugget=0.0,
        model_types=[1],
        model_cc=[1.0],
        model_aa=[100.0],
    )
    assert abs(result.estimates[0] - 10.0) < 0.01


def test_declus_uniform():
    import gslib_cpp

    n = 25
    x = [float(i % 5) * 10 for i in range(n)]
    y = [float(i // 5) * 10 for i in range(n)]
    z = [0.0] * n
    vr = [5.0] * n

    result = gslib_cpp.declus(x, y, z, vr, ncell=5, cmin=5.0, cmax=30.0, noff=3)
    total = sum(result.weights)
    assert abs(total - n) < 1e-4, f"Weights should sum to {n}, got {total}"


if __name__ == "__main__":
    test_gamv_constant()
    print("PASS: test_gamv_constant")
    test_kt3d_constant()
    print("PASS: test_kt3d_constant")
    test_declus_uniform()
    print("PASS: test_declus_uniform")
    print("All tests passed!")
