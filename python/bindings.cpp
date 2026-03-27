#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "gamv.h"
#include "kt3d.h"
#include "cova3.h"
#include "setrot.h"
#include "sqdist.h"
#include "declus.h"
#include "sort.h"

namespace py = pybind11;

PYBIND11_MODULE(gslib_cpp, m) {
    m.doc() = "GSLIB (Geostatistical Software Library) - C++17 implementation";

    // VariogramResult
    py::class_<gslib::VariogramResult>(m, "VariogramResult")
        .def_readonly("lags", &gslib::VariogramResult::lags)
        .def_readonly("semivariance", &gslib::VariogramResult::semivariance)
        .def_readonly("pair_counts", &gslib::VariogramResult::pair_counts)
        .def_readonly("head_mean", &gslib::VariogramResult::head_mean)
        .def_readonly("tail_mean", &gslib::VariogramResult::tail_mean)
        .def_readonly("head_var", &gslib::VariogramResult::head_var)
        .def_readonly("tail_var", &gslib::VariogramResult::tail_var);

    // gamv
    m.def("gamv", &gslib::gamv,
        "Experimental variogram for irregularly spaced 3D data",
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("values"),
        py::arg("n_lags") = 15,
        py::arg("xlag") = 10.0,
        py::arg("xltol") = -1.0,
        py::arg("azm") = 0.0,
        py::arg("atol") = 90.0,
        py::arg("bandwh") = 1.0e10,
        py::arg("dip") = 0.0,
        py::arg("dtol") = 90.0,
        py::arg("bandwd") = 1.0e10,
        py::arg("bhid") = std::vector<int>{},
        py::arg("ivtype") = 1);

    // KrigingResult
    py::class_<gslib::KrigingResult>(m, "KrigingResult")
        .def_readonly("estimates", &gslib::KrigingResult::estimates)
        .def_readonly("variances", &gslib::KrigingResult::variances);

    // kt3d
    m.def("kt3d", &gslib::kt3d,
        "Ordinary/Simple Kriging estimation in 3D",
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("values"),
        py::arg("xout"), py::arg("yout"), py::arg("zout"),
        py::arg("search_radius") = 200.0,
        py::arg("search_radius2") = -1.0,
        py::arg("search_radius3") = -1.0,
        py::arg("sang1") = 0.0,
        py::arg("sang2") = 0.0,
        py::arg("sang3") = 0.0,
        py::arg("ndmax") = 16,
        py::arg("ndmin") = 4,
        py::arg("noct") = 0,
        py::arg("nugget") = 0.0,
        py::arg("model_types") = std::vector<int>{1},
        py::arg("model_cc") = std::vector<double>{1.0},
        py::arg("model_aa") = std::vector<double>{100.0},
        py::arg("model_aa1") = std::vector<double>{},
        py::arg("model_aa2") = std::vector<double>{},
        py::arg("model_ang1") = std::vector<double>{},
        py::arg("model_ang2") = std::vector<double>{},
        py::arg("model_ang3") = std::vector<double>{},
        py::arg("ktype") = 1,
        py::arg("skmean") = 0.0,
        py::arg("unest") = -999.0);

    // DeclusResult
    py::class_<gslib::DeclusResult>(m, "DeclusResult")
        .def_readonly("weights", &gslib::DeclusResult::weights)
        .def_readonly("declustered_mean", &gslib::DeclusResult::declustered_mean)
        .def_readonly("weight_min", &gslib::DeclusResult::weight_min)
        .def_readonly("weight_max", &gslib::DeclusResult::weight_max)
        .def_readonly("cell_sizes", &gslib::DeclusResult::cell_sizes)
        .def_readonly("cell_means", &gslib::DeclusResult::cell_means);

    // declus
    m.def("declus", &gslib::declus,
        "Cell declustering (Deutsch, 1989)",
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("vr"),
        py::arg("anisy") = 1.0,
        py::arg("anisz") = 1.0,
        py::arg("minmax") = 0,
        py::arg("ncell") = 24,
        py::arg("cmin") = 1.0,
        py::arg("cmax") = 100.0,
        py::arg("noff") = 8,
        py::arg("maxcel") = 0);

    // CovResult
    py::class_<gslib::CovResult>(m, "CovResult")
        .def_readonly("cmax", &gslib::CovResult::cmax)
        .def_readonly("cova", &gslib::CovResult::cova);
}
