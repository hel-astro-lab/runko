#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


namespace tools {
  void bind_tools(py::module& m);
}

namespace emf {
  void bind_emf(py::module& m);
}

namespace vlv {
  void bind_vlv(py::module& m);
}

namespace pic {
  void bind_pic(py::module& m);
}

namespace qed {
  void bind_qed(py::module& m);
}

namespace ffe {
  void bind_ffe(py::module& m);
}

namespace cpl {
  void bind_cpl(py::module& m);
}
