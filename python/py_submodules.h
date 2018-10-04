#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


namespace tools {
  void bind_tools(py::module& m);
}

namespace fields {
  void bind_fields(py::module& m);
}

namespace vlv {
  void bind_vlv(py::module& m);
}

namespace pic {
  void bind_pic(py::module& m);
}

namespace rad {
  void bind_rad(py::module& m);
}



