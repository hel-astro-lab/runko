#include "corgi/corgi.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "runko/actions/env.h"
#include "tyvi/actions_ast.h"
#include "tyvi/actions_eval.h"

#include <chrono>
#include <exception>
#include <pika/execution.hpp>
#include <pika/init.hpp>
#include <pika/thread.hpp>
#include <print>
#include <ranges>
#include <thread>
#include <vector>

namespace {

namespace py = pybind11;
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace rn = std::ranges;
namespace rv = std::views;

ta::sexpr
  parse_element(const py::handle &obj)
{
  if(py::isinstance<py::int_>(obj)) {
    return obj.cast<long>();
  } else if(py::isinstance<py::str>(obj)) {
    return obj.cast<std::string>();
  } else if(py::isinstance<ta::intrinsic>(obj)) {
    return obj.cast<ta::intrinsic>();
  } else if(py::isinstance<runko::symbol>(obj)) {
    return obj.cast<runko::symbol>();
  } else if(py::isinstance<py::tuple>(obj)) {
    const auto tup = obj.cast<py::tuple>();

    if(rn::empty(tup)) { return ta::cons(); }

    const auto n = rn::size(tup);
    auto tail    = ta::sexpr { ta::cons(parse_element(tup[n - 1uz]), ta::null) };
    for(const auto i: rv::iota(0uz, n) | rv::reverse | rv::drop(1)) {
      tail = ta::cons(parse_element(tup[i]), std::move(tail));
    }

    return tail;
  } else {
    return "Trying to parse unsupported type.";
  }
}

void
  empty_context_eval(const py::handle &body_py)
{
  const auto body = parse_element(body_py);

  pika::start(0, nullptr);

  try {
    tyvi::this_thread::sync_wait(ta::eval<runko::symbol>(body, runko::build_stdenv()));
  } catch(const std::exception &e) {
    std::println("Exception while evaluation in print_context_eval: {}", e.what());
  }

  pika::finalize();
  pika::stop();
}
}  // namespace


namespace actions {

void
  bind_actions(py::module &m_sub)
{
  py::enum_<ta::intrinsic>(m_sub, "intrinsic")
    .value("car", ta::intrinsic::car)
    .value("cdr", ta::intrinsic::cdr)
    .value("quote", ta::intrinsic::quote)
    .export_values();

  py::enum_<runko::symbol>(m_sub, "symbol")
    .value("print", runko::symbol::print)
    .value("println", runko::symbol::println)
    .value("version", runko::symbol::version)
    .value("mt_showcase", runko::symbol::mt_showcase)
    .export_values();

  m_sub.def("empty_context_eval", &::empty_context_eval);
}
}  // namespace actions
