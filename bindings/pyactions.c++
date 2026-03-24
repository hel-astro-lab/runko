#include "corgi/corgi.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
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

enum class symbol : std::uint8_t { print, println, version, mt_showcase };

ta::sexpr
  parse_element(const py::handle &obj)
{
  if(py::isinstance<py::int_>(obj)) {
    return obj.cast<long>();
  } else if(py::isinstance<py::str>(obj)) {
    return obj.cast<std::string>();
  } else if(py::isinstance<ta::intrinsic>(obj)) {
    return obj.cast<ta::intrinsic>();
  } else if(py::isinstance<symbol>(obj)) {
    return obj.cast<symbol>();
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
  print_context_eval(const py::handle &body_py)
{
  auto print = [](const ta::sexpr &args) -> ta::sexpr_sender {
    const auto arg_list = std::get<ta::cons>(args);
    const auto arg0     = std::get<ta::atom>(arg_list.car());
    const auto str      = ta::atom_cast<std::string>(arg0);

    return te::just(str.value()) | te::then([](const auto &str) {
             std::print("{}", str);
             return ta::null;
           });
  };

  auto println = [=](const ta::sexpr &args) -> ta::sexpr_sender {
    return print(args) | te::then([](auto &&) -> ta::sexpr {
             std::println();
             return ta::null;
           });
  };

  auto version = ta::atom(std::string { "runko v6.x" });

  auto mt_showcase = [](const ta::sexpr &args) -> ta::sexpr_sender {
    const auto arg_list = std::get<ta::cons>(args);
    const auto arg0     = std::get<ta::atom>(arg_list.car());
    const auto n        = ta::atom_cast<long>(arg0);

    auto senders = std::vector<te::unique_any_sender<>> {};

    for(const auto i: rv::iota(0l, n.value())) {
      senders.push_back(
        te::just(i) | te::continues_on(te::thread_pool_scheduler {}) |
        te::then([](const auto x) {
          std::this_thread::sleep_for(std::chrono::milliseconds { 500 * x });
          std::println("Hello from thread: {}", x);
        }));
    }

    return te::when_all_vector(std::move(senders)) | te::drop_value() |
           te::then([] { return ta::null; });
  };

  const auto body = parse_element(body_py);
  const auto env  = ta::list(
    ta::cons(symbol::print, ta::procedure { print }),
    ta::cons(symbol::println, ta::procedure { println }),
    ta::cons(symbol::version, version),
    ta::cons(symbol::mt_showcase, ta::procedure { mt_showcase }));

  pika::start(0, nullptr);

  try {
    tyvi::this_thread::sync_wait(ta::eval<symbol>(body, env));
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

  py::enum_<symbol>(m_sub, "symbol")
    .value("print", symbol::print)
    .value("println", symbol::println)
    .value("version", symbol::version)
    .value("mt_showcase", symbol::mt_showcase)
    .export_values();

  m_sub.def("print_context_eval", &::print_context_eval);
}
}  // namespace actions
