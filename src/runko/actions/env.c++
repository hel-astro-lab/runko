// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/env.h"

#include "tyvi/actions_list.h"

#include <pika/execution.hpp>
#include <print>
#include <ranges>
#include <string>
#include <thread>
#include <tuple>
#include <vector>


namespace runko {
namespace ta = tyvi::actions;
namespace te = tyvi::exec;
namespace rn = std::ranges;
namespace rv = std::views;

[[nodiscard]]
tyvi::actions::sexpr
  build_stdenv()
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

  return ta::list(
    ta::cons(symbol::print, ta::procedure { print }),
    ta::cons(symbol::println, ta::procedure { println }),
    ta::cons(symbol::version, version),
    ta::cons(symbol::mt_showcase, ta::procedure { mt_showcase }));
}
}  // namespace runko
