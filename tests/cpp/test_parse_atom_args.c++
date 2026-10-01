// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/actions/args.h"

#include <boost/ut.hpp>


namespace {
using namespace boost::ut;
namespace ta = tyvi::actions;
namespace te = tyvi::exec;

const auto s = [] {
  "one arg"_test = [] {
    const auto n =
      tyvi::this_thread::sync_wait(runko::parse_atom_args<int>(ta::list(42)));
    expect(n == 42);
  };

  "two args"_test = [] {
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, char const*>(ta::list(42, "43")) |
      te::then([](const int n, const std::string_view m) {
        expect(n == 42);
        expect(m == "43");
      }));
  };

  "one optional arg"_test = [] {
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>>(ta::list(42, 43)) |
      te::then([](const int n, const std::optional<int> m) {
        expect(n == 42);
        expect(m.value() == 43);
      }));

    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>>(ta::list(42, "43")) |
      te::then([](const int n, const std::optional<int> m) {
        expect(n == 42);
        expect(not m);
      }));

    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>>(ta::list(42)) |
      te::then([](const int n, const std::optional<int> m) {
        expect(n == 42);
        expect(not m);
      }));
  };

  "one opt arg"_test = [] {
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<runko::opt_arg<int>>(ta::list(42)) |
      te::then([](const std::optional<int> n) { expect(n.value() == 42); }));
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<runko::opt_arg<int>>(ta::list()) |
      te::then([](const std::optional<int> n) { expect(not n); }));
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<runko::opt_arg<int>>(ta::list("42")) |
      te::then([](const std::optional<int> n) { expect(not n); }));
  };

  "two optional args"_test = [] {
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>, runko::opt_arg<int>>(
        ta::list(42, 43, 44)) |
      te::then([](const int n, const std::optional<int> m, const std::optional<int> l) {
        expect(n == 42);
        expect(m.value() == 43);
        expect(l.value() == 44);
      }));
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>, runko::opt_arg<int>>(
        ta::list(42, "43", 44)) |
      te::then([](const int n, const std::optional<int> m, const std::optional<int> l) {
        expect(n == 42);
        expect(not m);
        expect(l.value() == 44);
      }));
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>, runko::opt_arg<int>>(
        ta::list(42, 43)) |
      te::then([](const int n, const std::optional<int> m, const std::optional<int> l) {
        expect(n == 42);
        expect(m.value() == 43);
        expect(not l);
      }));
    tyvi::this_thread::sync_wait(
      runko::parse_atom_args<int, runko::opt_arg<int>, runko::opt_arg<int>>(
        ta::list(42)) |
      te::then([](const int n, const std::optional<int> m, const std::optional<int> l) {
        expect(n == 42);
        expect(not m);
        expect(not l);
      }));
  };
};
}  // namespace

int
  main(int argc, const char** argv)
{
  [[maybe_unused]]
  const suite<"parse_atom_args"> _ = s;
  return static_cast<int>(cfg<override>.run(run_cfg { .argc = argc, .argv = argv }));
}
