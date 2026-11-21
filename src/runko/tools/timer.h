// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <algorithm>
#include <chrono>
#include <map>
#include <mutex>
#include <print>
#include <string>
#include <tuple>
#include <vector>

namespace runko {

struct timer_duration {
  std::chrono::time_point<std::chrono::system_clock> begin, end;
};

/// Thread safe container for storing key-timer_duration pairs.
class duration_container {
  std::vector<std::string> keys_;
  std::vector<timer_duration> durations_;
  std::map<std::string, std::size_t> counts_;
  std::vector<std::string> metadata_;
  std::mutex mut_;

public:
  duration_container() = default;
  duration_container(duration_container&& other)
  {
    [[maybe_unused]] auto _ = std::lock_guard<std::mutex>(other.mut_);
    this->keys_             = std::move(other.keys_);
    this->durations_        = std::move(other.durations_);
    this->counts_           = std::move(other.counts_);
    this->metadata_         = std::move(other.metadata_);
  }
  duration_container& operator=(duration_container&&)      = delete;
  duration_container(const duration_container&)            = delete;
  duration_container& operator=(const duration_container&) = delete;

  void add(
    const std::string_view key,
    const timer_duration& dur,
    std::string metadata = "")
  {
    [[maybe_unused]] auto _ = std::lock_guard<std::mutex>(this->mut_);

    this->keys_.push_back(
      std::format("{} #{}", key, this->counts_[std::string { key }]++));
    this->durations_.push_back(dur);
    this->metadata_.push_back(std::move(metadata));
  }

  /// Returns the contained keys and durations and clears
  std::tuple<
    std::vector<std::string>,
    std::vector<timer_duration>,
    std::vector<std::string>>
    get_all()
  {
    [[maybe_unused]] auto _ = std::lock_guard<std::mutex>(this->mut_);

    auto ret = std::tuple { std::move(this->keys_),
                            std::move(this->durations_),
                            std::move(this->metadata_) };

    this->keys_.clear();
    this->durations_.clear();
    this->counts_.clear();
    this->metadata_.clear();
    return ret;
  }
};

}  // namespace runko
