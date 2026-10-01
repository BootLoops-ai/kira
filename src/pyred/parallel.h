/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)

This file is part of pyRed.

pyRed is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyRed is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyRed.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PYRED_PARALLEL_H
#define PYRED_PARALLEL_H

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pyred/defs.h"

namespace pyred {

/************************
 * Parallel distributor *
 ************************/

class cpu_count_error : public std::exception {
private:
  std::string msg;

public:
  inline cpu_count_error(const std::string &s) : msg(s) {}
  virtual inline const char *what() const noexcept { return msg.c_str(); }
};

inline int cpu_count() {
  // Get the number of physical CPU cores from /proc/cpuinfo.
  // If this failes (e.g. because /proc/cpuinfo does not exist on the platform),
  // return std::thread::hardware_concurrency().
  int n_physcores{0};
  try {
    std::ifstream instream{"/proc/cpuinfo"};
    if (!instream.good()) {
      throw cpu_count_error("Failed opening /proc/cpuinfo for reading.");
    }
    std::unordered_map<int, int> physid_cores;
    int physical_id{-1};
    std::string line;
    while (pyred::getline_normalized(instream, line)) {
      if (line.size() > 11 && line.substr(0, 11) == "physical id") {
        auto pos = line.find(':');
        if (pos >= line.size() - 1) throw cpu_count_error("cpu_count() error");
        physical_id = string_to_int(line.substr(pos + 1));
      }
      else if (line.size() > 9 && line.substr(0, 9) == "cpu cores") {
        auto pos = line.find(':');
        if (pos >= line.size() - 1) throw cpu_count_error("cpu_count() error");
        physid_cores.insert({physical_id, string_to_int(line.substr(pos + 1))});
      }
    }
    for (const auto &kv : physid_cores) {
      n_physcores += kv.second;
    }
  }
  catch (const std::exception &) {
    std::cout << "Failed determining physical cpu core count "
              << "from /proc/cpuinfo. "
              << "Using thread::hardware_concurrency() instead." << std::endl;
    n_physcores = std::thread::hardware_concurrency();
  }
  return n_physcores;
}

// Apply the function f concurrently to the elements of a vector.
// Accumulate the results into 'res' by calling acc(elem, res)
// on each processed element 'elem' in the order of the initial vector.
// parallel=0 chooses the number of threads automatically.
template <typename ITYPE, typename OTYPE, typename ACCTYPE>
class Distributor {
public:
  static ACCTYPE create(const std::function<OTYPE(const ITYPE &)> &,
                        const std::vector<ITYPE> &,
                        const std::function<void(OTYPE &&, ACCTYPE &)> &,
                        unsigned int = cpu_count());

private:
  // constructors
  Distributor(const std::function<OTYPE(const ITYPE &)> &,
              const std::vector<ITYPE> &, unsigned int,
              const std::function<void(OTYPE &&, ACCTYPE &)> &);
  // members
  std::function<OTYPE(const ITYPE &)> m_f;
  std::vector<ITYPE> m_in_items;
  std::vector<std::promise<OTYPE>> m_promises;
  std::function<void(OTYPE &&, ACCTYPE &)> m_acc;
  ACCTYPE m_out_data;
  unsigned int m_parallel;
  std::vector<std::thread> m_threads;
  std::mutex m_mtx;
  std::size_t m_pos;
  // methods
  void consume();
};

template <typename ITYPE, typename OTYPE, typename ACCTYPE>
Distributor<ITYPE, OTYPE, ACCTYPE>::Distributor(
    const std::function<OTYPE(const ITYPE &)> &f,
    const std::vector<ITYPE> &in_items, unsigned int parallel,
    const std::function<void(OTYPE &&, ACCTYPE &)> &acc)
    : m_f{f}, m_in_items{in_items}, m_acc{acc}, m_parallel{parallel}, m_pos{0} {
  if (!m_parallel) {
    m_parallel = cpu_count();
  }
  m_promises.resize(m_in_items.size());
  m_threads.reserve(m_parallel);
  for (unsigned int k = 0; k != m_parallel; ++k) {
    m_threads.push_back(std::thread(&Distributor::consume, this));
  }
  for (auto &prom : m_promises) {
    m_acc(std::move(prom.get_future().get()), m_out_data);
  }
  for (auto &t : m_threads) {
    t.join();
  }
}

template <typename ITYPE, typename OTYPE, typename ACCTYPE>
ACCTYPE Distributor<ITYPE, OTYPE, ACCTYPE>::create(
    const std::function<OTYPE(const ITYPE &)> &f,
    const std::vector<ITYPE> &in_items,
    const std::function<void(OTYPE &&, ACCTYPE &)> &acc,
    unsigned int parallel) {
  Distributor<ITYPE, OTYPE, ACCTYPE> dist{f, in_items, parallel, acc};
  return std::move(dist.m_out_data);
}

template <typename ITYPE, typename OTYPE, typename ACCTYPE>
void Distributor<ITYPE, OTYPE, ACCTYPE>::consume() {
  while (true) {
    std::size_t pos;
    {
      std::lock_guard<std::mutex> lock{m_mtx};
      if (m_pos == m_in_items.size()) {
        return;
      }
      pos = m_pos;
      ++m_pos;
    }
    m_promises[pos].set_value(m_f(m_in_items[pos]));
  }
}

// ============================================================= //
// New distributor for equation generator with on-the-fly solver //
// ============================================================= //

// parallel=0 chooses the number of threads automatically.
template <typename ITYPE, typename OTYPE, typename MASTER, typename AUX>
class Distributor2 {
public:
  // constructors
  Distributor2(const std::function<OTYPE(const ITYPE &, int, int, AUX &)> &,
               MASTER &, unsigned int);

private:
  // members
  std::function<OTYPE(const ITYPE &, int, int, AUX &)> m_f;
  std::vector<ITYPE> m_in_items;
  std::vector<std::promise<OTYPE>> m_promises;
  // m_master is a reference, so that the data in the callable object
  // can be accessed after when the distributor finished.
  MASTER &m_master;
  unsigned int m_parallel;
  std::vector<std::thread> m_threads;
  std::mutex m_mtx;
  std::condition_variable m_cv;
  std::size_t m_pos;
  bool m_finished;
  // methods
  void consume();
};

template <typename ITYPE, typename OTYPE, typename MASTER, typename AUX>
Distributor2<ITYPE, OTYPE, MASTER, AUX>::Distributor2(
    const std::function<OTYPE(const ITYPE &, int, int, AUX &)> &f,
    MASTER &master, unsigned int parallel)
    : m_f{f},
      m_master(master) // use () instead of {} because of a compiler bug
      ,
      m_parallel{parallel},
      m_pos{0},
      m_finished{false} {
  if (!m_parallel) {
    m_parallel = cpu_count();
  }
  m_promises.resize(m_in_items.size());
  m_threads.reserve(m_parallel);
  for (unsigned int k = 0; k != m_parallel; ++k) {
    m_threads.push_back(std::thread(&Distributor2::consume, this));
  }
  m_master(m_in_items, m_promises, m_mtx, m_cv, m_finished);
  for (auto &t : m_threads) {
    t.join();
  }
}

template <typename ITYPE, typename OTYPE, typename MASTER, typename AUX>
void Distributor2<ITYPE, OTYPE, MASTER, AUX>::consume() {
  AUX cache{};
  while (true) {
    std::size_t pos;
    {
      std::unique_lock<std::mutex> lock{m_mtx};
      // Wait until there are unprocessed items or m_finished has been set.
      m_cv.wait(lock,
                [this]() { return m_pos < m_in_items.size() || m_finished; });
      if (m_pos < m_in_items.size()) {
        // There are elements to process.
        pos = m_pos;
        ++m_pos;
      }
      else {
        // All elements have been processed and there will be no new ones.
        return;
      }
    }
    m_promises[pos].set_value(
        m_f(m_in_items[pos], pos + 1, m_promises.size(), cache));
  }
}

} // namespace pyred

#endif
