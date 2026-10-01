/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)

This file is part of Kira.

Kira is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Kira is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Kira.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef KIRA_BLACK_BOX_NUMERICAL_H_
#define KIRA_BLACK_BOX_NUMERICAL_H_

#include <chrono>
#include <mutex>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <firefly/BlackBoxBase.hpp>
#include <firefly/FFInt.hpp>
#include <firefly/ShuntingYardParser.hpp>

#include "pyred/defs.h"
#include "pyred/gauss.h"
#include "kira/exceptions.h"
#include "kira/tools.h"

// TODO factors? -> later, throw error for now

class BlackBoxNumerical : public firefly::BlackBoxBase<BlackBoxNumerical> {
public:
  BlackBoxNumerical(const std::vector<pyred::Weight>& mandatory_vec) {
    for (const auto &id : mandatory_vec) {
      mandatory.emplace(id);
    }
  }

  pyred::SystemOfEqs<firefly::FFInt> operator()(const std::vector<firefly::FFInt>& values) {
    std::vector<double> times{};

    pyred::SystemOfEqs<firefly::FFInt> result = solve(values, times);

    {
      std::lock_guard<std::mutex> lock(mut);

      if (iteration == 0) {
        parse_average = times[0];
        forward_average = times[1];
        back_average = times[2];
      }
      else {
        parse_average =
            (parse_average * iteration + times[0]) / (iteration + 1);
        forward_average =
            (forward_average * iteration + times[1]) / (iteration + 1);
        back_average = (back_average * iteration + times[2]) / (iteration + 1);
      }

      ++iteration;
    }

    return result;
  }

  template <typename FFIntTemp>
  pyred::SystemOfEqs<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
    std::vector<double> times{};
    std::size_t size = 1;

    pyred::SystemOfEqs<FFIntTemp> result = solve(values, times);

    if (!result.empty()) {
      size = result.begin()->front().second.size();
    }

    {
      std::lock_guard<std::mutex> lock(mut);

      if (iteration == 0) {
        parse_average = times[0] / size;
        forward_average = times[1] / size;
        back_average = times[2] / size;
      }
      else {
        parse_average =
            (parse_average * iteration + times[0]) / (iteration + size);
        forward_average =
            (forward_average * iteration + times[1]) / (iteration + size);
        back_average =
            (back_average * iteration + times[2]) / (iteration + size);
      }

      iteration += size;
    }

    return result;
  }

  inline void prime_changed() {
    parser.precompute_tokens();
  }

  inline void reserve(const std::size_t size) {
    system.reserve(size);
  }

  inline void add_eqn(std::vector<std::pair<pyred::Weight, std::size_t>>& eqn) {
    system.emplace_back(std::move(eqn));
  }

  inline void sort_system() {
    std::sort(system.begin(), system.end(), cmp_eqn);

    if (!system.empty()) {
      while (system.back().empty()) {
        system.pop_back();
      }
    }
  }

  inline void set_parser(firefly::ShuntingYardParser& par) {
    parser = std::move(par);
  }

  inline void force_precompute() {
    parser.precompute_tokens(true);
  }

private:
  friend class Kira;

  // system of equations
  // first entry in pair is the integral ID and the second the position of the
  // coefficient in parser
  std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> system {};
  firefly::ShuntingYardParser parser;
  std::unordered_set<pyred::Weight> mandatory;
  std::mutex mut;
  std::size_t iteration = 0;
  double parse_average = 0;
  double forward_average = 0;
  double back_average = 0;

  // cmp_eqn adapted from pyRed
  static bool cmp_eqn(const std::vector<std::pair<pyred::Weight, std::size_t>>& a, const std::vector<std::pair<pyred::Weight, std::size_t>>& b) {
    /*
    Compare Equations (like operator<) by
    * highest integral: lower first
    * length: shorter first
    * if same highest integral and same length:
      first if lower integrals following
    * if all integrals the same:
      lower equation number first.
    * if an equation is empty, place it last.
    Makes std::sort place lower equations first.
    */
    if (a.empty()) return false;
    if (b.empty()) return true;
    if (a.front().first != b.front().first) {
      return a.front().first < b.front().first;
    }
    if (a.size() != b.size()) return a.size() < b.size();
    for (std::size_t i = 1; i != a.size(); ++i) {
      if (a[i].first != b[i].first) {
        return a[i].first < b[i].first;
      }
    }
    return false;
  }

  template <typename FFIntTemp>
  pyred::SystemOfEqs<FFIntTemp> solve(
      const std::vector<FFIntTemp>& values,
      std::vector<double>& times) {
    auto time0 = std::chrono::high_resolution_clock::now();

    auto funs = parser.evaluate_pre(values);

    auto time1 = std::chrono::high_resolution_clock::now();

    auto numsys = pyred::SystemOfEqs<FFIntTemp>(system, funs, /*solve_otf*/ true);

    numsys.sort();

    auto time2 = std::chrono::high_resolution_clock::now();

    backward(numsys);

    auto time3 = std::chrono::high_resolution_clock::now();

    times = {std::chrono::duration<double>(time1 - time0).count(),
             std::chrono::duration<double>(time2 - time1).count(),
             std::chrono::duration<double>(time3 - time2).count()};

    return numsys;
  }

  template <typename FFIntTemp>
  void backward(pyred::SystemOfEqs<FFIntTemp>& numsys) {
    if (!mandatory.empty()) {
      std::unordered_set<pyred::Weight> mandatory_tmp = mandatory;
      std::vector<pyred::Equation<FFIntTemp>> new_sys_reverse;

      for (auto it = numsys.sys.rbegin(); it != numsys.sys.rend(); ++it) {
        auto found = mandatory_tmp.find(it->front().first);

        if (found != mandatory_tmp.end()) {
          for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
            mandatory_tmp.insert(itt->first);
          }

          new_sys_reverse.emplace_back(std::move(*it));
        }
      }

      numsys.sys.clear();
      numsys.sys.shrink_to_fit();
      numsys.sys.reserve(new_sys_reverse.size());

      for (auto it = new_sys_reverse.rbegin(); it != new_sys_reverse.rend();
           ++it) {
        numsys.sys.emplace_back(std::move(*it));
      }
    }

    numsys.solve(/*dosort*/ false);
  }
};

#endif // KIRA_BLACK_BOX_NUMERICAL_H_
