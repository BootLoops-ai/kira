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

#ifndef KIRA_BLACK_BOX_H_
#define KIRA_BLACK_BOX_H_

#include <algorithm>
#include <chrono>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <firefly/BaseReconst.hpp>
#include <firefly/BlackBoxBase.hpp>
#include <firefly/FFInt.hpp>
#include <firefly/ShuntingYardParser.hpp>

#include "pyred/defs.h"
#include "pyred/gauss.h"
#include "kira/exceptions.h"
#include "kira/tools.h"

static Loginfo& loggerbb = Loginfo::instance();

class BlackBoxKira : public firefly::BlackBoxBase<BlackBoxKira> {
public:
  BlackBoxKira(const std::vector<pyred::Weight>& mandatory_vec,
               const uint32_t mode_ = 0) : mode(mode_) {
    for (const auto &id : mandatory_vec) {
      mandatory.emplace(id);
    }
  }

  std::vector<firefly::FFInt> operator()(const std::vector<firefly::FFInt>& values) {
    std::vector<std::pair<pyred::Weight, pyred::Weight>> assignment_tmp{};
    std::unordered_map<pyred::Weight, int> equation_lengths_tmp{};
    std::vector<double> times{};

    std::vector<firefly::FFInt> result =
        solve(values, assignment_tmp, equation_lengths_tmp, times);

    {
      std::lock_guard<std::mutex> lock(mut);

      if (iteration == 0) {
        assignment = std::move(assignment_tmp);
        equation_lengths = std::move(equation_lengths_tmp);

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
  std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
    std::vector<std::pair<pyred::Weight, pyred::Weight>> assignment_tmp{};
    std::unordered_map<pyred::Weight, int> equation_lengths_tmp{};
    std::vector<double> times{};

    std::vector<FFIntTemp> result =
        solve(values, assignment_tmp, equation_lengths_tmp, times);

    std::size_t size = 1;

    if (!result.empty()) {
      size = result.front().size();
    }

    {
      std::lock_guard<std::mutex> lock(mut);

      if (iteration == 0) {
        assignment = std::move(assignment_tmp);
        equation_lengths = std::move(equation_lengths_tmp);

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
    loggerbb << "\nAverage times for probe solution steps:\n";
    loggerbb << "Coefficient evaluation: " << std::setprecision(6) << parse_average << " s\n";
    loggerbb << "Forward elimination:    " << forward_average << " s\n";
    loggerbb << "Back substitution:      " << back_average << " s\n\n";
    loggerbb << std::setprecision(1);

    parser.precompute_tokens();
    parser_factors.precompute_tokens();
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

  inline void set_system(const std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>>& system_) {
    system = system_;
  }

  inline std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> get_system() const {
    return system;
  }

  inline void set_parser(firefly::ShuntingYardParser& par) {
    parser = std::move(par);
  }

  inline void force_precompute() {
    loggerbb << "\033[1;34mFireFly info:\033[0m Precomputing tokens\n";
    parser.precompute_tokens(true);
  }

  inline void set_factors(std::unordered_map<pyred::Weight, std::unordered_map<pyred::Weight, std::size_t>>& system, firefly::ShuntingYardParser& par) {
    multiply_factors = true;
    system_factors = std::move(system);
    parser_factors = std::move(par);
    parser_factors.precompute_tokens(true);
  }

  std::pair<std::vector<firefly::FFInt>, std::vector<std::vector<std::string>>> post_select(const std::vector<std::string>& symbols, const std::unordered_set<pyred::Weight>& masters_set_to_zero) {
    std::vector<firefly::FFInt> values;
    values.reserve(symbols.size());

    for (std::size_t i = 0; i != symbols.size(); ++i) {
      values.emplace_back(firefly::BaseReconst().get_rand_32());
    }

    std::vector<firefly::FFInt> funs = parser.evaluate(values);

    std::size_t old_size = system.size();

    pyred::SystemOfEqs<firefly::FFInt>::post_select(system, funs, masters_set_to_zero);

    loggerbb << "Post selection: Trimmed system from " << old_size << " to " << system.size() << " equations\n";

    std::vector<std::vector<std::string>> rpn_functions;
    parser.move_rpn(rpn_functions);
    trim_parser(symbols, rpn_functions);

    return std::make_pair(std::move(funs), std::move(rpn_functions));
  }

  void post_select(const std::vector<std::string>& symbols, const std::unordered_set<pyred::Weight>& masters_set_to_zero, const std::vector<firefly::FFInt>& funs, const std::vector<std::vector<std::string>>& rpn_functions) {
    std::size_t old_size = system.size();

    auto zeroweights = pyred::SystemOfEqs<firefly::FFInt>::post_select(system, funs, masters_set_to_zero);

    loggerbb << "Post selection: Trimmed system from " << old_size << " to " << system.size() << " equations\n";

    trim_parser(symbols, rpn_functions);
  }

  void prepare_backward() {
    loggerbb << "Prepare back substitution with FireFly\n";

    std::size_t old_size = system.size();

    // check if the system is in triangular form
    // the system should be sorted at this point
    if (old_size != 0) {
      for (std::size_t i = 0; i != old_size - 1; ++i) {
        if (system[i][0].first == system[i + 1][0].first) {
          throw ExceptionInternal("BlackBoxKira::prepare_backward(): System is not triangular");
        }
      }
    }

    if (!mandatory.empty()) {
      std::unordered_set<pyred::Weight> mandatory_tmp = mandatory;
      std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> new_sys_reverse;

      for (auto it = system.rbegin(); it != system.rend(); ++it) {
        auto found = mandatory_tmp.find(it->front().first);

        if (found != mandatory_tmp.end()) {
          for (auto itt = ++(it->begin()); itt != it->end(); ++itt) {
            mandatory_tmp.insert(itt->first);
          }

          new_sys_reverse.emplace_back(std::move(*it));
        }
      }

      system.clear();
      system.shrink_to_fit();
      system.reserve(new_sys_reverse.size());

      for (auto it = new_sys_reverse.rbegin(); it != new_sys_reverse.rend();
           ++it) {
        system.emplace_back(std::move(*it));
      }
    }

    if (old_size != system.size()) {
      loggerbb << "Trimmed system from " << old_size << " to " << system.size() << " equations\n";

      std::unordered_set<std::size_t> required_functions;

      for (const auto& eq : system) {
        for (const auto& term : eq) {
          required_functions.emplace(term.second);
        }
      }

      std::unordered_map<size_t, size_t> new_positions = parser.trim(required_functions);

      // Insert new positions into the system
      for (auto& eq : system) {
        for (auto& term : eq) {
          term.second = new_positions[term.second];
        }
      }
    }
  }

private:
  friend class Kira;

  // system of equations
  // first entry in pair is the integral ID and the second the position of the
  // coeffcient in parser
  std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> system {};
  firefly::ShuntingYardParser parser;
  std::unordered_map<pyred::Weight, std::unordered_map<pyred::Weight, std::size_t>> system_factors {};
  firefly::ShuntingYardParser parser_factors;
  std::unordered_set<pyred::Weight> mandatory;
  std::vector<std::pair<pyred::Weight, pyred::Weight>> assignment;
  std::unordered_map<pyred::Weight, int> equation_lengths;
  std::mutex mut;
  bool multiply_factors = false;
  uint32_t mode = 0;
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

  void trim_parser(const std::vector<std::string>& symbols, const std::vector<std::vector<std::string>>& rpn_functions) {
    std::unordered_set<std::size_t> required_functions;

    for (const auto& eq : system) {
      for (const auto& term : eq) {
        required_functions.emplace(term.second);
      }
    }

    parser = firefly::ShuntingYardParser(std::vector<std::string>{}, symbols);
    parser.reserve(required_functions.size());

    for (const auto& el : required_functions) {
      parser.add_otf_precompute(rpn_functions[el]);
    }

    parser.precompute_tokens();

    // Get new positions
    std::size_t counter = 0;
    std::unordered_map<size_t, size_t> new_positions {};
    for (const auto& el : required_functions) {
      new_positions.emplace(std::make_pair(el, counter));
      ++counter;
    }

    // Insert new positions into the system
    for (auto& eq : system) {
      for (auto& term : eq) {
        term.second = new_positions[term.second];
      }
    }

    loggerbb << "\033[1;34mFireFly info:\033[0m Trimmed functions from " << rpn_functions.size() << " to " << required_functions.size() << "\n";
  }

  template <typename FFIntTemp>
  std::vector<FFIntTemp> solve(
      const std::vector<FFIntTemp>& values,
      std::vector<std::pair<pyred::Weight, pyred::Weight>>& assignment,
      std::unordered_map<pyred::Weight, int>& equation_lengths,
      std::vector<double>& times) {
    auto time0 = std::chrono::high_resolution_clock::now();

    auto funs = parser.evaluate_pre(values);

    auto time1 = std::chrono::high_resolution_clock::now();

    auto numsys = pyred::SystemOfEqs<FFIntTemp>(system, funs, /*solve_otf*/ true);

    if (mode == 0) {
      numsys.sort();
    }

    auto time2 = std::chrono::high_resolution_clock::now();

    std::vector<FFIntTemp> result =
        backward(numsys, assignment, equation_lengths, values);

    auto time3 = std::chrono::high_resolution_clock::now();

    if (mode == 0) {
      times = {std::chrono::duration<double>(time1 - time0).count(),
               std::chrono::duration<double>(time2 - time1).count(),
               std::chrono::duration<double>(time3 - time2).count()};
    }
    else {
      times = {std::chrono::duration<double>(time1 - time0).count(),
               std::chrono::duration<double>(0).count(),
               std::chrono::duration<double>(time3 - time1).count()};
    }

    return result;
  }

  template <typename FFIntTemp>
  std::vector<FFIntTemp> backward(
      pyred::SystemOfEqs<FFIntTemp>& numsys,
      std::vector<std::pair<pyred::Weight, pyred::Weight>>& assignment,
      std::unordered_map<pyred::Weight, int>& equation_lengths,
      const std::vector<FFIntTemp>& values) {
    if (mode == 0 && !mandatory.empty()) {
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

    if (mode == 0) {
      numsys.solve(/*dosort*/ false);
    }

    std::vector<FFIntTemp> result{};

    if (!multiply_factors) {
      if (!mandatory.empty()) {
        for (auto it = numsys.sys.begin(); it != numsys.sys.end(); ++it) {
          auto found = mandatory.find(it->front().first);

          if (found != mandatory.end()) {
            equation_lengths.emplace(
                std::make_pair(it->front().first, it->size()));

            for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
              assignment.emplace_back(
                  std::make_pair(it->front().first, itt->first));
              result.emplace_back(itt->second);
            }
          }

          it->eq.clear();
        }
      }
      else {
        for (auto it = numsys.sys.begin(); it != numsys.sys.end(); ++it) {
          equation_lengths.emplace(
              std::make_pair(it->front().first, it->size()));

          for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
            assignment.emplace_back(
                std::make_pair(it->front().first, itt->first));
            result.emplace_back(itt->second);
          }

          it->eq.clear();
        }
      }
    }
    else {
      std::vector<FFIntTemp> factors = parser_factors.evaluate_pre(values);

      if (!mandatory.empty()) {
        for (auto it = numsys.sys.begin(); it != numsys.sys.end(); ++it) {
          auto found = mandatory.find(it->front().first);

          if (found != mandatory.end()) {
            equation_lengths.emplace(
                std::make_pair(it->front().first, it->size()));

            for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
              assignment.emplace_back(
                  std::make_pair(it->front().first, itt->first));

              auto it_fac = system_factors.find(it->front().first);

              if (it_fac != system_factors.end()) {
                auto itt_fac = it_fac->second.find(itt->first);

                if (itt_fac != it_fac->second.end()) {
                  result.emplace_back(itt->second / factors[itt_fac->second]);
                }
                else {
                  result.emplace_back(itt->second);
                }
              }
              else {
                result.emplace_back(itt->second);
              }
            }
          }

          it->eq.clear();
        }
      }
      else {
        for (auto it = numsys.sys.begin(); it != numsys.sys.end(); ++it) {
          equation_lengths.emplace(
              std::make_pair(it->front().first, it->size()));

          for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
            assignment.emplace_back(
                std::make_pair(it->front().first, itt->first));

            auto it_fac = system_factors.find(it->front().first);

            if (it_fac != system_factors.end()) {
              auto itt_fac = it_fac->second.find(itt->first);

              if (itt_fac != it_fac->second.end()) {
                result.emplace_back(itt->second / factors[itt_fac->second]);
              }
              else {
                result.emplace_back(itt->second);
              }
            }
            else {
              result.emplace_back(itt->second);
            }
          }

          it->eq.clear();
        }
      }
    }

    numsys.sys.clear();

    return result;
  }
};

#endif // KIRA_BLACK_BOX_H_
