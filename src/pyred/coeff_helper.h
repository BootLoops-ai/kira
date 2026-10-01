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

#ifndef PYRED_COEFFICIENTS_H
#define PYRED_COEFFICIENTS_H

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "pyred/defs.h"

namespace pyred {

class Coeff_int;
class Coeff_vec;
template <int N>
class Coeff_arr;

class CoeffHelper {
  friend class Coeff_int;
  friend class Coeff_vec;
  template <int N>
  friend class Coeff_arr;

public:
  static bool is_valid_varname(const std::string &s);
  static void random_seed(int);
  static uint64_t random_minmax(uint64_t, uint64_t);
  static uint64_t mod_mul(uint64_t, uint64_t, uint64_t);
  static uint64_t mod_pow(uint64_t, uint64_t, uint64_t);
  static uint64_t mod_inv(uint64_t, uint64_t);
  static uint64_t parse_longint(const std::string &, uint64_t);
  static void add_invariant(const std::string &inv);
  static void settoone(const std::string &set1);
  static std::string settoone();
  static std::vector<std::string> invariants();
  static void clear_invariants() { s_invariants.clear(); }
  static void coeff_n(uint32_t);
  static uint32_t coeff_n();
  // defined in coeff_vec.cpp
  static std::vector<std::pair<int, int>> coeff_array_lengths();

private:
  static const std::vector<uint64_t> s_primes;
  static uint64_t randvar_lbound;
  static uint32_t s_coeff_n;
  static std::vector<std::string> s_invariants;
  static std::string s_settoone;
  static const std::vector<uint64_t> &primes();
};

} // namespace pyred

#endif
