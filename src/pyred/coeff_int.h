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

#ifndef PYRED_COEFF_INT_H
#define PYRED_COEFF_INT_H

#include <iostream>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pyred/defs.h"

namespace pyred {

class Coeff_int {
  friend Coeff_int operator-(const Coeff_int&);
  friend Coeff_int operator+(const Coeff_int&, const Coeff_int&);
  friend Coeff_int operator-(const Coeff_int&, const Coeff_int&);
  friend Coeff_int operator*(const Coeff_int&, const Coeff_int&);
  friend Coeff_int operator/(const Coeff_int&, const Coeff_int&);
  friend Coeff_int pow(const Coeff_int&, const Coeff_int&);
  friend bool operator==(const Coeff_int&, const Coeff_int&);
  friend bool operator!=(const Coeff_int&, const Coeff_int&);
  friend std::ostream& operator<<(std::ostream&, const Coeff_int&);

private:
  uint64_t c;
  static uint64_t prime;
  static bool auto_vars;
  static std::string set_to_one;
  static std::vector<std::pair<std::string, uint64_t>> randvars;
  static std::mutex s_mtx;

public:
  static void set_prime(uint64_t tempPrime);
  static uint64_t get_prime();
  static std::string auxName;
  static void init();
  Coeff_int();
  Coeff_int(const Coeff_int&);
  Coeff_int(const std::string&);
  template <typename T, typename = typename std::enable_if<
                            std::is_integral<T>::value>::type>
  Coeff_int(T);
  Coeff_int& operator=(const Coeff_int&) = default;
  Coeff_int& operator=(Coeff_int&&) = default;
  Coeff_int& operator+=(const Coeff_int&);
  Coeff_int& operator-=(const Coeff_int&);
  Coeff_int& operator*=(const Coeff_int&);
  Coeff_int& operator/=(const Coeff_int&);
  operator bool() const;
  std::size_t hash() const { return c; }
};

Coeff_int operator-(const Coeff_int&);
Coeff_int operator+(const Coeff_int&, const Coeff_int&);
Coeff_int operator-(const Coeff_int&, const Coeff_int&);
Coeff_int operator*(const Coeff_int&, const Coeff_int&);
Coeff_int operator/(const Coeff_int&, const Coeff_int&);
Coeff_int pow(const Coeff_int&, const Coeff_int&);
bool operator==(const Coeff_int&, const Coeff_int&);
bool operator!=(const Coeff_int&, const Coeff_int&);
std::ostream& operator<<(std::ostream&, const Coeff_int&);

template <typename T, typename>
inline Coeff_int::Coeff_int(T a) {
  if (a >= 0) {
    c = uint64_t(a);
  }
  else {
    c = prime - uint64_t(-a);
  }
}

inline Coeff_int::operator bool() const { return c != 0; }

} // namespace pyred

namespace std {
template <>
struct hash<pyred::Coeff_int> {
  size_t operator()(const pyred::Coeff_int& cf) const { return cf.hash(); }
};
} // namespace std

#endif
