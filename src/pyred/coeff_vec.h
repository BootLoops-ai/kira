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

#ifndef PYRED_COEFF_VEC_H
#define PYRED_COEFF_VEC_H

#include <array>
#include <iostream>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pyred/defs.h"
#include "pyred/ppmacros.h"

namespace pyred {

// ========= //
// Coeff_vec //
// ========= //

class Coeff_vec {
  friend Coeff_vec operator-(const Coeff_vec &);
  friend Coeff_vec operator+(const Coeff_vec &, const Coeff_vec &);
  friend Coeff_vec operator-(const Coeff_vec &, const Coeff_vec &);
  friend Coeff_vec operator*(const Coeff_vec &, const Coeff_vec &);
  friend Coeff_vec operator/(const Coeff_vec &, const Coeff_vec &);
  friend Coeff_vec pow(Coeff_vec &&, const Coeff_vec &);
  friend bool operator==(const Coeff_vec &, const Coeff_vec &);
  friend bool operator!=(const Coeff_vec &, const Coeff_vec &);
  friend std::ostream &operator<<(std::ostream &, const Coeff_vec &);

public:
  using coeff_t = std::vector<uint64_t>;

private:
  coeff_t c;
  static uint32_t ncoeffs; // number of coefficients in coefficient vector
  static coeff_t prime;
  static bool auto_vars;
  static std::string set_to_one;
  static std::vector<std::pair<std::string, coeff_t>> randvars;
  static std::mutex s_mtx;

public:
  static void init();
  Coeff_vec();
  Coeff_vec(const Coeff_vec &);
  Coeff_vec(Coeff_vec &&);
  Coeff_vec(const std::string &);
  template <typename T, typename = typename std::enable_if<
                            std::is_integral<T>::value>::type>
  Coeff_vec(T);
  Coeff_vec &operator=(const Coeff_vec &) = default;
  Coeff_vec &operator=(Coeff_vec &&) = default;
  Coeff_vec &operator+=(const Coeff_vec &);
  Coeff_vec &operator-=(const Coeff_vec &);
  Coeff_vec &operator*=(const Coeff_vec &);
  Coeff_vec &operator/=(const Coeff_vec &);
  operator bool() const;
  std::size_t hash() const { return c.front(); }
};

Coeff_vec operator-(const Coeff_vec &);
Coeff_vec operator+(const Coeff_vec &, const Coeff_vec &);
Coeff_vec operator-(const Coeff_vec &, const Coeff_vec &);
Coeff_vec operator*(const Coeff_vec &, const Coeff_vec &);
Coeff_vec operator/(const Coeff_vec &, const Coeff_vec &);
Coeff_vec pow(Coeff_vec &&, const Coeff_vec &);
inline Coeff_vec pow(const Coeff_vec &base, const Coeff_vec &exp) {
  return pow(Coeff_vec{base}, exp);
}
bool operator==(const Coeff_vec &, const Coeff_vec &);
bool operator!=(const Coeff_vec &, const Coeff_vec &);
std::ostream &operator<<(std::ostream &, const Coeff_vec &);

template <typename T, typename>
inline Coeff_vec::Coeff_vec(T a) {
  if (a >= 0) {
    c.resize(ncoeffs, a);
  }
  else {
    auto cc = uint64_t(-a);
    c.reserve(ncoeffs);
    for (const auto &p : prime) {
      c.emplace_back(p - cc);
    }
  }
}

inline Coeff_vec::operator bool() const {
  for (const auto &cc : c) {
    if (cc) return true;
  }
  return false;
}

// ========= //
// Coeff_arr //
// ========= //

template <int N>
class Coeff_arr {
  template <int NN>
  friend Coeff_arr<NN> operator-(const Coeff_arr<NN> &);
  template <int NN>
  friend Coeff_arr<NN> operator+(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend Coeff_arr<NN> operator-(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend Coeff_arr<NN> operator*(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend Coeff_arr<NN> operator/(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend Coeff_arr<NN> pow(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend bool operator==(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend bool operator!=(const Coeff_arr<NN> &, const Coeff_arr<NN> &);
  template <int NN>
  friend std::ostream &operator<<(std::ostream &, const Coeff_arr<NN> &);

public:
  using coeff_t = std::array<uint64_t, N>;

private:
  coeff_t c;
  static coeff_t prime;
  static bool auto_vars;
  static std::string set_to_one;
  static std::vector<std::pair<std::string, coeff_t>> randvars;
  static std::mutex s_mtx;

public:
  static void init();
  Coeff_arr() {}
  Coeff_arr(const Coeff_arr &a) : c{a.c} {}
  // Move constructor is pointless for std::array.
  Coeff_arr(const std::string &);
  template <typename T, typename = typename std::enable_if<
                            std::is_integral<T>::value>::type>
  Coeff_arr(T);
  Coeff_arr &operator=(const Coeff_arr &) = default;
  Coeff_arr &operator=(Coeff_arr &&) = default;
  Coeff_arr &operator+=(const Coeff_arr &);
  Coeff_arr &operator-=(const Coeff_arr &);
  Coeff_arr &operator*=(const Coeff_arr &);
  Coeff_arr &operator/=(const Coeff_arr &);
  operator bool() const;
  std::size_t hash() const { return c.front(); }
};

// static members
template <int N>
bool Coeff_arr<N>::auto_vars{true};
template <int N>
std::string Coeff_arr<N>::set_to_one{};
template <int N>
typename Coeff_arr<N>::coeff_t Coeff_arr<N>::prime{};
template <int N>
std::vector<std::pair<std::string, typename Coeff_arr<N>::coeff_t>>
    Coeff_arr<N>::randvars{};
template <int N>
std::mutex Coeff_arr<N>::s_mtx{};

template <int N>
template <typename T, typename>
inline Coeff_arr<N>::Coeff_arr(T a) {
  if (a >= 0) {
    c.fill(a);
  }
  else {
    auto cc = uint64_t(-a);
    for (int k = 0; k != N; ++k) {
      c[k] = prime[k] - cc;
    }
  }
}

template <int N>
inline Coeff_arr<N>::operator bool() const {
  for (const auto &cc : c) {
    if (cc) return true;
  }
  return false;
}

template <int N>
Coeff_arr<N> operator-(const Coeff_arr<N> &);

template <int N>
Coeff_arr<N> operator+(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
Coeff_arr<N> operator-(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
Coeff_arr<N> operator*(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
Coeff_arr<N> operator/(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
bool operator==(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
bool operator!=(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
Coeff_arr<N> pow(const Coeff_arr<N> &, const Coeff_arr<N> &);

template <int N>
std::ostream &operator<<(std::ostream &, const Coeff_arr<N> &);

// clang-format off
// Explicit template instantiation declarations for Coeff_arr
// extern template class Coeff_arr<2>;
// extern template Coeff_arr<2> operator-(const Coeff_arr<2> &);
// extern template Coeff_arr<2> operator+(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template Coeff_arr<2> operator-(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template Coeff_arr<2> operator*(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template Coeff_arr<2> operator/(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template bool operator==(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template bool operator!=(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template Coeff_arr<2> pow(const Coeff_arr<2> &, const Coeff_arr<2> &);
// extern template std::ostream & operator<<(std::ostream &, const Coeff_arr<2> &);
# define PYRED_PP_COEFFARR_EXPLICIT_DECL(k) \
extern template class PYRED_PP_COEFFCLASS(k); \
extern template PYRED_PP_COEFFCLASS(k) operator-(const PYRED_PP_COEFFCLASS(k) &); \
extern template PYRED_PP_COEFFCLASS(k) operator+(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template PYRED_PP_COEFFCLASS(k) operator-(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template PYRED_PP_COEFFCLASS(k) operator*(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template PYRED_PP_COEFFCLASS(k) operator/(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template bool operator==(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template bool operator!=(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template PYRED_PP_COEFFCLASS(k) pow(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
extern template std::ostream & operator<<(std::ostream &, const PYRED_PP_COEFFCLASS(k) &);
PYRED_PP_REPEATCARR(PYRED_PP_NCOEFFCLASSES, PYRED_PP_COEFFARR_EXPLICIT_DECL)
// clang-format on

} // namespace pyred

namespace std {

template <>
struct hash<pyred::Coeff_vec> {
  size_t operator()(const pyred::Coeff_vec &cf) const { return cf.hash(); }
};

template <int N>
struct hash<pyred::Coeff_arr<N>> {
  size_t operator()(const pyred::Coeff_arr<N> &cf) const { return cf.hash(); }
};

} // namespace std

#endif
