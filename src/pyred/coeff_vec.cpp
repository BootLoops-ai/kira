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

#include <algorithm>
#include <sstream>

#include "pyred/coeff_helper.h" // mod_mul, mod_inv
#include "pyred/coeff_vec.h"

namespace pyred {

// ========= //
// Coeff_vec //
// ========= //

// static members
bool Coeff_vec::auto_vars{true};
std::string Coeff_vec::set_to_one{};
uint32_t Coeff_vec::ncoeffs{2};
Coeff_vec::coeff_t Coeff_vec::prime{CoeffHelper::primes().at(0),
                                    CoeffHelper::primes().at(1)};
std::vector<std::pair<std::string, Coeff_vec::coeff_t>> Coeff_vec::randvars{};
std::mutex Coeff_vec::s_mtx{};

void Coeff_vec::init() {
  auto n = CoeffHelper::s_coeff_n;
  if (n < 1) {
    std::ostringstream ss;
    ss << "Invalid number of coefficients in vector: " << n << std::endl;
    throw input_error(ss.str());
  }
  if (!CoeffHelper::s_invariants.empty()) {
    // If at lease one variable was set manually,
    // turn off auto variable declaration.
    auto_vars = false;
    if (!CoeffHelper::s_settoone.empty() &&
        std::find(CoeffHelper::s_invariants.cbegin(),
                  CoeffHelper::s_invariants.cend(), CoeffHelper::s_settoone) ==
            CoeffHelper::s_invariants.cend()) {
      // If CoeffHelper::s_settoone is set and auto varible declaration is off,
      // CoeffHelper::s_settoone must be a declared variable.
      std::ostringstream ss;
      ss << "'" << CoeffHelper::s_settoone << "' is not declared as a variable."
         << std::endl;
      throw input_error(ss.str());
    }
  }
  set_to_one = CoeffHelper::s_settoone;
  ncoeffs = n;
  prime.clear();
  prime.reserve(n);
  auto pit = CoeffHelper::primes().cbegin();
  for (uint32_t k = 0; k != n; ++k) {
    prime.emplace_back(*pit);
    // If more coefficients are requested as primes are available,
    // use primes multiply with different random numbers.
    if (++pit == CoeffHelper::primes().cend()) {
      pit = CoeffHelper::primes().cbegin();
    }
  }
  randvars.clear();
  randvars.reserve(CoeffHelper::s_invariants.size());
  for (const auto& var : CoeffHelper::s_invariants) {
    randvars.push_back({var, coeff_t()});
    randvars.back().second.reserve(n);
  }
  for (const auto& p : prime) {
    for (auto& vv : randvars) {
      if (vv.first == CoeffHelper::s_settoone) {
        vv.second.push_back(1);
      }
      else {
        vv.second.push_back(
            CoeffHelper::random_minmax(CoeffHelper::randvar_lbound, p));
      }
    }
  }
  return;
}

Coeff_vec::Coeff_vec() : c(coeff_t()) {}
Coeff_vec::Coeff_vec(const Coeff_vec& a) : c{a.c} {}
Coeff_vec::Coeff_vec(Coeff_vec&& a) : c{std::move(a.c)} {}

Coeff_vec::Coeff_vec(const std::string& s) {
  if (auto_vars) s_mtx.lock();
  for (const auto& vc : randvars) {
    if (vc.first == s) {
      c = vc.second;
      if (auto_vars) s_mtx.unlock();
      return;
    }
  }
  if (auto_vars) s_mtx.unlock();
  std::istringstream ss{s};
  uint64_t cc;
  c.reserve(ncoeffs);
  auto success = static_cast<bool>(ss >> cc);
  if (!(success && ss.rdbuf()->in_avail() == 0)) {
    if (s.empty()) {
      // (ss >> c) fails if ss is empty
      throw parser_error("Coeff_vec(): empty argument");
    }
    else if (std::isalpha(s.front())) {
      if (!auto_vars) {
        throw parser_error(std::string("invalid coefficient string \"") + s +
                           "\"");
      }
      // auto_vars: add new symbols automatically and assign values to them.
      if (!CoeffHelper::is_valid_varname(s)) {
        throw parser_error(s + std::string("; symbols are expected to be of "
                                           "the form [A-Za-z_][A-Za-z0-9_]*."));
      }
      for (const auto& p : prime) {
        if (s == set_to_one) {
          c.push_back(1);
        }
        else {
          c.push_back(
              CoeffHelper::random_minmax(CoeffHelper::randvar_lbound, p));
        }
      }
      std::lock_guard<std::mutex> lock(s_mtx);
      CoeffHelper::s_invariants.push_back(s);
      randvars.push_back({s, c});
      // std::cout << "define symbol " << s << std::endl;
    }
    else {
      for (const auto& p : prime) {
        c.push_back(CoeffHelper::parse_longint(s, p));
      }
    }
  }
  else {
    for (const auto& p : prime) {
      if (cc >= p) {
        uint64_t ccc = cc % p;
        c.push_back(ccc);
      }
      else {
        c.push_back(cc);
      }
    }
  }
}

Coeff_vec& Coeff_vec::operator+=(const Coeff_vec& a) {
  auto ait = a.c.cbegin();
  auto pit = prime.cbegin();
  for (auto& cc : c) {
    cc += *ait++;
    if (cc >= *pit) cc -= *pit;
    ++pit;
  }
  return *this;
}
Coeff_vec& Coeff_vec::operator-=(const Coeff_vec& a) {
  auto ait = a.c.cbegin();
  auto cit = c.begin();
  for (const auto& p : prime) {
    if (*ait > *cit) *cit += p;
    *cit++ -= *ait++;
  }
  return *this;
}
Coeff_vec& Coeff_vec::operator*=(const Coeff_vec& a) {
  auto ait = a.c.cbegin();
  auto pit = prime.cbegin();
  for (auto& cc : c) {
    cc = CoeffHelper::mod_mul(cc, *ait++, *pit++);
  }
  return *this;
}
Coeff_vec& Coeff_vec::operator/=(const Coeff_vec& a) {
  auto ait = a.c.cbegin();
  auto pit = prime.cbegin();
  for (auto& cc : c) {
    if (*ait == static_cast<uint64_t>(0)) {
      throw zero_modular_division();
    }
    cc = CoeffHelper::mod_mul(cc, CoeffHelper::mod_inv(*ait++, *pit), *pit);
    ++pit;
  }
  return *this;
}
Coeff_vec pow(Coeff_vec&& b, const Coeff_vec& a) {
  Coeff_vec result = std::move(b);
  auto ait = a.c.cbegin();
  auto pit = Coeff_vec::prime.cbegin();
  for (auto& cc : result.c) {
    if (*ait == 2u) {
      cc = CoeffHelper::mod_mul(cc, cc, *pit);
    }
    else {
      uint64_t exp;
      uint64_t base;
      if (*ait < (*pit >> 1)) {
        // treat as positive exponent
        exp = *ait;
        base = cc;
      }
      else {
        // treat as negative exponent
        exp = *pit - *ait;
        base = CoeffHelper::mod_inv(cc, *pit); // =1/cc
      }
      cc = CoeffHelper::mod_pow(base, exp, *pit);
    }
    ++ait;
    ++pit;
  }
  return result;
}

Coeff_vec operator-(const Coeff_vec& a) {
  Coeff_vec cv;
  cv.c.reserve(Coeff_vec::ncoeffs);
  auto pit = Coeff_vec::prime.cbegin();
  for (const auto& aa : a.c) {
    if (aa == 0) {
      cv.c.emplace_back(aa);
      ++pit;
    }
    else {
      cv.c.emplace_back(*pit++ - aa);
    }
  }
  return cv;
}
Coeff_vec operator+(const Coeff_vec& a, const Coeff_vec& b) {
  Coeff_vec cv;
  cv.c.reserve(Coeff_vec::ncoeffs);
  auto ait = a.c.cbegin();
  auto bit = b.c.cbegin();
  for (const auto& p : Coeff_vec::prime) {
    auto sum = *ait++ + *bit++;
    if (sum >= p) sum -= p;
    cv.c.emplace_back(sum);
  }
  return cv;
}
Coeff_vec operator-(const Coeff_vec& a, const Coeff_vec& b) {
  Coeff_vec cv;
  cv.c.reserve(Coeff_vec::ncoeffs);
  auto ait = a.c.cbegin();
  auto pit = Coeff_vec::prime.cbegin();
  for (const auto& bb : b.c) {
    auto diff = *ait++;
    if (bb > diff) diff += *pit++;
    cv.c.emplace_back(diff - bb);
  }
  return cv;
}
Coeff_vec operator*(const Coeff_vec& a, const Coeff_vec& b) {
  Coeff_vec cv;
  cv.c.reserve(Coeff_vec::ncoeffs);
  auto ait = a.c.cbegin();
  auto bit = b.c.cbegin();
  for (const auto& p : Coeff_vec::prime) {
    cv.c.emplace_back(CoeffHelper::mod_mul(*ait++, *bit++, p));
  }
  return cv;
}
Coeff_vec operator/(const Coeff_vec& a, const Coeff_vec& b) {
  Coeff_vec cv;
  cv.c.reserve(Coeff_vec::ncoeffs);
  auto ait = a.c.cbegin();
  auto bit = b.c.cbegin();
  for (const auto& p : Coeff_vec::prime) {
    if (*bit == uint64_t(0)) {
      throw zero_modular_division();
    }
    cv.c.emplace_back(
        CoeffHelper::mod_mul(*ait++, CoeffHelper::mod_inv(*bit++, p), p));
  }
  return cv;
}

bool operator==(const Coeff_vec& a, const Coeff_vec& b) { return (a.c == b.c); }

bool operator!=(const Coeff_vec& a, const Coeff_vec& b) { return (a.c != b.c); }

std::ostream& operator<<(std::ostream& out, const Coeff_vec& a) {
  out << "{";
  if (!a.c.empty()) {
    out << a.c.front();
    for (auto it = a.c.cbegin() + 1; it != a.c.cend(); ++it) {
      out << "," << *it;
    }
  }
  out << "}";
  return out;
}

// ========= //
// Coeff_arr //
// ========= //

template <int N>
void Coeff_arr<N>::init() {
  if (!CoeffHelper::s_invariants.empty()) {
    // If at lease one variable was set manually,
    // turn off auto variable declaration.
    auto_vars = false;
    if (!CoeffHelper::s_settoone.empty() &&
        std::find(CoeffHelper::s_invariants.cbegin(),
                  CoeffHelper::s_invariants.cend(), CoeffHelper::s_settoone) ==
            CoeffHelper::s_invariants.cend()) {
      // If CoeffHelper::s_settoone is set and auto varible declaration is off,
      // CoeffHelper::s_settoone must be a declared variable.
      std::ostringstream ss;
      ss << "'" << CoeffHelper::s_settoone << "' is not declared as a variable."
         << std::endl;
      throw input_error(ss.str());
    }
  }
  set_to_one = CoeffHelper::s_settoone;
  auto pit = CoeffHelper::primes().cbegin();
  for (int k = 0; k != N; ++k) {
    prime[k] = *pit;
    // If more coefficients are requested as primes are available,
    // use primes multiply with different random numbers.
    if (++pit == CoeffHelper::primes().cend()) {
      pit = CoeffHelper::primes().cbegin();
    }
  }
  randvars.clear();
  randvars.reserve(CoeffHelper::s_invariants.size());
  for (const auto& var : CoeffHelper::s_invariants) {
    randvars.push_back({var, coeff_t()});
  }
  for (int k = 0; k != N; ++k) {
    for (auto& vv : randvars) {
      if (vv.first == CoeffHelper::s_settoone) {
        vv.second[k] = 1;
      }
      else {
        vv.second[k] =
            CoeffHelper::random_minmax(CoeffHelper::randvar_lbound, prime[k]);
      }
    }
  }
}

template <int N>
Coeff_arr<N>::Coeff_arr(const std::string& s) {
  if (auto_vars) s_mtx.lock();
  for (const auto& vc : randvars) {
    if (vc.first == s) {
      c = vc.second;
      if (auto_vars) s_mtx.unlock();
      return;
    }
  }
  if (auto_vars) s_mtx.unlock();
  std::istringstream ss{s};
  uint64_t cc;
  auto success = static_cast<bool>(ss >> cc);
  if (!(success && ss.rdbuf()->in_avail() == 0)) {
    if (s.empty()) {
      // (ss >> c) fails if ss is empty
      throw parser_error("Coeff_arr(): empty argument");
    }
    else if (std::isalpha(s.front())) {
      if (!auto_vars) {
        throw parser_error(std::string("invalid coefficient string \"") + s +
                           "\"");
      }
      // auto_vars: add new symbols automatically and assign values to them.
      if (!CoeffHelper::is_valid_varname(s)) {
        throw parser_error(s + std::string("; symbols are expected to be of "
                                           "the form [A-Za-z_][A-Za-z0-9_]*."));
      }
      if (s == set_to_one) {
        c.fill(1u);
      }
      else {
        for (int k = 0; k != N; ++k) {
          c[k] =
              CoeffHelper::random_minmax(CoeffHelper::randvar_lbound, prime[k]);
        }
      }
      std::lock_guard<std::mutex> lock(s_mtx);
      CoeffHelper::s_invariants.push_back(s);
      randvars.push_back({s, c});
      // std::cout << "define symbol " << s << std::endl;
    }
    else {
      for (int k = 0; k != N; ++k) {
        c[k] = CoeffHelper::parse_longint(s, prime[k]);
      }
    }
  }
  else {
    for (int k = 0; k != N; ++k) {
      c[k] = (cc % prime[k]);
    }
  }
}

template <int N>
Coeff_arr<N>& Coeff_arr<N>::operator+=(const Coeff_arr<N>& a) {
  for (int k = 0; k != N; ++k) {
    c[k] += a.c[k];
    if (c[k] >= prime[k]) c[k] -= prime[k];
  }
  return *this;
}

template <int N>
Coeff_arr<N>& Coeff_arr<N>::operator-=(const Coeff_arr<N>& a) {
  for (int k = 0; k != N; ++k) {
    if (a.c[k] > c[k]) c[k] += prime[k];
    c[k] -= a.c[k];
  }
  return *this;
}

template <int N>
Coeff_arr<N>& Coeff_arr<N>::operator*=(const Coeff_arr<N>& a) {
  for (int k = 0; k != N; ++k) {
    c[k] = CoeffHelper::mod_mul(c[k], a.c[k], prime[k]);
  }
  return *this;
}

template <int N>
Coeff_arr<N>& Coeff_arr<N>::operator/=(const Coeff_arr<N>& a) {
  for (int k = 0; k != N; ++k) {
    if (a.c[k] == 0u) {
      throw zero_modular_division();
    }
    c[k] = CoeffHelper::mod_mul(c[k], CoeffHelper::mod_inv(a.c[k], prime[k]),
                                prime[k]);
  }
  return *this;
}

template <int N>
Coeff_arr<N> pow(const Coeff_arr<N>& b, const Coeff_arr<N>& e) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    if (e.c[k] == 2u) {
      result.c[k] =
          CoeffHelper::mod_mul(b.c[k], b.c[k], Coeff_arr<N>::prime[k]);
    }
    else {
      if (e.c[k] < (Coeff_arr<N>::prime[k] >> 1)) {
        // treat as positive exponent
        result.c[k] =
            CoeffHelper::mod_pow(b.c[k], e.c[k], Coeff_arr<N>::prime[k]);
      }
      else {
        // treat as negative exponent
        result.c[k] = CoeffHelper::mod_pow(
            CoeffHelper::mod_inv(b.c[k], Coeff_arr<N>::prime[k]),
            Coeff_arr<N>::prime[k] - e.c[k], Coeff_arr<N>::prime[k]);
      }
    }
  }
  return result;
}

template <int N>
Coeff_arr<N> operator-(const Coeff_arr<N>& a) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    if (a.c[k] == 0u) {
      result.c[k] = 0u;
    }
    else {
      result.c[k] = Coeff_arr<N>::prime[k] - a.c[k];
    }
  }
  return result;
}

template <int N>
Coeff_arr<N> operator+(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    result = a.c[k] + b.c[k];
    if (result.c[k] >= Coeff_arr<N>::prime[k]) {
      result.c[k] -= Coeff_arr<N>::prime[k];
    }
  }
  return result;
}

template <int N>
Coeff_arr<N> operator-(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    if (b.c[k] <= a.c[k]) {
      result.c[k] = a.c[k] - b.c[k];
    }
    else {
      result.c[k] = Coeff_arr<N>::prime[k] + a.c[k] - b.c[k];
    }
  }
  return result;
}

template <int N>
Coeff_arr<N> operator*(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    result.c[k] = CoeffHelper::mod_mul(a.c[k], b.c[k], Coeff_arr<N>::prime[k]);
  }
  return result;
}

template <int N>
Coeff_arr<N> operator/(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  Coeff_arr<N> result;
  for (int k = 0; k != N; ++k) {
    if (b.c[k] == 0u) {
      throw zero_modular_division();
    }
    result.c[k] = CoeffHelper::mod_mul(
        a.c[k], CoeffHelper::mod_inv(b.c[k], Coeff_arr<N>::prime[k]),
        Coeff_arr<N>::prime[k]);
  }
  return result;
}

template <int N>
bool operator==(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  return (a.c == b.c);
}

template <int N>
bool operator!=(const Coeff_arr<N>& a, const Coeff_arr<N>& b) {
  return (a.c != b.c);
}

template <int N>
std::ostream& operator<<(std::ostream& out, const Coeff_arr<N>& a) {
  out << "{" << a.c.front();
  for (auto it = a.c.cbegin() + 1; it != a.c.cend(); ++it) {
    out << "," << *it;
  }
  out << "}";
  return out;
}

// clang-format off
// Explicit template instantiation definitions for Coeff_arr
// template class Coeff_arr<2>;
// template Coeff_arr<2> operator-(const Coeff_arr<2> &);
// template Coeff_arr<2> operator+(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template Coeff_arr<2> operator-(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template Coeff_arr<2> operator*(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template Coeff_arr<2> operator/(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template bool operator==(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template bool operator!=(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template Coeff_arr<2> pow(const Coeff_arr<2> &, const Coeff_arr<2> &);
// template std::ostream & operator<<(std::ostream &, const Coeff_arr<2> &);
# define PYRED_PP_COEFFARR_EXPLICIT_DEF(k) \
template class PYRED_PP_COEFFCLASS(k); \
template PYRED_PP_COEFFCLASS(k) operator-(const PYRED_PP_COEFFCLASS(k) &); \
template PYRED_PP_COEFFCLASS(k) operator+(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template PYRED_PP_COEFFCLASS(k) operator-(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template PYRED_PP_COEFFCLASS(k) operator*(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template PYRED_PP_COEFFCLASS(k) operator/(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template bool operator==(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template bool operator!=(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template PYRED_PP_COEFFCLASS(k) pow(const PYRED_PP_COEFFCLASS(k) &, const PYRED_PP_COEFFCLASS(k) &); \
template std::ostream & operator<<(std::ostream &, const PYRED_PP_COEFFCLASS(k) &);
PYRED_PP_REPEATCARR(PYRED_PP_NCOEFFCLASSES, PYRED_PP_COEFFARR_EXPLICIT_DEF)
// clang-format on

std::vector<std::pair<int, int>> CoeffHelper::coeff_array_lengths() {
  // For each available array coefficient type a pair {coeff_cls, length}
  // of the coefficient class number and the array length.
  std::vector<std::pair<int, int>> cal;
#define PYRED_PP_COEFFARRSIZES(k) \
  cal.push_back(                  \
      {k + 1, std::tuple_size<PYRED_PP_COEFFCLASS(k)::coeff_t>::value});
  PYRED_PP_REPEATCARR(PYRED_PP_NCOEFFCLASSES, PYRED_PP_COEFFARRSIZES)
  return cal;
}

} // namespace pyred
