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
#include <cctype> // isdigit, isalpha, isalnum
#include <random>
#include <sstream>

#include "pyred/coeff_helper.h"
#include "pyred/ppmacros.h"

namespace pyred {

namespace {
std::mt19937_64 random_gen{};
}

// Lower bound for random values for invariants. Must be larger than any root
// of a coefficient polynomial in the system of equations.
uint64_t CoeffHelper::randvar_lbound{2 << 15};
std::vector<std::string> CoeffHelper::s_invariants{};
std::string CoeffHelper::s_settoone{};
uint32_t CoeffHelper::s_coeff_n{1};

void CoeffHelper::random_seed(int s) {
  random_gen.seed(s);
  return;
}

uint64_t CoeffHelper::random_minmax(uint64_t min, uint64_t max) {
  // Generate random numbers until one is in [min,max).
  // May take a long time in small intervals.
  uint64_t r;
  do {
    r = random_gen();
  } while (r >= max || r < min);
  return r;
}

#ifndef VALGRIND
uint64_t CoeffHelper::mod_mul(uint64_t a, uint64_t b, uint64_t m) {
  // result = a*b mod m
  // m must be at most 63 bit.
  // Need floating point type with at least 64 bit mantissa (80 bit float).
  // "long double trick" from https://en.wikipedia.org/wiki/Modular_arithmetic.
  // NOTE: This does not work in valgrind.
  long double x = a;
  uint64_t c = x * b / m;
  int64_t r = static_cast<int64_t>(a * b - c * m) % static_cast<int64_t>(m);
  return r < 0 ? r + m : r;
}
#else
uint64_t CoeffHelper::mod_mul(uint64_t a, uint64_t b, uint64_t m) {
  // result = a*b mod m, safe version to use in valgrind.
  // m must be at most 63 bit.
  uint64_t d = 0;
  uint64_t mp2 = m >> 1;
  for (int i = 0; i != 64; ++i) {
    d = (d > mp2) ? (d << 1) - m : d << 1;
    if (a & 0x8000000000000000ULL) d += b;
    if (d > m) d -= m;
    a <<= 1;
  }
  return d;
}
#endif

uint64_t CoeffHelper::mod_pow(uint64_t base, uint64_t exp, uint64_t modulus) {
  // result = base^exp mod modulus
  // modulus must be at most 63 bit
  uint64_t res = 1;
  while (exp) {
    if (exp & 1) res = mod_mul(res, base, modulus);
    base = mod_mul(base, base, modulus);
    exp >>= 1;
  }
  return res;
}

// Extended Euclidian algorithm to calculate the multiplicative inverse.
// mod_inv(a,m) solves a*t = 1 mod m for t.
// (= mod_pow(a, m-2, m))
uint64_t CoeffHelper::mod_inv(uint64_t a, uint64_t m) {
  int64_t t{0};
  int64_t newt{1};
  int64_t tmpt;
  uint64_t r{m};
  uint64_t newr{a};
  uint64_t tmpr;
  uint64_t q;
  while (newr) {
    q = r / newr;
    tmpt = t;
    t = newt;
    newt = tmpt - q * newt;
    tmpr = r;
    r = newr;
    newr = tmpr - q * newr;
  }
  // if (r > 1) throw init_error("mod_inv: not invertible");
  return t < 0 ? t + m : t;
}

bool CoeffHelper::is_valid_varname(const std::string& s) {
  if (s.empty()) {
    return false;
  }
  else if (!(std::isalpha(s[0]) || s[0] == '_')) {
    return false;
  }
  for (const auto ch : s) {
    if (!(std::isalnum(ch) || ch == '_')) {
      // ok to check the first char again
      return false;
    }
  }
  return true;
}

uint64_t CoeffHelper::parse_longint(const std::string& s, uint64_t prime) {
  // Parse a long integer, passed as a string, take the modulus wrt. prime
  // and return it. The string is split into chunks of at most 18 digits
  // which are then put together via modular arithmetic.
  // Return zero if the string is empty.
  //
  // Make sure the input is an unsigned integer without whitespace padding.
  for (const auto ch : s) {
    if (!std::isdigit(ch)) {
      throw parser_error(
          std::string("parse_longint(): invalid number string \"") + s + "\"");
    }
  }
  uint64_t result = 0;
  std::size_t pos = 0;
  std::size_t len = ((s.size() - 1) % 18) + 1;
  while (pos < s.size()) {
    std::string strchunk = s.substr(pos, len);
    pos += len;
    len = 18;
    uint64_t intchunk;
    std::istringstream ss{strchunk};
    ss >> intchunk;
    // result=0 in the first pass or when the string is zero padded
    // on the left so that the first (few) chunks give zero.
    if (result) result = mod_mul(result, 1000000000000000000uLL, prime);
    result += intchunk;
    result %= prime;
  }
  return result;
}

void CoeffHelper::add_invariant(const std::string& inv) {
  if (Config::finished()) {
    throw init_error("Config::add_invariant(): cannot add new invariants "
                     "after config has been finalised.");
  }
  if (std::find(s_invariants.begin(), s_invariants.end(), inv) ==
      s_invariants.end()) {
    s_invariants.push_back(inv);
  }
}

void CoeffHelper::settoone(const std::string& set1) {
  if (Config::finished()) {
    throw init_error("Config::settoone(): cannot add set 'settoone' "
                     "after config has been finalised.");
  }
  s_settoone = set1;
}

std::string CoeffHelper::settoone() { return s_settoone; }

std::vector<std::string> CoeffHelper::invariants() { return s_invariants; }

void CoeffHelper::coeff_n(uint32_t n) {
  if (Config::finished()) {
    throw init_error("Config::coeff_n(): cannot be changed "
                     "after config has been finalised.");
  }
  s_coeff_n = n;
}

uint32_t CoeffHelper::coeff_n() { return s_coeff_n; }

const std::vector<uint64_t>& CoeffHelper::primes() {
  // the 100 largest 63-bit primes
  static const std::vector<uint64_t> s_primes = {
      9223372036854775783uLL, 9223372036854775643uLL, 9223372036854775549uLL,
      9223372036854775507uLL, 9223372036854775433uLL, 9223372036854775421uLL,
      9223372036854775417uLL, 9223372036854775399uLL, 9223372036854775351uLL,
      9223372036854775337uLL, 9223372036854775291uLL, 9223372036854775279uLL,
      9223372036854775259uLL, 9223372036854775181uLL, 9223372036854775159uLL,
      9223372036854775139uLL, 9223372036854775097uLL, 9223372036854775073uLL,
      9223372036854775057uLL, 9223372036854774959uLL, 9223372036854774937uLL,
      9223372036854774917uLL, 9223372036854774893uLL, 9223372036854774797uLL,
      9223372036854774739uLL, 9223372036854774713uLL, 9223372036854774679uLL,
      9223372036854774629uLL, 9223372036854774587uLL, 9223372036854774571uLL,
      9223372036854774559uLL, 9223372036854774511uLL, 9223372036854774509uLL,
      9223372036854774499uLL, 9223372036854774451uLL, 9223372036854774413uLL,
      9223372036854774341uLL, 9223372036854774319uLL, 9223372036854774307uLL,
      9223372036854774277uLL, 9223372036854774257uLL, 9223372036854774247uLL,
      9223372036854774233uLL, 9223372036854774199uLL, 9223372036854774179uLL,
      9223372036854774173uLL, 9223372036854774053uLL, 9223372036854773999uLL,
      9223372036854773977uLL, 9223372036854773953uLL, 9223372036854773899uLL,
      9223372036854773867uLL, 9223372036854773783uLL, 9223372036854773639uLL,
      9223372036854773561uLL, 9223372036854773557uLL, 9223372036854773519uLL,
      9223372036854773507uLL, 9223372036854773489uLL, 9223372036854773477uLL,
      9223372036854773443uLL, 9223372036854773429uLL, 9223372036854773407uLL,
      9223372036854773353uLL, 9223372036854773293uLL, 9223372036854773173uLL,
      9223372036854773069uLL, 9223372036854773047uLL, 9223372036854772961uLL,
      9223372036854772957uLL, 9223372036854772949uLL, 9223372036854772903uLL,
      9223372036854772847uLL, 9223372036854772801uLL, 9223372036854772733uLL,
      9223372036854772681uLL, 9223372036854772547uLL, 9223372036854772469uLL,
      9223372036854772429uLL, 9223372036854772367uLL, 9223372036854772289uLL,
      9223372036854772241uLL, 9223372036854772169uLL, 9223372036854772141uLL,
      9223372036854772061uLL, 9223372036854772051uLL, 9223372036854772039uLL,
      9223372036854771989uLL, 9223372036854771977uLL, 9223372036854771973uLL,
      9223372036854771953uLL, 9223372036854771869uLL, 9223372036854771841uLL,
      9223372036854771833uLL, 9223372036854771797uLL, 9223372036854771749uLL,
      9223372036854771737uLL, 9223372036854771727uLL, 9223372036854771703uLL,
      9223372036854771689uLL};
  return s_primes;
}

} // namespace pyred
