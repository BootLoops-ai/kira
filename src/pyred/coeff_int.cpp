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
#include <cctype> // isalpha
#include <sstream>

#include "pyred/coeff_helper.h" // mod_mul, mod_inv
#include "pyred/coeff_int.h"
#include "pyred/parser.h"
#include "yaml-cpp/yaml.h"

namespace pyred {

// static members
bool Coeff_int::auto_vars{true};
std::string Coeff_int::set_to_one{};
uint64_t Coeff_int::prime = CoeffHelper::primes().front();
std::vector<std::pair<std::string, uint64_t>> Coeff_int::randvars{};
std::mutex Coeff_int::s_mtx{};
std::string Coeff_int::auxName;

void Coeff_int::set_prime(uint64_t tempPrime){
  Coeff_int::prime = tempPrime;
}
uint64_t Coeff_int::get_prime(){
  return Coeff_int::prime;
}

void Coeff_int::init() {
  auto n = CoeffHelper::s_coeff_n;
  if (n < 1 || n > CoeffHelper::primes().size()) {
    std::ostringstream ss;
    ss << "Invalid choice of prime: " << n << std::endl
       << "Must be 1<=n<=" << CoeffHelper::primes().size() << std::endl;
    throw input_error(ss.str());
  }
  std::cout << "***********************************************************" << std::endl;
  std::cout << "Prepare numerical values for the generation of the system" << std::endl;
  std::cout << "***********************************************************" << std::endl;
  //std::cout << "The system is solved modulo prime: " << CoeffHelper::primes().at(n - 1) << std::endl;
  //std::cout << "which is the nth " << n << " prime in the list of hard coded primes in Kira." << std::endl;
  //std::cout << "To see the list of hard coded primes in Kira type kira --help." << std::endl;
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
  prime = CoeffHelper::primes().at(n - 1);
  randvars.clear();
  randvars.reserve(CoeffHelper::s_invariants.size());
  
  // import auxiliary.yaml which is automaticaly generated
  std::string auxiliary_file = auxName;
  YAML::Node auxiliary;
  std::ifstream ifile(auxiliary_file);
  if(static_cast<bool>(ifile)){
    try {
      auxiliary = YAML::LoadFile(auxiliary_file);
    }
    catch (const YAML::BadFile& e) {
      std::cerr << "Error reading file \"" << auxiliary_file << "\""
                << std::endl;
      throw;
    }
  }

  if(auxiliary["auxiliary"]){
    auto auxiliary_node = auxiliary["auxiliary"];
    std::cout << "set numerics" << std::endl;
    for (const auto& auxdef : auxiliary_node) {
      if(auxdef["user_prime"]){
        //uint64_t tmpNumber;
        std::stringstream temp;
        temp << auxdef["user_prime"].as<std::string>();
        temp >> prime;
        //prime = CoeffHelper::primes().at((auxdef["user_prime"].as<int>()) - 1);
        std::cout << "Prime used to aid the system generation: " << Coeff_int::prime << std::endl;
      }
      auto tmp_node = auxdef["user_numerics"];
      if (tmp_node) {
	      for(auto ix: tmp_node){
          std::cout << "User value: " << ix[1].as<std::string>() << std::endl;
          std::cout << "Is converted to a uint64_t number modulo prime "<< prime <<": " << pyred::parse_coeff<Coeff_int>(ix[1].as<std::string>()) << std::endl;
          //std::cout << ix[1].as<uint64_t>() << std::endl;
          uint64_t tmpNumber;
          std::stringstream temp;
          temp << pyred::parse_coeff<Coeff_int>(ix[1].as<std::string>());
          temp >> tmpNumber;
	        randvars.push_back({ix[0].as<std::string>(), tmpNumber});
          //randvars.push_back({ix[0].as<std::string>(), ix[1].as<uint64_t>()});
	      }
      }
    }
  }
  
//  std::cout << randvars[0].first << " " << randvars[0].second << std::endl;
  for (const auto& var : CoeffHelper::s_invariants) {
    if (var == CoeffHelper::s_settoone) {
      randvars.push_back({var, 1});
    }
    else {
      int flagC=0;
      for (const auto& vc : randvars) {
        if (vc.first == var) {
          flagC=1;
        }
      }
      if(flagC==1)
        continue;
//      std::cout << var << " " << CoeffHelper::random_minmax(
//                                   CoeffHelper::randvar_lbound, prime) << std::endl;
      randvars.push_back({var, CoeffHelper::random_minmax(
                                   CoeffHelper::randvar_lbound, prime)});
    }
  }
  return;
}

Coeff_int::Coeff_int() {}
Coeff_int::Coeff_int(const Coeff_int& a) : c(a.c) {}

Coeff_int::Coeff_int(const std::string& s) {
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
  auto success = static_cast<bool>(ss >> c);
  if (!(success && ss.rdbuf()->in_avail() == 0)) {
    if (s.empty()) {
      // (ss >> c) fails if ss is empty
      throw parser_error("Coeff_int(): empty argument");
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
        c = 1;
      }
      else {
        c = CoeffHelper::random_minmax(CoeffHelper::randvar_lbound, prime);
      }
      std::lock_guard<std::mutex> lock(s_mtx);
      CoeffHelper::s_invariants.push_back(s);
      randvars.push_back({s, c});
    }
    else {
      c = CoeffHelper::parse_longint(s, prime);
    }
  }
  else if (c >= prime) {
    // special case: the parsed value fits into c, but is >= prime
    c %= prime;
  }
}

Coeff_int& Coeff_int::operator+=(const Coeff_int& a) {
  c += a.c;
  if (c >= prime) c -= prime;
  return *this;
}
Coeff_int& Coeff_int::operator-=(const Coeff_int& a) {
  if (a.c > c) c += prime;
  c -= a.c;
  return *this;
}
Coeff_int& Coeff_int::operator*=(const Coeff_int& a) {
  c = CoeffHelper::mod_mul(c, a.c, prime);
  return *this;
}
Coeff_int& Coeff_int::operator/=(const Coeff_int& a) {
  c = CoeffHelper::mod_mul(c, CoeffHelper::mod_inv(a.c, prime), prime);
  return *this;
}
Coeff_int pow(const Coeff_int& b, const Coeff_int& a) {
  Coeff_int result;
  if (a.c == 2u) {
    // Fast-track the most common case
    result.c = CoeffHelper::mod_mul(b.c, b.c, Coeff_int::prime);
  }
  else {
    uint64_t exp;
    uint64_t base;
    if (a.c < (Coeff_int::prime >> 1)) {
      // treat as positive exponent
      exp = a.c;
      base = b.c;
    }
    else {
      // treat as negative exponent
      exp = Coeff_int::prime - a.c;
      base = CoeffHelper::mod_inv(b.c, Coeff_int::prime); // =1/b.c
    }
    result.c = CoeffHelper::mod_pow(base, exp, Coeff_int::prime);
  }
  return result;
}

Coeff_int operator-(const Coeff_int& a) {
  if (a.c == 0) {
    return Coeff_int(a);
  }
  else {
    return Coeff_int(Coeff_int::prime - a.c);
  }
}
Coeff_int operator+(const Coeff_int& a, const Coeff_int& b) {
  auto sum = a.c + b.c;
  if (sum >= Coeff_int::prime) sum -= Coeff_int::prime;
  return Coeff_int(sum);
}
Coeff_int operator-(const Coeff_int& a, const Coeff_int& b) {
  auto diff = a.c;
  if (b.c > diff) diff += Coeff_int::prime;
  diff -= b.c;
  return Coeff_int(diff);
}
Coeff_int operator*(const Coeff_int& a, const Coeff_int& b) {
  return Coeff_int(CoeffHelper::mod_mul(a.c, b.c, Coeff_int::prime));
}
Coeff_int operator/(const Coeff_int& a, const Coeff_int& b) {
  return Coeff_int(CoeffHelper::mod_mul(
      a.c, CoeffHelper::mod_inv(b.c, Coeff_int::prime), Coeff_int::prime));
}

bool operator==(const Coeff_int& a, const Coeff_int& b) { return (a.c == b.c); }

bool operator!=(const Coeff_int& a, const Coeff_int& b) { return (a.c != b.c); }

std::ostream& operator<<(std::ostream& out, const Coeff_int& a) {
  out << a.c;
  return out;
}

} // namespace pyred
