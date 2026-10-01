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

#ifndef TOOLS_H_
#define TOOLS_H_
#include <sys/stat.h>
#include <sys/time.h>

#include <fstream>
#include <sstream>
#include <streambuf>
#include <string>
#include <type_traits>

#include "pyred/defs.h"
#include "ginac/ginac.h"

void get_properties(pyred::Weight id,
                    std::tuple<std::string, unsigned, unsigned> &integral);

std::string masters_list_str(const std::vector<pyred::Weight> &);

std::string masters_list_str_trivial(const std::vector<pyred::Integral> &);


class Clock {
public:
  Clock() { gettimeofday(&tstart, 0); };
  double eval_time();

private:
  /*static*/ timeval tstart;
};

bool file_exists(const char *filename);
void load_bar(unsigned int x, unsigned int n, unsigned int w,
              unsigned int steps);
const GiNaC::possymbol &get_symbol(const std::string &s);
void generate_symbols(GiNaC::possymbol unknown[], std::string str1, int bound);
std::size_t binomial_coeff(std::size_t n, std::size_t k);

template <class SomeType>
GiNaC::ex fs(GiNaC::ex &koef, const SomeType &var) {
  koef = koef.expand().subs(var, GiNaC::subs_options::algebraic);
  return koef;
}

template <typename T>
std::string something_string(T &t) {
  std::ostringstream ss;
  ss << t;
  return ss.str();
}

template <typename T>
std::string int_string(T t) {
  std::ostringstream ss;
  ss << t;
  return ss.str();
}

template <typename T>
int something_int(T t) {
  int token;
  std::stringstream temp;
  temp << t;
  temp >> token;
  return token;
}

class Loginfo {
public:
  static inline Loginfo &instance() {
    static Loginfo loginfo;
    return loginfo;
  }
  void set_logfile(const std::string &s) {
    coutBuf = std::cout.rdbuf(); // save old buf
    outfile.open(s.c_str());
  };
  void silent_cout_buf() {
    silentBuf = std::cout.rdbuf(); // save old buf
    silentfile.open("/dev/null");
    std::cout.rdbuf(silentfile.rdbuf());
    coutBuf = std::cout.rdbuf();
  };
  void restore_silence() {
    std::cout.rdbuf(silentBuf);
    coutBuf = std::cout.rdbuf();
  };
  void set_level(int level_) { level = level_; };

  template <typename T>
  inline Loginfo &operator<<(const T &out) {
    using pyred::operator<<;
    switch (level) {
      case 0:
        std::cout.rdbuf(coutBuf); // restore original std::cout
        std::cout << out << std::flush;
        break;
      case 1:
        std::cout.rdbuf(outfile.rdbuf()); // redirect std::out outfile
        std::cout << out << std::flush;
        std::cout.rdbuf(coutBuf); // restore original std::cout
        std::cout << out << std::flush;
        break;
      case 2:
        std::cout.rdbuf(outfile.rdbuf()); // redirect std::out outfile
        std::cout << out << std::flush;
        std::cout.rdbuf(coutBuf); // restore original std::cout
        break;
      default:
        std::cout << out << std::flush;
    }
    return *this;
  }
  std::streambuf *coutBuf, *silentBuf;
  std::string logfile;
  std::ofstream outfile, silentfile;

private:
  Loginfo() : level(1){};
  int level;
};

unsigned powerINT(unsigned base, unsigned degree);

int mkpath(char *file_path, mode_t mode);

std::size_t get_next_num(std::vector<std::size_t>& v);
 
std::vector<std::size_t> generate_random(std::size_t n);

#endif
