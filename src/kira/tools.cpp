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

#include "kira/tools.h"

#include <sys/stat.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>

#include "pyred/integrals.h"

using namespace std;
using namespace GiNaC;

bool file_exists(const char* filename) {
  ifstream ifile(filename);
  return static_cast<bool>(ifile);
}

int mkpath(char* file_path, mode_t mode) {
  assert(file_path && *file_path);
  char* p;
  for (p = strchr(file_path + 1, '/'); p; p = strchr(p + 1, '/')) {
    *p = '\0';
    if (mkdir(file_path, mode) == -1) {
      if (errno != EEXIST) {
        *p = '/';
        return -1;
      }
    }
    *p = '/';
  }
  return 0;
}

std::size_t binomial_coeff(std::size_t n, std::size_t k) {
  std::size_t res = 1;
  if (n < k) return 0;
  if (k > n - k) k = n - k;
  for (std::size_t i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  return res;
}

unsigned powerINT(unsigned base, unsigned degree) {
  unsigned result = 1;
  unsigned term = base;
  while (degree) {
    if (degree & 1) result *= term;
    term *= term;
    degree = degree >> 1;
  }
  return result;
}

void load_bar(unsigned int x, unsigned int n, unsigned int w = 50,
              unsigned int steps = 100) {
  if ((x != n) && (x % (n / steps + 1) != 0)) return;

  float ratio = x / (float)n;
  unsigned int c = ratio * w;

  cout << std::fixed << std::setprecision(1) << setw(5) << (double)(ratio * 100)
       << "% [";
  for (unsigned int i = 0; i < c; i++)
    cout << "=";
  cout << "]\r" << flush;
}

const GiNaC::possymbol& get_symbol(const string& s) {
  static map<string, GiNaC::possymbol> directory;
  map<string, GiNaC::possymbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, GiNaC::possymbol(s))).first->second;
}

void generate_symbols(possymbol var[], string str1, int bound) {
  string str2;
  for (int i = 0; i < bound; i++) {
    str2 = str1 + something_string(i);
    var[i] = get_symbol(str2);
  }
}

void get_properties(pyred::Weight id,
                    std::tuple<std::string, unsigned, unsigned>& integral) {
  auto iglback = pyred::Integral(id);
  auto property = iglback.properties(id);

//   string strIndices;
//   for(size_t itt = 0; itt < iglback.m_powers.size(); itt++){
//     strIndices += to_string(iglback.m_powers[itt]);
//     if( itt != iglback.m_powers.size() -1 )
//       strIndices += ",";
//   }
  get<0>(integral) = iglback.to_string(pyred::Integral::StringFormat::indices);
  get<1>(integral) = property.topology;
  get<2>(integral) = property.sector;
}

std::string masters_list_str(const std::vector<pyred::Weight> &ids) {
  // Output format:
  // T[1,1,1,1,1]  # 31
  // T[1,-1,0,0,0] #  1
  // ...
  std::vector<std::pair<std::string,std::string>> iglstr_secstr_vec;
  for (const auto &iid: ids) {
    iglstr_secstr_vec.push_back({
      pyred::Integral(iid).to_string(),
      std::to_string(pyred::Integral::properties(iid).sector)});
  }
  std::size_t iglstr_len = 0u;
  std::size_t secstr_len = 0u;
  for (const auto &iglstr_secstr: iglstr_secstr_vec) {
    iglstr_len = std::max(iglstr_len, iglstr_secstr.first.size());
    secstr_len = std::max(secstr_len, iglstr_secstr.second.size());
  }
  std::string out = "";
  for (const auto &iglstr_secstr: iglstr_secstr_vec) {
    out += iglstr_secstr.first;
    for (std::size_t co = 0u; co < iglstr_len-iglstr_secstr.first.size(); ++co)
      out.push_back(' ');
    out += " # ";
    for (std::size_t co = 0u; co < secstr_len-iglstr_secstr.second.size(); ++co)
      out.push_back(' ');
    out += iglstr_secstr.second;
    out.push_back('\n');
  }
  return out;
}

std::string masters_list_str_trivial(const std::vector<pyred::Integral> &ints) {
  // TODO sector cannot easily be extracted for trivial sectors
  std::vector<std::string> iglstr_vec;
  for (const auto &inte: ints) {
    iglstr_vec.emplace_back(inte.to_string());
  }
  std::string out = "";
  for (const auto &iglstr_secstr: iglstr_vec) {
    out += iglstr_secstr;
    out.push_back('\n');
  }
  return out;
}

double Clock::eval_time() {
  timeval tend;
  gettimeofday(&tend, 0);
  return (tend.tv_sec - tstart.tv_sec +
          1.e-6 * (tend.tv_usec - tstart.tv_usec));
}

size_t get_next_num(vector<size_t>& v){
  size_t n = v.size();
  srand(time(NULL));
  size_t index = rand() % n;
  size_t num = v[index];
  swap(v[index], v[n - 1]);
  v.pop_back();
  return num;
}
 
vector<size_t> generate_random(size_t n){
  vector<size_t> v(n), v2(n);
  for (size_t i = 0; i < n; i++)
    v[i] = i;
  size_t count=0;
  while (v.size()) {
    v2[count++]=get_next_num(v);
  }
  return v2;
}

