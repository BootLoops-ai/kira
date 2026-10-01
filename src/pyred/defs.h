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

#ifndef PYRED_DEFS_H
#define PYRED_DEFS_H

#include <cstdint>
#include <exception>
#include <fstream>
#include <functional>
#include <istream>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pyred/config.h"

namespace pyred {
  inline std::size_t hash_combine(std::size_t h1, std::size_t h2) {
    h1 ^= (h2 + 0x9e3779b97f4a7c15uLL + (h1 << 6) + (h1 >> 2));
    // Alternative implementation from Boost:
    // const std::size_t fac = 0xc6a4a7935bd1e995uLL;
    // h2 *= fac;
    // h2 ^= h2 >> 47;
    // h2 *= fac;
    // h1 ^= h2;
    // h1 *= fac;
    // h1 += 0xe6546b64;
    return h1;
  }
} // namespace pyred

#include "pyred/weight.h"


// hash functions

namespace std {
  template <typename T>
  struct hash<vector<T>> {
    size_t operator()(const vector<T>& vec) const {
      size_t h{0u};
      for (const auto& elem : vec) {
        h = pyred::hash_combine(h, hash<T>{}(elem));
      }
      return h;
    }
  };
  template <typename T1, typename T2>
  struct hash<pair<T1, T2>> {
    size_t operator()(const pair<T1, T2>& p) const {
      size_t h{0u};
      h = pyred::hash_combine(h, hash<T1>{}(p.first));
      h = pyred::hash_combine(h, hash<T2>{}(p.second));
      return h;
    }
  };
} // namespace std


namespace pyred {

using std::uint32_t;
using std::uint64_t;

class init_error : public std::exception {
private:
  std::string msg;
public:
  inline init_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

class parser_error : public std::exception {
private:
  std::string msg;
public:
  inline parser_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
  std::string message() const { return msg; }
};

class zero_modular_division : public std::exception {
public:
  inline zero_modular_division() {}
  virtual inline const char* what() const noexcept {
    return "Modular arithmetic lead to a spurious zero division.";
  }
};

class input_error : public std::exception {
private:
  std::string msg;
public:
  inline input_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

class runtime_error : public std::exception {
private:
  std::string msg;

public:
  inline runtime_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

class uf_error : public std::exception {
private:
  std::string msg;
public:
  inline uf_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

class permutation_error : public std::exception {
private:
  std::string msg;
public:
  inline permutation_error(const std::string& s) : msg(s) {}
  virtual inline const char* what() const noexcept { return msg.c_str(); }
};

// little helper functions

inline int string_to_int(const std::string& s) {
  // Error checking conversion of a string into an integer.
  // Ignores leading and trailing spaces.
  std::size_t pos;
  int i;
  try {
    i = std::stoi(s, &pos);
  }
  catch (const std::invalid_argument&) {
    throw input_error(std::string("Invalid integer string '") + s + "'.");
  }
  while (pos != s.size()) {
    if (s[pos] != ' ') {
      throw input_error(std::string("Invalid integer string '") + s + "'.");
    }
    ++pos;
  }
  return i;
}

inline std::istream& getline_normalized(std::istream& in, std::string& out, char delim = '\n') {
  std::getline(in, out, delim);
  if (!out.empty() && out.back() == '\r')
    out.pop_back();
  return in;
}

// stream insertion operators

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream&, const std::pair<T1, T2>&);

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
  out << "[";
  for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
    if (it != vec.cbegin()) out << ",";
    out << *it;
  }
  out << "]";
  return out;
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& p) {
  out << "<" << p.first << "," << p.second << ">";
  return out;
}

class Integral;
using pow_type = int;
using is_equation = std::vector<std::pair<Integral, std::string>>;
template <typename Coeff>
using wc_pair = std::pair<Weight, Coeff>;
template <typename Coeff>
using wc_equation = std::vector<wc_pair<Coeff>>;
using eqdata = wc_equation<std::string>;

static_assert(!(sizeof(Weight) % sizeof(pow_type)),
              "sizeof(Weight) must be a multiple of sizeof(pow_type)");
static_assert(!(sizeof(Weight) % sizeof(uint64_t)),
              "sizeof(Weight) must be a multiple of sizeof(uint64_t)");

template <typename KEY, typename VALUE>
class Cache {
public:
  Cache() : m_cache{}, m_treatval{nullptr} {}
  bool insert(const KEY& key, const VALUE& val) {
    auto ins = m_cache.insert({key, val});
    if (m_treatval) {
      if (ins.second) {
        // New key, i.e. insertion took place.
        ins.first->second = m_treatval(key, val);
      }
    }
    return ins.second;
  }
  std::pair<VALUE, bool> lookup(const KEY& key) const {
    auto it = m_cache.find(key);
    if (it == m_cache.cend()) {
      return {VALUE{}, false};
    }
    return {it->second, true};
  }
  std::size_t size() const { return m_cache.size(); }
  bool empty() const { return m_cache.empty(); }
  void set_if_empty(const Cache<KEY, VALUE>& init) {
    if (m_cache.empty() && !init.empty()) {
      m_cache = init.get();
    }
  }
  void clear() { m_cache.clear(); }
  std::unordered_map<KEY, VALUE> get() const { return m_cache; }
  std::unordered_map<KEY, VALUE>& get_unsafe() { return m_cache; }
  const std::unordered_map<KEY, VALUE>& get_unsafe() const { return m_cache; }
  // A workaround for a compiler bug in gcc: call this manually to
  // call constructors when a cache that depends on a template parameter
  // is thread local. Works for the constructors, but the destructors still
  // won't be called. Therefore the code has been reworked to not use
  // thread local templated caches.
  //   void ensure_init() {
  //     thread_local static bool first{true};
  //     if (first) {
  //       m_cache = std::unordered_map<KEY,VALUE>{};
  //       m_treatval = nullptr;
  //       first = false;
  //     }
  //   }
  void treatval_func(
      const std::function<VALUE(const KEY&, const VALUE&)>& treatval) {
    m_treatval = treatval;
  }
  bool treatval_set() const { return bool(m_treatval); }

private:
  std::unordered_map<KEY, VALUE> m_cache;
  std::function<VALUE(const KEY&, const VALUE&)> m_treatval;
};

template <typename KEY, typename VALUE>
class LockedCache {
public:
  bool insert(const KEY& key, const VALUE& val) {
    std::lock_guard<std::mutex> lock(m_mtx);
    return m_cache.insert({key, val}).second;
  }
  std::pair<VALUE, bool> lookup(const KEY& key) const {
    std::lock_guard<std::mutex> lock(m_mtx);
    auto it = m_cache.find(key);
    if (it == m_cache.cend()) {
      return {VALUE{}, false};
    }
    return {it->second, true};
  }
  std::size_t size() const {
    std::lock_guard<std::mutex> lock(m_mtx);
    return m_cache.size();
  }
  void clear() {
    std::lock_guard<std::mutex> lock(m_mtx);
    m_cache.clear();
  }
  std::unordered_map<KEY, VALUE> get() const {
    std::lock_guard<std::mutex> lock(m_mtx);
    return m_cache;
  }
  std::unordered_map<KEY, VALUE>& get_unsafe() { return m_cache; }

private:
  std::unordered_map<KEY, VALUE> m_cache;
  mutable std::mutex m_mtx;
};

/*************************************************************
 * Logger and LockedLogger classes: defined in interface.cpp *
 *************************************************************/

enum class FileOpenMode { create, truncate, append };

class LockedLogger;

class Logger {
  // Logger log{cout_verbosity=1, file_verbosity=-2};
  //   if file_verbosity=-2 (default), set file_verbosity=cout_verbosity.
  // log(cout_level=1, file_level=-1) << out;
  //   if file_level=-1 (default), set file_level=cout_level, then
  //   write 'out' to cout if cout_level <= cout_verbosity,
  //   write 'out' to file if file_level <= file_verbosity.
  //   Locks a mutex to wite to cout exclusively as long as the LockedLogger
  //   object exists (don't keep it alive during long calculation!).
  //   TODO: buffer all data and only lock the mutex and write out all data
  //         when the LockedLogger is destructed.
  // attach_logfile(filename, mode)
  //   mode: FileOpenMode::create, FileOpenMode::truncate, FileOpenMode::append.
  // detach_logfile().
  friend class LockedLogger;

private:
  // m_verbosity.first: verbosity for cout,
  // m_verbosity.second: verbosity for file output.
  std::pair<int, int> m_verbosity;
  std::ofstream m_logfilestream;
  // The same thread may lock the mutex multiple times. This could happen
  // if a reference to a LockedLogger is kept while a new LockedLogger
  // is created. Locking a std::mutex multiple times from the same thread
  // is undefined behaviour (and does not necessarily block).
  std::recursive_mutex m_mtx;

public:
  Logger(int cout_verb = 1, int file_verb = -2);
  std::pair<int, int> verbosity(int cout_verb = 1, int file_verb = -2);
  void attach_logfile(const std::string&,
                      const FileOpenMode = FileOpenMode::truncate);
  void detach_logfile();
  LockedLogger operator()(int cout_lev = 1, int file_lev = -1);
};

class LockedLogger {
  friend class Logger;

private:
  std::pair<int, int> m_verbosity;
  std::pair<int, int> m_level;
  std::ofstream& m_logfilestream;
  std::unique_lock<std::recursive_mutex> m_lck;
  LockedLogger(Logger&, int cout_lev, int file_lev);

public:
  template <typename T>
  LockedLogger& operator<<(T&& out) {
    if (m_level.first <= m_verbosity.first) std::cout << std::forward<T>(out);
    if (m_level.second <= m_verbosity.second && m_logfilestream.is_open()) {
      m_logfilestream << std::forward<T>(out);
    }
    return *this;
  }
  // Specialisation for I/O manipulators std::endl and std::flush.
  LockedLogger& operator<<(std::ostream& (*)(std::ostream&));
};

/******************************************
 * Config class: defined in interface.cpp *
 *****************************************/

class Config {
public:
  static void coeff_cls(int);
  static int coeff_cls() { return s_coeff_cls_val; }
  static void backward(bool);
  static bool backward() { return s_backward_val; }
  static std::pair<int, int> verbosity(int = 1, int = -2);
  static void attach_logfile(const std::string& fn,
                             const FileOpenMode mode = FileOpenMode::truncate) {
    s_log_val.attach_logfile(fn, mode);
  }
  static void detach_logfile() { s_log_val.detach_logfile(); }
  static LockedLogger log(int cout_lev = 1, int file_lev = -1) {
    return s_log_val(cout_lev, file_lev);
  }
  static void parallel(int);
  static int parallel() { return s_parallel; }
  static void auto_symseed(const bool);
  static bool auto_symseed() { return s_auto_symseed; }
  static void symlimits(const int, const int);
  static std::pair<int, int> symlimits() { return s_symlimits; }
  static int lookahead() { return s_lookahead; }
  static void lookahead(int);
  static int insertion_tracer() { return s_insertion_tracer; }
  static void insertion_tracer(int);
  static void database_file(const std::string& fn, bool overwrite = false) {
    s_database_file = {fn, overwrite};
  }
  static std::pair<std::string, bool> database_file() {
    return s_database_file;
  }
  static std::pair<std::string, std::string> parse(const std::string& = "");
  static std::pair<std::string, std::string> parse(
      const std::vector<std::string>&);
  static void finish();
  static bool finished() { return s_finished; }
  static bool johanntrick() { return s_johanntrick; }
  static void johanntrick(bool jt) { s_johanntrick = jt; }

private:
  static Logger s_log_val;
  static int s_coeff_cls_val; // 1,2,... as defined in ppmacros.h
  static bool s_backward_val;
  static int s_parallel; // 0 means automatic
  static bool s_auto_symseed;
  // Limits for automatic symmetry seed generation:
  // use seeds up to s_symlimits=(maxdots,maxsps).
  // Do not symmetrise if maxdots<0 and/or maxsps<0
  static std::pair<int, int> s_symlimits;
  // lookahead in OTF generator/solver:
  // -2: Generate sectors in sector order;
  //     sort equations within each sector seed,
  //     first sorted IBP equations, then sorted symmetry equations.
  //     Solve the next sector as soon as it has been generated.
  // -1: Generate sectors in sector order;
  //     sort equations within each sector seed.
  //     Solve the next sector as soon as it has been generated.
  //  0: Schedule the sectors to be solved right after their bunch
  //     has been generated.
  //  1: Schedule the sectors to be solved after the next bunch
  //     has been generated.
  //  2: Solve after all sectors have been generated.
  static int s_lookahead;
  // Insertion tracer:
  //   0: off,
  //   1: in memory (default),
  //   2: in SQLite database,
  //   3: in Kyoto Cabinet database,
  //   4: Kyoto Cabinet if available, otherwise SQLite.
  // Only effective if compiled with PYRED_EXTERNAL_TRACER
  static int s_insertion_tracer;
  // Temporary file without file extension used for the insertion database.
  // bool: if true, overwrite existing file.
  static std::pair<std::string, bool> s_database_file;
  static bool s_finished;
  static bool s_johanntrick;
};

// Defined in interface.cpp.
class FileSystem {
public:
  enum class FileType {regular, directory, other};
  static FileType get_filetype(const std::string &/*fname*/);
  static std::vector<std::string>
    get_filenames(const std::string &/*file_or_dir*/,
                  const std::vector<std::string> &/*fileexts*/ = {});
};

} // namespace pyred

#endif
