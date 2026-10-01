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

#ifndef KEYVALUEDB_H
#define KEYVALUEDB_H

#include <cstdint>
#include <exception>
#include <string>
#include <vector>

#ifdef PYRED_KCDB
#include <kcpolydb.h>
#endif

#include "pyred/config.h"
#include "defs.h"
#include "sqlite3/sqlite3.h"

namespace keyvaluedb {

using std::uint64_t;
using pyred::Weight;
using pyred::Config;

class database_error : public std::exception {
private:
  std::string msg;

public:
  inline database_error(const std::string &s) : msg(s) {}
  virtual inline const char *what() const noexcept { return msg.c_str(); }
};

// Base class for database interfaces.
class KeyValueDB {
public:
  virtual ~KeyValueDB(){};
  virtual void put(uint64_t, const std::vector<Weight> &) = 0;
  // get() must be non-const, because the Kyoto Cabinet get() is not const.
  virtual std::vector<Weight> get(uint64_t, bool = true /*fatal*/) = 0;
  virtual void remove(uint64_t, bool = true /*fatal*/) = 0;
  virtual std::string filename() const { return ""; }
  static void tmpdir(const std::string &dir) { s_tmpdir = dir; }
  static std::string tmpdir() { return s_tmpdir; }
  static bool file_exists(const std::string &filename);

private:
  static std::string s_tmpdir;
};

class KeyValueDiscard : public KeyValueDB {
public:
  KeyValueDiscard(const std::string & = "", std::size_t = 0,
                  bool /*overwrite*/ = false) {}
  void put(uint64_t, const std::vector<Weight> &) override {}
  std::vector<Weight> get(uint64_t, bool = true) override { return {}; }
  void remove(uint64_t, bool = true /*fatal*/) override {}
};

class KeyValueVector : public KeyValueDB {
public:
  KeyValueVector(const std::string & = "", std::size_t sz = 0,
                 bool /*overwrite*/ = false) {
    m_db.resize(sz);
  }
  void put(uint64_t key, const std::vector<Weight> &val) override {
    if (m_db.size() <= key) {
      auto newsz = m_db.size();
      while (newsz <= key) {
        newsz *= 1.5;
      }
      m_db.reserve(newsz);
      m_db.resize(newsz);
    }
    m_db[key] = val;
  }
  std::vector<Weight> get(uint64_t key, bool fatal = true) override {
    if (m_db.size() <= key) {
      if (fatal) {
        throw database_error(std::string("Error retrieving key ") +
                             std::to_string(key) +
                             " from vector database: overflow");
      }
      else {
        return {};
      }
    }
    return m_db[key];
  }
  void remove(uint64_t key, bool fatal = true) override {
    if (m_db.size() >= key && fatal) {
      throw database_error(std::string("Error deleting key ") +
                           std::to_string(key) +
                           " from vector database: overflow");
    }
    m_db[key].clear();
  }

private:
  std::vector<std::vector<Weight>> m_db;
};

class KeyValueSQLite : public KeyValueDB {
  // SQLite3 wrapper to provide key-value database functionality
  // for uint64_t keys and std::vector<Weight> values.
  // No transactions.
public:
  KeyValueSQLite(const std::string & = "", std::size_t = 0,
                 bool overwrite = false);
  ~KeyValueSQLite() override;
  void put(uint64_t, const std::vector<Weight> &) override;
  std::vector<Weight> get(uint64_t, bool = true /*fatal*/) override;
  void remove(uint64_t key, bool fatal = true) override;
  std::string filename() const override;

private:
  sqlite3 *m_db;
  sqlite3_stmt *m_put_stmt;
  sqlite3_stmt *m_get_stmt;
  sqlite3_stmt *m_rm_stmt;
  void execute(const std::string &);
};

#ifdef PYRED_KCDB
class KeyValueKC : public KeyValueDB {
  // Kyoto Cabinet wrapper to mimic the interface of the SQLite3 wrapper.
public:
  KeyValueKC(const std::string & = ":", std::size_t = 0,
             bool overwrite = false);
  ~KeyValueKC() override;
  void put(uint64_t, const std::vector<Weight> &) override;
  std::vector<Weight> get(uint64_t, bool = true /*fatal*/) override;
  void remove(uint64_t key, bool fatal = true) override;
  std::string filename() const override { return m_filename; };

private:
  kyotocabinet::PolyDB m_db;
  std::string m_filename;
};
#endif // PYRED_KCDB

} // namespace keyvaluedb

#endif // KEYVALUEDB_H
