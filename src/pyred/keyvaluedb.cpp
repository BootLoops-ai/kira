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

#include <cstdio>
#include <iostream>

#include "sys/stat.h"

#include "pyred/keyvaluedb.h"

namespace keyvaluedb {

bool KeyValueDB::file_exists(const std::string &filename) {
  struct stat st;
  return (stat(filename.c_str(), &st) == 0);
}

/******************
 * KeyValueSQLite *
 ******************/

std::string KeyValueDB::s_tmpdir{};

KeyValueSQLite::KeyValueSQLite(const std::string &filename, std::size_t,
                               bool overwrite) {
  if (file_exists(filename)) {
    if (overwrite) {
      std::remove(filename.c_str());
    }
    else {
      throw database_error(std::string("Database file '") + filename +
                           "' already exists.");
    }
  }
  auto ret = sqlite3_open(filename.c_str(), &m_db);
  if (ret != SQLITE_OK) {
    throw database_error(std::string("Error opening database '") + filename +
                         "': " + sqlite3_errmsg(m_db) + ".");
  }
  // TEMP databases always use synchronous=off.
  execute("PRAGMA synchronous = off;");
  execute("PRAGMA journal_mode = off;");
  // TEMP databases always use locking_mode=exclusive.
  execute("PRAGMA locking_mode = exclusive;");
  if (!tmpdir().empty()) {
    execute(std::string("PRAGMA temp_store_directory = \"") + tmpdir() + "\"");
  }
  // BLOB KEY
  // execute("CREATE TABLE kv(k BLOB PRIMARY KEY, v BLOB);");
  execute("CREATE TABLE kv(k INTEGER PRIMARY KEY, v BLOB);");
  std::string putcmd = "REPLACE INTO kv VALUES(?, ?)";
  ret = sqlite3_prepare_v2(m_db, putcmd.c_str(), -1, &m_put_stmt, nullptr);
  if (ret != SQLITE_OK) {
    throw database_error(std::string("Database '") + filename +
                         "': Error compiling put statement.");
  }
  std::string getcmd = "SELECT v FROM kv WHERE k = ?";
  ret = sqlite3_prepare_v2(m_db, getcmd.c_str(), -1, &m_get_stmt, nullptr);
  if (ret != SQLITE_OK) {
    throw database_error(std::string("Database '") + filename +
                         "': Error compiling get statement.");
  }
  std::string rmcmd = "DELETE FROM kv WHERE k = ?";
  ret = sqlite3_prepare_v2(m_db, rmcmd.c_str(), -1, &m_rm_stmt, nullptr);
  if (ret != SQLITE_OK) {
    throw database_error(std::string("Database '") + filename +
                         "': Error compiling remove statement.");
  }
}

KeyValueSQLite::~KeyValueSQLite() {
  auto fn = filename();
  sqlite3_finalize(m_put_stmt);
  sqlite3_finalize(m_get_stmt);
  sqlite3_finalize(m_rm_stmt);
  // Close even if open failed.
  auto ret = sqlite3_close(m_db);
  if (ret != SQLITE_OK) {
    Config::log(1) << "Warning: error closing database '" << fn << "'."
                   << std::endl;
  }
  if (!fn.empty() && fn.front() != ':') {
    std::remove(fn.c_str());
  }
}

void KeyValueSQLite::put(uint64_t key, const std::vector<Weight> &val) {
  // BLOB KEY
  // sqlite3_bind_blob(m_put_stmt, 1, &key, static_cast<int>(sizeof key),
  //                   SQLITE_STATIC);
  sqlite3_bind_int64(m_put_stmt, 1, static_cast<sqlite3_int64>(key));
  sqlite3_bind_blob(m_put_stmt, 2, val.data(),
                    static_cast<int>(val.size() * (sizeof val.front())),
                    SQLITE_STATIC);
  auto ret = sqlite3_step(m_put_stmt);
  if (ret != SQLITE_DONE) {
    throw database_error(std::string("Error inserting key ") +
                         std::to_string(key) + " into database '" + filename() +
                         "'.");
  }
  sqlite3_reset(m_put_stmt);
}

std::vector<Weight> KeyValueSQLite::get(uint64_t key, bool fatal) {
  std::vector<Weight> val{};
  // BLOB KEY
  // sqlite3_bind_blob(m_get_stmt, 1, &key, static_cast<int>(sizeof key),
  //                   SQLITE_STATIC);
  sqlite3_bind_int64(m_get_stmt, 1, static_cast<sqlite3_int64>(key));
  auto ret = sqlite3_step(m_get_stmt);
  if (ret == SQLITE_ROW) {
    auto elems = static_cast<std::size_t>(sqlite3_column_bytes(m_get_stmt, 0) /
                                          (sizeof val.front()));
    auto data =
        reinterpret_cast<const Weight *>(sqlite3_column_blob(m_get_stmt, 0));
    val.assign(data, data + elems);
  }
  else if (fatal) {
    throw database_error(std::string("Error retrieving key ") +
                         std::to_string(key) + " from database '" + filename() +
                         "'.");
  }
  sqlite3_reset(m_get_stmt);
  return val;
}

void KeyValueSQLite::remove(uint64_t key, bool fatal) {
  sqlite3_bind_int64(m_rm_stmt, 1, static_cast<sqlite3_int64>(key));
  auto ret = sqlite3_step(m_rm_stmt);
  if (ret != SQLITE_DONE && fatal) {
    throw database_error(std::string("Error removing key ") +
                         std::to_string(key) + " from database '" + filename() +
                         "'.");
  }
  sqlite3_reset(m_rm_stmt);
}

std::string KeyValueSQLite::filename() const {
  auto ptr = sqlite3_db_filename(m_db, "main");
  if (ptr) {
    return std::string(ptr);
  }
  return ":nofile:";
}

void KeyValueSQLite::execute(const std::string &cmd) {
  auto ret = sqlite3_exec(m_db, cmd.c_str(), nullptr, nullptr, nullptr);
  if (ret != SQLITE_OK) {
    throw database_error(std::string("Database '") + filename() +
                         "': " + "Error executing command '" + cmd + "'.");
  }
}

/**************
 * KeyValueKC *
 **************/

#ifdef PYRED_KCDB
KeyValueKC::KeyValueKC(const std::string &filename, std::size_t sz,
                       bool overwrite)
    : m_filename{filename} {
  using kyotocabinet::PolyDB;
  if (file_exists(filename)) {
    if (overwrite) {
      std::remove(filename.c_str());
    }
    else {
      throw database_error(std::string("Database file '") + filename +
                           "' already exists.");
    }
  }
  if (!m_db.open(filename + "#bnum=" + std::to_string(sz),
                 PolyDB::OWRITER | PolyDB::OCREATE | PolyDB::OTRUNCATE)) {
    throw database_error(std::string("Error opening database '") + filename +
                         "': " + m_db.error().name() + ".");
  }
}

KeyValueKC::~KeyValueKC() {
  if (!m_db.close()) {
    Config::log(1) << "Warning: error closing database '" << m_filename << "'."
                   << std::endl;
  }
  if (!m_filename.empty() && m_filename.front() != ':') {
    std::remove(m_filename.c_str());
  }
}

void KeyValueKC::put(uint64_t key, const std::vector<Weight> &val) {
  auto success = m_db.set(reinterpret_cast<char *>(&key), sizeof key,
                          reinterpret_cast<const char *>(val.data()),
                          val.size() * (sizeof val.front()));
  if (!success) {
    throw database_error(std::string("Error inserting key ") +
                         std::to_string(key) + " into database '" + filename() +
                         "'.");
  }
}

std::vector<Weight> KeyValueKC::get(uint64_t key, bool fatal) {
  std::vector<Weight> val{};
  std::size_t buf_size;
  auto buf = m_db.get(reinterpret_cast<char *>(&key), sizeof key, &buf_size);
  if (buf) {
    auto elems = buf_size / (sizeof val.front());
    auto data = reinterpret_cast<const Weight *>(buf);
    val.assign(data, data + elems);
    delete[] buf;
  }
  else if (fatal) {
    throw database_error(std::string("Error retrieving key ") +
                         std::to_string(key) + " from database '" + filename() +
                         "'.");
  }
  return val;
}

void KeyValueKC::remove(uint64_t key, bool fatal) {
  if (!m_db.remove(reinterpret_cast<char *>(&key), sizeof key) && fatal) {
    throw database_error(std::string("Error removing key ") +
                         std::to_string(key) + " from database '" + filename() +
                         "'.");
  }
}
#endif // PYRED_KCDB

} // namespace keyvaluedb
