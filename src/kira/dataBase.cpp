/*
Copyright (C) 2017-2025 The Kira Developers (see AUTHORS file)
Modified 2026 for BootLoops by Matthew D. Schwartz; distributed under the GNU GPL v3 or later, see COPYING. The changes are described in PATCHES.md.

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

#include "kira/dataBase.h"
#include "kira/exceptions.h"
#include "kira/integral.h"
#include "kira/kira.h"
#include "kira/tools.h"

using namespace std;
using namespace GiNaC;

static Loginfo& logger = Loginfo::instance();

bool sortDBintegral(const DBintegral& l, const DBintegral& r) {
  if (l.id > r.id)
    return true;
  else if (l.id < r.id)
    return false;

  return false;
}

DataBase::DataBase(string name) {
  errorMessage = 0;

  // check if the database .db file exists
  //   rc = sqlite3_open_v2(name.c_str(), &db,SQLITE_OPEN_READONLY,NULL);

  rc = sqlite3_open(name.c_str(), &db);

  if (rc) {
    sqlite3_close(db);
    throw ExceptionDisk("Cannot open database + '" + std::string(sqlite3_errmsg(db)) + "'");
  }

  string equation("EQUATION");
  string version("VERSION");

  int checkReduction = checkTable(equation);
  int checkVersion = checkTable(version);

  if(checkReduction == 0 && checkVersion == 0){
    create_version_table();
    save_version_number();
    logger << "Kira database version number: " << get_version_number() << "\n";
  }

//   else {
//     fprintf(stderr, "\nOpened database successfully\n");
//   }
}

void DataBase::open(string name) {
  errorMessage = 0;

  rc = sqlite3_open(name.c_str(), &db);

  if (rc) {
    sqlite3_close(db);
    throw ExceptionDisk("Cannot open database + '" + std::string(sqlite3_errmsg(db)) + "'");
  }

  string equation("EQUATION");
  string version("VERSION");

  int checkReduction = checkTable(equation);
  int checkVersion = checkTable(version);

  if(checkReduction == 0 && checkVersion == 0){
    create_version_table();
    save_version_number();
    logger << "Kira database version number: " << get_version_number() << "\n";
  }

//   else {
//     fprintf(stderr, "\nOpened database successfully\n");
//   }
}

DataBase::~DataBase() { sqlite3_close_v2(db); }

int DataBase::checkTable(string& name) {
  string sql = "select 1 from sqlite_master where type='table' and tbl_name='" +
               name + "';";

  sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL);

  sqlite3_step(stmt);

  auto column_text = sqlite3_column_text(stmt, 0);

  if (column_text == nullptr) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 0;
  }

  int token = something_int(column_text);

  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);

  return token;
}

void DataBase::merge_databases(std::vector<std::string>& files) {
  create_equation_table();
  create_integral_ordering_table();
  create_weight_bits_table();

  for (size_t it = 0; it < files.size(); it++) {
    string sql =
        "attach '" + files[it] + "' as toMerge" + to_string(it) + ";\n";
    sql += "BEGIN;\n";

    sql += "insert or ignore into EQUATION select * from toMerge" +
           to_string(it) + ".equation;\n";

    sql += "insert or ignore into INTEGRALORDERING select * from toMerge" +
           to_string(it) + ".INTEGRALORDERING;\n";

    sql += "insert or ignore into WEIGHTBITS select * from toMerge" +
           to_string(it) + ".WEIGHTBITS;\n";

    sql += "commit;\n";
    execute_statement(sql, "");
    if (!sqlite3_get_autocommit(db)) {
      logger << "This happens if you list results/kira.db in your current "
                "directory as a file to merge. This causes a lock. This lock "
                "will now be lifted. So the next merge will be ok. Kira merges "
                "all .db files into results/kira.db.\n";
      sql = "rollback;\n";
      execute_statement(sql, "Lock issue resolved: " + files[it] + "\n");
    }

    sql = "detach toMerge" + to_string(it) + ";\n";
    execute_statement(sql, "Processed: " + files[it] + "\n");
  }
}

void DataBase::create_integral_table(int /*lengthIndices*/) {
  string sql = "CREATE TABLE INTEGRAL("
               "ID            TEXT   NOT NULL,";
  sql += "INDICES       TEXT   NOT NULL,"
         "SECTOR        INT    NOT NULL,"
         "TOPOLOGY      INT    NOT NULL,"
         "DENCOUNT      INT    NOT NULL,"
         "DOTS          INT    NOT NULL,"
         "NUMS          INT    NOT NULL,"
         "AUXILIARY     INT    NOT NULL,"
         "UNIQUE (ID),"
         "PRIMARY KEY (INDICES,TOPOLOGY));";

  columnName.push_back("ID");
  columnName.push_back("INDICES");
  columnName.push_back("SECTOR");
  columnName.push_back("TOPOLOGY");
  columnName.push_back("DENCOUNT");
  columnName.push_back("DOTS");
  columnName.push_back("NUMS");
  columnName.push_back("AUXILIARY");

  //   totalColumn=NDBINT;
  //   if(totalColumn!=columnName.size()){
  //     logger << __FILE__ << "something is wrong\n";
  //     exit(-1);
  //   }

  execute_statement(sql, "Integral table created successfully\n");
}

void DataBase::create_skipid_table() {
  string sql = "CREATE TABLE SKIPID("
               "ID              TEXT    NOT NULL,"
               "PRIMARY KEY (ID));";

  execute_statement(sql, "SKIPID table created successfully\n");
}

void DataBase::create_equation_table() {
  string sql = "CREATE TABLE EQUATION("
               "ID              TEXT    NOT NULL,"
               "ARRAYID         TEXT    NOT NULL,"
               "AUXILIARY       INT    NOT NULL,"
               "SECTOR          INT    NOT NULL,"
               "LENGTH          INT    NOT NULL,"
               "COEFFICIENT     TEXT   NOT NULL,"
               "PRIMARY KEY (ID,ARRAYID));";

  execute_statement(sql, "Equation table created successfully\n");
}

void DataBase::create_weight_bits_table() {
  string sql = "CREATE TABLE WEIGHTBITS("
               "A              INT    NOT NULL,"
               "B              INT    NOT NULL,"
               "C              INT    NOT NULL,"
               "D              INT    NOT NULL,"
               "PRIMARY KEY (A));";

  execute_statement(sql, "WEIGHTBITS table created successfully\n");
}

void DataBase::save_weight_bits(std::vector<std::uint32_t>& weightBits) {

  string input = "INSERT INTO WEIGHTBITS VALUES (?1,?2,?3,?4)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  if (weightBits.size() != 4) {
    throw ExceptionInternal("DataBase::save_weight_bits(): Weight-bits length is not 4");
  }

  sqlite3_bind_int(stmt, 1, weightBits[0]);
  sqlite3_bind_int(stmt, 2, weightBits[1]);
  sqlite3_bind_int(stmt, 3, weightBits[2]);
  sqlite3_bind_int(stmt, 4, weightBits[3]);

  sqlite3_step(stmt);
  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
}

int DataBase::table_weight_bits_empty() {
  string input = "SELECT * FROM WEIGHTBITS";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  int res = sqlite3_step(stmt);

  if (res == SQLITE_DONE) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 1;
  }
  else if (res == SQLITE_ROW) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 0;
  }
  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
  return 1;
}

int DataBase::get_weight_bits(std::vector<std::uint32_t>& weightBits) {

  weightBits.clear();
  string input = "SELECT * FROM WEIGHTBITS";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  int res = sqlite3_step(stmt);

  if (res == SQLITE_ROW) {
    for (int ju = 0; ju < 4; ju++) {
      int number = sqlite3_column_int(stmt, ju);
      weightBits.push_back(number);
    }

    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 1;
  }
  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 0;
  }

  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
  return 0;
}

void DataBase::create_integral_ordering_table() {
  string sql = "CREATE TABLE INTEGRALORDERING("
               "ID              INT    NOT NULL,"
               "PRIMARY KEY (ID));";

  execute_statement(sql, "INTEGRALORDERING table created successfully\n");
}

void DataBase::save_integral_ordering(int integralOrdering) {
  string input = "INSERT INTO INTEGRALORDERING VALUES (?1)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  sqlite3_bind_int(stmt, 1, integralOrdering);

  sqlite3_step(stmt);
  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
}

int DataBase::get_integral_ordering() {
  string input = "SELECT * FROM INTEGRALORDERING";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  int res = sqlite3_step(stmt);
  if (res == SQLITE_ROW) {
    int integralOrdering = sqlite3_column_int(stmt, 0);
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return integralOrdering;
  }
  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 0;
  }

  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
  return 0;
}

void DataBase::create_firefly_table() {
  string sql = "CREATE TABLE FIREFLY("
               "ID              INT    NOT NULL,"
               "PRIMARY KEY (ID));";

  execute_statement(sql, "FIREFLY table created successfully\n");
}

void DataBase::save_firefly_flag(int flagFiryfly) {
  string input = "INSERT INTO FIREFLY VALUES (?1)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  sqlite3_bind_int(stmt, 1, flagFiryfly);

  sqlite3_step(stmt);
  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
}

int DataBase::get_firefly_flag() {
  string input = "SELECT * FROM FIREFLY";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  int res = sqlite3_step(stmt);
  if (res == SQLITE_ROW) {
    int fireflyFlag = sqlite3_column_int(stmt, 0);
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return fireflyFlag;
  }
  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return 0;
  }

  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
  return 0;
}


void DataBase::create_version_table() {
  string sql = "CREATE TABLE VERSION("
               "number              TEXT    NOT NULL,"
               "PRIMARY KEY (number));";

  execute_statement(sql, "VERSION table created successfully\n");
}

void DataBase::save_version_number() {
  string input = "INSERT INTO VERSION VALUES (?1)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  string version("2.1");

  sqlite3_bind_text(stmt, 1, version.c_str(), version.size(), SQLITE_TRANSIENT);

  sqlite3_step(stmt);
  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
}

string DataBase::get_version_number() {
  string input = "SELECT * FROM VERSION";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);

  string version = "0";

  int res = sqlite3_step(stmt);

  if (res == SQLITE_ROW) {
    version = (char*)sqlite3_column_text(stmt, 0);
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return version;
  }
  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    sqlite3_finalize(stmt);
    return version;
  }

  sqlite3_reset(stmt);
  sqlite3_finalize(stmt);
  return version;
}


void DataBase::attach_table(std::string& attachName) {
  string sql = " attach '" + attachName + "' as toMerge;";

  execute_statement(sql, "1. ");

  sql = " begin;";

  execute_statement(sql, "2. ");

  sql = " insert into equation select * from toMerge.equation;";

  execute_statement(sql, "3. ");

  sql = " commit;";

  execute_statement(sql, "4. ");

  sql = " detach toMerge;";

  execute_statement(sql, "Old results successfully loaded.\n");
}

void DataBase::execute_statement(string& sql, string message) {
  rc = sqlite3_exec(db, sql.c_str(), NULL, NULL, &errorMessage);

  if (rc != SQLITE_OK) {
    fprintf(stderr, "SQLite3: %s\n", errorMessage);
    sqlite3_free(errorMessage);
  }
  else {
    logger << message;
  }
}

void DataBase::begin_transaction() {
  // Several call sites issue begin/commit in tight loops; if a previous
  // BEGIN is still open (e.g. an inner step returned SQLITE_OK but the
  // outer loop re-enters), a second BEGIN would error.  Only open a new
  // transaction when SQLite reports autocommit (i.e. no txn active).
  if (!sqlite3_get_autocommit(db)) {
    return;
  }
  string sql = "BEGIN TRANSACTION";
  rc = sqlite3_exec(db, sql.c_str(), NULL, NULL, &errorMessage);

  if (rc != SQLITE_OK) {
    fprintf(stderr, "SQL error: %s\n", errorMessage);
    sqlite3_free(errorMessage);
  }
  //   else{
  //     logger<< message;
  //   }
}

void DataBase::commit_transaction() {
  // Guard the COMMIT with the connection's actual transaction state.
  // sqlite3_get_autocommit() == 1 means no explicit transaction is open,
  // so a COMMIT here would emit "cannot commit - no transaction is
  // active" (printed on essentially every commit during long reductions
  // when an intermediate sqlite3_step() has already auto-rolled-back on
  // an error).  In that case we have nothing to commit and should only
  // finalize the prepared statement.  Non-fatal but pure log noise; the
  // guard makes the call idempotent.
  if (!sqlite3_get_autocommit(db)) {
    string sql = "COMMIT TRANSACTION";
    rc = sqlite3_exec(db, sql.c_str(), NULL, NULL, &errorMessage);

    if (rc != SQLITE_OK) {
      fprintf(stderr, "SQL error: %s\n", errorMessage);
      sqlite3_free(errorMessage);
    }
    //   else{
    //     logger<< message;
    //   }
  }

  sqlite3_finalize(stmt);
}

void DataBase::prepare_integral() {
  string input = "INSERT INTO INTEGRAL VALUES (";

  int jj = 1;
  for (int i = 0; i < (NDBINT - 1); i++) {
    input += "?" + int_string(jj++) + ",";
  }

  input += "?" + int_string(jj++) + ")";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::prepare_pyred() {
  string input = "SELECT ID FROM INTEGRAL WHERE ";

  //   input+=columnName[0];//id
  //   input+="=?"+int_string(1);

  //   input+=" AND ";

  input += columnName[3]; // topology
  input += "=?" + int_string(1);

  input += " AND ";

  //   input+=columnName[2];//sector
  //   input+="&?"+int_string(2) + "="+columnName[2];

  input += columnName[2]; // sector
  input += " & ~?" + int_string(2) + "==0";

  input += " AND ";

  input += columnName[5]; // r
  //   input+="-";
  //   input+=columnName[4];//t
  input += "<=?" + int_string(3);

  input += " AND ";

  input += columnName[6]; // s
  input += "<=?" + int_string(4);

  input += " AND ";

  input += columnName[7]; // auxiliary
  input += "=?" + int_string(5);

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::prepare_id() {
  string input = "SELECT * FROM INTEGRAL WHERE ";

  input += columnName[1];
  input += "=?" + int_string(1);
  input += " AND ";
  input += columnName[3];
  input += "=?" + int_string(2);

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::select_skipid() {
  string input = "SELECT * FROM SKIPID WHERE ID=?1";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

// void DataBase::prepare_lookup_id(){
//
//   string input = "SELECT * FROM EQUATION INNER JOIN INTEGRAL ON
//   EQUATION.ARRAYID=INTEGRAL.ID WHERE EQUATION.ID=?";
//
//   sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
// }

void DataBase::prepare_lookup_id() {
  string input = "SELECT * FROM EQUATION WHERE ID=?";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::prepare_backsubstitution() {
  string input = "INSERT INTO EQUATION VALUES (?1, ?2, ?3, ?4, ?5, ?6)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::prepare_skipid() {
  string input = "INSERT INTO SKIPID VALUES (?1)";

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

void DataBase::bind_integral(SEEDIntegral& integral, int lengthIndices) {
  sqlite3_bind_text(stmt, 1, integral.id.to_string().c_str(), -1,
                    SQLITE_TRANSIENT);

  string indices;
  for (int i = 0; i < lengthIndices; i++) {
    indices += int_string(integral.indices[i]);
    if (i != lengthIndices - 1) indices += ",";
  }
  sqlite3_bind_text(stmt, 2, indices.c_str(), indices.size(), SQLITE_TRANSIENT);

  sqlite3_bind_int(stmt, 3, integral.characteristics[SECTOR]);
  sqlite3_bind_int(stmt, 4, integral.characteristics[TOPOLOGY]);
  sqlite3_bind_int(stmt, 5, integral.characteristics[DENCOUNT]);
  sqlite3_bind_int(stmt, 6, integral.characteristics[DOTS]);
  sqlite3_bind_int(stmt, 7, integral.characteristics[NUM]);
  sqlite3_bind_int(stmt, 8, integral.flag2);

  sqlite3_step(stmt);
  sqlite3_reset(stmt);
}

void DataBase::bind_denominators(pyred::Weight id, std::string& coefficient,
                                 int length, pyred::Weight ID) {
  tuple<string, unsigned, unsigned> properties;

  get_properties(id, properties);

  sqlite3_bind_text(stmt, 1, ID.to_string().c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, id.to_string().c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt, 3, 0);
  sqlite3_bind_int(stmt, 4, get<2>(properties));
  sqlite3_bind_int(stmt, 5, length);
  sqlite3_bind_text(stmt, 6, coefficient.c_str(), coefficient.size(),
                    SQLITE_STATIC);

  sqlite3_step(stmt);

  if (sqlite3_step(stmt) != SQLITE_DONE) {
    sqlite3_free(errorMessage);
  }

  sqlite3_reset(stmt);
}

void DataBase::bind_equation(
    BaseIntegral& integral, int length, pyred::Weight ID,
    std::unordered_map<pyred::Weight, pyred::Weight> /*mastersReMap*/) {
  //   auto mapContent = mastersReMap.find(ID);
  //
  //   if(mapContent != mastersReMap.end()){
  //
  //     sqlite3_bind_text(stmt, 1, to_string(mapContent->second).c_str(), -1,
  //     SQLITE_TRANSIENT);
  //   }
  //   else{
  //     sqlite3_bind_text(stmt, 1, to_string(ID).c_str(), -1,
  //     SQLITE_TRANSIENT);
  //   }
  sqlite3_bind_text(stmt, 1, ID.to_string().c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, integral.id.to_string().c_str(), -1,
                    SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt, 3, integral.flag2);
  sqlite3_bind_int(stmt, 4, integral.characteristics[SECTOR]);
  sqlite3_bind_int(stmt, 5, length);
  sqlite3_bind_text(stmt, 6, integral.coefficientString.c_str(),
                    integral.coefficientString.size(), SQLITE_STATIC);

  sqlite3_step(stmt);

  if (sqlite3_step(stmt) != SQLITE_DONE) {
    sqlite3_free(errorMessage);
  }

  sqlite3_reset(stmt);
}

void DataBase::bind_zero_equation(pyred::Weight ID, std::uint32_t sector) {
  std::string coef = "1";
  sqlite3_bind_text(stmt, 1, ID.to_string().c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, ID.to_string().c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt, 3, 1);
  sqlite3_bind_int(stmt, 4, sector);
  sqlite3_bind_int(stmt, 5, 1);
  sqlite3_bind_text(stmt, 6, coef.c_str(),
                    coef.size(), SQLITE_STATIC);

  sqlite3_step(stmt);

  if (sqlite3_step(stmt) != SQLITE_DONE) {
    sqlite3_free(errorMessage);
  }

  sqlite3_reset(stmt);
}

void DataBase::bind_skipid(pyred::Weight ID) {
  sqlite3_bind_text(stmt, 1, ID.to_string().c_str(), -1, SQLITE_TRANSIENT);

  sqlite3_step(stmt);

  if (sqlite3_step(stmt) != SQLITE_DONE) {
    sqlite3_free(errorMessage);
  }

  sqlite3_reset(stmt);
}

int DataBase::number_of_columns() {
  string sql = "SELECT * FROM INTEGRAL;";

  sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL);

  totalColumn = sqlite3_column_count(stmt);

  sqlite3_step(stmt);

  for (unsigned int iterator = 0; iterator < totalColumn; iterator++) {
    columnName.push_back(sqlite3_column_name(stmt, iterator));
    sqlite3_reset(stmt);
  }

  sqlite3_finalize(stmt);
  if (totalColumn == 0) logger << "The database is yet empty, skip it\n";
  return totalColumn;
}

void DataBase::prepare_find_master() {
  string input = "SELECT * FROM INTEGRAL WHERE ";

  input += columnName[0]; // id
  input += "=?" + int_string(1);

  sqlite3_prepare_v2(db, input.c_str(), input.size(), &stmt, NULL);
}

int DataBase::bind_id_get_integral(
    pyred::Weight id, tuple<string, unsigned, unsigned>& integral) {
  sqlite3_bind_text(stmt, 1, id.to_string().c_str(), -1, SQLITE_STATIC);

  int res = sqlite3_step(stmt);
  if (res == SQLITE_ROW) {
    get<0>(integral) = (char*)sqlite3_column_text(stmt, 1);

    unsigned token;
    std::string unsignedLong = (char*)sqlite3_column_text(stmt, 3);
    std::istringstream ss(unsignedLong);
    ss >> token;
    get<1>(integral) = token;
    sqlite3_reset(stmt);
    return 1;
  }
  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    return 0;
  }
  sqlite3_reset(stmt);
  return 0;
}

int DataBase::bind_id_get_BSequation(
    pyred::Weight id, BaseIntegral*& integral,
    std::vector<pyred::Weight>& masterVectorSkip,
    std::unordered_map<pyred::Weight, unsigned>& occurrence, int flagOcc) {
  sqlite3_bind_text(stmt, 1, id.to_string().c_str(), -1, SQLITE_TRANSIENT);
  int res = sqlite3_step(stmt);

  if (res == SQLITE_ROW) {

    uint32_t eqLength = 0;

    vector<BaseIntegral> vecBaseIntegral;
    while (1) {
      BaseIntegral TMPintegral;
      //       integral[j].id = std::stoul((char*)sqlite3_column_text(stmt,1));
      TMPintegral.id = pyred::Weight((char*)sqlite3_column_text(stmt, 1));

      tuple<string, unsigned, unsigned> properties;
      //       get_properties(integral[j].id,properties);
      get_properties(TMPintegral.id, properties);

      //       integral[j].characteristics[TOPOLOGY] = get<1>(properties);
      //       integral[j].characteristics[SECTOR] =  get<2>(properties);
      //       //sector
      //
      //       integral[j].length = eqLength;
      //       integral[j].coefficientString =
      //       (char*)sqlite3_column_text(stmt,5); integral[j].id =
      //       std::stoul((char*)sqlite3_column_text(stmt,1)); integral[j].flag2
      //       = sqlite3_column_int(stmt,2);

      TMPintegral.characteristics[TOPOLOGY] = get<1>(properties);
      TMPintegral.characteristics[SECTOR] = get<2>(properties); // sector

      //       TMPintegral.length = eqLength;
      TMPintegral.coefficientString = (char*)sqlite3_column_text(stmt, 5);
      TMPintegral.id = pyred::Weight((char*)sqlite3_column_text(stmt, 1));
      TMPintegral.flag2 = sqlite3_column_int(stmt, 2);

      //       auto findSkipMaster = find(masterVectorSkip.begin(),
      //       masterVectorSkip.end(), integral[j].id);
      auto findSkipMaster = find(masterVectorSkip.begin(),
                                 masterVectorSkip.end(), TMPintegral.id);

      //       if ( findSkipMaster != masterVectorSkip.end() ) {
      // 	integral[j].coefficientString = "0";
      //       }
      if (findSkipMaster != masterVectorSkip.end()) {
        TMPintegral.coefficientString = "0";
      }

      if (flagOcc) {
        auto occContent = occurrence.find(TMPintegral.id);

        if (occContent != occurrence.end()) {
          (occContent->second)++;
        }
        else {
          occurrence.insert(pair<pyred::Weight, unsigned>(TMPintegral.id, 1));
        }
      }

      //       j--;

      res = sqlite3_step(stmt);
      eqLength++;
      vecBaseIntegral.push_back(TMPintegral);
      if (res == SQLITE_DONE || res == SQLITE_ERROR) {
        break;
      }
    }
    sqlite3_reset(stmt);

    integral = new BaseIntegral[eqLength];

    for (int it = eqLength - 1; it >= 0; it--) {
      integral[it] = vecBaseIntegral[it];
      integral[it].length = eqLength;
    }

    if (eqLength > 0) {
      for (unsigned it = 0; it < eqLength - 1; it++) {
        for (unsigned itt = it; itt < eqLength; itt++) {
          if ((integral[itt].id) > (integral[it].id)) {
            swap(integral[it], integral[itt]);
          }
        }
      }
    }

    return 1;
  }
  else {
    sqlite3_reset(stmt);
    return 0;
  }
}

void DataBase::bind_id_get_equation(pyred::Weight id,
                                    vector<DBintegral>& equation) {
  sqlite3_bind_text(stmt, 1, id.to_string().c_str(), -1, SQLITE_TRANSIENT);

  while (1) {
    int res = sqlite3_step(stmt);

    if (res == SQLITE_ROW) {
      DBintegral integral;

      integral.id = pyred::Weight((char*)sqlite3_column_text(stmt, 1));

      tuple<string, unsigned, unsigned> properties;
      get_properties(integral.id, properties);

      integral.topology = get<1>(properties);
      integral.indices = get<0>(properties);
      //     get<2>(properties); //sector

      int massDimension = 0;
      size_t found;
      int indices;
      string line = integral.indices;
      while ((found = line.find_first_of(",")) != string::npos) {
        istringstream(line.substr(0, found)) >> indices;
        line = line.substr(found + 1);
        massDimension += indices;
      }
      istringstream(line) >> indices;
      massDimension += indices;

      //     integral.id = std::stoul((char*)sqlite3_column_text(stmt, 1));

      integral.massDimension = massDimension;
      integral.coefficient = (char*)sqlite3_column_text(stmt, 5);
      if (integral.id == id)
        equation.insert(equation.begin(), integral);
      else
        equation.push_back(integral);
    }
    if (res == SQLITE_DONE || res == SQLITE_ERROR) {
      break;
    }
  }
  if (equation.size() > 1)
    sort((equation.begin() + 1), equation.end(), sortDBintegral);

  sqlite3_reset(stmt);
}

int DataBase::find_skipid(pyred::Weight id) {
  sqlite3_bind_text(stmt, 1, id.to_string().c_str(), -1, SQLITE_STATIC);

  int res = sqlite3_step(stmt);
  if (res == SQLITE_ROW) {
    sqlite3_reset(stmt);

    return 1;
  }

  if (res == SQLITE_DONE || res == SQLITE_ERROR) {
    sqlite3_reset(stmt);
    return 0;
  }
  return 0;
}

void DataBase::finalize() {sqlite3_finalize(stmt); /*sqlite3_close(db);*/}

void DataBase::bind_get_answer(int topology, int sector, int dots, int nums,
                               std::vector<pyred::Weight>& mandatory) {
  sqlite3_bind_int(stmt, 1, topology);
  sqlite3_bind_int(stmt, 2, sector);
  sqlite3_bind_int(stmt, 3, dots);
  sqlite3_bind_int(stmt, 4, nums);
  sqlite3_bind_int(stmt, 5, 0);

  while (1) {
    int res = sqlite3_step(stmt);

    if (res == SQLITE_ROW) {
      pyred::Weight token;
      std::string unsignedLong = (char*)sqlite3_column_text(stmt, 0);
      std::istringstream ss(unsignedLong);
      ss >> token;
      mandatory.push_back(token);
    }
    if (res == SQLITE_DONE || res == SQLITE_ERROR) {
      break;
    }
  }

  sqlite3_reset(stmt);
}

pyred::Weight DataBase::bind_get_id(vector<int>& indices, int topologyNumber) {
  string indicesSTRING;
  for (unsigned int iterator = 0; iterator < indices.size(); iterator++) {
    indicesSTRING += int_string(indices[iterator]);
    if (iterator != indices.size() - 1) indicesSTRING += ",";
  }
  sqlite3_bind_text(stmt, 1, indicesSTRING.c_str(), indicesSTRING.size(),
                    SQLITE_STATIC);

  sqlite3_bind_int(stmt, 2, topologyNumber);

  int row = sqlite3_step(stmt);

  if (row == SQLITE_DONE || row == SQLITE_ERROR) {
    logger << "This integral is not inside the range of the reduction "
              "parameters (S,r,s): \n";
    for (unsigned int iterator = 0; iterator < indices.size(); iterator++) {
      logger << indices[iterator] << " ";
    }
    logger << "\n";
    sqlite3_reset(stmt);
    return pyred::Weight(0u);
  }

  if (something_int(sqlite3_column_text(stmt, 7)) != 0) {
    logger << "This integral is not inside the range of the reduction "
              "parameters (S,r,s) and \n";
    for (unsigned int iterator = 0; iterator < indices.size(); iterator++) {
      logger << indices[iterator] << " ";
    }
    logger << "\n";
    sqlite3_reset(stmt);
    return pyred::Weight(0u);
  }

  pyred::Weight id;
  std::string unsignedLong = (char*)sqlite3_column_text(stmt, 0);
  std::istringstream ss(unsignedLong);
  ss >> id;

  sqlite3_reset(stmt);
  return id;
}
