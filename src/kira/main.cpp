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

#include <chrono>
#include <ctime>
#include <getopt.h>

#include "pyred/config.h"
#include "pyred/defs.h"
#include "pyred/parallel.h"
#include "kira/exceptions.h"
#include "kira/kira.h"
#include "kira/tools.h"

#ifdef KIRAFIREFLY
#include "firefly/config.hpp" // for WITH_MPI, FLINT
#ifdef WITH_MPI
#include <mpi.h>
#include "firefly/MPIWorker.hpp"
#include "kira/mpi_work.h"
#endif
#endif


using namespace std;

static Loginfo& logger = Loginfo::instance();
std::time_t timestamp;

class ParseOptions {
private:
  std::string version_info() {
    std::string buildtype;
    bool weight128 = false;
    bool with_firefly = false;
    bool firefly_subproject = false;
    bool with_flint = false;
    bool with_mpi = false;
    int with_jemalloc = -1;
#ifdef BUILDTYPE
    buildtype = BUILDTYPE;
    with_jemalloc = 0;
#else
    buildtype = "unknown (use the Meson build system!)";
#endif
#ifdef STATICBUILD
    buildtype = buildtype + " (static)";
#endif
#ifdef WITH_MPI
    with_mpi = true;
#endif
#ifdef WEIGHT128
    weight128 = true;
#endif
#ifdef JEMALLOC
    with_jemalloc = 1;
#endif
#ifdef KIRAFIREFLY
    with_firefly = true;
#endif
#ifdef FIREFLYSUBPROJECT
    firefly_subproject = true;
#endif
#ifdef FLINT
    with_flint = true;
#endif
    std::string s = "Kira version: " VERSION;
#ifdef GITVERSION
    s = s + " (Git: " + GITVERSION + ")";
#endif
    s = s + "\nBuild configuration: " + buildtype + "; " +
    (with_mpi ? "MPI; " : "non-MPI; ") +
    (weight128 ? "128 bit weights; " : "64 bit weights; ");
    s = s + (with_jemalloc >= 0 ?
      (with_jemalloc ? "with jemalloc;" : "without jemalloc;") : "") + "\n";
    if (with_firefly) {
      s = s + "with " +
      (firefly_subproject ? "FireFly as subproject" : "external FireFly") +
      (with_flint ? " with FLINT\n" : " without FLINT\n");
    }
    else {
      s = s + "without FireFly\n";
    }
    return s;
  }

public:
  ParseOptions(int argc, char* argv[]) {
    silent_flag = 0;
    pyred_config = "";
    parallel = 1;
    integralOrdering = 9;
    int c;

    static struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        {"silent", no_argument, &silent_flag, 1},
        {"set_value", required_argument, nullptr, 's'},
        {"set_sector", required_argument, nullptr, 'S'},
        {"force_database_format", required_argument, nullptr, 'd'},
        {"trim", required_argument, nullptr, 't'},
        {"dir_num", required_argument, nullptr, 'w'},
        {"auxiliary_name", required_argument, nullptr, 'a'},
        {"parallel", optional_argument, nullptr, 'p'},
        //{"prime_n_number", required_argument, nullptr, 'P'},
        {"prime", required_argument, nullptr, 'P'},
        {"integral_ordering", required_argument, nullptr, 'i'},
        {"pyred_config", required_argument, nullptr, 'c'},
        {"log_time_stamp", no_argument, nullptr, 'l'},
#ifdef KIRAFIREFLY
        {"bunch_size", required_argument, nullptr, 'b'},
#endif
        {nullptr, 0, nullptr, 0} // mark end of array
    };

    bool lFound = false;
    // Step 1: Manually scan for -l first
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--log_time_stamp") == 0) {
            char timeString[20];
            strftime(timeString, 20, "%F_%T", gmtime(&timestamp));
            logger.set_logfile("kira_" + string(timeString) + ".log");
            lFound = true;
            break;
        }
    }
    if (!lFound) {
        logger.set_logfile("kira.log");
    }

    while (1) {
      int option_index = 0;
#ifndef KIRAFIREFLY
      c = getopt_long_only(argc, argv, "hvla:p::i:c:s:t:w:a:S:d:P:", long_options,
                            &option_index);
#else
      c = getopt_long_only(argc, argv, "hvla:p::i:c:s:t:w:a:S:d:b:P:", long_options,
                            &option_index);
#endif
      if (c == 'l') {
        continue;
      }

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c) {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0) break;
          printf("option %s", long_options[option_index].name);
          if (optarg) printf(" with arg %s", optarg);
          printf("\n");
          break;
        case 'k':
          break;
        case 'h':
          PrintUsage();
          throw ExceptionQuit();
        case 'v':
          logger << version_info();
          throw ExceptionQuit();
        case 'p': {
          string optstr = "physical";
          if (optarg) {
            optstr = string(optarg);
          }
          parallel = -1;
          if (optstr == "physical") {
            // TODO: in case of MPI runs, set this per node.
            parallel = pyred::cpu_count();
          }
          else if (optstr == "logical") {
            parallel = thread::hardware_concurrency();
          }
          else {
            try {
              parallel = pyred::string_to_int(optstr);
            }
            catch (const pyred::input_error &) {}
          }
          if (parallel > 0)
            logger << parallel
                    << " Fermat/Blackbox instances will be run in parallel.\n";
          else {
            throw ExceptionCommandLine("parallel (p): Invalid argument");
          }
        } break;
        case 'i':
          integralOrdering = atoi(optarg);
          if (1 > integralOrdering || integralOrdering > 8) {
            throw ExceptionCommandLine("integral_ordering (i): Invalid argument");
          }
          break;
        case 's':
          set_variables.push_back(optarg);
          logger <<"Set the symbol: " << optarg << "\n";
          break;
        case 'S':
          set_sector = optarg;
          logger <<"Force Kira to do the reduction of sector: " << optarg << "\n";
          break;
        case 'P':
          prime_user = optarg;
          logger <<"Force Kira to generate the system modulo prime number nth: " << optarg << "\n";
          break;
        case 'd':
          database_format = optarg;
          if(database_format != "fermat" && database_format != "firefly"){
            throw ExceptionCommandLine("force_database_format (d): Invalid argument");
          }
          logger <<"Set database format to: " << optarg << "\n";
          break;
        case 't':
          trim.push_back(optarg);
          break;
        case 'w':
          pathNum = optarg;
          break;
        case 'a':
          auxName = optarg;
          break;
        case 'c':
          pyred_config += optarg;
          break;
#ifdef KIRAFIREFLY
        case 'b':
          bunch_size = atoi(optarg);

          if ((bunch_size & (bunch_size - 1)) == 0 && bunch_size <= 128) {
            logger << "FireFly will use " << bunch_size
                    << " as maximum bunch size.\n";
          }
          else {
            throw ExceptionCommandLine("The option bunch_size accepts only numbers of powers of 2 up to 128\n");
          }
          break;
#endif
        case '?':
        default:
          PrintUsage();
          abort();
      }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind != (argc - 1)) {
      throw ExceptionCommandLine("Unexpected number of command line arguments");
    }

#if defined(KIRAFIREFLY) && defined(WITH_MPI)
    int process;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);

    if (0 == file_exists(argv[optind]) && process == firefly::master) {
#else
    if (0 == file_exists(argv[optind])) {
#endif
      throw ExceptionCommandLine("Cannot find '" + job_name + "'");
    }

    job_name = argv[optind];

    if (silent_flag) {
      logger.silent_cout_buf();
    }
  }

  ~ParseOptions() {
    if (silent_flag) {
      logger.restore_silence(); // restore the original stream buffer
    }
  }

  int getCores() { return (parallel); }

  int get_integral_ordering() { return (integralOrdering); }

  std::string getPyred() { return (pyred_config); }

  string getJobname() { return (job_name); }

#ifdef KIRAFIREFLY
  uint32_t get_bunch_size() { return static_cast<uint32_t>(bunch_size); }
#endif

  vector<string> get_set_variables() { return (set_variables); }

  string get_set_sector() { return (set_sector); }

  string get_prime_user() { return (prime_user); }

  string get_database_format() { return (database_format); }

  vector<string> get_trim() { return (trim); }
  string get_path_num() { return (pathNum); }
  string get_aux_name() { return (auxName);}

  int get_silent() {return silent_flag;}

  void PrintUsage() {
    logger << "\n";
    logger << "Usage:\n";
    logger << "  kira FILE        Run the job specified in the job file FILE."
              "\n";
    logger << "                   The project directory is set"
              " to the working directory.\n";
    logger << "Options:\n";
    logger << "   --version       Print the Kira version and exit.\n";
    logger << "   --help          Print this help and exit.\n";
    logger << "   --silent        Switch to silent mode.\n";
    logger << "   --parallel=n    Run n instances of Fermat rsp. FireFly \n"
              "                   black-box evaluators in parallel. \n"
              "                   If the value is omitted or =physical, use \n"
              "                   all CPU cores. \n"
              "                   With =logical, all logical cores are used.\n";
#ifdef KIRAFIREFLY
    logger << "   --bunch_size=n  Set FireFly's maximum bunch size to n.\n";
#endif
    logger << "   --integral_ordering=n    Use the n-th integral ordering \n"
              "                   (possible values: 1,...,8).\n";
    logger << "   --set_value=\"var=val\"    Set the value of the symbol "
              "\"var\" to \"val\".\n";
    logger << "   --set_sector=n  Perform only the reduction of "
              "sector n.\n";
    logger << "   --prime_user=n  User makes free choice for a prime number\n"
              "                   maximum prime number supported is 9223372036854775783\n";

    logger << "\n";
    logger << "   --force_database_format=fermat|firefly  Updates the database to version 2.1\n"
              "                   and the correct format.\n";
    logger << "\n";
  }

private:
  int silent_flag;
  vector<string> set_variables;
  string set_sector;
  string prime_user;
  vector<string> trim;
  string pathNum, auxName;
  //   , algebra_flag;
  string job_name;
  string database_format;
  ofstream fout;
  int parallel;
  int integralOrdering;
  std::string pyred_config;
#ifdef KIRAFIREFLY
  int bunch_size = 1;
#endif
  //   std::string algebra_config;
  ParseOptions(){};
};

int main(int argc, char* argv[]) {
#if defined(KIRAFIREFLY) && defined(WITH_MPI)
  // Initialization of MPI processes
  // ---------------------------------------------------------------------------
  // int provided =
  MPIHelper::init();
  int process;
  MPI_Comm_rank(MPI_COMM_WORLD, &process);
  // ---------------------------------------------------------------------------

  if (process == firefly::master) {
#endif
    Clock clock;

    try {
      timestamp = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      ParseOptions parseoptions(argc, argv);

#ifndef KIRAFIREFLY
    Kira kira(parseoptions.getJobname(),
              parseoptions.getCores() /*,parseoptions.getAlgebra()*/,
              parseoptions.getPyred(), parseoptions.get_integral_ordering(),
              parseoptions.get_set_variables(), parseoptions.get_trim(),
              parseoptions.get_set_sector(),
              parseoptions.get_prime_user(),parseoptions.get_database_format(), parseoptions.get_path_num(), parseoptions.get_aux_name());
#else
  Kira kira(parseoptions.getJobname(),
            parseoptions.getCores() /*,parseoptions.getAlgebra()*/,
            parseoptions.getPyred(), parseoptions.get_integral_ordering(),
            parseoptions.get_set_variables(), parseoptions.get_trim(),
            parseoptions.get_set_sector(), parseoptions.get_prime_user(),
            parseoptions.get_bunch_size(), parseoptions.get_silent(),
            parseoptions.get_database_format(), parseoptions.get_path_num(),
            parseoptions.get_aux_name());
#endif
      kira.read_config();
      kira.execute_jobs();
    }
    catch (pyred::init_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::parser_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::zero_modular_division& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::input_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::runtime_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::uf_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (pyred::permutation_error& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionCommandLine& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionConfigFiles& e) {
      logger << e.what() << "\n\0";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionDisk& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionFermat& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionInternal& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionPipe& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionResume& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionRuntime& e) {
      logger << e.what() << "\n";
      logger << "Kira run aborted\n";
      return 1;
    }
    catch (ExceptionQuit& e) {
      return 0;
    }

#if defined(KIRAFIREFLY) && defined(WITH_MPI)
    bool cont = false;
    MPI_Bcast(&cont, 1, MPI_CXX_BOOL, firefly::master, MPI_COMM_WORLD);
#endif

    logger << "Total time:\n( " << clock.eval_time() << " s )\n";
#if defined(KIRAFIREFLY) && defined(WITH_MPI)
  }
  else {
    ParseOptions parseoptions(argc, argv);
    mpi_work(static_cast<uint32_t>(parseoptions.getCores()),
             parseoptions.get_bunch_size());
  }

  // MPI is finalised and the MPI weight type freed by ~MPIHelper().
#endif
  return 0;
}
