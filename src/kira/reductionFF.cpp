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

#include "kira/kira.h"

#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <iomanip>
#include <mutex>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <firefly/BaseReconst.hpp>
#include "firefly/config.hpp" // for WITH_MPI
#include <firefly/FFInt.hpp>
#include <firefly/FFIntVec.hpp>
#include <firefly/RationalFunction.hpp>
#include <firefly/RatReconst.hpp>
#include <firefly/ShuntingYardParser.hpp>
#include <firefly/utils.hpp>

#include "pyred/config.h" // for GZSTREAM_NAMESPACE
#include "pyred/coeff_int.h"
#include "pyred/defs.h"
#include "kira/tools.h"
#include "gzstream/gzstream.h"

#ifdef WITH_MPI
#include <mpi.h>

#include "kira/mpi_work.h" // MPIHelper
#include <firefly/MPIWorker.hpp>
#endif

static Loginfo& logger = Loginfo::instance();

static std::mutex reconst_done;
static std::condition_variable cond_reconst_done;

static std::mutex num_done;
static std::condition_variable cond_num_done;

void Kira::run_firefly(/*const*/ std::vector<pyred::Weight>& mandatory_vec,
    const uint32_t mode_, const std::string topology_name,
    const int factor_scan) { // mandatory not const due to MPI_Bcast
#ifdef WITH_MPI
  MPI_Datatype &mpi_weight_type = MPIHelper::weight_type();
#endif
  if (mode_ == 0) {
    reconstruction_mode = 0; // full reconstruction
  }
  else if (mode_ == 1) {
    reconstruction_mode = 1; // only back substitution
  }
  else {
    throw ExceptionInternal("Kira::run_firefly(): Unkown mode '" + std::to_string(mode_) + "'");
  }

  if (symbols.size() == 0) {
    // Note that the MPI_Bcast() transmitting mpi_symbols below relies
    // on this check, so that it is guaranteed not to be empty.
    throw ExceptionConfigFiles("'run_firefly' can only be run if there variables in the system");
  }

  // Create main director for saved states
  mkdir("firefly_saves", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("firefly_saves/logs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  size_t runs = 1;

  if (!mastersSetZero.empty()) {
    runs = mastersSetZero.size();
  }

  std::vector<std::vector<std::string>> rpn_functions {};
  std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> system_copy {};
  std::vector<firefly::FFInt> funs_copy {};

  for (size_t run = 0; run != runs; ++run) {
    if (runs != 1) {
      logger << "\nReduction " << run + 1 << " / " << runs << "\n";
    }

#ifdef WITH_MPI
    // Send configuration to the workers
    bool cont = true;
    MPI_Bcast(&cont, 1, MPI_CXX_BOOL, firefly::master, MPI_COMM_WORLD);

    // Symbols
    std::string mpi_symbols = "";
    for (const auto& el : symbols) {
      mpi_symbols += el + "!";
    }
    mpi_symbols.pop_back();

    int symbol_size = static_cast<int>(mpi_symbols.size());
    MPI_Bcast(&symbol_size, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
    // Before C++17 the is no non-const string::data() method.
    // Note that it is guaranteed here that mpi_symbols is not empty.
    MPI_Bcast(&mpi_symbols[0], symbol_size, MPI_CHAR, firefly::master,
              MPI_COMM_WORLD);

    // Reconstruction mode
    MPI_Bcast(&reconstruction_mode, 1, MPI_UINT32_T, firefly::master,
              MPI_COMM_WORLD);

    // Mandatory
    int mandatory_size = static_cast<int>(mandatory_vec.size());
    MPI_Bcast(&mandatory_size, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
    MPI_Bcast(mandatory_vec.data(), mandatory_size, mpi_weight_type, firefly::master,
              MPI_COMM_WORLD);
#endif

    bb = new BlackBoxKira(mandatory_vec, reconstruction_mode);

    std::unordered_set<pyred::Weight> masters_set_to_zero;
    std::string file_postfix = "";

    if (mastersSetZero.empty()) {
#ifdef WITH_MPI
      // Zero masters
      int zero_masters_size = 0;
      std::vector<pyred::Weight> tmp {};
      MPI_Bcast(&zero_masters_size, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
      MPI_Bcast(tmp.data(), zero_masters_size, mpi_weight_type, firefly::master,
                MPI_COMM_WORLD);
#endif
      masters_set_to_zero.reserve(0);
    }
    else {
#ifdef WITH_MPI
      // Zero masters
      int zero_masters_size = static_cast<int>(mastersSetZero[run].size());
      MPI_Bcast(&zero_masters_size, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
      MPI_Bcast(mastersSetZero[run].data(), zero_masters_size, mpi_weight_type, firefly::master,
                MPI_COMM_WORLD);
#endif
      file_postfix = "_" + std::to_string(run);
      masters_set_to_zero = std::unordered_set<pyred::Weight> (mastersSetZero[run].begin(), mastersSetZero[run].end());
    }

    if (run == 0) {
#ifdef WITH_MPI
      // 1: normal mode, 2: store system, 0: use stored system
      int load_system = 1;

      if (!mastersSetZero.empty()) {
        load_system = 2;
      }

      MPI_Bcast(&load_system, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
#endif

      load_ff_system(bb, symbols);

      if (mode_ == 1) {
        bb->prepare_backward();
      }

      //TODO make this an option
      if (!mastersSetZero.empty()) {
        system_copy = bb->get_system();
      }
    } else {
#ifdef WITH_MPI
      int load_system = 0;
      MPI_Bcast(&load_system, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
#endif
      bb->set_system(system_copy);
    }

    if (!mastersSetZero.empty()) {
      if (run == 0) {
        auto tmp = bb->post_select(symbols, masters_set_to_zero);
        funs_copy = std::move(tmp.first);
        rpn_functions = std::move(tmp.second);
      } else {
        bb->post_select(symbols, masters_set_to_zero, funs_copy, rpn_functions);
      }
    }
    else {
      bb->force_precompute();
    }

    if (!multiply_factors) {
#ifdef WITH_MPI
      int send_factors = 0;
      MPI_Bcast(&send_factors, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
#endif
    }
    else {
#ifdef WITH_MPI
      int send_factors = 1;
      MPI_Bcast(&send_factors, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
#endif
      prepare_factors(mandatory_vec, masters_set_to_zero);
    }

    // disable the insertion tracer
    pyred::Config::insertion_tracer(0);

    logger << "Reconstructing in: ";
    logger << symbols[0];

    for (auto it = symbols.begin() + 1; it != symbols.end(); ++it) {
      logger << ", " << *it;
    }

    logger << "\n";

    database_reconst = new DataBase(outputDir + "/results/kira.db");
    database_reconst->create_equation_table();

    if (silent_flag == 0) {
      reconst = new firefly::Reconstructor<BlackBoxKira>(
          static_cast<uint32_t>(symbols.size()), static_cast<uint32_t>(coreNumber),
          max_bunch_size, *bb /*, firefly::Reconstructor<BlackBoxKira>::CHATTY*/);
    } else {
      reconst = new firefly::Reconstructor<BlackBoxKira>(
          static_cast<uint32_t>(symbols.size()), static_cast<uint32_t>(coreNumber),
          max_bunch_size, *bb, firefly::Reconstructor<BlackBoxKira>::SILENT);
    }

    //  reconst->set_safe_interpolation();
    if (symbols.size() > 1) {
      reconst->enable_shift_scan();
    }

    if ((symbols.size() > 2 && factor_scan == -1) || factor_scan == 1) {
      reconst->enable_factor_scan();
    }

    bool save = true;
    bool ff_save_exists = false;

    if (save) { // save intermediate results, TODO: use previous pyRed run instead
                // of running it again?
      std::string copied_dir = "firefly_saves/ff_save_" + topology_name + file_postfix;
      std::string save_dir = "ff_save/states";

      if (file_exists(copied_dir.c_str())) {
        logger << "Reconstruction for " << topology_name + file_postfix << " already done\n";
        logger << "Loading results to write them into the database\n";

	// Check whether the ff_save directory exists already and rename it.
        if (file_exists("ff_save")) {
	  char old_file_name[] = "ff_save";
	  char tmp_file_name[] = "ff_save_tmp";
	  std::rename(old_file_name, tmp_file_name);
	  ff_save_exists = true;
	}

        // rename FireFly log file
        char old_log[] = "firefly.log";
        std::string new_log = "firefly_saves/logs/firefly_" + topology_name + file_postfix + ".log";
        std::rename(new_log.c_str(), old_log);

        // rename ff_save
        char old_name[] = "ff_save";
        std::string tmp = "firefly_saves/ff_save_" + topology_name + file_postfix;
        std::rename(tmp.c_str(), old_name);

        // Call the black box once so that equation_lengths is set
        // TODO remove
        std::vector<firefly::FFInt> values;
        for (size_t i = 0; i != symbols.size(); ++i) {
          values.emplace_back(firefly::BaseReconst().get_rand_32());
        }

        (*bb)(values);

        reconst->resume_from_saved_state();
      }
      else if (file_exists(save_dir.c_str())) {
        reconst->resume_from_saved_state();
      }
      else {
        mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::vector<std::string> tags;

        std::vector<firefly::FFInt> values;
        for (size_t i = 0; i != symbols.size(); ++i) {
          values.emplace_back(firefly::BaseReconst().get_rand_32());
        }

        (*bb)(values);

        std::string filename = "ff_save/tags";
        std::ofstream filestream;
        filestream.open(filename.c_str());

        for (const auto& ids : bb->assignment) {
          std::string tag =
              ids.first.to_string() + "_" + ids.second.to_string();
          tags.emplace_back(tag);
          filestream << tag << "\n";
        }

        filestream.close();

        reconst->set_tags(tags);
      }
    }

    bool done = false;
    std::unordered_set<pyred::Weight> masters_set {};

    ThreadPool tp;
    tp.initialize(1);
    tp.enqueue([this, &done]() {
      reconst->reconstruct();

      std::lock_guard<std::mutex> lock(reconst_done);
      done = true;
      cond_reconst_done.notify_one();
    });

    {
      std::unique_lock<std::mutex> lock(reconst_done);

      if (!done) {
        while (cond_reconst_done.wait_for(lock, std::chrono::minutes(10)) ==
               std::cv_status::timeout) {
          write_to_database(masters_set);
        }
      }
    }

    write_to_database(masters_set);

    logger << "\nAverage times for probe solution steps:\n";
    logger << "Coefficient evaluation: " << std::setprecision(6) << bb->parse_average << " s\n";
    logger << "Forward elimination:    " << bb->forward_average << " s\n";
    logger << "Back substitution:      " << bb->back_average << " s\n";
    logger << std::setprecision(1);

    std::vector<pyred::Weight> masters_vec;
    masters_vec.reserve(masters_set.size());

    for (const auto& master : masters_set) {
      masters_vec.emplace_back(master);
    }

    std::sort(masters_vec.begin(), masters_vec.end());

    std::string masters_str = masters_list_str(masters_vec);
    std::ofstream masters_file(
        (outputDir + "/results/" + integralfamily.name + "/masters.final").c_str());

    logger << "\n*****Master integrals***********************************\n";
    logger << masters_str;
    masters_file << masters_str;
    logger << "\nTotal number of master integrals: " << masters_vec.size() << "\n\n";

    masters_file.close();

    firefly::RatReconst::reset();

    // reset the insertion tracer
    pyred::Config::insertion_tracer(1);

    // rename FireFly log file
    char old_log[] = "firefly.log";
    std::string new_log = "firefly_saves/logs/firefly_" + topology_name + file_postfix + ".log";
    std::rename(old_log, new_log.c_str());

    if (save) {
      char old_name[] = "ff_save";
      std::string tmp = "firefly_saves/ff_save_" + topology_name + file_postfix;

      std::rename(old_name, tmp.c_str());

      if (ff_save_exists) {
        char tmp_file_name[] = "ff_save_tmp";
        std::rename(tmp_file_name, old_name);
      }
    }

    first = true;
    delete database_reconst;
    delete bb;
    delete reconst;
  }

  integral_data.clear();
}

template <typename BlackBoxType>
void Kira::load_ff_system(BlackBoxType* bb_, const std::vector<std::string>& variables) {
  std::string inputName;

  if (reconstruction_mode == 0) {
    inputName =
        (outputDir + "/tmp/" + integralfamily.name + "/SYSTEMconfig").c_str();
  }
  else if (reconstruction_mode == 1) {
    inputName =
        (outputDir + "/tmp/" + integralfamily.name + "/VERconfig").c_str();
  }
  else {
    throw ExceptionInternal("Kira::run_firefly(): Unkown mode '" + std::to_string(reconstruction_mode) + "'");
  }

  std::ifstream inputConfig;
  std::size_t number_of_equations = 0;

  if (file_exists(inputName.c_str())) {
    inputConfig.open(inputName.c_str());
    inputConfig >> number_of_equations;
    inputConfig.close();
  }
  else {
    throw ExceptionInternal("Kira::run_firefly(): Could not read the number of equations from '" + inputName + "'");
  }

  logger << "Loading " << number_of_equations << " equations.\n\n";

  bb_->reserve(number_of_equations);

  int eqLength;
  pyred::Weight ID;
  size_t parsed_functions = 0;
  size_t parsed_equations = 0;
  size_t curr_percentage = 1;
  std::string line;

  std::vector<std::string> add;
  add.push_back("");
  add.push_back("sym");

  firefly::ShuntingYardParser parser(std::vector<std::string>{}, variables, true);
  auto time0 = std::chrono::high_resolution_clock::now();

  for (unsigned k = 0; k != collectReductions.size(); ++k) {
    for (int sectorNumber = 0; sectorNumber != (1 << integralfamily.jule);
         ++sectorNumber) {

      for (int ll = 0; ll != 2; ++ll) {
        /*backward compability with Kira versions prior to 2.1*/
        int itF = -1;

        while (true) {
          if (reconstruction_mode == 0) {
            if (itF == -1) {
              inputName = outputDir + "/tmp/" + integralfamily.name + "/SYSTEM" +
                          add[ll] + "_" + collectReductions[k] + "_" +
                          something_string(sectorNumber) + ".gz";
            }
            else {
              inputName = outputDir + "/tmp/" + integralfamily.name + "/SYSTEM" + add[ll] + "_" + collectReductions[k] + "_" + something_string(sectorNumber)+ "_" + something_string(itF) + ".gz";
            }
          }
          else if (reconstruction_mode == 1) {
            if (itF == -1) {
              inputName = outputDir + "/tmp/" + integralfamily.name + "/VER" +
                          add[ll] + "_" + collectReductions[k] + "_" +
                          something_string(sectorNumber) + ".gz";
            }
            else {
              inputName = outputDir + "/tmp/" + integralfamily.name + "/VER" +
                          add[ll] + "_" + collectReductions[k] + "_" +
                          something_string(sectorNumber) + "_" + something_string(itF) +  ".gz";
            }
          }

          if (file_exists(inputName.c_str())) {
            bool unexpected_end = true;
            GZSTREAM_NAMESPACE::igzstream input;
            input.open(inputName.c_str());

            while (true) {
              std::vector<std::pair<pyred::Weight, std::size_t>> tmp_eq;
#ifdef WITH_MPI
              std::string eqn_s = "";
#endif

              // the loop has to exit here, otherwise the file is corrupted
              unexpected_end = false;
              if (!(input >> line)) break;

              unexpected_end = true;
              if (!(input >> ID)) break;
              if (!(input >> eqLength)) break;

              for (int j = 0; j != eqLength; ++j) {
                int skip;
                std::string coef;
                pyred::Weight id;
                int sector;
                int flag2;

                if (!(input >> skip)) break;
                if (!(input >> coef)) break;
                if (!(input >> id)) break;
                if (!(input >> sector)) break; // sector
                if (!(input >> skip)) break;   // topology
                if (!(input >> flag2)) break;  // flag2

                if (reconstruction_mode == 1 && j == 0) {
                  coef += "*(-1)";
                }

                // TODO Horner form?
#ifdef WITH_MPI
                eqn_s += coef + "!" + id.to_string() + "|";
#endif

                std::size_t pos = parser.add_otf(coef, true);
                ++parsed_functions;
                tmp_eq.emplace_back(std::make_pair(id, pos));

                integral_data.emplace(
                    std::make_pair(id, std::make_pair(sector, flag2)));
              }

              if (parsed_equations + 1 < number_of_equations + 1 && parsed_equations + 1 > curr_percentage*number_of_equations/10) {
                ++curr_percentage;
                if (silent_flag == 0) {
                  std::cerr << "\033[1;34mFireFly info:\033[0m " << parsed_equations + 1 << " / " << number_of_equations << "\r";
                }
              }
              ++parsed_equations;

#ifdef WITH_MPI
              eqn_s.pop_back();

              int batch_size = 2147483645; // 2147483647 is the largest signed 32-bit integer
              int split = static_cast<int>(eqn_s.size()) / batch_size;

              if (split > 1) {
                int tmp;
                if (static_cast<int>(eqn_s.size()) % batch_size == 0) {
                  tmp = -split;
                }
                else {
                  tmp = -1 - split;
                }
                MPI_Bcast(&tmp, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
              }
              else if (split == 1 &&
                        static_cast<int>(eqn_s.size()) != batch_size) {
                int tmp = -1 - split;
                MPI_Bcast(&tmp, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
              }

              for (int i = 0; i != 1 + split; ++i) {
                if (i == split) {
                  int amount = static_cast<int>(eqn_s.size()) - i * batch_size;
                  if (amount != 0) {
                    MPI_Bcast(&amount, 1, MPI_INT, firefly::master,
                              MPI_COMM_WORLD);
                    MPI_Bcast(&eqn_s[i * batch_size], amount, MPI_CHAR,
                              firefly::master, MPI_COMM_WORLD);
                  }
                }
                else {
                  int amount = batch_size;
                  MPI_Bcast(&amount, 1, MPI_INT, firefly::master,
                            MPI_COMM_WORLD);
                  MPI_Bcast(&eqn_s[i * batch_size], amount, MPI_CHAR,
                            firefly::master, MPI_COMM_WORLD);
                }
              }
#endif
              bb_->add_eqn(tmp_eq);
            }

            input.close();

            if (unexpected_end) {
              throw ExceptionInternal("Kira::run_firefly(): File '" + inputName + "' corrupted");
            }

            ++itF;
          }
          else if (itF == 0)
            ++itF;
          else
            break;
        }
      }
    }
  }

  bb_->sort_system();

  auto time1 = std::chrono::high_resolution_clock::now();
  std::cout << "\033[1;34mFireFly info:\033[0m " << number_of_equations << " / " << number_of_equations << "\n";
  logger << "\033[1;34mFireFly info:\033[0m Parsed " << parsed_functions
      << " functions of which " << parser.get_size() << " are distinct\n";

#ifdef WITH_MPI
  int end = -1;
  MPI_Bcast(&end, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  bb_->set_parser(parser);

  logger << "\nSystem loaded in " << std::to_string(std::chrono::duration<double>(time1 - time0).count()) << " s\n";
}

void Kira::prepare_factors(const std::vector<pyred::Weight>& mandatory_vec, const std::unordered_set<pyred::Weight>& masters_set_to_zero) {
#ifdef WITH_MPI
  MPI_Datatype &mpi_weight_type = MPIHelper::weight_type();
#endif
  logger << "Preparing factors\n";

  std::unordered_map<pyred::Weight, std::unordered_map<pyred::Weight, std::size_t>> system {};
  system.reserve(mandatory_vec.size());
  firefly::ShuntingYardParser parser(std::vector<std::string>{}, symbols, true);

  for (const auto &id : mandatory_vec) {
    auto it = prefactorEquations.find(id);

    if (it != prefactorEquations.end()) {
      std::unordered_map<pyred::Weight, std::size_t> equation {};
      equation.reserve(it->second.size());

#ifdef WITH_MPI
      int next = 1;
      MPI_Bcast(&next, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

      pyred::Weight eqid = it->first;
      MPI_Bcast(&eqid, 1, mpi_weight_type, firefly::master, MPI_COMM_WORLD);
#endif

      for (const auto& integral : it->second) {
        if (integral.first == it->first) {
          continue;
        }
        else if (masters_set_to_zero.find(integral.first) == masters_set_to_zero.end()) {
          // copy required due to SYP deleting string at Initialization
          std::string coef = integral.second;

#ifdef WITH_MPI
          int next = 2;
          MPI_Bcast(&next, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

          pyred::Weight iid = integral.first;
          MPI_Bcast(&iid, 1, mpi_weight_type, firefly::master, MPI_COMM_WORLD);

          int batch_size =
              2147483645; // 2147483647 is the largest signed 32-bit integer
          int split = static_cast<int>(coef.size()) / batch_size;

          if (split > 1) {
            int tmp;
            if (static_cast<int>(coef.size()) % batch_size == 0) {
              tmp = -split;
            }
            else {
              tmp = -1 - split;
            }
            MPI_Bcast(&tmp, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
          }
          else if (split == 1 &&
                   static_cast<int>(coef.size()) != batch_size) {
            int tmp = -1 - split;
            MPI_Bcast(&tmp, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
          }

          for (int i = 0; i != 1 + split; ++i) {
            if (i == split) {
              int amount =
                  static_cast<int>(coef.size()) - i * batch_size;
              if (amount != 0) {
                MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
                MPI_Bcast(&coef[i * batch_size], amount, MPI_CHAR,
                          firefly::master, MPI_COMM_WORLD);
              }
            }
            else {
              int amount = batch_size;
              MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
              MPI_Bcast(&coef[i * batch_size], amount, MPI_CHAR,
                        firefly::master, MPI_COMM_WORLD);
            }
          }
#endif

          std::size_t pos = parser.add_otf(coef, true);
          equation.emplace(std::make_pair(integral.first, pos));
        }
      }

      system.emplace(std::make_pair(id, std::move(equation)));
    }
  }

#ifdef WITH_MPI
  int end = -1;
  MPI_Bcast(&end, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
#endif

  bb->set_factors(system, parser);

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Kira::write_to_database(std::unordered_set<pyred::Weight>& masters) {
  std::vector<std::pair<std::string, firefly::RationalFunction>> results =
      reconst->get_early_results();

  if (!results.empty()) {
    database_reconst->begin_transaction();
    database_reconst->prepare_backsubstitution();

    if (first) {
      std::unique_lock<std::mutex> lock(bb->mut);

      if (!(bb->equation_lengths.empty())) {
        lock.unlock();

        for (const auto& equation : bb->equation_lengths) {
          BaseIntegral integral;
          integral.id = equation.first;
          integral.coefficientString = "1";
          auto values = integral_data.at(equation.first);
          integral.flag2 = values.second;
          integral.characteristics[SECTOR] = values.first;

          database_reconst->bind_equation(integral, equation.second,
                                          equation.first, mastersReMap);
        }

        first = false;
      }
    }

    for (const auto& res : results) {
      std::pair<pyred::Weight, pyred::Weight> ids_eq_int = std::make_pair(
          pyred::Weight(res.first.substr(0, res.first.find("_"))),
          pyred::Weight(res.first.substr(res.first.find("_") + 1)));

      BaseIntegral integral;
      integral.id = ids_eq_int.second;
      auto values = integral_data.at(ids_eq_int.second);
      integral.flag2 = values.second;
      integral.characteristics[SECTOR] = values.first;

      masters.insert(ids_eq_int.second);

      std::string factor = "";

      if (multiply_factors) {
        auto it = prefactorEquations.find(ids_eq_int.first);

        if (it != prefactorEquations.end()) {
          auto itt = it->second.find(ids_eq_int.second);

          if (itt != it->second.end()) {
            factor = "*prefactor[" + itt->second + "]";
          }
        }
      }

      integral.coefficientString =  res.second.generate_horner(symbols) + factor;

      database_reconst->bind_equation(integral,
                                      bb->equation_lengths.at(ids_eq_int.first),
                                      ids_eq_int.first, mastersReMap);
    }

    database_reconst->commit_transaction();
  }
}

void Kira::run_numerical(const std::vector<pyred::Weight>& mandatory_vec, const std::string& topology_name, const std::string& input_file_name) {
  bbNum = new BlackBoxNumerical(mandatory_vec);
  std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>> results {};
  std::vector<std::string> variables = std::get<1>(numericalPointsAux.front());

  // disable the insertion tracer
  pyred::Config::insertion_tracer(0);

  load_ff_system(bbNum, variables);
  bbNum->force_precompute();

  ThreadPool tp;
  tp.initialize(static_cast<uint32_t>(coreNumber));

  for (const auto & set_of_points : numericalPointsAux) {
    if (std::get<1>(set_of_points) != variables) {
      throw ExceptionConfigFiles("'numerical_points': the order of variables in '" + input_file_name + "' has to be the same for all sets of sample points");
    }

    firefly::FFInt::set_new_prime(std::stoul(std::get<0>(set_of_points)));
    bbNum->prime_changed();
    logger << "Prime set to " << firefly::FFInt::p << "\n";

    if (max_bunch_size == 1) {
      for (std::size_t i = 0; i != std::get<2>(set_of_points).size(); ++i) {
        tp.enqueue([this, &results, &set_of_points, i]() {
          numerical_compute(results, set_of_points, i);
        });
      }
    }
    else {
      std::size_t next_bunch_size = firefly::compute_bunch_size(std::get<2>(set_of_points).size(), static_cast<uint32_t>(coreNumber), max_bunch_size);
      std::size_t i = 0;

      while(i != std::get<2>(set_of_points).size()) {
        tp.enqueue([this, &results, &set_of_points, i, next_bunch_size]() {
          switch(next_bunch_size) {
            case 1:
              numerical_compute(results, set_of_points, i);
              break;
            case 2:
              numerical_compute_bunch<2>(results, set_of_points, i);
              break;
            case 4:
              numerical_compute_bunch<4>(results, set_of_points, i);
              break;
            case 8:
              numerical_compute_bunch<8>(results, set_of_points, i);
              break;
            case 16:
              numerical_compute_bunch<16>(results, set_of_points, i);
              break;
            case 32:
              numerical_compute_bunch<32>(results, set_of_points, i);
              break;
            case 64:
              numerical_compute_bunch<64>(results, set_of_points, i);
              break;
            case 128:
              numerical_compute_bunch<128>(results, set_of_points, i);
              break;
          }
        });

        i += next_bunch_size;
        next_bunch_size = firefly::compute_bunch_size(std::get<2>(set_of_points).size() - i, static_cast<uint32_t>(coreNumber), max_bunch_size);
      }
    }

    numerical_write_results(results, mandatory_vec, topology_name, input_file_name, std::get<2>(set_of_points).size());

    logger << "Finished prime\n" << "Average times:\n";
    logger << "Parsing: " << std::setprecision(3) << bbNum->parse_average << " s\n";
    logger << "Forward: " << bbNum->forward_average << " s\n";
    logger << "Backsubs: " << bbNum->back_average << " s\n\n";
    logger << std::setprecision(1);
  }

  // TODO write masters.final?

  // reset the insertion tracer
  pyred::Config::insertion_tracer(1);

  delete bbNum;
}

// TODO move most of the functionality to the black box?

void Kira::numerical_compute(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::tuple<std::string, std::vector<std::string>, std::vector<std::vector<std::string>>>& set_of_points, std::size_t pos) {
  // TODO the shunting-yard-parser is an overkill, but probably irrelevant for the final performance ...
  // TODO we need a silent mode ...
  firefly::ShuntingYardParser parser = firefly::ShuntingYardParser(std::get<2>(set_of_points)[pos], {}, false, false, true);
  std::vector<firefly::FFInt> values = parser.evaluate_pre(std::vector<firefly::FFInt> {});
  pyred::SystemOfEqs<firefly::FFInt> result = (*bbNum)(values);

  std::lock_guard<std::mutex> lock(num_done);
  results.emplace(std::move(values), std::move(result));
  cond_num_done.notify_one();
}

template <std::size_t N>
void Kira::numerical_compute_bunch(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::tuple<std::string, std::vector<std::string>, std::vector<std::vector<std::string>>>& set_of_points, std::size_t pos) {
  // TODO the shunting-yard-parser is an overkill, but probably irrelevant for the final performance ...
  // TODO we need a silent mode ...
  firefly::ShuntingYardParser parser;
  std::array<std::vector<firefly::FFInt>, N> values_arr;
  std::vector<firefly::FFIntVec<N>> bunch_values_vec(std::get<2>(set_of_points)[pos].size());

  for (std::size_t i = 0; i != N; ++i) {
    parser = firefly::ShuntingYardParser(std::get<2>(set_of_points)[pos + i], {}, false, false, true);
    values_arr[i] = parser.evaluate_pre(std::vector<firefly::FFInt> {});

    for (std::size_t j = 0; j != std::get<2>(set_of_points)[pos + i].size(); ++j) {
      bunch_values_vec[j][i] = values_arr[i][j];
    }
  }

  pyred::SystemOfEqs<firefly::FFIntVec<N>> bunch_result = (*bbNum)(bunch_values_vec);

  // TODO this whole conversion back to the scalar format looks like a huge waste of time, but we have to got to the drawing board about the format

  std::array<pyred::SystemOfEqs<firefly::FFInt>, N> result_arr;

  for (auto it = bunch_result.sys.begin(); it != bunch_result.sys.end(); ++it) {
    // TODO array would be optimal, but pyred::Equation has no default constructor
    std::vector<pyred::Equation<firefly::FFInt>> equation_vec(N, it->eqnum);

    for (auto itt = it->begin(); itt != it->end(); ++itt) {
      for (std::size_t i = 0; i != N; ++i) {
        equation_vec[i].eq.emplace_back(std::make_pair(itt->first, itt->second[i]));
      }
    }

    for (std::size_t i = 0; i != N; ++i) {
      equation_vec[i].eq.shrink_to_fit();
      result_arr[i].sys.emplace_back(std::move(equation_vec[i]));
    }
  }

  std::lock_guard<std::mutex> lock(num_done);
  for (std::size_t i = 0; i != N; ++i) {
    results.emplace(std::move(values_arr[i]), std::move(result_arr[i]));
  }
  cond_num_done.notify_one();
}

// TODO record trivial integrals -> no, write them already out during the generation -> Johann

void Kira::numerical_write_results(std::queue<std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>>>& results, const std::vector<pyred::Weight>& mandatory_vec, const std::string& topology_name, const std::string& input_file_name, const std::size_t points) {
  std::pair<std::vector<firefly::FFInt>, pyred::SystemOfEqs<firefly::FFInt>> result;
  std::string filename;
  std::ofstream filestream;

  std::unique_lock<std::mutex> lock(num_done);

  for (std::size_t counter = 0; counter != points; ++counter) {
    if (counter % pointsPerFile == 0) {
      // TODO also put in filename of selected integrals
      filename = inputDir + "results/" + topology_name + "/" + input_file_name + "_" + std::to_string(firefly::FFInt::p) + "_" + std::to_string(counter/pointsPerFile) + ".m";
      filestream.open(filename.c_str());
      filestream << "{" << firefly::FFInt::p << ",\n";
    }

    if (results.empty()) {
      cond_num_done.wait(lock);
    }

    result = std::move(results.front());
    results.pop();

    lock.unlock();

    filestream << "{{";

    for (auto it = result.first.begin(); it != result.first.end(); ++it) {
      filestream << *it;

      if (it != result.first.end() - 1) {
        filestream << ", ";
      }
    }

    filestream << "},\n{";

    if (!mandatory_vec.empty()) {
      auto it_mandatory = mandatory_vec.begin();

      for (auto it = result.second.sys.begin(); it != result.second.sys.end(); ++it) {
        while (it->front().first > *it_mandatory) {
          ++it_mandatory;
        }

        if (it->front().first == *it_mandatory) {
          numerical_write_equation(it, filestream);

          ++it_mandatory;

          if (it != result.second.sys.end() - 1) {
            filestream << ",\n";
          }

          if (it_mandatory == mandatory_vec.end()) {
            break;
          }
        }
      }
    }
    else {
      for (auto it = result.second.sys.begin(); it != result.second.sys.end(); ++it) {
        numerical_write_equation(it, filestream);

        if (it != result.second.sys.end() - 1) {
          filestream << ",\n";
        }
      }
    }

    if ((counter + 1) % pointsPerFile != 0 && counter + 1 != points) {
      filestream << "}},\n";
    }
    else {
      filestream << "}}}\n";
      filestream.close();
    }

    lock.lock();
  }
}

void Kira::numerical_write_equation(const std::vector<pyred::Equation<firefly::FFInt>>::iterator it, std::ofstream& filestream) {
  filestream << pyred::Integral(it->eq.begin()->first).to_string() << " ->\n";

  if (it->eq.size() > 1) {
    for (auto itt = ++(it->eq.begin()); itt != it->eq.end(); ++itt) {
      filestream << " + " << pyred::Integral(itt->first).to_string() << "*" << std::to_string(itt->second.n) << "\n";
    }
  }
  else {
    filestream << "0\n";
  }
}
