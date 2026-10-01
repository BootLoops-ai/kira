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

#include <firefly/config.hpp> // for WITH_MPI

#ifdef WITH_MPI

#include "kira/mpi_work.h"

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <mpi.h>

#include <firefly/MPIWorker.hpp>
#include <firefly/RatReconst.hpp>
#include <firefly/ShuntingYardParser.hpp>

#include <pyred/defs.h>
#include <kira/black_box.h>
#include "kira/exceptions.h"

void mpi_work(const uint32_t threads, const uint32_t bunch_size) {
  std::vector<std::vector<std::string>> rpn_functions {};
  std::vector<std::vector<std::pair<pyred::Weight, std::size_t>>> system_copy {};
  std::vector<firefly::FFInt> funs_copy {};
  MPI_Datatype &mpi_weight_type = MPIHelper::weight_type();

  while (true) {
    bool cont;
    MPI_Bcast(&cont, 1, MPI_CXX_BOOL, firefly::master, MPI_COMM_WORLD);

    if (!cont) {
      break;
    }

    // Receive configuration

    // Symbols
    uint32_t n_vars;
    std::vector<std::string> symbols;
    int amount;
    MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

    char* symbols_c = new char[amount];
    MPI_Bcast(symbols_c, amount, MPI_CHAR, firefly::master, MPI_COMM_WORLD);

    std::string var = "";
    for (int i = 0; i != amount; ++i) {
      char el = symbols_c[i];
      if (el != '!')
        var += el;
      else {
        symbols.emplace_back(var);
        var = "";
      }
    }

    symbols.emplace_back(var);
    n_vars = symbols.size();

    delete[] symbols_c;

    // Reconstruction mode
    uint32_t mode;
    MPI_Bcast(&mode, 1, MPI_UINT32_T, firefly::master, MPI_COMM_WORLD);

    // Mandatory
    MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
    std::vector<pyred::Weight> mandatory_vec;
    mandatory_vec.resize(amount);
    MPI_Bcast(mandatory_vec.data(), amount, mpi_weight_type, firefly::master,
              MPI_COMM_WORLD);

    // Black box and system of equations
    BlackBoxKira bb(mandatory_vec, mode);

    // Zero masters
    MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
    std::vector<pyred::Weight> masters_set_to_zero_vec;
    masters_set_to_zero_vec.resize(amount);
    MPI_Bcast(masters_set_to_zero_vec.data(), amount, mpi_weight_type, firefly::master,
              MPI_COMM_WORLD);

    std::unordered_set<pyred::Weight> masters_set_to_zero(masters_set_to_zero_vec.begin(), masters_set_to_zero_vec.end());
    masters_set_to_zero_vec.clear();
    masters_set_to_zero_vec.resize(0);

    // 1: normal mode, 2: store system, 0: use stored system
    int load_system;
    MPI_Bcast(&load_system, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

    if (load_system) {
      // Receive system of equations
      firefly::ShuntingYardParser parser(std::vector<std::string>{}, symbols,
                                         true);
      std::size_t position = 0;

      while (true) {
        int amount;
        int multiple = 1;
        MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

        if (amount == -1) {
          break;
        }
        else if (amount < -1) {
          multiple = -amount;
        }

        std::string tmp_s = "";
        std::vector<std::pair<pyred::Weight, std::size_t>> tmp_eqn;

        for (int j = 0; j != multiple; ++j) {
          if (multiple != 1) {
            MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
          }

          char* eqn_c = new char[amount];
          MPI_Bcast(eqn_c, amount, MPI_CHAR, firefly::master, MPI_COMM_WORLD);

          for (int i = 0; i != amount; ++i) {
            char c = eqn_c[i];
            switch (c) {
              case '!': {
                position = parser.add_otf(tmp_s, true);
                tmp_s = "";
                break;
              }
              case '|': {
                tmp_eqn.emplace_back(
                    std::make_pair(pyred::Weight(tmp_s), position));
                tmp_s = "";
                break;
              }
              default: {
                tmp_s += c;
                break;
              }
            }
          }

          delete[] eqn_c;
        }

        tmp_eqn.emplace_back(std::make_pair(pyred::Weight(tmp_s), position));
        bb.add_eqn(tmp_eqn);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      bb.set_parser(parser);

      if (mode == 1) {
        bb.prepare_backward();
      }

      if (load_system == 2) {
        system_copy = bb.get_system();
      }
    }
    else {
      bb.set_system(system_copy);
    }

    if (load_system == 2) {
      auto tmp = bb.post_select(symbols, masters_set_to_zero);
      funs_copy = std::move(tmp.first);
      rpn_functions = std::move(tmp.second);
    }
    else if (load_system == 0) {
      bb.post_select(symbols, masters_set_to_zero, funs_copy, rpn_functions);
    }
    else {
      bb.force_precompute();
    }

    int receive_factors;
    MPI_Bcast(&receive_factors, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

    if (receive_factors) {
      // Receive factors
      std::unordered_map<pyred::Weight, std::unordered_map<pyred::Weight, std::size_t>> system {};
      system.reserve(mandatory_vec.size());
      firefly::ShuntingYardParser parser(std::vector<std::string>{}, symbols, true);
      pyred::Weight equation_id;
      std::unordered_map<pyred::Weight, std::size_t> equation {};

      while (true) {
        int next;
        MPI_Bcast(&next, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

        if (next == -1) {
          if (!equation.empty()) {
            system.emplace(std::make_pair(equation_id, std::move(equation)));
            equation = std::unordered_map<pyred::Weight, std::size_t> {};
          }

          break;
        }
        else if (next == 1) {
          if (!equation.empty()) {
            system.emplace(std::make_pair(equation_id, std::move(equation)));
            equation = std::unordered_map<pyred::Weight, std::size_t> {};
          }

          MPI_Bcast(&equation_id, 1, mpi_weight_type, firefly::master, MPI_COMM_WORLD);
        }
        else if (next == 2) {
          pyred::Weight integral_id;
          MPI_Bcast(&integral_id, 1, mpi_weight_type, firefly::master, MPI_COMM_WORLD);

          int amount;
          int multiple = 1;
          MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);

          if (amount < -1) {
            multiple = -amount;
          }

          std::string tmp_s = "";

          for (int j = 0; j != multiple; ++j) {
            if (multiple != 1) {
              MPI_Bcast(&amount, 1, MPI_INT, firefly::master, MPI_COMM_WORLD);
            }

            char* tmp_c = new char[amount];
            MPI_Bcast(tmp_c, amount, MPI_CHAR, firefly::master, MPI_COMM_WORLD);

            tmp_s += std::string(tmp_c, amount);

            delete[] tmp_c;
          }

          std::size_t pos = parser.add_otf(tmp_s, true);
          equation.emplace(std::make_pair(integral_id, pos));
        }
        else {
          throw ExceptionInternal("mpi_work(): Unkown value received while preparing factors");
        }
      }

      bb.set_factors(system, parser);

      MPI_Barrier(MPI_COMM_WORLD);
    }

    firefly::MPIWorker<BlackBoxKira>(n_vars, threads, bunch_size, bb);

    firefly::RatReconst::reset();
  }
}

#endif // WITH_MPI
