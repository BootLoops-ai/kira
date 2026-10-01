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

#ifndef CREATESEEDS_H
#define CREATESEEDS_H
#include <tuple>

#include "kira/kira.h"

class SeedsKira {
public:
  SeedsKira(Kira& kira, int nOfTopology);
  ~SeedsKira();
  void gen_int(int n, int** array, int maxI, int option);
  void fusion(int j, int iS2, int topology, int noIBP, int symDOTS,
              int symNUMS);
  void prepare_seeds();
  void generate_dots();
  void generate_nums();
  void create_seeds();
  void create_symmetry_seeds();
  int save();
  Kira* kira;
  int den_max;
  int num_max, num_max_sym;
  int nOfTopology;
  std::vector<std::uint64_t>** aIBPint;
  int ***dotsMaster, ***scalMaster;
  int *nSector, **mappingSector, **mappingSectorReverse;
  std::set<std::vector<int> > testSet;
  std::unordered_map<std::uint64_t, std::tuple<int, int, int> >
      integralMap; // sector, topology, flag2
};

void initiate_seeds(Kira& kira, std::size_t nOfTopology);
void run_relations(Kira& kira_, std::uint32_t nOfTopology_);
void get_properties(std::uint64_t id,
                    std::tuple<std::string, unsigned, unsigned>& integral);
#endif
