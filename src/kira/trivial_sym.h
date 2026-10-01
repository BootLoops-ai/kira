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

#ifndef TRIVIAL_SYM_H_
#define TRIVIAL_SYM_H_
#include <sys/stat.h>

#include <algorithm>
#include <tuple>

#include "ginac/ginac.h"
#include "pyred/defs.h"
#include "kira/ReadYamlFiles.h"
#include "kira/kira.h"
#include "kira/tools.h"

class Kira;

class TermsGiNaC {
public:
  TermsGiNaC(int size_) { size = size_; }
  ~TermsGiNaC() {}

  //   friend std::ostream& operator<< (std::ostream& out, const TermsGiNaC&
  //   term);
  std::string coef;
  uint64_t coefN;
  std::vector<int> elem;
  int size;
  int powerSum;
  int cutPropagator;
  std::vector<std::pair<pyred::Weight, std::string> > equation;
  GiNaC::ex coefEx;
};

class Terms {
public:
  Terms(int size_) { size = size_; }
  ~Terms() {}

  friend std::ostream& operator<<(std::ostream& out, const Terms& term);
  std::string coef;
  uint64_t coefN;
  std::vector<int> elem;
  int size;
  int powerSum;
  int cutPropagator;
  std::vector<std::pair<pyred::Weight, std::string> > equation;
};

// pyred::Coeff_int::set_prime(9223372036854775783) should be used before to
// ensure that users did not set a weird prime
class Pak {
public:
  Pak(){};

  ~Pak(){
    if(canonica.size()==0){
      for(auto itCanonica: canonicaX){
          delete itCanonica;
      }
    }
    for(auto itCanonica: canonica){
      for(auto it: itCanonica)
        delete it;
    }
  };

  // symmetry
  Pak(GiNaC::ex listOfTerms, uint32_t jule, std::vector<int>& seed,
      std::vector<int>& holes,
      std::vector<std::tuple<GiNaC::lst, GiNaC::lst, int, GiNaC::ex,
                             std::vector<std::string> > >& externalTransf,
      GiNaC::lst& invariantsReplacementRev, unsigned maskCut);

  std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> >
  generate_combinatorics(
      std::vector<int>& seed, std::vector<int>& holes,
      std::vector<std::tuple<GiNaC::lst, GiNaC::lst, int, GiNaC::ex,
                             std::vector<std::string> > >& externalTransf,
      std::vector<std::string>& invarSol);

  void collect(std::vector<Terms*>& canonicaX);

  // trivial
  // needs pyred configured:
  // pyred::Config::verbosity(0);
  // pyred::Config::insertion_tracer(1);
  Pak(GiNaC::ex listOfTerms, uint32_t jule, std::vector<int>& seed,
      std::vector<int>& holes);

  void check_trivial(std::vector<int>& seed, std::vector<int>& holes,
                     int& check);

  void collect_trivial(std::vector<TermsGiNaC*>&);

  //   bool compare(Pak& pak);

  int l_Indices;
  uint32_t jule;
  std::vector<Terms*> canonicaX;
  std::vector<std::vector<Terms*> > canonica;
};

struct compMatrix {
  compMatrix(int level_) { level = level_; }

  bool operator()(const Terms* s1, const Terms* s2) const {
    if (level == 0) {
      if(s1->cutPropagator < s2->cutPropagator)
        return true;
      if(s1->cutPropagator > s2->cutPropagator)
        return false;

      if(s1->powerSum < s2->powerSum)
        return true;
      if(s1->powerSum > s2->powerSum)
        return false;

      return
        s1->coefN < s2->coefN;
    }
    else {
      if(s1->cutPropagator < s2->cutPropagator)
        return true;
      if(s1->cutPropagator > s2->cutPropagator)
        return false;

      if(s1->powerSum < s2->powerSum)
        return true;
      if(s1->powerSum > s2->powerSum)
        return false;

      if (s1->coefN < s2->coefN)
        return true;
      else if (s1->coefN > s2->coefN)
        return false;

      for (int i = 0; i < level; i++) {
        if (s1->elem[i] < s2->elem[i]) return true;
        if (s1->elem[i] > s2->elem[i]) return false;
      }
      return false;
    }
    return false;
  }

  int level;
};

struct compMatrix2 {
  compMatrix2(int level_, int row_) {
    level = level_;
    row = row_;
  }

  bool operator()(const Terms* s1, const Terms* s2) const {
    if (level == 0) {
      if(s1->cutPropagator< s2->cutPropagator)
        return true;
      if(s1->cutPropagator > s2->cutPropagator)
        return false;

      if(s1->powerSum < s2->powerSum)
        return true;
      if(s1->powerSum > s2->powerSum)
        return false;

      return
        s1->coefN < s2->coefN;
    }
    else {
      if(s1->cutPropagator< s2->cutPropagator)
        return true;
      if(s1->cutPropagator > s2->cutPropagator)
        return false;

      if(s1->powerSum < s2->powerSum)
        return true;
      if(s1->powerSum > s2->powerSum)
        return false;

      if (s1->coefN < s2->coefN)
        return true;
      else if (s1->coefN > s2->coefN)
        return false;

      for (int i = 0; i < level - 1; i++) {
        if (s1->elem[i] < s2->elem[i]) return true;
        if (s1->elem[i] > s2->elem[i]) return false;
      }

      if (s1->elem[row] < s2->elem[row]) return true;
      if (s1->elem[row] > s2->elem[row]) return false;
    }
    return false;
  }

  int row;
  int level;
};

struct compMatrixV {
  compMatrixV(int level_) { level = level_; }

  bool operator()(
      const std::tuple<std::vector<Terms*>, std::vector<int> >& s1,
      const std::tuple<std::vector<Terms*>, std::vector<int> >& s2) const {
    for (size_t i = 0; i < std::get<0>(s1).size(); i++) {
      if (std::get<0>(s1)[i]->elem[level - 1] <
          std::get<0>(s2)[i]->elem[level - 1])
        return true;
      if (std::get<0>(s1)[i]->elem[level - 1] >
          std::get<0>(s2)[i]->elem[level - 1])
        return false;
    }
    return false;
  }

  int level;
};

struct equalMatrix {
  equalMatrix(int level_) { level = level_; }

  bool operator()(const Terms* s1, const Terms* s2) const {
    if (s1->elem[level - 1] == s2->elem[level - 1]) return true;
    return false;
  }

  int level;
};

struct equalAllMatrix {
  equalAllMatrix(int size_) { size = size_; }

  bool operator()(const Terms* s1, const Terms* s2) const {
    if (s1->coefN == s2->coefN)
      return true;
    else
      return false;

    for (int i = 0; i < size; i++) {
      if (s1->elem[i] == s2->elem[i])
        return true;
      else
        return false;
    }
    return false;
  }
  int size;
};

class combProps {
public:
  combProps(int j_, int nn_);
  ~combProps();
  std::vector<int*> vecComb;
  int nn, j;
};

#endif
