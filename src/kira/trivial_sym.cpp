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

// jule is the number of propagators

#include <sys/stat.h>

#include <algorithm>
#include <tuple>

#include "pyred/coeff_int.h"
#include "pyred/interface.h"
#include "pyred/parser.h"
#include "kira/kira.h"
#include "kira/tools.h"
#include "kira/trivial_sym.h"

using namespace std;
using namespace GiNaC;

static Loginfo& logger = Loginfo::instance();

Pak::Pak(ex listOfTerms_, uint32_t jule_, vector<int>& seed, vector<int>& holes) {
  l_Indices = seed.size();
  jule = jule_;
  possymbol* feynman;
  feynman = new possymbol[jule];
  generate_symbols(feynman, "b", jule);

  std::vector<TermsGiNaC*> canonicaGiNaC;

  for (size_t itTerms = 0; itTerms < listOfTerms_.nops(); itTerms++) {
    ex termEx = listOfTerms_.op(itTerms);

    TermsGiNaC* term = new TermsGiNaC(l_Indices);

    for (auto itF : seed) {
      int good = 0;
      for (int itDegree = 1; itDegree < 100; itDegree++) {
        if (termEx.coeff(feynman[itF], itDegree) != 0) {
          termEx = termEx.coeff(feynman[itF], itDegree);
          term->elem.push_back(itDegree);
          good++;
        }
      }
      if (!good) {
        term->elem.push_back(0);
      }
    }
    int good = 1;
    for (auto itF : holes) {
      for (int itDegree = 1; itDegree < 100; itDegree++) {
        if (termEx.coeff(feynman[itF], itDegree) != 0) {
          good = 0;
        }
      }
    }

    if (good) {
      term->coef = something_string(termEx);
      term->coefEx = termEx;
      canonicaGiNaC.push_back(term);
    }
    else
      delete term;
  }

  collect_trivial(canonicaGiNaC);

  for (auto itR : canonicaGiNaC) {
    Terms* term = new Terms(l_Indices);

    term->coef = itR->coef;
    term->coefN = pyred::parse_coeff<pyred::Coeff_int>(itR->coef).hash();

    for (auto itB : itR->elem) {
      term->elem.push_back(itB);
    }
    term->equation = itR->equation;

    term->size = itR->size;

    canonicaX.push_back(term);
  }

  for (auto itR : canonicaGiNaC)
    delete itR;

  delete[] feynman;
}

void Pak::collect_trivial(std::vector<TermsGiNaC*>& canonicaGiNaC) {
  unsigned end = canonicaGiNaC.size();

  for (size_t itX = 0; itX != end; itX++) {
    for (size_t itY = itX + 1; itY != end; itY++) {
      int count = 0;
      for (int i = 0; i < l_Indices; i++) {
        if (canonicaGiNaC[itX]->elem[i] == canonicaGiNaC[itY]->elem[i]) count++;
      }
      if (count == l_Indices) {
        canonicaGiNaC[itX]->coef += ("+" + canonicaGiNaC[itY]->coef);
        canonicaGiNaC[itX]->coefEx += canonicaGiNaC[itY]->coefEx;
        delete canonicaGiNaC[itY];
        canonicaGiNaC.erase(canonicaGiNaC.begin() + itY);
        end--;
        itY--;
      }
    }
  }

  possymbol* unknown = new possymbol[jule];
  generate_symbols(unknown, "k", jule);

  for (size_t itX = 0; itX != canonicaGiNaC.size(); itX++) {
    for (uint32_t i = 0; i < jule; i++) {
      ex term = canonicaGiNaC[itX]->coefEx.coeff(unknown[i], 1);
      if (term != 0) {
        canonicaGiNaC[itX]->equation.push_back(
            make_pair(pyred::Weight(i + 1), something_string(term)));
        canonicaGiNaC[itX]->coefEx =
            canonicaGiNaC[itX]->coefEx.subs(unknown[i] == 0);
      }
    }

    canonicaGiNaC[itX]->equation.push_back(
        make_pair(pyred::Weight(0u), something_string(canonicaGiNaC[itX]->coefEx)));
  }
  delete [] unknown;
}

void Pak::check_trivial(vector<int>& /*seed*/, vector<int>& holes, int& check) {
  auto sys = pyred::System();

  int count = 0;
  for (auto itX : canonicaX) {
    int good = 1;
    for (auto itF : holes) {
      if (itX->elem[itF]) {
        good = 0;
        break;
      }
    }
    if (!good) continue;

    sys.add(itX->equation);
    count++;
  }
  if (count != 0) {
    //     {
    //       std::lock_guard<std::mutex> lock(m);
    sys.solve();
    //     }

    auto content = sys.reduction_content();

    if (content.find(pyred::Weight(0u)) == content.end()) {
      check = 1;
    }
    else {
      check = 0;
    }
  }
  else {
    check = 1;
  }
}

int Kira::find_zero_sectors() {
  logger << "\n*****Search trivial sectors********************************\n";
  Clock clock;

  int numberOfNonTrivialSectors = 0;

  int nn = integralfamily.loopVar.size();
  vector<ex> vec;
  possymbol token("token");

  /*generate Feynman parameters*/

  possymbol* feynman = new possymbol[integralfamily.jule];
  generate_symbols(feynman, "b", integralfamily.jule);

  ex schwinger;
  for (uint32_t i = 0; i < integralfamily.jule; i++) {
    schwinger += -feynman[i] * integralfamily.props[i].expand();
  }

  lst kinematic4sym;

  for (size_t itE = 0; itE < kinematic.nops(); itE++)
    kinematic4sym.append(kinematic[itE].subs(invariants4sym));

  fs<lst>(schwinger, kinematic4sym);

  /*Create mZS, collect coefficients for loop_mom O(2)*/
  lst mZS;
  for (int i = 0; i < nn; i++) {
    ex superSchwinger = schwinger;
    for (int j = 0; j < nn; j++) {
      fs<relational>(
          superSchwinger,
          (integralfamily.loopVar[i] * integralfamily.loopVar[j] == token));
      if (i == j) {
        mZS.append(-superSchwinger.coeff(token, 1));
      }
      else {
        mZS.append(-superSchwinger.coeff(token, 1) / 2);
      }
      superSchwinger = superSchwinger.subs(token == 0);
    }
  }

  for (int i = 0; i < nn; i++) {
    for (int j = i; j < nn; j++) {
      fs<relational>(
          schwinger,
          (integralfamily.loopVar[i] * integralfamily.loopVar[j] == token));
      fs<relational>(schwinger, token == 0);
    }
  }

  matrix mZS2(integralfamily.loopVar.size(), integralfamily.loopVar.size(),
              mZS);
  ex mZS2det = mZS2.determinant();

  /*Create vector, collect coefficients for loop_mom O(1)*/
  for (int i = 0; i < nn; i++) {
    fs<relational>(schwinger, integralfamily.loopVar[i] == token);
    vec.push_back(schwinger.coeff(token, 1) / 2);
    fs<relational>(schwinger, token == 0);
  }

  ex J = schwinger;
  integralfamily.U = mZS2det.expand();

  lst mZSAdj;

  for (int k = 0; k < nn; k++) {
    for (int l = 0; l < nn; l++) {
      lst adj;
      int countAdj = 0;

      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          if (i == k && j == l) {
            adj.append(1);
            countAdj++;
          }
          else if (i == k || j == l) {
            adj.append(0);
            countAdj++;
          }
          else
            adj.append(mZS[countAdj++]);
        }
      }
      matrix adjMatr(integralfamily.loopVar.size(),
                     integralfamily.loopVar.size(), adj);

      mZSAdj.append(adjMatr.determinant());
    }
  }

  ex ZScalar = 0;
  int countM = 0;
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      ZScalar += mZSAdj[countM++] * vec[i] * vec[j];
    }
  }

  ex F = (ZScalar + J * mZS2det).expand();

  fs<lst>(F, kinematic4sym);

  F = expand(F.subs(kinematicReverse));

  integralfamily.F = F;

  //   logger << integralfamily.U << "\n\n";
  //   logger << integralfamily.F << "\n\n";

  integralfamily.Gsym = integralfamily.U + integralfamily.F;
  integralfamily.G = expand(integralfamily.Gsym.subs(invariants4symRev));

  ifstream input;
  string inputName = outputDir + "/sectormappings/" + integralfamily.name +
                     "/nonTrivialSector";
  if (file_exists(inputName.c_str())) {
    input.open(inputName.c_str());
    logger << "read from: " << inputName << "\n";
    while (1) {
      int sector;
      int nOfProps;
      if (!(input >> sector)) break;
      if (!(input >> nOfProps)) break;
      numberOfNonTrivialSectors++;
      integralfamily.mask[nOfProps].push_back(sector);
    }
    input.close();

    logger << "\nNon trivial sectors in total: " << numberOfNonTrivialSectors
           << "\n";
    logger << "Trivial sectors in total: "
           << ((1u << integralfamily.jule) - numberOfNonTrivialSectors) << "\n";
    logger << "( " << clock.eval_time() << " s )\n\n";

    delete[] feynman;

    if (numberOfNonTrivialSectors == 0) {
      logger << "The system you try to reduce contains only zero integrals.\n";
      return 1;
    }
    else
      return 0;
  }

  /*initialize Lee's zero sector criterium*/
  possymbol* unknown = new possymbol[integralfamily.jule];
  generate_symbols(unknown, "k", integralfamily.jule);
  ex kriterium = 0, kriterium2 = 0, kritU = 0, kritF = 0;
  for (uint32_t i = 0; i < integralfamily.jule; i++) {
    kriterium += unknown[i] * feynman[i] * integralfamily.G.diff(feynman[i]);
    kritU += unknown[i] * feynman[i] * integralfamily.U.diff(feynman[i]);
    kritF += unknown[i] * feynman[i] * integralfamily.F.diff(feynman[i]);
  }
  kriterium = expand(kriterium - integralfamily.G);
  kritU = expand(kritU - integralfamily.U);
  kritF = expand(kritF - integralfamily.F);

  vector<int> seedInit;
  vector<int> holesInit;

  for (uint32_t y = 0; y < integralfamily.jule; y++) {
    if (((1u << y) & ((1u << integralfamily.jule) - 1))) {
      seedInit.push_back(y);
    }
    else {
      holesInit.push_back(y);
    }
  }

  // prepare pyRed to be used by Pak
  pyred::Config::verbosity(0);
  pyred::Config::insertion_tracer(1);
  // store prime set by user to restore it later
  uint64_t saveUserPrime = pyred::Coeff_int::get_prime();
  pyred::Coeff_int::set_prime(9223372036854775783);

  Pak pak(kritF + kritU, integralfamily.jule, seedInit, holesInit);

  uint32_t countLoop = 0;
  uint32_t loopSIZE = (1 << integralfamily.jule);
  std::mutex m;
  std::condition_variable cond_var;
  bool processed = false;

  if (loopSIZE != 0) {
    for (uint32_t sector_number = 0; sector_number < (1u << integralfamily.jule);
         sector_number++) {
      pool->enqueue([this, sector_number, &pak, &numberOfNonTrivialSectors,
                     &countLoop, &processed, &loopSIZE, &cond_var, &m]() {
        int userSetZero = 0;
        for (uint32_t jIt : integralfamily.userZeroSectors) {
          if (sector_number == jIt) {
            std::lock_guard<std::mutex> lock(m);
            integralfamily.zeroSector.push_back(sector_number);
            userSetZero = 1;
          }
        }

        if (!userSetZero) {
          unsigned forceTrivial;
          if (integralfamily.cutProps.size())
            forceTrivial = integralfamily.cutProps.size();
          else
            forceTrivial = 0;

          for (unsigned j = 0; j < integralfamily.cutProps.size(); j++) {
            if ((1 << (integralfamily.cutProps[j] - 1)) & sector_number) {
              forceTrivial--;
            }
          }

          if (forceTrivial > 0) {
            {
              std::lock_guard<std::mutex> lock(m);
              integralfamily.zeroSector.push_back(sector_number);
            }
          }
          else {
            vector<int> seedCheck;
            vector<int> holesCheck;
            for (uint32_t y = 0; y < integralfamily.jule; y++) {
              if (((1 << y) & sector_number)) {
                seedCheck.push_back(y);
              }
              else {
                holesCheck.push_back(y);
              }
            }
            int check = 0;
            pak.check_trivial(seedCheck, holesCheck, check);

            /*symbolic IBP no trivial sectors*/ //???
            for (auto itSymb : integralfamily.symbolicIBP) {
              if (check == 1 && !((1 << itSymb) & sector_number)) {
                check = 0;
              }
            }

            if (check) {
              {
                std::lock_guard<std::mutex> lock(m);
                integralfamily.zeroSector.push_back(sector_number);
              }
            }
            else {
              int count = 0;

              for (uint32_t j = 0; j < integralfamily.jule; j++) {
                if ((1 << j) & sector_number) count++;
              }
              {
                std::lock_guard<std::mutex> lock(m);
                numberOfNonTrivialSectors++;
                integralfamily.mask[count].push_back(sector_number);
              }
            }
          }
        }
        {
          std::lock_guard<std::mutex> lock(m);
          countLoop++;
          load_bar(sector_number + 1, (1u << integralfamily.jule), 50, 100);
          if (countLoop == loopSIZE) {
            processed = true;
            cond_var.notify_one();
          }
        }
      });
    }

    {
      std::unique_lock<std::mutex> lock(m);
      cond_var.wait(lock, [&processed] { return processed; });
    }
  }

  /*Write results*/
  for (uint32_t i = 0; i < integralfamily.jule + 1; i++) {
    std::sort(integralfamily.mask[i].begin(), integralfamily.mask[i].end());

    for(size_t ix = 0; ix < integralfamily.mask[i].size(); ix++){
      for(size_t ix2 = 0; ix2 < integralfamily.mask[i].size(); ix2++){
        if(integralfamily.sectorweight[integralfamily.mask[i][ix]].first <
          integralfamily.sectorweight[integralfamily.mask[i][ix2]].first){
          swap(integralfamily.mask[i][ix],integralfamily.mask[i][ix2]);
        }
      }
    }

    integralfamily.mask[i].resize(distance(
        integralfamily.mask[i].begin(),
        unique(integralfamily.mask[i].begin(), integralfamily.mask[i].end())));
  }
  ItVint endZero = unique(integralfamily.zeroSector.begin(),
                          integralfamily.zeroSector.end());
  integralfamily.zeroSector.resize(
      distance(integralfamily.zeroSector.begin(), endZero));
  sort(integralfamily.zeroSector.begin(), integralfamily.zeroSector.end());
  logger << "\nNon trivial sectors in total: " << numberOfNonTrivialSectors
         << "\n";

  ofstream fileZeroSector;
  fileZeroSector.open(
      (outputDir + "/sectormappings/" + integralfamily.name + "/trivialsector")
          .c_str());

  for (unsigned it = 0, end1 = integralfamily.zeroSector.size(); it < end1;
       it++) {
    if (it < end1 - 1)
      fileZeroSector << integralfamily.zeroSector[it] << ",";
    else
      fileZeroSector << integralfamily.zeroSector[it];
  }

  fileZeroSector.close();

  ofstream nonTrivialSector;
  nonTrivialSector.open((outputDir + "/sectormappings/" + integralfamily.name +
                         "/nonTrivialSector.back")
                            .c_str());

  for (uint32_t i = 0; i < integralfamily.jule + 1; i++) {
    for (unsigned it = 0; it < integralfamily.mask[i].size(); it++) {
      nonTrivialSector << integralfamily.mask[i][it] << " " << i << endl;
    }
  }

  nonTrivialSector.close();

  rename((outputDir + "/sectormappings/" + integralfamily.name +
          "/nonTrivialSector.back")
             .c_str(),
         (outputDir + "/sectormappings/" + integralfamily.name +
          "/nonTrivialSector")
             .c_str());

  remove((outputDir + "/sectormappings/" + integralfamily.name +
          "/nonTrivialSector.back")
             .c_str());

  logger << "Trivial sectors in total: " << integralfamily.zeroSector.size()
         << "\n";
  logger << "( " << clock.eval_time() << " s )\n\n";

  delete[] feynman;
  delete[] unknown;

  // reset to prime set by user
  pyred::Coeff_int::set_prime(saveUserPrime);

  if (numberOfNonTrivialSectors == 0) {
    logger << "The system you try to reduce contains only zero integrals.\n";
    return 1;
  }
  else
    return 0;
}

int Kira::testProps(lst& matr, string& topoName) {
  lst mZS;
  for (size_t i = 0; i < matr.nops(); i++) {
    for (size_t j = 0; j < matr.nops(); j++) {
      mZS.append(diff(matr[i], ex_to<symbol>(topology[topoName].loopVarList[j])));
    }
  }
  matrix mZS2(topology[topoName].loopVarList.nops(),
              topology[topoName].loopVarList.nops(), mZS);

  int det = 0;
  det = something_int(mZS2.determinant());

  return det;
}

int Kira::test(lst& matr, string& topoName) {
  lst mZS;
  for (size_t i = 0; i < matr.nops(); i++) {
    for (size_t j = 0; j < matr.nops(); j++) {
      mZS.append(diff(integralfamily.loopVarList[i].subs(matr[i]),
                      ex_to<symbol>(topology[topoName].loopVarList[j])));
    }
  }
  matrix mZS2(topology[topoName].loopVarList.nops(),
              topology[topoName].loopVarList.nops(), mZS);

  int det = 0;
  det = something_int(mZS2.determinant());

  return det;
}

bool sortVariables(const tuple<size_t, size_t, size_t>& i,
                   const tuple<size_t, size_t, size_t>& j) {
  if (get<0>(i) < get<0>(j)) return 1;
  if (get<0>(i) > get<0>(j)) return 0;

  return false;
}

bool myfunction(const symmetries& i, const symmetries& j) {
  if (i.topology < j.topology) return 1;
  if (i.topology > j.topology) return 0;

  if (i.sector < j.sector) return 1;
  if (i.sector > j.sector) return 0;

  return (i.symDOTS < j.symDOTS);
}

bool isEqual(const symmetries& i, const symmetries& j) {
  return (i.topology == j.topology && i.sector == j.sector &&
          i.symDOTS == j.symDOTS);
}

ostream& operator<<(ostream& out, const Terms& term) {
  out << std::setw(8) << term.coef << " ";
  out << std::setw(5) << term.coefN << " ";
  for (int i = 0; i < term.size; i++) {
    out << term.elem[i] << " ";
  }
  out << endl;

  return out;
}

// bool Pak::compare(Pak& pak){
//
//   if(canonica.size() != pak.canonica.size())
//     return false;
//
//   if(equal(canonica.begin(), canonica.end(), pak.canonica.begin(),
//   equalAllMatrix(l_Indices)))
//     return true;
//   else
//     return false;
//
// }

void Pak::collect(std::vector<Terms*>& matrix) {
  unsigned end = matrix.size();
  for (size_t itX = 0; itX != end; itX++) {
    for (size_t itY = itX + 1; itY != end; itY++) {
      int count = 0;
      for (int i = 0; i < l_Indices; i++) {
        if (matrix[itX]->elem[i] == matrix[itY]->elem[i]) count++;
      }
      if (count == l_Indices) {
        matrix[itX]->coef += ("+" + matrix[itY]->coef);
        delete matrix[itY];
        matrix.erase(matrix.begin() + itY);
        // 	swap(matrix[itY], matrix[end-1]);
        end--;
        itY--;
      }
    }
  }
  for (size_t itX = 0; itX != matrix.size(); itX++) {
    matrix[itX]->coefN =
        pyred::parse_coeff<pyred::Coeff_int>((matrix[itX]->coef)).hash();
  }
}

Pak::Pak(ex listOfTerms_, uint32_t jule, vector<int>& seed, vector<int>& holes,
         std::vector<std::tuple<GiNaC::lst, GiNaC::lst, int, GiNaC::ex,
                                std::vector<std::string> > >& externalTransf,
         GiNaC::lst& invariantsReplacementRev, unsigned maskCut) {
  l_Indices = seed.size();
  possymbol* feynman;
  feynman = new possymbol[jule];
  generate_symbols(feynman, "b", jule);

  for (auto itExt : externalTransf) {
    canonicaX.clear();

    for (size_t itTerms = 0; itTerms < listOfTerms_.nops(); itTerms++) {
      ex termEx = listOfTerms_.op(itTerms);

      Terms* term = new Terms(l_Indices);
      term->powerSum = 0;
      term->cutPropagator = 0;
      for (auto itF : seed) {
        int good = 0;
        for (int itDegree = 1; itDegree < 100; itDegree++) {
          if (termEx.coeff(feynman[itF], itDegree) != 0) {
            termEx = termEx.coeff(feynman[itF], itDegree);
            term->elem.push_back(itDegree);
            term->powerSum += itDegree;
            if(((1 << itF) & maskCut)){
              term->cutPropagator = 1;
            }
            good++;
          }
        }
        if (!good) {
          term->elem.push_back(0);
        }
      }
      int good = 1;
      for (auto itF : holes) {
        for (int itDegree = 1; itDegree < 100; itDegree++) {
          if (termEx.coeff(feynman[itF], itDegree) != 0) {
            good = 0;
          }
        }
      }

      if (good) {
        termEx = termEx.subs(get<3>(itExt), subs_options::algebraic)
                     .subs(invariantsReplacementRev, subs_options::algebraic);

        term->coef = something_string(termEx);

        canonicaX.push_back(term);
      }
      else
        delete term;
    }

    /*generates coefN as well*/
    collect(canonicaX);

    std::sort(canonicaX.begin(), canonicaX.end(), compMatrix(0));

    canonica.push_back(canonicaX);
  }

  delete[] feynman;
}

int check_variable(std::string coef, uint64_t coefN,
                   std::vector<std::string>& invarSol,
                   std::vector<std::string>& replaceString, int mask) {
  coef = "(" + coef + ")";

  std::vector<std::tuple<size_t, size_t, size_t> > posVarVec;

  vector<string> signs = {"*", "+", "-", "("};

  for (size_t itY = 0; itY < invarSol.size(); itY++) {
    if (((1 << itY) & mask)) {
      string foundStr1;

      vector<string> testString = {"*" + invarSol[itY], "+" + invarSol[itY],
                                   "-" + invarSol[itY], "(" + invarSol[itY]};

      for (size_t it = 0; it < testString.size(); it++) {
        size_t posVar = coef.find(testString[it]);

        if (posVar != std::string::npos) {
          posVarVec.push_back(make_tuple(posVar, itY, it));
        }
      }
    }
  }

  sort(posVarVec.begin(), posVarVec.end(), sortVariables);

  for (size_t itX = 0; itX < posVarVec.size(); itX++) {
    string foundStr1 = coef.substr(get<0>(posVarVec[itX]) + 1 +
                                   (invarSol[get<1>(posVarVec[itX])].size()));

    coef = coef.substr(0, get<0>(posVarVec[itX])) +
           signs[get<2>(posVarVec[itX])] + "(" +
           replaceString[get<1>(posVarVec[itX])] + ")" +
           coef.substr(get<0>(posVarVec[itX]) + 1 +
                       (invarSol[get<1>(posVarVec[itX])].size()));

    for (size_t itZ = itX + 1; itZ < posVarVec.size(); itZ++) {
      get<0>(posVarVec[itZ]) += replaceString[get<1>(posVarVec[itX])].size() -
                                (invarSol[get<1>(posVarVec[itX])].size()) + 2;
    }
  }
  uint64_t coefC = pyred::parse_coeff<pyred::Coeff_int>(coef).hash();

  if ((coefN - coefC) != 0) {
    return 1;
  }

  //   {
  //     string foundStr1;
  //     vector<string> testString =
  //     {"*"+var+"^","+"+var+"^","-"+var+"^","("+var+"^"};
  //
  //     for(auto it: testString){
  //
  //       size_t posVar = coef.find(it);
  //
  //       if (posVar!=std::string::npos){
  //
  // 	return 1;
  //       }
  //     }
  //   }
  return 0;
}

std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> >
Pak::generate_combinatorics(
    vector<int>& seed, vector<int>& holes,
    std::vector<std::tuple<GiNaC::lst, GiNaC::lst, int, GiNaC::ex,
                           std::vector<std::string> > >& externalTransf,
    std::vector<std::string>& invarSol) {

  l_Indices = seed.size();

  std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> > symmetriesVec;
  std::vector<std::tuple<size_t, std::vector<int*>, int, int> >
      symmetriesCrossedVec;

  unsigned countCanonica = 0;

  std::vector<uint64_t> originalCanonic;

  for (auto itCanonica : canonica) {

    std::vector<Terms*> matrix;
    std::vector<std::vector<int> > permutation;
    std::vector<std::tuple<std::vector<Terms*>, std::vector<int> > > matrixV;

    for (auto itX : itCanonica) {
      int good = 1;
      for (auto itF : holes) {
        if (itX->elem[itF]) {
          good = 0;
          break;
        }
      }
      if (!good) continue;

      Terms* term = new Terms(l_Indices);

      for (auto itF : seed) {
        term->elem.push_back(itX->elem[itF]);
      }
      term->coefN = itX->coefN;
      term->coef = itX->coef;
      term->powerSum = itX->powerSum;
      term->cutPropagator = itX->cutPropagator;

      matrix.push_back(term);
    }

    std::sort(matrix.begin(), matrix.end(), compMatrix(0));

    vector<int> tmp;
    matrixV.push_back(make_tuple(matrix, tmp));

    for (int i = 0; i < l_Indices - 1; i++) {
      vector<tuple<vector<int>, int, int> > collumnVec;

      for (size_t itM = 0; itM != matrixV.size(); itM++) {
        for (int it = i; it < l_Indices; it++) {
          std::sort(get<0>(matrixV[itM]).begin(), get<0>(matrixV[itM]).end(),
                    compMatrix2(i + 1, it));

          vector<int> collumn;

          for (auto itX : get<0>(matrixV[itM])) {
            collumn.push_back((*itX).elem[it]);
          }

          collumnVec.push_back(make_tuple(collumn, itM, it));
        }
      }

      std::sort(collumnVec.begin(), collumnVec.end());

      vector<int> eraseElements2;

      for (size_t j = 0; j < collumnVec.size() - 1; j++) {
        if (!equal(get<0>(collumnVec.back()).begin(),
                   get<0>(collumnVec.back()).end(),
                   get<0>(collumnVec[j]).begin())) {
          eraseElements2.push_back(j);
        }
      }

      for (auto it = eraseElements2.rbegin(); it != eraseElements2.rend();
           ++it) {
        collumnVec.erase(collumnVec.begin() + (*it));
      }

      vector<tuple<vector<Terms*>, vector<int> > > matrixVtmp2 = matrixV;

      matrixV.clear();

      for (auto itM : collumnVec) {
        vector<Terms*> matrix;

        for (size_t itX = 0; itX < get<0>(matrixVtmp2[get<1>(itM)]).size();
             itX++) {
          Terms* term = new Terms(l_Indices);
          (*term).coefN = (*get<0>(matrixVtmp2[get<1>(itM)])[itX]).coefN;
          (*term).coef = (*get<0>(matrixVtmp2[get<1>(itM)])[itX]).coef;
          (*term).elem = (*get<0>(matrixVtmp2[get<1>(itM)])[itX]).elem;
          (*term).powerSum = (*get<0>(matrixVtmp2[get<1>(itM)])[itX]).powerSum;
          (*term).cutPropagator = (*get<0>(matrixVtmp2[get<1>(itM)])[itX]).cutPropagator;
          swap((*term).elem[i], (*term).elem[get<2>(itM)]);

          matrix.push_back(term);
        }
        vector<int> tmp = get<1>(matrixVtmp2[get<1>(itM)]);
        tmp.push_back(get<2>(itM));
        matrixV.push_back(make_tuple(matrix, tmp));
      }

      for (auto itM = matrixV.begin(); itM != matrixV.end(); itM++) {
        std::sort(get<0>(*itM).begin(), get<0>(*itM).end(), compMatrix(i + 1));
      }

      for (auto itX : matrixVtmp2) {
        for (auto itY : get<0>(itX)) {
          delete itY;
        }
      }
    }

    for (auto itM = matrixV.begin(); itM != matrixV.end(); itM++) {
      std::sort(get<0>(*itM).begin(), get<0>(*itM).end(),
                compMatrix(l_Indices + 1));
    }

    for (size_t itM = 0; itM < matrixV.size(); itM++) {
      vector<int> perArray;

      for (int it = 0; it < l_Indices; it++) {
        perArray.push_back(seed[it]);
      }

      int count = 0;
      for (auto itX : get<1>(matrixV[itM])) {
        iter_swap(perArray.begin() + count++, perArray.begin() + itX);
      }

      permutation.push_back(perArray);
    }
    stringstream inputStr;

    int crossedFlag = 0;
    for (size_t itX = 0; itX < get<0>(matrixV.front()).size(); itX++) {

      if (check_variable(get<0>(matrixV.front())[itX]->coef,
                         get<0>(matrixV.front())[itX]->coefN, invarSol,
                         get<4>(externalTransf[countCanonica]),
                         get<2>(externalTransf[countCanonica]))) {
        crossedFlag = 1;
      }

      inputStr << get<0>(matrixV.front())[itX]->coefN << " ";
    }
    for (int i = 0; i < l_Indices; i++) {
      for (auto itX : get<0>(matrixV.front())){
        inputStr << itX->elem[i] << " ";
      }
    }

    std::hash<std::string> hash_fn;
    size_t str_hash = hash_fn(inputStr.str());

    for (size_t itM = 0; itM != matrixV.size(); itM++) {
      for (auto itV : std::get<0>(matrixV[itM])) {
        delete itV;
      }
    }
    matrixV.clear();
    // check if cross or normal symmetry
    if (!crossedFlag) {
      symmetriesVec.push_back(
          make_tuple(str_hash, permutation, l_Indices, countCanonica));
    }
    countCanonica++;
  }

  return symmetriesVec;
}

void Kira::search_symmetry_relations() {
  logger << "Search symmetry relations: \n";

  // store prime set by user to restore it later
  auto saveUserPrime = pyred::Coeff_int::get_prime();
  pyred::Coeff_int::set_prime(9223372036854775783);

  Clock clock;

  ifstream input;
  string inputName =
      (outputDir + "/sectormappings/" + integralfamily.name + "/canonica");
  int countLines = 0;

  if (file_exists(inputName.c_str())) {
    input.open(inputName.c_str());
    logger << "read from: " << inputName << "\n";

    while (1) {
      int sector;
      if (!(input >> sector)) break;

      size_t hashNumber;
      if (!(input >> hashNumber)) break;

      int vectorSize;
      if (!(input >> vectorSize)) break;

      int sizeArray;
      if (!(input >> sizeArray)) break;

//       vector<int*> vecInts;
      vector<vector<int> > vecInts;
      for (int itB = 0; itB < vectorSize; itB++) {
//         int* arrayIn = new int[sizeArray];
        vector<int> arrayIn;
        for (int itA = 0; itA < sizeArray; itA++) {
          arrayIn.push_back(0);
        }
        for (int itA = 0; itA < sizeArray; itA++) {
          if (!(input >> arrayIn[itA])) break;
        }
        vecInts.push_back(arrayIn);
      }

      countLines++;

      std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> > tokenVec;

      tokenVec.push_back(make_tuple(hashNumber, vecInts, sizeArray, 0));

      integralfamily.symmetries.insert(
          std::pair<
              int,
              std::vector<std::tuple<size_t, std::vector<std::vector<int> >, int, int> > >(
              sector, tokenVec));
    }
    logger << "Number of symmetries read out: " << countLines << "\n";
    input.close();
  }
  else {
    vector<int> seed1;
    vector<int> holes1;
    for (uint32_t y = 0; y < integralfamily.jule; y++) {
      if (((1 << y) & ((1 << integralfamily.jule) - 1))) {
        seed1.push_back(y);
      }
      else {
        holes1.push_back(y);
      }
    }

    for (auto itE : invarSol)
      pyred::parse_coeff<pyred::Coeff_int>(itE);

    Pak pak(integralfamily.Gsym, integralfamily.jule, seed1, holes1,
            externalTransf, invariantsReplacementRev, integralfamily.maskCut);

    for (int j = 0; j < integralfamily.biggestBound + 1; ++j) {
      uint32_t countLoop = 0;
      uint32_t loopSIZE = integralfamily.mask[j].size();
      std::mutex m;
      std::condition_variable cond_var;
      bool processed = false;

      if (loopSIZE != 0) {
        for (size_t it = 0; it < integralfamily.mask[j].size(); ++it) {
          pool->enqueue([this, j, it, &pak, &countLoop, &processed, &loopSIZE,
                         &cond_var, &m]() {

            vector<int> seed1;
            vector<int> holes1;
            for (uint32_t y = 0; y < integralfamily.jule; y++) {
              if (((1 << y) & integralfamily.mask[j][it])) {
                seed1.push_back(y);
              }
              else {
                holes1.push_back(y);
              }
            }
            auto symmetries = pak.generate_combinatorics(
                seed1, holes1, externalTransf, invarSol);

            {
              std::lock_guard<std::mutex> lock(m);

              integralfamily.symmetries.insert(
                  std::pair<int, std::vector<std::tuple<
                                     size_t, std::vector<vector<int> >, int, int> > >(
                      integralfamily.mask[j][it], symmetries));

              countLoop++;
              if (countLoop == loopSIZE) {
                processed = true;
                cond_var.notify_one();
              }
            }
          });
        }
        {
          std::unique_lock<std::mutex> lock(m);
          cond_var.wait(lock, [&processed] { return processed; });
        }
      }
      load_bar(j, integralfamily.biggestBound, 50, 100);
    }

    ofstream fileSymVec;
    fileSymVec.open(
        (outputDir + "/sectormappings/" + integralfamily.name + "/canonica"));

    for (auto itX : integralfamily.symmetries) {
      if ((itX.second).size()) {
        fileSymVec << (itX.first) << " ";

        auto symmetry = (itX.second).front();
        fileSymVec << get<0>(symmetry) << " ";
        fileSymVec << get<1>(symmetry).size() << " ";
        fileSymVec << get<2>(symmetry) << "\n";

        for (auto itX2 : get<1>(symmetry)) {
          for (int itA = 0; itA < get<2>(symmetry); itA++)
            fileSymVec << itX2[itA] << " ";
          fileSymVec << "\n";
        }
        fileSymVec << "\n";
      }
    }
    fileSymVec.close();
  }

  // reset to prime set by user
  pyred::Coeff_int::set_prime(saveUserPrime);

  logger << "\n( " << clock.eval_time() << " s )\n";
}

combProps::combProps(int j_, int nn_) {

  j = j_;
  nn = nn_;

  int *arrayTest, *arrayBU, *arraySave;
  arrayBU = new int[nn];

  arrayTest = new int[nn];

  for (int it = 0; it < nn; it++)
    arrayTest[it] = it;

  int end = j - 1;

  for (int it = 0; it < nn; it++)
    arrayBU[it] = arrayTest[it];

  int go = 1;

  if (nn == 1) {
    arraySave = new int[nn];

    for (int it = 0; it < nn; it++) {
      arraySave[it] = arrayTest[it];
    }

    vecComb.push_back(arraySave);
  }
  else {
    while (go == 1) {
      for (int it = 0; it < nn; it++)
        arrayTest[it] = arrayBU[it];

      int inIt = arrayTest[nn - 1];

      arraySave = new int[nn];

      for (int it = 0; it < nn; it++) {
        arraySave[it] = arrayTest[it];
      }

      vecComb.push_back(arraySave);

      while (inIt < end) {
        inIt = ++arrayTest[nn - 1];

        arraySave = new int[nn];

        for (int it = 0; it < nn; it++) {
          arraySave[it] = arrayTest[it];
        }

        vecComb.push_back(arraySave);
      }

      for (int it = nn - 2; it >= 0; it--) {
        if (arrayBU[it] < end - ((nn - it) - 1)) {
          arrayBU[it]++;

          for (int itt = it + 1; itt < nn; itt++)
            arrayBU[itt] = arrayBU[itt - 1] + 1;

          break;
        }
        if (it == 0) go = 0;
      }
    }
  }
  delete[] arrayTest;
  delete[] arrayBU;
}

combProps::~combProps() {
  for (auto j : vecComb) {
    delete[] j;
  }
}

int Kira::prepare_symmetry() {
  for (int j = 0; j < integralfamily.biggestBound + 1; ++j) {
    for (size_t itt = 0; itt != integralfamily.mask[j].size();
         ++itt) { // all possible non trivial sectors

      int sectorNumber = integralfamily.mask[j][itt];
      int num_ones = 0;
      unsigned testSector = sectorNumber;
      for (size_t i = 0; i < SEEDSIZE; ++i, testSector >>= 1) {
        if ((testSector & 1) == 1) ++num_ones;
      }
      if (num_ones < SEEDSIZE) {
        int flagy = 0;
        for (size_t it = 0; it < integralfamily.lowestSectors.size(); it++) {
          if ((integralfamily.lowestSectors[it] & sectorNumber) ==
              integralfamily.lowestSectors[it]) {
            flagy = 1;
            break;
          }
        }
        if (!flagy) {
          integralfamily.lowestSectors.push_back(sectorNumber);
        }
      }
    }
  }

  if (file_exists(
          (outputDir + "/sectormappings/" + integralfamily.name + "/magicFlag")
              .c_str())) {
    ifstream input;
    input.open(
        (outputDir + "/sectormappings/" + integralfamily.name + "/magicFlag")
            .c_str());

    string magicFlag;

    input >> magicFlag;
    if (integralfamily.magic_relations != magicFlag) {
      std::string msg = "";
      msg += "You have ones generated the symmetries with the\n";
      msg += "option magic_relations switched ";
      if (magicFlag == "true")
        msg += "on.\n";
      else
        msg += "off.\n";

      msg += "If you change the option magic_relations compared\n";
      msg += "to the previous run you need to delete the directory:\n";
      msg += outputDir + "/sectormappings/" + integralfamily.name;
      throw ExceptionResume(msg);
    }
  }
  else {
    ofstream output;
    output.open(
        (outputDir + "/sectormappings/" + integralfamily.name + "/magicFlag")
            .c_str());

    string magicFlag;

    output << integralfamily.magic_relations;
  }

  int flagSYM = read_symmetries((outputDir + "/sectormappings/" +
                                 integralfamily.name + "/sectorRelations"),
                                integralfamily.symVec);
  int flagREL = read_symmetries((outputDir + "/sectormappings/" +
                                 integralfamily.name + "/sectorSymmetries"),
                                integralfamily.relVec);

  if (flagSYM && flagREL) return 1;
  return 0;
}

int Kira::symmetry_finder(int it, int itt, string topoName, int j, std::vector<int> array2,
                          std::vector<int> array3, int flag, unsigned klop) {
  // check if cut propagators are mixed with uncut

  for (int g = 0; g < j; g++) {
    if (!(!(1 << array2[g] & integralfamily.maskCut) ==
          !(1 << array3[g] & topology[topoName].maskCut))) {
      return 0;
    }
  }

  // checkSymbolic IBP
  size_t theBig =
      integralfamily.symbolicIBP.size() > topology[topoName].symbolicIBP.size()
          ? integralfamily.symbolicIBP.size()
          : topology[topoName].symbolicIBP.size();

  if (theBig > 0) {
    size_t checkV = 0;

    for (int g = 0; g < j; g++) { //???

      auto itSY = find(integralfamily.symbolicIBP.begin(),
                       integralfamily.symbolicIBP.end(), array2[g]);

      if (itSY != integralfamily.symbolicIBP.end()) checkV++;
    }

    if (checkV != theBig)
      return 0;

    checkV = 0;

    for (int g = 0; g < j; g++) { //???

      auto itSY = find(topology[topoName].symbolicIBP.begin(),
                       topology[topoName].symbolicIBP.end(), array3[g]);

      if (itSY != topology[topoName].symbolicIBP.end()) checkV++;
    }

    if (checkV != theBig)
      return 0;
  }

  symmetries halloB;

  halloB.symDOTS = 0;

  for (uint32_t g = 0; g < integralfamily.jule; g++) {
    halloB.ing[g] = -1;
  }
  for (int g = 0; g < j; g++) {
    halloB.ing[array2[g]] = array3[g];
  }

  int nn = integralfamily.loopVar.size();

  //   Check if momenta Flow is available
  for (int g = 0; g < j; g++) {
    if (!((1 << array2[g] & integralfamily.propsMomentaFlowMask) &&
          (1 << array3[g] & topology[topoName].propsMomentaFlowMask))) {
      halloB.det = 1;
      halloB.sector = topology[topoName].mask[j][itt];
      halloB.nOfProps = j;
      halloB.topology = topology[topoName].topology;
      halloB.externalSymmetry = klop;
      halloB.symDOTS = 1;

      symbol jj("placeholder");
      for (int kk = 0; kk < nn; kk++)
        halloB.subst.append(jj == jj);

      if (flag) {
        integralfamily.symVec[integralfamily.mask[j][it]].push_back(halloB);
        integralfamily.skipSector[j].insert(integralfamily.mask[j][it]);
      }
      else {
        integralfamily.relVec[integralfamily.mask[j][it]].push_back(halloB);
      }

      return 1;
    }
  }

  if (integralfamily.symbolicIBP.size() > 0) //???
    halloB.symDOTS = 1;

  // choose linearly independent set of propagators
  int chosenMask = 0;
  vector<int> chosenProps;

  combProps combinatorics(j, nn);

  for (auto itComb : combinatorics.vecComb) {
    lst check1, check2;

    for (int g = 0; g < nn; g++) {
      check1.append(integralfamily.propsMomFlowA[array2[itComb[g]]]/*.subs(get<0>(externalTransf[klop]), subs_options::algebraic).subs(get<1>(externalTransf[klop]), subs_options::algebraic)*/);

      check2.append(
          topology[topoName].propsMomFlowB[klop /*0*/][array3[itComb[g]]].subs(
              topology[topoName].loop2loop));
    }

    if ((testProps(check1,integralfamily.name) != 0) && (testProps(check2,topology[topoName].name) != 0)) {
      for (int g = 0; g < nn; g++) {
        chosenMask |= 1 << itComb[g];
        chosenProps.push_back(itComb[g]);
      }
      break;
    }
  }

  //     Go through the sign ambiguity
  for (int gi = 0; gi < (1 << nn); gi++) {
    lst testsolve;
    for (size_t g = 0; g < chosenProps.size(); g++) {
      if (((1 << g) & gi)) {
        testsolve.append(
            (-integralfamily.propsMomFlowA
                  [array2[chosenProps[g]]]) /*.subs(get<0>(externalTransf[klop]),
                                               subs_options::algebraic).subs(get<1>(externalTransf[klop]),
                                               subs_options::algebraic)*/
            == (topology[topoName]
                    .propsMomFlowB[klop /*0*/][array3[chosenProps[g]]]));
      }
      if (!((1 << g) & gi)) {
        testsolve.append(
            (integralfamily.propsMomFlowA
                 [array2[chosenProps[g]]]) /*.subs(get<0>(externalTransf[klop]),
                                              subs_options::algebraic).subs(get<1>(externalTransf[klop]),
                                              subs_options::algebraic)*/
            == (topology[topoName]
                    .propsMomFlowB[klop /*0*/][array3[chosenProps[g]]]));
      }
    }

    lst sT2;
    sT2.append(lsolve(testsolve, integralfamily.loopVarList)
                   .subs(topology[topoName].loop2loop));



    lst testsolve1, testsolve2;

    for (int g = 0; g < j; g++) {
      if (!((1 << g) & chosenMask)) {

        testsolve1.append(
            (integralfamily
                 .propsMomFlowA[array2[g]] /*.subs(get<0>(externalTransf[klop]),
                                              subs_options::algebraic).subs(get<1>(externalTransf[klop]),
                                              subs_options::algebraic)*/
                 .subs(sT2[0], subs_options::algebraic)));

        testsolve2.append(
            (topology[topoName].propsMomFlowB[klop /*0*/][array3[g]].subs(
                topology[topoName].loop2loop)));

      }
    }

    int countZero = 0;
    int countSign = 0;

    for (size_t itS = 0; itS < testsolve1.nops(); itS++) {
      ex testEx = (testsolve1[itS] - testsolve2[itS]).normal();

      if (!testEx.is_zero()) {
        if ((testEx / 2 - testsolve1[itS]).is_zero()) {
          countSign++;
        }
      }
      else
        countZero++;
    }

    if ((countZero + countSign) == (j - nn)) {

      int re;
      halloB.det = 0;
      halloB.subst.remove_all();

      for (size_t hi = 0; hi < sT2[0].nops(); hi++) {
        halloB.subst.append(sT2[0][hi]);
      }

      if (sT2[0].nops()) {
        re = test(halloB.subst, topoName);

        if (re != 0) {
          halloB.det = re;
          halloB.sector = topology[topoName].mask[j][itt];
          halloB.nOfProps = j;
          halloB.topology = topology[topoName].topology;
          halloB.externalSymmetry = klop;

          if (klop != 0) {
            halloB.symDOTS = 1;
          }

          if (flag) {
            integralfamily.symVec[integralfamily.mask[j][it]].push_back(halloB);
            integralfamily.skipSector[j].insert(integralfamily.mask[j][it]);
          }
          else {
            integralfamily.relVec[integralfamily.mask[j][it]].push_back(halloB);
          }

          return 1;
        }
        else {
          halloB.det = 0;
          halloB.subst.remove_all();
        }
      }

      break;
    }
  }

  halloB.det = 1;
  halloB.sector = topology[topoName].mask[j][itt];
  halloB.nOfProps = j;
  halloB.topology = topology[topoName].topology;
  halloB.externalSymmetry = klop;
  halloB.symDOTS = 1;
  symbol jj("placeholder");
  for (int kk = 0; kk < nn; kk++)
    halloB.subst.append(jj == jj);

  if (flag) {
    integralfamily.symVec[integralfamily.mask[j][it]].push_back(halloB);
    integralfamily.skipSector[j].insert(integralfamily.mask[j][it]);
  }
  else {
    integralfamily.relVec[integralfamily.mask[j][it]].push_back(halloB);
  }

  return 1;
}

int Kira::skip_symmetry_topology(unsigned op) {
  string topoName = collectReductions[op];

  unsigned countMappings = 0;

  for (vector<uint32_t>::iterator topIt = topology[topoName].topLevelSectors.begin();
       topIt != topology[topoName].topLevelSectors.end(); topIt++) {
    auto itF1 = topology[topoName].symmetries.find((*topIt));

    if (itF1 == topology[topoName].symmetries.end()) {
      continue;
    }
    unsigned tokenCount = countMappings;

    for (unsigned topoNumber = 0; topoNumber < op; topoNumber++) {
      string topoNameLower = collectReductions[topoNumber];

      for (vector<uint32_t>::iterator topItLower =
               topology[topoNameLower].topLevelSectors.begin();
           topItLower != topology[topoNameLower].topLevelSectors.end();
           topItLower++) {
        auto itF2 = topology[topoNameLower].symmetries.find((*topItLower));

        if (itF2 == topology[topoNameLower].symmetries.end()) {
          continue;
        }

        auto symmetry1 = (itF1->second).front();
        auto symmetry2 = (itF2->second).front();

        if (get<0>(symmetry1) == get<0>(symmetry2)) {
          countMappings++;
          break;
        }
      }
      if (tokenCount < countMappings)
        break;
    }
  }

  if (topology[topoName].topLevelSectors.size() == countMappings)
    return 0;

  countMappings = 0;

  for (vector<uint32_t>::iterator topIt = integralfamily.topLevelSectors.begin();
       topIt != integralfamily.topLevelSectors.end(); topIt++) {

    auto itF1 = integralfamily.symmetries.find((*topIt));

    if (itF1 == integralfamily.symmetries.end()) {
      continue;
    }

    unsigned tokenCount = countMappings;

    for (unsigned topoNumber = 0; topoNumber < op; topoNumber++) {

      string topoNameLower = collectReductions[topoNumber];

      for (vector<uint32_t>::iterator topItLower =
               topology[topoNameLower].topLevelSectors.begin();
           topItLower != topology[topoNameLower].topLevelSectors.end();
           topItLower++) {

        if (pyred::count_set_bits(*topItLower) > pyred::count_set_bits(*topIt)) {

          int nOfProps = pyred::count_set_bits(*topIt);

          for (size_t it = 0;
               it < topology[topoNameLower].mask[nOfProps].size(); ++it) {
            uint32_t sec = topology[topoNameLower].mask[nOfProps][it];

            if ((sec & (*topItLower)) == sec) {
              auto itF2 = topology[topoNameLower].symmetries.find(sec);

              if (itF2 == topology[topoNameLower].symmetries.end()) {
                continue;
              }

              auto symmetry1 = (itF1->second).front();
              auto symmetry2 = (itF2->second).front();

              if (get<0>(symmetry1) == get<0>(symmetry2)) {
                countMappings++;
                break;
              }
            }
          }
        }
        else {
          auto itF2 = topology[topoNameLower].symmetries.find((*topItLower));

          if (itF2 == topology[topoNameLower].symmetries.end()) {
            continue;
          }

          auto symmetry1 = (itF1->second).front();
          auto symmetry2 = (itF2->second).front();

          if (get<0>(symmetry1) == get<0>(symmetry2)) {
            countMappings++;
            break;
          }
        }
      }
      if (tokenCount < countMappings) break;
    }
  }

  if (topology[topoName].topLevelSectors.size() == countMappings)
    return 0;

  return 1;
}

int Kira::skip_symmetry(string topoName, int j, int it, int itt) {
  auto itSkip = integralfamily.skipSector[j].find(integralfamily.mask[j][it]);

  if (itSkip != integralfamily.skipSector[j].end()) return 1;

  int skipTopSym = 1;

  for (vector<uint32_t>::iterator topIt = topology[topoName].topLevelSectors.begin();
       topIt != topology[topoName].topLevelSectors.end(); topIt++) {
    if ((topology[topoName].mask[j][itt] & (*topIt)) ==
        topology[topoName].mask[j][itt]) {
      skipTopSym = 0;
    }
  }

  if (skipTopSym == 0) {
    skipTopSym = 1;
    for (auto topIt = integralfamily.topLevelSectors.begin();
         topIt != integralfamily.topLevelSectors.end(); topIt++) {
      if ((integralfamily.mask[j][it] & (*topIt)) ==
          integralfamily.mask[j][it]) {
        skipTopSym = 0;
      }
    }
  }

  if (topology[topoName].magic_relations == "true") {
    for (vector<uint32_t>::iterator topIt =
             topology[topoName].topLevelSectors.begin();
         topIt != topology[topoName].topLevelSectors.end(); topIt++) {
      if (((*topIt) & topology[topoName].mask[j][itt])) {
        skipTopSym = 0;
      }
    }
  }

  if (skipTopSym == 1) return 1;

  return 0;
}

int Kira::symmetry_relations(unsigned op) {
  Clock clock;

  if (!skip_symmetry_topology(op)) {
    logger << "Skip mappings because all top level sectors are mapped\n";
    logger << "already to a lower topology.\n";
    logger << "( " << clock.eval_time() << " s )\n";
    return 0;
  }

  string topoName = collectReductions[op];

  for (int j = 0; j < integralfamily.biggestBound + 1; j++) {
    for (size_t itt = 0; itt < topology[topoName].mask[j].size();
         ++itt) { // all possible non trivial sectors

      if (integralfamily.mask[j].size() == 0) continue;

      size_t it2 = 0;

      if (integralfamily.name == topoName)
        it2 = itt;

      for (size_t it = it2; it < integralfamily.mask[j].size(); it++) {

        if (skip_symmetry(topoName, j, it, itt) == 1) {
          continue;
        }

        auto itF1 = integralfamily.symmetries.find(integralfamily.mask[j][it]);

        auto itF2 =
            topology[topoName].symmetries.find(topology[topoName].mask[j][itt]);

        if (!((itF1 != integralfamily.symmetries.end()) &&
              (itF2 != topology[topoName].symmetries.end()))) {
          continue;
        }


        auto symmetry1 = (itF1->second).front();
        auto symmetry2 = (itF2->second).front();

        if (get<0>(symmetry1) != get<0>(symmetry2)) {
          continue;
        }

        if (integralfamily.name == topoName) {
          if (integralfamily.mask[j][it] != topology[topoName].mask[j][itt]) {
            for (size_t permIt = 0; permIt < get<1>(symmetry2).size();
                 permIt++) {
              symmetry_finder(it, itt, topoName, j, get<1>(symmetry1).front(),
                              get<1>(symmetry2)[permIt], 1, get<3>(symmetry2));
            }
          }
          else if (integralfamily.mask[j][it] == integralfamily.mask[j][itt]) {
            for (size_t permIt = 1; permIt < get<1>(symmetry2).size();
                 permIt++) {
              symmetry_finder(it, itt, topoName, j, get<1>(symmetry1).front(),
                              get<1>(symmetry2)[permIt], 0, get<3>(symmetry2));
            }
          }
        }
        else {
          for (size_t permIt = 0; permIt < get<1>(symmetry2).size(); permIt++) {
            symmetry_finder(it, itt, topoName, j, get<1>(symmetry1).front(),
                            get<1>(symmetry2)[permIt], 1, get<3>(symmetry2));
          }
        }
      }
    }

    load_bar(j + 1, integralfamily.biggestBound + 1, 50, 100);
  }

  int countR = 0;
  for (uint32_t it = 0; it < (1u << integralfamily.jule) + 1; it++) {
    sort(integralfamily.symVec[it].begin(), integralfamily.symVec[it].end(),
         myfunction);
    countR += integralfamily.symVec[it].size();
  }

  int countS = 0;
  for (uint32_t it = 0; it < (1u << integralfamily.jule) + 1; it++) {
    sort(integralfamily.relVec[it].begin(), integralfamily.relVec[it].end(),
         myfunction);
    countS += integralfamily.relVec[it].size();
  }

  write_symmetries((outputDir + "/sectormappings/" + integralfamily.name +
                    "/sectorRelations"),
                   integralfamily.symVec);
  write_symmetries((outputDir + "/sectormappings/" + integralfamily.name +
                    "/sectorSymmetries"),
                   integralfamily.relVec);

  logger << "\nSector relations: " << countR << "\n";
  logger << "Sector symmetries: " << countS << "\n";
  logger << "( " << clock.eval_time() << " s )\n";

  return 1;
}

void Kira::write_symmetries(const string otputName, SYM symVec[]) {
  ofstream fileSymVec;

  string otputNametmp = otputName + ".back";

  fileSymVec.open(otputNametmp.c_str());

  for (uint32_t i = 0; i < (1u << integralfamily.jule) + 1; i++) {
    for (unsigned it = 0; it < symVec[i].size(); it++) {
      symmetries halloB = symVec[i][it];
      fileSymVec << i << " ";

      for (unsigned itt = 0; itt < halloB.subst.nops(); itt++) {
        fileSymVec << halloB.subst.op(itt).lhs() << " "
                   << halloB.subst.op(itt).rhs() << " ";
      }
      fileSymVec << halloB.det << " ";
      fileSymVec << halloB.sector << " ";
      fileSymVec << halloB.nOfProps << " ";
      fileSymVec << halloB.externalSymmetry << " ";
      fileSymVec << get<0>(externalTransf[halloB.externalSymmetry]) << " ";
      fileSymVec << get<1>(externalTransf[halloB.externalSymmetry]) << " ";
      for (uint32_t g = 0; g < integralfamily.jule; g++) {
        fileSymVec << halloB.ing[g] << " ";
      }
      fileSymVec << halloB.symDOTS << " ";
      fileSymVec << halloB.topology << endl;
    }
  }

  fileSymVec.close();

  rename(otputNametmp.c_str(), otputName.c_str());

  remove(otputNametmp.c_str());
}

int Kira::read_symmetries(const string inputName, SYM symVec[]) {
  ifstream input;
  int count = 0;
  if (file_exists(inputName.c_str())) {
    input.open(inputName.c_str());
    logger << "read from: " << inputName << "\n";
    parser symbolReader(GiNaCSymbols);
    while (1) {
      symmetries halloB;

      int sector;
      string skipRead;

      if (!(input >> sector)) break;
      for (unsigned it = 0; it < integralfamily.loopVar.size(); it++) {
        string substleft;
        string substright;

        if (!(input >> substleft)) break;
        if (!(input >> substright)) break;
        halloB.subst.append(symbolReader(substleft) ==
                            symbolReader(substright));
      }
      if (!(input >> halloB.det)) break;
      if (!(input >> halloB.sector)) break;
      if (!(input >> halloB.nOfProps)) break;
      if (!(input >> halloB.externalSymmetry)) break;
      if (!(input >> skipRead)) break;
      if (!(input >> skipRead)) break;
      for (uint32_t g = 0; g < integralfamily.jule; g++) {
        if (!(input >> halloB.ing[g])) break;
      }
      if (!(input >> halloB.symDOTS)) break;
      if (!(input >> halloB.topology)) break;

      symVec[sector].push_back(halloB);

      if (sector != halloB.sector) {
        integralfamily.skipSector[halloB.nOfProps].insert(sector);
      }
      count++;
    }
    logger << "Number of symmetries read out: " << count << "\n";
    input.close();

    return 1;
  }
  return 0;
}
